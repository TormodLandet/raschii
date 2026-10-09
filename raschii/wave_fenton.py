import math

from numpy import (
    arange,
    array,
    asarray,
    atleast_1d,
    broadcast_arrays,
    cos,
    cosh,
    isfinite,
    linspace,
    newaxis,
    pi,
    sin,
    sinh,
    stack,
    sum,
    zeros,
)
from numpy.linalg import solve
from numpy.typing import NDArray

from .base_classes import AirPhaseModel, WaveModel
from .common import (
    Frame,
    NonConvergenceError,
    RaschiiError,
    blend_air_and_wave_velocities,
    cosh_by_cosh,
    cosh_ratio,
    sinh_by_cosh,
    trapezoid_integration,
)
from .cpp import FentonCppGenerator


class FentonWave(WaveModel):
    required_input = {"height", "depth", "length", "N"}
    optional_input = {"air": None, "g": 9.81, "relax": 0.5}

    def __init__(
        self,
        height: float,
        depth: float,
        length: float | None = None,
        *,
        N: int = 5,
        period: float | None = None,
        air: AirPhaseModel | None = None,
        g: float = 9.81,
        relax: float = 0.5,
        maxiter: int = 500,
        num_steps: int | None = None,
        current: float = 0.0,
        current_criterion: str = "euler",
    ):
        """
        Implement stream function waves based on the paper by Rienecker and
        Fenton (1981)

        * height: wave height above still water level
        * depth: still water distance from the flat sea bottom to the free surface
          in meters, but you can give -1.0 for infinite depth
        * length: the periodic length of the wave (optional, if not given then period is used)
        * N: the number of coefficients in the truncated Fourier series
        * period: the wave period (optional, if not given then length is used).
          This is the period observed at a fixed point, so a current Doppler shifts it
        * current: the mean current in m/s, default 0. The wave shape is found in the frame
          where the Eulerian mean current is zero (as in Fenton's Fourier program) and the
          wave field is then moved by the current. See ``current_criterion``
        * current_criterion: how ``current`` is defined, following Fenton's Fourier program:
          ``"euler"`` (default) is the mean velocity at a fixed point below the wave troughs,
          ``"stokes"`` is the mean mass transport velocity (Eulerian current plus the
          volume flux due to the waves divided by the depth)
        """
        current_criterion = _check_current_criterion(current_criterion)
        if current != 0.0 and air is not None:
            raise RaschiiError("Air-phase models do not support a current")

        if length is None:
            if period is None:
                raise RaschiiError("Either length or period must be given, both are None!")
            length = compute_length_from_period(
                height=height,
                depth=depth,
                period=period,
                N=N,
                g=g,
                relax=relax,
                current=current,
                current_criterion=current_criterion,
            )

        #: The mean current in [m/s] as specified, see :attr:`current_criterion`
        self.current_input: float = float(current)

        #: How the current is specified: "euler" or "stokes"
        self.current_criterion: str = current_criterion

        #: The wave height
        self.height: float = height

        #: The water depth
        self.depth: float = depth

        #: The wave length
        self.length: float = length

        #: The approximation order
        self.order: int = N

        #: The optional air-phase model
        self.air: AirPhaseModel | None = air

        #: The acceleration of gravity
        self.g: float = g

        #: The numerical relaxation in the optimization loop
        self.relax: float = relax

        #: Warnings raised when generating this wave
        self.warnings: str = ""

        if N < 1:
            self.warnings = "Fenton order must be at least 1, using order 1"
            self.order = 1

        # Find the coeffients through optimization
        data = fenton_coefficients(
            height, depth, length, N, g, relax=relax, maxiter=maxiter, num_steps=num_steps
        )
        self.set_data(data)

        # For evaluating velocities close to the free surface
        self.eta_eps: float = self.height / 1e5

        # Provide velocities also in the air phase
        if self.air is not None:
            self.air.set_wave(self)

        # Niche use use case: Provide a C++ code generator for this wave model
        self.cpp = FentonCppGenerator(self)

    def set_data(self, data):
        """
        Update the coefficients defining this stream-function wave
        """
        self.data = data

        #: Wave elevation at colocation points
        self.eta = data["eta"]

        #: Positions of colocation points
        self.x = data["x"]

        #: Wave number (2 pi / wavelength) in [1/m]
        self.k = data["k"]

        #: Wave celerity in [m/s] relative to the water, i.e., in the frame where the
        #: Eulerian mean current is zero
        self.c = data["c"]

        # Cosine series coefficients for the elevation
        N = len(self.eta) - 1
        self.E = zeros(N + 1, float)
        J = arange(0, N + 1)
        self.E = trapezoid_integration(self.eta * cos(J * J[:, newaxis] * pi / N))

        #: The value *Q* from the Fenton stream-function wave
        self.Q = data["Q"]

        # Volume flux due to the waves q = c d - Q, where the mean mass transport
        # (Stokes) current is the Eulerian current plus q / d, see Fenton's Fourier program
        d = 25 * self.length if self.depth < 0 else self.depth
        q = self.c * d - self.Q

        #: The mean Eulerian current in [m/s]
        self.current: float = self.current_input
        if self.current_criterion == "stokes":
            self.current = self.current_input - q / d

        #: The mean mass transport (Stokes) current in [m/s]
        self.current_stokes: float = self.current + q / d

        #: Wave phase speed in [m/s] as seen from a fixed point
        self.phase_speed: float = self.c + self.current
        if self.phase_speed <= 0.0:
            raise RaschiiError("The current is too strong against the waves, they cannot propagate")

        #: Wave period in [s] observed at a fixed point
        self.period = self.length / self.phase_speed

        #: Wave frequency in [rad/s] observed at a fixed point
        self.omega = self.phase_speed * self.k

    def stream_function(self, x, z, t=0, frame=Frame.EARTH):
        """
        Compute the stream function at time t for position(s) x.

        * frame: :class:`~raschii.Frame` – ``Frame.EARTH`` (default) includes
          the constant base-flow term; ``Frame.WAVE`` returns only the
          oscillatory part.
        """
        if isinstance(x, (float, int)):
            x, z = [x], [z]
        x2 = asarray(x, dtype=float) - self.phase_speed * t
        z2 = asarray(z, dtype=float)
        x2, z2 = x2[:, newaxis], z2[:, newaxis]

        N = len(self.eta) - 1
        B = self.data["B"]
        k = self.k
        J = arange(1, N + 1)

        psi = (sinh(J * k * z2) / cosh(J * k * self.depth) * cos(J * k * x2)).dot(B[1:])

        if frame == Frame.EARTH:
            return B[0] * z + psi
        elif frame == Frame.WAVE:
            return psi
        else:
            raise ValueError(f"Unknown frame {frame!r}; use Frame.EARTH or Frame.WAVE")

    def _surface_elevation(self, x: NDArray, t: NDArray, include_depth: bool) -> NDArray:
        N = len(self.eta) - 1
        J = arange(0, N + 1)
        k, c = self.k, self.phase_speed
        x_r = x[newaxis, :, newaxis]  # (1, n_x, 1)
        t_r = t[:, newaxis, newaxis]  # (T, 1, 1)
        J_r = J[newaxis, newaxis, :]  # (1, 1, N+1)
        eta = 2 * trapezoid_integration(self.E * cos(J_r * k * (x_r - c * t_r)), axis=-1) / N
        # eta shape: (T, n_x)
        if include_depth:
            if self.depth < 0:
                raise RaschiiError("Cannot include depth in elevation for infinite depth")
            subtract = 0.0
        else:
            subtract = 25 * self.length if self.depth < 0 else self.depth
        return eta - subtract

    def surface_slope(self, x, t=0):
        """
        Compute the x derivative of the surface elevation at time t
        """
        if isinstance(x, (float, int)):
            x = array([x], float)
        x = asarray(x)

        # Cosine transformation of the elevation
        N = len(self.eta) - 1
        J = arange(0, N + 1)
        k, c = self.k, self.phase_speed
        return -2 * trapezoid_integration(self.E * J * k * sin(J * k * (x[:, newaxis] - c * t))) / N

    def _velocity(self, x: NDArray, z: NDArray, t: NDArray, all_points_wet: bool) -> NDArray:
        x_1d, z_1d, t_1d = x, z, t  # save (N,), (N,), (T,) before reshaping
        x = x[newaxis, :]  # (1, n_points)
        z = z[newaxis, :]  # (1, n_points)
        t = t[:, newaxis]  # (n_times, 1)

        N = len(self.eta) - 1
        B = self.data["B"]
        k = self.k
        c = self.phase_speed
        depth = self.depth
        J = arange(1, N + 1)

        phase = J * k * (x - c * t)[..., newaxis]  # (n_times, n_points, N)
        kz = J * k * z[..., newaxis]  # (1, n_points, N)
        kh = J * k * depth  # (N,)

        vel_x = k * sum(J * B[1:] * cos(phase) * cosh(kz) / cosh(kh), axis=-1)
        vel_z = k * sum(J * B[1:] * sin(phase) * sinh(kz) / cosh(kh), axis=-1)
        vel_x = vel_x + self.current  # the current moves the whole wave field
        vel = stack([vel_x, vel_z], axis=-1)  # (n_times, n_points, 2)

        if not all_points_wet:
            for i, ti in enumerate(t_1d):
                blend_air_and_wave_velocities(
                    x_1d, z_1d, float(ti), self, self.air, vel[i], self.eta_eps
                )

        return vel

    def acceleration(
        self,
        x: float | NDArray,
        z: float | NDArray,
        t: float | NDArray = 0.0,
        all_points_wet: bool = False,
    ):
        """
        Compute the horizontal and vertical fluid acceleration at each position
        ``(x, z)`` at each time ``t``.

        .. note::

           This method is **water-phase only**.  Acceleration above the free
           surface is set to zero when ``all_points_wet=False`` (the default).
           Air-phase blending for accelerations is not yet implemented.  If an
           air model is attached and you query points above the free surface with
           ``all_points_wet=False``, a :exc:`~raschii.RaschiiError` is raised.

        Parameters
        ----------
        x : float | array
            Horizontal position(s).
        z : float | array
            Vertical position(s) where z = 0 at the sea floor and
            z = depth at the free surface.
        t : float | array, optional
            Time(s) at which to compute the acceleration (default 0).
        all_points_wet : bool, optional
            If ``True``, evaluate the wave formula at all points regardless of
            whether they are above the free surface.  Useful for testing.

        Returns
        -------
        ndarray
            Shape ``(2,)`` for scalar inputs, ``(N, 2)`` for array points and
            scalar time, ``(T, 2)`` for scalar point and array time,
            ``(T, N, 2)`` for array points and array time.
        """
        time_is_scalar = asarray(t).ndim == 0
        point_is_scalar = asarray(x).ndim == 0 and asarray(z).ndim == 0

        x_1d = atleast_1d(asarray(x, dtype=float))
        z_1d = atleast_1d(asarray(z, dtype=float))
        t_1d = atleast_1d(asarray(t, dtype=float))
        x_1d, z_1d = broadcast_arrays(x_1d, z_1d)

        x = x_1d[newaxis, :]  # shape (1, n_points)
        z = z_1d[newaxis, :]  # shape (1, n_points)
        t = t_1d[:, newaxis]  # shape (n_times, 1)

        N = len(self.eta) - 1
        B = self.data["B"]
        k = self.k
        c = self.phase_speed
        depth = self.depth
        J = arange(1, N + 1)

        phase = J * k * (x - c * t)[..., newaxis]  # shape (n_times, n_points, N)
        kz = J * k * z[..., newaxis]  # shape (1, n_points, N)
        kh = J * k * depth  # shape (N,)

        acc_x = k * sum(J * B[1:] * J * k * c * sin(phase) * cosh(kz) / cosh(kh), axis=-1)
        acc_z = k * sum(J * B[1:] * J * k * -c * cos(phase) * sinh(kz) / cosh(kh), axis=-1)

        acc = stack([acc_x, acc_z], axis=-1)  # shape (n_times, n_points, 2)

        if not all_points_wet:
            if self.air is not None:
                eta = asarray(
                    [self.surface_elevation(x_1d, float(ti)) for ti in t_1d], dtype=float
                )  # (n_times, n_points)
                above = z_1d[newaxis, :] > eta + self.eta_eps  # (n_times, n_points)
                if above.any():
                    raise RaschiiError(
                        "acceleration() does not support air-phase blending. "
                        "Points above the free surface were detected. "
                        "Use all_points_wet=True to evaluate the raw wave formula, "
                        "or only query points below the free surface."
                    )
            else:
                # Zero out accelerations above the free surface (no air model)
                for i, ti in enumerate(t_1d):
                    eta_i = asarray(self.surface_elevation(x_1d, float(ti)), dtype=float)
                    above_i = z_1d > eta_i + self.eta_eps
                    if above_i.any():
                        acc[i, above_i] = 0.0

        if time_is_scalar and point_is_scalar:
            return acc.squeeze()  # shape (2,)
        elif time_is_scalar:
            return acc[0]  # shape (n_points, 2)
        elif point_is_scalar:
            return acc[:, 0]  # shape (n_times, 2)

        return acc  # shape (n_times, n_points, 2)

    def _velocity_potential(self, x: NDArray, z: NDArray, t: NDArray) -> NDArray:
        depth = 25.0 * self.length if self.depth < 0 else self.depth
        N = len(self.eta) - 1
        B = self.data["B"]
        k = self.k
        c = self.phase_speed
        J = arange(1, N + 1)
        x_r = x[newaxis, :, newaxis]  # (1, M, 1)
        t_r = t[:, newaxis, newaxis]  # (T, 1, 1)
        z_r = z[newaxis, :, newaxis]  # (1, M, 1)
        J_r = J[newaxis, newaxis, :]  # (1, 1, N)
        phi = (
            B[1:] * sin(J_r * k * (x_r - c * t_r)) * cosh_ratio(J_r * k * z_r, J_r * k * depth)
        ).sum(axis=-1)  # (T, M)
        return phi + self.current * x[newaxis, :]

    def write_swd(self, path, dt, tmax=None, nperiods=None, amp: int = 1):
        """
        Write a SWD-file of the wave field.

        * path:     Full path of the new SWD file
        * dt:       The temporal sampling spacing in the SWD file
        * tmax:     The temporal sampling range in the SWD file is [0, tmax]
        * nperiods: Alternative specification: tmax = nperiods * wave_period
        * amp:      SWD amp flag (1, 2, or 3).  Default is 1.

                    - 1: store potential coefficients at z=0 (calm surface)
                    - 2: store potential coefficients on the actual wavy free surface
                    - 3: store elevation only (no potential data)

        See the SWD documentation for the details:

        * https://spectral-wave-data.readthedocs.io/en/latest/shape_2.html
        * https://spectral-wave-data.readthedocs.io/en/latest/swd_format.html
        """
        if self.current != 0.0:
            raise RaschiiError("Writing SWD files is not supported for waves with a current")

        from .swd import SwdWriterFenton

        SwdWriterFenton(self).write(path, dt, tmax=tmax, nperiods=nperiods, amp=amp)


def fenton_coefficients(
    height, depth, length, N, g=9.8, maxiter=500, tolerance=1e-8, relax=1.0, num_steps=None
):
    """
    Find B, Q and R by Newton-Raphson following Rienecker and Fenton (1981)

    Using relaxation can help in some difficult cases, try a value less than 1
    to decrease convergence speed, but increase chances of converging.
    """
    if depth < 0:
        depth = 25 * length

    # Non dimensionalised input
    H = height / depth
    lam = length / depth
    k = 2 * pi / lam
    c = (math.tanh(k) / k) ** 0.5
    D = 1
    N_unknowns = 2 * (N + 1) + 2

    # Input data arrays
    J = arange(1, N + 1)
    M = arange(0, N + 1)
    x = M * lam / (2 * N)

    def initial_guess(H):
        """
        Initial guesses for the unknowns (linear wave)
        """
        B = zeros(N + 1, float)
        B[0] = c
        B[1] = -H / (4 * c * k)
        eta = 1 + H / 2 * cos(k * x)
        Q = c
        R = 1 + 0.5 * c**2
        return B, Q, R, eta

    def optimize(B, Q, R, eta, H):
        """
        Find B, Q and R by Newton iterations starting from the given initial
        guesses. According to Rienecker and Fenton (1981) a linear theory
        initial guess should work unless H close to breaking, then an initial
        guess from the optimization routine run with a slightly lower H should
        be used instead.
        """
        # Insert initial guesses into coefficient vector
        coeffs = zeros(N_unknowns, float)
        coeffs[: N + 1] = B
        coeffs[N + 1 : 2 * N + 2] = eta
        coeffs[2 * N + 2] = Q
        coeffs[2 * N + 3] = R
        f = func(coeffs, H, k, D, J, M)

        for it in range(1, maxiter + 1):
            jac = fprime(coeffs, H, k, D, J, M)
            delta = solve(jac, -f)
            coeffs += delta * relax
            f = func(coeffs, H, k, D, J, M)

            # Check the progress
            error = abs(f).max()
            eta_max = coeffs[N + 1 : 2 * N + 2].max()
            eta_min = coeffs[N + 1 : 2 * N + 2].min()
            if eta_max > 2:
                raise NonConvergenceError(
                    "Optimization did not converge. Got "
                    "max(eta)/depth = %r in iteration %d" % (eta_max, it)
                )
            elif eta_min < 0.1:
                raise NonConvergenceError(
                    "Optimization did not converge. Got "
                    "min(eta)/depth = %r in iteration %d" % (eta_min, it)
                )
            elif not isfinite(error):
                raise NonConvergenceError(
                    "Optimization did not converge. Got error %r in iteration %d" % (error, it)
                )
            elif error < tolerance:
                B = coeffs[: N + 1]
                eta = coeffs[N + 1 : 2 * N + 2]
                Q = coeffs[2 * N + 2]
                R = coeffs[2 * N + 3]
                return B, Q, R, eta, error, it
        raise NonConvergenceError(
            "Optimization did not converge after %d iterations, error = %r" % (it, error)
        )

    # Perform the optimization, optionally in steps gradually increasing H
    steps = wave_height_steps(num_steps, D, lam, H)
    B, Q, R, eta = initial_guess(steps[0])
    for Hi in steps:
        B, Q, R, eta, error, niter = optimize(B, Q, R, eta, Hi)

    # Scale back to physical space
    B[0] *= (g * depth) ** 0.5
    B[1:] *= (g * depth**3) ** 0.5
    return {
        "x": x * depth,
        "eta": eta * depth,
        "B": B,
        "Q": Q * (g * depth**3) ** 0.5,
        "R": R * g * depth,
        "k": k / depth,
        "c": B[0],
        "error": error,
        "niter": niter,
    }


def wave_height_steps(num_steps, D, lam, H):
    """
    Compute the breaking height and use this to select how many steps take when
    gradually increasing the wave height to improve convergence on high waves
    """
    # Breaking height
    Hb = 0.142 * math.tanh(2 * pi * D / lam) * lam

    # Try with progressively higher waves to get better initial conditions
    if num_steps is not None:
        pass
    if H > 0.75 * Hb:
        num_steps = 10
    elif H > 0.65 * Hb:
        num_steps = 5
    else:
        num_steps = 3

    if num_steps == 1:
        return [H]
    else:
        return linspace(H / num_steps, H, num_steps)


def func(coeffs, H, k, D, J, M):
    "The function to minimize"
    N_unknowns = coeffs.size
    N = J.size

    B0 = coeffs[0]
    B = coeffs[1 : N + 1]
    eta = coeffs[N + 1 : 2 * N + 2]
    Q = coeffs[2 * N + 2]
    R = coeffs[2 * N + 3]

    # The function to me minimized
    f = zeros(N_unknowns, float)

    # Loop over the N + 1 points along the half wave
    for m in M:
        S1 = sinh_by_cosh(J * k * eta[m], J * k * D)
        C1 = cosh_by_cosh(J * k * eta[m], J * k * D)
        S2 = sin(J * m * pi / N)
        C2 = cos(J * m * pi / N)

        # Velocity at the free surface
        # The sign of B0 is swapped from what is in the paper
        um = -B0 + k * J.dot(B * C1 * C2)
        vm = 0 + k * J.dot(B * S1 * S2)

        # Enforce a streamline along the free surface
        # The sign of B0 is swapped from what is in the paper
        f[m] = -B0 * eta[m] + B.dot(S1 * C2) + Q

        # Enforce the dynamic free surface boundary condition
        f[N + 1 + m] = (um**2 + vm**2) / 2 + eta[m] - R

    # Enforce mean(eta) = D
    f[-2] = trapezoid_integration(eta) / N - 1

    # Enforce eta_0 - eta_N = H, the wave height criterion
    f[-1] = eta[0] - eta[-1] - H

    return f


def fprime_num(coeffs, H, k, D, J, M):
    "The Jacobian of the function to minimize (numerical version)"
    N_unknowns = coeffs.size
    dc = 1e-10
    jac = zeros((N_unknowns, N_unknowns), float)
    f0 = func(coeffs, H, k, D, J, M)
    for i in range(N_unknowns):
        incr = zeros(N_unknowns, float)
        incr[i] = dc
        f1 = func(coeffs + incr, H, k, D, J, M)
        jac[:, i] = (f1 - f0) / dc
    return jac


def fprime(coeffs, H, k, D, J, M):
    "The Jacobian of the function to minimize"
    N_unknowns = coeffs.size
    N = J.size

    jac = zeros((N_unknowns, N_unknowns), float)
    B0 = coeffs[0]
    B = coeffs[1 : N + 1]
    eta = coeffs[N + 1 : 2 * N + 2]

    for m in range(N + 1):
        S1 = sinh_by_cosh(J * k * eta[m], J * k * D)
        C1 = cosh_by_cosh(J * k * eta[m], J * k * D)
        S2 = sin(J * m * pi / N)
        C2 = cos(J * m * pi / N)

        SC = S1 * C2
        SS = S1 * S2
        CC = C1 * C2
        CS = C1 * S2

        # Velocity at the free surface
        um = -B0 + k * J.dot(B * CC)
        vm = 0 + k * J.dot(B * SS)

        # Derivatives of the eq. for the streamline along the free surface
        jac[m, N + 1 + m] = um
        jac[0 : N + 1, 0] = -eta
        jac[m, 1 : N + 1] = SC
        jac[m, -2] = 1

        # Derivatives of the dynamic free surface boundary condition
        jac[N + 1 + m, N + 1 + m] = 1 + (
            um * k**2 * B.dot(J**2 * SC) + vm * k**2 * B.dot(J**2 * CS)
        )
        jac[N + 1 + m, -1] = -1
        jac[N + 1 + m, 0] = -um
        jac[N + 1 + m, 1 : N + 1] = k * um * J * CC + k * vm * J * SS

    # Derivative of mean(eta) = 1
    jac[-2, N + 1 : 2 * N + 2] = M * 0 + 1 / N
    jac[-2, N + 1] = 1 / (2 * N)
    jac[-2, 2 * N + 1] = 1 / (2 * N)

    # Derivative of the wave height criterion
    jac[-1, N + 1] = 1
    jac[-1, 2 * N + 1] = -1

    return jac


def _check_current_criterion(current_criterion: str) -> str:
    crit = str(current_criterion).lower()
    if crit not in ("euler", "stokes"):
        raise RaschiiError(f"current_criterion must be 'euler' or 'stokes', not {current_criterion!r}")
    return crit


def compute_length_from_period(
    height: float,
    depth: float,
    period: float,
    N: int = 5,
    g: float = 9.81,
    relax: float = 0.5,
    current: float = 0.0,
    current_criterion: str = "euler",
):
    """
    Compute the wave length from the wave period using the Fenton wave theory

    The period is the one observed at a fixed point, so a current gives a Doppler
    shift and a different wave length than without a current.

    This would be much faster if we had an implementation of the Fenton wave
    theory dispersion relation for arbitrary order N
    """
    from .wave_airy import compute_length_from_period as airy_compute_length_from_period

    # Initial guess is based on the linear dispersion relation for deep water waves
    length = airy_compute_length_from_period(depth=depth, period=period, g=g)
    if current != 0.0:
        # The linear dispersion relation with a Doppler shifted period: L = (c(L) + U) T
        for _ in range(50):
            k = 2 * pi / length
            c_airy = math.sqrt(g / k * (1.0 if depth < 0 else math.tanh(k * depth)))
            length_new = (c_airy + current) * period
            if length_new <= 0.0:
                raise RaschiiError("The current is too strong against the waves for this period")
            if abs(length_new - length) < 1e-6 * length:
                break
            length = 0.5 * (length + length_new)

    # Find the length by Newton iterations
    kw = dict(
        height=height,
        depth=depth,
        N=N,
        g=g,
        relax=relax,
        current=current,
        current_criterion=current_criterion,
    )
    wave1 = FentonWave(length=length * 0.95, **kw)
    wave2 = FentonWave(length=length * 1.05, **kw)

    length_N = 0.0
    iter = 0
    while abs(length_N - length) > 1e-4:
        # Store the previous length
        length = length_N

        # New guess for the wave length by interpolation
        f = (period - wave1.period) / (wave2.period - wave1.period)
        length_N = wave1.length + (wave2.length - wave1.length) * f

        # Resulting wave period for the new length from the dispersion relation
        waveN = FentonWave(length=length_N, **kw)

        # Update the two points used for the interpolation in the next iteration
        if waveN.period < period:
            wave1 = waveN
        else:
            wave2 = waveN

        iter += 1
        if iter > 100:
            raise NonConvergenceError(
                "Failed to converge when computing wave length from period for Fenton waves"
            )

    return length_N
