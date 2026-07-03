from numpy import pi, cos, sin, zeros, array, asarray, sinh, cosh, tanh
from .common import (
    blend_air_and_wave_velocities,
    blend_air_and_wave_velocity_cpp,
    cosh_ratio,
    np2py,
    RasciiError,
    NonConvergenceError,
)
from .base_classes import WaveModel


class AiryWave(WaveModel):
    required_input = ("height", "depth", "length")
    optional_input = {"air": None, "g": 9.81}

    def __init__(
        self,
        height: float,
        depth: float,
        length: float | None = None,
        period: float | None = None,
        air=None,
        g: float = 9.81,
    ):
        """
        Linear Airy waves

        * height: wave height above still water level
        * depth: still water distance from the flat sea bottom to the surface
          in meters, but you can give -1.0 for infinite depth
        * length: the periodic length of the wave (optional, if not given then period is used)
        * period: the wave period (optional, if not given then length is used)
        """
        if length is None:
            if period is None:
                raise RasciiError("Either length or period must be given, both are None!")
            length = compute_length_from_period(depth=depth, period=period, g=g)

        self.height: float = height  #: The wave height
        self.depth: float = depth  #: The water depth
        self.length: float = length  #: The wave length
        self.air = air  #: The optional air-phase model
        self.g: float = g  #: The acceleration of gravity
        self.warnings: str = ""  #: Warnings raised when generating this wave

        self.k = 2 * pi / length
        if self.depth < 0:
            # Infinite depth
            self.omega = (self.k * g) ** 0.5
        else:
            # Finite depth
            self.omega = (self.k * g * tanh(self.k * depth)) ** 0.5
        self.c = self.omega / self.k
        self.T = self.length / self.c  # Wave period

        # For evaluating velocities close to the free surface
        self.eta_eps = self.height / 1e5

        # Provide velocities also in the air phase
        if self.air is not None:
            self.air.set_wave(self)

    def surface_elevation(self, x: float | list[float], t: float = 0.0, include_depth: bool = True):
        """
        Compute the surface elavation at time t for position(s) x
        """
        if isinstance(x, (float, int)):
            x = array([x], float)
        x = asarray(x)

        if include_depth:
            if self.depth < 0:
                raise RasciiError("Cannot include depth in elevation for infinite depth")
            offset = self.depth
        else:
            offset = 0.0

        return offset + self.height / 2 * cos(self.k * x - self.omega * t)

    def velocity(
        self,
        x: float | list[float],
        z: float | list[float],
        t: float = 0,
        all_points_wet: bool = False,
    ):
        """
        Compute the fluid velocity at time t for position(s) (x, z)
        where z is 0 at the bottom and equal to depth at the free surface
        """
        if self.depth < 0:
            raise RasciiError("Cannot currently compute velocity for infinite depth waves")
        if isinstance(x, (float, int)):
            x, z = [x], [z]
        x = asarray(x, dtype=float)
        z = asarray(z, dtype=float)

        H = self.height
        k = self.k
        d = self.depth
        w = self.omega

        vel = zeros((x.size, 2), float) + 1
        vel[:, 0] = w * H / 2 * cosh(k * z) / sinh(k * d) * cos(k * x - w * t)
        vel[:, 1] = w * H / 2 * sinh(k * z) / sinh(k * d) * sin(k * x - w * t)

        if not all_points_wet:
            blend_air_and_wave_velocities(x, z, t, self, self.air, vel, self.eta_eps)

        return vel

    def elevation_cpp(self):
        """
        Return C++ code for evaluating the elevation of this specific wave.
        The positive traveling direction is x[0]
        """
        # Repr of np.float64(42.0) is "np.float64(42.0)" and not "42.0"
        # We use repr to make Python output a "smart" amount of digits
        depth = np2py(self.depth)
        height = np2py(self.height)
        k = np2py(self.k)
        c = np2py(self.c)

        return f"{depth!r} + {height!r} / 2.0 * cos({k!r} * (x[0] - {c!r} * t))"

    def velocity_cpp(self, all_points_wet=False):
        """
        Return C++ code for evaluating the particle velocities of this specific
        wave. Returns the x and z components only with z positive upwards. The
        positive traveling direction is x[0] and the vertical coordinate is x[2]
        which is zero at the bottom and equal to +depth at the mean water level.
        """
        # Repr of np.float64(42.0) is "np.float64(42.0)" and not "42.0"
        # We use repr to make Python output a "smart" amount of digits
        H = float(self.height)
        k = float(self.k)
        d = float(self.depth)
        w = float(self.omega)
        a = float(w * H / (2 * sinh(k * d)))

        cpp_x = f"{a!r} * cosh({k!r} * x[2]) * cos({k!r} * x[0] - {w!r} * t)"
        cpp_z = f"{a!r} * sinh({k!r} * x[2]) * sin({k!r} * x[0] - {w!r} * t)"

        if all_points_wet:
            return cpp_x, cpp_z

        # Handle velocities above the free surface
        e_cpp = self.elevation_cpp()
        cpp_ax = cpp_az = None
        cpp_psiw = cpp_psia = cpp_slope = None
        if self.air is not None:
            cpp_ax, cpp_az = self.air.velocity_cpp()
            cpp_psiw = self.stream_function_cpp(frame="c")
            cpp_psia = self.air.stream_function_cpp(frame="c")
            cpp_slope = self.slope_cpp()

        cpp_x = blend_air_and_wave_velocity_cpp(
            cpp_x, cpp_ax, e_cpp, "x", self.eta_eps, self.air, cpp_psiw, cpp_psia, cpp_slope
        )
        cpp_z = blend_air_and_wave_velocity_cpp(
            cpp_z, cpp_az, e_cpp, "z", self.eta_eps, self.air, cpp_psiw, cpp_psia, cpp_slope
        )

        return cpp_x, cpp_z

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

        Implemented via a Stokes N=1 wave (analytically identical to Airy theory).

        See the SWD documentation for the details:

        * https://spectral-wave-data.readthedocs.io/en/latest/shape_2.html
        * https://spectral-wave-data.readthedocs.io/en/latest/swd_format.html
        """
        from .swd import SwdWriterAiry

        SwdWriterAiry(self).write(path, dt, tmax=tmax, nperiods=nperiods, amp=amp)

    def velocity_potential(
        self,
        x: float | list[float],
        z: float | list[float],
        t: float = 0,
    ):
        """
        Compute the earth-frame velocity potential φ at time t for position(s) (x, z).

        z is measured from the sea floor (z=0 at bottom, z≈depth at calm surface).
        The gradient ∇φ equals the oscillatory fluid velocity as returned by
        :meth:`velocity`; the formula is

        .. math::

           \\phi(x, z, t) = \\frac{\\omega H}{2k}
               \\frac{\\cosh(kz)}{\\sinh(kd)}\\sin(kx - \\omega t)

        Works for both finite and infinite depth.  For infinite-depth waves
        (``depth=-1``) an effective depth of 25 × wave_length is used internally;
        z should be supplied in the same coordinate system (z=0 at the effective
        sea floor).
        """
        if isinstance(x, (float, int)):
            x, z = [x], [z]
        x = asarray(x, dtype=float)
        z = asarray(z, dtype=float)

        H = self.height
        k = self.k
        w = self.omega

        d = 25.0 * self.length if self.depth < 0 else self.depth

        # phi = (w H) / (2 k) * cosh(kz) / sinh(kd) * sin(kx - wt)
        #     = (w H) / (2 k) * cosh_ratio(kz, kd) / tanh(kd) * sin(kx - wt)
        # tanh(kd) -> 1 for deep water and never overflows.
        return (w * H) / (2 * k) * cosh_ratio(k * z, k * d) / tanh(k * d) * sin(k * x - w * t)


def compute_length_from_period(depth: float, period: float, g: float = 9.81) -> float:
    """
    Compute the wave length from the wave height, depth and period
    using the linear dispersion relation.
    """
    # Infinite depth approximation
    length = g * period**2 / (2 * pi)

    if depth < 0:
        # No need to compute the length via iterations for infinite depth
        return length

    # Find the length by Newton iterations
    length_1 = length * 0.95
    length_2 = length * 1.05
    T_1 = (length_1 * 2 * pi / (g * tanh(2 * pi * depth / length_1))) ** 0.5
    T_2 = (length_2 * 2 * pi / (g * tanh(2 * pi * depth / length_2))) ** 0.5
    length_N = 0.0
    iter = 0
    while abs(length_N - length) > 1e-4:
        # Store the previous length
        length = length_N

        # New guess for the wave length by interpolation
        f = (period - T_1) / (T_2 - T_1)
        length_N = length_1 + (length_2 - length_1) * f

        # Resulting period for the new length from the dispersion relation
        T_N = (length_N * 2 * pi / (g * tanh(2 * pi * depth / length_N))) ** 0.5

        # Update the two points used for the interpolation in the next iteration
        if T_N < period:
            length_1 = length_N
            T_1 = T_N
        else:
            length_2 = length_N
            T_2 = T_N

        iter += 1
        if iter > 100:
            raise NonConvergenceError(
                "Failed to converge when computing wave length from period for Airy waves"
            )

    return length_N
