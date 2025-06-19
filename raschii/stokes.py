from math import pi, sinh, tanh, exp, sqrt, pow
import numpy as np
from .common import blend_air_and_wave_velocities, RaschiiError, NonConvergenceError
from .swd_tools import SwdShape1and2
from .base_classes import WaveModel


class StokesWave(WaveModel):
    required_input = {"height", "depth", "length", "N"}
    optional_input = {"air": None, "g": 9.81}

    def __init__(
        self,
        height: float,
        depth: float,
        length: float | None = None,
        N: int = 5,
        period: float | None = None,
        air=None,
        g: float = 9.81,
    ):
        """
        Implement Stokes waves based on the paper by J. D. Fenton (1985),
        "A Fifth-Order Stokes Theory for Steady Waves".

        * height: wave height above still water level
        * depth: still water distance from the flat sea bottom to the surface
          in meters, but you can give -1.0 for infinite depth
        * length: the periodic length of the wave (optional, if not given then period is used)
        * N: the number of coefficients in the truncated Fourier series
        * period: the wave period (optional, if not given then length is used)
        """
        if length is None:
            if period is None:
                raise RaschiiError("Either length or period must be given, both are None!")
            length = compute_length_from_period(height=height, depth=depth, period=period, N=N, g=g)

        self.height: float = height  #: The wave height
        self.depth: float = depth  #: The water depth
        self.length: float = length  #: The wave length
        self.order: int = N  #: The approximation order
        self.air = air  #: The optional air-phase model
        self.g: float = g  #: The acceleration of gravity
        self.warnings: str = ""  #: Warnings raised when generating this wave

        if N < 1:
            self.warnings = "Stokes order must be at least 1, using order 1"
            self.order = 1
        elif N > 5:
            self.warnings = "Stokes order is maximum 5, using order 5"
            self.order = 4

        # Find the coeffients through explicit formulas
        self.k = 2 * pi / length  # The wave number
        data = stokes_coefficients(self.k * depth, self.order)
        self.set_data(data)

        # For evaluating velocities close to the free surface
        self.eta_eps = self.height / 1e5

        # Provide velocities also in the air phase
        if self.air is not None:
            self.air.set_wave(self)

    def set_data(self, data):
        """
        Update the coefficients defining this Stokes wave
        """
        self.data = data
        k = 2 * pi / self.length
        eps = k * self.height / 2
        d = self.depth

        # Apply consistent water depth limitation with 'stokes_coefficients': kd <= 50 * pi
        if d < 0:
            d = 50 * pi / self.k
        d = min(d, 50 * pi / self.k)

        c = (data["C0"] + pow(eps, 2) * data["C2"] + pow(eps, 4) * data["C4"]) * sqrt(self.g / k)
        Q = (c * d * sqrt(k**3 / self.g) + data["D2"] * eps**2 + data["D4"] * eps**4) * sqrt(
            self.g / k**3
        )

        self.c = c  # Phase speed
        self.cs = c - Q  # Mean Stokes drift speed (TODO: verify this!)
        self.T = self.length / self.c  # Wave period
        self.omega = self.c * self.k  # Wave frequency

    def surface_elevation(self, x: float | list[float], t: float = 0.0, include_depth: bool = True):
        """
        Compute the surface elavation at time t for position(s) x
        """
        if isinstance(x, (float, int)):
            x = np.array([x], float)
        x = np.asarray(x)
        x2 = x - self.c * t

        k = self.k
        eps = k * self.height / 2
        D = self.data
        cos = np.cos
        eta = (
            eps * cos(k * x2)
            + eps**2 * D["B22"] * cos(2 * k * x2)
            + eps**3 * D["B31"] * (cos(k * x2) - cos(3 * k * x2))
            + eps**4 * (D["B42"] * cos(2 * k * x2) + D["B44"] * cos(4 * k * x2))
            + eps**5
            * (
                -(D["B53"] + D["B55"]) * cos(k * x2)
                + D["B53"] * cos(3 * k * x2)
                + D["B55"] * cos(5 * k * x2)
            )
        ) / k

        if include_depth:
            if self.depth < 0:
                raise RaschiiError("Cannot include depth in elevation for infinite depth")
            offset = self.depth
        else:
            offset = 0.0

        return offset + eta

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
            raise RaschiiError("Cannot currently compute velocity for infinite depth waves")
        if isinstance(x, (float, int)):
            x, z = [x], [z]
        x = np.asarray(x, dtype=float)
        z = np.asarray(z, dtype=float)
        x2 = x - self.c * t

        def my_cosh_cos(i, j):
            n = "A%d%d" % (i, j)
            if self.data[n] == 0.0:
                return 0.0
            else:
                return (
                    pow(eps, i)
                    * self.data[n]
                    * j
                    * self.k
                    * np.cosh(j * self.k * z)
                    * np.cos(j * self.k * x2)
                )

        def my_sinh_sin(i, j):
            n = "A%d%d" % (i, j)
            if self.data[n] == 0.0:
                return 0.0
            else:
                return (
                    pow(eps, i)
                    * self.data[n]
                    * j
                    * self.k
                    * np.sinh(j * self.k * z)
                    * np.sin(j * self.k * x2)
                )

        eps = self.k * self.height / 2
        vel = np.zeros((x.size, 2), float)
        vel[:, 0] = (
            my_cosh_cos(1, 1)
            + my_cosh_cos(2, 2)
            + my_cosh_cos(3, 1)
            + my_cosh_cos(3, 3)
            + my_cosh_cos(4, 2)
            + my_cosh_cos(4, 4)
            + my_cosh_cos(5, 1)
            + my_cosh_cos(5, 3)
            + my_cosh_cos(5, 5)
        )
        vel[:, 1] = (
            my_sinh_sin(1, 1)
            + my_sinh_sin(2, 2)
            + my_sinh_sin(3, 1)
            + my_sinh_sin(3, 3)
            + my_sinh_sin(4, 2)
            + my_sinh_sin(4, 4)
            + my_sinh_sin(5, 1)
            + my_sinh_sin(5, 3)
            + my_sinh_sin(5, 5)
        )
        vel *= self.data["C0"] * sqrt(self.g / self.k**3)

        if not all_points_wet:
            blend_air_and_wave_velocities(x, z, t, self, self.air, vel, self.eta_eps)

        return vel

    def write_swd(self, path, dt, tmax=None, nperiods=None):
        """
        Write a SWD-file of the wave field according to the file
        specification in the Github repository spectral-wave-data ....

        * path:     Full path of the new SWD file
        * dt:       The temporal sampling spacing in the SWD file
        * tmax:     The temporal sampling range in the SWD file is [0, tmax]
        * nperiods: Alternative specification: tmax = nperiods * wave_period
        """
        if tmax is None:
            assert nperiods is not None
            tmax = nperiods * self.T
        assert tmax > dt > 0.0

        # Apply consistent water depth limitation with 'stokes_coefficients': kd <= 50 * pi
        if self.depth < 0:
            depth = 50 * pi / self.k
        else:
            depth = min(self.depth, 50 * pi / self.k)

        kd = self.k * depth
        eps = self.k * self.height / 2
        D = self.data

        # The swd coordinate system is earth fixed with zswd=0 in the calm
        # surface. Hence:
        #   xswd = x + t * self.c    and    zswd = z - depth

        # zeta(x) = sum[ecs[j] * cos(j * self.k * x) for j in range(nc)]

        # Note that ecs[j] * cos(j * self.k * x) == Re{h[j, t] * exp(-I * j * self.k * xswd)}
        # where h[j, t] = ecs[j] * exp(-I * j * self.k * (-self.c * t)),   j>=0, I=sqrt(-1)
        # Hence dh[j, t]/dt = I * j * self.k * self.c * h[j, t]

        nc = self.order + 1  # Include Bias term
        ecs = np.empty(nc, complex)
        ecs[0] = 0.0  # zero at calm water line
        ecs[1] = eps + eps**3 * D["B31"] - eps**5 * (D["B53"] + D["B55"])
        if self.order > 1:
            ecs[2] = eps**2 * D["B22"] + eps**4 * D["B42"]
        if self.order > 2:
            ecs[3] = -(eps**3) * D["B31"] + eps**5 * D["B53"]
        if self.order > 3:
            ecs[4] = eps**4 * D["B44"]
        if self.order > 4:
            ecs[5] = eps**5 * D["B55"]
        ecs /= self.k

        # From particle velocities we construct the SWD velocity potential...

        # Note:  vcs[j] * cos(j * self.k * x) == Im{c[j, t] * exp(-I * j * self.k * xswd)}
        # where c[j, t] = I * vcs[j] * exp(-I * j * self.k * (-self.c * t)),     j>0
        # Hence dc[j, t]/dt = I * j * self.k * self.c * c[j, t]

        vcs = np.zeros(nc, complex)
        vcs[0] = 0.0
        vcs[1] = (eps * D["A11"] + eps**3 * D["A31"] + eps**5 * D["A51"]) * np.cosh(kd)

        if self.order > 1:
            multiplier = eps**2 * D["A22"] + eps**4 * D["A42"]
            if multiplier != 0.0:
                vcs[2] = multiplier * np.cosh(2 * kd)

        if self.order > 2:
            multiplier = eps**3 * D["A33"] + eps**5 * D["A53"]
            if multiplier != 0.0:
                vcs[3] = multiplier * np.cosh(3 * kd)

        if self.order > 3:
            multiplier = eps**4 * D["A44"]
            if multiplier != 0.0:
                vcs[4] = multiplier * np.cosh(4 * kd)

        if self.order > 4:
            multiplier = eps**5 * D["A55"]
            if multiplier != 0.0:
                vcs[5] = multiplier * np.cosh(5 * kd)

        vcs *= 1.0j * D["C0"] * sqrt(self.g / self.k**3)

        input_data = {
            "model": "Stokes",
            "T": self.T,
            "length": self.length,
            "height": self.height,
            "depth": self.depth,
            "depth_actual": depth,
            "N": self.order,
            "air": self.air.__class__.__name__,
            "g": self.g,
            "c": self.c,
        }

        swd = SwdShape1and2(
            self.T, self.length, self.depth, vcs, ecs, input_data, self.g, order_zpos=self.order
        )
        swd.write(path, dt, tmax=tmax)


def sech(x):
    "Hyperbolic secant"
    return 2 * exp(x) / (exp(2 * x) + 1)


def csch(x):
    "Hyperbolic cosecant"
    return 2 * exp(x) / (exp(2 * x) - 1)


def cotanh(x):
    "Hyperbolic cotangent"
    return (1 + exp(-2 * x)) / (1 - exp(-2 * x))


def stokes_coefficients(kd, N):
    """
    Define the Stokes expansion coefficients based on "A Fifth‐Order Stokes
    Theory for Steady Waves" (J. D. Fenton, 1985)

    The code uses pow instead of ** to be compatible with Dart
    """
    # Limit depth to 25 wave lengths to avoid overflow
    if kd > 50 * pi or kd < 0:
        kd = 50 * pi

    S = sech(2 * kd)
    Sh = sinh(kd)
    Th = tanh(kd)
    CTh = cotanh(kd)

    # Parameters are zero if not used in linear theory
    data = {
        "A11": csch(kd),
        "A22": 0.0,
        "A31": 0.0,
        "A33": 0.0,
        "A42": 0.0,
        "A44": 0.0,
        "A51": 0.0,
        "A53": 0.0,
        "A55": 0.0,
        "B22": 0.0,
        "B31": 0.0,
        "B42": 0.0,
        "B44": 0.0,
        "B53": 0.0,
        "B55": 0.0,
        "C0": sqrt(Th),
        "C2": 0.0,
        "C4": 0.0,
        "D2": 0.0,
        "D4": 0.0,
        "E2": 0.0,
        "E4": 0.0,
    }

    if N == 1:
        return data

    # Define additional constants needed for second order Stokes waves
    data["A22"] = 3 * pow(S, 2) / (2 * pow(1 - S, 2))
    data["B22"] = CTh * (1 + 2 * S) / (2 * (1 - S))
    data["C2"] = sqrt(Th) * (2 + 7 * pow(S, 2)) / (4 * pow(1 - S, 2))
    data["D2"] = -sqrt(CTh) / 2
    data["E2"] = Th * (2 + 2 * S + 5 * pow(S, 2)) / (4 * pow(1 - S, 2))

    if N == 2:
        return data

    # Define additional constants needed for third order Stokes waves
    data["A31"] = (-4 - 20 * S + 10 * pow(S, 2) - 13 * pow(S, 3)) / (8 * Sh * pow(1 - S, 3))
    data["A33"] = (-2 * pow(S, 2) + 11 * pow(S, 3)) / (8 * Sh * pow(1 - S, 3))
    data["B31"] = -3 * (1 + 3 * S + 3 * pow(S, 2) + 2 * pow(S, 3)) / (8 * pow(1 - S, 3))

    if N == 3:
        return data

    # Define additional constants needed for forth order Stokes waves
    data["A42"] = (12 * S - 14 * pow(S, 2) - 264 * pow(S, 3) - 45 * pow(S, 4) - 13 * pow(S, 5)) / (
        24 * pow(1 - S, 5)
    )
    data["A44"] = (10 * pow(S, 3) - 174 * pow(S, 4) + 291 * pow(S, 5) + 278 * pow(S, 6)) / (
        48 * (3 + 2 * S) * pow(1 - S, 5)
    )
    data["B42"] = (
        CTh
        * (6 - 26 * S - 182 * pow(S, 2) - 204 * pow(S, 3) - 25 * pow(S, 4) + 26 * pow(S, 5))
        / (6 * (3 + 2 * S) * pow(1 - S, 4))
    )
    data["B44"] = (
        CTh
        * (24 + 92 * S + 122 * pow(S, 2) + 66 * pow(S, 3) + 67 * pow(S, 4) + 34 * pow(S, 5))
        / (24 * (3 + 2 * S) * pow(1 - S, 4))
    )
    data["C4"] = (
        sqrt(Th)
        * (4 + 32 * S - 116 * pow(S, 2) - 400 * pow(S, 3) - 71 * pow(S, 4) + 146 * pow(S, 5))
        / (32 * pow(1 - S, 5))
    )
    data["D4"] = sqrt(CTh) * (2 + 4 * S + pow(S, 2) + 2 * pow(S, 3)) / (8 * pow(1 - S, 3))
    data["E4"] = (
        Th
        * (8 + 12 * S - 152 * pow(S, 2) - 308 * pow(S, 3) - 42 * pow(S, 4) + 77 * pow(S, 5))
        / (32 * pow(1 - S, 5))
    )

    if N == 4:
        return data

    # Define additional constants needed for fifth order Stokes waves
    data["A51"] = (
        -1184
        + 32 * S
        + 13232 * pow(S, 2)
        + 21712 * pow(S, 3)
        + 20940 * pow(S, 4)
        + 12554 * pow(S, 5)
        - 500 * pow(S, 6)
        - 3341 * pow(S, 7)
        - 670 * pow(S, 8)
    ) / (64 * Sh * (3 + 2 * S) * (4 + S) * pow(1 - S, 6))
    data["A53"] = (
        4 * S
        + 105 * pow(S, 2)
        + 198 * pow(S, 3)
        - 1376 * pow(S, 4)
        - 1302 * pow(S, 5)
        - 117 * pow(S, 6)
        + 58 * pow(S, 7)
    ) / (32 * Sh * (3 + 2 * S) * pow(1 - S, 6))
    data["A55"] = (
        -6 * pow(S, 3)
        + 272 * pow(S, 4)
        - 1552 * pow(S, 5)
        + 852 * pow(S, 6)
        + 2029 * pow(S, 7)
        + 430 * pow(S, 8)
    ) / (64 * Sh * (3 + 2 * S) * (4 + S) * pow(1 - S, 6))
    data["B53"] = (
        9
        * (
            132
            + 17 * S
            - 2216 * pow(S, 2)
            - 5897 * pow(S, 3)
            - 6292 * pow(S, 4)
            - 2687 * pow(S, 5)
            + 194 * pow(S, 6)
            + 467 * pow(S, 7)
            + 82 * pow(S, 8)
        )
        / (128 * (3 + 2 * S) * (4 + S) * pow(1 - S, 6))
    )
    data["B55"] = (
        5
        * (
            300
            + 1579 * S
            + 3176 * pow(S, 2)
            + 2949 * pow(S, 3)
            + 1188 * pow(S, 4)
            + 675 * pow(S, 5)
            + 1326 * pow(S, 6)
            + 827 * pow(S, 7)
            + 130 * pow(S, 8)
        )
        / (384 * (3 + 2 * S) * (4 + S) * pow(1 - S, 6))
    )

    return data


def compute_length_from_period(
    height: float,
    depth: float,
    period: float,
    N: int = 5,
    g: float = 9.81,
):
    """
    Compute the wave length from the wave period using the Fenton wave theory

    This would be much faster if we had an implementation of the Stokes wave
    theory dispersion relation, which should be not too hard to implement ...
    """
    from .airy import compute_length_from_period as airy_compute_length_from_period

    # Initial guess is based on the linear dispersion relation for deep water waves
    length = airy_compute_length_from_period(depth=depth, period=period, g=g)

    # Find the length by Newton iterations
    wave1 = StokesWave(height=height, depth=depth, length=length * 0.95, N=N, g=g)
    wave2 = StokesWave(height=height, depth=depth, length=length * 1.05, N=N, g=g)

    length_N = 0.0
    iter = 0
    while abs(length_N - length) > 1e-4:
        # Store the previous length
        length = length_N

        # New guess for the wave length by interpolation
        f = (period - wave1.T) / (wave2.T - wave1.T)
        length_N = wave1.length + (wave2.length - wave1.length) * f

        # Resulting wave period for the new length from the dispersion relation
        waveN = StokesWave(height=height, depth=depth, length=length_N, N=N, g=g)

        # Update the two points used for the interpolation in the next iteration
        if waveN.T < period:
            wave1 = waveN
        else:
            wave2 = waveN

        iter += 1
        if iter > 100:
            raise NonConvergenceError(
                "Failed to converge when computing wave length from period for Stokes waves"
            )

    return length_N
