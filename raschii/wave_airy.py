from numpy import (
    cos,
    cosh,
    newaxis,
    pi,
    sin,
    sinh,
    stack,
    tanh,
)
from numpy.typing import NDArray

from .base_classes import AirPhaseModel, WaveModel
from .common import (
    NonConvergenceError,
    RaschiiError,
    blend_air_and_wave_velocities,
    cosh_ratio,
)
from .cpp import AiryCppGenerator


class AiryWave(WaveModel):
    required_input = {"height", "depth", "length"}
    optional_input = {"air": None, "g": 9.81}

    def __init__(
        self,
        height: float,
        depth: float,
        length: float | None = None,
        *,
        period: float | None = None,
        air: AirPhaseModel | None = None,
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
                raise RaschiiError("Either length or period must be given, both are None!")
            length = compute_length_from_period(depth=depth, period=period, g=g)

        #: The wave height
        self.height: float = height

        #: The water depth
        self.depth: float = depth

        #: The wave length
        self.length: float = length

        #: The optional air-phase model
        self.air: AirPhaseModel | None = air

        #: The acceleration of gravity
        self.g: float = g

        self.warnings: str = ""  #: Warnings raised when generating this wave

        #: Wave number (2 * pi / wavelength) in [1/m]
        self.k: float = 2 * pi / length

        # Wave angular frequency from the linear dispersion relation
        if self.depth < 0:
            # Infinite depth
            omega: float = (self.k * g) ** 0.5
        else:
            # Finite depth
            omega = (self.k * g * tanh(self.k * depth)) ** 0.5

        #: Wave angular frequency in [rad/s]
        self.omega: float = omega

        #: Wave celerity (phase speed) in [m/s]
        self.c: float = self.omega / self.k

        #: Wave period in [s]
        self.period: float = self.length / self.c

        # For evaluating velocities close to the free surface
        self.eta_eps: float = self.height / 1e5

        # Provide velocities also in the air phase
        if self.air is not None:
            self.air.set_wave(self)

        self.cpp = AiryCppGenerator(self)

    def _surface_elevation(self, x: NDArray, t: NDArray, include_depth: bool) -> NDArray:
        if include_depth:
            if self.depth < 0:
                raise RaschiiError("Cannot include depth in elevation for infinite depth")
            offset = self.depth
        else:
            offset = 0.0
        return offset + self.height / 2 * cos(self.k * x[newaxis, :] - self.omega * t[:, newaxis])

    def _velocity(self, x: NDArray, z: NDArray, t: NDArray, all_points_wet: bool) -> NDArray:
        if self.depth < 0:
            raise RaschiiError("Cannot currently compute velocity for infinite depth waves")
        x_1d, z_1d, t_1d = x, z, t  # save (N,), (N,), (T,) before reshaping
        H = self.height
        k = self.k
        d = self.depth
        w = self.omega
        phase = k * x[newaxis, :] - w * t[:, newaxis]  # (T, N)
        vel_x = w * H / 2 * cosh(k * z[newaxis, :]) / sinh(k * d) * cos(phase)
        vel_z = w * H / 2 * sinh(k * z[newaxis, :]) / sinh(k * d) * sin(phase)
        vel = stack([vel_x, vel_z], axis=-1)  # (T, N, 2)
        if not all_points_wet:
            for i, ti in enumerate(t_1d):
                blend_air_and_wave_velocities(
                    x_1d, z_1d, float(ti), self, self.air, vel[i], self.eta_eps
                )
        return vel

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

    def _velocity_potential(self, x: NDArray, z: NDArray, t: NDArray) -> NDArray:
        H = self.height
        k = self.k
        w = self.omega
        d = 25.0 * self.length if self.depth < 0 else self.depth
        return (
            (w * H)
            / (2 * k)
            * cosh_ratio(k * z[newaxis, :], k * d)
            / tanh(k * d)
            * sin(k * x[newaxis, :] - w * t[:, newaxis])
        )


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
