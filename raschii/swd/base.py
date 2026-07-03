import numpy as np

from .swd_file import SwdShape1and2
from raschii.common import RasciiError


class SwdWriter:
    """
    Base class for writing SWD files from Raschii wave models.

    Subclasses implement the model-specific hooks; this base class provides
    the orchestration, amp validation, and the generic amp=2 surface-potential
    extraction that works for any wave model exposing ``velocity_potential`` and
    ``surface_elevation``.
    """

    def __init__(self, wave):
        self.wave = wave

    def write(self, path, dt, tmax=None, nperiods=None, amp: int = 1):
        """
        Write the SWD file.

        * path:     Full path of the SWD file to create
        * dt:       Temporal sampling spacing [s]
        * tmax:     Duration of the time series [s]
        * nperiods: Alternative to tmax: tmax = nperiods * wave_period
        * amp:      SWD amp flag (1, 2, or 3):

          - 1: store velocity-potential coefficients evaluated at z=0 (the calm
               free surface).  Full kinematics available everywhere via the
               vertical shape functions Z_j(z).  *Default.*
          - 2: store velocity-potential coefficients evaluated on the actual
               wavy free surface ψ(x) = φ(x, η(x)).  The Z_j dependency is
               removed; reconstruction of kinematics at arbitrary depth requires
               the H2-operator in the experimental SWD reader.
          - 3: store elevation only.  No potential data; only surface-elevation
               queries are possible with the resulting file.
        """
        if amp not in (1, 2, 3):
            raise RasciiError(f"SWD amp must be 1, 2, or 3, got {amp!r}")

        wave = self.wave
        if tmax is None:
            if nperiods is None:
                raise RasciiError("Either tmax or nperiods must be given")
            tmax = nperiods * wave.T
        if not (tmax > dt > 0.0):
            raise RasciiError(f"Must have tmax > dt > 0, got tmax={tmax!r}, dt={dt!r}")

        depth = self._effective_depth()
        ecs = self._elevation_coefficients(depth)
        nc = len(ecs)

        if amp == 1:
            vcs = self._potential_coefficients_z0(depth, nc)
        elif amp == 2:
            # The surface potential φ(x, η(x)) generates harmonics up to ~2N due to
            # cosh(jk·η)·sin(jk·x) cross-products.  Storing only N harmonics gives
            # ~2% error for nonlinear shallow-water waves.  We therefore expand the
            # spectral resolution to AMP2_NC_FACTOR × N for amp=2 files, zero-padding
            # the elevation coefficients beyond N so the elevation reconstruction is
            # unchanged while the potential is represented accurately.
            AMP2_NC_FACTOR = 4
            nc_amp2 = nc * AMP2_NC_FACTOR
            ecs_amp2 = np.zeros(nc_amp2, complex)
            ecs_amp2[:nc] = ecs
            vcs = self._potential_coefficients_surface(depth, nc_amp2)
            ecs = ecs_amp2
        else:  # amp == 3
            vcs = np.zeros(nc, complex)  # not written, but SwdShape1and2 requires same length

        swd = SwdShape1and2(
            wave.T,
            wave.length,
            wave.depth,
            vcs,
            ecs,
            self._input_data(depth),
            wave.g,
            order_zpos=self._order_zpos(),
        )
        swd.write(path, dt, tmax=tmax, amp=amp)

    # ------------------------------------------------------------------
    # Generic amp=2 extraction
    # ------------------------------------------------------------------

    def _potential_coefficients_surface(self, depth: float, nc: int) -> np.ndarray:
        """
        Compute SWD velocity-potential coefficients for amp=2 by sampling
        φ(x, η(x)) over one wavelength with oversample=8, then extracting the
        sine-series amplitudes via an inverse DFT.

        The stored c_j coefficients satisfy:
            c_j(t=0) = i * p_j
        where φ(x) = Σ p_j * sin(j k x)  (pure sine series, x measured from trough).
        Time evolution follows c_j(t) = c_j(0) * exp(i j ω t) as for amp=1.
        """
        wave = self.wave
        M = 8 * nc  # oversampling factor = 8

        x_samp = np.linspace(0, wave.length, M, endpoint=False)

        # Surface elevation in z-from-bottom coordinates.
        # surface_elevation(include_depth=True) raises for infinite depth, so
        # we fall back to the zero-mean elevation and add the effective depth.
        if wave.depth < 0:
            z_surf = wave.surface_elevation(x_samp, t=0, include_depth=False) + depth
        else:
            z_surf = wave.surface_elevation(x_samp, t=0, include_depth=True)

        phi_samp = wave.velocity_potential(x_samp, z_surf, t=0)

        # DFT: Im{rfft[phi][j]} = -sum_m phi_m * sin(2π j m / M)
        # so p_j = -2/M * Im{rfft[phi][j]}
        X = np.fft.rfft(phi_samp)
        vcs = np.zeros(nc, complex)
        for j in range(1, nc):
            p_j = -2.0 / M * X[j].imag
            vcs[j] = 1j * p_j
        return vcs

    # ------------------------------------------------------------------
    # Abstract hooks — implement in subclasses
    # ------------------------------------------------------------------

    def _effective_depth(self) -> float:
        """Return the depth used in the SWD file (finite positive number)."""
        raise NotImplementedError

    def _elevation_coefficients(self, depth: float) -> np.ndarray:
        """Return the real elevation Fourier amplitudes ecs[j]."""
        raise NotImplementedError

    def _potential_coefficients_z0(self, depth: float, nc: int) -> np.ndarray:
        """Return complex SWD potential coefficients for amp=1 (evaluated at z=0)."""
        raise NotImplementedError

    def _input_data(self, depth: float) -> dict:
        """Return the metadata dict stored in the SWD file header."""
        raise NotImplementedError

    def _order_zpos(self) -> int:
        """Return the SWD order_zpos flag."""
        raise NotImplementedError
