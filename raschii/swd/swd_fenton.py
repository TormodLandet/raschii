import numpy as np
from numpy import array, empty
from numpy.fft import irfft

from .base import SwdWriter


class SwdWriterFenton(SwdWriter):
    """SWD writer for Fenton stream-function waves."""

    def _effective_depth(self) -> float:
        wave = self.wave
        return 25.0 * wave.length if wave.depth < 0 else wave.depth

    def _elevation_coefficients(self, depth: float) -> np.ndarray:
        wave = self.wave
        nc = len(wave.eta)

        # Verify equally-spaced collocation points
        dx0 = abs(wave.x[1] - wave.x[0])
        assert all(
            abs(abs(wave.x[i + 2] - wave.x[i + 1]) - abs(wave.x[i + 1] - wave.x[i])) < 1e-4 * dx0
            for i in range(nc - 2)
        )

        # Reconstruct elevation Fourier cosine coefficients via irfft of the
        # half-wave collocation values mirrored to a full wavelength.
        nc2 = 2 * nc - 1
        etas = array([wave.eta[i] if i < nc else wave.eta[nc2 - i - 1] for i in range(nc2)])
        res = irfft(etas)

        ecs = empty(nc)
        ecs[0] = res[0].real - depth  # shift origin to calm free surface
        for j in range(1, nc - 1):
            ecs[j] = 2.0 * res[2 * j].real
        ecs[-1] = res[nc2 - 1].real
        return ecs

    def _potential_coefficients_z0(self, depth: float, nc: int) -> np.ndarray:
        # For Fenton, vcs[j] = i * B_j because:
        #   B_j * sin(j k (x - ct)) = Re{ i B_j exp(i j omega t) exp(-i j k x) }
        wave = self.wave
        B = wave.data["B"]
        vcs = empty(nc, complex)
        vcs[0] = 0.0
        for j in range(1, nc):
            vcs[j] = 1.0j * B[j]
        return vcs

    def _input_data(self, depth: float) -> dict:
        wave = self.wave
        return {
            "model": "Fenton",
            "T": wave.T,
            "height": wave.height,
            "depth": wave.depth,
            "depth_actual": depth,
            "N": wave.order,
            "air": wave.air.__class__.__name__,
            "g": wave.g,
            "c": wave.c,
            "relax": wave.relax,
        }

    def _order_zpos(self) -> int:
        return -1  # fully nonlinear
