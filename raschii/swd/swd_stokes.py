from math import pi, sqrt

import numpy as np

from .base import SwdWriter


class SwdWriterStokes(SwdWriter):
    """SWD writer for Stokes waves."""

    def _effective_depth(self) -> float:
        wave = self.wave
        if wave.depth < 0:
            return 50.0 * pi / wave.k
        return min(wave.depth, 50.0 * pi / wave.k)

    def _elevation_coefficients(self, depth: float) -> np.ndarray:
        wave = self.wave
        k = wave.k
        eps = k * wave.height / 2
        D = wave.data
        nc = wave.order + 1

        ecs = np.zeros(nc, complex)
        ecs[0] = 0.0  # zero at calm water line
        ecs[1] = eps + eps**3 * D["B31"] - eps**5 * (D["B53"] + D["B55"])
        if wave.order > 1:
            ecs[2] = eps**2 * D["B22"] + eps**4 * D["B42"]
        if wave.order > 2:
            ecs[3] = -(eps**3) * D["B31"] + eps**5 * D["B53"]
        if wave.order > 3:
            ecs[4] = eps**4 * D["B44"]
        if wave.order > 4:
            ecs[5] = eps**5 * D["B55"]
        ecs /= k
        return ecs

    def _potential_coefficients_z0(self, depth: float, nc: int) -> np.ndarray:
        wave = self.wave
        kd = wave.k * depth
        eps = wave.k * wave.height / 2
        D = wave.data

        vcs = np.zeros(nc, complex)
        vcs[1] = (eps * D["A11"] + eps**3 * D["A31"] + eps**5 * D["A51"]) * np.cosh(kd)
        if wave.order > 1:
            m = eps**2 * D["A22"] + eps**4 * D["A42"]
            if m != 0.0:
                vcs[2] = m * np.cosh(2 * kd)
        if wave.order > 2:
            m = eps**3 * D["A33"] + eps**5 * D["A53"]
            if m != 0.0:
                vcs[3] = m * np.cosh(3 * kd)
        if wave.order > 3:
            m = eps**4 * D["A44"]
            if m != 0.0:
                vcs[4] = m * np.cosh(4 * kd)
        if wave.order > 4:
            m = eps**5 * D["A55"]
            if m != 0.0:
                vcs[5] = m * np.cosh(5 * kd)
        vcs *= 1.0j * D["C0"] * sqrt(wave.g / wave.k**3)
        return vcs

    def _input_data(self, depth: float) -> dict:
        wave = self.wave
        return {
            "model": "Stokes",
            "T": wave.T,
            "length": wave.length,
            "height": wave.height,
            "depth": wave.depth,
            "depth_actual": depth,
            "N": wave.order,
            "air": wave.air.__class__.__name__,
            "g": wave.g,
            "c": wave.c,
        }

    def _order_zpos(self) -> int:
        return self.wave.order
