import numpy as np
from .common import AIR_BLENDING_HEIGHT_FACTOR, np2py


class ConstantAirPhase:
    def __init__(self, height, blending_height=None):
        """
        Constant horizontal velocity equal to the phase speed
        """
        self.height = height
        self.blending_height = blending_height

    def set_wave(self, wave):
        """
        Connect this air phase with the wave in the water phase
        """
        self.c = wave.c
        self.depth_water = wave.depth

        if self.blending_height is None:
            self.blending_height = AIR_BLENDING_HEIGHT_FACTOR * wave.height

        from .cpp import ConstantAirCppGenerator
        self.cpp = ConstantAirCppGenerator(self)

    def stream_function(self, x, z, t=0, frame="b"):
        """
        Compute the stream function at time t for position(s) x
        """
        if isinstance(x, (float, int)):
            x, z = [x], [z]
        z = np.asarray(z, dtype=float)

        if frame == "e":
            return self.c * z
        elif frame == "c":
            return 0.0 * z

    def velocity(self, x, z, t=0):
        """
        Compute the air phase particle velocity at time t for position(s) (x, z)
        where z is 0 at the bottom and equal to depth_water at the free surface
        and equal to depth_water + depth air at the top free slip lid above the
        air phase
        """
        if isinstance(x, (float, int)):
            x, z = [x], [z]
        x = np.asarray(x, dtype=float)
        z = np.asarray(z, dtype=float)

        return np.zeros((x.size, 2), float)

    def __repr__(self):
        return f"ConstantAirPhase(height={self.height}, blending_height={self.blending_height})"
