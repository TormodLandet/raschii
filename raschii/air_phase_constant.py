import numpy as np

from .common import AIR_BLENDING_HEIGHT_FACTOR, Frame
from .base_classes import AirPhaseModel


class ConstantAirPhase(AirPhaseModel):
    def __init__(self, height, blending_height=None):
        """
        Air phase model with zero velocity in the earth frame (still air).

        The wave-frame stream function contains a backward drift at the wave
        phase speed, which allows divergence-free blending with the water-phase
        velocity field across the free surface.
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

    def stream_function(self, x, z, t=0, frame=Frame.EARTH):
        """
        Compute the stream function at time t for position(s) x.

        * frame: :class:`~raschii.Frame` – ``Frame.EARTH`` (default) returns
          zero (still air in the earth frame); ``Frame.WAVE`` returns ``-c*z``
          (air appears to move backward at the wave phase speed in the
          co-moving frame).
        """
        if isinstance(x, (float, int)):
            x, z = [x], [z]
        z = np.asarray(z, dtype=float)

        if frame == Frame.EARTH:
            return np.zeros_like(z)
        elif frame == Frame.WAVE:
            return -self.c * z
        else:
            raise ValueError(f"Unknown frame {frame!r}; use Frame.EARTH or Frame.WAVE")

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
