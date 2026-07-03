from .swd_stokes import SwdWriterStokes


class SwdWriterAiry:
    """
    SWD writer for Airy (linear) waves.

    Implemented via a Stokes N=1 wave, which is analytically identical to
    Airy theory.  Supports amp=1, 2, and 3.
    """

    def __init__(self, wave):
        from raschii.stokes import StokesWave

        # Build a Stokes N=1 wave from the Airy parameters.
        # The air model does not affect SWD output.
        self._stokes_wave = StokesWave(
            height=wave.height,
            depth=wave.depth,
            length=wave.length,
            N=1,
            g=wave.g,
        )
        self._writer = SwdWriterStokes(self._stokes_wave)

    def write(self, path, dt, tmax=None, nperiods=None, amp: int = 1):
        self._writer.write(path, dt, tmax=tmax, nperiods=nperiods, amp=amp)
