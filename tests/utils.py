import sys
import pytest

try:
    from spectral_wave_data import SpectralWaveData  # NOQA

    skip_swd = False
except ImportError:
    skip_swd = True

skip_on_windows = pytest.mark.skipif(
    sys.platform.startswith("win"), reason="Skipping CMake/C++compilation on windows"
)

skip_swd_uninstalled = pytest.mark.skipif(
    skip_swd, reason="Could not import spectral_wave_data package"
)
