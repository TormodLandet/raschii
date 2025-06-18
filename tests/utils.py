import sys
import pytest
import os

try:
    from spectral_wave_data import SpectralWaveData  # NOQA

    skip_swd = False
except ImportError:
    skip_swd = True

if "SKIP_RASCHII_COMPILATION_TESTS" in os.environ:
    skipif_no_compile = pytest.mark.skip(
        reason="Skipping CMake/C++compilation (SKIP_RASCHII_COMPILATION_TESTS)",
    )
else:
    skipif_no_compile = pytest.mark.skipif(
        sys.platform.startswith("win"),
        reason="Skipping CMake/C++compilation on windows",
    )


skip_swd_uninstalled = pytest.mark.skipif(
    skip_swd, reason="Could not import spectral_wave_data package"
)
