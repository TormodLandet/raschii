import sys
import pytest

skip_on_windows = pytest.mark.skipif(
    sys.platform.startswith("win"), reason="Skipping CMake/C++compilation on windows"
)
