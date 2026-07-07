"""
Integration tests for the raschii command line interface.

These tests run the CLI tools in a subprocess using the same Python
interpreter that is running the test suite, so no path juggling is needed.
Plot tests are skipped when matplotlib is not installed; install
``raschii[plot]`` to run them.
"""
import os
import subprocess
import sys

import pytest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _run(*args):
    """Run a command as a subprocess and return CompletedProcess."""
    return subprocess.run(
        [sys.executable, *args],
        capture_output=True,
        text=True,
    )


try:
    import matplotlib  # noqa: F401

    has_matplotlib = True
except ImportError:
    has_matplotlib = False

skip_no_matplotlib = pytest.mark.skipif(
    not has_matplotlib,
    reason="matplotlib not installed — pip install raschii[plot] to run plot tests",
)


# ---------------------------------------------------------------------------
# raschii.cmd.swd  (no matplotlib dependency)
# ---------------------------------------------------------------------------


def test_swd_cmd_by_length(tmp_path):
    """SWD writer: positional wave_length argument."""
    swd = str(tmp_path / "wave.swd")
    result = _run("-m", "raschii.cmd.swd", "-N", "5", swd, "Fenton", "0.2", "1.5", "2.0")
    assert result.returncode == 0, f"Command failed:\n{result.stderr}"
    assert os.path.isfile(swd), "SWD file was not created"


def test_swd_cmd_by_period(tmp_path):
    """SWD writer: --period flag instead of positional wave_length."""
    swd = str(tmp_path / "wave.swd")
    result = _run(
        "-m", "raschii.cmd.swd", "-N", "5", "--period", "1.2",
        swd, "Fenton", "0.2", "1.5",
    )
    assert result.returncode == 0, f"Command failed:\n{result.stderr}"
    assert os.path.isfile(swd), "SWD file was not created"


def test_swd_cmd_stokes(tmp_path):
    """SWD writer: Stokes wave model."""
    swd = str(tmp_path / "stokes.swd")
    result = _run("-m", "raschii.cmd.swd", "-N", "3", swd, "Stokes", "0.2", "1.5", "2.0")
    assert result.returncode == 0, f"Command failed:\n{result.stderr}"
    assert os.path.isfile(swd)


def test_swd_cmd_missing_length_and_period(tmp_path):
    """SWD writer: exits non-zero when neither wave_length nor --period is given."""
    swd = str(tmp_path / "wave.swd")
    result = _run("-m", "raschii.cmd.swd", swd, "Fenton", "0.2", "1.5")
    assert result.returncode != 0


# ---------------------------------------------------------------------------
# raschii.cmd.plot  (require matplotlib)
# ---------------------------------------------------------------------------


@skip_no_matplotlib
def test_plot_cmd_saves_png_by_length(tmp_path):
    """Plot command saves PNG files when --output-prefix is given (wave_length)."""
    prefix = str(tmp_path / "plot")
    result = _run(
        "-m", "raschii.cmd.plot", "-N", "5",
        "--output-prefix", prefix,
        "Fenton", "0.2", "1.5", "2.0",
    )
    assert result.returncode == 0, f"Command failed:\n{result.stderr}"
    assert os.path.isfile(prefix + "_elevation.png"), "Elevation PNG not created"
    assert os.path.isfile(prefix + "_velocities.png"), "Velocities PNG not created"


@skip_no_matplotlib
def test_plot_cmd_saves_png_by_period(tmp_path):
    """Plot command saves PNG files when --period is used instead of wave_length."""
    prefix = str(tmp_path / "plot_period")
    result = _run(
        "-m", "raschii.cmd.plot", "-N", "5",
        "--period", "1.2", "--output-prefix", prefix,
        "Fenton", "0.2", "1.5",
    )
    assert result.returncode == 0, f"Command failed:\n{result.stderr}"
    assert os.path.isfile(prefix + "_elevation.png")
    assert os.path.isfile(prefix + "_velocities.png")


@skip_no_matplotlib
def test_plot_cmd_multiple_models(tmp_path):
    """Plot command handles a comma-separated list of wave models."""
    prefix = str(tmp_path / "multi")
    result = _run(
        "-m", "raschii.cmd.plot", "-N", "5",
        "--output-prefix", prefix,
        "Fenton,Stokes", "0.2", "1.5", "2.0",
    )
    assert result.returncode == 0, f"Command failed:\n{result.stderr}"
    assert os.path.isfile(prefix + "_elevation.png")
    assert os.path.isfile(prefix + "_velocities.png")
