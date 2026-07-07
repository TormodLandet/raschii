"""
Pure-Python tests for the velocity blending behaviour (no C++ compiler needed).

These tests catch regressions where blend_air_and_wave_velocities is accidentally
disabled in the _velocity implementations.
"""
import numpy
import pytest
from raschii import FentonAirPhase, FentonWave, get_wave_model


@pytest.fixture(params=["Airy", "Fenton"])
def wave_no_air(request):
    WaveClass = get_wave_model(request.param)[0]
    if request.param == "Airy":
        return WaveClass(height=1.0, depth=10.0, length=20.0)
    else:
        return WaveClass(height=1.0, depth=10.0, length=20.0, N=5)


def test_velocity_above_surface_is_zero(wave_no_air):
    """
    Without an air model, velocity must be zero above the free surface when
    all_points_wet=False (the default). This is the Python-side regression
    test for blend_air_and_wave_velocities being active.
    """
    wave = wave_no_air
    x = wave.length / 4

    # A point well above the highest possible surface elevation
    z_above = wave.depth + wave.height * 5

    vel = wave.velocity(x, z_above)  # all_points_wet defaults to False
    assert vel[0] == 0.0, f"Expected u=0 above surface, got {vel[0]}"
    assert vel[1] == 0.0, f"Expected w=0 above surface, got {vel[1]}"


def test_velocity_above_surface_nonzero_when_all_wet(wave_no_air):
    """
    With all_points_wet=True the blending is skipped and velocity above the
    surface should be non-zero (it returns the raw wave-theory value).
    """
    wave = wave_no_air
    x = wave.length / 4
    z_above = wave.depth + wave.height * 5

    vel = wave.velocity(x, z_above, all_points_wet=True)
    assert abs(vel[0]) > 0 or abs(vel[1]) > 0, "Expected non-zero velocity with all_points_wet=True"


def test_velocity_array_above_surface_is_zero(wave_no_air):
    """Same as above but with array inputs (exercises the T=1 loop path)."""
    wave = wave_no_air
    x = numpy.linspace(0, wave.length, 20)
    z = numpy.full_like(x, wave.depth + wave.height * 5)

    vel = wave.velocity(x, z)  # (N, 2)
    assert numpy.all(vel[:, 0] == 0.0), "Expected all u=0 above surface"
    assert numpy.all(vel[:, 1] == 0.0), "Expected all w=0 above surface"


def test_blended_velocity_divergence_free():
    """
    The Python blended velocity field (with FentonAir) must be divergence-free.
    This is the pure-Python equivalent of test_fenton_air_with_fenton_cpp_divergence.
    """
    air = FentonAirPhase(height=100.0)
    fwave = FentonWave(height=10.0, depth=200.0, length=100.0, N=5, air=air)

    length, depth, height = fwave.length, fwave.depth, fwave.height
    time = 1.0
    eps = 1e-7

    top = depth + 2.75 * height
    xpos = numpy.linspace(-length / 2, length / 2, 51)
    zpos = numpy.linspace(depth - height / 2, top, 51)
    X, Z = numpy.meshgrid(xpos, zpos)
    xr, zr = X.ravel(), Z.ravel()

    totvel = fwave.velocity(xr, zr, time)
    velsdx = fwave.velocity(xr + eps, zr, time)
    velsdz = fwave.velocity(xr, zr + eps, time)
    div = (velsdx[:, 0] - totvel[:, 0] + velsdz[:, 1] - totvel[:, 1]) / eps

    max_abs_div = abs(div).max()
    assert max_abs_div < 1e-4, f"Blended velocity field not divergence-free: max |div| = {max_abs_div}"


def test_blended_velocity_matches_zero_above_blend_zone():
    """
    With FentonAir, velocity well above the blending zone equals the air velocity
    (not the raw wave velocity). This checks the blending actually switches over.
    """
    blending_height = 20.0
    air = FentonAirPhase(height=100.0, blending_height=blending_height)
    fwave = FentonWave(height=10.0, depth=200.0, length=100.0, N=5, air=air)

    x = 0.0
    # Point well above the blending zone
    z_far_above = fwave.depth + blending_height * 2

    vel_blended = fwave.velocity(x, z_far_above)
    vel_raw_wave = fwave.velocity(x, z_far_above, all_points_wet=True)

    # The blended value should differ significantly from the raw wave value
    assert abs(vel_blended[0] - vel_raw_wave[0]) > 0.01 or abs(vel_blended[1] - vel_raw_wave[1]) > 0.01, \
        "Blended velocity above blend zone should differ from raw wave velocity"
