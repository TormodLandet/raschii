import numpy as np
import pytest

import raschii


def make_wave(depth=10.0, **kw):
    return raschii.FentonWave(height=3.0, depth=depth, length=40.0, N=10, **kw)


def points(wave):
    x = np.linspace(0, wave.length, 7)
    z = wave.depth - 0.5 - np.linspace(0, 1, 7) * 3
    return x, z


def test_local_acceleration_matches_finite_difference():
    wave = make_wave()
    x, z = points(wave)
    dt = 1e-5
    fd = (wave.velocity(x, z, dt, all_points_wet=True) - wave.velocity(x, z, -dt, all_points_wet=True)) / (2 * dt)
    acc = wave.acceleration(x, z, 0.0, all_points_wet=True, kind="local")
    assert np.allclose(acc, fd, atol=1e-6)


def test_material_acceleration_matches_finite_difference_along_path():
    wave = make_wave()
    x, z = points(wave)
    dt = 1e-5
    vel = wave.velocity(x, z, 0.0, all_points_wet=True)
    # Follow the particles a small step forward and backward in time
    vp = wave.velocity(x + vel[:, 0] * dt, z + vel[:, 1] * dt, dt, all_points_wet=True)
    vm = wave.velocity(x - vel[:, 0] * dt, z - vel[:, 1] * dt, -dt, all_points_wet=True)
    fd = (vp - vm) / (2 * dt)
    acc = wave.acceleration(x, z, 0.0, all_points_wet=True)
    assert np.allclose(acc, fd, atol=1e-5)


def test_material_is_default_and_differs_from_local():
    wave = make_wave()
    x, z = points(wave)
    default = wave.acceleration(x, z, all_points_wet=True)
    material = wave.acceleration(x, z, all_points_wet=True, kind="material")
    local = wave.acceleration(x, z, all_points_wet=True, kind="local")
    assert np.array_equal(default, material)
    assert not np.allclose(material, local)


def test_invalid_kind():
    with pytest.raises(ValueError):
        make_wave().acceleration(0.0, 5.0, kind="bad")


@pytest.mark.parametrize("kind", ["material", "local"])
def test_deep_water_matches_finite_depth(kind):
    wave = raschii.FentonWave(height=3.0, depth=-1, length=100.0, N=10)
    deep = raschii.FentonWave(height=3.0, depth=5000.0, length=100.0, N=10)
    z_surf = 25 * wave.length  # infinite depth is modelled as 25 wave lengths
    x = np.array([0.0, 10.0, 30.0])
    z = np.array([z_surf - 1.0, z_surf - 5.0, z_surf - 20.0])
    acc = wave.acceleration(x, z, kind=kind, all_points_wet=True)
    ref = deep.acceleration(x, deep.depth - (z_surf - z), kind=kind, all_points_wet=True)
    assert np.all(np.isfinite(acc))
    assert np.allclose(acc, ref, atol=1e-6)


def test_above_surface_is_zero_without_air():
    wave = make_wave()
    acc = wave.acceleration(0.0, wave.depth + 10.0)
    assert np.allclose(acc, 0.0)
