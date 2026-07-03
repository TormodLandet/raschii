import pytest

from raschii import WAVE_MODELS, WaveModel, get_wave_model


def make_example_wave(model: str) -> WaveModel:
    wave_model, _ = get_wave_model(model)

    height = 10.0
    depth = 200.0
    length = 100.0
    N = 5 if model in ["Fenton", "Stokes"] else 1
    return wave_model(height, depth, length, N)


@pytest.mark.parametrize("model", WAVE_MODELS.keys())
def test_elevation_shape(model: str):
    wave = make_example_wave(model)

    x = wave.length / 10
    elev = wave.surface_elevation(x)
    assert elev.shape == (1,)

    x2 = [x, x + 1, x + 2, x + 3]
    elev2 = wave.surface_elevation(x2)
    assert elev2.shape == (4,)

    elev3 = wave.surface_elevation(x, t=[0, 1, 2])
    assert elev3.shape == (3,)

    elev4 = wave.surface_elevation(x2, t=[0, 1, 2])
    assert elev4.shape == (3, 4)


@pytest.mark.parametrize("model", WAVE_MODELS.keys())
def test_velocity_shape(model: str):
    wave = make_example_wave(model)

    x = wave.length / 10
    z = wave.   depth - 1
    vel = wave.velocity(x, z, all_points_wet=True)
    assert vel.shape == (2,)

    x2 = [x, x + 1, x + 2, x + 3]
    z2 = [z, z, z, z]
    vel2 = wave.velocity(x2, z2, all_points_wet=True)
    assert vel2.shape == (4, 2)

    vel3 = wave.velocity(x, z, t=[0, 1, 2], all_points_wet=True)
    assert vel3.shape == (3, 2)

    vel4 = wave.velocity(x2, z2, t=[0, 1, 2], all_points_wet=True)
    assert vel4.shape == (3, 4, 2)


@pytest.mark.parametrize("model", WAVE_MODELS.keys())
def test_velocity_potential_shape(model: str):
    wave = make_example_wave(model)

    x = wave.length / 10
    z = wave.depth - 1
    vp = wave.velocity_potential(x, z)
    assert vp.shape == (1,)

    x2 = [x, x + 1, x + 2, x + 3]
    z2 = [z, z, z, z]
    vp2 = wave.velocity_potential(x2, z2)
    assert vp2.shape == (4,)

    vp3 = wave.velocity_potential(x, z, t=[0, 1, 2])
    assert vp3.shape == (3,)

    vp4 = wave.velocity_potential(x2, z2, t=[0, 1, 2])
    assert vp4.shape == (3, 4)
