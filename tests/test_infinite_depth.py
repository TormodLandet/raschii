import pytest


@pytest.mark.parametrize("give_period", [True, False])
def test_infinite_depth_airy(give_period: bool):
    import raschii

    WaveModel, _AirModel = raschii.get_wave_model("Airy")

    if give_period:
        inp = dict(height=12, period=15)
    else:
        inp = dict(height=12, length=350)

    assert raschii.check_breaking_criteria(depth=5000, **inp) == ("", "")
    assert raschii.check_breaking_criteria(depth=-1, **inp) == ("", "")

    wave1 = WaveModel(depth=5000, **inp)
    wave2 = WaveModel(depth=-1, **inp)

    assert abs(wave1.length - wave2.length) < 1e-4
    assert abs(wave1.period - wave2.period) < 1e-4
    assert abs(wave1.c - wave2.c) < 1e-4


@pytest.mark.parametrize("order", [1, 2, 3, 5])
@pytest.mark.parametrize("give_period", [True, False])
def test_infinite_depth_stokes(order: int, give_period: bool):
    import raschii

    WaveModel, _AirModel = raschii.get_wave_model("Stokes")

    if give_period:
        inp = dict(height=12, period=15)
    else:
        inp = dict(height=12, length=350)

    assert raschii.check_breaking_criteria(depth=5000, **inp) == ("", "")
    assert raschii.check_breaking_criteria(depth=-1, **inp) == ("", "")

    inp["N"] = order

    wave1 = WaveModel(depth=5000, **inp)
    wave2 = WaveModel(depth=-1, **inp)

    assert abs(wave1.length - wave2.length) < 1e-4
    assert abs(wave1.period - wave2.period) < 1e-4
    assert abs(wave1.c - wave2.c) < 1e-4


@pytest.mark.parametrize("order", [1, 5, 10])
@pytest.mark.parametrize("give_period", [True, False])
def test_infinite_depth_fenton(order: int, give_period: bool):
    import raschii

    WaveModel, _AirModel = raschii.get_wave_model("Fenton")

    if give_period:
        inp = dict(height=12, period=15)
    else:
        inp = dict(height=12, length=350)

    assert raschii.check_breaking_criteria(depth=5000, **inp) == ("", "")
    assert raschii.check_breaking_criteria(depth=-1, **inp) == ("", "")

    inp["N"] = order

    wave1 = WaveModel(depth=5000, **inp)
    wave2 = WaveModel(depth=-1, **inp)

    assert abs(wave1.length - wave2.length) < 1e-4
    assert abs(wave1.period - wave2.period) < 1e-4
    assert abs(wave1.c - wave2.c) < 1e-4


def test_fenton_velocity_deep_water_is_finite():
    import numpy as np

    import raschii

    wave = raschii.FentonWave(height=3.0, depth=-1, length=100.0, N=10)
    z_surf = 25 * wave.length  # infinite depth is modelled as 25 wave lengths
    z = np.array([z_surf - 1.0, z_surf - 20.0])
    vel = wave.velocity(np.array([0.0, 10.0]), z, all_points_wet=True)
    assert np.all(np.isfinite(vel))

    # Matches a deep but finite depth
    deep = raschii.FentonWave(height=3.0, depth=5000.0, length=100.0, N=10)
    v1 = wave.velocity(0.0, z_surf - 1.0, all_points_wet=True)
    v2 = deep.velocity(0.0, deep.depth - 1.0, all_points_wet=True)
    assert np.allclose(v1, v2, atol=1e-6)
