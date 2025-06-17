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
    assert abs(wave1.T - wave2.T) < 1e-4
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
    assert abs(wave1.T - wave2.T) < 1e-4
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
    assert abs(wave1.T - wave2.T) < 1e-4
    assert abs(wave1.c - wave2.c) < 1e-4
