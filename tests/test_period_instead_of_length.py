import pytest


def test_period_instead_of_length_airy():
    import raschii

    WaveModel, _AirModel = raschii.get_wave_model("Airy")
    wave = WaveModel(height=12, depth=200, period=15)

    assert abs(wave.length - 350.7521) < 1e-2
    assert abs(wave.T - 15.0) < 1e-2


@pytest.mark.parametrize("order", [1, 2, 3, 5])
def test_period_instead_of_length_stokes(order: int):
    import raschii

    WaveModel, _AirModel = raschii.get_wave_model("Stokes")
    wave = WaveModel(height=12, depth=200, period=15, N=5)

    assert abs(wave.length - 354.7048) < 1e-2
    assert abs(wave.T - 15.0) < 1e-2


@pytest.mark.parametrize("order", [1, 2, 3, 5, 10, 20])
def test_period_instead_of_length_fenton(order: int):
    import raschii

    WaveModel, _AirModel = raschii.get_wave_model("Fenton")
    wave = WaveModel(height=12, depth=200, period=15, N=5)

    assert abs(wave.length - 354.7048) < 1e-2
    assert abs(wave.T - 15.0) < 1e-2
