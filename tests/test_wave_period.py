import numpy as np
from raschii import get_wave_model, WAVE_MODELS
import pytest

WAVE_MODEL_NAMES = list(WAVE_MODELS.keys())


def create_test_wave(wave_model_name):
    model_class = get_wave_model(wave_model_name)[0]

    wave_inp = dict(height=1, depth=10, length=20)
    if "N" in model_class.required_input:
        wave_inp["N"] = 5

    return model_class(**wave_inp)


@pytest.mark.parametrize("wave_model_name", WAVE_MODEL_NAMES)
def test_wave_period(wave_model_name):
    wave_model = create_test_wave(wave_model_name)
    depth = wave_model.depth

    N = 3000
    t = np.linspace(0, 100, N)
    x = 0
    eta = np.array([wave_model.surface_elevation(x, ti) for ti in t])

    # Find zero up-crossings
    up_crossings = []
    for i in range(N):
        if eta[i] > depth and eta[i - 1] <= depth:
            up_crossings.append(i)

    # Find periods
    periods = []
    for i in range(1, len(up_crossings)):
        it0 = up_crossings[i - 1]
        it1 = up_crossings[i]
        periods.append(t[it1] - t[it0])

    # Find the average wave period and check that it is accurate
    period = np.mean(periods)
    stdev_pst = np.std(periods) / period * 100
    print(wave_model_name, wave_model.T, period, stdev_pst)
    assert stdev_pst < 0.6, "The standard dev should small compared to the mean"

    # Check period of eta vs analytical period
    error_pst = abs(period - wave_model.T) / wave_model.T * 100
    print(error_pst)
    assert error_pst < 0.05, "The period should be close to the analytical"

    # Check that omega is defined and equal to 2pi/T
    omega = 2 * np.pi / period
    err_omega_pst = abs(omega - wave_model.omega) / wave_model.omega * 100
    print(err_omega_pst)
    assert err_omega_pst < 0.05, "The omegas do not match"


if __name__ == "__main__":
    for wmn in WAVE_MODEL_NAMES:
        test_wave_period(wmn)
