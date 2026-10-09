import numpy as np
import pytest

import raschii
from raschii import FentonWave, RaschiiError


def test_zero_current_is_unchanged():
    w0 = FentonWave(height=3.0, depth=10.0, period=8.0, N=10)
    w1 = FentonWave(height=3.0, depth=10.0, period=8.0, N=10, current=0.0, current_criterion="stokes")
    assert w0.current == 0.0
    assert abs(w0.period - 8.0) < 1e-3
    assert abs(w0.length - w1.length) < 1e-3 * w0.length or w1.current_stokes == 0.0


def test_against_fourier_program_no_current():
    # Fenton's Fourier program, Data.dat test wave: H/d=0.3, T sqrt(g/d)=10, N=20, no current
    w = FentonWave(height=0.3, depth=1.0, period=10.0, N=20, g=1.0)
    assert abs(w.length - 9.6897) < 1e-3
    assert abs(w.c - 0.9690) < 1e-4
    assert abs(w.current_stokes - 0.0111) < 1e-4  # "Stokes current" in Solution.res
    assert abs(w.current) < 1e-12


@pytest.mark.parametrize("current", [-0.5, 0.5, 1.5])
@pytest.mark.parametrize("criterion", ["euler", "stokes"])
def test_period_is_observed_period_with_current(current, criterion):
    w = FentonWave(
        height=2.0, depth=20.0, period=10.0, N=10, current=current, current_criterion=criterion
    )
    assert abs(w.period - 10.0) < 1e-3
    assert abs(w.length / w.period - w.phase_speed) < 1e-3
    assert abs(w.current_input - current) < 1e-12
    if criterion == "euler":
        assert abs(w.current - current) < 1e-12
    else:
        assert abs(w.current_stokes - current) < 1e-3

    # A following current lengthens the wave, an opposing one shortens it
    w0 = FentonWave(height=2.0, depth=20.0, period=10.0, N=10)
    assert (w.length > w0.length) == (current > 0)


def test_euler_current_is_a_galilean_shift():
    U = 0.7
    w0 = FentonWave(height=2.0, depth=10.0, length=60.0, N=10)
    w1 = FentonWave(height=2.0, depth=10.0, length=60.0, N=10, current=U)
    x = np.linspace(0, 60, 5)
    z = np.full(5, 9.0)
    t = 1.3
    v0 = w0.velocity(x - U * t, z, t, all_points_wet=True)
    v1 = w1.velocity(x, z, t, all_points_wet=True)
    assert np.allclose(v1[:, 0], v0[:, 0] + U, atol=1e-10)
    assert np.allclose(v1[:, 1], v0[:, 1], atol=1e-10)
    assert np.allclose(w1.surface_elevation(x, t), w0.surface_elevation(x - U * t, t))


def test_invalid_criterion_and_air():
    with pytest.raises(RaschiiError):
        FentonWave(height=1.0, depth=10.0, length=50.0, current_criterion="bad")
    with pytest.raises(RaschiiError):
        air = raschii.FentonAirPhase(10.0, 0.5)
        FentonWave(height=1.0, depth=10.0, length=50.0, current=1.0, air=air)


def test_current_too_strong_against_waves():
    with pytest.raises(RaschiiError):
        FentonWave(height=1.0, depth=10.0, period=8.0, current=-20.0)
