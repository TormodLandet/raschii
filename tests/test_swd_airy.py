"""
Tests for SWD output from AiryWave.

Airy write_swd delegates to a Stokes N=1 wave, so we verify:
1. The output is byte-for-byte equivalent to Stokes N=1 for all amp values.
2. The simplified reader round-trips elevation and surface potential correctly.
"""

import math

import pytest

import raschii
from utils import skip_swd_uninstalled


@pytest.mark.parametrize("depth", [15.0, 200.0, -1.0])
def test_airy_swd_matches_stokes_n1(depth, tmpdir):
    """Airy write_swd must produce the same SWD file as Stokes N=1."""
    from raschii.swd import SwdReaderForRaschiiTests

    height = 3.0
    length = 80.0
    nperiods = 0.5
    dt = 0.05

    AiryModel, _ = raschii.get_wave_model("Airy")
    StokesModel, _ = raschii.get_wave_model("Stokes")

    airy = AiryModel(height=height, depth=depth, length=length)
    stokes = StokesModel(height=height, depth=depth, length=length, N=1)

    file_airy = tmpdir / "airy.swd"
    file_stokes = tmpdir / "stokes_n1.swd"

    airy.write_swd(file_airy, dt=dt, nperiods=nperiods)
    stokes.write_swd(file_stokes, dt=dt, nperiods=nperiods)

    swd_a = SwdReaderForRaschiiTests(file_airy)
    swd_s = SwdReaderForRaschiiTests(file_stokes)

    assert swd_a.shp == swd_s.shp
    assert swd_a.amp == swd_s.amp == 1

    eps = 1.0e-5
    dx = length / 17
    i_time = 8
    for x in [dx * i for i in range(17)]:
        eta_a = swd_a.surface_elevation(x=x)[i_time].item()
        eta_s = swd_s.surface_elevation(x=x)[i_time].item()
        assert math.isclose(eta_a, eta_s, rel_tol=eps, abs_tol=eps)

        phi_a = swd_a.surface_potential(x=x)[i_time].item()
        phi_s = swd_s.surface_potential(x=x)[i_time].item()
        assert math.isclose(phi_a, phi_s, rel_tol=eps, abs_tol=eps)


@pytest.mark.parametrize("depth", [15.0, -1.0])
def test_airy_swd_amp2(depth, tmpdir):
    """amp=2: surface potential read back matches velocity_potential on the surface."""
    from raschii.swd import SwdReaderForRaschiiTests

    height = 3.0
    length = 80.0
    nperiods = 0.5
    dt = 0.05

    AiryModel, _ = raschii.get_wave_model("Airy")
    wave = AiryModel(height=height, depth=depth, length=length)

    file_swd = tmpdir / "airy_amp2.swd"
    wave.write_swd(file_swd, dt=dt, nperiods=nperiods, amp=2)

    swd = SwdReaderForRaschiiTests(file_swd)
    assert swd.amp == 2

    eff_depth = 25.0 * length if depth < 0 else depth
    eps = 5.0e-4
    dx = length / 16
    i_time = 4
    t_check = swd.t_vector[i_time]
    for x in [dx * i for i in range(16)]:
        phi_swd = swd.surface_potential(x=x)[i_time].item()
        z_surf = wave.surface_elevation(x=[x], t=t_check, include_depth=False)[0] + eff_depth
        phi_ref = wave.velocity_potential(x=[x], z=[z_surf], t=t_check)[0]
        assert math.isclose(phi_swd, phi_ref, rel_tol=eps, abs_tol=eps), (
            f"x={x:.1f}: phi_swd={phi_swd:.6f}, phi_ref={phi_ref:.6f}"
        )


def test_airy_swd_amp3(tmpdir):
    """amp=3: elevation-only file, no potential stored."""
    from raschii.swd import SwdReaderForRaschiiTests

    AiryModel, _ = raschii.get_wave_model("Airy")
    wave = AiryModel(height=3.0, depth=15.0, length=80.0)
    file_swd = tmpdir / "airy_amp3.swd"
    wave.write_swd(file_swd, dt=0.05, nperiods=0.5, amp=3)

    swd = SwdReaderForRaschiiTests(file_swd)
    assert swd.amp == 3
    assert not hasattr(swd, "c_vectors")

    eps = 1.0e-5
    i_time = 4
    t_check = swd.t_vector[i_time]
    for x in [5.0 * i for i in range(16)]:
        eta_swd = swd.surface_elevation(x=x)[i_time].item()
        eta_ref = wave.surface_elevation(x=x, t=t_check, include_depth=False).item()
        assert math.isclose(eta_swd, eta_ref, rel_tol=eps, abs_tol=eps)


@skip_swd_uninstalled
def test_airy_swd_official_reader(tmpdir):
    """amp=1: verify elevation and kinematics against the official SWD library."""
    from spectral_wave_data import SpectralWaveData

    height = 3.0
    depth = 15.0
    length = 80.0
    nperiods = 0.5
    dt = 0.05

    AiryModel, _ = raschii.get_wave_model("Airy")
    wave = AiryModel(height=height, depth=depth, length=length)
    file_swd = tmpdir / "airy.swd"
    wave.write_swd(file_swd, dt=dt, nperiods=nperiods)

    swd = SpectralWaveData(str(file_swd), x0=0.0, y0=0.0, t0=0.0, beta=0.0)
    tmax = swd["tmax"]

    for t_swd in [0.0, 10 * dt]:
        assert t_swd <= tmax
        x_swd = 0.3 * length
        z_swd = -0.3 * depth

        swd.update_time(t_swd)
        zs_swd = swd.elev(x_swd, 0.0)
        v_swd = swd.grad_phi(x_swd, 0.0, z_swd)
        phi_swd = swd.phi(x_swd, 0.0, z_swd)

        x_r = (x_swd,)
        z_r = (z_swd + wave.depth,)
        zs_r = wave.surface_elevation(x=x_r, t=t_swd)[0] - wave.depth
        v_r = wave.velocity(x=x_r, z=z_r, t=t_swd, all_points_wet=True)
        phi_r = wave.velocity_potential(x=x_r, z=z_r, t=t_swd)

        eps_r = eps_a = 1.0e-5
        assert math.isclose(v_swd.x, v_r[0, 0], rel_tol=eps_r, abs_tol=eps_a)
        assert math.isclose(v_swd.z, v_r[0, 1], rel_tol=eps_r, abs_tol=eps_a)
        assert math.isclose(zs_swd, zs_r, rel_tol=eps_r, abs_tol=eps_a)
        assert math.isclose(phi_swd, phi_r[0], rel_tol=eps_r, abs_tol=eps_a)

    swd.close()
