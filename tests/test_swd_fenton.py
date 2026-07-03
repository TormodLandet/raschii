import math

import pytest

import raschii
from utils import skip_swd_uninstalled


@skip_swd_uninstalled
def test_swd_fenton(tmpdir):
    from spectral_wave_data import SpectralWaveData

    file_swd = tmpdir / "fenton.swd"
    height = 5.0
    depth = 15.0
    length = 200
    nperiods = 0.4
    dt = 0.1
    ncomp = 30

    WaveModel, AirModel = raschii.get_wave_model("Fenton")
    wave = WaveModel(height=height, depth=depth, length=length, N=ncomp)
    wave.write_swd(file_swd, dt=dt, nperiods=nperiods)
    # raschii and swd apply the same definitions for initial x, t and wave propagation direction...
    swd = SpectralWaveData(str(file_swd), x0=0.0, y0=0.0, t0=0.0, beta=0.0)

    assert swd["prog"].startswith("raschii-")
    assert swd["n"] == ncomp
    tmax = swd["tmax"]

    for t_swd in [0.0, 15 * dt]:  # No interpolation round-off in swd at exact time steps
        assert t_swd <= tmax

        x_swd = 0.4 * length
        y_swd = 0.0
        z_swd = -0.3 * depth

        swd.update_time(t_swd)
        zs_swd = swd.elev(x_swd, y_swd)
        v_swd = swd.grad_phi(x_swd, y_swd, z_swd)
        phi_swd = swd.phi(x_swd, y_swd, z_swd)

        t_raschii = t_swd
        x_raschii = (x_swd,)
        z_raschii = (z_swd + wave.depth,)
        zs_raschii = wave.surface_elevation(x=x_raschii, t=t_raschii)
        zs_raschii_wl = zs_raschii[0] - wave.depth
        v_raschii = wave.velocity(x=x_raschii, z=z_raschii, t=t_raschii, all_points_wet=True)
        phi_raschii = wave.velocity_potential(x=x_raschii, z=z_raschii, t=t_raschii)

        eps_r = 1.0e-6  # Swd files store data in float (single) precision
        eps_a = 1.0e-6

        assert math.isclose(v_swd.x, v_raschii[0, 0], rel_tol=eps_r, abs_tol=eps_a)
        assert math.isclose(v_swd.y, 0.0, rel_tol=eps_r, abs_tol=eps_a)
        assert math.isclose(v_swd.z, v_raschii[0, 1], rel_tol=eps_r, abs_tol=eps_a)
        assert math.isclose(zs_swd, zs_raschii_wl, rel_tol=eps_r, abs_tol=eps_a)
        assert math.isclose(phi_swd, phi_raschii[0], rel_tol=eps_r, abs_tol=eps_a)

    swd.close()


@pytest.mark.parametrize("depth", [15.0, 200, -1.0])
def test_swd_fenton_simplified_reader(depth, tmpdir):
    from raschii.swd import SwdReaderForRaschiiTests

    file_swd = tmpdir / "fenton.swd"
    height = 5.0
    length = 200
    nperiods = 0.4
    dt = 0.1
    norder = 5

    WaveModel, AirModel = raschii.get_wave_model("Fenton")
    wave = WaveModel(height=height, depth=depth, length=length, N=norder)
    wave.write_swd(file_swd, dt=dt, nperiods=nperiods)

    eps_r = 1.0e-6  # Swd files store data in float (single) precision
    eps_a = 1.0e-6

    swd = SwdReaderForRaschiiTests(file_swd)
    assert math.isclose(swd.dk, wave.k, rel_tol=eps_r, abs_tol=eps_a)
    assert math.isclose(swd.depth, wave.depth, rel_tol=eps_r, abs_tol=eps_a)

    if depth < 1:
        assert swd.shp == 1, "Expected infinite depth SWD file"
    else:
        assert swd.shp == 2, "Expected finite depth SWD file"

    # Compare surface elevations
    eps_r = eps_a = 1.0e-5
    dx = length / 23
    i_time = 13
    t_check = swd.t_vector[i_time]
    for x in [dx * i for i in range(30)]:
        swdfile_eta = swd.surface_elevation(x=x)[i_time].item()
        raschii_eta = wave.surface_elevation(x=x, t=t_check, include_depth=False).item()

        assert math.isclose(swdfile_eta, raschii_eta, rel_tol=eps_r, abs_tol=eps_a)


@pytest.mark.parametrize("depth", [15.0, -1.0])
def test_swd_fenton_amp2(depth, tmpdir):
    """amp=2: surface potential coefficients round-trip via the simplified reader."""
    from raschii.swd import SwdReaderForRaschiiTests

    file_swd = tmpdir / "fenton_amp2.swd"
    height = 5.0
    length = 200.0
    nperiods = 0.4
    dt = 0.1
    norder = 5

    WaveModel, _ = raschii.get_wave_model("Fenton")
    wave = WaveModel(height=height, depth=depth, length=length, N=norder)
    wave.write_swd(file_swd, dt=dt, nperiods=nperiods, amp=2)

    swd = SwdReaderForRaschiiTests(file_swd)
    assert swd.amp == 2

    # Compare surface potential: SWD reconstruction vs wave.velocity_potential on surface
    eps_r = eps_a = 5.0e-4  # oversampling=8 and N=5 harmonics limits precision
    dx = length / 20
    i_time = 5
    t_check = swd.t_vector[i_time]

    eff_depth = 25.0 * length if depth < 0 else depth
    for x in [dx * i for i in range(20)]:
        phi_swd = swd.surface_potential(x=x)[i_time].item()
        z_surf = wave.surface_elevation(x=[x], t=t_check, include_depth=False)[0] + eff_depth
        phi_ref = wave.velocity_potential(x=[x], z=[z_surf], t=t_check)[0]
        assert math.isclose(phi_swd, phi_ref, rel_tol=eps_r, abs_tol=eps_a), (
            f"x={x:.1f}: phi_swd={phi_swd:.6f}, phi_ref={phi_ref:.6f}"
        )


@pytest.mark.parametrize("depth", [15.0, -1.0])
def test_swd_fenton_amp3(depth, tmpdir):
    """amp=3: elevation-only file — c data absent, elevation correct."""
    from raschii.swd import SwdReaderForRaschiiTests

    file_swd = tmpdir / "fenton_amp3.swd"
    height = 5.0
    length = 200.0
    nperiods = 0.4
    dt = 0.1

    WaveModel, _ = raschii.get_wave_model("Fenton")
    wave = WaveModel(height=height, depth=depth, length=length, N=5)
    wave.write_swd(file_swd, dt=dt, nperiods=nperiods, amp=3)

    swd = SwdReaderForRaschiiTests(file_swd)
    assert swd.amp == 3
    assert not hasattr(swd, "c_vectors"), "c_vectors should not be stored for amp=3"

    eps_r = eps_a = 1.0e-5
    dx = length / 23
    i_time = 7
    t_check = swd.t_vector[i_time]
    for x in [dx * i for i in range(20)]:
        swdfile_eta = swd.surface_elevation(x=x)[i_time].item()
        raschii_eta = wave.surface_elevation(x=x, t=t_check, include_depth=False).item()
        assert math.isclose(swdfile_eta, raschii_eta, rel_tol=eps_r, abs_tol=eps_a)
