import os
import math

import raschii
from utils import skip_swd_uninstalled


@skip_swd_uninstalled
def test_swd_stokes(tmpdir):
    from spectral_wave_data import SpectralWaveData

    dir_swd = str(tmpdir.mkdir("swd"))
    file_swd = os.path.join(dir_swd, "stokes.swd")

    height = 5.0
    depth = 15.0
    length = 200
    nperiods = 0.4
    dt = 0.1
    norder = 5

    WaveModel, AirModel = raschii.get_wave_model("Stokes")
    wave = WaveModel(height=height, depth=depth, length=length, N=norder)
    wave.write_swd(file_swd, dt=dt, nperiods=nperiods)
    # raschii and swd apply the same definitions for initial x, t and wave propagation direction...
    swd = SpectralWaveData(file_swd, x0=0.0, y0=0.0, t0=0.0, beta=0.0)

    assert swd["prog"].startswith("raschii-")
    assert swd["n"] == norder
    tmax = swd["tmax"]

    for t_swd in [0.0, 15 * dt]:  # No interpolation round-off in swd at exact time steps
        assert t_swd <= tmax

        x_swd = 0.4 * length
        y_swd = 0.0
        z_swd = -0.3 * depth

        swd.update_time(t_swd)
        zs_swd = swd.elev(x_swd, y_swd)
        v_swd = swd.grad_phi(x_swd, y_swd, z_swd)

        t_raschii = t_swd
        x_raschii = (x_swd,)
        z_raschii = (z_swd + wave.depth,)
        zs_raschii = wave.surface_elevation(x=x_raschii, t=t_raschii)
        zs_raschii_wl = zs_raschii[0] - wave.depth
        v_raschii = wave.velocity(x=x_raschii, z=z_raschii, t=t_raschii, all_points_wet=True)

        eps_r = 1.0e-6  # Swd files store data in float (single) precision
        eps_a = 1.0e-6

        assert math.isclose(v_swd.x, v_raschii[0, 0], rel_tol=eps_r, abs_tol=eps_a)
        assert math.isclose(v_swd.y, 0.0, rel_tol=eps_r, abs_tol=eps_a)
        assert math.isclose(v_swd.z, v_raschii[0, 1], rel_tol=eps_r, abs_tol=eps_a)
        assert math.isclose(zs_swd, zs_raschii_wl, rel_tol=eps_r, abs_tol=eps_a)

    swd.close()
