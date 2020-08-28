import os
import pytest


def test_readme_example():
    import raschii

    fwave = raschii.FentonWave(height=0.25, depth=0.5, length=2.0, N=20)
    print(fwave.surface_elevation(x=0))
    print(fwave.surface_elevation(x=[0, 0.1, 0.2, 0.3]))
    print(fwave.velocity(x=0, z=0.2))


@pytest.mark.parametrize("wave_model_name", ["Fenton", "Stokes"])
def test_readme_swd_waves(wave_model_name, tmpdir):
    from raschii.cmd.swd import write_swd

    swd_file = os.path.join(tmpdir.mkdir("swd"), "readme_test_%s.swd" % wave_model_name)

    write_swd(swd_file, wave_model_name, 0.2, 1.5, 2, N=5, dt=0.01, tmax=10)


def test_user_doc():
    import raschii

    WaveModel, AirModel = raschii.get_wave_model("Fenton")
    wave = WaveModel(height=12, depth=200, length=100, N=5)

    elev = wave.surface_elevation([0.0, 10.0, 20.0])
    print(elev)

    vel = wave.velocity(0.0, elev[0])
    print(vel)

    WaveModel, AirModel = raschii.get_wave_model("Fenton", "FentonAir")
    air = AirModel(height=100, blending_height=20)
    wave = WaveModel(height=12, depth=200, length=100, N=5, air=air)

    vel = wave.velocity(0.0, 208.0)
    print(vel)
