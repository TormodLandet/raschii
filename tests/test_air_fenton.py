import numpy


def test_fenton_air_with_fenton():
    from raschii import FentonAirPhase, FentonWave
    height = 10.0
    depth = 200.0
    length = 100.0
    N = 5
    air = FentonAirPhase(100)
    FentonWave(height, depth, length, N, air=air)
    
    # Compare velocities with numerical differentiation of the stream function
    eps = 1e-7
    for x in numpy.linspace(0, length, 21):
        z = depth + height / 2
        vel = air.velocity(x, z)
        sf0 = air.stream_function(x, z, frame='c')
        sfX = air.stream_function(x + eps, z, frame='c')
        sfZ = air.stream_function(x, z + eps, frame='c')
        assert vel.shape == (1, 2) and sf0.shape == (1,) and sfX.shape == (1,)
        sfvel_x = (sfZ[0] - sf0) / eps
        sfvel_z = -(sfX[0] - sf0) / eps
        print('x: %r, z: %r, vel: %r, vel_num: %r' %
              (x, z, vel, (sfvel_x, sfvel_z)))
        assert abs(vel[0, 0] - sfvel_x) < 1e-5
        assert abs(vel[0, 1] - sfvel_z) < 1e-5
