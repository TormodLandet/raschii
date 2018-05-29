import numpy


def test_fenton_air_with_fenton(plot=False):
    from raschii import FentonAirPhase, FentonWave
    height = 10.0
    depth = 200.0
    length = 100.0
    height_air = 100.0
    N = 5
    time = 0.0
    air = FentonAirPhase(height_air)
    fwave = FentonWave(height, depth, length, N, air=air)
    
    # Locations to check
    xpos = numpy.linspace(-length / 2, length / 2, 101)
    zpos = numpy.linspace(depth - height / 2, depth + 2.75 * height, 101)
    X, Z = numpy.meshgrid(xpos, zpos)
    xr = X.ravel()
    zr = Z.ravel()
    eps = 1e-7
    
    # Compare velocities with numerical differentiation of the stream function
    avel = air.velocity(xr, zr, time)
    sf0 = air.stream_function(xr, zr, time, frame='c')
    sfX = air.stream_function(xr + eps, zr, time, frame='c')
    sfZ = air.stream_function(xr, zr + eps, time, frame='c')
    sfvel_x = (sfZ - sf0) / eps
    sfvel_z = -(sfX - sf0) / eps
    err = abs(avel[:, 0] - sfvel_x) + abs(avel[:, 1] - sfvel_z)
    
    maxi = err.argmax()
    xmax = xr[maxi]
    zmax = zr[maxi]
    max_vel_err = err[maxi]
    print('\nThe maximum velocity error is', max_vel_err)
    print('The location is x/lambda = %.5f and (z - D)/H = %.5f' %
          (xmax / length, (zmax - depth) / height))
    print('The expected velocity is %r, got %r' % (tuple(avel[maxi]),
                                                   (sfvel_x[maxi],
                                                    sfvel_z[maxi])))
    assert max_vel_err < 1e-5
    
    # Check that the blended velocity field is divergence free
    totvel = fwave.velocity(xr, zr, time, all_points_wet=False)
    velsdx = fwave.velocity(xr + eps, zr, time, all_points_wet=False)
    velsdz = fwave.velocity(xr, zr + eps, time, all_points_wet=False)
    div = (velsdx[:, 0] - totvel[:, 0] + velsdz[:, 1] - totvel[:, 1]) / eps
    adiv = abs(div)
    
    if plot:
        from matplotlib import pyplot
        c = pyplot.contourf(X, Z, adiv.reshape(X.shape))
        pyplot.colorbar(c)
        pyplot.plot(xpos, fwave.surface_elevation(xpos, time))
        pyplot.show()
    
    maxi = adiv.argmax()
    xmax = xr[maxi]
    zmax = zr[maxi]
    max_abs_div = adiv[maxi]
    print('\nThe maximum absolute divergence is', max_abs_div)
    print('The location is x/lambda = %.5f and (z - D)/H = %.5f' %
          (xmax / length, (zmax - depth) / height))
    print('The velocity at the location is %r' % (tuple(totvel[maxi]), ))
    assert max_abs_div < 1e-5
