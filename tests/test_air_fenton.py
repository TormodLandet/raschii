import numpy
import pytest
from jit_helper import jit_compile


@pytest.fixture(params=['deep1', 'shallow1'])
def wave_with_air_model(request):
    N = 5
    time = 0.0
    blending_height = None
    plot = False  # for debugging only
    
    if request.param == 'deep1':
        height = 10.0
        depth = 200.0
        length = 100.0
        height_air = 100.0
    elif request.param == 'shallow1':
        height = 0.5
        depth = 7.0
        length = 20
        height_air = 3.0
        blending_height = 2
    
    from raschii import FentonAirPhase, FentonWave
    air = FentonAirPhase(height_air, blending_height)
    fwave = FentonWave(height, depth, length, N, air=air)
    return fwave, air, time, plot


def test_fenton_air_with_fenton(wave_with_air_model):
    # Get the wave from the fixture
    fwave, air, time, plot = wave_with_air_model
    length, depth, height = fwave.length, fwave.depth, fwave.height
    
    # Locations to check
    top = depth + 2.75 * height
    if air.blending_height < 5 * height:
        top = max(depth + air.blending_height * 1.2, top)
    xpos = numpy.linspace(-length / 2, length / 2, 101)
    zpos = numpy.linspace(depth - height / 2, top, 101)
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


def test_fenton_air_with_fenton_cpp_divergence(tmpdir, wave_with_air_model):
    # Get the wave from the fixture
    fwave, air, time, plot = wave_with_air_model
    length, depth, height = fwave.length, fwave.depth, fwave.height
    
    cpp_wrapper = """
    #define _USE_MATH_DEFINES
    #include <vector>
    #include <cmath>
    #include <pybind11/pybind11.h>
    #include <pybind11/stl.h>
    
    using namespace std;
    const double pi = M_PI;
    
    vector<double> vel(const vector<double> x, const vector<double> z,
                       const double t=0.0) {
        vector<double> velocity(x.size() * 2);
        for (int i = 0; i < x.size(); i++) {
            double xpos = x[i];
            double zpos = z[i];
            velocity[i * 2 + 0] = CODE_X_GOES_HERE;
            velocity[i * 2 + 1] = CODE_Z_GOES_HERE;
        }
        return velocity;
    }
    namespace py = pybind11;
    PYBIND11_MODULE(MODNAME, m) {
        m.def("vel", &vel, py::arg("xr"), py::arg("zr"), py::arg("t")=0.0);
    }
    """
    cache_dir = tmpdir.ensure('jit_cache', dir=True)
    
    # Check that the wave model produces the same results in C++ and Python
    cppx, cppz = fwave.velocity_cpp()
    
    cpp = cpp_wrapper.replace('CODE_X_GOES_HERE', cppx)\
                     .replace('CODE_Z_GOES_HERE', cppz)\
                     .replace('x[0]', 'xpos')\
                     .replace('x[2]', 'zpos')
    mod = jit_compile(cpp, cache_dir)
    
    # Locations to check
    top = depth + 2.75 * height
    if air.blending_height < 5 * height:
        top = max(depth + air.blending_height * 1.2, top)
    xpos = numpy.linspace(-length / 2, length / 2, 101)
    zpos = numpy.linspace(depth - height / 2, top, 101)
    X, Z = numpy.meshgrid(xpos, zpos)
    xr = X.ravel()
    zr = Z.ravel()
    eps = 1e-7
    
    # Check that the blended velocity field is divergence free
    totvel = numpy.asarray(mod.vel(xr, zr, time)).reshape((xr.size, 2))
    velsdx = numpy.asarray(mod.vel(xr + eps, zr, time)).reshape((xr.size, 2))
    velsdz = numpy.asarray(mod.vel(xr, zr + eps, time)).reshape((xr.size, 2))
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
    assert max_abs_div < 1e-4
