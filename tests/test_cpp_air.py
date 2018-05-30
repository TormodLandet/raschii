import numpy
import pytest
from raschii import get_wave_model
from jit_helper import jit_compile


@pytest.fixture(params=['ConstantAir', 'FentonAir'])
def air_model(request):
    WaveClass, AirClass = get_wave_model('Fenton', request.param)
    air = AirClass(3.0, 2.0)
    height = 0.5
    depth = 7.0
    length = 20.0
    WaveClass(height, depth, length, N=5, air=air)
    return air, height, depth, length


def air_locations(height, depth, length, blending_height):
    # Locations to check
    top = depth + 2.75 * height
    if blending_height < 5 * height:
        top = max(depth + blending_height * 1.2, top)
    xpos = numpy.linspace(-length / 2, length / 2, 101)
    zpos = numpy.linspace(depth - height / 2, top, 101)
    X, Z = numpy.meshgrid(xpos, zpos)
    xr = X.ravel()
    zr = Z.ravel()
    return xr, zr


def check_results(xr, zr, expected, computed, name, tolerance):
    aerr = abs(expected - computed)
    maxi = aerr.argmax()
    xmax = xr[maxi]
    zmax = zr[maxi]
    max_abs_err = aerr[maxi]
    print('\nThe maximum error for %s is %r' % (name, max_abs_err))
    print('The location is x = %.5f and z = %.5f' % (xmax, zmax))
    print('The expected value here is %r' % expected[maxi])
    print('The computed value here is %r' % computed[maxi])
    assert max_abs_err < tolerance


def test_cpp_vs_py_air_stream_function(tmpdir, air_model):
    cpp_wrapper = """
    #define _USE_MATH_DEFINES
    #include <vector>
    #include <cmath>
    #include <pybind11/pybind11.h>
    #include <pybind11/stl.h>
    
    using namespace std;
    const double pi = M_PI;
    
    double stream_function(vector<double> x, double t=0.0) {
        double value = CODE_GOES_HERE;
        return value;
    }
    
    namespace py = pybind11;
    PYBIND11_MODULE(MODNAME, m) {
        m.def("sfunc", &stream_function, py::arg("x"), py::arg("t")=0.0);
    }
    """
    cache_dir = tmpdir.ensure('jit_cache', dir=True)
    air_model, height, depth, length = air_model
    
    # Check that the wave model produces the same results in C++ and Python
    cpp = air_model.stream_function_cpp(frame='c')
    cpp = cpp_wrapper.replace('CODE_GOES_HERE', cpp)\
                     .replace('x[2]', 'x[1]')
    mod = jit_compile(cpp, cache_dir)
    
    # The input values
    xr, zr = air_locations(height, depth, length, air_model.blending_height)
    t = -23.9
    
    # Compute the stream function using both the C++ and the Python versions
    sf_cpp = numpy.zeros_like(xr)
    for i in range(xr.size):
        sf_cpp[i] = mod.sfunc([xr[i], zr[i]], t)
    sf_py = air_model.stream_function(xr, zr, t, frame='c')
    
    # Check the results
    test_name = '%s C++ stream function' % air_model.__class__.__name__
    check_results(xr, zr, sf_py, sf_cpp, test_name, 1e-4)


def test_cpp_vs_py_air_velocity(tmpdir, air_model):
    cpp_wrapper = """
    #define _USE_MATH_DEFINES
    #include <vector>
    #include <cmath>
    #include <pybind11/pybind11.h>
    #include <pybind11/stl.h>
    
    using namespace std;
    const double pi = M_PI;
    
    double vel_x(vector<double> x, double t=0.0) {
        double value = CODE_X_GOES_HERE;
        return value;
    }
    
    double vel_z(vector<double> x, double t=0.0) {
        double value = CODE_Z_GOES_HERE;
        return value;
    }
    
    namespace py = pybind11;
    PYBIND11_MODULE(MODNAME, m) {
        m.def("vel_x", &vel_x, py::arg("x"), py::arg("t")=0.0);
        m.def("vel_z", &vel_z, py::arg("x"), py::arg("t")=0.0);
    }
    """
    cache_dir = tmpdir.ensure('jit_cache', dir=True)
    air_model, height, depth, length = air_model
    
    # Check that the wave model produces the same results in C++ and Python
    cppx, cppz = air_model.velocity_cpp()
    cpp = cpp_wrapper.replace('CODE_X_GOES_HERE', cppx)\
                     .replace('CODE_Z_GOES_HERE', cppz)\
                     .replace('x[2]', 'x[1]')
    mod = jit_compile(cpp, cache_dir)
    
    # The input values
    xr, zr = air_locations(height, depth, length, air_model.blending_height)
    t = 4.2
    
    # Compute the velocities using both the C++ and the Python versions
    vx_cpp = numpy.zeros_like(xr)
    vz_cpp = numpy.zeros_like(xr)
    for i in range(xr.size):
        vx_cpp[i] = mod.vel_x([xr[i], zr[i]], t)
        vz_cpp[i] = mod.vel_z([xr[i], zr[i]], t)
    vel_py = air_model.velocity(xr, zr, t)

    # Check the results
    test_name_x = '%s C++ x-velocity' % air_model.__class__.__name__
    test_name_z = '%s C++ z-velocity' % air_model.__class__.__name__
    check_results(xr, zr, vel_py[:, 0], vx_cpp, test_name_x, 1e-5)
    check_results(xr, zr, vel_py[:, 1], vz_cpp, test_name_z, 1e-5)
