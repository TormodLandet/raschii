import numpy
import pytest
from raschii import get_wave_model
from jit_helper import jit_compile


@pytest.fixture(params=['Airy', 'Fenton'])
def wave_model(request):
    if request.param == 'Airy':
        model = get_wave_model('Airy')[0](height=1, depth=10, length=20)
    elif request.param == 'Fenton':
        model = get_wave_model('Fenton')[0](height=1, depth=10, length=20, N=5)
        
    if ('stream_function' in request.node.name and
            not hasattr(model, 'stream_function_cpp')):
        raise pytest.skip('Stream function C++ missing from %s' % request.param)
    
    if ('slope' in request.node.name and
            not hasattr(model, 'slope_cpp')):
        raise pytest.skip('Surface slope C++ missing from %s' % request.param)
    
    return model


def test_cpp_jit(tmpdir):
    # The example from the pybind11 docs
    cpp_code = """
    #include <pybind11/pybind11.h>
    
    int add(int i, int j) {
        return i + j;
    }
    
    PYBIND11_MODULE(MODNAME, m) {
        m.doc() = "pybind11 example plugin"; // optional module docstring
    
        m.def("add", &add, "A function which adds two numbers");
    }
    """
    cache_dir = tmpdir.ensure('jit_cache', dir=True)
    mod = jit_compile(cpp_code, cache_dir)
    assert mod.add(5, 37) == 42


def test_cpp_vs_py_elevation(tmpdir, wave_model):
    cpp_wrapper = """
    #define _USE_MATH_DEFINES
    #include <vector>
    #include <cmath>
    #include <pybind11/pybind11.h>
    #include <pybind11/stl.h>
    
    using namespace std;
    const double pi = M_PI;
    
    double elevation(vector<double> x, double t=0.0) {
        double value = CODE_GOES_HERE;
        return value;
    }
    
    namespace py = pybind11;
    PYBIND11_MODULE(MODNAME, m) {
        m.def("elevation", &elevation, py::arg("x"), py::arg("t")=0.0);
    }
    """
    cache_dir = tmpdir.ensure('jit_cache', dir=True)
    
    # Check that the wave model produces the same results in C++ and Python
    cpp = wave_model.elevation_cpp()
    mod = jit_compile(cpp_wrapper.replace('CODE_GOES_HERE', cpp), cache_dir)
    
    for pos in ([0.0, 0.0], [4.0, 7.0]):
        e_cpp = mod.elevation(pos, t=1.3)
        e_py = wave_model.surface_elevation(pos[0], t=1.3)[0]
        err = abs(e_cpp - e_py)
        print(wave_model.__class__.__name__, pos, e_cpp, e_py, err)
        assert err < 1e-14


def test_cpp_vs_py_velocity(tmpdir, wave_model):
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
    
    # Check that the wave model produces the same results in C++ and Python
    cppx, cppz = wave_model.velocity_cpp()
    cpp = cpp_wrapper.replace('CODE_X_GOES_HERE', cppx)\
                     .replace('CODE_Z_GOES_HERE', cppz)\
                     .replace('x[2]', 'x[1]')
    mod = jit_compile(cpp, cache_dir)
    
    t = 4.2
    for pos in ([-10.0, 5.0], [5.0, 6.0]):
        vx_cpp = mod.vel_x(pos, t)
        vz_cpp = mod.vel_z(pos, t)
        vx_py, vz_py = wave_model.velocity(pos[0], pos[1], t)[0]
        
        # Compute the relative error
        rerr = abs((vx_cpp - vx_py) / vx_py) + abs((vz_cpp - vz_py) / vz_py)
        print(wave_model.__class__.__name__, pos, (vx_cpp, vz_cpp),
              (vx_py, vz_py), rerr)
        assert rerr < 1e-2


def test_cpp_vs_py_stream_function(tmpdir, wave_model):
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
    
    # Check that the wave model produces the same results in C++ and Python
    cpp = wave_model.stream_function_cpp(frame='c')
    cpp = cpp_wrapper.replace('CODE_GOES_HERE', cpp)\
                     .replace('x[2]', 'x[1]')
    mod = jit_compile(cpp, cache_dir)
    
    t = -23.9
    for pos in ([-10.0, 5.0], [5.0, 6.0]):
        sf_cpp = mod.sfunc(pos, t)
        sf_py = wave_model.stream_function(pos[0], pos[1], t, frame='c')[0]
        
        # Compute the relative error
        rerr = abs((sf_cpp - sf_py) / sf_py)
        print(wave_model.__class__.__name__, pos, sf_cpp, sf_py, rerr)
        assert rerr < 1e-2


def test_cpp_vs_py_slope(tmpdir, wave_model):
    cpp_wrapper = """
    #define _USE_MATH_DEFINES
    #include <vector>
    #include <cmath>
    #include <pybind11/pybind11.h>
    #include <pybind11/stl.h>
    
    using namespace std;
    const double pi = M_PI;
    
    double slope(vector<double> x, double t=0.0) {
        double value = CODE_GOES_HERE;
        return value;
    }
    
    namespace py = pybind11;
    PYBIND11_MODULE(MODNAME, m) {
        m.def("slope", &slope, py::arg("x"), py::arg("t")=0.0);
    }
    """
    cache_dir = tmpdir.ensure('jit_cache', dir=True)
    
    # Check that the wave model produces the same results in C++ and Python
    cpp = wave_model.slope_cpp()
    cpp = cpp_wrapper.replace('CODE_GOES_HERE', cpp)\
                     .replace('x[2]', 'x[1]')
    mod = jit_compile(cpp, cache_dir)
    
    t = 100.0
    for x in numpy.linspace(0, 20, 201):
        slope_cpp = mod.slope([x, 0], t)
        slope_py = wave_model.surface_slope(x, t)[0]
        
        # Compute the relative error
        rerr = abs((slope_cpp - slope_py) / slope_py)
        print(wave_model.__class__.__name__, x, slope_cpp, slope_py, rerr)
        assert rerr < 1e-14
