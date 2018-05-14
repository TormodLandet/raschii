from raschii import get_wave_model
from jit_helper import jit_compile


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


def test_cpp_vs_py_elevation(tmpdir):
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
    
    # Create wave models
    airy = get_wave_model('Airy')(height=1, depth=10, length=20)
    fenton = get_wave_model('Fenton')(height=1, depth=10, length=20, N=5)
    
    # Check that each model produces the same results in C++ and Python
    for model in [airy, fenton]:
        cpp = model.elevation_cpp()
        mod = jit_compile(cpp_wrapper.replace('CODE_GOES_HERE', cpp), cache_dir)
        
        for pos in ([0.0, 0.0], [4.0, 7.0]):
            e_cpp = mod.elevation(pos, t=1.3)
            e_py = model.surface_elevation(pos[0], t=1.3)[0]
            err = abs(e_cpp - e_py)
            print(model.__class__.__name__, pos, e_cpp, e_py, err)
            assert err < 1e-16
