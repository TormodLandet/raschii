import numpy
import pytest
from raschii import get_wave_model
from jit_helper import jit_compile
from utils import skipif_no_compile


@pytest.fixture(params=["Airy", "Stokes", "Fenton"])
def wave_model(request):
    if request.param == "Airy":
        model = get_wave_model("Airy")[0](height=1, depth=10, length=20)
    elif request.param == "Stokes":
        model = get_wave_model("Stokes")[0](height=1, depth=10, length=20, N=5)
    elif request.param == "Fenton":
        model = get_wave_model("Fenton")[0](height=1, depth=10, length=20, N=5)

    for tname, mname in [
        ("elevation", "elevation_cpp"),
        ("velocity", "velocity_cpp"),
        ("stream_function", "stream_function_cpp"),
        ("slope", "slope_cpp"),
    ]:
        if tname in request.node.name and not hasattr(model, mname):
            pytest.xfail("Missing %sWave.%s method" % (request.param, mname))

    return model


def wave_locations(wave_model):
    """
    Define locations to check in the tests
    """
    height, depth = wave_model.height, wave_model.depth
    xpos = numpy.linspace(-wave_model.length / 2, wave_model.length / 2, 101)
    zpos = numpy.linspace(depth - height * 2, depth + height, 101)
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
    print("\nThe maximum error for %s is %r" % (name, max_abs_err))
    print("The location is x = %.5f and z = %.5f" % (xmax, zmax))
    print("The expected value here is %r" % expected[maxi])
    print("The computed value here is %r" % computed[maxi])
    assert max_abs_err < tolerance


@skipif_no_compile
def test_cpp_jit(tmpdir):
    """
    Test the example from the pybind11 docs. If this fails then the build tool
    chain is at fault and not Raschii itself
    """
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
    cache_dir = tmpdir.ensure("jit_cache", dir=True)
    mod = jit_compile(cpp_code, cache_dir)
    assert mod.add(5, 37) == 42


@skipif_no_compile
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
    cache_dir = tmpdir.ensure("jit_cache", dir=True)

    # Check that the wave model produces the same results in C++ and Python
    cpp = wave_model.elevation_cpp()
    assert "x[2]" not in cpp
    mod = jit_compile(cpp_wrapper.replace("CODE_GOES_HERE", cpp), cache_dir)

    # The input values
    xr, zr = wave_locations(wave_model)
    t = 1.3

    # Compute the elevation using both the C++ and the Python versions
    e_cpp = numpy.zeros_like(xr)
    for i in range(xr.size):
        e_cpp[i] = mod.elevation([xr[i], zr[i]], t)
    e_py = wave_model.surface_elevation(xr, t)

    # Check the results
    test_name = "%s C++ elevation" % wave_model.__class__.__name__
    check_results(xr, zr, e_py, e_cpp, test_name, 1e-14)


@skipif_no_compile
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
    cache_dir = tmpdir.ensure("jit_cache", dir=True)

    # Check that the wave model produces the same results in C++ and Python
    cppx, cppz = wave_model.velocity_cpp()
    cpp = (
        cpp_wrapper.replace("CODE_X_GOES_HERE", cppx)
        .replace("CODE_Z_GOES_HERE", cppz)
        .replace("x[2]", "x[1]")
    )
    mod = jit_compile(cpp, cache_dir)

    # The input values
    xr, zr = wave_locations(wave_model)
    t = 4.2

    # Compute the velocities using both the C++ and the Python versions
    vx_cpp = numpy.zeros_like(xr)
    vz_cpp = numpy.zeros_like(xr)
    for i in range(xr.size):
        vx_cpp[i] = mod.vel_x([xr[i], zr[i]], t)
        vz_cpp[i] = mod.vel_z([xr[i], zr[i]], t)
    vel_py = wave_model.velocity(xr, zr, t)

    # Check the results
    test_name_x = "%s C++ x-velocity" % wave_model.__class__.__name__
    test_name_z = "%s C++ z-velocity" % wave_model.__class__.__name__
    check_results(xr, zr, vel_py[:, 0], vx_cpp, test_name_x, 1e-5)
    check_results(xr, zr, vel_py[:, 1], vz_cpp, test_name_z, 1e-5)


@skipif_no_compile
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
    cache_dir = tmpdir.ensure("jit_cache", dir=True)

    # Check that the wave model produces the same results in C++ and Python
    cpp = wave_model.stream_function_cpp(frame="c")
    cpp = cpp_wrapper.replace("CODE_GOES_HERE", cpp).replace("x[2]", "x[1]")
    mod = jit_compile(cpp, cache_dir)

    # The input values
    xr, zr = wave_locations(wave_model)
    t = -23.9

    # Compute the stream function using both the C++ and the Python versions
    sf_cpp = numpy.zeros_like(xr)
    for i in range(xr.size):
        sf_cpp[i] = mod.sfunc([xr[i], zr[i]], t)
    sf_py = wave_model.stream_function(xr, zr, t, frame="c")

    # Check the results
    test_name = "%s C++ stream function" % wave_model.__class__.__name__
    check_results(xr, zr, sf_py, sf_cpp, test_name, 1e-3)


@skipif_no_compile
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
    cache_dir = tmpdir.ensure("jit_cache", dir=True)

    # Check that the wave model produces the same results in C++ and Python
    cpp = wave_model.slope_cpp()
    assert "x[2]" not in cpp
    cpp = cpp_wrapper.replace("CODE_GOES_HERE", cpp)
    mod = jit_compile(cpp, cache_dir)

    # The input values
    xr, zr = wave_locations(wave_model)
    t = 100.0

    # Compute the surface slope using both the C++ and the Python versions
    slope_cpp = numpy.zeros_like(xr)
    for i in range(xr.size):
        slope_cpp[i] = mod.slope([xr[i], zr[i]], t)
    slope_py = wave_model.surface_slope(xr, t)

    # Check the results
    test_name = "%s C++ surface slope" % wave_model.__class__.__name__
    check_results(xr, zr, slope_py, slope_cpp, test_name, 1e-16)


def test_cpp_elevation_nocompile(wave_model):
    # Check that the wave model produces valid C++ code
    cpp = wave_model.elevation_cpp()
    assert "x[2]" not in cpp
    if 'np.' in cpp or 'numpy' in cpp:
        print(cpp)
        assert False


def test_cpp_velocity_nocompile(wave_model):
    # Check that the wave model produces valid C++ code
    cpp = wave_model.velocity_cpp()
    if 'np.' in cpp or 'numpy' in cpp:
        print(cpp)
        assert False


def test_cpp_stream_function_nocompile(wave_model):
    # Check that the wave model produces valid C++ code
    cpp = wave_model.stream_function_cpp()
    if 'np.' in cpp or 'numpy' in cpp:
        print(cpp)
        assert False


def test_cpp_slope_nocompile(wave_model):
    # Check that the wave model produces valid C++ code
    cpp = wave_model.slope_cpp()
    if 'np.' in cpp or 'numpy' in cpp:
        print(cpp)
        assert False
