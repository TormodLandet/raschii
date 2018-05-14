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
