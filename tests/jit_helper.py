import sys
import os
import subprocess
import importlib
import hashlib


# CMake code from pybind11 docs
CMAKE_TEMPLATE = """
cmake_minimum_required(VERSION 2.8.12)
project(jit_compile_cpp_module)
find_package(pybind11 REQUIRED)
pybind11_add_module(MODNAME jitmod.cpp)
"""


def jit_compile(cpp_code, cache_dir, verbose=True):
    """
    A simple JIT for C++ code using pybind11, included to test the generated C++
    code from the tests
    """
    # Hash the code object to make module name and check for cached version
    cpp_hash = hashlib.sha1(cpp_code.encode("utf-8", "replace")).hexdigest()

    # Create a unique Python module name
    assert "PYBIND11_MODULE(MODNAME, m)" in cpp_code
    modname = "jitmod_%s" % cpp_hash
    cpp_code = cpp_code.replace("PYBIND11_MODULE(MODNAME, m)", "PYBIND11_MODULE(%s, m)" % modname)

    # Ensure that the JIT directory exists
    dirname = os.path.join(cache_dir, "pybind11_jit_%s" % cpp_hash)
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

    mod = import_mod(dirname, modname)
    if mod is not None:
        # Module exists in cache
        return mod

    # Write necessary files to compile the module
    cpp_file = os.path.join(dirname, "jitmod.cpp")
    cmake_file = os.path.join(dirname, "CMakeLists.txt")
    with open(cpp_file, "wt") as f:
        f.write(cpp_code)
    with open(cmake_file, "wt") as f:
        f.write(CMAKE_TEMPLATE.replace("MODNAME", modname))

    # Compile the module
    if not compile_mod(dirname, verbose):
        raise ValueError("Could not compile module in %s" % dirname)

    # Import the compiled module
    mod = import_mod(dirname, modname, verbose)

    if mod is None:
        raise ValueError("Could not import the compiled module in %s" % dirname)
    else:
        return mod


def import_mod(dirname, modname, verbose=False):
    try:
        sys.path.insert(0, dirname)
        return importlib.import_module(modname)
    except ImportError as err:
        if verbose:
            print(err)
        return None
    finally:
        del sys.path[0]


def compile_mod(dirname, verbose):
    # Compile the module
    stdoutf = os.path.join(dirname, "stdout.txt")
    stderrf = os.path.join(dirname, "stderr.txt")
    shellscript = os.path.join(dirname, "runme.sh")
    with open(shellscript, "wt") as sh:
        sh.write("#!/bin/sh\n")
        sh.write("set -eu\n")
        sh.write("cmake .\n")
        sh.write("make\n")
    with open(stdoutf, "w") as out, open(stderrf, "w") as err:
        try:
            p = subprocess.Popen(["sh", shellscript], cwd=dirname, stdout=out, stderr=err)
            ok = p.wait() == 0
        except Exception as e:
            err.write("\n\nGOT Python EXCEPTION\n%r" % repr(e))
            ok = False

    if not ok:
        print("#############################################################")
        print("Compiler stdout:\n" + open(stdoutf, "rt").read())
        print("#############################################################")
        print("Compiler stderr:\n" + open(stderrf, "rt").read())
        print("#############################################################")
    return ok
