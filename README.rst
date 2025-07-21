Raschii
=======

Raschii is a Python library for constructing non-linear regular waves.
Supported wave models are:

- Stream function waves (M. M. Rienecker and J. D. Fenton, 1981)
- Stokes second- through fifth-order waves (based on J. D. Fenton, 1985) 
- Airy waves, the standard linear wave theory

Raschii directly computes the series expansion coefficients (for Stokes waves),
or solves a minimization problem for the coefficients (for stream function
waves). Once the expansion coefficients of the wave has been found you can

1) Compute the free-surface elevation
2) Compute the particle velocities at any point in the wave field (and in the
   air phase by use of the potential-flow solution above the free-surface).
3) Export the wave field to a file in the SWD_ file format (Spectral Wave Data, see
   more below). The SWD library can then be used to compute many other properties
   such as the velocity potential, the particle accelerations, the stream function,
   the Bernoulli pressure, `and several other properties
   <https://spectral-wave-data.readthedocs.io/en/latest/theory.html#spectral-kinematics>`_.

.. figure:: http://raschii.readthedocs.io/en/latest/_static/fenton_stokes.png
   :alt: A comparison of Stokes and Fenton waves of fifth order

   A comparison of fifth-order Stokes waves and fifth-order Fenton stream
   function waves. Deep water, wave height 12 m, wave length 100 m.


Using the Raschii waves
-----------------------

Raschii can output waves in the SWD_ (spectral wave data) standard file format for
use as the incoming incident waves in flow analysis programs such as potential-flow
or CFD (Euler and Navier-Stokes equations etc). The SWD file format is used by other
linear and non-linear wave-generation tools, such as the non-linear irregular-wave
HOSM solver DNV WAMOD. By supporting the SWD file format you avoid having to implement
your own wave models for all types of waves, and you can instead use a uniform wave
API for all types of linear and non-linear waves and let the SWD library handle the
details of the wave kinematics.

Raschii also includes a command line program to plot regular waves from the supported
wave models and C++ code generation for using the results in other programs, such as
in `FEniCS <https://www.fenicsproject.org/>`_ expressions for initial and boundary
conditions in a FEM solver (if you do not want to link to the SWD library for some
reason).

.. _SWD: https://github.com/SpectralWaveData/spectral_wave_data


Installation and running
------------------------

Raschii can be installed using pip, uv, or any other Python package manager that
supports the `Python Package Index (PyPI) <https://pypi.org/project/raschii>`_.
Example:

.. code:: bash

    pip install raschii


Using Raschii from Python
.........................

An example of using Raschii from Python:

.. code:: python

    import raschii
    
    fwave = raschii.FentonWave(height=0.25, depth=0.5, length=2.0, N=20)

    print(fwave.surface_elevation(x=0))
    print(fwave.surface_elevation(x=[0, 0.1, 0.2, 0.3]))
    print(fwave.velocity(x=0, z=0.2))

This will output:

.. code:: output

    [0.67352456]
    [0.67352456 0.61795882 0.57230232 0.53352878]
    [[0.27263788 0.        ]]

See the `documentation <https://raschii.readthedocs.io/en/latest/usage.html>`_ for more
information on the available parameters, methods and attributes of the wave classes.


Using Raschii from the command line
...................................

You can also use Raschii from the command line. You can plot the wave
elevation and particle velocities, and also write SWD files. See the 
help for the command line programs to get detailed usage info.

.. code:: bash

  python -m raschii.cmd.plot -h
  python -m raschii.cmd.swd -h

Substitute ``python`` with ``python3`` or ``uv run``  as appropriate to your installation.
You must install the ``matplotlib`` Python package to be able to use the plot command.

An example of using Raschii from the command line:

.. code:: bash

  # Plot a 0.2 m high wave that is 2 meters long in 1.5 meters water depth
  # Some information about the wave is also shown
  python -m raschii.cmd.plot -N 5 Fenton 0.2 1.5 2

  # Save the same stream function wave to a SWD file
  python -m raschii.cmd.swd -N 5 fenton.swd Fenton 0.2 1.5 2  

The plot tool allows comparing multiple waves, the SWD file writer only
supports a single wave at a time and does currently not support Airy waves.


Documentation
-------------

.. TOC_STARTS_HERE  - in the Sphinx documentation a table of contents will be inserted here 

The documentation can be found on `Raschii's Read-the-Docs pages
<https://raschii.readthedocs.io/en/latest/index.html#documentation>`_.

.. TOC_ENDS_HERE


Development
-----------

Raschii is developed in Python on `GitHub <https://github.com/TormodLandet/raschii>`_
using the Git version control system.

Raschii is automatically tested using pytest and GitHub Actions and the current CI build status is
|ci_status|.

.. |ci_status| image:: https://github.com/TormodLandet/raschii/actions/workflows/pytest.yml/badge.svg
  :target: https://github.com/TormodLandet/raschii/actions/workflows/pytest.yml


Releases
--------

Version 1.1.1 - Unreleased
............................

- Use PyScript to run Raschii in the browser (instead of the incomplete Dart port).
- Increase the speed of the SWD file writer

Version 1.1.0 - Jun 18. 2025
.............................

- Support for giving the wave period instead of the wave length
- Support for infinite depth waves. This is not fully complete, but should be
  sufficient to export proper SWD files for infinite depth waves.
  Set depth=-1 to use infinite depth waves.
- Better testing of the SWD file exporter when the SpectralWaveData package is not installed
  by including a simplified SWD file reader for the tests.
- Move repository and CI to GitHub. Tested on Python 3.10 (Ubuntu 22.04), and Python 3.12 (uv).

Version 1.0.7 - Sep 30. 2024
.............................

- Support for numpy 2.1
- Drop support for Python 3.9 and older (`following numpy <https://numpy.org/neps/nep-0029-deprecation_policy.html>`_)
- Added testing with latest Python available via uv (currently CPython 3.12)

Version 1.0.6 - Jun 28. 2024
.............................

- Support for numpy 2.0
- Add type annotations
- Add API docs for public API functions

Version 1.0.5 - Jan 25. 2024
............................

- Update the documentation
- Unbreak the read-the-docs builder
- Switch to pyproject.toml from setup.py (replace setuptools with hatchling)
- No new code or functionality added or removed, just housekeeping

Version 1.0.4 - Aug 28. 2020
............................

- Add the ``raschii.cmd.plot`` and ``raschii.cmd.swd`` command line programs

Version 1.0.3 - Aug 28. 2020
............................

- Fix missing time dependency in Stokes surface elevation
- Ensure all wave models implement ``T`` and ``omega`` attributes
- Test that the surface elevation has the correct period for all wave models
- Include `SWD <https://github.com/SpectralWaveData/spectral_wave_data>`_ file 
  format support for writing generated waves to files for interchange with other
  tools.

Version 1.0.2 - Jun 4. 2018
............................

Some more work on air-phase / water phase velocity blending 

- Change the air blending zone to be horizontal at the top (still follows the
  wave profile at the bottom). The air phase blending still has no influence on
  the wave profile or water-phase velocities, but the transition from blended to
  pure air-phase velocities is now a bit smoother for steep waves and the 
  divergence of the resulting field is lower when projected into a FEM function
  space (analytically the divergence is always zero).  

Version 1.0.1 - May 31. 2018
............................

Small bugfix release

- Fix bug related to sign of x component of FentonAir C++ velocity
- Improve unit testing suite
- Improve FEM interpolation demo

Version 1.0.0 - May 29. 2018
............................

The initial release of Raschii

- Support for Fenton stream functions (Rienecker and Fenton, 1981)
- Support for Stokes 1st - 5th order waves (Fenton, 1985)
- Support for Airy waves
- Support for C++ code generation (for FEniCS expressions etc)
- Command line program for plotting waves
- Command line demo for converting fields to FEniCS
- Unit tests for most things
- Documentation and (currently non-complete online demo)
- Support for computing a combined wave and air velocity field which is
  continuous across the free surface and divergence free (currently only works
  with the Fenton stream function wave model).


Copyright and license
---------------------

Raschii is copyright Tormod Landet (2018--).

Raschii is licensed under the Apache 2.0 license,
a permissive free software license compatible with version 3 of the GNU GPL.
See the file ``LICENSE`` for the details.
