Raschii
=======

Raschii is a Python library for constructing non-linear regular waves and is
named after `Thysanoessa raschii
<https://en.wikipedia.org/wiki/Thysanoessa_raschii>`_, the Arctic Krill.

Supported wave models are currently:

- Stream function waves (M. M. Rienecker and J. D. Fenton, 1981)
- Stokes second through fifth order waves (based on J. D. Fenton, 1985) 
- Airy waves, standard linar waves cosine waves

Raschii includes a command line program to plot regular waves from the suported
wave models and C++ code generation for using the results in other programs, 
such as in `FEniCS <https://www.fenicsproject.org/>`_ expressions for initial
and boundary conditions in a FEM solver. There is also a limited `Dart port
<https://bitbucket.org/trlandet/raschiidart>`_ which is used in the `online demo
<https://raschii.readthedocs.io/en/latest/raschii_dart.html>`_.

.. figure:: http://raschii.readthedocs.io/en/latest/_static/fenton_stokes.png
   :alt: A comparison of Stokes and Fenton waves of fifth order

   A comparison of fifth order Stokes waves and fifth order Fenton stream
   function waves. Deep water, wave heigh 12 m, wave length 100 m.

As of version 1.0.3, Raschii can output waves on the SWD_ (spectral wave data)
standard file format for use in analysis programs (CFD, boundary element etc).

.. _SWD: https://github.com/SpectralWaveData/spectral_wave_data

Installation and running
------------------------

Raschii can be installed by running::

    python -m pip install raschii

Substitute ``python`` with ``python3`` as appropriate to your installation.
The command will also install any dependencies (numpy).

Using Raschii from Python
.........................

An example of using Raschii from Python::

    import raschii
    
    fwave = raschii.FentonWave(height=0.25, depth=0.5, length=2.0, N=20)
    print(fwave.surface_elevation(x=0))
    print(fwave.surface_elevation(x=[0, 0.1, 0.2, 0.3]))
    print(fwave.velocity(x=0, z=0.2))

This will output::

    [0.67352456]
    [0.67352456 0.61795882 0.57230232 0.53352878]
    [[0.27263788 0.        ]]

Using Raschii from the command line
...................................

You can also use Raschii from the command line. You can plot the wave
elevation and particle velocities, and also write swd files. See the 
help for the command line programs to get detailed usage info::

  python -m raschii.cmd.plot -h
  python -m raschii.cmd.swd -h

Substitute ``python`` with ``python3`` as appropriate to your installation.
You must install the ``matplotlib`` Python package to be able to use the
plot command.

An example of using Rascii from the command line::

  # Plot a 0.2 m high wave that is 2 meters long in 1.5 meters water depth
  # Some information abot the wave is also shown
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

Raschii is developed in Python on `Bitbucket <https://bitbucket.org/trlandet/raschii>`_
by use of the Git version control system. If you are reading this on github,
please be aware that you are seeing a mirror that could potentially be months
out of date. All pull requests and issues should go to the Bitbucket repository.

Raschii is automatically tested on `CircleCI <https://circleci.com/bb/trlandet/raschii/tree/master>`_  
and the current CI build status is |circleci_status|.

.. |circleci_status| image:: https://circleci.com/bb/trlandet/raschii.svg?style=svg&circle-token=d0d6c55654d1c7ba49a9679d7dd1623e1b52b748
  :target: https://circleci.com/bb/trlandet/raschii/tree/master


Releases
--------

Version 1.0.6 - June 28. 2024
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

- Fix missing time dependency in Stokes surface surface elevation
- Ensure all wave models implement ``T`` and ``omega`` attributes
- Test that the surface elevation has the correct period for all wave models
- Include `SWD <https://github.com/SpectralWaveData/spectral_wave_data>`_ file 
  format support for writing generated waves to files for interchange with other
  tools.

Version 1.0.2 - Jun 4. 2018
............................

Some more work on air-phase / water phase velocity blending 

- Change the air blending zone to be horizontal at the top (still follows the
  wave profile at the bottom). The air phase blanding still has no influenece on
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

Raschii is copyright Tormod Landet, 2018. Raschii is licensed under the Apache
2.0 license, a  permissive free software license compatible with version 3 of
the GNU GPL. See the file ``LICENSE`` for the details.
