=========================
Using Raschii from Python
=========================


Basic usage
===========

Most of the interaciton with Raschii will be through a WaveModel object. To get
such an object, first get the class and then instantiate the class to get the
wave model

.. code-block:: python

    import raschii

    WaveModel, AirModel = raschii.get_wave_model("Fenton")

    wave = WaveModel(height=12, depth=200, length=100, N=5)

You can ignore the air model class if you are just interested in the wave. To
show that the maximum elevation for this wave is 207.45 meters you can run

.. code-block:: python

    elev = wave.surface_elevation([0.0, 10.0, 20.0])
    print(elev)

You can get the crest velocity by running

.. code-block:: python

    vel = wave.velocity(0.0, elev[0])
    print(vel)

This will show that the crest velocity is approximately 7.6 m/s.

You can give period instead of wave length if you want (added in version 1.1.0).
This is slightly slower, but still relatively instant unless you are generating
a ton of waves. The period is used to iteratively try to find a wave length that
results in the correct period, which is easy for deep-water Airy waves and not
so easy for finite depth and higher-order waves.


Most common parameters
----------------------

Wave *height*:
  given in meters.
  This is the double-amplitude for linear waves.
  The crests are higher than the troughs are deep for non-linear waves,
  so the height is no-longer the double of the wave crest height.

Water *depth*:
 given in meters.

Wave *length*:
  given in meters.
  If you give *None* then the period is used instead (slower).

Wave *period*:
  given in seconds.
  *Only used when the wave length is None*.

Wave order, *N*:
  the order of the wave model to use.
  For Airy waves this is always 1, for Stokes waves it can be 1, 2, or 5,
  and for Fenton waves it can be any integer greater than or equal to 1.
  The higher the order, the more accurate the wave model is, but also the
  slower it is to compute (and for Fenton, the more it is likely to give
  non-sensical results for steep waves, such as finding an irregular
  two-peak solution).

See the API docs for :class:`~raschii.AiryWave`, :class:`~raschii.StokesWave`, and
:class:`~raschii.FentonWave` for more information.


Air-phase model
===============

Asking for the velocity above the free surface will result in zero. To get velocities
above the free surface you need to specify a method to compute the velocities in the
air phase, see :ref:`sec_blending` and the description of the air-phase models above
that  section to understand how Raschii handles this.

The code to compute velocities with an air-phase model is

.. code-block:: python

    import raschii

    WaveModel, AirModel = raschii.get_wave_model('Fenton', 'FentonAir')
    
    air = AirModel(height=100, blending_height=20)
    wave = WaveModel(height=12, depth=200, length=100, N=5, air=air)
    
    vel0 = wave.velocity(0.0, 207.0)  # Slightly below the free surface
    vel1 = wave.velocity(0.0, 208.0)  # Slightly above the free surface
    vel2 = wave.velocity(0.0, 220.0)  # Significantly above the free surface
    print(vel0, vel1, vel2)

This computes the velocities in the air above the crest. In this blended model
the velocities will increase slightly above the free surface before they reduce,
change direction, and then reduce to zero (in the vertical direction) at a
distance ``blending_height`` above the mean free surface. The ``height`` of the
air domain should be at least as large as the ``blending_height``.


SWD: Spectral Wave Data format
==============================

To write the wave elevation and kinematics to the SWD (Spectral Wave Data) file
format, e.g. for use as an incident wave field in a CFD or potential flow simulation,
use the `write_swd` method on the wave class

.. code-block:: python

    import raschii

    WaveModel, _AirModel = raschii.get_wave_model('Fenton')

    wave = WaveModel(height=12, depth=200, length=100, N=5)

    wave.write_swd("my_fenton_wave.swd", tmax=200.0, dt=0.01)

More infomation about SWD can
be found at the GitHub repo at https://github.com/SpectralWaveData/spectral_wave_data
and in the documentation at https://spectral-wave-data.readthedocs.io/ where the
underlying SWD wave description is also described. Raschii waves are stored as SWD
shape-class 2 (long-crested waves in constant water depth with constant spacing
:math:`\Delta k`)

The air model is not a part of the SWD file format and the kinematics above the free
surface are hence decided by the SWD library you use and how your program choses to
use the SWD data. Some versions of OpenFOAM will query the wave model to get the
elevation and only look up kinematics below the free surface, treating the air-phase
totally separate. Adapters for using SWD-files in OpenFOAM, Star CCM+, DNV Wasim and
other wave-simulation programs exist, but currently none that are open source as far
as we know. Writing a custom adapter is relatively straight forward since the SWD
library itself is open source. Interfacing with Raschii waves using the SWD file
format is a recommended way to integrate other programs with Raschii.


===================================
Using Raschii from the command line
===================================

You can generate SWD-files from the command line, run

.. code:: bash

    python -m raschii.cmd.swd --help
    
to see the options. You can also plot the waves with similar input parameters

.. code:: bash

    python -m raschii.cmd.plot Fenton 10 100 100 --velocities --ymin 50

Use ``--help`` to see all options. Two plot windows should pop up on your screen.
