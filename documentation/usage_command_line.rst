.. _command_line_interface:

======================
Command line interface
======================

Below you will find a short introduction to using Raschii from the command line.
The main way of using Raschii is as a library in your Python code, see
the :ref:`user_manual` for more information on that.

.. contents::
  :local:


Creating SWD files
==================

Raschii comes with a command line interface for generating SWD files from the non-linear wave models.
To learn how to generate SWD files from the command line, run

.. code:: bash

    python -m raschii.cmd.swd --help
    
to see the options. At the time of writing, the output is

.. code:: text

  usage: raschii.cmd.swd [-h] [-T PERIOD] [-N N] [--dt DT] [--tmax TMAX]
                         [--swd-amp {1,2,3}] [-f]
                         swd_file wave_type wave_height water_depth [wave_length]

  Write a Raschii wave to file (SWD format)

  positional arguments:
    swd_file              Name of the SWD file to write.
    wave_type             Name of the wave model.
    wave_height           Wave height
    water_depth           The still water depth
    wave_length           Distance between wave crests. Mutually exclusive with
                          --period / -T.

  options:
    -h, --help            show this help message and exit
    -T PERIOD, --period PERIOD
                          Wave period in seconds (alternative to the positional
                          wave_length argument).
    -N N                  Approximation order
    --dt DT               Timestep
    --tmax TMAX           Duration
    --swd-amp {1,2,3}     SWD amp flag.
                          1 (default): store potential at z=0 (calm surface);
                          2: store potential on the wavy free surface;
                          3: store elevation only (no potential).
    -f, --force           Allow exceeding breaking criteria

Examples::

    # Write by wave length
    python -m raschii.cmd.swd -N 5 fenton.swd Fenton 0.2 1.5 2

    # Write by wave period
    python -m raschii.cmd.swd -N 5 --period 1.2 fenton.swd Fenton 0.2 1.5


Basic plots
===========

Raschii can create basic plots of the wave elevation and kinematics.
This only works if you separately install the `matplotlib` package which is not a required
dependency of Raschii since it is not needed for the main functionality.

To plot a wave by wave length, run

.. code:: bash

    python -m raschii.cmd.plot Fenton 10 100 100 --velocities --ymin 50

To plot a wave by wave period instead, use ``--period`` / ``-T``::

    python -m raschii.cmd.plot -T 8 Fenton 10 100 --velocities

Use ``--help`` to see all options. Two plot windows should pop up on your screen.
