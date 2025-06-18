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

  usage: raschii.cmd.swd [-h] [-N N] [--dt DT] [--tmax TMAX] [-f]
                         swd_file     wave_type
                         wave_height  water_depth  wave_length

  Write a Raschii wave to file (SWD format)

  positional arguments:
    swd_file     Name of the SWD file to write.
    wave_type    Name of the wave model.
    wave_height  Wave height
    water_depth  The still water depth
    wave_length  Distance between peaks

  options:
    -h, --help   show this help message and exit
    -N N         Approximation order
    --dt DT      Timestep
    --tmax TMAX  Duration
    -f, --force  Allow exceeding breaking criteria


Basic plots
===========

Raschii can create basic plots of the wave elevation and kinematics.
This only works if you separately install the `matplotlib` package which is not a required
dependency of Raschii since it is not needed for the main functionality.

To plot a wave, run

.. code:: bash

    python -m raschii.cmd.plot Fenton 10 100 100 --velocities --ymin 50

Use ``--help`` to see all options. Two plot windows should pop up on your screen.
