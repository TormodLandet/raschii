=========================
Raschii Python API
=========================

Documentation of the Raschii Python API, automatically generated from the source-code comments.

.. contents::
   :local:


Main functions
==============

.. autofunction:: raschii.get_wave_model

.. autofunction:: raschii.check_breaking_criteria

.. .. autodata:: raschii.WAVE_MODELS

.. .. autodata:: raschii.AIR_MODELS

.. .. autodata:: raschii.__version__


Wave model classes
==================


Airy waves
----------

Raschii linear waves, see :ref:`the Airy wave model <theory_airy>`.

.. autoclass:: raschii.AiryWave
    :class-doc-from: init
    :members:


Stokes waves
------------

Raschii implements the Stokes 1st, 2nd, and 5th order wave models, see :ref:`the Stokes wave model <theory_stokes>`.

.. autoclass:: raschii.StokesWave
    :class-doc-from: init
    :members:


Fenton stream-function waves
----------------------------

Raschii implements the Fenton stream-function wave model as described in :ref:`the Fenton wave model <theory_fenton>`.

.. autoclass:: raschii.FentonWave
    :class-doc-from: init
    :members:


Air model classes
=================

Raschii implements special support for kinematics above the free surface,
see :ref:`theory_air-phase_models` for details.
You can use these to construct a fully divergence free velocity field for a computational domain
with both water and air phases.
This is normally not done in lower-order methods such as the typicall finite-volume solvers
(OpenFOAM etc.), but has been used in a higher-order fully divergence free DG-FEM solver to
construct consistent initial and boundary conditions.

.. autoclass:: raschii.FentonAirPhase
    :class-doc-from: init
    :members:


.. autoclass:: raschii.ConstantAirPhase
    :class-doc-from: init
    :members:


Exceptions
=================

.. autoexception:: raschii.RasciiError

.. autoexception:: raschii.NonConvergenceError
