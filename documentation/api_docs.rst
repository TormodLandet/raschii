=========================
Raschii Python API
=========================

Documentation of the Raschii Python API, automatically generated from the source-code comments.

Main functions
==============

.. autofunction:: raschii.check_breaking_criteria

.. autofunction:: raschii.get_wave_model

.. autodata:: raschii.WAVE_MODELS

.. autodata:: raschii.AIR_MODELS

.. autodata:: raschii.__version__


Wave model classes
==================

.. autoclass:: raschii.AiryWave
    :members:

.. autoclass:: raschii.StokesWave
    :members:

.. autoclass:: raschii.FentonWave
    :members:


Air model classes
=================

.. autoclass:: raschii.FentonAirPhase
    :members:

.. autoclass:: raschii.ConstantAirPhase
    :members:


Exceptions
=================

.. autoexception:: raschii.RasciiError

.. autoexception:: raschii.NonConvergenceError
