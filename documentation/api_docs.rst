=========================
Raschii Python API
=========================

Documentation of the Raschii Python API, automatically generated from the source code comments.

.. contents::
   :local:


Main functions
==============

.. autofunction:: raschii.get_wave_model

.. autofunction:: raschii.check_breaking_criteria

.. .. autodata:: raschii.WAVE_MODELS

.. .. autodata:: raschii.AIR_MODELS

.. .. autodata:: raschii.__version__


Wave model API summary
======================

All three wave model classes (:class:`~raschii.AiryWave`,
:class:`~raschii.StokesWave`, and :class:`~raschii.FentonWave`) share the
following attributes and methods.

.. list-table:: Attributes and methods available on all wave models
   :widths: 45 55
   :header-rows: 1

   * - Attribute / Method
     - Description
   * - ``height``, ``depth``, ``length``, ``period``, ``g``
     - Constructor input parameters, stored as instance attributes.
   * - ``omega``
     - Angular frequency [rad/s] (attribute).
   * - ``k``
     - Wave number [1/m]  (attribute).
   * - ``c``
     - Phase speed [m/s]  (attribute).
   * - ``warnings``
     - String with any warnings raised during construction (empty if none, attribute).
   * - ``surface_elevation(x, t=0, include_depth=True)``
     - Free-surface elevation.  Returns a scalar when both *x* and *t* are
       scalar, an ndarray otherwise (NumPy scalar-in / scalar-out convention).
   * - ``velocity(x, z, t=0, all_points_wet=False)``
     - Fluid velocity at each *(x, z)* point.  Shape ``(2,)`` for scalar inputs,
       ``(N, 2)`` or ``(T, N, 2)`` for array inputs.
   * - ``velocity_potential(x, z, t=0)``
     - Earth-frame velocity potential φ (mean flow excluded).  Same scalar/array
       convention as ``surface_elevation``.
   * - ``write_swd(path, dt, tmax=None, nperiods=None, amp=1)``
     - Write the wave field to a file in the
       `SWD format <https://spectral-wave-data.readthedocs.io/>`_.
   * - ``cpp``
     - C++ code generator (experimental — see :ref:`cpp-code-generation`).

The following methods are only available on specific wave model classes.

.. list-table:: Methods available on specific wave models only
   :widths: 38 42 20
   :header-rows: 1

   * - Method
     - Description
     - Available on
   * - ``stream_function(x, z, t=0, frame=Frame.EARTH)``
     - Stream function ψ.  Use :class:`~raschii.Frame` to choose the reference
       frame (``Frame.EARTH`` or ``Frame.WAVE``).
     - :class:`~raschii.FentonWave`
   * - ``surface_slope(x, t=0)``
     - Horizontal derivative dη/dx of the free-surface elevation.
     - :class:`~raschii.FentonWave`
   * - ``acceleration(x, z, t=0, all_points_wet=False)``
     - Fluid acceleration.  Water-phase only; returns zero above the free
       surface (air blending is not yet implemented for accelerations).
     - :class:`~raschii.FentonWave`


Wave model classes
==================


Airy waves
----------

Raschii linear waves, see :ref:`the Airy wave model <theory_airy>`.

.. autoclass:: raschii.AiryWave
    :class-doc-from: init
    :members:
    :inherited-members: raschii.WaveModel


Stokes waves
------------

Raschii implements the Stokes 1st through 5th order wave models, see :ref:`the Stokes wave model <theory_stokes>`.

.. autoclass:: raschii.StokesWave
    :class-doc-from: init
    :members:
    :inherited-members: raschii.WaveModel


Fenton stream-function waves
----------------------------

Raschii implements the Fenton stream-function wave model as described in :ref:`the Fenton wave model <theory_fenton>`.

.. autoclass:: raschii.FentonWave
    :class-doc-from: init
    :members:
    :inherited-members: raschii.WaveModel


Air model classes
=================

Raschii implements special support for kinematics above the free surface,
see :ref:`theory_air-phase_models` for details.
You can use these to construct a fully divergence-free velocity field for a computational domain
with both water and air phases.
This is normally not done in lower-order methods such as the typical finite-volume solvers
(OpenFOAM etc.), but has been used in a higher-order fully divergence-free DG-FEM solver to
construct consistent initial and boundary conditions.

.. autoclass:: raschii.FentonAirPhase
    :class-doc-from: init
    :members:
    :inherited-members: raschii.AirPhaseModel

.. autoclass:: raschii.ConstantAirPhase
    :class-doc-from: init
    :members:
    :inherited-members: raschii.AirPhaseModel


Exceptions
==========

.. autoexception:: raschii.RaschiiError

.. autoexception:: raschii.NonConvergenceError


Enumerations
============

.. autoclass:: raschii.Frame
    :members:


Other
=====


.. _cpp-code-generation:

C++ code generators
-------------------

Each wave model exposes a ``wave.cpp`` attribute that can generate C++ code strings for use in e.g.
FEniCS boundary conditions.  The available methods are

* ``wave.cpp.elevation()``
* ``wave.cpp.velocity(all_points_wet=False)``
* ``wave.cpp.stream_function(frame=Frame.EARTH)``
* ``wave.cpp.slope()``

(the last two are implemented only for :class:`~raschii.FentonWave`).

The air-phase model classes expose a ``air.cpp`` attribute that can generate C++ code strings for
the air-phase velocity field.  The available methods are

* ``air.cpp.velocity()``
* ``air.cpp.stream_function(frame=Frame.EARTH)``

which are available in both air-phase models.

.. warning::

   The C++ code generation interfaces (``wave.cpp.*`` and ``air.cpp.*``) are **experimental**.
   The code-generation API may change or be removed in a future release without following the normal
   deprecation cycle. They are currently not used as far as the author of Raschii is aware.
