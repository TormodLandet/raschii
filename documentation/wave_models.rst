===========
Wave models
===========

All wave models are implemented in such a way that there is a crest at ``x=0``
at ``t = 0`` and the waves are traveling in the positive x-direction. The origin
of the coordinate system is on the sea floor, which is assumed to be horizontal.
The default acceleration of gravity is ``9.81``, so the input quantities should
be given in the SI system.


Airy
====

This is the simplest wave theory, a single cosine wave. It is included mostly
for testing the other models and is only valid in the limit of very low waves
in very deep water. 


Stokes
======

Stokes waves of 1st to 5th order are implemented following John D. Fenton's 1985
paper *A Fifth‐Order Stokes Theory for Steady Waves*. By increasing the 
requested order ``N``, the code will include more and more expansion
coefficients, starting from linear Airy waves at 1st order and properly
replicating 2nd, 3rd, 4th and finally 5th order Stokes waves. Any higher order
wave requested will issue a warning and return the fift order solution.

Stokes waves are good approximation in the deep water limit. No time consuming 
calculations are required to generate these waves and the computations will not
diverge. Both of these issues can be problematic for stream function waves that
need to optimize a non-linear function of many parameters.

Further details and analytical expression for all the coefficients in the 
perturbation expansion can be found in the original paper, available on `John D.
Fentons web pages <http://johndfenton.com/Papers/Papers-John%20Fenton.html>`_.


Fenton stream function waves
=============================

High order regular wave theories based on truncated Fourier series approximating
the stream function was pioneered by Dean (1965). Our implementation is based on
Rienecker and Fenton's 1981 paper, *A Fourier approximation method for steady
water waves*, which is often refered to as "Fenton" stream function wave theory
to differentiate it from the original "Dean" stream function wave theory.

The method is based on collocation (solving the non-linear equations exactly in
N + 1 points) and is based on Newton–Raphson iterations to tackle the
non-linearities. The unknowns are the expansion coefficients :math:`B_j`, the
wave elevation :math:`\eta(x_m)`, the stream functions value at the free surface
:math:`Q` and the Bernoulli constant at the free surface :math:`R`.

The stream function a-priori satisfies the bottom boundary condition at z=0 and
also the Laplace equation :math:`\nabla^2\Psi=0`. It is defined as

.. math::

    \Psi = B_0 z + \sum_{j=1}^{N} B_j \frac{\sinh jkz}{\cosh jkD} \cos(j k x),

which is non-linear in :math:`\eta` on the free surface where `z=\eta`. To find
the unknowns the following conditions are requested to be met:

- The free surface is a stream line, such that
  :math:`\Psi(x, \eta) = -Q`.
- The pressure is constant at the free surface
  (Bernoulli constant :math:`R`).
- The wave height is :math:`H`, such that
  :math:`\eta(0) - \eta(\lambda/2) = H`
- The mean wave height is :math:`D`, such that
  :math:`\int_0^{\lambda/2} \eta\,\mathrm d x = D \lambda / 2`. 

Further details and analytical expression for all terms of the Jacobian matrix
used in the Newton-iterations can be found in the original paper, available on
`John D. Fentons web pages
<http://johndfenton.com/Papers/Papers-John%20Fenton.html>`_.


Air phase
=========

Raschii was origially written to provide good initial conditions for a two-phase
Navier-Stokes solver. In order to initialise the domain with a divergence free
velocity field it is important to also compute the velocities in the air phase.
The Fourier series stream function from Rienecker and Fenton (1981) is used also
for the air phase. Using a stream function ensures an exactly divergence free
velocity field.

The "Fenton" stream function is linear in the unknown parameters when the wave
elevation is known. The air-phase velocities are computed after the water wave
has been established, so this means that the expansion coefficients
:math:`B_{1..N}` can be found by a simple linear solve to satisfy that the free
surface is a stream function also in the air phase. In order to use the "Fenton"
stream function the z-coordinate is flipped such that the velocity is purely
horizontal a distance ``depth_air`` above the free surface.

.. figure:: figures/air_vel_compare.png
   :alt: comparison of blended and unblended stream function velocities near the
         free surface

   Comparison of blended and unblended stream function velocities near the free
   surface

The resulting velocity field can be seen on the left in the above figure. As can
be seen the velocities in the water and air phases are reasonably continuous,
but some discrepancies are highlighted. To combat these problems a blending is
done between the two stream functions. The stream function in the water domain
is left entirely undisturbed, but from the free surface and a distance 
``air_blend_distance`` up—by default the same as the wave height—the stream 
function in the water is smoothly transitioned into the the stream function in
the air. The results can be seen in the right image. Since the blending is done
to create a new stream function, :math:`\Psi_\text{blended} = (1 - f)
\Psi_\text{water} + f \Psi_\text{air}`, the resulting velocity field is exactly
divergence free.    

 