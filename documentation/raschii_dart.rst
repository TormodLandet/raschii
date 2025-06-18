===========================
Online wave calculator demo
===========================


If you are reading this in a web browser then a demo should be shown below.
Please note that the Python to Javascript transpilers are currently not advanced
enough to convert the Raschii Python code to a format runnable on the web, so 
what you are seeing is a `port of Raschii to the web-friendly Dart language 
<https://bitbucket.org/trlandet/raschiidart>`_.

The web ecosystem is currently missing powerful libraries like NumPy, so parts
of the necessary linear algebra routines had to be re-implemented in Dart. For
this reason the online calculator may not give exactly the same results as the
Python version. The Python version is continuously tested with unit tests
and every change is compared against known solutions, so that is the most
trusted version should there be any discrepancies.


.. raw:: html
    :file: raschii_dart.html
    :encoding: utf-8

Hint: click the plot of the generated wave to see the particle velocities! The
Airy and Stokes wave models are fast, the Fenton stream function wave model may
take some time to generate a wave. When the wave is generated it is very fast to
compute particle velocities for all the implemented wave models.

This web visualization was made after being inspired by the `original online
wave calculators <http://www.coastal.udel.edu/faculty/rad/>`_ from the 1990s by 
Robert A. Dalrymple, specifically the `Dean stream function wave theory 
calculator <http://www.coastal.udel.edu/faculty/rad/streamless.html>`_.
Unfortunately Dalrymple's calculators are based on Java applet technology which
does not work well in modern web browsers (though, the JAR files can be 
downloaded and run locally by a sufficiently technologically proficient user).