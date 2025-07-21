===========================
Online wave calculator demo
===========================

If you are reading this in a web browser then a demo should be shown below.
What you are then seeing is a Raschii running in the browser using
`PyScript <https://pyscript.net/>`_ and `Pyodide <https://pyodide.org/>`_.
Due to limitations in the floating point precision of JavaScript, the
wave models implemented in Raschii are not as accurate as you typically
get when running in a normal Python interpreter, but they are still
sufficiently accurate for most cases.

.. raw:: html
    
    <iframe src="./_static/raschii_pyscript.html" height="1000px" width="100%"></iframe>

Hint: click the plot of the generated wave to see the particle velocities!
It takes some time to download Python, numpy, and Raschii the first time
you run this demo, so please be patient.

This web visualization was made after being inspired by the `original online
wave calculators <http://www.coastal.udel.edu/faculty/rad/>`_ from the 1990s by 
Robert A. Dalrymple, specifically the `Dean stream function wave theory 
calculator <http://www.coastal.udel.edu/faculty/rad/streamless.html>`_.
Unfortunately Dalrymple's calculators are based on Java applet technology which
does not work well in modern web browsers (though, the JAR files can be 
downloaded and run locally by a sufficiently technologically proficient user).

There was once a port of Raschii to Dart (as an alternative to Python), but now
that PyScript (Pyodide) is available, it is no longer maintained.
