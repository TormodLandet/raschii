.. _raschii_pyscript:

===========================
Online demo
===========================

If you are reading this in a web browser then a demo should be shown below.
What you are then seeing is a Raschii running in the browser using
`PyScript <https://pyscript.net/>`_ and `Pyodide <https://pyodide.org/>`_.


Visualization
=============

You can click the figure with the wave profile to see the particle velocities.
It takes some time to download Python, numpy, and Raschii the first time you
run this demo, so please be patient.

.. raw:: html
    
    <iframe src="./_static/raschii_pyscript.html" height="1000px" width="100%"></iframe>

This web visualization was made after being inspired by the `original online
wave calculators <http://www.coastal.udel.edu/faculty/rad/>`_ from the 1990s by 
Robert A. Dalrymple, specifically the `Dean stream function wave theory 
calculator <http://www.coastal.udel.edu/faculty/rad/streamless.html>`_.
Unfortunately Dalrymple's calculators are based on Java applet technology which
does not work well in modern web browsers (though, the JAR files can be 
downloaded and run locally by a sufficiently technologically proficient user).

There was once a `port of Raschii to Dart <https://bitbucket.org/trlandet/raschiidart/>`_
(as an alternative to Python), but now that PyScript (Pyodide) is available, it is no
longer maintained.


Online code editor
==================

You can also write your own code using Raschii (and numpy) in the below editor.
Press the run (play) button to run the code.

.. raw:: html
    
    <iframe src="./_static/raschii_pyscript_editor.html" height="400px" width="100%"></iframe>
