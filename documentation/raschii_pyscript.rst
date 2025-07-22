.. _raschii_pyscript:

===========================
Online demo
===========================

If you are reading this in a web browser then a demo should be shown below.
This is Raschii running in your browser using `PyScript <https://pyscript.net/>`_
and `Pyodide <https://pyodide.org/>`_.
Note that some numerical calculations may give slightly different results than in
"normal" Python running on your computer.

.. contents::
   :local:


Wave visualization
==================

You can click the figure with the wave profile to see the particle velocities
after the wave profile has been calculated.
You can also `open in a full window <./_static/raschii_pyscript.html>`_.

.. raw:: html
    
    <iframe src="./_static/raschii_pyscript.html" height="1000px" width="100%"></iframe>

This web visualization was made after being inspired by the `original online
wave calculators <http://www.coastal.udel.edu/faculty/rad/>`_ from the 1990s by 
Robert A. Dalrymple, specifically the `Dean stream function wave theory 
calculator <http://www.coastal.udel.edu/faculty/rad/streamless.html>`_.
Unfortunately, Dalrymple's calculators are based on Java applet technology which
does not work well in modern web browsers (though the JAR files can be 
downloaded and run locally by a sufficiently technologically proficient user).

There was once a `port of Raschii to Dart <https://bitbucket.org/trlandet/raschiidart/>`_
(as an alternative to Python), but now that PyScript (Pyodide) is available, it is no
longer maintained.


Online code editor
==================

You can write and run your own code using Raschii, numpy, and matplotlib in the online code editor.
Press the run (play) button to run the code; this button may not appear until you start editing.
It will take a while to run the code for the first time, but after that it is relatively fast.

`Open the online code editor <./_static/raschii_pyscript_editor.html>`_!
