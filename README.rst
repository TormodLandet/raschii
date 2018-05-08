Raschii
=======

Raschii is a Python implementation of stream function waves theory and is named
after `Thysanoessa raschii <https://en.wikipedia.org/wiki/Thysanoessa_raschii>`_,
the Arctic Krill.


Installation and running
------------------------

Raschii can be installed by running::

    python3 -m pip install raschii
    
This will also install dependencies like numpy.

An example of using Raschii from Python::

    import raschii
    
    fwave = raschii.FentonWave(height=0.25, depth=0.5, length=2.0, N=20)
    print(fwave.surface_elevation(x=0))
    print(fwave.surface_elevation(x=[0, 0.1, 0.2, 0.3]))
    print(fwave.velocity(x=0, z=0.2))


Development
-----------

Raschii is developed in Python on `Bitbucket <https://bitbucket.org/trlandet/raschii>`_
by use of the Git version control system. If you are reading this on github,
please be aware that you are seeing a mirror that could potentially be months
out of date. All pull requests and issues should go to the Bitbucket repository.

Raschii is automatically tested on `CircleCI <https://circleci.com/bb/trlandet/raschii/tree/master>`_  
and the current CI build status is |circleci_status|.

.. |circleci_status| image:: https://circleci.com/bb/trlandet/raschii.svg?style=svg&circle-token=d0d6c55654d1c7ba49a9679d7dd1623e1b52b748
  :target: https://circleci.com/bb/trlandet/raschii/tree/master


Releases
--------

Version 1.0.0 - May XX. 2018
............................

The initial release of Raschii

- Support for Fenton stream functions, Rienecker and Fenton (1981)


Copyright and license
---------------------

Raschii is copyright Tormod Landet, 2018. Raschii is licensed under the Apache
2.0 license, a  permissive free software license compatible with version 3 of
the GNU GPL. See the file ``LICENSE`` for the details.
