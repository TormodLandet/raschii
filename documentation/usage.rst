=========================
Using Raschii from Python
=========================

Most of the interaciton with Raschii will be through a WaveModel object. To get
such an object, first get the class and then instantiate the class to get the
wave model::

    import raschii
    WaveModel, AirModel = raschii.get_wave_model('Fenton')
    wave = WaveModel(height=12, depth=200, length=100, N=5)

You can ignore the air model class if you are just interested in the wave. To
show that the maximum elevation for this wave is 207.45 meters you can run::

    elev = wave.surface_elevation([0.0, 10.0, 20.0])
    print(elev)

You can get the crest velocity by running::

    vel = wave.velocity(0.0, elev[0])
    print(vel)

This will show that the crest velocity is approximately 7.6 m/s. Note that
asking for the velocity at ``(0.0, 207.45)`` will result in 0.0 since this point
is ever so slightly above the free surface. To get velocities above the free
surface you need to specify a method to compute the velocities in the air phase,
see :ref:`sec_blending` and the description of the air-phase models above that 
section to understand how Raschii handles this.

The code to compute velocities with an air-phase model is::

    WaveModel, AirModel = raschii.get_wave_model('Fenton', 'FentonAir')
    air = AirModel(height=100, blending_height=20)
    wave = WaveModel(height=12, depth=200, length=100, N=5, air=air)
    
    vel = wave.velocity(0.0, 208.0)
    print(vel)

This computes the velocities in the air above the crest. In this blended model
the velocities will increase slightly above the free surface before they reduce,
change direction, and then reduce to zero (in the vertical direction) at a
distance ``blending_height`` above the mean free surface. The ``height`` of the
air domain should be at least as large as the ``blending_height``.