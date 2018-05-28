from numpy import asarray, zeros


class ConstantAirPhase:
    def __init__(self, height, blending_height=None):
        """
        Constant horizontal velocity equal to the phase speed
        """
        self.height = height
        self.blending_height = blending_height
    
    def set_wave(self, wave):
        """
        Connect this air phase with the wave in the water phase
        """
        self.c = wave.c
        
        if self.blending_height is None:
            self.blending_height = 2 * wave.height
    
    def stream_function(self, x, z, t=0, frame='b'):
        """
        Compute the stream function at time t for position(s) x
        """
        if isinstance(x, (float, int)):
            x, z = [x], [z]
        z = asarray(z, dtype=float)
        
        if frame == 'e':
            return self.c * z
        elif frame == 'c':
            return 0.0 * z
    
    def velocity(self, x, z, t=0):
        """
        Compute the air phase particle velocity at time t for position(s) (x, z)
        where z is 0 at the bottom and equal to depth_water at the free surface
        and equal to depth_water + depth air at the top free slip lid above the
        air phase
        """
        if isinstance(x, (float, int)):
            x, z = [x], [z]
        x = asarray(x, dtype=float)
        z = asarray(z, dtype=float)
        
        return zeros((x.size, 2), float)
    
    def stream_function_cpp(self, frame='b'):
        """
        Return C++ code for evaluating the stream function of this specific
        wave. The positive traveling direction is x[0] and the vertical
        coordinate is x[2] which is zero at the bottom and equal to +depth at
        the mean water level.
        """
        if frame == 'b':
            return '%r * x[2]' % self.c
        elif frame == 'c':
            return '0.0'
    
    def velocity_cpp(self):
        """
        Return C++ code for evaluating the particle velocities of this specific
        wave. Returns the x and z components only with z positive upwards. The
        positive traveling direction is x[0] and the vertical coordinate is x[2]
        which is zero at the bottom and equal to +depth at the mean water level.
        """
        cpp_x = '0.0'
        cpp_z = '0.0'
        return (cpp_x, cpp_z)
    
    def __repr__(self):
        return ('ConstantAirPhase(height={s.height}, blending_height='
                '{s.blending_height})').format(s=self)
