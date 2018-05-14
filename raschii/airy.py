from numpy import pi, cos, sin, zeros, array, asarray


class AiryWave:
    required_input = ('height', 'depth', 'length')
    
    def __init__(self, height, depth, length, g=9.81):
        """
        Linear Airy waves
        
        * height: wave height above still water level
        * depth: still water distance from the flat sea bottom to the surface
        * length: the periodic length of the wave (distance between peaks)
        """
        self.height = height
        self.depth = depth
        self.length = length
        self.g = g
    
    def surface_elevation(self, x, t=0):
        """
        Compute the surface elavation at time t for position(s) x
        """
        if isinstance(x, (float, int)):
            x = array([x], float)
        x = asarray(x)
        
        return self.depth + self.height / 2 * cos(2 * pi * x / self.length)
    
    def velocity(self, x, z, t=0):
        """
        Compute the fluid velocity at time t for position(s) (x, z)
        """
        pass
    
    def elevation_cpp(self):
        return '%r + %r / 2.0 * cos(2 * pi * x[0] / %r)' % (self.depth,
                                                            self.height,
                                                            self.length)
