from numpy import pi, cos, sin, zeros, array, asarray, sinh, cosh, tanh


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
        self.k = 2 * pi / length
        self.omega = (self.k * g * tanh(self.k * depth))**0.5
        self.c = self.omega / self.k
    
    def surface_elevation(self, x, t=0):
        """
        Compute the surface elavation at time t for position(s) x
        """
        if isinstance(x, (float, int)):
            x = array([x], float)
        x = asarray(x)
        x = x + self.c * t
        return self.depth + self.height / 2 * cos(self.k * x)
    
    def velocity(self, x, z, t=0):
        """
        Compute the fluid velocity at time t for position(s) (x, z)
        """
        if isinstance(x, (float, int)):
            x, z = [x], [z]
        x = asarray(x, dtype=float)
        z = asarray(z, dtype=float)
        
        H = self.height
        c = self.c
        k = self.k
        d = self.depth
        w = self.omega
        x = x + c * t
        
        vel = zeros((x.size, 2), float) + 1
        vel[:, 0] = w * H / 2 * cosh(k * (z + d)) / sinh(k * d) * cos(k * x)
        vel[:, 1] = w * H / 2 * sinh(k * (z + d)) / sinh(k * d) * sin(k * x)
        zmax = self.surface_elevation(x, t)
        vel[z > zmax] = 0
        
        return vel
    
    def elevation_cpp(self):
        return '%r + %r / 2.0 * cos(%r * (x[0] + %r * t))' % \
            (self.depth, self.height, self.k, self.c)
