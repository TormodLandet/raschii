from numpy import pi, cos, sin, zeros, array, asarray, sinh, cosh, tanh
from .air_phase import StreamFunctionAirPhase


class AiryWave:
    required_input = ('height', 'depth', 'length')
    optional_input = {'g': 9.81, 'depth_air': 0}
    
    def __init__(self, height, depth, length, g=9.81, depth_air=0):
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
        self.depth_air = depth_air
        self.include_air_phase = (depth_air > 0)
        self.warnings = ''
        
        self.k = 2 * pi / length
        self.omega = (self.k * g * tanh(self.k * depth))**0.5
        self.c = self.omega / self.k
        
        # Provide velocities also in the air phase
        if self.include_air_phase:
            self.air = StreamFunctionAirPhase(self, 1, length, depth, depth_air)
        
        # For evaluating velocities close to the free surface
        self.eta_eps = self.height / 1e5
    
    def surface_elevation(self, x, t=0):
        """
        Compute the surface elavation at time t for position(s) x
        """
        if isinstance(x, (float, int)):
            x = array([x], float)
        x = asarray(x)
        return self.depth + self.height / 2 * cos(self.k * x - self.omega * t)
    
    def velocity(self, x, z, t=0):
        """
        Compute the fluid velocity at time t for position(s) (x, z)
        where z is 0 at the bottom and equal to depth at the free surface
        """
        if isinstance(x, (float, int)):
            x, z = [x], [z]
        x = asarray(x, dtype=float)
        z = asarray(z, dtype=float)
        
        H = self.height
        k = self.k
        d = self.depth
        w = self.omega
        
        vel = zeros((x.size, 2), float) + 1
        vel[:, 0] = w * H / 2 * cosh(k * z) / sinh(k * d) * cos(k * x - w * t)
        vel[:, 1] = w * H / 2 * sinh(k * z) / sinh(k * d) * sin(k * x - w * t)
        zmax = self.surface_elevation(x, t)
        
        above = z > zmax + self.eta_eps
        if self.include_air_phase:
            vel_air = self.air.velocity(x[above], z[above], t)
            vel[above] = vel_air
        else:
            vel[above] = 0
        
        return vel
    
    def elevation_cpp(self):
        """
        Return C++ code for evaluating the elevation of this specific wave.
        The positive traveling direction is x[0]
        """
        return '%r + %r / 2.0 * cos(%r * (x[0] - %r * t))' % \
            (self.depth, self.height, self.k, self.c)
    
    def velocity_cpp(self):
        """
        Return C++ code for evaluating the particle velocities of this specific
        wave. Returns the x and z components only with z positive upwards. The
        positive traveling direction is x[0] and the vertical coordinate is x[2]
        which is zero at the bottom and equal to +depth at the mean water level.
        """
        H = self.height
        k = self.k
        d = self.depth
        w = self.omega
        
        cpp_x = '%r * cosh(%r * x[2]) * cos(%r * x[0] - %r * t)' %\
                (w * H / (2 * sinh(k * d)), k, k, w)
        cpp_z = '%r * sinh(%r * x[2]) * sin(%r * x[0] - %r * t)' %\
                (w * H / (2 * sinh(k * d)), k, k, w)
        e_cpp = self.elevation_cpp()
        
        return ('x[2] < (%s) + %r ? (%s) : 0.0' % (e_cpp, self.eta_eps, cpp_x),
                'x[2] < (%s) + %r ? (%s) : 0.0' % (e_cpp, self.eta_eps, cpp_z))
