from numpy import zeros, asarray, arange, sin, cos, sinh, cosh, newaxis
from numpy.linalg import solve
from .common import sinh_by_cosh


class FentonAirPhase:
    def __init__(self, height, blending_height=None):
        """
        Given a set of colocation points with precomputed surface elevations
        obtained from a wave model in the water phase, produce a stream function
        approximation of the velocities in the air phase.
        """
        self.height = height
        self.blending_height = blending_height
    
    def set_wave(self, wave):
        """
        Connect this air phase with the wave in the water phase
        """
        N = wave.order
        self.x = arange(N + 1) * wave.length / (2 * N)
        self.eta = wave.surface_elevation(self.x)
        self.c = wave.c
        self.k = wave.k
        self.g = wave.g
        BQ = air_velocity_coefficients(self.x, self.eta, self.c, self.k,
                                       wave.depth, self.height)
        self.B, self.Q = BQ
        self.depth_water = wave.depth
        
        if self.blending_height is None:
            self.blending_height = 2 * wave.height
    
    def stream_function(self, x, z, t=0, frame='b'):
        """
        Compute the stream function at time t for position(s) x
        """
        if isinstance(x, (float, int)):
            x, z = [x], [z]
        x2 = asarray(x, dtype=float) - self.c * t
        z2 = self.depth_water + self.height - asarray(z, dtype=float)
        x2, z2 = x2[:, newaxis], z2[:, newaxis]
        
        N = len(self.eta) - 1
        B0 = self.c
        B = self.B
        k = self.k
        J = arange(1, N + 1)
        
        psi = (sinh(J * k * z2) / cosh(J * k * self.height) *
               cos(J * k * x2)).dot(B)
        
        if frame == 'e':
            return B0 * z + psi
        elif frame == 'c':
            return psi
    
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
        
        B = self.B
        k = self.k
        c = self.c
        top = self.depth_water + self.height
        J = arange(1, B.size + 1)
        D = self.height
        x2 = x[:, newaxis] - c * t
        z2 = top - z[:, newaxis]
        
        vel = zeros((x.size, 2), float)
        vel[:, 0] = (k * B * cos(J * k * x2) * cosh(J * k * z2) /
                     cosh(J * k * D)).dot(J) * -1
        vel[:, 1] = (k * B * sin(J * k * x2) * sinh(J * k * z2) /
                     cosh(J * k * D)).dot(J)
        
        return vel
    
    def stream_function_cpp(self, frame='b'):
        """
        Return C++ code for evaluating the stream function of this specific
        wave. The positive traveling direction is x[0] and the vertical
        coordinate is x[2] which is zero at the bottom and equal to +depth at
        the mean water level.
        """
        N = len(self.eta) - 1
        J = arange(1, N + 1)
        k = self.k
        c = self.c
        
        Jk = J * k
        facs = self.B / cosh(Jk * self.height)
        
        z2 = '(%r - x[2])' % (self.depth_water + self.height, )
        cpp = ' + '.join('%r * cos(%f * (x[0] - %r * t)) * sinh(%r * %s)' %
                         (facs[i], Jk[i], c, Jk[i], z2) for i in range(N))
        
        if frame == 'b':
            B0 = self.c
            return '%r * x[2] + %s' % (B0, cpp)
        elif frame == 'c':
            return cpp
    
    def velocity_cpp(self):
        """
        Return C++ code for evaluating the particle velocities of this specific
        wave. Returns the x and z components only with z positive upwards. The
        positive traveling direction is x[0] and the vertical coordinate is x[2]
        which is zero at the bottom and equal to +depth at the mean water level.
        """
        N = len(self.eta) - 1
        J = arange(1, N + 1)
        k = self.k
        c = self.c
        
        Jk = J * k
        facs = J * self.B * k / cosh(Jk * self.height)
        
        z2 = '(%r - x[2])' % (self.depth_water + self.height, )
        cpp_x = ' + '.join('%r * cos(%f * (x[0] - %r * t)) * cosh(%r * %s)' %
                           (facs[i], Jk[i], c, Jk[i], z2) for i in range(N))
        cpp_z = ' + '.join('%r * sin(%f * (x[0] - %r * t)) * sinh(%r * %s)' %
                           (facs[i], Jk[i], c, Jk[i], z2) for i in range(N))
        return (cpp_x, cpp_z)
    
    def __repr__(self):
        return ('FentonAirPhase(height={s.height}, blending_height='
                '{s.blending_height})').format(s=self)


def air_velocity_coefficients(x, eta, c, k, depth_water, height_air):
    """
    This uses the same method as in M.M.Rienecker and J. D. Fenton (1981), but
    since the surface elvation and phase speed is known the problem is now
    linear in the unknowns B1..BN and Q
    """
    Nm = len(eta)
    Nj = Nm - 1
    Neq = Nm
    Nuk = Nj + 1
    J = arange(1, Nj + 1)
    D = height_air
    z = depth_water + height_air - eta
    
    lhs = zeros((Neq, Nuk), float)
    rhs = zeros(Neq, float)
    for m in range(Nm):
        S1 = sinh_by_cosh(J * k * z[m], J * k * D)
        C2 = cos(J * k * x[m])
        
        # The free surface is a stream line (stream func = const Q)
        lhs[m, :Nj] = S1 * C2
        lhs[m, -1] = 1
        rhs[m] = - c * z[m]
    
    BQ = solve(lhs, rhs)
    B = BQ[:-1]
    Q = BQ[-1]
    
    return B, Q
