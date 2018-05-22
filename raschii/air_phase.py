from numpy import pi, zeros, asarray, arange, sin, cos, sinh, cosh, newaxis
from numpy.linalg import solve
from .common import sinh_by_cosh


class StreamFunctionAirPhase:
    def __init__(self, wave, N, length, depth_water, depth_air):
        """
        Given a set of colocation points with precomputed surface elevations
        obtained from a wave model in the water phase, produce a stream function
        approximation of the velocities in the air phase.
        """
        self.x = arange(N + 1) * length / (2 * N)
        self.eta = wave.surface_elevation(self.x)
        self.c = wave.c
        self.k = wave.k
        self.depth_water = depth_water
        self.depth_air = depth_air
        B, Q = air_velocity_coefficients(self.x, self.eta, self.c, self.k,
                                         depth_water, depth_air)
        self.B = B
        self.Q = Q
    
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
        
        N = len(self.eta) - 1
        B = self.B
        k = self.k
        c = self.c
        top = self.depth_water + self.depth_air
        z = top - z
        J = arange(1, N + 1)
        
        vel = zeros((x.size, 2), float)
        vel[:, 0] = k * (B * cos(J * k * x[:, newaxis] - c * t) *
                         cosh(J * k * z[:, newaxis]) /
                         cosh(J * k * self.depth_air)).dot(J)
        vel[:, 1] = k * (B * sin(J * k * x[:, newaxis] - c * t) *
                         sinh(J * k * z[:, newaxis]) /
                         cosh(J * k * self.depth_air)).dot(J)
        return vel
    
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
        facs = J * self.B * k / cosh(Jk * self.depth_air)
        
        x2 = '(%r - x[2])' % (self.depth_water + self.depth_air)
        cpp_x = ' + '.join('%r * cos(%f * (x[0] - %r * t)) * cosh(%r * %s)' %
                           (facs[i], Jk[i], c, Jk[i], x2) for i in range(N))
        cpp_z = ' + '.join('%r * sin(%f * (x[0] - %r * t)) * sinh(%r * %s)' %
                           (facs[i], Jk[i], c, Jk[i], x2) for i in range(N))
        return (cpp_x, cpp_z)


def air_velocity_coefficients(x, eta, c, k, depth_water, depth_air):
    """
    This uses the same method as in M.M.Rienecker and J. D. Fenton (1981), but
    since the surface elvation and phase speed is known the problem is now
    linear in the unknowns B1..BN and Q
    """
    N = len(eta) - 1
    J = arange(1, N + 1)
    D = depth_air
    
    eta = depth_water + depth_air - eta
    eta += (eta.max() - eta.min()) / 10
    
    lhs = zeros((N + 1, N + 1), float)
    rhs = zeros(N + 1, float)
    for m in range(0, N + 1):
        S1 = sinh_by_cosh(J * k * eta[m], J * k * D)
        C2 = cos(J * m * pi / N)
        
        lhs[m, :-1] = S1 * C2
        lhs[m, -1] = 1
        rhs[m] = - c * eta[m]
    
    BQ = solve(lhs, rhs)
    B = BQ[:-1]
    Q = BQ[-1]
    return B, Q
