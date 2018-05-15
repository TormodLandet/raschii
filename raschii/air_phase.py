from numpy import pi, zeros, asarray, arange, sin, cos, sinh, cosh, newaxis
from numpy.linalg import solve
from .fenton import sinh_by_cosh


class StreamFunctionAirPhase:
    def __init__(self, x, eta, c, k, depth_water, depth_air):
        """
        Given a set of colocation points with precomputed surface elevations
        obtained from a wave model in the water phase, produce a stream function
        approximation of the velocities in the air phase.
        """
        self.x = x
        self.eta = eta
        self.c = c
        self.k = k
        self.depth_water = depth_water
        self.depth_air = depth_air
        B, Q = air_velocity_coefficients(x, eta, c, k, depth_water, depth_air)
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
        x = x - c * t
        z = top - z
        J = arange(1, N + 1)
        
        vel = zeros((x.size, 2), float)
        vel[:, 0] = k * (B * cos(J * k * x[:, newaxis]) *
                         cosh(J * k * z[:, newaxis]) /
                         cosh(J * k * self.depth_air)).dot(J)
        vel[:, 1] = k * (B * sin(J * k * x[:, newaxis]) *
                         sinh(J * k * z[:, newaxis]) /
                         cosh(J * k * self.depth_air)).dot(J)
        return vel


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
