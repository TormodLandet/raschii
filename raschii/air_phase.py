from numpy import pi, zeros, asarray, arange, sin, cos, sinh, cosh, newaxis
from numpy.linalg import solve, lstsq
from .common import sinh_by_cosh, cosh_by_cosh


class StreamFunctionAirPhase:
    def __init__(self, wave, N, length, depth_water, depth_air):
        """
        Given a set of colocation points with precomputed surface elevations
        obtained from a wave model in the water phase, produce a stream function
        approximation of the velocities in the air phase.
        """
        self.depth_water = depth_water
        self.depth_air = depth_air
        self.N = N
        self.x = arange(N + 1) * length / (2 * N)
        self.eta = wave.surface_elevation(self.x)
        self.vel = wave.velocity(self.x, self.eta)
        self.c = wave.c
        self.k = wave.k
        BPQ = air_velocity_coefficients(self.x, self.eta, self.vel, self.c,
                                        self.k, depth_water, depth_air)
        self.B, self.Px, self.Pz, self.Q = BPQ
        
        vel2 = self.velocity(self.x, self.eta)
        print('vel:\n%r' % self.vel)
        print('vel2:\n%r' % vel2)
        print('vel - vel2:\n%r' % (self.vel - vel2,))
        print('vel2 / vel:\n%r' % (vel2 / self.vel,))
    
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
        
        Nm = len(x)
        J = arange(1, len(self.B) + 1)
        D = self.depth_air
        x = x - self.c * t
        z = self.depth_water + self.depth_air - z
        k = self.k
        
        vel = zeros((x.size, 2), float)
        for m in range(Nm):
            S1 = sinh_by_cosh(J * k * z[m], J * k * D)
            C1 = cosh_by_cosh(J * k * z[m], J * k * D)
            S2 = sin(J * k * x[m])
            C2 = cos(J * k * x[m])
            
            Psi = S1 * C2
            u = J * k * C1 * C2
            v = J * k * S1 * S2
            vel[m, 0] = u.dot(self.B) + u.dot(self.Px) * x[m] + (u * z[m] + Psi).dot(self.Pz)
            vel[m, 1] = v.dot(self.B) + v.dot(self.Pz) * z[m] + (v * x[m] - Psi).dot(self.Px)
        
        return vel
        
        
        B = self.B
        k = self.k
        c = self.c
        top = self.depth_water + self.depth_air
        J = arange(1, B.size + 1)
        D = self.depth_air
        x2 = x[:, newaxis] - c * t
        z2 = top - z[:, newaxis]
        
        B2 = self.B + self.Px * x2 + self.Pz * z2
        Psi = B2 * (sinh(J * k * z2) / cosh(J * k * D) * cos(J * k * x2))
        
        vel = zeros((x.size, 2), float)
        vel[:, 0] = (k * B2 * cos(J * k * x2) * cosh(J * k * z2) /
                     cosh(J * k * D)).dot(J) + Psi.dot(self.Pz)
        vel[:, 1] = (k * B2 * sin(J * k * x2) * sinh(J * k * z2) /
                     cosh(J * k * D)).dot(J) - Psi.dot(self.Px)
        
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


def air_velocity_coefficients(x, eta, vel, c, k, depth_water, depth_air):
    """
    This uses the same method as in M.M.Rienecker and J. D. Fenton (1981), but
    since the surface elvation and phase speed is known the problem is now
    linear in the unknowns B1..BN and Q
    """
    Nm = len(eta)
    Nj = Nm - 1
    Neq = Nm * 3
    Nuk = 3*Nj + 1 
    J = arange(1, Nj + 1)
    D = depth_air
    z = depth_water + depth_air - eta
    
    lhs = zeros((Neq, Nuk), float)
    rhs = zeros(Neq, float)
    for m in range(Nm):
        S1 = sinh_by_cosh(J * k * z[m], J * k * D)
        C1 = cosh_by_cosh(J * k * z[m], J * k * D)
        S2 = sin(J * k * x[m])
        C2 = cos(J * k * x[m])
        
        Psi = S1 * C2
        u = J * k * C1 * C2
        v = J * k * S1 * S2
        
        # The free surface is a stream line (stream func = const Q)
        lhs[m, :Nj] = Psi
        lhs[m, -1] = 1
        rhs[m] = - c * z[m]
        #lhs[m, 1 * Nj:2 * Nj] = Psi * x[m]
        #lhs[m, 2 * Nj:3 * Nj] = Psi * z[m]
        
        # Continuity of the velocity in the x-direction
        m2 = m + Nm
        lhs[m2, :Nj] = u
        rhs[m2] = vel[m, 0]
        #lhs[m, 1 * Nj:2 * Nj] = u * x[m]
        #lhs[m, 2 * Nj:3 * Nj] = Psi + u * z[m]
        
        # Continuity of the velocity in the y-direction
        m3 = m2 + Nm
        lhs[m3, :Nj] = v
        rhs[m3] = -vel[m, 1]
        #lhs[m, 1 * Nj:2 * Nj] = - Psi + v * x[m]
        #lhs[m, 2 * Nj:3 * Nj] = v * z[m]
    
    BPQ, *_  = lstsq(lhs, rhs, rcond=None)
    B = BPQ[:Nj]
    Px = BPQ[Nj:2 * Nj]
    Pz = BPQ[2 * Nj: 3 * Nj]
    Q = BPQ[-1]
    print('B', B, B.shape)
    print('Px', Px, Px.shape)
    print('Pz', Pz, Pz.shape)
    print('Q', Q)
    
    if False:
        lhs2 = lhs[:Nm,:Nm].copy()
        lhs2[:,-1] = 1
        BQ = solve(lhs2, rhs[:Nm])
        B = BQ[:-1]
        Q = BQ[-1]
        Px = Pz = B * 0
    
    return B, Px, Pz, Q
