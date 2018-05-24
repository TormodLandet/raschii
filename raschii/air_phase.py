from numpy import zeros, asarray, arange, sin, cos, sinh, cosh, newaxis
from numpy.linalg import solve
from .common import sinh_by_cosh


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
        self.c = wave.c
        self.k = wave.k
        BQ = air_velocity_coefficients(self.x, self.eta, self.c, self.k,
                                       depth_water, depth_air)
        self.B, self.Q = BQ
    
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
        top = self.depth_water + self.depth_air
        J = arange(1, B.size + 1)
        D = self.depth_air
        x2 = x[:, newaxis] - c * t
        z2 = top - z[:, newaxis]
        
        vel = zeros((x.size, 2), float)
        vel[:, 0] = (k * B * cos(J * k * x2) * cosh(J * k * z2) /
                     cosh(J * k * D)).dot(J)
        vel[:, 1] = (k * B * sin(J * k * x2) * sinh(J * k * z2) /
                     cosh(J * k * D)).dot(J)
        
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
    Nm = len(eta)
    Nj = Nm - 1
    Neq = Nm
    Nuk = Nj + 1
    J = arange(1, Nj + 1)
    D = depth_air
    z = depth_water + depth_air - eta
    
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


def blend_air_and_wave_velocities(x, z, t, wave, air, vel, eta_eps,
                                  air_blend_distance, include_air_phase):
    """
    Compute velocities in the air phase and blend the water and air velocities
    in a divergence free manner from the free surface and a distance
    ``air_blend_distance`` up.
    
    The blending is done as follows. Introduce a new coordinate Z which is zero
    on the free surface and 1 at air_blend_distance above the free surface. Then
    the blending function is ``f = 3Z² - 2Z³`` which is used on the stream
    functions in the water and in the air. After taking the derivatives of this
    blended stream function the velocities are as can be seen in the code below.
    """
    zmax = wave.surface_elevation(x, t)
    above = z > zmax + eta_eps
    if include_air_phase and above.any():
        xa = x[above]
        za = z[above]
        ea = zmax[above]
        vel_air = air.velocity(xa, za, t)
        
        blend = za < ea + air_blend_distance
        if air_blend_distance > 0 and blend.any():
            zb = za[blend]
            eb = ea[blend]
            Z = (zb - eb) / air_blend_distance
            fw = 2 * Z**3 - 3 * Z**2 + 1
            fa = - Z**2 * (2 * Z - 3)
            vel_water = vel[above]
            vel_air[blend, 0] = fw * vel_water[blend, 0] + fa * vel_air[blend, 0]
            vel_air[blend, 1] = fw * vel_water[blend, 1] + fa * vel_air[blend, 1]
        vel[above] = vel_air
    else:
        vel[above] = 0


def blend_air_and_wave_velocity_cpp(wave_cpp, air_cpp, elevation_cpp, eta_eps,
                                    air_blend_distance, include_air_phase):
    """
    Return C++ code which blends the velocities in the water into the velocities
    in the air in such a way that the C++ code replicates the Python results
    from the blend_air_and_wave_velocities() function
    """
    if not include_air_phase:
        return 'x[2] < (%s) + %r ? (%s) : (%s)' % (elevation_cpp, eta_eps,
                                                   wave_cpp, '0.0')
    
    return """[&]() {
        const double elev = (%s);
        
        const double val_water = (%s);
        if (x[2] < elev + %r) {
            // The point is below the free surface
            return val_water;
        }
        
        const double dist_blend = %r;
        const double val_air = (%s);
        if (x[2] < elev + dist_blend) {
            // The point is in the blending zone
            const double Z = (x[2] - elev) / dist_blend;
            const double fw = 2*Z*Z*Z  - 3*Z*Z + 1;
            const double fa = - Z*Z * (2*Z - 3);
            return fw * val_water + fa * val_air;
        }
        
        // The point is above the blending zone
        return val_air;
    }()""" % (elevation_cpp, wave_cpp, eta_eps, air_blend_distance, air_cpp)
