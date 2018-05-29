import math
from math import pi, tanh
import numpy


class RasciiError(Exception):
    pass


class NonConvergenceError(RasciiError):
    pass


def check_breaking_criteria(height, depth, length):
    """
    Return two empty strings if everything is OK, else a string with
    warnings about breaking criteria and a string with warnings about
    being close to a breaking criterion
    """
    h1 = 0.14 * length
    h2 = 0.78 * depth
    h3 = 0.142 * tanh(2 * pi * depth / length) * length
    
    err = warn = ''
    for name, hmax in (('Length criterion', h1),
                       ('Depth criterion', h2),
                       ('Combined criterion', h3)):
        if height > hmax:
            err += '%s is exceeded, %.2f > %.2f\n' % (name, height, hmax)
        elif height > hmax * 0.9:
            warn += '%s is close to exceeded, %.2f = %.2f * %.3f\n' % \
                (name, height, hmax, height / hmax)
    return err, warn


def sinh_by_cosh(a, b):
    """
    A version of sinh(a)/cosh(b) where "b = a * f" and f is close
    to 1. This can then be written exp(a * (1 - f)) for large a
    """
    ans = numpy.zeros(a.size, float)
    for i, (ai, bi) in enumerate(zip(a, b)):
        if ai == 0:
            continue
        f = bi / ai
        if ((ai > 30 and 0.5 < f < 1.5) or (ai > 200 and 0.1 < f < 1.9)):
            ans[i] = math.exp(ai * (1 - f))
        else:
            sa = math.sinh(ai)
            cb = math.cosh(bi)
            ans[i] = sa / cb
    return ans


def cosh_by_cosh(a, b):
    """
    A version of cosh(a)/cosh(b) where "b = a * f" and f is close
    to 1. This can then be written exp(a * (1 - f)) for large a
    """
    ans = numpy.ones(a.size, float)
    for i, (ai, bi) in enumerate(zip(a, b)):
        if ai == 0:
            continue
        f = bi / ai
        if ((ai > 30 and 0.5 < f < 1.5) or (ai > 200 and 0.1 < f < 1.9)):
            ans[i] = math.exp(ai * (1 - f))
        else:
            ca = math.cosh(ai)
            cb = math.cosh(bi)
            ans[i] = ca / cb
    return ans


def blend_air_and_wave_velocities(x, z, t, wave, air, vel, eta_eps):
    """
    Compute velocities in the air phase and blend the water and air velocities
    in a divergence free manner from the free surface and a distance
    ``air.blending_height`` up. If this is ``None`` then blend all the way up to
    ``air.height``.
    
    The blending is done as follows. Introduce a new coordinate Z which is zero
    on the free surface and 1 at air_blend_distance above the free surface. Then
    the blending function is ``f = 3Z² - 2Z³`` which is used on the stream
    functions in the water and in the air. After taking the derivatives of this
    blended stream function the velocities are as can be seen in the code below.
    """
    eta = wave.surface_elevation(x, t)
    above = z > eta + eta_eps
    if air is not None and above.any():
        xa = x[above]
        za = z[above]
        ea = eta[above]
        vel_air = air.velocity(xa, za, t)
        
        blend = za < ea + air.blending_height
        if air.blending_height > 0 and blend.any():
            xb = xa[blend]
            zb = za[blend]
            eb = ea[blend]
            Z = (zb - eb) / air.blending_height
            vel_water = vel[above]
            psi_wave = wave.stream_function(xb, zb, t, frame='c')
            psi_air = air.stream_function(xb, zb, t, frame='c')
            detadx = wave.surface_slope(xb, t)
            
            if False:
                # Cubic smoothstep
                f = Z * Z * (3 - 2 * Z)
                dfdZ = 6 * Z - 6 * Z * Z
            else:
                # Fift order smootherstep
                f = Z * Z * Z * (Z * (Z * 6 - 15) + 10)
                dfdZ = Z * Z * (30 + Z * (30 * Z - 60))
            
            dZdx = -1 / air.blending_height * detadx
            dZdz = 1 / air.blending_height
            dfdx = dfdZ * dZdx
            dfdz = dfdZ * dZdz
            
            vel_air[blend, 0] = (1 - f) * vel_water[blend, 0] + f * vel_air[blend, 0]
            vel_air[blend, 1] = (1 - f) * vel_water[blend, 1] + f * vel_air[blend, 1]
            vel_air[blend, 0] += - dfdz * psi_wave + dfdz * psi_air
            vel_air[blend, 1] -= - dfdx * psi_wave + dfdx * psi_air
        vel[above] = vel_air
    else:
        vel[above] = 0


def blend_air_and_wave_velocity_cpp(wave_cpp, air_cpp, elevation_cpp, direction,
                                    eta_eps, air=None, psi_wave_cpp=None,
                                    psi_air_cpp=None, slope_cpp=None):
    """
    Return C++ code which blends the velocities in the water into the velocities
    in the air in such a way that the C++ code replicates the Python results
    from the blend_air_and_wave_velocities() function
    """
    if air is None:
        return 'x[2] < (%s) + %r ? (%s) : (%s)' % (elevation_cpp, eta_eps,
                                                   wave_cpp, '0.0')
    
    cpp = """[&]() {{
        const double elev = ({ecpp});
        
        const double val_water = ({uw_cpp});
        if (x[2] < elev + {eps!r}) {{
            // The point is below the free surface
            return val_water;
        }}
        
        const double dist_blend = {dist_blend!r};
        const double val_air = ({ua_cpp});
        if (x[2] < elev + dist_blend) {{
            // The point is in the blending zone
            const double Z = (x[2] - elev) / dist_blend;
            
            // Cubic smoothstep
            //const double f = Z * Z * (3 - 2 * Z);
            //const double dfdZ = Z * (6 - 6 * Z);
            
            // Fift order smootherstep
            const double f = Z * Z * Z * (Z * (Z * 6 - 15) + 10);
            const double dfdZ = Z * Z * (30 + Z * (30 * Z - 60));
            
            // Derivatives needed in the product rules
            const double dZdx = - 1.0 / dist_blend * ({slope_cpp});
            const double dZdz = 1.0 / dist_blend;
            const double dfdx = dfdZ * dZdx;
            const double dfdz = dfdZ * dZdz;
            
            return (1 - f) * val_water + f * val_air + %r * (%s);
        }}
        
        // The point is above the blending zone
        return val_air;
    }}()""".format(ecpp=elevation_cpp, uw_cpp=wave_cpp, ua_cpp=air_cpp,
                   eps=eta_eps, dist_blend=air.blending_height,
                   slope_cpp=slope_cpp)
    
    # Compute the product rule code
    if direction == 'x':
        sign = -1
        prcode = 'dfdz * (%s) - dfdz * (%s)' % (psi_wave_cpp, psi_air_cpp)
    elif direction == 'z':
        sign = +1
        prcode = 'dfdx * (%s) - dfdx * (%s)' % (psi_wave_cpp, psi_air_cpp)
    
    return cpp % (sign, prcode)
