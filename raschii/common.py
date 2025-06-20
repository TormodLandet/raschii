import math
from math import pi, tanh
import numpy as np


# If the air phase blending_height is None then the wave height times this
# default factor will be used
AIR_BLENDING_HEIGHT_FACTOR = 2


class RaschiiError(Exception):
    pass


class NonConvergenceError(RaschiiError):
    pass


def check_breaking_criteria(
    height: float, depth: float, length: float | None = None, period: float | None = None
):
    """
    Return two empty strings if everything is OK, else a string with
    warnings about breaking criteria and a string with warnings about
    being close to a breaking criterion

    * height: wave height above still water level
    * depth: still water distance from the flat sea bottom to the free surface
        in meters, but you can give -1.0 for infinite depth
    * length: the periodic length of the wave (optional, if not given then period is used)
    * period: the wave period (optional, if not given then length is used)
      Since we need the wave length we assume Airy to convert period to length!
    """
    if length is None:
        if period is None:
            raise RaschiiError("Either length or period must be given")

        # We do not know the wave model, so we assume it is Airy
        # which should give a ballpark OK-ish answer
        from .airy import compute_length_from_period

        length = compute_length_from_period(depth=depth, period=period)
    
    if depth < 0.0:
        # Use a large depth for infinite depth waves, same as is used
        # in the stokes_coefficients and fenton_coefficients functions
        depth = 25 * length

    h1 = 0.14 * length
    h2 = 0.78 * depth
    h3 = 0.142 * tanh(2 * pi * depth / length) * length

    err = warn = ""
    for name, hmax in (
        ("Length criterion", h1),
        ("Depth criterion", h2),
        ("Combined criterion", h3),
    ):
        if height > hmax:
            err += "%s is exceeded, %.2f > %.2f\n" % (name, height, hmax)
        elif height > hmax * 0.9:
            warn += "%s is close to exceeded, %.2f = %.2f * %.3f\n" % (
                name,
                height,
                hmax,
                height / hmax,
            )
    return err, warn


def sinh_by_cosh(a, b):
    """
    A version of sinh(a)/cosh(b) where "b = a * f" and f is close
    to 1. This can then be written exp(a * (1 - f)) for large a
    """
    ans = np.zeros(a.size, float)
    for i, (ai, bi) in enumerate(zip(a, b)):
        if ai == 0:
            continue
        f = bi / ai
        if (ai > 30 and 0.5 < f < 1.5) or (ai > 200 and 0.1 < f < 1.9):
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
    ans = np.zeros(a.size, float)
    for i, (ai, bi) in enumerate(zip(a, b)):
        if ai == 0:
            ans[i] = 1.0 / math.cosh(bi)
            continue
        f = bi / ai
        if (ai > 30 and 0.5 < f < 1.5) or (ai > 200 and 0.1 < f < 1.9):
            ans[i] = math.exp(ai * (1 - f))
        else:
            ca = math.cosh(ai)
            cb = math.cosh(bi)
            ans[i] = ca / cb
    return ans


def blend_air_and_wave_velocities(x, z, t, wave, air, vel, eta_eps):
    """
    Compute velocities in the air phase and blend the water and air velocities
    in a divergence free manner up a distance ``air.blending_height - eta(x)``.
    If this is ``air.blending_height `` is ``None`` then blend all the way up to
    ``air.height``.

    The blending is done as follows. Introduce a new coordinate Z which is zero
    on the free surface and 1 at air_blend_distance. Then the blending function
    a smooth step function of a coordinate Z which is zero at the free surface
    and 1 at the blending heigh. After taking the derivatives of this blended
    stream function the velocities are as can be seen in the code below.
    """
    eta = wave.surface_elevation(x, t)
    above = z > eta + eta_eps
    if air is not None and above.any():
        d = air.blending_height + wave.depth
        xa = x[above]
        za = z[above]
        ea = eta[above]
        vel_air = air.velocity(xa, za, t)

        blend = za < d
        if d > 0 and blend.any():
            xb = xa[blend]
            zb = za[blend]
            eb = ea[blend]
            Z = (zb - eb) / (d - eb)
            vel_water = vel[above]
            psi_wave = wave.stream_function(xb, zb, t, frame="c")
            psi_air = air.stream_function(xb, zb, t, frame="c")
            detadx = wave.surface_slope(xb, t)

            if False:
                # Cubic smoothstep
                f = Z * Z * (3 - 2 * Z)
                dfdZ = 6 * Z - 6 * Z * Z
            else:
                # Fift order smootherstep
                f = Z * Z * Z * (Z * (Z * 6 - 15) + 10)
                dfdZ = Z * Z * (30 + Z * (30 * Z - 60))

            dZdx = (zb - d) / (d - eb) ** 2 * detadx
            dZdz = 1 / (d - eb)
            dfdx = dfdZ * dZdx
            dfdz = dfdZ * dZdz

            vel_air[blend, 0] = (1 - f) * vel_water[blend, 0] + f * vel_air[blend, 0]
            vel_air[blend, 1] = (1 - f) * vel_water[blend, 1] + f * vel_air[blend, 1]
            vel_air[blend, 0] += -dfdz * psi_wave + dfdz * psi_air
            vel_air[blend, 1] -= -dfdx * psi_wave + dfdx * psi_air
        vel[above] = vel_air
    else:
        vel[above] = 0


def blend_air_and_wave_velocity_cpp(
    wave_cpp,
    air_cpp,
    elevation_cpp,
    direction,
    eta_eps,
    air=None,
    psi_wave_cpp=None,
    psi_air_cpp=None,
    slope_cpp=None,
):
    """
    Return C++ code which blends the velocities in the water into the velocities
    in the air in such a way that the C++ code replicates the Python results
    from the blend_air_and_wave_velocities() function
    """
    if air is None:
        return "x[2] < (%s) + %r ? (%s) : (%s)" % (elevation_cpp, eta_eps, wave_cpp, "0.0")

    cpp = """[&]() {{
        const double elev = ({ecpp});

        const double val_water = ({uw_cpp});
        if (x[2] < elev + {eps!r}) {{
            // The point is below the free surface
            return val_water;
        }}

        const double d = {d!r};
        const double val_air = ({ua_cpp});
        if (x[2] < d) {{
            // The point is in the blending zone
            const double Z = (x[2] - elev) / (d - elev);

            // Cubic smoothstep
            //const double f = Z * Z * (3 - 2 * Z);
            //const double dfdZ = Z * (6 - 6 * Z);

            // Fift order smootherstep
            const double f = Z * Z * Z * (Z * (Z * 6 - 15) + 10);
            const double dfdZ = Z * Z * (30 + Z * (30 * Z - 60));

            // Derivatives needed in the product rules
            const double dZdx = (x[2] - d) / pow(d - elev, 2) * ({slope_cpp});
            const double dZdz = 1.0 / (d - elev);
            const double dfdx = dfdZ * dZdx;
            const double dfdz = dfdZ * dZdz;

            return (1 - f) * val_water + f * val_air + %r * (%s);
        }}

        // The point is above the blending zone
        return val_air;
    }}()""".format(
        ecpp=elevation_cpp,
        uw_cpp=wave_cpp,
        ua_cpp=air_cpp,
        eps=eta_eps,
        d=air.blending_height + air.depth_water,
        slope_cpp=slope_cpp,
    )

    # Compute the product rule code
    if direction == "x":
        sign = -1
        prcode = "dfdz * (%s) - dfdz * (%s)" % (psi_wave_cpp, psi_air_cpp)
    elif direction == "z":
        sign = +1
        prcode = "dfdx * (%s) - dfdx * (%s)" % (psi_wave_cpp, psi_air_cpp)

    return cpp % (sign, prcode)


def trapezoid_integration(*argv, **kwargs):
    """
    Compatibility for numpy 2.0 rename of np.trapz to np.trapezoid
    """
    if hasattr(np, "trapezoid"):
        return np.trapezoid(*argv, **kwargs)
    else:
        return np.trapz(*argv, **kwargs)


def np2py(val):
    """
    Convert a numpy array or numpy number into a list of Python floats
    or a single Python float

    We want the base types so that we can call repr() and get something
    that is pure Python (works in C++ after repr() conversion to code
    string) and does not include the string "np.float64" or similar
    (which repr() will in numpy 2.0)
    """
    if hasattr(val, "tolist"):
        return val.tolist()  # Convert numpy array to Python list of float
    elif hasattr(val, "item"):
        return val.item()  # Convert numpy float64 to Python float
    else:
        return val  # Assume this is allready a Python float
