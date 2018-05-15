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
