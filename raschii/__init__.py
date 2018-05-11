__version__ = '1.0.0.dev0'


class RasciiError(Exception):
    pass


class NonConvergenceError(RasciiError):
    pass


# The wave models
WAVE_MODELS = {}
from .fenton import FentonWave
WAVE_MODELS['Fenton'] = FentonWave


def get_wave_model(wave_type):
    return WAVE_MODELS[wave_type]


def check_breaking_criteria(height, depth, length):
    """
    Return two empty strings if everything is OK, else a string with
    warnings about breaking criteria and a string with warnings about
    being close to a breaking criterion
    """
    from math import pi, tanh
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
