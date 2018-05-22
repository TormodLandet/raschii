"""
Raschii is a library for constructing non-linear regular waves

Raschii is written by Tormod Landet (c) 2018- and is named after Thysanoessa
raschii, the Arctic Krill.

SPDX-License-Identifier: Apache-2.0
"""
__version__ = '1.0.0.dev0'
from .common import check_breaking_criteria, RasciiError, NonConvergenceError
from .air_phase import StreamFunctionAirPhase
from .airy import AiryWave
from .fenton import FentonWave
from .stokes import StokesWave


# The available wave models
WAVE_MODELS = {
    'Airy': AiryWave,
    'Fenton': FentonWave,
    'Stokes': StokesWave
}


def get_wave_model(model_name):
    """
    Get a Raschii wave model by name
    """
    if model_name in WAVE_MODELS:
        return WAVE_MODELS[model_name]
    else:
        raise RasciiError('Wave model %r is not supported, supported wave '
                          'models are %s' % (model_name,
                                             ', '.join(WAVE_MODELS.keys())))
