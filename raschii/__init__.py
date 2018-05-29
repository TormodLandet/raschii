"""
Raschii is a library for constructing non-linear regular waves

Raschii is written by Tormod Landet (c) 2018- and is named after Thysanoessa
raschii, the Arctic Krill.

SPDX-License-Identifier: Apache-2.0
"""
__version__ = '1.0.0'
from .common import check_breaking_criteria, RasciiError, NonConvergenceError
from .airy import AiryWave
from .fenton import FentonWave
from .stokes import StokesWave
from .air_phase_fenton import FentonAirPhase
from .air_phase_constant import ConstantAirPhase


# The available wave models
WAVE_MODELS = {
    'Airy': AiryWave,
    'Fenton': FentonWave,
    'Stokes': StokesWave
}

# Air phase models
AIR_MODELS = {'FentonAir': FentonAirPhase,
              'ConstantAir': ConstantAirPhase}


def get_wave_model(model_name, air_model_name=None):
    """
    Get a Raschii wave model by name
    """
    if '+' in model_name:
        model_name, air_model_name = model_name.split('+')
    
    if model_name not in WAVE_MODELS:
        raise RasciiError('Wave model %r is not supported, supported wave '
                          'models are %s' % (model_name,
                                             ', '.join(WAVE_MODELS.keys())))
    wave = WAVE_MODELS[model_name]
    
    if air_model_name is None:
        return wave, None
    
    if air_model_name not in AIR_MODELS:
        raise RasciiError('Air model %r is not supported, supported air phase '
                          'models are %s' % (air_model_name,
                                             ', '.join(AIR_MODELS.keys())))
    air = AIR_MODELS[air_model_name]
    
    return wave, air
