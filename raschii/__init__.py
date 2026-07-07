"""
Raschii is a library for constructing non-linear regular waves

Raschii is written by Tormod Landet (c) 2018- and is named after Thysanoessa
raschii, the Arctic Krill.

SPDX-License-Identifier: Apache-2.0
"""

from typing import Literal, TypeAlias, overload

#: The three-digit version number of Raschii
__version__ = "2.0.0"
version = __version__

from .air_phase_constant import ConstantAirPhase
from .air_phase_fenton import FentonAirPhase
from .base_classes import AirPhaseModel, Frame, WaveModel
from .common import NonConvergenceError, RaschiiError, check_breaking_criteria
from .wave_airy import AiryWave
from .wave_fenton import FentonWave
from .wave_stokes import StokesWave

# Type aliases for the wave and air-phase model names
WaveModelName: TypeAlias = Literal["Airy", "Fenton", "Stokes"]
AirModelName: TypeAlias = Literal["FentonAir", "ConstantAir"]

#: The available wave models. A dictionary mapping model name (str) to model class.
WAVE_MODELS: dict[WaveModelName, type[WaveModel]] = {
    "Airy": AiryWave,
    "Fenton": FentonWave,
    "Stokes": StokesWave,
}
#: The available air-phase models. A dictionary mapping model name (str) to model class.
AIR_MODELS: dict[AirModelName, type[AirPhaseModel]] = {
    "FentonAir": FentonAirPhase,
    "ConstantAir": ConstantAirPhase,
}


# Return value when both wave model and air model are requested
@overload
def get_wave_model(
    model_name: WaveModelName, air_model_name: AirModelName
) -> tuple[type[WaveModel], type[AirPhaseModel]]: ...


# Return value when only wave model is requested
@overload
def get_wave_model(
    model_name: WaveModelName, air_model_name: None
) -> tuple[type[WaveModel], None]: ...


def get_wave_model(
    model_name: WaveModelName, air_model_name: AirModelName | None = None
) -> tuple[
    type[WaveModel],
    type[AirPhaseModel] | None,
]:
    """
    Get a Raschii wave model by name (returns the class, not an instance)

    If you want to use an air-phase model, you can specify it with the `air_model_name` argument
    and the air-model class will be returned as the second element of the tuple, otherwise the
    second returned element will be None.

    Returns a tuple of (wave_model_class, air_model_class | None)
    """
    # Old style: "Fenton+FentonAir" is equivalent to model_name="Fenton", air_model_name="FentonAir"
    if "+" in model_name:
        assert air_model_name is None
        model_name, air_model_name = model_name.split("+")

    if model_name not in WAVE_MODELS:
        raise RaschiiError(
            "Wave model %r is not supported, supported wave "
            "models are %s" % (model_name, ", ".join(WAVE_MODELS.keys()))
        )
    wave = WAVE_MODELS[model_name]

    if air_model_name is None:
        return wave, None

    if air_model_name not in AIR_MODELS:
        raise RaschiiError(
            "Air model %r is not supported, supported air phase "
            "models are %s" % (air_model_name, ", ".join(AIR_MODELS.keys()))
        )
    air = AIR_MODELS[air_model_name]

    return wave, air
