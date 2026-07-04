from typing import TypeAlias
from .cpp_airy import AiryCppGenerator
from .cpp_fenton import FentonCppGenerator
from .cpp_stokes import StokesCppGenerator
from .cpp_air import ConstantAirCppGenerator, FentonAirCppGenerator

WaveModelCppGenerator: TypeAlias = AiryCppGenerator | FentonCppGenerator | StokesCppGenerator
AirPhaseModelCppGenerator: TypeAlias = ConstantAirCppGenerator | FentonAirCppGenerator
