"""
SWD (Spectral Wave Data) file writing and testing utilities.

Public API: there is no public API for raschii.swd. The classes and functions in this module are for
internal use by raschii only.
"""

from .swd_file import SwdShape1and2, SwdReaderForRaschiiTests
from .swd_fenton import SwdWriterFenton
from .swd_stokes import SwdWriterStokes
from .swd_airy import SwdWriterAiry
