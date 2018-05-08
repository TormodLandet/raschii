__version__ = '1.0.0.dev0'


class RasciiError(Exception):
    pass


class NonConvergenceError(RasciiError):
    pass


from .fenton import FentonWave
