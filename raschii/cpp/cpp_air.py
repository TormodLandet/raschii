from numpy import arange, cosh

from ..common import Frame, np2py


class ConstantAirCppGenerator:
    """C++ code generator for ConstantAirPhase."""

    def __init__(self, air):
        self.air = air

    def stream_function(self, frame=Frame.EARTH):
        """
        Return C++ code for evaluating the stream function of this specific
        air phase model.
        """
        c = np2py(self.air.c)
        if frame == Frame.EARTH:
            return "0.0"
        elif frame == Frame.WAVE:
            return f"-{c!r} * x[2]"

    def velocity(self):
        """
        Return C++ code for evaluating the particle velocities of this specific
        air phase model.
        """
        return ("0.0", "0.0")


class FentonAirCppGenerator:
    """C++ code generator for FentonAirPhase."""

    def __init__(self, air):
        self.air = air

    def stream_function(self, frame=Frame.EARTH):
        """
        Return C++ code for evaluating the stream function of this specific
        air phase model. The positive traveling direction is x[0] and the
        vertical coordinate is x[2] which is zero at the bottom and equal to
        +depth at the mean water level.
        """
        air = self.air
        N = len(air.eta) - 1
        J = arange(1, N + 1)
        k = air.k

        Jk = J * k
        facs = air.B / cosh(Jk * air.height)

        # Repr of np.float64(42.0) is "np.float64(42.0)" and not "42.0"
        # We use repr to make Python output a "smart" amount of digits
        z2 = np2py(air.depth_water + air.height)
        c = np2py(air.c)
        Jk = np2py(Jk)
        facs = np2py(facs)

        z2_cpp = f"({z2!r} - x[2])"
        cpp = " + ".join(
            f"{facs[i]!r} * cos({Jk[i]!r} * (x[0] - {c!r} * t)) * sinh({Jk[i]!r} * {z2_cpp})"
            for i in range(N)
        )

        if frame == Frame.EARTH:
            B0 = np2py(air.c)
            return f"{B0!r} * x[2] + {cpp}"
        elif frame == Frame.WAVE:
            return cpp

    def velocity(self):
        """
        Return C++ code for evaluating the particle velocities of this specific
        air phase model. Returns the x and z components only with z positive
        upwards. The positive traveling direction is x[0] and the vertical
        coordinate is x[2] which is zero at the bottom and equal to +depth at
        the mean water level.
        """
        air = self.air
        N = len(air.eta) - 1
        J = arange(1, N + 1)
        k = air.k

        Jk = J * k
        facs = J * air.B * k / cosh(Jk * air.height)

        # Repr of np.float64(42.0) is "np.float64(42.0)" and not "42.0"
        # We use repr to make Python output a "smart" amount of digits
        z2 = np2py(air.depth_water + air.height)
        c = np2py(air.c)
        Jk = np2py(Jk)
        facs = np2py(facs)

        z2_cpp = f"({z2!r} - x[2])"
        cpp_x = " + ".join(
            f"{-facs[i]!r} * cos({Jk[i]!r} * (x[0] - {c!r} * t)) * cosh({Jk[i]!r} * {z2_cpp})"
            for i in range(N)
        )
        cpp_z = " + ".join(
            f"{facs[i]!r} * sin({Jk[i]!r} * (x[0] - {c!r} * t)) * sinh({Jk[i]!r} * {z2_cpp})"
            for i in range(N)
        )
        return (cpp_x, cpp_z)
