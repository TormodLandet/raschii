from numpy import arange, cosh

from ..common import RaschiiError, blend_air_and_wave_velocity_cpp, np2py


class FentonCppGenerator:
    """C++ code generator for FentonWave."""

    def __init__(self, wave):
        self.wave = wave

    def stream_function(self, frame="b"):
        """
        Return C++ code for evaluating the stream function of this specific
        wave. The positive traveling direction is x[0] and the vertical
        coordinate is x[2] which is zero at the bottom and equal to +depth at
        the mean water level.
        """
        wave = self.wave
        if wave.depth < 0:
            raise RaschiiError("Cannot currently generate C++ code for infinite depth waves")
        N = len(wave.eta) - 1
        J = arange(1, N + 1)
        B = wave.data["B"]
        k = wave.k

        Jk = J * k
        facs = B[1:] / cosh(Jk * wave.depth)

        # Repr of np.float64(42.0) is "np.float64(42.0)" and not "42.0"
        # We use repr to make Python output a "smart" amount of digits
        c = np2py(wave.c)
        Jk = np2py(Jk)
        facs = np2py(facs)

        cpp = " + ".join(
            f"{facs[i]!r} * cos({Jk[i]!r} * (x[0] - {c!r}* t)) * sinh({Jk[i]!r} * x[2])"
            for i in range(N)
        )

        if frame == "b":
            return f"{np2py(B[0])!r} * x[2] + {cpp}"
        elif frame == "c":
            return cpp

    def elevation(self):
        """
        Return C++ code for evaluating the elevation of this specific wave.
        The positive traveling direction is x[0]
        """
        wave = self.wave
        if wave.depth < 0:
            raise RaschiiError("Cannot currently generate C++ code for infinite depth waves")
        N = wave.E.size - 1
        facs = wave.E * 2 / N
        facs[0] *= 0.5
        facs[-1] *= 0.5

        # Repr of np.float64(42.0) is "np.float64(42.0)" and not "42.0"
        # We use repr to make Python output a "smart" amount of digits
        k = np2py(wave.k)
        c = np2py(wave.c)
        facs = np2py(facs)

        code = " + ".join(
            f"{facs[j]!r} * cos({j:d} * {k!r} * (x[0] - {c!r} * t))" for j in range(0, N + 1)
        )
        return code

    def slope(self):
        """
        Return C++ code for evaluating the surface slope of this specific wave.
        The positive traveling direction is x[0]
        """
        wave = self.wave
        if wave.depth < 0:
            raise RaschiiError("Cannot currently generate C++ code for infinite depth waves")
        N = wave.E.size - 1
        facs = wave.E * 2 / N * wave.k * -1.0
        facs[0] *= 0.5
        facs[-1] *= 0.5

        # Repr of np.float64(42.0) is "np.float64(42.0)" and not "42.0"
        # We use repr to make Python output a "smart" amount of digits
        k = np2py(wave.k)
        c = np2py(wave.c)
        facs = np2py(facs)

        code = " + ".join(
            f"{facs[j]!r} * {j:d} * sin({j:d} * {k!r} * (x[0] - {c!r} * t))"
            for j in range(0, N + 1)
        )
        return code

    def velocity(self, all_points_wet=False):
        """
        Return C++ code for evaluating the particle velocities of this specific
        wave. Returns the x and z components only with z positive upwards. The
        positive traveling direction is x[0] and the vertical coordinate is x[2]
        which is zero at the bottom and equal to +depth at the mean water level.
        """
        wave = self.wave
        if wave.depth < 0:
            raise RaschiiError("Cannot currently generate C++ code for infinite depth waves")
        N = len(wave.eta) - 1
        J = arange(1, N + 1)
        B = wave.data["B"]
        k = wave.k

        Jk = J * k
        facs = J * B[1:] * k / cosh(Jk * wave.depth)

        # Repr of np.float64(42.0) is "np.float64(42.0)" and not "42.0"
        # We use repr to make Python output a "smart" amount of digits
        c = np2py(wave.c)
        Jk = np2py(Jk)
        facs = np2py(facs)

        cpp_x = " + ".join(
            f"{facs[i]!r} * cos({Jk[i]!r} * (x[0] - {c!r} * t)) * cosh({Jk[i]!r} * x[2])"
            for i in range(N)
        )
        cpp_z = " + ".join(
            f"{facs[i]!r} * sin({Jk[i]!r} * (x[0] - {c!r} * t)) * sinh({Jk[i]!r} * x[2])"
            for i in range(N)
        )

        if all_points_wet:
            return cpp_x, cpp_z

        # Handle velocities above the free surface
        e_cpp = self.elevation()
        cpp_ax = cpp_az = None
        cpp_psiw = cpp_psia = cpp_slope = None
        if wave.air is not None:
            cpp_ax, cpp_az = wave.air.cpp.velocity()
            cpp_psiw = self.stream_function(frame="c")
            cpp_psia = wave.air.cpp.stream_function(frame="c")
            cpp_slope = self.slope()

        cpp_x = blend_air_and_wave_velocity_cpp(
            cpp_x, cpp_ax, e_cpp, "x", wave.eta_eps, wave.air, cpp_psiw, cpp_psia, cpp_slope
        )
        cpp_z = blend_air_and_wave_velocity_cpp(
            cpp_z, cpp_az, e_cpp, "z", wave.eta_eps, wave.air, cpp_psiw, cpp_psia, cpp_slope
        )

        return cpp_x, cpp_z
