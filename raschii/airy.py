from numpy import pi, cos, sin, zeros, array, asarray, sinh, cosh, tanh
from .common import blend_air_and_wave_velocities, blend_air_and_wave_velocity_cpp


class AiryWave:
    required_input = ("height", "depth", "length")
    optional_input = {"air": None, "g": 9.81}

    def __init__(self, height, depth, length, air=None, g=9.81):
        """
        Linear Airy waves

        * height: wave height above still water level
        * depth: still water distance from the flat sea bottom to the surface
        * length: the periodic length of the wave (distance between peaks)
        """
        self.height = height
        self.depth = depth
        self.length = length
        self.air = air
        self.g = g
        self.order = 1
        self.warnings = ""

        self.k = 2 * pi / length
        self.omega = (self.k * g * tanh(self.k * depth)) ** 0.5
        self.c = self.omega / self.k
        self.T = self.length / self.c  # Wave period

        # For evaluating velocities close to the free surface
        self.eta_eps = self.height / 1e5

        # Provide velocities also in the air phase
        if self.air is not None:
            self.air.set_wave(self)

    def surface_elevation(self, x, t=0):
        """
        Compute the surface elavation at time t for position(s) x
        """
        if isinstance(x, (float, int)):
            x = array([x], float)
        x = asarray(x)
        return self.depth + self.height / 2 * cos(self.k * x - self.omega * t)

    def velocity(self, x, z, t=0, all_points_wet=False):
        """
        Compute the fluid velocity at time t for position(s) (x, z)
        where z is 0 at the bottom and equal to depth at the free surface
        """
        if isinstance(x, (float, int)):
            x, z = [x], [z]
        x = asarray(x, dtype=float)
        z = asarray(z, dtype=float)

        H = self.height
        k = self.k
        d = self.depth
        w = self.omega

        vel = zeros((x.size, 2), float) + 1
        vel[:, 0] = w * H / 2 * cosh(k * z) / sinh(k * d) * cos(k * x - w * t)
        vel[:, 1] = w * H / 2 * sinh(k * z) / sinh(k * d) * sin(k * x - w * t)

        if not all_points_wet:
            blend_air_and_wave_velocities(x, z, t, self, self.air, vel, self.eta_eps)

        return vel

    def elevation_cpp(self):
        """
        Return C++ code for evaluating the elevation of this specific wave.
        The positive traveling direction is x[0]
        """
        # Repr of np.float64(42.0) is "np.float64(42.0)" and not "42.0"
        # We use repr to make Python output a "smart" amount of digits
        depth = float(self.depth)
        height = float(self.height)
        k = float(self.k)
        c = float(self.c)

        return f"{depth!r} + {height!r} / 2.0 * cos({k!r} * (x[0] - {c!r} * t))"

    def velocity_cpp(self, all_points_wet=False):
        """
        Return C++ code for evaluating the particle velocities of this specific
        wave. Returns the x and z components only with z positive upwards. The
        positive traveling direction is x[0] and the vertical coordinate is x[2]
        which is zero at the bottom and equal to +depth at the mean water level.
        """
        # Repr of np.float64(42.0) is "np.float64(42.0)" and not "42.0"
        # We use repr to make Python output a "smart" amount of digits
        H = float(self.height)
        k = float(self.k)
        d = float(self.depth)
        w = float(self.omega)
        a = float(w * H / (2 * sinh(k * d)))

        cpp_x = f"{a!r} * cosh({k!r} * x[2]) * cos({k!r} * x[0] - {w!r} * t)"
        cpp_z = f"{a!r} * sinh({k!r} * x[2]) * sin({k!r} * x[0] - {w!r} * t)"

        if all_points_wet:
            return cpp_x, cpp_z

        # Handle velocities above the free surface
        e_cpp = self.elevation_cpp()
        cpp_ax = cpp_az = None
        cpp_psiw = cpp_psia = cpp_slope = None
        if self.air is not None:
            cpp_ax, cpp_az = self.air.velocity_cpp()
            cpp_psiw = self.stream_function_cpp(frame="c")
            cpp_psia = self.air.stream_function_cpp(frame="c")
            cpp_slope = self.slope_cpp()

        cpp_x = blend_air_and_wave_velocity_cpp(
            cpp_x, cpp_ax, e_cpp, "x", self.eta_eps, self.air, cpp_psiw, cpp_psia, cpp_slope
        )
        cpp_z = blend_air_and_wave_velocity_cpp(
            cpp_z, cpp_az, e_cpp, "z", self.eta_eps, self.air, cpp_psiw, cpp_psia, cpp_slope
        )

        return cpp_x, cpp_z

    def write_swd(self, path, dt, tmax=None, nperiods=None):
        """
        NOT IMPLEMENTED FOR AIRY!

        Write a SWD-file of the wave field according to the file
        specification in the Github repository spectral-wave-data ....

        * path:     Full path of the new SWD file
        * dt:       The temporal sampling spacing in the SWD file
        * tmax:     The temporal sampling range in the SWD file is [0, tmax]
        * nperiods: Alternative specification: tmax = nperiods * wave_period
        """

        print("Airy does not currently support writing SWD files.")
        print("Please use the Stokes wave model with order N=1")
        return NotImplemented
