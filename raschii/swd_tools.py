import os
from datetime import datetime
from struct import pack

import numpy as np


class SwdShape2:
    required_input = {"t_wave", "length", "depth", "c_cofs", "h_cofs", "input_data"}
    optional_input = {"g": 9.81, "order_zpos": -1}

    def __init__(self, t_wave, length, depth, c_cofs, h_cofs, input_data, g=9.81, order_zpos=-1):
        """
        For SWD output of stationary 2d-waves in finite depth.
        (shape=2 is a SWD class using constant spaced wave numbers in finite depth)

        * t_wave: wave period [sec]
        * length: the periodic length of the wave (distance between peaks) [m]
        * depth: distance from the flat sea bottom to the calm free surface [m]
        * h_cofs: the complex sequence of SWD elevation coefficients including
                  the bias term (j=0) is = h_cofs[j] * exp(I j * omega * t),  j=0,1,...
        * c_cofs: the complex sequence of SWD potential coefficients including
                  the bias term (j=0) is = c_cofs[j] * exp(I j * omega * t),  j=0,1,...
        * input_data: dictionary of input variables to reconstruct the waves (stored in SWD file)
        * order_zpos: default perturbation order when evaluating fields above calm surface
                      (-1 if fully nonlinear exponential terms apply above z=0)
                      (this value can be replaced by the SWD application constructor.)
        """

        assert len(c_cofs) == len(h_cofs)
        self.n_swd = len(c_cofs) - 1

        self.t_wave = t_wave
        self.omega = 2.0 * np.pi / t_wave
        self.kwave = 2.0 * np.pi / length
        self.depth = depth
        self.order_zpos = order_zpos
        self.c_cofs = c_cofs
        self.h_cofs = h_cofs
        self.g = g
        self.lscale = 1.0
        self.input_data = str(input_data).encode("utf-8")

    def write(self, path, dt, tmax=None, nperiods=None):
        """
        Write the actual SWD file

        * path: Name of SWD file
        * dt: Sampling spacing of time in SWD file
        * tmax: Approximate length of time series to be stored in the SWD file
        * nperiods: Alternative specification of tmax in terms of number of oscillation periods.
        """
        if tmax is None:
            assert nperiods is not None
            tmax = nperiods * self.t_wave
        assert tmax > dt > 0.0

        # Get the raschii version number
        version = "x.y.z"
        for line in open(os.path.join(os.path.dirname(__file__), "__init__.py"), encoding="utf-8"):
            if line.startswith("__version__"):
                version = line.split("=")[1].strip()[1:-1]
        prog = "raschii-" + version

        h_swd = np.empty(self.n_swd + 1, complex)
        ht_swd = np.empty(self.n_swd + 1, complex)
        c_swd = np.empty(self.n_swd + 1, complex)
        ct_swd = np.empty(self.n_swd + 1, complex)
        dtime = datetime.now().strftime("%Y:%m:%d %H:%M:%S")
        nsteps = int((tmax + 0.0001 * dt) / dt + 1)

        out = open(path, "wb")
        out.write(pack("<f", 37.0221))  # Magic number
        out.write(pack("<i", 100))  # fmt
        out.write(pack("<i", 2))  # shp for long-crested constant finite depth
        out.write(pack("<i", 1))  # amp (both elevation and field data)
        out.write(pack("<30s", prog.ljust(30).encode("utf-8")))  # prog name
        out.write(pack("<20s", dtime.ljust(20).encode("utf-8")))  # date
        nid = len(self.input_data)
        out.write(pack("<i", nid))  # length of input file
        out.write(pack("<{0}s".format(nid), self.input_data))
        out.write(pack("<f", self.g))  # acc. of gravity
        out.write(pack("<f", self.lscale))
        out.write(pack("<i", 0))  # nstrip
        out.write(pack("<i", nsteps))
        out.write(pack("<f", dt))
        out.write(pack("<i", self.order_zpos))
        out.write(pack("<i", self.n_swd))
        out.write(pack("<f", self.kwave))  # delta_k in swd
        out.write(pack("<f", self.depth))

        def dump_cofs(vals):
            for j in range(self.n_swd + 1):
                r = vals[j]
                out.write(pack("<f", r.real))
                out.write(pack("<f", r.imag))

        for istep in range(nsteps):
            t = istep * dt
            for j in range(self.n_swd + 1):
                fac = complex(0.0, j * self.omega)
                h_swd[j] = self.h_cofs[j] * np.exp(fac * t)
                ht_swd[j] = fac * h_swd[j]
                c_swd[j] = self.c_cofs[j] * np.exp(fac * t)
                ct_swd[j] = fac * c_swd[j]
            dump_cofs(h_swd)
            dump_cofs(ht_swd)
            dump_cofs(c_swd)
            dump_cofs(ct_swd)
        out.close()
