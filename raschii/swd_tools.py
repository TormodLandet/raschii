from datetime import datetime
from struct import pack, unpack

import numpy as np

from raschii import version as raschii_version


class SwdShape1and2:
    def __init__(
        self,
        t_wave: float,
        length: float,
        depth: float,
        c_cofs: list[complex],
        h_cofs: list[complex],
        input_data: str,
        g: float = 9.81,
        order_zpos: int = -1,
    ):
        """
        For SWD output of stationary 2d-waves in finite depth.

        SWD shape=1:
            SWD class using constant spaced wave numbers in infinite finite depth

        SWD shape=2:
            SWD class using constant spaced wave numbers in finite depth

        Parameters:

        * t_wave: wave period [sec]
        * length: the periodic length of the wave (distance between peaks) [m]
        * depth: distance from the flat sea bottom to the calm free surface [m]
          or -1.0 for infinite depth
        * h_cofs: the complex sequence of SWD elevation coefficients including
                  the bias term (j=0) is = h_cofs[j] * exp(I j * omega * t),  j=0,1,...
        * c_cofs: the complex sequence of SWD potential coefficients including
                  the bias term (j=0) is = c_cofs[j] * exp(I j * omega * t),  j=0,1,...
        * input_data: dictionary of input variables to reconstruct the waves (stored in SWD file)
        * order_zpos: default perturbation order when evaluating fields above calm surface
                      (-1 if fully nonlinear exponential terms apply above z=0)
                      (this value can be replaced by the SWD application constructor.)

        See:

        * https://spectral-wave-data.readthedocs.io/en/latest/shape_1.html
        * https://spectral-wave-data.readthedocs.io/en/latest/shape_2.html
        * https://spectral-wave-data.readthedocs.io/en/latest/swd_format.html
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

        See: https://spectral-wave-data.readthedocs.io/en/latest/swd_format.html
        """
        if tmax is None:
            assert nperiods is not None
            tmax = nperiods * self.t_wave
        assert tmax > dt > 0.0

        prog = f"raschii-{raschii_version}"
        dtime = datetime.now().strftime("%Y:%m:%d %H:%M:%S")
        nsteps = int((tmax + 0.0001 * dt) / dt + 1)

        with open(path, "wb") as out:
            out.write(pack("<f", 37.0221))  # Magic number
            out.write(pack("<i", 100))  # fmt
            if self.depth < 0:
                out.write(pack("<i", 1))  # shp for long-crested infinite depth
            else:
                out.write(pack("<i", 2))  # shp for long-crested constant finite depth
            out.write(pack("<i", 1))  # amp (both elevation and field data)
            out.write(pack("<30s", prog.ljust(30).encode("utf-8")))  # prog name
            out.write(pack("<20s", dtime.ljust(20).encode("utf-8")))  # date
            nid = len(self.input_data)
            out.write(pack("<i", nid))  # length of input file
            out.write(pack(f"<{nid}s", self.input_data))
            out.write(pack("<f", self.g))  # acc. of gravity
            out.write(pack("<f", self.lscale))
            out.write(pack("<i", 0))  # nstrip
            out.write(pack("<i", nsteps))
            out.write(pack("<f", dt))
            out.write(pack("<i", self.order_zpos))
            out.write(pack("<i", self.n_swd))
            out.write(pack("<f", self.kwave))  # delta_k in swd
            if self.depth >= 0:
                out.write(pack("<f", self.depth))

            # Wrte the h, ht, c, and ct coefficient arrays for each time step
            N2 = self.n_swd + 1
            fac = np.array([1j * j * self.omega for j in range(N2)], dtype=np.complex64)
            h_swd = np.empty(N2, np.complex64)
            ht_swd = np.empty(N2, np.complex64)
            c_swd = np.empty(N2, np.complex64)
            ct_swd = np.empty(N2, np.complex64)
            for istep in range(nsteps):
                t = istep * dt
                fac2 = np.exp(fac * t)

                h_swd[:] = self.h_cofs * fac2
                ht_swd[:] = h_swd * fac
                c_swd[:] = self.c_cofs * fac2
                ct_swd[:] = c_swd * fac
                
                out.write(h_swd.tobytes())
                out.write(ht_swd.tobytes())
                out.write(c_swd.tobytes())
                out.write(ct_swd.tobytes())


class SwdReaderForRaschiiTests:
    def __init__(self, swd_file: str):
        """
        Read an SWD file of shape 1 or 2.

        Only used in tests to verify the SWD files we write from Raschii.

        See: https://spectral-wave-data.readthedocs.io/en/latest/swd_format.html
        """
        with open(swd_file, "rb") as f:
            magic, fmt, shp, amp = unpack("<fiii", f.read(4 * 4))
            assert abs(magic - 37.0221) < 1e-4, "Invalid SWD file magic number"
            assert fmt == 100, "Unsupported SWD file format"

            self.shp: int = shp
            self.amp: int = amp

            prog, date, nid = unpack("<30s20si", f.read(30 + 20 + 4))
            self.prog: str = prog.decode("utf-8")
            self.date: str = date.decode("utf-8")
            self.input_data: str = unpack(f"<{nid}s", f.read(nid))[0].decode("utf-8")

            self.g: float = unpack("<f", f.read(4))[0]  # acceleration of gravity
            self.lscale: float = unpack("<f", f.read(4))[0]  # length scale
            self.nstrip: int = unpack("<i", f.read(4))[
                0
            ]  # number of initial time steps stripped off
            self.nsteps: int = unpack("<i", f.read(4))[0]  # number of time steps
            self.dt: float = unpack("<f", f.read(4))[0]  # time step
            self.order: int = unpack("<i", f.read(4))[
                0
            ]  # order of wave kinematics above the still water level
            self.nx: int = unpack("<i", f.read(4))[0]  # number of spectral components
            self.dk: float = unpack("<f", f.read(4))[0]  # wave number spacing
            if self.shp == 1:
                self.depth: float = -1.0
            elif self.shp == 2:
                self.depth: float = unpack("<f", f.read(4))[0]
            else:
                raise ValueError(f"Unsupported SWD shape: {self.shp}")

            # The timesteps
            self.t_vector = np.arange(self.nsteps) * self.dt

            # The wave numbers
            nsig = self.nx + 1
            self.k_vector = np.arange(nsig) * self.dk

            # Read the h(t, k) array by looping through the time steps
            self.h_vectors = np.empty((self.nsteps, nsig), dtype=np.complex64)
            for it in range(self.nsteps):
                self.h_vectors[it, :] = np.fromfile(f, dtype=np.complex64, count=nsig)
                pos = f.tell()
                spool = 8 * nsig  # spool past ht(t,k)
                if self.amp < 3:
                    spool += 8 * nsig  # spool past c(t,k)
                    spool += 8 * nsig  # spool past ct(t,k)
                f.seek(pos + spool)

    def surface_elevation(self, x: float = 0.0) -> np.ndarray:
        """
        Compute the SWD surface elevation at position (x, y) and time t.
        For more details see the documentation at

        Shape 1:
            https://spectral-wave-data.readthedocs.io/en/latest/shape_1.html

        Shape 2:
            https://spectral-wave-data.readthedocs.io/en/latest/shape_2.html
        """
        X = np.exp(-1j * self.k_vector[1:] * x)
        h = self.h_vectors[:, 1:]
        return np.real(h @ X)
