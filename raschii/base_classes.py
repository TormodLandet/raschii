import numpy as np
from numpy.typing import NDArray


class WaveModel:
    """
    Base class for Raschii wave models.
    """

    def __init__(self, height: float, depth: float, length: float, g: float = 9.81):
        self.height: float = height  #: The wave height
        self.depth: float = depth  #: The water depth
        self.length: float = length  #: The wave length
        self.g: float = g  #: The acceleration of gravity

        self.T: float  #: The wave period [t], to be defined in subclasses
        self.omega: float  #: The wave angular frequency [rad/s], to be defined in subclasses
        self.k: float  #: The wave number [1/m], to be defined in subclasses
        self.c: float  #: The wave celerity [m/s], to be defined in subclasses

    def surface_elevation(
        self,
        x: float | list[float] | NDArray,
        t: float | list[float] | NDArray = 0.0,
        include_depth: bool = True,
    ) -> NDArray | float:
        """
        Compute the surface elevation at time t for position(s) x.

        Returns a scalar when both x and t are scalar, an ndarray otherwise.
        Output shape follows the NumPy convention (scalar in → scalar out):

        - scalar x, scalar t  → scalar
        - array x (N), scalar t → (N,)
        - scalar x, array t (T) → (T,)
        - array x (N), array t (T) → (T, N)
        """
        x_was_scalar = np.asarray(x).ndim == 0
        t_was_scalar = np.asarray(t).ndim == 0
        x_arr = np.atleast_1d(np.asarray(x, dtype=float))
        t_arr = np.atleast_1d(np.asarray(t, dtype=float))
        result = self._surface_elevation(x_arr, t_arr, include_depth)  # (T, N)
        if x_was_scalar and t_was_scalar:
            return result[0, 0]
        elif t_was_scalar:
            return result[0]
        elif x_was_scalar:
            return result[:, 0]
        return result

    def velocity(
        self,
        x: float | list[float] | NDArray,
        z: float | list[float] | NDArray,
        t: float | list[float] | NDArray = 0.0,
        all_points_wet: bool = False,
    ) -> NDArray:
        """
        Compute the fluid velocity at time t for position(s) (x, z)
        where z is 0 at the bottom and equal to depth at the free surface.

        Output shape:

        - (2,) if x, z, and t are scalar
        - (N, 2) if x/z are arrays and t is scalar
        - (T, 2) if x/z are scalar and t is an array
        - (T, N, 2) if x/z are arrays and t is an array
        """
        point_is_scalar = np.asarray(x).ndim == 0 and np.asarray(z).ndim == 0
        time_is_scalar = np.asarray(t).ndim == 0
        x_arr = np.atleast_1d(np.asarray(x, dtype=float))
        z_arr = np.atleast_1d(np.asarray(z, dtype=float))
        t_arr = np.atleast_1d(np.asarray(t, dtype=float))
        x_arr, z_arr = np.broadcast_arrays(x_arr, z_arr)
        result = self._velocity(x_arr, z_arr, t_arr, all_points_wet)  # (T, N, 2)
        if time_is_scalar and point_is_scalar:
            return result.squeeze()
        elif time_is_scalar:
            return result[0]
        elif point_is_scalar:
            return result[:, 0]
        return result

    def velocity_potential(
        self,
        x: float | list[float] | NDArray,
        z: float | list[float] | NDArray,
        t: float | list[float] | NDArray = 0.0,
    ) -> NDArray | float:
        """
        Compute the earth-frame velocity potential φ at time t for position(s) (x, z).

        z is measured from the sea floor (z=0 at bottom, z≈depth at calm surface).
        The gradient of φ equals the oscillatory fluid velocity as returned by
        :meth:`velocity`; the mean-flow current term is excluded.

        Returns a scalar when both x/z and t are scalar, an ndarray otherwise.
        Output shape:

        - scalar x/z, scalar t → scalar
        - array x/z (N), scalar t → (N,)
        - scalar x/z, array t (T) → (T,)
        - array x/z (N), array t (T) → (T, N)
        """
        point_is_scalar = np.asarray(x).ndim == 0 and np.asarray(z).ndim == 0
        t_was_scalar = np.asarray(t).ndim == 0
        x_arr = np.atleast_1d(np.asarray(x, dtype=float))
        z_arr = np.atleast_1d(np.asarray(z, dtype=float))
        t_arr = np.atleast_1d(np.asarray(t, dtype=float))
        x_arr, z_arr = np.broadcast_arrays(x_arr, z_arr)
        result = self._velocity_potential(x_arr, z_arr, t_arr)  # (T, N)
        if point_is_scalar and t_was_scalar:
            return result[0, 0]
        elif t_was_scalar:
            return result[0]
        elif point_is_scalar:
            return result[:, 0]
        return result

    def _surface_elevation(self, x: NDArray, t: NDArray, include_depth: bool) -> NDArray:
        """Compute surface elevation. x: (N,), t: (T,) → returns (T, N)."""
        raise NotImplementedError("This method should be implemented in subclasses.")

    def _velocity(self, x: NDArray, z: NDArray, t: NDArray, all_points_wet: bool) -> NDArray:
        """Compute fluid velocity. x, z: (N,), t: (T,) → returns (T, N, 2)."""
        raise NotImplementedError("This method should be implemented in subclasses.")

    def _velocity_potential(self, x: NDArray, z: NDArray, t: NDArray) -> NDArray:
        """Compute velocity potential. x, z: (N,), t: (T,) → returns (T, N)."""
        raise NotImplementedError("This method should be implemented in subclasses.")
