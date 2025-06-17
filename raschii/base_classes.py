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
        self.c: float  #: The wave celery [m/s], to be defined in subclasses

    def surface_elevation(self, x: float | list[float], t: float = 0.0, include_depth: bool = True):
        """
        Compute the surface elavation at time t for position(s) x
        """
        raise NotImplementedError("This method should be implemented in subclasses.")

    def velocity(
        self,
        x: float | list[float],
        z: float | list[float],
        t: float = 0,
        all_points_wet: bool = False,
    ):
        """
        Compute the fluid velocity at time t for position(s) (x, z)
        where z is 0 at the bottom and equal to depth at the free surface
        """
        raise NotImplementedError("This method should be implemented in subclasses.")
