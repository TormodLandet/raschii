class FentonWave:
    def __init__(self, height, depth, length, order):
        """
        Implement stream function waves based on the paper by Rienecker and
        Fenton (1981)
        
        * height: wave height above still water level
        * depth: still water distance from the flat sea bottom to the surface
        * length: the periodic length of the wave (distance between peaks)
        * order: the number of coefficients in the truncated Fourier series
        """
        self.height = height
        self.depth = depth
        self.length = length
        self.order = order
        self.coefficients = fenton_coefficients(height, depth, length, order)
    
    def surface_elevation(self, x):
        return None
    
    def velocity(self, x, z):
        return None


def fenton_coefficients(height, depth, length, order):
    return None
