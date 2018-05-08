from numpy import pi, sinh, cosh, tanh, cos, sin, zeros, arange, trapz, isfinite
from numpy.linalg import solve
from . import NonConvergenceError


class FentonWave:
    def __init__(self, height, depth, length, N, g=9.81):
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
        self.N = N
        
        self.coeffs = fenton_coefficients(height, depth, length, N, g)
        self.B = self.coeffs['B']
        self.Q = self.coeffs['Q']
        self.R = self.coeffs['R']
        self.eta = self.coeffs['eta'] * depth
        self.x = self.coeffs['x'] * depth
        self.k = self.coeffs['k'] / depth
    
    @property
    def coefficients(self):
        return (self.B, self.Q, self.R)
    
    def surface_elevation(self, x):
        return None
    
    def velocity(self, x, z):
        return None


def fenton_coefficients(height, depth, length, N, g):
    """
    Find B, Q and R by Newton-Raphson following Rienecker and Fenton (1981)
    """
    # Non dimensionalised input
    H = height / depth
    lam = length / depth
    k = 2 * pi / lam
    c = (tanh(k) / k) ** 0.5
    D = 1
    N_unknowns = 2 * N + 5
    
    # Input data arrays
    J = arange(0, N)
    M = arange(0, N + 1)
    x = M * lam / (2 * N)
    
    # Initial guesses for the unknowns (linear wave)
    B = zeros(N + 1, float)
    B[0] = -c
    B[1] = -H / (4 * c * k)
    eta = 1 + H / 2 * cos(k * x)
    Q = c
    R = 1 + 0.5 * c**2
    
    def optimize(B, Q, R, eta, k, maxiter=100, tol=1e-6):
        """
        Find B, Q and R by Newton iterations starting from the given initial
        guesses. According to Rienecker and Fenton (1981) a linear theory
        initial guess should work unless H close to breaking, then an initial
        guess from the optimization routine run with a slightly lower H should
        be used instead
        """
        # Initial guess
        coeffs = zeros(N_unknowns, float)
        coeffs[:N + 1] = B
        coeffs[N + 1:2 * N + 2] = eta
        coeffs[2 * N + 2] = Q
        coeffs[2 * N + 3] = R
        coeffs[2 * N + 4] = k
        f = func(coeffs)
        
        for _ in range(maxiter):
            jac = fprime(coeffs)
            delta = solve(jac, -f)
            coeffs += delta
            f = func(coeffs)
            
            assert all(isfinite(f))
            error = abs(f).max()
            if error < tol:
                return B, eta, Q, R, k, error
        raise NonConvergenceError('Optimization did not converge, err = %r'
                                  % error)
    
    def func(coeffs):
        "The function to minimize"
        B0 = coeffs[0]
        B = coeffs[1:N + 1]
        eta = coeffs[N + 1:2 * N + 2]
        Q = coeffs[2 * N + 2]
        R = coeffs[2 * N + 3]
        k = coeffs[2 * N + 4]
        
        f = zeros(N_unknowns, float)
        for m in range(N + 1):
            S = sinh(J * k * eta[m])
            C = cosh(J * k * eta[m])
            c = cos(J * m * pi / N)
            s = sin(J * m * pi / N)
            CD = cosh(J * k * D)
            
            # Velocity at the free surface
            um = B0 + k * J.dot(B * C / CD * c)
            vm = 0 + k * J.dot(B * S / CD * s)
            
            # Enforce a streamline along the free surface
            f[m] = B0 * eta[m] + B.dot(S / CD * c + Q)
            
            # Enforce the dynamic free surface boundary condition
            f[N + 1 + m] = um**2 + vm**2 + eta[m] - R
            
        # Enforce mean(eta) = D
        f[2 * N + 2] = 1 / (2 * N) * trapz(eta) - 1
            
        # Enforce eta_0 - eta_N = H
        f[2 * N + 3] = eta[0] - eta[-1] - H
        
        # Enforce k = 2*pi/lambda
        f[2 * N + 4] = k - 2 * pi / lam
        
        return f
    
    def fprime(coeffs):
        "The Jacobian of the function to minimize"
        dc = 1e-10
        jac = zeros((N_unknowns, N_unknowns), float)
        for i in range(N_unknowns):
            cpdc = coeffs.copy()
            cpdc[i] += dc
            jac[i] = (func(cpdc) - func(coeffs)) / dc
        return jac
    
    # Perform optimization
    B, eta, Q, R, k, error = optimize(B, Q, R, eta, k)
    
    return dict(B=B, eta=eta, Q=Q, R=R, k=k, c=c, error=error, x=x)
