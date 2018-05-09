from numpy import pi, sinh, cosh, tanh, cos, sin, zeros, arange, trapz, isfinite
from numpy.linalg import solve, cond
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
        
        self.coeffs = fenton_coefficients(height, depth, length, N)
        self.B = self.coeffs['B']
        self.Q = self.coeffs['Q']
        self.R = self.coeffs['R']
        self.eta = self.coeffs['eta'] * depth
        self.x = self.coeffs['x'] * depth
        self.k = self.coeffs['k'] / depth
        self.c = self.coeffs['B'][0] * (g * depth)**0.5
        self.T = length / self.c
    
    @property
    def coefficients(self):
        return (self.B, self.Q, self.R)
    
    def surface_elevation(self, x):
        return None
    
    def velocity(self, x, z):
        return None


def fenton_coefficients(height, depth, length, N, relax=1.0):
    """
    Find B, Q and R by Newton-Raphson following Rienecker and Fenton (1981)
    """
    # Non dimensionalised input
    H = height / depth
    lam = length / depth
    k = 2 * pi / lam
    c = (tanh(k) / k)**0.5
    D = 1
    N_unknowns = 2 * (N + 1) + 2
    
    # Input data arrays
    J = arange(1, N + 1)
    M = arange(0, N + 1)
    x = M * lam / (2 * N)
    
    # Initial guesses for the unknowns (linear wave)
    B = zeros(N + 1, float)
    B[0] = -c
    B[1] = -H / (4 * c * k)
    eta = 1 + H / 2 * cos(k * x)
    Q = c
    R = 1 + 0.5 * c**2
    
    def optimize(B, Q, R, eta, maxiter=100, tol=1e-6, relax=1.0):
        """
        Find B, Q and R by Newton iterations starting from the given initial
        guesses. According to Rienecker and Fenton (1981) a linear theory
        initial guess should work unless H close to breaking, then an initial
        guess from the optimization routine run with a slightly lower H should
        be used instead
        """
        # Insert initial guesses into coefficient vector
        coeffs = zeros(N_unknowns, float)
        coeffs[:N + 1] = B
        coeffs[N + 1:2 * N + 2] = eta
        coeffs[2 * N + 2] = Q
        coeffs[2 * N + 3] = R
        f = func(coeffs, H, k, D, J, M)
        
        allc = [coeffs.copy()]
        
        def pl():
            from matplotlib import pyplot
            pyplot.clf()
            for i, ci in enumerate(allc):
                pyplot.plot(ci[N + 1:2 * N + 2], label=str(i))
                pyplot.plot(ci[:N + 1], ls=':', label=str(i) + 'b')
            pyplot.legend()
            pyplot.show()
        
        for i in range(maxiter):
            jac = fprime(coeffs, H, k, D, J, M)
            print('Iteration %2d has condition number %.3e' % (i + 1, cond(jac)))
            delta = solve(jac, -f)
            coeffs += delta * relax
            f = func(coeffs, H, k, D, J, M)
            
            allc.append(coeffs.copy())
            # DEBUG: pl()
            
            error = abs(f).max()
            if not isfinite(error):
                raise NonConvergenceError('Optimization did not converge. Got '
                                          '%r in iteration %d' % (error, i + 1))
            elif error < tol:
                B = coeffs[:N + 1]
                eta = coeffs[N + 1:2 * N + 2]
                Q = coeffs[2 * N + 2]
                R = coeffs[2 * N + 3]
                return B, eta, Q, R, error
        raise NonConvergenceError('Optimization did not converge after %d'
                                  'iterations, error = %r' % (i + 1, error))
    
    # Perform optimization
    B, eta, Q, R, error = optimize(B, Q, R, eta, relax=relax)
    
    return dict(B=B, eta=eta, Q=Q, R=R, k=k, c=c, error=error, x=x)


def func(coeffs, H, k, D, J, M):
    "The function to minimize"
    N_unknowns = coeffs.size
    N = (coeffs.size - 4) // 2
    assert N == 10
    
    B0 = coeffs[0]
    B = coeffs[1:N + 1]
    eta = coeffs[N + 1:2 * N + 2]
    Q = coeffs[2 * N + 2]
    R = coeffs[2 * N + 3]
    
    # The function to me minimized
    f = zeros(N_unknowns, float)
    
    # Loop over the N + 1 points along the half wave
    for m in M:
        S1 = sinh(J * k * eta[m])
        C1 = cosh(J * k * eta[m])
        S2 = sin(J * m * pi / N)
        C2 = cos(J * m * pi / N)
        CD = cosh(J * k * D)
        
        # Velocity at the free surface
        um = B0 + k * J.dot(B * C1 / CD * C2)
        vm = 0 + k * J.dot(B * S1 / CD * S2)
        
        # Enforce a streamline along the free surface
        f[m] = B0 * eta[m] + B.dot(S1 / CD * C2) + Q
        
        # Enforce the dynamic free surface boundary condition
        f[N + 1 + m] = (um**2 + vm**2) / 2 + eta[m] - R
        
    # Enforce mean(eta) = D
    f[-2] = trapz(eta) / N - 1
        
    # Enforce eta_0 - eta_N = H, the wave height criterion
    f[-1] = eta[0] - eta[-1] - H
    
    return f


def fprime_num(coeffs, H, k, D, J, M):
    "The Jacobian of the function to minimize (numerical version)"
    N_unknowns = coeffs.size
    N = (coeffs.size - 4) // 2
    assert N == 10
    
    dc = 1e-10
    jac = zeros((N_unknowns, N_unknowns), float)
    f0 = func(coeffs, H, k, D, J, M)
    for i in range(N_unknowns):
        incr = zeros(N_unknowns, float)
        incr[i] = dc
        f1 = func(coeffs + incr, H, k, D, J, M)
        jac[i] = (f1 - f0) / dc
    return jac


def fprime(coeffs, H, k, D, J, M):
    "The Jacobian of the function to minimize"
    N_unknowns = coeffs.size
    N = (coeffs.size - 4) // 2
    assert N == 10
    
    jac = zeros((N_unknowns, N_unknowns), float)
    B0 = coeffs[0]
    B = coeffs[1:N + 1]
    eta = coeffs[N + 1:2 * N + 2]
    
    for m in range(N + 1):
        S1 = sinh(J * k * eta[m])
        C1 = cosh(J * k * eta[m])
        S2 = sin(J * m * pi / N)
        C2 = cos(J * m * pi / N)
        CD = cosh(J * k * D)
        
        SC = S1 / CD * C2
        SS = S1 / CD * S2
        CC = C1 / CD * C2
        CS = C1 / CD * S2
        
        # Velocity at the free surface
        um = B0 + k * J.dot(B * CC)
        vm = 0 + k * J.dot(B * SS)
        
        # Derivatives of the eq. for the streamline along the free surface
        jac[N + 1 + m, m] = um
        jac[0, 0:N + 1] = eta  # FIXME: this is different sign than in the paper <--------------------------
        jac[1:N + 1, m] = SC
        jac[-2, m] = 1
        
        # Derivatives of the dynamic free surface boundary condition
        jac[N + 1 + m, N + 1 + m] = 1 + (um * k**2 * B.dot(J**2 * SC) +
                                         vm * k**2 * B.dot(J**2 * CS))
        jac[-1, N + 1 + m] = -1
        jac[0, N + 1 + m] = um  # FIXME: this is different sign than in the paper <--------------------------
        jac[1:N + 1, N + 1 + m] = k * um * J * CC + k * vm * J * SS
    
    # Derivative of mean(eta) = 1
    jac[N + 1:2 * N + 2, -2] = M * 0 + 1 / N
    jac[N + 1, -2] = 1 / (2 * N)
    jac[2 * N + 1, -2] = 1 / (2 * N)
    
    # Derivative of the wave height criterion
    jac[N + 1, -1] = 1
    jac[2 * N + 1, -1] = -1
    
    # For debugging and comparing to fprime_num
    # print('coeffs = [%s]' % ', '.join(repr(ci) for ci in coeffs))
    # print('H, k, D, J, M =', ', '.join(repr(ci) for ci in (H, k, D, J, M)))
    
    return jac
