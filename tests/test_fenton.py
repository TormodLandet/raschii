from numpy import array


def test_fenton_jacobian():
    # From case with height=0.1, depth=0.5, length=2, N=10
    coeffs = [-0.7641186499221805, -0.04165712827663715, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 1.1, 1.0951056516295155,
              1.0809016994374947, 1.0587785252292474, 1.0309016994374947, 1.0,
              0.9690983005625052, 0.9412214747707527, 0.9190983005625053,
              0.9048943483704847, 0.9, 0.7641186499221805, 1.2919386555794479]
    coeffs = array(coeffs)
    H, k, D = 0.2, 1.5707963267948966, 1
    J = array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    M = array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    
    # Compute the jacobian using the two methods
    from raschii.fenton import fprime, fprime_num
    jacA = fprime(coeffs, H, k, D, J, M)
    jacN = fprime_num(coeffs, H, k, D, J, M)
    
    # Verify that the jacobians are close
    Nc = coeffs.size
    err = False
    for i in range(Nc):
        for j in range(Nc):
            if abs(jacA[i, j] - jacN[i, j]) > 1e-5:
                print('Expected equal elements at [%d, %d], found %r and %r '
                      'with diff %r' % (i, j, jacA[i, j], jacN[i, j],
                                        jacA[i, j] - jacN[i, j]))
                err = True
    assert not err


def _test_long_flat_wave():
    from raschii import FentonWave
    fwave = FentonWave(height=1, depth=5000, length=100, N=10)
    
    h = fwave.surface_elevation(0)
    assert abs(h - 1) < 1e-6


def test_compare_fenton_m_01():
    """
    Compare with results obtained by
    https://github.com/roenby/fentonWave/blob/master/tests/fenton.m
    """
    from raschii import FentonWave
    
    fwave = FentonWave(height=0.1, depth=0.5, length=2, N=10)
    
    ml_eta = [625.633134481387e-3, 608.838129402900e-3, 574.190216147621e-3,
              538.287549720044e-3, 506.621746098552e-3, 480.480195361792e-3,
              459.890477455227e-3, 444.500105760016e-3, 433.883580207027e-3,
              427.674865588164e-3, 425.633134481387e-003]
    ml_B = [288.233415443277e-3, 19.0836854437944e-3, 821.414465583978e-6,
            131.328047779987e-6, 34.6376174890207e-6, 6.11862403532633e-6,
            1.09123607695571e-6, 238.738040611394e-6, 55.1743593679826e-6,
            11.0463615856899e-6]
    ml_Q = 799.771104212840e-003
    ml_R = 6.26369379030432e+000
    ml_k = 4.18879020478639e+000
    
    assert abs(fwave.eta - ml_eta).max() < 1e-15
    assert abs(fwave.B - ml_B).max() < 1e-15
    assert abs(fwave.Q - ml_Q) < 1e-15
    assert abs(fwave.R - ml_R) < 1e-15
    assert abs(fwave.coeffs['k'] - ml_k) < 1e-15


if __name__ == '__main__':
    # For testing the tests
    import numpy
    numpy.set_printoptions(linewidth=200, precision=3, suppress=True)
    test_fenton_jacobian()
    test_compare_fenton_m_01()
