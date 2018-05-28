import numpy
from numpy import array


def test_sinh_by_cosh():
    from raschii.fenton import sinh_by_cosh
    
    end = 45
    worst_err1 = worst_err2 = 0
    for f in numpy.linspace(0.001, 2, 100):
        # Compute the two approximations
        a = numpy.linspace(0, end, 1001)
        b = numpy.linspace(0, end, 1001) * f
        f1 = numpy.sinh(a) / numpy.cosh(b)
        f2 = sinh_by_cosh(a, b)
        
        # Compute the absolute and the relative error
        err = abs(f1 - f2)
        imax = numpy.argmax(err)
        if err.max() == 0:
            e1 = e2 = 0
        else:
            e1 = err.max()
            e2 = err[imax] / f1[imax]
        
        worst_err1 = max(worst_err1, e1)
        worst_err2 = max(worst_err2, e2)
    print('Worst error is', worst_err1, worst_err2)
    assert worst_err1 < 1e-6
    assert worst_err2 < 1e-12


def test_cosh_by_cosh():
    from raschii.fenton import cosh_by_cosh
    
    end = 45
    worst_err1 = worst_err2 = 0
    for f in numpy.linspace(0.001, 2, 100):
        # Compute the two approximations
        a = numpy.linspace(0, end, 1001)
        b = numpy.linspace(0, end, 1001) * f
        f1 = numpy.cosh(a) / numpy.cosh(b)
        f2 = cosh_by_cosh(a, b)
        
        # Compute the absolute and the relative error
        err = abs(f1 - f2)
        imax = numpy.argmax(err)
        if err.max() == 0:
            e1 = e2 = 0
        else:
            e1 = err.max()
            e2 = err[imax] / f1[imax]
        
        worst_err1 = max(worst_err1, e1)
        worst_err2 = max(worst_err2, e2)
    print('Worst error is', worst_err1, worst_err2)
    assert worst_err1 < 1e-6
    assert worst_err2 < 1e-12


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
    assert abs(h - 0.5) < 1e-6


def test_compare_fenton_m_01():
    """
    Compare with results obtained by
    https://github.com/roenby/fentonWave/blob/master/tests/fenton.m
    """
    from raschii import FentonWave, check_breaking_criteria
    height = 0.2
    depth = 0.5
    length = 2
    N = 30
    print(*check_breaking_criteria(height, depth, length))
    fwave = FentonWave(height, depth, length, N)
    py_res = fwave.data
    
    ml_eta = array([6.256332118992537e-01, 6.235402803532166e-01,
                    6.176398294183506e-01, 6.088383657880129e-01,
                    5.981433292625251e-01, 5.863925484803192e-01,
                    5.741905043147202e-01, 5.619444604418909e-01,
                    5.499201111324211e-01, 5.382877714262410e-01,
                    5.271547620101914e-01, 5.165868059915903e-01,
                    5.066219276102671e-01, 4.972795106660977e-01,
                    4.885662821396637e-01, 4.804803408120494e-01,
                    4.730139306009313e-01, 4.661553962008353e-01,
                    4.598905967179053e-01, 4.542039528250898e-01,
                    4.490792405884328e-01, 4.445002058220727e-01,
                    4.404510478467331e-01, 4.369168054556746e-01,
                    4.338836674511877e-01, 4.313392232541929e-01,
                    4.292726644934076e-01, 4.276749454052453e-01,
                    4.265389076728026e-01, 4.258593738023461e-01,
                    4.256332118992525e-01])

    ml_x = array([0.000000000000000e+00, 3.333333333333333e-02,
                  6.666666666666667e-02, 1.000000000000000e-01,
                  1.333333333333333e-01, 1.666666666666667e-01,
                  2.000000000000000e-01, 2.333333333333333e-01,
                  2.666666666666667e-01, 3.000000000000000e-01,
                  3.333333333333333e-01, 3.666666666666666e-01,
                  4.000000000000000e-01, 4.333333333333334e-01,
                  4.666666666666667e-01, 5.000000000000000e-01,
                  5.333333333333333e-01, 5.666666666666667e-01,
                  6.000000000000000e-01, 6.333333333333333e-01,
                  6.666666666666666e-01, 7.000000000000000e-01,
                  7.333333333333333e-01, 7.666666666666667e-01,
                  8.000000000000000e-01, 8.333333333333334e-01,
                  8.666666666666668e-01, 9.000000000000000e-01,
                  9.333333333333333e-01, 9.666666666666667e-01,
                  1.000000000000000e+00])

    ml_B = array([1.799024665755609e+00, 2.882335299753824e-01,
                  1.908361026480882e-02, 8.213767815747559e-04,
                  1.313097979767291e-04, 3.462697008046542e-05,
                  6.111647882537242e-06, 1.086244274356666e-06,
                  2.348957533655237e-07, 5.188765616282414e-08,
                  1.118382204680194e-08, 2.478538275444786e-09,
                  5.663408583688279e-10, 1.314393502678408e-10,
                  3.145814127549186e-11, 8.003417970111735e-12,
                  3.322723890828176e-12, 3.846369798481719e-13,
                  1.962769427566047e-12, -4.233208440136547e-13,
                  1.977390134075710e-12, -1.407633578386493e-12,
                  2.704910823483152e-12, -1.465002188938580e-12,
                  1.163053015764832e-12, 8.764635112300480e-13,
                  -1.550809821892457e-12, 3.324853049346685e-12,
                  -3.355716762751551e-12, 4.576015189520344e-12,
                  -2.400972953139157e-12])

    ml_res = {'eta': ml_eta, 'x': ml_x, 'B': ml_B,
              'Q': 0.873795415021738, 'R': 6.53355446469293,
              'k': 3.14159265358979, 'c': ml_B[0]}

    # Scale the fenton.m results
    ml_B[1:] *= (9.81 / ml_res['k']**3)**0.5
    
    if False:
        from matplotlib import pyplot
        pyplot.plot(fwave.x, fwave.eta)
        pyplot.plot(ml_x, ml_eta)
        pyplot.show()
    
    has_err = False
    for name in ml_res:
        ml = ml_res[name]
        py = py_res[name]
        print('%s ml:\n%r\n%s py:\n%r' % (name, ml, name, py))
        
        if hasattr(ml, 'size'):
            err = 1e100 if ml.size != py.size else abs(py - ml).max()
        else:
            err = abs(py - ml)
        
        if err > 1e-7:
            print('ERROR %s: %r' % (name, err))
            has_err = True
        else:
            print('diff %s: %r' % (name, err))
    
    assert not has_err


def test_fenton_stream_function_and_slope():
    from raschii import FentonWave
    height = 10.0
    depth = 200.0
    length = 100.0
    N = 5
    fwave = FentonWave(height, depth, length, N)
    
    # Compare velocities with numerical differentiation of the stream function
    eps = 1e-7
    for x in numpy.linspace(0, length, 21):
        z = depth + height / 2
        vel = fwave.velocity(x, z, all_points_wet=True)
        sf0 = fwave.stream_function(x, z, frame='c')
        sfX = fwave.stream_function(x + eps, z, frame='c')
        sfZ = fwave.stream_function(x, z + eps, frame='c')
        assert vel.shape == (1, 2) and sf0.shape == (1,) and sfX.shape == (1,)
        sfvel_x = (sfZ[0] - sf0) / eps
        sfvel_z = -(sfX[0] - sf0) / eps
        print('x: %r, z: %r, vel: %r, vel_num: %r' %
              (x, z, vel, (sfvel_x, sfvel_z)))
        assert abs(vel[0, 0] - sfvel_x) < 1e-5
        assert abs(vel[0, 1] - sfvel_z) < 1e-5
    
    # Compare slope with numerical differentiation of the elevation
    for x in numpy.linspace(0, length, 21):
        e0 = fwave.surface_elevation(x)
        slope = fwave.surface_slope(x)
        assert e0.shape == (1, ) and slope.shape == (1, )
        
        e1 = fwave.surface_elevation(x + eps)
        slope_num = (e1 - e0) / eps
        print('x: %r, eta: %r, slope: %r, slope_num: %r' %
              (x, e0, slope, slope_num))
        assert abs(slope[0] - slope_num) < 1e-5
