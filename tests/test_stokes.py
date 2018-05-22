from raschii import stokes


def test_coefficients():
    kd = 0.753982
    
    # First check that the higher order coefficients are zero for lower order
    for N in range(1, 6):
        data = stokes.stokes_coefficients(kd, N)
        for k, v in data.items():
            if int(k[1]) > N:
                assert v == 0.0
            else:
                assert v != 0.0
    
    # Test against the values in "A Fifth‐Order Stokes Theory for Steady Waves"
    # by John D. Fenton (1985) given for kd = 0.753982
    print('Tests for kd = %r' % kd)
    refvals = {'A11': 1.208490,
               'A22': 0.799840,
               'A31': -9.105340,
               'A33': 0.368275,
               'A42': -12.196150,
               'A44': 0.058723,
               'A51': 108.467921,
               'A53': -6.941756,
               'A55': -0.074979,
               'B22': 2.502414,
               'B31': -5.731666,
               'B42': -32.407508,
               'B44': 14.033758,
               'B53': -103.445042,
               'B55': 37.200027,
               'C0': 0.798448,
               'C2': 1.940215,
               'C4': -12.970403,
               'D2': -0.626215,
               'D4': 3.257104,
               'E2': 1.781926,
               'E4': -11.573657}
    for name, ref in refvals.items():
        calc = data[name]
        rerr = abs(calc - ref) / abs(ref)
        print('%3s: %15.5e' % (name, rerr))
        assert rerr < 1e-5, 'Discrepancy found for %s: %e' % (name, rerr)
    
    # Check for overflow prevention
    for kd in range(1, 1000):
        stokes.stokes_coefficients(kd, N)
