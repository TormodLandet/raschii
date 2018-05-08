def _test_fenton_non_steep():
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
    fwave = FentonWave(height=0.2, depth=0.5, length=1.5, N=10)
    
    ml_B = [323.589155595150e-003, 22.0140805523597e-003, 3.52575820413638e-003,
            1.07573388366430e-003, 364.710791280593e-006, 136.215131403687e-006,
            54.6358173146445e-006, 22.9669854313857e-006, 9.96497612562210e-006,
            3.19234961117461e-006]
    ml_Q = 799.771104212840e-003
    ml_R = 6.26369379030432e+000
    ml_k = 4.18879020478639e+000
    
    assert abs(fwave.B - ml_B).max() < 1e-15
    assert abs(fwave.Q - ml_Q) < 1e-15
    assert abs(fwave.R - ml_R) < 1e-15
    assert abs(fwave.coeffs['k'] - ml_k) < 1e-15
