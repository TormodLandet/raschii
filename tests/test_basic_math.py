import numpy

from raschii.common import sinh_by_cosh, cosh_by_cosh, cosh_ratio


def test_sinh_by_cosh():
    end = 45
    for f in numpy.linspace(0.001, 2, 100):
        # Compute the two approximations
        a = numpy.linspace(0, end, 1001)
        b = numpy.linspace(0, end, 1001) * f
        f1 = numpy.sinh(a) / numpy.cosh(b)
        f2 = sinh_by_cosh(a, b)
        check_arrays(f1, f2, 1e-5, 1e-12)

    # Some handpicked tests
    a = numpy.array([0.0, 0.0, 1.0, 1.0, 199.0], float)
    b = numpy.array([0.0, 1.0, 0.0, 1.0, 199.0], float)
    f1 = numpy.sinh(a) / numpy.cosh(b)
    f2 = sinh_by_cosh(a, b)
    check_arrays(f1, f2, 1e-5, 1e-12)


def test_cosh_by_cosh():
    end = 45
    for f in numpy.linspace(0.001, 2, 100):
        # Compute the two approximations
        a = numpy.linspace(0, end, 1001)
        b = numpy.linspace(0, end, 1001) * f
        f1 = numpy.cosh(a) / numpy.cosh(b)
        f2 = cosh_by_cosh(a, b)
        check_arrays(f1, f2, 1e-5, 1e-12)

    # Some handpicked tests
    a = numpy.array([0.0, 0.0, 1.0, 1.0, 199.0], float)
    b = numpy.array([0.0, 1.0, 0.0, 1.0, 199.0], float)
    f1 = numpy.cosh(a) / numpy.cosh(b)
    f2 = cosh_by_cosh(a, b)
    check_arrays(f1, f2, 1e-5, 1e-12)


def test_cosh_ratio():
    """cosh(a)/cosh(b) stably for 0 <= a <= b, including overflow of cosh(b)."""
    # --- Normal regime: compare against direct float64 computation ---
    for b_end in [1.0, 10.0, 45.0]:
        b = numpy.linspace(0.0, b_end, 501)
        for frac in numpy.linspace(0.0, 1.0, 11):
            a = b * frac
            expected = numpy.cosh(a) / numpy.cosh(b)
            got = cosh_ratio(a, b)
            check_arrays(expected, got, 1e-10, 1e-12)

    # a == b should give exactly 1.0
    a = numpy.array([0.0, 0.5, 5.0, 30.0])
    b = numpy.array([0.0, 0.5, 5.0, 30.0])
    numpy.testing.assert_allclose(cosh_ratio(a, b), 1.0, rtol=1e-12)

    # a == 0 (bottom boundary) in non-overflow range
    b = numpy.linspace(0.0, 45.0, 501)
    a = numpy.zeros_like(b)
    expected = 1.0 / numpy.cosh(b)
    check_arrays(expected, cosh_ratio(a, b), 1e-10, 1e-12)

    # --- Overflow regime: cosh(b) overflows float64 for b > ~709 ---
    # Exact: cosh(a)/cosh(b) = exp(a-b) * (1+exp(-2a))/(1+exp(-2b))
    # For a > 30: exp(-2a) < 1e-26, so the approximation exp(a-b) is exact
    # to within double precision.
    # Start a at b-600 so that exp(a-b) > exp(-600) ~ 1e-261, well above
    # the subnormal range (where relative-error comparison is meaningless).
    for b_val in [800.0, 1500.0]:
        a_start = b_val - 600.0  # exp(a_start - b_val) = exp(-600) ~ 1e-261
        a_vals = numpy.linspace(a_start, b_val, 50)  # a > 30 throughout
        b_vals = numpy.full(50, b_val)
        expected = numpy.exp(a_vals - b_vals)
        check_arrays(expected, cosh_ratio(a_vals, b_vals), 1e-10, 1e-12)

    # Scalar overflow: single large-b call must not raise
    got_scalar = cosh_ratio(numpy.array([700.0]), numpy.array([800.0]))
    numpy.testing.assert_allclose(got_scalar, numpy.exp(700.0 - 800.0), rtol=1e-12)

    # --- Arbitrary-shape (2-D) arrays ---
    a2 = numpy.array([[1.0, 5.0], [10.0, 20.0]])
    b2 = numpy.array([[2.0, 7.0], [15.0, 25.0]])
    expected2 = numpy.cosh(a2) / numpy.cosh(b2)
    numpy.testing.assert_allclose(cosh_ratio(a2, b2), expected2, rtol=1e-12, atol=1e-12)


def check_arrays(f1, f2, atol, rtol, atol2=1e5, atol2_lim=1e10):
    """
    Compute the absolute and the relative error
    """
    assert len(f1) == len(f2)
    err = abs(f1 - f2)
    for i, e1 in enumerate(err):
        if e1 == 0:
            continue

        # Relative error
        e2 = e1 / f1[i] if f1[i] != 0 else err

        # Change atol if f1 is VERY large
        at = atol if abs(f1[i]) < atol2_lim else atol2

        # Check for errors
        if e1 > at or e2 > rtol:
            print("Found abserr %r (tol: %r) and relerr %r (tol: %r)" % (e1, at, e2, rtol))
            print("i = %r, f1[imax] = %r, f2[imax] = %r" % (i, f1[i], f2[i]))
        assert e1 < at
        assert e2 < rtol
