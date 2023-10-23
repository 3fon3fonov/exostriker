from __future__ import division, print_function
import numpy
from transitleastsquares_ES import resample

if __name__ == "__main__":
    print("Starting test resample...", end="")

    # Resampling
    testtime = numpy.linspace(0, 1, 1000)
    testflux = numpy.linspace(0.99, 1.01, 1000)
    a, b = resample(time=testtime, flux=testflux, factor=100)
    expected_result_a = (
        0.0,
        0.11111111,
        0.22222222,
        0.33333333,
        0.44444444,
        0.55555556,
        0.66666667,
        0.77777778,
        0.88888889,
        1.0,
    )
    numpy.testing.assert_almost_equal(a, expected_result_a)
    expected_result_b = (
        0.99,
        0.99222222,
        0.99444444,
        0.99666667,
        0.99888889,
        1.00111111,
        1.00333333,
        1.00555556,
        1.00777778,
        1.01,
    )
    numpy.testing.assert_almost_equal(b, expected_result_b)
    numpy.testing.assert_equal(len(a), 10)
    numpy.testing.assert_almost_equal(min(a), 0)
    numpy.testing.assert_almost_equal(max(a), 1)
    numpy.testing.assert_almost_equal(numpy.mean(a), 0.5)
    numpy.testing.assert_equal(len(b), 10)
    numpy.testing.assert_almost_equal(min(b), 0.99)
    numpy.testing.assert_almost_equal(max(b), 1.01)
    numpy.testing.assert_almost_equal(numpy.mean(b), 1)
    print("passed")
