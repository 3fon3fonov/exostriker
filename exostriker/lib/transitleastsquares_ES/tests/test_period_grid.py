from __future__ import division, print_function
import numpy
from transitleastsquares_ES import period_grid

if __name__ == "__main__":
    print("Starting test: period_grid...", end="")

    periods = period_grid(R_star=1, M_star=1, time_span=0.1)
    numpy.testing.assert_almost_equal(max(periods), 2.4999999999999987)
    numpy.testing.assert_almost_equal(min(periods), 0.6002621413799498)
    numpy.testing.assert_equal(len(periods), 268)

    periods = period_grid(R_star=1, M_star=1, time_span=20)
    numpy.testing.assert_almost_equal(max(periods), 10)
    numpy.testing.assert_almost_equal(min(periods), 0.6015575922909607)
    numpy.testing.assert_equal(len(periods), 1716)

    periods = period_grid(
        R_star=5,  # R_sun
        M_star=1,  # M_sun
        time_span=20,  # days
        period_min=0,
        period_max=999,
        oversampling_factor=3,
    )
    numpy.testing.assert_almost_equal(max(periods), 10)
    numpy.testing.assert_almost_equal(min(periods), 0.6015575922909607)
    numpy.testing.assert_equal(len(periods), 1716)

    periods = period_grid(
        R_star=1,  # R_sun
        M_star=1,  # M_sun
        time_span=20,  # days
        period_min=0,
        period_max=999,
        oversampling_factor=3,
    )
    numpy.testing.assert_almost_equal(max(periods), 10)
    numpy.testing.assert_almost_equal(min(periods), 0.60155759)
    numpy.testing.assert_equal(len(periods), 1716)

    periods = period_grid(
        R_star=0.1,  # R_sun
        M_star=1,  # M_sun
        time_span=1000,  # days
        period_min=0,
        period_max=999,
        oversampling_factor=3,
    )
    numpy.testing.assert_equal(len(periods), 4308558)
    print("passed")
