from __future__ import division, print_function
import numpy
import batman
from transitleastsquares_ES import transitleastsquares_ES


if __name__ == "__main__":
    print("Starting test: uncertainties...", end="")
    numpy.random.seed(seed=0)  # reproducibility

    # Create test data
    start = 48
    days = 365.25 * 3
    samples_per_day = 12  # 48
    samples = int(days * samples_per_day)  # 48
    t = numpy.linspace(start, start + days, samples)

    # Use batman to create transits
    ma = batman.TransitParams()
    ma.t0 = (
        start + 20
    )  # time of inferior conjunction; first transit is X days after start
    ma.per = 365.25  # orbital period
    ma.rp = 6371 / 696342  # 6371 planet radius (in units of stellar radii)
    ma.a = 217  # semi-major axis (in units of stellar radii)
    ma.inc = 90  # orbital inclination (in degrees)
    ma.ecc = 0  # eccentricity
    ma.w = 90  # longitude of periastron (in degrees)
    ma.u = [0.5]  # limb darkening coefficients
    ma.limb_dark = "linear"  # limb darkening model
    m = batman.TransitModel(ma, t)  # initializes model
    original_flux = m.light_curve(ma)  # calculates light curve

    # Create noise and merge with flux
    ppm = 5
    stdev = 10 ** -6 * ppm
    noise = numpy.random.normal(0, stdev, int(samples))
    y = original_flux + noise

    # Inject excess noise near the end of the time series.
    # When searching without uncertainties, the SDE is 4.827
    # When searching with uncertainties, the SDE is larger, 5.254. Test passed!
    noise = numpy.random.normal(0, 10 * stdev, 3149)
    y[10000:] = y[10000:] + noise

    dy = numpy.full(len(y), stdev)
    dy[10000:] = 10 * stdev

    model = transitleastsquares_ES(t, y, dy)
    results = model.power(
        period_min=360,
        period_max=370,
        oversampling_factor=3,
        duration_grid_step=1.05,
        T0_fit_margin=0.2,
    )
    numpy.testing.assert_almost_equal(results.SDE, 5.292594615900944, decimal=5)
    print("passed")
