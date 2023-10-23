from __future__ import division, print_function
import numpy
import scipy
import scipy.signal
from transitleastsquares_ES import transitleastsquares_ES, transit_mask, cleaned_array


def loadfile(filename):
    data = numpy.genfromtxt(filename, delimiter=",", dtype="f8, f8", names=["t", "y"])
    return data["t"], data["y"]


if __name__ == "__main__":
    print("Starting test: Multi-planet...", end="")
    t, y = loadfile("EPIC201367065.csv")
    trend = scipy.signal.medfilt(y, 25)
    y_filt = y / trend

    model = transitleastsquares_ES(t, y_filt)
    results = model.power()

    numpy.testing.assert_almost_equal(max(results.power), 45.49085809486116, decimal=3)
    numpy.testing.assert_almost_equal(
        max(results.power_raw), 42.93056655774114, decimal=3
    )
    numpy.testing.assert_almost_equal(min(results.power), -0.6175100139942546, decimal=3)
    numpy.testing.assert_almost_equal(
        min(results.power_raw), -0.3043720539933344, decimal=3
    )
    print("Detrending of power spectrum from power_raw passed")

    # Mask of the first planet
    intransit = transit_mask(t, results.period, 2 * results.duration, results.T0)
    y_second_run = y_filt[~intransit]
    t_second_run = t[~intransit]
    t_second_run, y_second_run = cleaned_array(t_second_run, y_second_run)

    # Search for second planet
    model_second_run = transitleastsquares_ES(t_second_run, y_second_run)
    results_second_run = model_second_run.power()
    numpy.testing.assert_almost_equal(
        results_second_run.duration, 0.15061016994013998, decimal=3
    )
    numpy.testing.assert_almost_equal(
        results_second_run.SDE, 34.9911304598618, decimal=3
    )
    numpy.testing.assert_almost_equal(
        results_second_run.rp_rs, 0.025852178872027086, decimal=3
    )

    print("Passed")
