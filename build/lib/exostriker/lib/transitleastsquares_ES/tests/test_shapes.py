from __future__ import division, print_function
import numpy
import scipy
import scipy.signal
from transitleastsquares_ES import transitleastsquares_ES


def loadfile(filename):
    data = numpy.genfromtxt(filename, delimiter=",", dtype="f8, f8", names=["t", "y"])
    return data["t"], data["y"]


if __name__ == "__main__":
    print("Starting test: transit shapes...", end="")

    # Testing transit shapes
    t, y = loadfile("EPIC206154641.csv")
    trend = scipy.signal.medfilt(y, 25)
    y_filt = y / trend

    # box
    model_box = transitleastsquares_ES(t, y_filt)
    results_box = model_box.power(transit_template="box")
    numpy.testing.assert_almost_equal(
        results_box.duration, 0.06111785726416931, decimal=5)
    numpy.testing.assert_almost_equal(results_box.rp_rs, 0.08836981203437415, decimal=5)

    # grazing
    model_grazing = transitleastsquares_ES(t, y_filt)
    results_grazing = model_grazing.power(transit_template="grazing")

    numpy.testing.assert_almost_equal(
        results_grazing.duration, 0.08948265482047034, decimal=5
    )
    numpy.testing.assert_almost_equal(
        min(results_grazing.chi2red), 0.06759475703796078, decimal=5)
    print("Test passed: Grazing-shaped")



    print("Test passed: Box-shaped")
    print("All tests passed")
