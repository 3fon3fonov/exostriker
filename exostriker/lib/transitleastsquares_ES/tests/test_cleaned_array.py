from __future__ import division, print_function
import numpy
from transitleastsquares_ES import cleaned_array

if __name__ == "__main__":
    print("Starting test: cleaned_array...", end="")

    dirty_array = numpy.ones(10, dtype=object)
    time_array = numpy.linspace(1, 10, 10)
    dy_array = numpy.ones(10, dtype=object)
    dirty_array[1] = None
    dirty_array[2] = numpy.inf
    dirty_array[3] = -numpy.inf
    dirty_array[4] = numpy.nan
    dirty_array[5] = -99
    time_array[8] = numpy.nan
    dy_array[9] = numpy.inf

    t, y, dy = cleaned_array(time_array, dirty_array, dy_array)
    numpy.testing.assert_equal(t, [1, 7, 8])
    numpy.testing.assert_equal(y, [1, 1, 1])
    numpy.testing.assert_equal(dy, [1, 1, 1])
    print("passed")
