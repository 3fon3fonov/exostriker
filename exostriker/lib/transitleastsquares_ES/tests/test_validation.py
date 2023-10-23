from __future__ import division, print_function
import numpy
from transitleastsquares_ES import transitleastsquares_ES

if __name__ == "__main__":

    t = y = numpy.linspace(0, 1, 2000)

    try:
        model = transitleastsquares_ES(t, y)
        results = model.power(use_threads=0)
    except:
        print("Test passed: use_threads=0")
    try:
        model = transitleastsquares_ES(t, y)
        results = model.power(use_threads="1")
    except:
        print("Test passed: use_threads='1'")
