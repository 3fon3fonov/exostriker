import numpy as np

def running_median(data, kernel):
    """Returns sliding median of width 'kernel' and same length as data """
    idx = numpy.arange(kernel) + numpy.arange(len(data) - kernel + 1)[:, None]
#    idx = idx.astype(numpy.int)  # needed if oversampling_factor is not int
    med = numpy.median(data[idx], axis=1)

    # Append the first/last value at the beginning/end to match the length of
    # data and returned median
    first_values = med[0]
    last_values = med[-1]
    missing_values = len(data) - len(med)
    values_front = int(missing_values * 0.5)
    values_end = missing_values - values_front
    med = numpy.append(numpy.full(values_front, first_values), med)
    med = numpy.append(med, numpy.full(values_end, last_values))
    return med


