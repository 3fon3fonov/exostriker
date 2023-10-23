from __future__ import division, print_function
import numpy
import numba
import transitleastsquares_ES.tls_constants as tls_constants
from transitleastsquares_ES.grid import T14
from transitleastsquares_ES.helpers import running_mean


@numba.jit(fastmath=True, parallel=False, nopython=True)
def fold(time, period, T0):
    """Normal phase folding"""
    return (time - T0) / period - numpy.floor((time - T0) / period)


@numba.jit(fastmath=True, parallel=False, nopython=True)
def foldfast(time, period):
    """Fast phase folding with T0=0 hardcoded"""
    return time / period - numpy.floor(time / period)


@numba.jit(fastmath=True, parallel=False, nopython=True)
def edge_effect_correction(flux, patched_data, dy, inverse_squared_patched_dy):
    regular = numpy.sum(((1 - flux) ** 2) * 1 / dy ** 2)
    patched = numpy.sum(((1 - patched_data) ** 2) * inverse_squared_patched_dy)
    return patched - regular


@numba.jit(fastmath=True, parallel=False, nopython=True)
def lowest_residuals_in_this_duration(
    mean,
    transit_depth_min,
    patched_data_arr,
    duration,
    signal,
    inverse_squared_patched_dy_arr,
    overshoot,
    ootr,
    summed_edge_effect_correction,
    chosen_transit_row,
    datapoints,
    T0_fit_margin,
):

    # if nothing is fit, we fit a straight line: signal=1. Then, at dy=1,
    # the squared sum of residuals equals the number of datapoints
    summed_residual_in_rows = datapoints
    best_row = 0
    best_depth = 0

    xth_point = 1  # How many cadences the template shifts forward in each step
    if T0_fit_margin > 0 and duration > T0_fit_margin:
        T0_fit_margin = 1 / T0_fit_margin
        xth_point = int(duration / T0_fit_margin)
        if xth_point < 1:
            xth_point = 1

    for i in range(len(mean)):
        if (mean[i] > transit_depth_min) and (i % xth_point == 0):
            data = patched_data_arr[i : i + duration]
            dy = inverse_squared_patched_dy_arr[i : i + duration]
            target_depth = mean[i] * overshoot
            scale = tls_constants.SIGNAL_DEPTH / target_depth
            reverse_scale = 1 / scale  # speed: one division now, many mults later

            # Scale model and calculate residuals
            intransit_residual = 0
            for j in range(len(signal)):
                sigi = (1 - signal[j]) * reverse_scale
                intransit_residual += ((data[j] - (1 - sigi)) ** 2) * dy[j]
            current_stat = intransit_residual + ootr[i] - summed_edge_effect_correction
            if current_stat < summed_residual_in_rows:
                summed_residual_in_rows = current_stat
                best_row = chosen_transit_row
                best_depth = 1 - target_depth

    return summed_residual_in_rows, best_row, best_depth


@numba.jit(fastmath=True, parallel=False, nopython=True)
def out_of_transit_residuals(data, width_signal, dy):
    chi2 = numpy.zeros(len(data) - width_signal + 1)
    fullsum = numpy.sum(((1 - data) ** 2) * dy)
    window = numpy.sum(((1 - data[:width_signal]) ** 2) * dy[:width_signal])
    chi2[0] = fullsum - window
    for i in range(1, len(data) - width_signal + 1):
        becomes_visible = i - 1
        becomes_invisible = i - 1 + width_signal
        add_visible_left = (1 - data[becomes_visible]) ** 2 * dy[becomes_visible]
        remove_invisible_right = (1 - data[becomes_invisible]) ** 2 * dy[
            becomes_invisible
        ]
        chi2[i] = chi2[i - 1] + add_visible_left - remove_invisible_right
    return chi2


def search_period(
    period,
    t,
    y,
    dy,
    transit_depth_min,
    R_star_min,
    R_star_max,
    M_star_min,
    M_star_max,
    lc_arr,
    lc_cache_overview,
    T0_fit_margin,
):
    """Core routine to search the flux data set 'injected' over all 'periods'"""

    # duration (in samples) of widest transit in lc_cache (axis 0: rows; axis 1: columns)
    durations = numpy.unique(lc_cache_overview["width_in_samples"])
    maxwidth_in_samples = int(max(durations))
    if maxwidth_in_samples % 2 != 0:
        maxwidth_in_samples = maxwidth_in_samples + 1

    # Phase fold
    phases = foldfast(t, period)
    sort_index = numpy.argsort(phases, kind="mergesort")  # 8% faster than Quicksort
    phases = phases[sort_index]
    flux = y[sort_index]
    dy = dy[sort_index]

    # faster to multiply than divide
    patched_dy = numpy.append(dy, dy[:maxwidth_in_samples])
    inverse_squared_patched_dy = 1 / patched_dy ** 2

    # Due to phase folding, the signal could start near the end of the data
    # and continue at the beginning. To avoid (slow) rolling,
    # we patch the beginning again to the end of the data
    patched_data = numpy.append(flux, flux[:maxwidth_in_samples])

    this_edge_effect_correction = edge_effect_correction(
        flux, patched_data, dy, inverse_squared_patched_dy
    )

    # Set "best of" counters to max, in order to find smaller residuals
    smallest_residuals_in_period = float("inf")
    summed_residual_in_rows = float("inf")

    # Make unique to avoid duplicates in dense grids
    duration_max = T14(R_s=R_star_max, M_s=M_star_max, P=period, small=False)
    duration_min = T14(R_s=R_star_min, M_s=M_star_min, P=period, small=True)

    # Fractional transit duration can be longer than this.
    # Example: Data length 11 days, 2 transits at 0.5 days and 10.5 days
    length = max(t) - min(t)
    no_of_transits_naive = length / period
    no_of_transits_worst = no_of_transits_naive + 1
    correction_factor = no_of_transits_worst / no_of_transits_naive

    duration_min_in_samples = int(numpy.floor(duration_min * len(y)))
    duration_max_in_samples = int(numpy.ceil(duration_max * len(y) * correction_factor))
    durations = durations[durations >= duration_min_in_samples]
    durations = durations[durations <= duration_max_in_samples]

    skipped_all = True
    best_row = 0  # shortest and shallowest transit
    best_depth = 0

    for duration in durations:
        chosen_transit_row = 0
        while lc_cache_overview["width_in_samples"][chosen_transit_row] != duration:
            chosen_transit_row += 1
        this_residual, this_row, this_depth = lowest_residuals_in_this_duration(
            mean=1 - running_mean(patched_data, duration),
            transit_depth_min=transit_depth_min,
            patched_data_arr=patched_data,
            duration=duration,
            signal=lc_arr[chosen_transit_row],
            inverse_squared_patched_dy_arr=inverse_squared_patched_dy,
            overshoot=lc_cache_overview["overshoot"][chosen_transit_row],
            ootr=out_of_transit_residuals(
                patched_data, duration, inverse_squared_patched_dy
            ),
            summed_edge_effect_correction=this_edge_effect_correction,
            chosen_transit_row=chosen_transit_row,
            datapoints=len(flux),
            T0_fit_margin=T0_fit_margin,
        )

        if this_residual < summed_residual_in_rows:
            summed_residual_in_rows = this_residual
            best_row = chosen_transit_row
            best_depth = this_depth

    return [period, summed_residual_in_rows, best_row, best_depth]
