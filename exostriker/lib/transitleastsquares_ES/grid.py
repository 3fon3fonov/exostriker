from __future__ import division, print_function
import numba
import transitleastsquares_ES.tls_constants as tls_constants
import numpy
from numpy import pi, sqrt, arccos, degrees, floor, ceil
import warnings


@numba.jit(fastmath=True, parallel=False, nopython=True)
def T14(
    R_s, M_s, P, upper_limit=tls_constants.FRACTIONAL_TRANSIT_DURATION_MAX, small=False
):
    """Input:  Stellar radius and mass; planetary period
               Units: Solar radius and mass; days
       Output: Maximum planetary transit duration T_14max
               Unit: Fraction of period P"""

    P = P * tls_constants.SECONDS_PER_DAY
    R_s = tls_constants.R_sun * R_s
    M_s = tls_constants.M_sun * M_s

    if small:  # small planet assumption
        T14max = R_s * ((4 * P) / (pi * tls_constants.G * M_s)) ** (1 / 3)
    else:  # planet size 2 R_jup
        T14max = (R_s + 2 * tls_constants.R_jup) * (
            (4 * P) / (pi * tls_constants.G * M_s)
        ) ** (1 / 3)

    result = T14max / P
    if result > upper_limit:
        result = upper_limit
    return result


def duration_grid(periods, shortest, log_step=tls_constants.DURATION_GRID_STEP):
    
    duration_max = T14(
        R_s=tls_constants.R_STAR_MAX,
        M_s=tls_constants.M_STAR_MAX,
        P=min(periods),
        small=False  # large planet for long transit duration
    )
    duration_min = T14(
        R_s=tls_constants.R_STAR_MIN,
        M_s=tls_constants.M_STAR_MIN,
        P=max(periods),
        small=True  # small planet for short transit duration
    )

    durations = [duration_min]
    current_depth = duration_min
    while current_depth * log_step < duration_max:
        current_depth = current_depth * log_step
        durations.append(current_depth)
    durations.append(duration_max)  # Append endpoint. Not perfectly spaced.
    return durations


def period_grid(
    R_star,
    M_star,
    time_span,
    period_min=0,
    period_max=float("inf"),
    oversampling_factor=tls_constants.OVERSAMPLING_FACTOR,
    n_transits_min=tls_constants.N_TRANSITS_MIN,
):
    """Returns array of optimal sampling periods for transit search in light curves
       Following Ofir (2014, A&A, 561, A138)"""

    if R_star < 0.01:
        text = (
            "Warning: R_star was set to 0.01 for period_grid (was unphysical: "
            + str(R_star)
            + ")"
        )
        warnings.warn(text)
        R_star = 0.1

    if R_star > 10000:
        text = (
            "Warning: R_star was set to 10000 for period_grid (was unphysical: "
            + str(R_star)
            + ")"
        )
        warnings.warn(text)
        R_star = 10000

    if M_star < 0.01:
        text = (
            "Warning: M_star was set to 0.01 for period_grid (was unphysical: "
            + str(M_star)
            + ")"
        )
        warnings.warn(text)
        M_star = 0.01

    if M_star > 1000:
        text = (
            "Warning: M_star was set to 1000 for period_grid (was unphysical: "
            + str(M_star)
            + ")"
        )
        warnings.warn(text)
        M_star = 1000

    R_star = R_star * tls_constants.R_sun
    M_star = M_star * tls_constants.M_sun
    time_span = time_span * tls_constants.SECONDS_PER_DAY  # seconds

    # boundary conditions
    f_min = n_transits_min / time_span
    f_max = 1.0 / (2 * pi) * sqrt(tls_constants.G * M_star / (3 * R_star) ** 3)

    # optimal frequency sampling, Equations (5), (6), (7)
    A = (
        (2 * pi) ** (2.0 / 3)
        / pi
        * R_star
        / (tls_constants.G * M_star) ** (1.0 / 3)
        / (time_span * oversampling_factor)
    )
    C = f_min ** (1.0 / 3) - A / 3.0
    N_opt = (f_max ** (1.0 / 3) - f_min ** (1.0 / 3) + A / 3) * 3 / A

    X = numpy.arange(N_opt) + 1
    f_x = (A / 3 * X + C) ** 3
    P_x = 1 / f_x

    # Cut to given (optional) selection of periods
    periods = P_x / tls_constants.SECONDS_PER_DAY
    selected_index = numpy.where(
        numpy.logical_and(periods > period_min, periods <= period_max)
    )

    number_of_periods = numpy.size(periods[selected_index])

    if number_of_periods > 10 ** 6:
        text = (
            "period_grid generates a very large grid ("
            + str(number_of_periods)
            + "). Recommend to check physical plausibility for stellar mass, radius, and time series duration."
        )
        warnings.warn(text)

    if number_of_periods < tls_constants.MINIMUM_PERIOD_GRID_SIZE:
        if time_span < 5 * tls_constants.SECONDS_PER_DAY:
            time_span = 5 * tls_constants.SECONDS_PER_DAY
        warnings.warn(
            "period_grid defaults to R_star=1 and M_star=1 as given density yielded grid with too few values"
        )
        return period_grid(
            R_star=1, M_star=1, time_span=time_span / tls_constants.SECONDS_PER_DAY
        )
    else:
        return periods[selected_index]  # periods in [days]
