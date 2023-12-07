from __future__ import division, print_function
import batman  # https://www.cfa.harvard.edu/~lkreidberg/batman/
import numpy
import transitleastsquares_ES.tls_constants as tls_constants
from transitleastsquares_ES.interpolation import interp1d


def reference_transit(samples, per, rp, a, inc, ecc, w, u, limb_dark):
    """Returns an Earth-like transit of width 1 and depth 1"""

    f = numpy.ones(tls_constants.SUPERSAMPLE_SIZE)
    duration = 1  # transit duration in days. Increase for exotic cases
    t = numpy.linspace(-duration * 0.5, duration * 0.5, tls_constants.SUPERSAMPLE_SIZE)
    ma = batman.TransitParams()
    ma.t0 = 0  # time of inferior conjunction
    ma.per = per  # orbital period, use Earth as a reference
    ma.rp = rp  # planet radius (in units of stellar radii)
    ma.a = a  # semi-major axis (in units of stellar radii)
    ma.inc = inc  # orbital inclination (in degrees)
    ma.ecc = ecc  # eccentricity
    ma.w = w  # longitude of periastron (in degrees)
    ma.u = u  # limb darkening coefficients
    ma.limb_dark = limb_dark  # limb darkening model
    m = batman.TransitModel(ma, t)  # initializes model
    flux = m.light_curve(ma)  # calculates light curve

    # Determine start of transit (first value < 1)
    idx_first = numpy.argmax(flux < 1)
    intransit_flux = flux[idx_first : -idx_first + 1]
    intransit_time = t[idx_first : -idx_first + 1]

    # Downsample (bin) to target sample size
    x_new = numpy.linspace(t[idx_first], t[-idx_first - 1], samples)
    f = interp1d(x_new, intransit_time)
    downsampled_intransit_flux = f(intransit_flux)

    # Rescale to height [0..1]
    rescaled = (numpy.min(downsampled_intransit_flux) - downsampled_intransit_flux) / (
        numpy.min(downsampled_intransit_flux) - 1
    )

    return rescaled


def fractional_transit(
    duration,
    maxwidth,
    depth,
    samples,
    per,
    rp,
    a,
    inc,
    ecc,
    w,
    u,
    limb_dark,
    cached_reference_transit=None,
):
    """Returns a scaled reference transit with fractional width and depth"""

    if cached_reference_transit is None:
        reference_flux = reference_transit(
            samples=samples,
            per=per,
            rp=rp,
            a=a,
            inc=inc,
            ecc=ecc,
            w=w,
            u=u,
            limb_dark=limb_dark,
        )
    else:
        reference_flux = cached_reference_transit

    # Interpolate to shorter interval - new method without scipy
    reference_time = numpy.linspace(-0.5, 0.5, samples)
    occupied_samples = int((duration / maxwidth) * samples)
    x_new = numpy.linspace(-0.5, 0.5, occupied_samples)
    f = interp1d(x_new, reference_time)
    y_new = f(reference_flux)

    # Patch ends with ones ("1")
    missing_samples = samples - occupied_samples
    emtpy_segment = numpy.ones(int(missing_samples * 0.5))
    result = numpy.append(emtpy_segment, y_new)
    result = numpy.append(result, emtpy_segment)
    if numpy.size(result) < samples:  # If odd number of samples
        result = numpy.append(result, numpy.ones(1))

    # Depth rescaling
    result = 1 - ((1 - result) * depth)

    return result


def get_cache(durations, maxwidth_in_samples, per, rp, a, inc, ecc, w, u,
              limb_dark, verbose=True):
    """Fetches (size(durations)*size(depths)) light curves of length 
        maxwidth_in_samples and returns these LCs in a 2D array, together with 
        their metadata in a separate array."""

    if verbose:
        print("Creating model cache for", str(len(durations)), "durations")
    lc_arr = []
    rows = numpy.size(durations)
    lc_cache_overview = numpy.zeros(
        rows,
        dtype=[("duration", "f8"), ("width_in_samples", "i8"), ("overshoot", "f8")],
    )
    cached_reference_transit = reference_transit(
        samples=maxwidth_in_samples,
        per=per,
        rp=rp,
        a=a,
        inc=inc,
        ecc=ecc,
        w=w,
        u=u,
        limb_dark=limb_dark,
    )

    row = 0
    for duration in durations:
        scaled_transit = fractional_transit(
            duration=duration,
            maxwidth=numpy.max(durations),
            depth=tls_constants.SIGNAL_DEPTH,
            samples=maxwidth_in_samples,
            per=per,
            rp=rp,
            a=a,
            inc=inc,
            ecc=ecc,
            w=w,
            u=u,
            limb_dark=limb_dark,
            cached_reference_transit=cached_reference_transit,
        )
        lc_cache_overview["duration"][row] = duration
        used_samples = int((duration / numpy.max(durations)) * maxwidth_in_samples)
        lc_cache_overview["width_in_samples"][row] = used_samples
        full_values = numpy.where(
            scaled_transit < (1 - tls_constants.NUMERICAL_STABILITY_CUTOFF)
        )
        first_sample = numpy.min(full_values)
        last_sample = numpy.max(full_values) + 1
        signal = scaled_transit[first_sample:last_sample]
        lc_arr.append(signal)

        # Fraction of transit bottom and mean flux
        overshoot = numpy.mean(signal) / numpy.min(signal)

        # Later, we multiply the inverse fraction ==> convert to inverse percentage
        lc_cache_overview["overshoot"][row] = 1 / (2 - overshoot)
        row += +1

    lc_arr = numpy.array(lc_arr, dtype=object)
    return lc_cache_overview, lc_arr
