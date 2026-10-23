"""Readable posterior values with precision set by the uncertainty."""
import numpy as np


def number(value, uncertainty=None, rounded=True):
    if not rounded:
        return str(float(value))
    if not np.isfinite(value):
        return 'undefined'
    if uncertainty is None or not np.isfinite(uncertainty) or uncertainty == 0:
        return f'{value:.6g}'
    exponent = int(np.floor(np.log10(abs(uncertainty))))
    decimals = 1-exponent  # Two significant digits in the uncertainty.
    if abs(value) >= 1e7 or (value != 0 and abs(value) < 1e-4):
        return f'{value:.{max(1, int(np.floor(np.log10(abs(value))))-exponent+1)}e}'
    return f'{round(float(value), decimals):.{max(0, decimals)}f}'


def interval_text(center, lower, upper, rounded=True):
    minus, plus = center-lower, upper-center
    scale = min((abs(x) for x in (minus, plus) if np.isfinite(x) and x != 0), default=0.)
    return f'{number(center, scale, rounded)} - {number(minus, scale, rounded)} + {number(plus, scale, rounded)}'


def print_summary(samples, labels, mode='median', reference=None, logl=None,
                  rounded=True, show_mad=False, stable=False):
    samples = np.asarray(samples)
    suffix = ' (stable)' if stable else ''
    print('\n'+{'mean': 'Means', 'median': 'Medians', 'best': 'Best-fit'}[mode]+' and their 68.3% intervals'+suffix)
    best = None
    if logl is not None and np.isfinite(logl).any():
        best = int(np.nanargmax(logl))
    for j, label in enumerate(labels):
        x = samples[:, j]
        x = x[np.isfinite(x)]
        if not len(x):
            print(label+' = undefined')
            continue
        low, high = np.percentile(x, [15.85, 84.15])
        center = np.mean(x) if mode == 'mean' else np.median(x)
        if mode == 'best':
            center = samples[best, j] if best is not None else reference[j] if reference is not None else center
        print(label+' = '+interval_text(center, low, high, rounded))
    if show_mad:
        print('\nMedian Absolute Deviation values'+suffix)
        for j, label in enumerate(labels):
            x = samples[:, j]; x = x[np.isfinite(x)]
            mad = np.median(np.abs(x-np.median(x))) if len(x) else np.nan
            print(label+'  MAD = '+number(mad, mad, rounded))
