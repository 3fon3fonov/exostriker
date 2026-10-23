"""Continuous plotting branches for periodic angles; never mutate posteriors."""
import re
import numpy as np


def is_periodic_label(label):
    label = str(label)
    text = label.lower()
    if any(word in text for word in ('ampl', 'semiamp', 'excursion', 'resultant', 'concentration')):
        return False
    if 'sin(' in text or 'cos(' in text or 'sin(' in text.replace(' ', '') or 'cos(' in text.replace(' ', ''):
        return False
    # Historical make_stab.py uses capital Delta for angular amplitudes.
    if label.strip().startswith((r'$\Delta\omega', r'$\Delta\theta')):
        return False
    if re.search(r'pl\d+.*_(?:omega|Omega|varpi|lambda|M|delta_omega|delta_varpi|theta\d+)_(?:mean)$', label):
        return True
    if 'MMR' in label and re.search(r'_theta\d+_mean$', label):
        return True
    return (any(symbol in label for symbol in (r'\theta', r'\omega', r'\Omega', r'\varpi', r'\lambda'))
            or label.startswith(('MA$', 'MA_', 'mean M$', 'M$')))


def fold_plot_angles(samples, labels, reference=None, summary=None):
    """Center periodic columns on a circular-mean branch, including references.

    summary='mean'/'median' recomputes angle markers after folding: an old
    arithmetic marker at 180 degrees cannot simply be wrapped back to zero.
    No row reordering or folding of amplitudes, inclinations, h/k or other data.
    """
    values = np.array(samples, dtype=float, copy=True)
    marker = None if reference is None else np.array(reference, dtype=float, copy=True)
    for column, label in enumerate(labels):
        if not is_periodic_label(label):
            continue
        good = np.isfinite(values[:, column])
        if not good.any():
            continue
        angles = values[good, column] % 360.
        vector = np.mean(np.exp(1j*np.radians(angles)))
        # A uniform/antipodal posterior has no meaningful circular center.
        if abs(vector) < 1e-12:
            continue
        center = float(np.degrees(np.angle(vector)))
        values[good, column] = center + (angles-center+180.) % 360. - 180.
        if marker is not None:
            if summary in ('mean', 'median'):
                marker[column] = (np.mean if summary == 'mean' else np.median)(values[good, column])
            else:
                marker[column] = center + (marker[column]-center+180.) % 360. - 180.
    return values, marker
