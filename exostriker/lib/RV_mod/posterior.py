"""Posterior summaries for linear parameters and orbital angles in degrees."""
import numpy as np


def angle_columns(flags, ndset, npl, hkl=False):
    """Inclination is bounded, not periodic; h/k are linear, lambda is periodic."""
    return [column for column, index in enumerate(flags)
            if 0 <= index - 2 * ndset < 7 * npl
            and ((index - 2 * ndset) % 7 in (4, 6)
                 or ((index - 2 * ndset) % 7 == 3 and not hkl))]


def circular_mean(angles):
    radians = np.radians(np.asarray(angles) % 360.)
    vector = np.mean(np.exp(1j * radians))
    # An antipodal/uniform posterior has no preferred direction. Use a
    # deterministic sample as the reference for its branch-cut summaries.
    center = angles[0] if abs(vector) < 1e-12 else np.degrees(np.angle(vector))
    return float(center % 360.)


def centered_samples(samples, flags, ndset, npl, hkl=False):
    """Return a copy with each periodic column on a common local branch."""
    values = np.array(samples, dtype=float, copy=True)
    for column in angle_columns(flags, ndset, npl, hkl):
        center = circular_mean(values[:, column])
        values[:, column] = center + (values[:, column] - center + 180.) % 360. - 180.
    return values


def posterior_summary(samples, flags, ndset, npl, hkl=False):
    values = centered_samples(samples, flags, ndset, npl, hkl)
    means = np.mean(values, axis=0)
    medians = np.median(values, axis=0)
    for column in angle_columns(flags, ndset, npl, hkl):
        means[column] = circular_mean(np.asarray(samples)[:, column])
        medians[column] %= 360.
    return means, medians


def posterior_errors(samples, best, level, flags, ndset, npl, hkl=False):
    values = centered_samples(samples, flags, ndset, npl, hkl)
    reference = np.array(best, dtype=float, copy=True)
    for column in angle_columns(flags, ndset, npl, hkl):
        center = circular_mean(np.asarray(samples)[:, column])
        reference[column] = center + (reference[column] - center + 180.) % 360. - 180.
    lower, upper = np.percentile(values, [level, 100. - level], axis=0)
    return np.column_stack((reference - lower, upper - reference)).tolist()
