"""Provenance and time averages for cached osculating orbital periods."""
import numpy as np


def orbit_signature(obj):
    slots = sorted((i for i in range(len(obj.use_planet)) if obj.use_planet[i]),
                   key=lambda i: obj.P[i])
    return (float(obj.params.stellar_mass), float(getattr(obj, 'epoch', 0.)),
            tuple((i, *(float(getattr(obj, name)[i])
                        for name in ('K', 'P', 'e', 'i', 'Node', 'w', 'M0')))
                  for i in slots))


def record_period_source(obj, slots):
    # Arbitrary-orbit adapter objects do not contain fitted RV parameters.
    if hasattr(obj, 'K'):
        obj.phase_period_source = dict(signature=orbit_signature(obj), slots=list(slots))


def mean_osculating_period(obj, slot):
    """Time-weighted mean P(t) in days over the available integration (years).

    Evolution outputs use period-ordered integration indices, not GUI slots.
    Reject stale/unattributed caches and arbitrary-planet runs.
    """
    source = getattr(obj, 'phase_period_source', None)
    if not source:
        return None, 'Run orb. evol for the fitted system to obtain osculating periods.'
    if source['signature'] != orbit_signature(obj):
        return None, 'Orbital parameters changed; rerun orb. evol for the current fit.'
    try:
        rank = source['slots'].index(slot)
        times = np.asarray(obj.evol_T[rank], dtype=float)
        periods = np.asarray(obj.evol_Per[rank], dtype=float)
    except (AttributeError, KeyError, IndexError, ValueError):
        return None, 'No integration periods available for this planet.'
    if times.ndim != 1 or periods.shape != times.shape or len(times) < 2:
        return None, 'At least two integration samples are required.'
    if not (np.all(np.isfinite(times)) and np.all(np.isfinite(periods)) and np.all(periods > 0)):
        return None, 'Integration contains invalid or unbound orbital periods.'
    order = np.argsort(times)
    times, periods = times[order], periods[order]
    if np.any(np.diff(times) <= 0):
        return None, 'Integration times must be distinct.'
    trapezoid = getattr(np, 'trapezoid', None) or np.trapz
    mean = float(trapezoid(periods, times) / (times[-1] - times[0]))
    return mean, f'Time-weighted over integration t = {times[0]:.8g} … {times[-1]:.8g} yr from its starting epoch.'
