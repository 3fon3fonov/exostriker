"""Shared stellar density uses stellar bounds/priors, never scaled-axis limits."""
import copy
import numpy as np


def ensure_density_settings(session):
    for name, default in (
        ('stellar_rho_bounds', [100., 5000.]),
        ('stellar_rho_norm_pr', [1000., 100., False]),
        ('stellar_rho_jeff_pr', [100., 5000., False]),
        ('stellar_rho_str', r'$\rho_{\rm \star}$ [kg/m^3]'),
    ):
        if not hasattr(session, name):
            setattr(session, name, copy.deepcopy(default))


def density_limits(session):
    ensure_density_settings(session)
    lower, upper = map(float, session.stellar_rho_bounds)
    if not np.all(np.isfinite([lower, upper])) or lower <= 0 or lower >= upper:
        raise ValueError('Stellar density bounds must be positive and ordered (kg/m^3), under Limits and Priors > Stellar.')
    normal = session.stellar_rho_norm_pr
    jeffreys = session.stellar_rho_jeff_pr
    if normal[2] and jeffreys[2]:
        raise ValueError('Select only one stellar density prior family under Limits and Priors > Stellar.')
    if normal[2] and (not np.all(np.isfinite(normal[:2])) or normal[1] <= 0):
        raise ValueError('Stellar density Gaussian prior requires a finite mean and positive sigma.')
    if jeffreys[2] and (not np.all(np.isfinite(jeffreys[:2])) or not 0 < jeffreys[0] < jeffreys[1]):
        raise ValueError('Stellar density Jeffreys prior requires positive ordered limits.')
    return session.stellar_rho_bounds


def set_density_mode(session, enabled):
    """Keep density and a/Rstar values distinct when switching GUI coordinates."""
    enabled = bool(enabled)
    if enabled:
        session.stellar_mass_radius = False
    if enabled == session.rho:
        return
    slots = list(np.flatnonzero(session.use_planet))
    if enabled:
        session._scaled_axes_before_density = copy.deepcopy(session.pl_a)
        for slot in slots:
            session.pl_a[slot] = float(session.stellar_rho)
    else:
        saved = getattr(session, '_scaled_axes_before_density', {})
        density = session.pl_a[slots[0]] if slots else session.stellar_rho
        for slot in slots:
            if slot in saved:
                session.pl_a[slot] = saved[slot]
            else:
                # Old density sessions have no saved scaled-axis coordinate.
                session.pl_a[slot] = (6.67408e-11*density*(session.P[slot]*86400.)**2/(3*np.pi))**(1/3)
    session.rho = enabled
