"""Alternative mass coordinates, with K retained for existing model interfaces.

Masses are in Jupiter masses. Dynamical masses are physical Jacobi masses;
Keplerian masses are exact edge-on minimum masses (usually called m sin i).
"""
import copy
import numpy as np
from RV_mod.jacobi import GMSUN, TWOPI, mass_a_from_jacobi_el

MSUN_TO_MJUP = 1047.5654817267318


def amplitude_root(session):
    """Name the independently fitted amplitude coordinate."""
    return 'mp' if getattr(session, 'mass_parameterization', False) else 'K'


def ensure_mass_settings(session):
    """Upgrade old sessions without changing their default K parameterization."""
    if not hasattr(session, 'mass_parameterization'):
        session.mass_parameterization = False
    defaults = {
        'mp': 0.,
        'mp_use': True,
        'mp_bound': [0., 1000.],
        'mp_norm_pr': [0., 1., False],
        'mp_jeff_pr': [0.001, 1000., False],
        'mp_err': [0., 0.],
    }
    for name, value in defaults.items():
        if not hasattr(session, name):
            setattr(session, name, {i: copy.deepcopy(value) for i in range(9)})
    session.mp_str = {i: r'$m_%s$ [M$_J$]' % chr(98+i) for i in range(9)}
    if not session.mod_dynamical:
        session.mp_str = {i: r'$m_{%s,\min}$ [M$_J$]' % chr(98+i) for i in range(9)}


def masses_to_k(rows, masses, stellar_mass, dynamical=True):
    """Convert a full period-ordered Jacobi system to equivalent RV amplitudes."""
    rows = np.asarray(rows, float).reshape(-1, 7)
    masses = np.asarray(masses, float) / MSUN_TO_MJUP
    if masses.shape != (len(rows),):
        raise ValueError('Provide one mass for each orbital row.')
    if not np.isfinite(stellar_mass) or stellar_mass <= 0:
        raise ValueError('Stellar mass must be finite and positive.')
    if not np.all(np.isfinite(rows)) or not np.all(np.isfinite(masses)) or np.any(masses < 0):
        raise ValueError('Orbital parameters and nonnegative masses must be finite.')
    if np.any(rows[:, 1] <= 0) or np.any(rows[:, 2] < 0) or np.any(rows[:, 2] >= 1):
        raise ValueError('Mass conversion needs positive periods and bound orbits.')
    amplitudes = np.zeros(len(rows))
    inner = stellar_mass
    for index in np.argsort(rows[:, 1], kind='stable'):
        mass = masses[index]
        projection = abs(np.sin(np.radians(rows[index, 5]))) if dynamical else 1.
        coefficient = (TWOPI * GMSUN / (rows[index, 1] * 86400.)) ** (1./3.)
        amplitudes[index] = coefficient * mass * projection / (
            (inner + mass) ** (2./3.) * np.sqrt(1.-rows[index, 2]**2))
        inner += mass
    return amplitudes


def session_rows(session, slots):
    """Build classical orbital rows, respecting h/k and coplanar modes."""
    names = ('K', 'P', 'e', 'w', 'M0', 'i', 'Node')
    rows = np.array([[getattr(session, name)[i] for name in names] for i in slots], float).reshape(-1, 7)
    if session.hkl and len(slots):
        h = np.array([session.e_sinw[i] for i in slots])
        k = np.array([session.e_cosw[i] for i in slots])
        rows[:, 2] = np.hypot(h, k)
    if session.mod_dynamical and session.copl_incl and len(slots):
        rows[:, 5] = rows[0, 5]
    return rows


def sync_amplitudes(session):
    """Update the derived coordinate, preserving the independently fitted one."""
    mass_mode = getattr(session, 'mass_parameterization', False)
    if not mass_mode and not hasattr(session, 'use_planet'):
        return
    ensure_mass_settings(session)
    mixed = bool(np.asarray(getattr(session, 'mixed_fit', [False])[0]).item())
    if mass_mode and mixed:
        raise ValueError('Mass parameterization is not yet supported for mixed dynamical subsets.')
    slots = list(np.flatnonzero(session.use_planet))
    rows = session_rows(session, slots)
    if mass_mode:
        amplitudes = masses_to_k(
            rows, [session.mp[i] for i in slots],
            session.params.stellar_mass, session.mod_dynamical)
        for index, slot in enumerate(slots):
            session.K[slot] = amplitudes[index]
    else:
        if not session.mod_dynamical:
            rows[:, 5] = 90.
        try:
            masses, _ = mass_a_from_jacobi_el(
                rows, len(slots), session.params.stellar_mass)
        except ValueError:
            # Invalid trial inputs must not change the existing fit's error handling.
            masses = np.full(len(slots), np.nan)
        for index, slot in enumerate(slots):
            session.mp[slot] = masses[index]


def set_mass_mode(session, enabled):
    """Switch coordinates while preserving the current physical orbit."""
    ensure_mass_settings(session)
    enabled = bool(enabled)
    if enabled and bool(np.asarray(getattr(session, 'mixed_fit', [False])[0]).item()):
        raise ValueError('Mass parameterization is not yet supported for mixed dynamical subsets.')
    if enabled == session.mass_parameterization:
        return
    slots = list(np.flatnonzero(session.use_planet))
    if enabled:
        rows = session_rows(session, slots)
        if np.any(rows[:, 0] < 0):
            raise ValueError('Switching to mass requires nonnegative RV amplitudes.')
        if not session.mod_dynamical:
            rows[:, 5] = 90.
        masses, _ = mass_a_from_jacobi_el(rows, len(slots), session.params.stellar_mass)
        for index, slot in enumerate(slots):
            session.mp[slot] = masses[index]
        session.mass_settings_initialized = True
    else:
        sync_amplitudes(session)
    # Preserve the planet's fit selection when changing the fitted coordinate.
    for slot in range(9):
        selected = session.mp_use[slot] if session.mass_parameterization else session.K_use[slot]
        session.K_use[slot] = bool(selected)
        session.mp_use[slot] = bool(selected)
    session.mass_parameterization = enabled


def configure_mass_options(session, options):
    """Carry coordinate conventions into legacy and native inference contexts."""
    if getattr(session, 'mass_parameterization', False):
        if bool(np.asarray(getattr(session, 'mixed_fit', [False])[0]).item()):
            raise ValueError('Mass parameterization is not yet supported for mixed dynamical subsets.')
    options['mass_parameterization'] = getattr(session, 'mass_parameterization', False)
    options['mass_dynamical'] = session.mod_dynamical


def vector_to_k(parameters, npl, ndset, options):
    """Return a model-only K vector; the sampler vector retains mass coordinates."""
    result = np.asarray(parameters, float).copy()
    rows = result[2*ndset:2*ndset+7*npl].reshape(npl, 7)
    classical = rows.copy()
    if options['hkl']:
        classical[:, 2] = np.hypot(rows[:, 2], rows[:, 3])
    if options.get('copl_incl') and options['mass_dynamical'] and npl:
        classical[:, 5] = classical[0, 5]
    rows[:, 0] = masses_to_k(classical, rows[:, 0], result[-1], options['mass_dynamical'])
    return result


def model_vector(session):
    """Return K-based parameters for legacy plotting functions."""
    if not getattr(session, 'mass_parameterization', False):
        return session.parameters
    ndset = len(getattr(session, 'mcmc_dataset_slots', ()))
    options = dict(hkl=session.hkl, copl_incl=session.copl_incl,
                   mass_dynamical=session.mod_dynamical)
    return vector_to_k(session.parameters, session.npl, ndset, options)


def mass_model_session(session):
    """Refresh derived amplitudes on a shallow copy for read-only model APIs."""
    if not getattr(session, 'mass_parameterization', False):
        return session
    model_session = copy.copy(session)
    model_session.K = copy.copy(session.K)
    sync_amplitudes(model_session)
    return model_session
