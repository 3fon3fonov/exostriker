"""Systemic radial-velocity settings for Cartesian astrometric fits."""
import numpy as np


def ensure_velocity_settings(session):
    """Add fit settings to older sessions while retaining their fixed RV."""
    defaults = {
        'ast_radial_velocity': 0.,
        'astrometric_cartesian': False,
        'ast_radial_velocity_use': False,
        'ast_radial_velocity_bound': [-1000., 1000.],
        'ast_radial_velocity_norm_pr': [0., 100., False],
        'ast_radial_velocity_jeff_pr': [1., 1000., False],
        'ast_radial_velocity_err': [0., 0.],
        'ast_radial_velocity_index': None,
    }
    for name, value in defaults.items():
        if not hasattr(session, name):
            setattr(session, name, value.copy() if isinstance(value, list) else value)


def pack_velocity(session, parameters, flags, labels, bounds, gaussian, jeffreys):
    """Pack systemic RV before the established six-entry astrometric tail."""
    ensure_velocity_settings(session)
    session.ast_radial_velocity_index = None
    if not session.astrometric_cartesian:
        return
    session.ast_radial_velocity_index = len(parameters)
    parameters.append(session.ast_radial_velocity)
    active = session.type_fit['AST'] and session.use_ast_hipp_gaia
    flags.append(bool(active and session.ast_radial_velocity_use))
    labels.append(r'$v_{rad}$ [km/s]')
    bounds.append(session.ast_radial_velocity_bound)
    gaussian.append(session.ast_radial_velocity_norm_pr)
    jeffreys.append(session.ast_radial_velocity_jeff_pr)


def velocity_from_parameters(parameters, options):
    """Use the trial velocity when packed, otherwise the fixed session value."""
    index = options.get('ast_radial_velocity_index')
    if index is None:
        return float(options.get('ast_radial_velocity', 0.))
    return float(parameters[index])


def publish_velocity(session, parameters, errors=None):
    """Return fitted systemic RV and its uncertainty to the session."""
    ensure_velocity_settings(session)
    index = session.ast_radial_velocity_index
    if index is not None:
        session.ast_radial_velocity = float(parameters[index])
        if errors is not None:
            session.ast_radial_velocity_err = np.asarray(errors[index], dtype=float).copy()
