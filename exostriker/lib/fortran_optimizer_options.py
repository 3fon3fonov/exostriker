"""Session-owned Fortran solver controls. Zero retains each backend's old default."""
import copy
import math
import numpy as np

FIELDS = ('maxiter', 'maxfev', 'ftol', 'xtol', 'gtol', 'step_scale',
          'initial_damping', 'diff_step_scale', 'restart_ftol')
DEFAULTS = {
    'simplex': dict(maxiter=0, maxfev=0, ftol=0., xtol=0., gtol=0.,
                    step_scale=1., initial_damping=.001, diff_step_scale=1.,
                    restart_ftol=1e-6, restarts=0),
    'lm': dict(maxiter=0, maxfev=0, ftol=0., xtol=0., gtol=0.,
               step_scale=1., initial_damping=.001, diff_step_scale=1.,
               restart_ftol=1e-6, restarts=0),
}


def ensure_options(session):
    """Fill old sessions without discarding configured values; reject invalid input."""
    stored = getattr(session, 'fortran_optimizer_options', {})
    if not isinstance(stored, dict):
        raise ValueError('fortran_optimizer_options must be a dictionary')
    result = copy.deepcopy(DEFAULTS)
    for method, defaults in DEFAULTS.items():
        values = stored.get(method, {})
        unknown = set(values) - set(defaults)
        if unknown:
            raise ValueError(f'Unknown {method} Fortran options: {sorted(unknown)}')
        for key, value in values.items():
            if isinstance(value, (bool, np.bool_)) or not math.isfinite(float(value)):
                raise ValueError(f'{method}.{key} must be finite and numeric')
            if key in ('maxiter', 'maxfev', 'restarts'):
                if int(value) != value or not 0 <= value <= 100000000:
                    raise ValueError(f'{method}.{key} must be an integer from 0 to 100000000')
                value = int(value)
            elif value < 0 or (key in ('step_scale', 'initial_damping', 'diff_step_scale', 'restart_ftol') and value == 0):
                raise ValueError(f'{method}.{key} must be positive (tolerances may be zero for automatic)')
            result[method][key] = value
    session.fortran_optimizer_options = result
    return result


def option_vector(session, method):
    values = ensure_options(session)[method]
    return np.array([values[key] for key in FIELDS], dtype=float)
