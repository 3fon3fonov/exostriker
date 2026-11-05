"""Transit dilution state, migration, and parameter validation.

D is the target's fraction of the total flux: F_observed = 1 + D*(F_transit-1).
The use flag selects fitting; a fixed value is applied even when it is false.
"""
import numpy as np


def ensure_dilution_parameters(fit):
    """Preserve legacy applied dilution, without turning Apply flags into fit flags."""
    legacy = getattr(fit, '_dilution_parameter_version', 0) < 1
    defaults = dict(tra_dil=1., tra_dil_err=[0., 0.], tra_dil_use=False,
                    tra_dil_bounds=[0., 1.], tra_dil_norm_pr=[1., .1, False],
                    tra_dil_jeff_pr=[.01, 1., False])
    for name, default in defaults.items():
        if not hasattr(fit, name):
            setattr(fit, name, {})
        values = getattr(fit, name)
        if not isinstance(values, dict):
            values = dict(enumerate(values))
            setattr(fit, name, values)
        for slot in range(60):
            if slot not in values:
                values[slot] = np.array(default) if isinstance(default, list) else default
    if not hasattr(fit, 'tra_dil_str'):
        fit.tra_dil_str = {i: r'tr. data dilution$_%s$' % (i+1) for i in range(60)}
    if legacy:
        datasets = getattr(fit, 'tra_data_sets', {})
        records = datasets.items() if hasattr(datasets, 'items') else enumerate(datasets)
        for slot, data in records:
            if len(data) > 8:
                fit.tra_dil[slot] = float(data[8])
        for slot in range(60):
            fit.tra_dil_use[slot] = False
        fit._dilution_parameter_version = 1


def validate_dilution(fit, slot):
    value = float(fit.tra_dil[slot])
    lo, hi = np.asarray(fit.tra_dil_bounds[slot], dtype=float)
    if not np.isfinite([value, lo, hi]).all() or not 0 <= lo < hi <= 1:
        raise ValueError(f'Transit dataset {slot+1}: dilution requires 0 <= lower bound < upper bound <= 1')
    if not lo <= value <= hi:
        raise ValueError(f'Transit dataset {slot+1}: dilution must lie within its bounds')
    mean, sigma, use_normal = fit.tra_dil_norm_pr[slot]
    alpha, beta, use_jeff = fit.tra_dil_jeff_pr[slot]
    if use_normal and use_jeff:
        raise ValueError(f'Transit dataset {slot+1}: select either a Gaussian or Jeffreys dilution prior')
    if use_normal and (not np.isfinite([mean, sigma]).all() or sigma <= 0):
        raise ValueError(f'Transit dataset {slot+1}: Gaussian dilution sigma must be positive')
    if use_jeff and (not np.isfinite([alpha, beta]).all() or not 0 < alpha < beta <= 1):
        raise ValueError(f'Transit dataset {slot+1}: Jeffreys dilution limits require 0 < alpha < beta <= 1')
    if use_jeff and (max(lo, alpha) >= min(hi, beta) or not alpha <= value <= beta):
        raise ValueError(f'Transit dataset {slot+1}: dilution must lie within the Jeffreys limits and bounds')
