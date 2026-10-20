"""Prepare Gaia along-scan residuals once, retaining the original measurements.

The input is the seven-row HIP-like array returned by ExTRA.gaia_init.
Exo-Striker's Gaia importer expects centroids that still include the catalogue
single-star solution. Only parallax and proper motion are removed here: the
catalogue RA/Dec specify the reference position, not offsets in mas.
"""
import numpy as np

_METADATA_KEY = '_exostriker_gaia_baseline_v1'
_STANDARD_KEYS = ('RAdeg', 'DEdeg', 'Plx', 'pm_RA', 'pm_DE')


def prepare_gaia_dataset(dataset):
    """Return catalogue residuals and migrate a legacy five-field record once.

    The record layout and filename at [-1] stay unchanged. Private metadata
    lives in the header dictionary and survives session serialization. Retain
    raw along-scan measurements so changing the catalogue header can safely
    recompute the residuals instead of subtracting twice.
    """
    data, header = dataset[0], dataset[3]
    solution = tuple(float(header[key]) for key in _STANDARD_KEYS)
    metadata = header.get(_METADATA_KEY)
    if metadata is not None and metadata['solution'] == solution:
        return data
    if not np.all(np.isfinite(solution)):
        raise ValueError('Gaia catalogue position, parallax and proper motion must be finite')

    data = np.asarray(data, dtype=float)
    if data.ndim != 2 or data.shape[0] != 7 or data.shape[1] == 0:
        raise ValueError('Gaia data must have seven rows and at least one observation')
    if not np.all(np.isfinite(data)) or np.any(data[-1] <= 0):
        raise ValueError('Gaia observations must be finite with positive uncertainties')
    raw_al = data[-2].copy() if metadata is None else np.array(metadata['raw_al'], copy=True)
    if raw_al.shape != data[-2].shape:
        raise ValueError('Gaia cached measurements do not match the observation count')

    # Same subtraction/order as abs_res(raw, catalogue, [RA, Dec, 0, 0, 0], data).
    residual = raw_al.copy()
    for row, parameter in zip(data[2:5], solution[2:]):
        residual -= row * parameter
    prepared = data.copy()
    prepared[-2] = residual
    prepared_header = header.copy()
    prepared_header[_METADATA_KEY] = {'solution': solution, 'raw_al': raw_al}
    # Publish only after validation and calculation, keeping a failed import
    # or header edit from leaving a partially prepared record.
    dataset[0], dataset[3] = prepared, prepared_header
    return prepared
