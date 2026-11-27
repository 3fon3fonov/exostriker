"""Validate Gaia measurements for rigorous propagation of the full solution."""
import numpy as np

_METADATA_KEY = '_exostriker_gaia_baseline_v1'
_STANDARD_KEYS = ('RAdeg', 'DEdeg', 'Plx', 'pm_RA', 'pm_DE')


def prepare_gaia_dataset(dataset, allow_residuals=False):
    """Return raw along-scan measurements, restoring old cached residuals.

    Rigorous ExTRA propagates the catalogue solution itself. Older sessions
    may store subtracted measurements; their cache retains the original data.
    """
    data, header = dataset[0], dataset[3]
    if header.get('ALModelRemoved', 0.) and not allow_residuals:
        raise ValueError('Rigorous Gaia fitting needs original centroids; export with -m n')
    solution = tuple(float(header[key]) for key in _STANDARD_KEYS)
    offsets = (
        float(header.get('delta_RA', header.get('ALOffsetRA', 0.))),
        float(header.get('deltaDE', header.get('ALOffsetDE', 0.))),
    )
    reference_epoch = float(header.get('RefEpoch', 2017.5))
    if not np.all(np.isfinite((*solution, *offsets, reference_epoch))):
        raise ValueError('Gaia solution and comparison coordinates must be finite')
    if not 0 <= solution[0] < 360 or not -90 <= solution[1] <= 90:
        raise ValueError('Gaia comparison coordinates must be in degrees')

    data = np.asarray(data, dtype=float)
    if data.ndim != 2 or data.shape[0] != 7 or data.shape[1] == 0:
        raise ValueError('Gaia data must have seven rows and at least one observation')
    if not np.all(np.isfinite(data)) or np.any(data[-1] <= 0):
        raise ValueError('Gaia observations must be finite with positive uncertainties')

    metadata = header.get(_METADATA_KEY)
    if metadata is not None:
        raw_al = np.asarray(metadata['raw_al'], dtype=float)
        if raw_al.shape != data[-2].shape or not np.all(np.isfinite(raw_al)):
            raise ValueError('Gaia cached measurements do not match the observations')
        data = data.copy()
        data[-2] = raw_al
        header = header.copy()
        del header[_METADATA_KEY]
        dataset[0], dataset[3] = data, header
    return data
