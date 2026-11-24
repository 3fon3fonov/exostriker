"""Restore Hipparcos catalogue motion for rigorous centroid likelihoods."""
import numpy as np

_CENTROID_KEY = '_exostriker_hip_centroids_v1'


def prepare_hipparcos_dataset(dataset):
    """Convert imported or saved catalogue residuals to centroids once.

    Hipparcos IAD gives residuals about its catalogue solution. Rigorous
    ExTRA subtracts the full propagated solution, so restore parallax and
    proper motion first. Absolute RA/Dec define the comparison point and
    must never be added as milliarcsecond offsets.
    """
    data, header = dataset[0], dataset[3]
    data = np.asarray(data, dtype=float)
    if data.ndim != 2 or data.shape[0] != 7 or data.shape[1] == 0:
        raise ValueError('Hipparcos data must have seven rows and at least one observation')
    if not np.all(np.isfinite(data)) or np.any(data[-1] <= 0.):
        raise ValueError('Hipparcos observations must be finite with positive uncertainties')
    catalogue = np.array([
        header['RAdeg'], header['DEdeg'], header['Plx'],
        header['pm_RA'], header['pm_DE'],
    ], dtype=float)
    if not np.all(np.isfinite(catalogue)):
        raise ValueError('Hipparcos catalogue solution must be finite')
    if header.get(_CENTROID_KEY, False):
        return data

    centroids = data.copy()
    for derivative, parameter in zip(data[2:5], catalogue[2:]):
        centroids[-2] += derivative * parameter
    prepared_header = header.copy()
    prepared_header[_CENTROID_KEY] = True
    # Publish the measurements and marker together, after validation.
    dataset[0], dataset[3] = centroids, prepared_header
    return centroids
