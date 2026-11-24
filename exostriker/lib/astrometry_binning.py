"""Plot-only binning of absolute astrometric measurements."""
import numpy as np


def bin_astrometry(coordinates, sky_data, along_scan_errors, width, origin=0.):
    """Average displayed points in fixed [origin+k*w, origin+(k+1)*w) bins.

    Use the along-scan inverse variances for both sky coordinates, preserving
    their common measurement weights. Propagate the projected uncertainties
    through the weighted mean rather than treating small projected errors as
    independent high-precision measurements. Input observations are unchanged.
    """
    coordinates = np.asarray(coordinates, dtype=float)
    sky_data = np.asarray(sky_data, dtype=float)
    errors = np.asarray(along_scan_errors, dtype=float)
    if not np.isfinite(width) or width <= 0. or not np.isfinite(origin):
        raise ValueError('Binning needs a finite positive width and origin')
    if (coordinates.ndim != 1 or sky_data.shape != (4, coordinates.size)
            or errors.shape != coordinates.shape):
        raise ValueError('Astrometric bin coordinates, errors and points must align')
    if (not np.all(np.isfinite(coordinates)) or not np.all(np.isfinite(sky_data))
            or not np.all(np.isfinite(errors)) or np.any(errors <= 0.)):
        raise ValueError('Astrometric bins need finite measurements and positive errors')
    if not coordinates.size:
        return np.empty(0), np.empty((4, 0)), np.empty(0, dtype=int)

    bin_ids = np.floor((coordinates - origin) / width)
    bin_coordinates = []
    bin_points = []
    counts = []
    for bin_id in np.unique(bin_ids):
        selected = bin_ids == bin_id
        selected_errors = errors[selected]
        weights = (np.min(selected_errors) / selected_errors) ** 2
        weights /= np.sum(weights)
        points = sky_data[:, selected]
        mean_x = np.sum(weights * points[0])
        mean_y = np.sum(weights * points[2])
        error_x = np.sqrt(np.sum((weights * points[1]) ** 2))
        error_y = np.sqrt(np.sum((weights * points[3]) ** 2))
        bin_coordinates.append(np.sum(weights * coordinates[selected]))
        bin_points.append([mean_x, error_x, mean_y, error_y])
        counts.append(np.count_nonzero(selected))
    return np.asarray(bin_coordinates), np.asarray(bin_points).T, np.asarray(counts)
