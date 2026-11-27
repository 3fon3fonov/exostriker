"""Catalogue residuals and plots for the classical astrometric likelihoods."""
import numpy as np


def standard_astrometry_iad(data, header, gaia=False):
    """Subtract the linear catalogue baseline from prepared absolute centroids.

    Work on a copy: rigorous likelihoods and future session saves still need
    the original centroids. Position offsets are in mas, never absolute degrees.
    """
    residual_iad = np.asarray(data, dtype=float).copy()
    if gaia and header.get('ALModelRemoved', 0.):
        return residual_iad
    alpha = header.get('delta_RA', header.get('ALOffsetRA', 0.)) if gaia else 0.
    delta = header.get('deltaDE', header.get('ALOffsetDE', 0.)) if gaia else 0.
    baseline = np.array([
        alpha, delta, header['Plx'], header['pm_RA'], header['pm_DE'],
    ], dtype=float)
    if gaia:
        # Gaia IAD time derivatives always refer to J2017.5.
        elapsed = 2017.5 - header.get('RefEpoch', 2017.5)
        baseline[:2] += elapsed * baseline[3:]
    residual_iad[-2] -= baseline @ residual_iad[:5]
    return residual_iad


def standard_astrometry_plot(ex, data, correction, planets, timestamps,
                             model_max, model_samp, plot_planet=None):
    """Publish sky-plane points using the same residuals as L_hip/L_gaia."""
    model_start = min(timestamps)
    model_end = max(timestamps) + model_max
    if plot_planet is not None:
        model_end = model_start + planets[plot_planet, 0]
    model_times = np.linspace(model_start, model_end, model_samp)
    orbit_x, orbit_y = ex.orbit_total(planets, timestamps)
    corrected = ex.abs_res(data[-2], correction, np.zeros(5), data)
    residuals = corrected - data[0] * orbit_x - data[1] * orbit_y

    plot_iad = data.copy()
    plot_iad[-2] = residuals
    residual_x, error_x, residual_y, error_y = ex.hip_2d(plot_iad)
    displayed_planets = planets
    if plot_planet is not None:
        displayed_planets = planets[plot_planet:plot_planet + 1]
        selected_x, selected_y = ex.orbit_total(displayed_planets, timestamps)
        corrected -= data[0] * (orbit_x - selected_x)
        corrected -= data[1] * (orbit_y - selected_y)
        orbit_x, orbit_y = selected_x, selected_y
    plot_iad[-2] = corrected
    observed = ex.hip_2d(plot_iad)
    observations_on_orbit = (
        residual_x + orbit_x, error_x, residual_y + orbit_y, error_y,
    )
    sky_residuals = (residual_x, error_x, residual_y, error_y)
    plot_data = np.empty(4, dtype=object)
    plot_data[:] = [observed, observations_on_orbit, timestamps, sky_residuals]
    model_x, model_y = ex.orbit_total(displayed_planets, model_times)
    model_x = np.broadcast_to(model_x, model_times.shape)
    model_y = np.broadcast_to(model_y, model_times.shape)
    return plot_data, np.array([model_x, model_y, model_times])
