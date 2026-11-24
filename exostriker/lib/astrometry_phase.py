"""Time coordinates for absolute-astrometry RA/Dec plots."""
import numpy as np


def astrometry_plot_times(data_times, model_times, period=None, epoch=0., phase_offset=0.,
                          phase_units="phase"):
    """Fold data and model onto [0, 1), sorting the connected model curve.

    Phase zero is the orbital element epoch plus the slider offset in cycles.
    Period units express the folded coordinates in days instead of cycles.
    Keep observations in their original order so uncertainties remain paired with their measurements.
    A missing period preserves the ordinary BJD view.
    """
    data_times = np.asarray(data_times)
    model_times = np.asarray(model_times)
    if period is None:
        return data_times, model_times, np.arange(len(model_times))
    if not np.isfinite(period) or period <= 0. or not np.isfinite(epoch):
        raise ValueError('Phase folding needs a finite positive period and epoch')

    data_phase = np.remainder(data_times - epoch - phase_offset * period, period) / period
    model_phase = np.remainder(model_times - epoch - phase_offset * period, period) / period
    if phase_units == 'period':
        data_phase = data_phase * period
        model_phase = model_phase * period
    elif phase_units != 'phase':
        raise ValueError('Folded coordinates must use phase or period units')
    model_order = np.argsort(model_phase, kind='stable')
    return data_phase, model_phase[model_order], model_order
