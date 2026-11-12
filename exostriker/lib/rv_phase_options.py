"""Display-only RV refolding from cached, absolute-time model samples."""
import numpy as np
from collections.abc import Mapping
from PyQt6 import QtWidgets


def raw_signal(fit, slot):
    """Read one planet component from the cached absolute-time RV model."""

    slots = sorted((i for i in fit.P if fit.use_planet[i]), key=fit.P.get)
    rank = slots.index(slot)
    result = fit.fit_results
    saved = getattr(result, 'res', None)
    if not isinstance(saved, Mapping) or saved.get('fit') is None:
        raise ValueError('Original RV model samples are unavailable.')
    times = np.asarray(result.model_jd)
    samples = np.asarray(saved['fit'])
    data = np.asarray(result.model_data)
    if (times.ndim != 1 or not times.size or samples.ndim != 2
            or samples.shape[0] != times.size or samples.shape[1] <= 2 + rank
            or data.ndim != 2 or data.shape[0] <= 6 + rank):
        raise ValueError('Original RV model samples are unavailable or incomplete.')
    return (times, samples[:, 2 + rank], data[6 + rank])


def fold_cached(fit, slot, period):
    """Keep RV decomposition, errors, GP subtraction, and residuals unchanged."""
    if not np.isfinite(period) or period <= 0:
        raise ValueError('Folding period must be positive.')
    times, model, component = raw_signal(fit, slot)
    result = fit.fit_results
    jd = np.asarray(result.rv_model.jd)
    if jd.ndim != 1 or not jd.size:
        raise ValueError('No RV observations are available for phase folding.')
    origin = jd[0]
    mx = (times - origin) % period
    dx = (jd - origin) % period
    ms, ds = (np.argsort(mx), np.argsort(dx))
    residual = np.array(result.rv_model.o_c, copy=True)
    if fit.doGP:
        residual -= np.asarray(fit.gp_model_data[0])
    return ([dx[ds], (component + residual)[ds], np.asarray(result.rv_model.rv_err)[ds],
             np.asarray(result.idset)[ds], residual[ds]], [mx[ms], model[ms]])


class RVPhaseOptions(QtWidgets.QDialog):
    """Keep display-only folding choices separate from fitted periods."""

    def __init__(self, gui, session):
        super().__init__(gui)
        self.gui, self.session = (gui, session)
        self.settings = {}
        from rv_phase_options_ui import Ui_RVPhaseOptions
        ui = Ui_RVPhaseOptions()
        ui.setupUi(self)
        self.info, self.mean, self.delta = (ui.info, ui.mean, ui.delta)
        self.points, self.status = (ui.points, ui.status)
        gui.extra_plot_RVphase_norm.hide()
        gui.extra_plot_RVphase_norm = ui.normalize
        ui.reset.clicked.connect(self.reset)
        for widget, signal in (
            (self.mean, 'toggled'),
            (self.delta, 'valueChanged'),
            (self.points, 'toggled'),
        ):
            getattr(widget, signal).connect(self.changed)
        self.slot = None
        self.signature = None

    def configure(self, slot):
        """Load one planet and reset choices when the fitted model changes."""

        from PyQt6 import QtCore
        fit = self.session()
        signature = (id(fit), id(fit.fit_results), tuple(fit.P.items()))
        if signature != self.signature:
            self.signature = signature
            self.settings.clear()
        self.slot = slot
        self.nominal = float(fit.P[slot])
        from integration_period import mean_osculating_period
        self.estimate, integration_info = mean_osculating_period(fit, slot)
        mean, delta, points = self.settings.get(slot, (False, 0., False))
        blockers = [QtCore.QSignalBlocker(w) for w in (self.mean, self.delta, self.points)]
        self.mean.setEnabled(self.estimate is not None)
        self.mean.setChecked(mean and self.estimate is not None)
        self.delta.setValue(delta)
        self.points.setChecked(points)
        del blockers
        estimate = f'{self.estimate:.9g} days' if self.estimate is not None else 'unavailable'
        self.info.setText(f'Planet {slot + 1}: fitted period {self.nominal:.9g} days.\nMean osculating period: {estimate}.\n{integration_info}\nDisplay only; uses cached model times and leaves the fit unchanged.')
        self.period()

    def period(self):
        """Return a positive display period without changing the fit."""

        value = (self.estimate if self.mean.isChecked() and self.estimate is not None else self.nominal) + self.delta.value()
        self.status.setText(f'Folding period: {value:.9g} days' if value > 0 else 'Period must be positive; showing fitted period.')
        return value if value > 0 else self.nominal

    def changed(self, *args):
        """Save the display settings and redraw the selected phase plot."""

        if self.slot is None:
            return
        self.settings[self.slot] = (self.mean.isChecked(), self.delta.value(), self.points.isChecked())
        self.gui.phase_plots(self.slot + 1)

    def reset(self):
        """Restore the fitted-period display for the current planet."""

        self.settings.pop(self.slot, None)
        self.configure(self.slot)
        self.gui.phase_plots(self.slot + 1)
