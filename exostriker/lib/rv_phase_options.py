"""Display-only RV refolding from cached, absolute-time model samples."""
import numpy as np
from PyQt6 import QtWidgets


def raw_signal(fit, slot):
    slots = sorted((i for i in fit.P if fit.use_planet[i]), key=fit.P.get)
    rank = slots.index(slot)
    result = fit.fit_results
    return (np.asarray(result.model_jd), np.asarray(result.res['fit']).T[2 + rank],
            np.asarray(result.model_data)[6 + rank])


def fold_cached(fit, slot, period):
    """Keep RV decomposition, errors, GP subtraction, and residuals unchanged."""
    if not np.isfinite(period) or period <= 0:
        raise ValueError('Folding period must be positive.')
    times, model, component = raw_signal(fit, slot)
    result = fit.fit_results
    jd = np.asarray(result.rv_model.jd)
    origin = jd[0]
    mx = (times - origin) % period
    dx = (jd - origin) % period
    ms, ds = np.argsort(mx), np.argsort(dx)
    residual = np.array(result.rv_model.o_c, copy=True)
    if fit.doGP:
        residual -= np.asarray(fit.gp_model_data[0])
    return ([dx[ds], (component + residual)[ds], np.asarray(result.rv_model.rv_err)[ds],
             np.asarray(result.idset)[ds], residual[ds]], [mx[ms], model[ms]])


class RVPhaseOptions(QtWidgets.QDialog):
    def __init__(self, gui, session):
        super().__init__(gui)
        self.gui, self.session = gui, session
        self.settings = {}
        self.setWindowTitle('RV phase options')
        form = QtWidgets.QFormLayout(self)
        self.info = QtWidgets.QLabel()
        self.info.setWordWrap(True)
        form.addRow(self.info)
        gui.extra_plot_RVphase_norm.setParent(self)
        gui.extra_plot_RVphase_norm.setText('Normalize phase to −0.5 … +0.5')
        gui.extra_plot_RVphase_norm.show()
        form.addRow(gui.extra_plot_RVphase_norm)
        self.mean = QtWidgets.QCheckBox('Use mean osculating period from integration')
        self.mean.setToolTip('Time-weighted mean of P(t) from the fitted system’s cached orbital integration.\nRun orb. evol after fitting; no integration is launched by this option.')
        form.addRow(self.mean)
        self.delta = QtWidgets.QDoubleSpinBox()
        self.delta.setDecimals(8)
        self.delta.setRange(-1.e7, 1.e7)
        self.delta.setSingleStep(1.)
        form.addRow('Period adjustment ± [days]', self.delta)
        self.points = QtWidgets.QCheckBox('Show model as points (useful for N-body spread)')
        form.addRow(self.points)
        self.status = QtWidgets.QLabel()
        form.addRow(self.status)
        reset = QtWidgets.QPushButton('Reset this planet to fitted period')
        form.addRow(reset)
        reset.clicked.connect(self.reset)
        for widget, signal in ((self.mean, 'toggled'), (self.delta, 'valueChanged'), (self.points, 'toggled')):
            getattr(widget, signal).connect(self.changed)
        self.slot = None
        self.signature = None

    def configure(self, slot):
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
        self.info.setText(f'Planet {slot+1}: fitted period {self.nominal:.9g} days.\nMean osculating period: {estimate}.\n{integration_info}\nDisplay only; uses cached model times and leaves the fit unchanged.')
        self.period()

    def period(self):
        value = (self.estimate if self.mean.isChecked() and self.estimate is not None else self.nominal) + self.delta.value()
        self.status.setText(f'Folding period: {value:.9g} days' if value > 0 else 'Period must be positive; showing fitted period.')
        return value if value > 0 else self.nominal

    def changed(self, *args):
        if self.slot is None:
            return
        self.settings[self.slot] = (self.mean.isChecked(), self.delta.value(), self.points.isChecked())
        self.gui.phase_plots(self.slot + 1)

    def reset(self):
        self.settings.pop(self.slot, None)
        self.configure(self.slot)
        self.gui.phase_plots(self.slot + 1)
