"""Independent phased transit extra plot using cached fit results."""
import numpy as np
from PyQt6 import QtCore, QtWidgets
import pyqtgraph as pg


def transit_phase(times, period, epoch, shift=0., normalized=False):
    """Center on the fitted transit epoch, preserving a common phase scale."""
    phase = (np.asarray(times, dtype=float) - epoch - shift + period/2) % period - period/2
    return phase / period if normalized else phase


def transit_series(fit, dataset, subtract_gp=True):
    """Return independent arrays; GP = residual before GP minus residual after GP."""
    sep = fit.transit_results[1]
    t, flux, error, model, residual = [np.array(sep[k][dataset], dtype=float, copy=True)
                                      for k in range(5)]
    if any(values.shape != t.shape for values in (flux, error, model, residual)):
        raise ValueError('Transit result arrays have inconsistent lengths')
    if subtract_gp and fit.tra_doGP:
        corrected = np.asarray(sep[5][dataset], dtype=float)
        if corrected.shape != residual.shape:
            raise ValueError('Transit GP residual arrays have inconsistent lengths')
        flux -= residual - corrected
        residual = corrected.copy()
    elif fit.tra_doGP:
        # Show the full transit+GP model when GP remains in the data.
        corrected = np.asarray(sep[5][dataset], dtype=float)
        model += residual - corrected
        residual = corrected.copy()
    return t, flux, error, model, residual


class TransitPhaseView(QtCore.QObject):
    def __init__(self, gui, session):
        super().__init__(gui)
        self.gui, self.session = gui, session
        self.page = gui.tab_phased_transits
        self.graphics = gui.transit_phase_graphics
        self.graphics.setBackground('w')
        self.top = self.graphics.addPlot(row=0, col=0)
        self.bottom = self.graphics.addPlot(row=1, col=0)
        self.top.hideAxis('bottom')
        self.bottom.setXLink(self.top)
        self.graphics.ci.layout.setRowStretchFactor(0, 3)
        self.graphics.ci.layout.setRowStretchFactor(1, 1)
        self.legend = self.top.addLegend()
        self.planet, self.dataset = gui.transit_phase_planet, gui.transit_phase_dataset
        self.cross, self.slider = gui.transit_phase_cross, gui.transit_phase_slider
        from transit_phase_options_ui import Ui_TransitPhaseOptions
        self.dialog = QtWidgets.QDialog(gui)
        ui = Ui_TransitPhaseOptions()
        ui.setupUi(self.dialog)
        self.info, self.norm, self.gp = ui.info, ui.normalize, ui.gp
        self.points, self.dense = ui.points, ui.dense
        self.show_legend, self.auto, self.window = ui.show_legend, ui.auto, ui.window
        for widget in (self.norm, self.gp, self.points, self.dense, self.show_legend, self.auto):
            widget.toggled.connect(self.draw)
        self.window.valueChanged.connect(self.draw)
        ui.reset.clicked.connect(self.reset)
        gui.transit_phase_options_button.clicked.connect(self.show_options)
        self.planet.currentIndexChanged.connect(self.draw)
        self.dataset.currentIndexChanged.connect(self.draw)
        self.slider.valueChanged.connect(self.draw)
        self.cross.toggled.connect(self.draw)
        self._window_key = None
        self.refresh()

    def show_options(self):
        self.dialog.show()
        self.dialog.raise_()

    def reset(self):
        widgets = (self.slider, self.norm, self.gp, self.points, self.dense, self.show_legend, self.auto, self.window)
        blockers = [QtCore.QSignalBlocker(w) for w in widgets]
        self.slider.setValue(0)
        self.norm.setChecked(False)
        self.gp.setChecked(True)
        self.points.setChecked(False)
        self.dense.setChecked(True)
        self.show_legend.setChecked(False)
        self.auto.setChecked(False)
        self.window.setValue(0)
        self._window_key = None
        del blockers
        self.draw()

    def refresh(self):
        fit = self.session()
        old_slot, old_dataset = self.planet.currentData(), self.dataset.currentData()
        blockers = [QtCore.QSignalBlocker(w) for w in (self.planet, self.dataset)]
        self.planet.clear()
        self.dataset.clear()
        self.dataset.addItem('All light curves', -1)
        if fit.type_fit.get('Transit', False):
            for slot in range(len(fit.use_planet)):
                if fit.use_planet[slot]:
                    self.planet.addItem(f'pl {slot+1}', slot)
            try:
                sep = fit.transit_results[1]
                for j in range(len(fit.tra_data_sets)):
                    data = fit.tra_data_sets[j]
                    if len(data) and len(sep[0][j]):
                        self.dataset.addItem(str(data[-1]), j)
            except (TypeError, IndexError):
                pass
        for combo, previous in ((self.planet, old_slot), (self.dataset, old_dataset)):
            index = combo.findData(previous)
            combo.setCurrentIndex(index if index >= 0 else (0 if combo.count() else -1))
        self.planet.setEnabled(self.planet.count() > 0)
        self.dataset.setEnabled(self.dataset.count() > 1)
        del blockers
        self.draw()

    def draw(self, *args):
        fit = self.session()
        self.top.clear()
        self.bottom.clear()
        self.legend.clear()
        self.legend.setVisible(self.show_legend.isChecked())
        style = self.gui.plot_font_style()
        for panel in (self.top, self.bottom):
            panel.setLogMode(False, False)
            for side in ('left', 'bottom', 'top', 'right'):
                axis = panel.getAxis(side)
                axis.setTickFont(self.gui.plot_font)
                axis.setPen('k')
                axis.setTextPen('k')
                axis.enableAutoSIPrefix(False)
            panel.getAxis('left').setWidth(round(60 * self.gui.plot_font.pointSizeF()/11))
            panel.showAxis('right')
            panel.showAxis('top')
        self.top.setLabel('left', 'Relative flux', **style)
        self.bottom.setLabel('left', 'o-c', **style)
        self.bottom.setLabel('bottom', 'phase' if self.norm.isChecked() else 'phase [days]', **style)
        for panel in (self.top, self.bottom):
            panel.getAxis('top').setStyle(showValues=False)
            panel.getAxis('right').setStyle(showValues=False)
            for side in ('left', 'bottom', 'top', 'right'):
                panel.getAxis(side).enableAutoSIPrefix(False)
        slot = self.planet.currentData()
        if slot is None or not fit.type_fit.get('Transit', False) or self.dataset.count() < 2:
            self.info.setText('Initialize or fit a transit model to display phased light curves.')
            self.top.setTitle('No fitted transit model', **style)
            return
        period, epoch = float(fit.P[slot]), float(fit.t0[slot])
        if not np.isfinite(period) or period <= 0 or not np.isfinite(epoch):
            return
        self.top.setTitle('')
        self.info.setText(f'Planet {slot+1}: P = {period:.9g} days; transit epoch = {epoch:.9f}.\n'
                          'Uses cached fitted light curves and residuals. In multiple-planet systems, '
                          'the cached model includes all transits, folded at the selected planet’s period.')
        shift = self.slider.value()/100 * period
        plotted = False
        selected = self.dataset.currentData()
        dense_model = None
        if self.dense.isChecked() and not fit.tra_doGP and selected == -1:
            try:
                rich = fit.transit_results[3]
                times, values = np.asarray(rich[0], dtype=float), np.asarray(rich[1], dtype=float)
                if times.ndim == 1 and len(times) and values.shape == times.shape:
                    dense_model = (transit_phase(times, period, epoch, shift, self.norm.isChecked()), values)
            except (TypeError, IndexError, ValueError):
                pass
        for index in range(1, self.dataset.count()):
            j = self.dataset.itemData(index)
            if selected != -1 and selected != j:
                continue
            try:
                t, flux, error, model, residual = transit_series(fit, j, self.gp.isChecked())
            except (TypeError, IndexError, ValueError):
                continue
            x = transit_phase(t, period, epoch, shift, self.norm.isChecked())
            color = pg.mkColor(fit.tra_colors[j])
            color.setAlpha(int(fit.pyqt_color_alpha_tra[j]))
            symbol, size = fit.pyqt_symbols_tra[j], fit.pyqt_symbols_size_tra[j]
            for panel, y in ((self.top, flux), (self.bottom, residual)):
                panel.plot(x, y, pen=None, symbol=symbol, symbolSize=size,
                           symbolPen=color, symbolBrush=color,
                           name=self.dataset.itemText(index) if panel is self.top else None)
                panel.addItem(pg.ErrorBarItem(x=x, y=y, top=error, bottom=error, beam=0, pen=color))
            if dense_model is None:
                self.plot_model(x, model, fit)
            plotted = True
        if dense_model is not None:
            self.plot_model(*dense_model, fit)
        self.bottom.addLine(y=0, pen=pg.mkPen('#ff9933', width=.8))
        if self.cross.isChecked():
            self.gui.cross_hair(self.top, log=False)
            self.gui.cross_hair(self.bottom, log=False)
        if self.auto.isChecked() and plotted:
            self.top.autoRange()
            self.bottom.autoRange()
        # An explicit phase-window change controls X independently of auto range.
        # Other redraws preserve manual zoom when auto range is off.
        window_key = (id(fit), slot, period, self.norm.isChecked(), self.window.value())
        if plotted and (self.auto.isChecked() or window_key != self._window_key):
            width = min(self.window.value() or period/2, period/2)
            if self.norm.isChecked():
                width /= period
            self.top.setXRange(-width, width, padding=.002)
            self._window_key = window_key

    def plot_model(self, x, model, fit):
        order = np.argsort(x)
        self.top.plot(x[order], model[order],
                      pen=None if self.points.isChecked() else pg.mkPen(fit.tra_colors[-1], width=self.gui.tra_model_width.value()+1),
                      symbol='o' if self.points.isChecked() else None, symbolSize=2,
                      symbolPen=None, symbolBrush=fit.tra_colors[-1]).setZValue(self.gui.tra_model_z.value())
