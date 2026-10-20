"""Sample correlation rendering and dialogs, independent of the fitting GUI."""
from html import escape
from math import erf, sqrt

import numpy as np
import pyqtgraph as pg
from PyQt6 import QtCore, QtGui, QtWidgets
from scipy.stats import spearmanr
from scipy.ndimage import gaussian_filter
import contourpy
from color_dialog import choose_color
from print_info_window import print_info
from sample_correlation_options_ui import Ui_SampleCorrelationOptions
from sample_corner_plot import CornerGrid
from sample_plot_labels import plot_label

DEFAULTS = dict(symbol='o', size=2., opacity=.25, max_points=50000,
                bins=80, font_size=9, summary='None', color='#0066ff', plot_points=True,
                plot_contours=True, fill_contours=False, smooth=1.,
                convention='Corner 2D sigma', level1=39.3469, level2=86.4665,
                level3=98.8891, quantile_lines=False, quantiles='0.16, 0.5, 0.84')

COMBO_FIELDS = ('symbol', 'summary', 'convention')
VALUE_FIELDS = ('size', 'opacity', 'max_points', 'bins', 'smooth',
                'level1', 'level2', 'level3', 'font_size')
CHECK_FIELDS = ('plot_points', 'plot_contours', 'fill_contours', 'quantile_lines')


def parse_quantiles(text):
    try:
        quantiles = np.asarray([float(value.strip()) for value in text.split(',')])
    except ValueError:
        raise ValueError('Enter comma-separated quantiles, for example 0.16, 0.5, 0.84.') from None
    if (not len(quantiles) or not np.isfinite(quantiles).all()
            or np.any(quantiles <= 0) or np.any(quantiles >= 1)
            or np.any(np.diff(quantiles) <= 0)):
        raise ValueError('Quantiles must be increasing values strictly between 0 and 1.')
    return quantiles


def sample_data(fit, nested=False):
    sampler = getattr(fit, 'ns_sampler' if nested else 'mcmc_sampler', None)
    samples = np.asarray(getattr(sampler, 'samples', []), dtype=float)
    if samples.ndim != 2 or not samples.size:
        return None
    labels = list(getattr(sampler, 'lbf', {}).get('OrigLabels', []))
    labels = [str(labels[i]) if i < len(labels) else 'Parameter %d' % (i + 1)
              for i in range(samples.shape[1])]
    return samples, labels


class OptionsDialog(QtWidgets.QDialog, Ui_SampleCorrelationOptions):
    def __init__(self, parent, options, labels=None, selected=None):
        super().__init__(parent)
        self.setupUi(self)
        self.parametersTab.setEnabled(labels is not None)
        for index, label in enumerate(labels or []):
            item = QtWidgets.QListWidgetItem(str(label), self.cornerParameters)
            item.setFlags(item.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(QtCore.Qt.CheckState.Checked if index in (selected or [])
                               else QtCore.Qt.CheckState.Unchecked)
        self._color = options['color']
        self.quantiles.setText(options['quantiles'])
        for name in COMBO_FIELDS:
            getattr(self, name).setCurrentText(options[name])
        for name in VALUE_FIELDS:
            getattr(self, name).setValue(options.get(name, DEFAULTS[name]))
        for name in CHECK_FIELDS:
            getattr(self, name).setChecked(options[name])
        self.convention.currentTextChanged.connect(self.apply_preset)
        for name in ('level1', 'level2', 'level3'):
            getattr(self, name).valueChanged.connect(self.custom_levels)
        self.color.clicked.connect(self.choose_color)
        self.update_color()
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)

    def apply_preset(self, preset):
        if preset == 'Custom percentages':
            return
        levels = (1 - np.exp(-.5*np.arange(1, 4)**2)) if preset == 'Corner 2D sigma' else np.array([erf(n/sqrt(2)) for n in (1, 2, 3)])
        for name, value in zip(('level1', 'level2', 'level3'), levels):
            control = getattr(self, name)
            with QtCore.QSignalBlocker(control):
                control.setValue(value*100)

    def custom_levels(self):
        self.convention.setCurrentText('Custom percentages')

    def selected_parameters(self):
        return [i for i in range(self.cornerParameters.count())
                if self.cornerParameters.item(i).checkState() == QtCore.Qt.CheckState.Checked]

    def accept(self):
        if self.parametersTab.isEnabled() and not 1 <= len(self.selected_parameters()) <= 12:
            QtWidgets.QMessageBox.warning(self, 'Corner parameters',
                                         'Select between 1 and 12 parameters for the corner plot.')
            return
        try:
            parse_quantiles(self.quantiles.text())
        except ValueError as error:
            QtWidgets.QMessageBox.warning(self, 'Marginal quantiles', str(error))
            return
        if not self.plot_points.isChecked() and not self.plot_contours.isChecked():
            QtWidgets.QMessageBox.warning(self, 'Plot options', 'Enable points or contours.')
            return
        if not self.level1.value() < self.level2.value() < self.level3.value():
            QtWidgets.QMessageBox.warning(self, 'Contour probabilities',
                                         'Enter three increasing percentages between 0 and 100.')
            return
        super().accept()

    def update_color(self):
        self.color.setText(self._color)
        self.color.setStyleSheet('border: 2px solid %s;' % self._color)

    def choose_color(self):
        color = choose_color(self, self._color, 'Sample correlation color')
        if color.isValid():
            self._color = color.name()
            self.update_color()

    def values(self):
        result = {name: getattr(self, name).currentText()
                  for name in COMBO_FIELDS}
        result.update({name: getattr(self, name).value()
                       for name in VALUE_FIELDS})
        result.update({name: getattr(self, name).isChecked() for name in CHECK_FIELDS})
        return dict(result, color=self._color, quantiles=self.quantiles.text())


def render(plot, samples, labels, xindex, yindex, options, font_size=9):
    font_size = options.get('font_size', font_size)
    plot.clear()
    font = QtGui.QFont()
    font.setPointSize(font_size)
    for axis in ('left', 'bottom'):
        plot.getAxis(axis).setStyle(tickFont=font)
    plot.setLabel('bottom', plot_label(labels[xindex]), **{'font-size': '%dpt' % font_size})
    plot.setLabel('left', plot_label(labels[yindex]), **{'font-size': '%dpt' % font_size})
    pair = samples[:, [xindex, yindex]]
    pair = pair[np.isfinite(pair).all(axis=1)]
    if not len(pair):
        return
    color = QtGui.QColor(options['color'])
    color.setAlphaF(options['opacity'])
    if options['plot_points']:
        maximum = options['max_points']
        shown = pair
        if maximum and len(pair) > maximum:
            shown = pair[np.random.default_rng(0).choice(len(pair), maximum, replace=False)]
        plot.plot(shown[:, 0], shown[:, 1], pen=None, symbol=options['symbol'],
                  symbolSize=options['size'], symbolPen=pg.mkPen(color),
                  symbolBrush=pg.mkBrush(color))
    if options['plot_contours'] and np.all(np.ptp(pair, axis=0) > 0) and len(pair) >= 4:
        draw_contours(plot, pair, options)
    if options['quantile_lines']:
        for axis, angle in enumerate((90, 0)):
            for value in np.quantile(pair[:, axis], parse_quantiles(options['quantiles'])):
                plot.addItem(pg.InfiniteLine(value, angle=angle,
                    pen=pg.mkPen(options['color'], style=QtCore.Qt.PenStyle.DashLine)))
    if options['summary'] != 'None':
        center = (np.median if options['summary'] == 'Median' else np.mean)(pair, axis=0)
        for angle, value in zip((90, 0), center):
            plot.addItem(pg.InfiniteLine(value, angle=angle,
                                        pen=pg.mkPen(options['color'], width=1.5)))
    plot.enableAutoRange()


def contour_histogram(pair, bins, smooth, probabilities):
    """Corner-style thresholds from sorted histogram mass, not density height percentiles."""
    hist, xe, ye = np.histogram2d(pair[:, 0], pair[:, 1], bins=bins)
    if smooth > 0:
        hist = gaussian_filter(hist, smooth)
    ordered = np.sort(hist.ravel())[::-1]
    cumulative = np.cumsum(ordered) / ordered.sum()
    thresholds = []
    for probability in probabilities:
        # Match corner.hist2d: last cell at or below the requested cumulative mass.
        below = ordered[cumulative <= probability]
        thresholds.append(below[-1] if len(below) else ordered[0])
    return hist, xe, ye, np.asarray(thresholds)


def draw_contours(plot, pair, options):
    probabilities = [options['level%d' % n]/100 for n in (1, 2, 3)]
    hist, xe, ye, thresholds = contour_histogram(pair, options['bins'], options['smooth'], probabilities)
    dx, dy = xe[1]-xe[0], ye[1]-ye[0]
    x = np.r_[xe[0]-dx/2, (xe[:-1]+xe[1:])/2, xe[-1]+dx/2]
    y = np.r_[ye[0]-dy/2, (ye[:-1]+ye[1:])/2, ye[-1]+dy/2]
    grid = contourpy.contour_generator(x=x, y=y, z=np.pad(hist.T, 1),
                                      line_type='Separate', fill_type='OuterOffset')
    # Draw outer regions first; tied histogram levels are rendered only once.
    for threshold in sorted(set(thresholds)):
        if options['fill_contours']:
            polygons, offsets = grid.filled(threshold, float(hist.max()+1))
            for polygon, rings in zip(polygons, offsets):
                path = QtGui.QPainterPath()
                path.setFillRule(QtCore.Qt.FillRule.OddEvenFill)
                for start, end in zip(rings[:-1], rings[1:]):
                    points = polygon[start:end]
                    path.moveTo(*points[0])
                    for point in points[1:]:
                        path.lineTo(*point)
                    path.closeSubpath()
                item = QtWidgets.QGraphicsPathItem(path)
                color = QtGui.QColor(options['color'])
                color.setAlphaF(.13)
                item.setBrush(pg.mkBrush(color))
                item.setPen(pg.mkPen(None))
                item.setZValue(2)
                plot.addItem(item)
        for points in grid.lines(threshold):
            curve = pg.PlotCurveItem(points[:, 0], points[:, 1],
                                     pen=pg.mkPen(options['color'], width=1.5))
            curve.setZValue(3)
            plot.addItem(curve)


def statistics_html(samples, labels, xindex, yindex, source, options=None):
    options = DEFAULTS if options is None else options
    rows = ['<h1>%s sample statistics</h1>' % escape(source),
            '<p>%d samples; %d parameters. All finite values are used independently '
            'of display sampling.<br>σ is the sample standard deviation (ddof=1). '
            'Angles are summarized linearly on their stored branch.</p>' % samples.shape]
    table = '<table border="1" cellspacing="0" cellpadding="5">'
    rows.append(table+'<tr><th>Selected pair</th><th>Finite pairs</th><th>Pearson r</th><th>Spearman ρ</th></tr>')
    pair = samples[:, [xindex, yindex]]
    pair = pair[np.isfinite(pair).all(axis=1)]
    pearson = spearman = 'undefined'
    if len(pair) > 1 and np.all(np.ptp(pair, axis=0) > 0):
        pearson = '%.6g' % np.corrcoef(pair.T)[0, 1]
        spearman = '%.6g' % spearmanr(pair[:, 0], pair[:, 1]).statistic
    rows.append('<tr><td>%s vs %s</td><td>%d</td><td>%s</td><td>%s</td></tr></table>' %
                (escape(labels[xindex]), escape(labels[yindex]), len(pair), pearson, spearman))
    rows.append('<h2>Parameter summaries</h2>'+table+'<tr>'+''.join(
        '<th>%s</th>' % text for text in ['Parameter', 'Finite / excluded', 'Mean', 'Median', 'σ', 'MAD', 'IQR', 'Minimum', 'Maximum'])+'</tr>')
    columns = []
    for column, label in enumerate(labels):
        values = samples[:, column]
        values = values[np.isfinite(values)]
        columns.append(values)
        cells = [escape(label), '%d / %d' % (len(values), len(samples)-len(values))]
        if len(values):
            cells += ['%.8g' % np.mean(values), '%.8g' % np.median(values),
                      '%.8g' % np.std(values, ddof=1) if len(values)>1 else 'undefined',
                      '%.8g' % np.median(np.abs(values-np.median(values))),
                      '%.8g' % np.subtract(*np.percentile(values, [75, 25])),
                      '%.8g' % values.min(), '%.8g' % values.max()]
        else:
            cells += ['undefined']*7
        rows.append('<tr>'+''.join('<td>%s</td>' % cell for cell in cells)+'</tr>')
    rows.append('</table><h2>Marginal equal-tail credible intervals</h2>'+table+
                '<tr><th>Parameter</th><th>1σ (68.2689%)</th><th>2σ (95.4500%)</th><th>3σ (99.7300%)</th></tr>')
    for label, values in zip(labels, columns):
        cells = [escape(label)]
        for n in (1, 2, 3):
            coverage = erf(n/sqrt(2))
            bounds = np.quantile(values, [(1-coverage)/2, (1+coverage)/2]) if len(values) else None
            cells.append('[%.8g, %.8g]' % tuple(bounds) if bounds is not None else 'undefined')
        rows.append('<tr>'+''.join('<td>%s</td>' % cell for cell in cells)+'</tr>')
    rows.append('</table><h2>Adopted 2D contour probabilities</h2>'+table+
                '<tr><th>Contour</th><th>Enclosed posterior probability</th></tr>')
    for n in (1, 2, 3):
        rows.append('<tr><td>%d</td><td>%.4f%%</td></tr>' % (n, options['level%d' % n]))
    rows.append('</table><p>Contour levels enclose probability in the joint 2D distribution, '
                'estimated from the smoothed histogram, as in corner.py.<br>They differ from '
                'the marginal equal-tail intervals above and from mean ± nσ. '
                '<br>These summaries do not measure convergence or effective sample size.</p>')
    return ''.join(rows)


class SampleCorrelationController:
    def __init__(self, parent, fit_provider, plot):
        self.parent, self.fit_provider, self.plot = parent, fit_provider, plot
        self.options = DEFAULTS.copy()
        self.info_window = None
        self.corner = CornerGrid(parent.sample_corner_grid)
        self.parameter_selections = {}
        parent.corner_mode.toggled.connect(self.render)
        parent.corr_plot_options.clicked.connect(self.show_options)
        parent.samp_corr_print_info.clicked.connect(self.show_statistics)
        parent.corr_plot_options.setEnabled(True)
        parent.samp_corr_print_info.setEnabled(True)

    def data(self):
        return sample_data(self.fit_provider(), self.parent.samp_corr_ns.isChecked())

    def indices(self, samples):
        x = self.parent.comboBox_samp_corr_1.currentIndex()
        y = self.parent.comboBox_samp_corr_2.currentIndex()
        return (x, y) if 0 <= x < samples.shape[1] and 0 <= y < samples.shape[1] else None

    def selection_key(self, labels):
        return self.parent.samp_corr_ns.isChecked(), tuple(labels)

    def selected_parameters(self, labels):
        key = self.selection_key(labels)
        return self.parameter_selections.setdefault(key, list(range(min(4, len(labels)))))

    def render(self, *args):
        corner_mode = self.parent.corner_mode.isChecked()
        self.parent.sample_plot_stack.setCurrentIndex(int(corner_mode))
        for combo in (self.parent.comboBox_samp_corr_1, self.parent.comboBox_samp_corr_2):
            combo.setEnabled(not corner_mode)
        data = self.data()
        if corner_mode:
            if data is None:
                self.corner.clear()
                self.parent.sample_corner_grid.setMinimumSize(280, 240)
                self.parent.sample_corner_grid.addLabel('Run or load MCMC / nested samples first.')
                return
            self.corner.render(*data, self.selected_parameters(data[1]), self.options,
                               render, parse_quantiles(self.options['quantiles']),
                               self.parent.plot_font.pointSize())
        else:
            self.corner.clear()
            if data is None or self.indices(data[0]) is None:
                self.plot.clear()
                return
            render(self.plot, *data, *self.indices(data[0]), self.options,
                   self.parent.plot_font.pointSize())

    def show_options(self):
        data = self.data()
        labels = data[1] if data is not None else None
        selected = self.selected_parameters(labels) if labels is not None else None
        dialog = OptionsDialog(self.parent, self.options, labels, selected)
        if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            self.options = dialog.values()
            if labels is not None:
                self.parameter_selections[self.selection_key(labels)] = dialog.selected_parameters()
            self.render()

    def show_statistics(self):
        data = self.data()
        if data is None or self.indices(data[0]) is None:
            QtWidgets.QMessageBox.information(self.parent, 'Sample statistics',
                                             'Run or load a sampler and select two parameters first.')
            return
        if self.info_window is None:
            self.info_window = print_info(self.parent)
            self.info_window.text.setLineWrapMode(QtWidgets.QTextEdit.LineWrapMode.NoWrap)
        self.info_window.setWindowTitle('Sample correlation statistics')
        source = 'Nested posterior (equal-weight resampled)' if self.parent.samp_corr_ns.isChecked() else 'MCMC'
        self.info_window.text.setHtml(statistics_html(*data, *self.indices(data[0]), source, self.options))
        self.info_window.resize(1200, 700)
        self.info_window.show()
        self.info_window.raise_()
