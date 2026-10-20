"""Scrollable, linked pyqtgraph corner grid, independent of the main GUI."""
from sample_plot_labels import plot_label

import numpy as np
import pyqtgraph as pg
from PyQt6 import QtCore, QtGui
from scipy.ndimage import gaussian_filter1d


class CornerGrid:
    def __init__(self, widget):
        self.widget = widget
        self.panels = {}
        self.connections = []
        self._linking = False

    def clear(self):
        for signal, callback in self.connections:
            signal.disconnect(callback)
        self.connections.clear()
        self.widget.clear()
        self.panels.clear()

    def connect(self, signal, callback):
        signal.connect(callback)
        self.connections.append((signal, callback))

    def sync_axes(self, target, axis, interval):
        if self._linking:
            return
        self._linking = True
        try:
            target.enableAutoRange(axis=axis, enable=False)
            target.setRange(**{'xRange' if axis == 0 else 'yRange': interval}, padding=0)
        finally:
            self._linking = False

    def render(self, samples, labels, selected, options, pair_renderer, quantiles, font_size=9):
        font_size = options.get('font_size', font_size)
        self.clear()
        n = len(selected)
        self.widget.setMinimumSize(max(280, n*190), max(240, n*180))
        if not n:
            self.widget.addLabel('Choose parameters in Options → Corner parameters.')
            return
        subset = samples[:, selected]
        names = [labels[i] for i in selected]
        font = QtGui.QFont()
        font.setPointSize(font_size)
        for row in range(n):
            for column in range(row+1):
                panel = self.widget.addPlot(row=row, col=column)
                self.panels[row, column] = panel
                panel.setMenuEnabled(True)
                panel.hideButtons()
                for axis in ('left', 'bottom'):
                    panel.getAxis(axis).setStyle(tickFont=font, autoExpandTextSpace=True)
                    panel.getAxis(axis).enableAutoSIPrefix(False)
                if row == column:
                    self.histogram(panel, subset[:, column], options, quantiles)
                    panel.setTitle(plot_label(names[column]), size='%dpt' % font_size)
                    if row == 0:
                        panel.setLabel('left', 'Count', **{'font-size': '%dpt' % font_size})
                else:
                    pair_renderer(panel, subset, names, column, row, options, font_size)
                # Interior plots keep ticks but omit repeated text and labels.
                panel.setLabel('bottom', plot_label(names[column]) if row == n-1 else '',
                               **{'font-size': '%dpt' % font_size})
                panel.getAxis('bottom').setStyle(showValues=row == n-1)
                if row != column:
                    panel.setLabel('left', plot_label(names[row]) if column == 0 else '',
                                   **{'font-size': '%dpt' % font_size})
                panel.getAxis('left').setStyle(showValues=column == 0)
                panel.getAxis('left').setWidth(max(65, font_size*5) if column == 0 else 8)
                panel.getAxis('bottom').setHeight(max(48, font_size*4) if row == n-1 else 10)
                self.widget.ci.layout.setColumnStretchFactor(column, 1)
                self.widget.ci.layout.setRowStretchFactor(row, 1)
        for parameter in range(n):
            diagonal = self.panels[parameter, parameter].getViewBox()
            finite = subset[:, parameter]
            finite = finite[np.isfinite(finite)]
            if len(finite):
                lower, upper = finite.min(), finite.max()
                margin = (upper-lower)*.04 if upper > lower else max(abs(lower)*.01, .5)
                diagonal.setXRange(lower-margin, upper+margin, padding=0)
            for row in range(parameter+1, n):
                self.panels[row, parameter].setXLink(self.panels[parameter, parameter])
            if parameter:
                row_anchor = self.panels[parameter, 0].getViewBox()
                row_anchor.setYRange(*diagonal.viewRange()[0], padding=0)
                for column in range(1, parameter):
                    self.panels[parameter, column].setYLink(self.panels[parameter, 0])
                # Link this parameter's histogram X axis to its joint-panel Y axes.
                self.connect(diagonal.sigXRangeChanged,
                             lambda view, interval, target=row_anchor: self.sync_axes(target, 1, interval))
                self.connect(row_anchor.sigYRangeChanged,
                             lambda view, interval, target=diagonal: self.sync_axes(target, 0, interval))
        # Histogram count axes stay independent from parameter axes.
        for index in range(n):
            self.panels[index, index].enableAutoRange(axis='y', enable=True)

    @staticmethod
    def histogram(panel, values, options, quantiles):
        values = values[np.isfinite(values)]
        if not len(values):
            return
        counts, edges = np.histogram(values, bins=options['bins'])
        counts = counts.astype(float)
        if options['smooth'] > 0:
            counts = gaussian_filter1d(counts, options['smooth'])
        panel.plot(edges, counts, stepMode='center',
                   pen=pg.mkPen(options['color'], width=1.5), fillLevel=0,
                   brush=pg.mkBrush(*QtGui.QColor(options['color']).getRgb()[:3], 35))
        if options['quantile_lines']:
            for value in np.quantile(values, quantiles):
                panel.addItem(pg.InfiniteLine(value, angle=90,
                    pen=pg.mkPen(options['color'], style=QtCore.Qt.PenStyle.DashLine)))
        if options['summary'] != 'None':
            center = (np.median if options['summary'] == 'Median' else np.mean)(values)
            panel.addItem(pg.InfiniteLine(center, angle=90,
                                         pen=pg.mkPen(options['color'], width=1.5)))
