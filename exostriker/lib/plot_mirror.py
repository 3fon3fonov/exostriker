"""Independent, live pyqtgraph display copies for the Extra Plots panel."""
import copy

import numpy as np
import pyqtgraph as pg
from PyQt6 import QtCore, QtGui, QtWidgets


def plot_item(plot):
    return plot.getPlotItem() if isinstance(plot, pg.PlotWidget) else plot


def plot_navigation(plot):
    """Find the source widget and its available enclosing tab pages."""
    widget = plot if isinstance(plot, QtWidgets.QWidget) else plot.getViewWidget()
    if widget is None:
        return None, []
    pages = []
    ancestor = widget.parentWidget()
    while ancestor is not None:
        if isinstance(ancestor, QtWidgets.QTabWidget):
            page = next((ancestor.widget(index) for index in range(ancestor.count())
                         if ancestor.widget(index) is widget
                         or ancestor.widget(index).isAncestorOf(widget)), None)
            # Removed tabs remain parented to their old stacked widget but
            # cannot be selected. Do not navigate to an unrelated outer page.
            if page is None:
                return None, []
            pages.append((ancestor, page))
        ancestor = ancestor.parentWidget()
    return widget, pages


def reveal_plot(plot):
    """Open every enclosing tab needed to show an original plot widget."""
    widget, pages = plot_navigation(plot)
    if widget is None:
        return False
    for tabs, page in pages:
        tabs.setCurrentWidget(page)
    widget.setFocus(QtCore.Qt.FocusReason.OtherFocusReason)
    return True


def display_copy(value):
    """Copy arrays/styles without attempting to pickle Qt graphics objects."""
    if isinstance(value, np.ndarray):
        return value.copy()
    if isinstance(value, (QtGui.QPen, QtGui.QBrush, QtGui.QColor, QtGui.QFont,
                          QtGui.QPainterPath, QtGui.QTransform)):
        return type(value)(value)
    if isinstance(value, dict):
        return {key: display_copy(item) for key, item in value.items()}
    if isinstance(value, list):
        return [display_copy(item) for item in value]
    if isinstance(value, tuple):
        return tuple(display_copy(item) for item in value)
    return value


def clone_item(source, copies):
    """Clone display data only; never reparent an item or copy its callbacks."""
    if isinstance(source, pg.FillBetweenItem):
        curves = [copies.get(curve) for curve in source.curves or ()]
        if len(curves) != 2 or any(curve is None for curve in curves):
            return None
        result = pg.FillBetweenItem(*curves, brush=QtGui.QBrush(source.brush()),
                                   pen=QtGui.QPen(source.pen()), fillRule=source.fillRule())
    elif isinstance(source, pg.PlotDataItem):
        options = display_copy(source.opts)
        result = pg.PlotDataItem(x=display_copy(source.xData), y=display_copy(source.yData), **options)
    elif isinstance(source, pg.PlotCurveItem):
        x, y = source.getData()
        result = pg.PlotCurveItem(x=display_copy(x), y=display_copy(y), **display_copy(source.opts))
    elif isinstance(source, pg.ErrorBarItem):
        result = pg.ErrorBarItem(**display_copy(source.opts))
    elif isinstance(source, pg.ScatterPlotItem):
        spots = [dict(pos=point.pos(), size=point.size(), symbol=display_copy(point.symbol()),
                      pen=point.pen(), brush=point.brush())
                 for point in source.points() if point.isVisible()]
        options = {key: display_copy(source.opts[key])
                   for key in ('pxMode', 'useCache', 'antialias', 'name')}
        result = pg.ScatterPlotItem(spots=spots, **options)
    elif isinstance(source, pg.InfiniteLine):
        result = pg.InfiniteLine(pos=source.pos(), angle=source.angle,
                                 pen=QtGui.QPen(source.pen), movable=False,
                                 span=source.span)
        if getattr(source, 'label', None) is not None:
            # Copy the currently displayed annotation, not an interactive link.
            label = source.label
            result.label = pg.InfLineLabel(result, text=label.format,
                                           position=label.orthoPos,
                                           color=label.color)
    elif isinstance(source, pg.TextItem):
        result = pg.TextItem(html=source.textItem.toHtml(), anchor=source.anchor,
                             border=display_copy(source.border), fill=display_copy(source.fill),
                             angle=source.angle, rotateAxis=source.rotateAxis)
        result.setFont(QtGui.QFont(source.textItem.font()))
        # QTextDocument HTML omits the QGraphicsTextItem's default text color.
        # GLS peak labels use that default, so HTML alone changes gray to black.
        result.setColor(QtGui.QColor(source.textItem.defaultTextColor()))
    elif isinstance(source, pg.ArrowItem):
        result = pg.ArrowItem(**display_copy(source.opts))
    elif isinstance(source, pg.LinearRegionItem):
        result = pg.LinearRegionItem(values=source.getRegion(), orientation=source.orientation,
                                     brush=QtGui.QBrush(source.brush),
                                     pen=QtGui.QPen(source.lines[0].pen), movable=False,
                                     span=source.span, swapMode=source.swapMode)
    elif isinstance(source, pg.ImageItem):
        result = pg.ImageItem(image=display_copy(source.image), axisOrder=source.axisOrder)
        result.setLevels(display_copy(source.getLevels()))
        if source.lut is not None:
            result.setLookupTable(display_copy(source.lut))
    else:
        return None
    if isinstance(result, pg.PlotDataItem):
        result.setPos(source.pos().x(), source.pos().y())
    else:
        result.setPos(source.pos())
    result.setTransform(QtGui.QTransform(source.transform()))
    result.setRotation(source.rotation())
    result.setScale(source.scale())
    result.setOpacity(source.opacity())
    result.setZValue(source.zValue())
    result.setVisible(source.isVisible())
    return result


class PlotMirror(QtCore.QObject):
    """Debounced one-way refresh; pan/zoom in the copy never edits its source."""
    statusChanged = QtCore.pyqtSignal(str)
    refreshed = QtCore.pyqtSignal()

    def __init__(self, target, parent=None):
        super().__init__(parent)
        self.target = target
        self.source = None
        self.bottom_axis = None
        self._scene = None
        self._last_range = None
        self.copies = {}
        self.timer = QtCore.QTimer(self)
        self.timer.setSingleShot(True)
        self.timer.setInterval(100)
        self.timer.timeout.connect(self.refresh)
        self.target.installEventFilter(self)

    def eventFilter(self, watched, event):
        if watched is self.target:
            if event.type() == QtCore.QEvent.Type.Show:
                self.schedule()
            elif event.type() == QtCore.QEvent.Type.Hide:
                self.timer.stop()
        return False

    def set_source(self, source, bottom_axis=None):
        self.timer.stop()
        if self._scene is not None:
            try:
                self._scene.changed.disconnect(self.schedule)
            except (RuntimeError, TypeError):
                pass
        self.source = plot_item(source)
        self.bottom_axis = bottom_axis
        self._last_range = None
        self._scene = self.source.scene() if self.source is not None else None
        self.target.clear()
        self.copies = {}
        if self._scene is not None:
            self._scene.changed.connect(self.schedule)
        self.schedule()

    def schedule(self, *args):
        if self.source is not None and self.target.isVisible() and not self.timer.isActive():
            self.timer.start()

    def refresh(self):
        self.timer.stop()
        if self.source is None:
            return
        source, target = self.source, self.target.getPlotItem()
        source_view = source.getViewBox()
        # Hidden tabs may not have painted since their data changed. Complete
        # their pending autorange so the mirror does not inherit a stale view.
        if source_view._autoRangeNeedsUpdate:
            source_view.updateAutoRange()
        target.clear()
        if target.legend is not None:
            target.legend.clear()
            target.legend.setVisible(False)
        # PlotItem applies its own log/FFT/display controls when adding curves.
        # Copy those controls before adding items so transformations occur once.
        target.setLogMode(source.ctrl.logXCheck.isChecked(), source.ctrl.logYCheck.isChecked())
        target.ctrl.fftCheck.setChecked(source.ctrl.fftCheck.isChecked())
        target.setDownsampling(*source.downsampleMode())
        target.setClipToView(source.clipToViewMode())
        target.ctrl.alphaGroup.setChecked(source.ctrl.alphaGroup.isChecked())
        target.ctrl.autoAlphaCheck.setChecked(source.ctrl.autoAlphaCheck.isChecked())
        target.ctrl.alphaSlider.setValue(source.ctrl.alphaSlider.value())
        target.ctrl.subtractMeanCheck.setChecked(source.ctrl.subtractMeanCheck.isChecked())
        for name in ('left', 'bottom', 'right', 'top'):
            original = source.getAxis(name)
            if name == 'bottom' and not original.labelText and self.bottom_axis is not None:
                original = self.bottom_axis
            axis = target.getAxis(name)
            axis.setLabel(original.labelText, units=original.labelUnits,
                          **display_copy(original.labelStyle))
            style = display_copy(original.style)
            style['tickTextOffset'] = int(style['tickTextOffset'][0 if name in ('left', 'right') else 1])
            axis.setStyle(**style)
            axis.setGrid(original.grid)
            axis.setPen(QtGui.QPen(original.pen()))
            axis.setTextPen(QtGui.QPen(original.textPen()))
            # Disabling SI prefixes in this bundled pg version still computes
            # a prefix once from the old range. Copy the actual axis scale so
            # switching from BJD to small astrometric offsets keeps its units.
            axis.autoSIPrefix = original.autoSIPrefix
            axis.autoSIPrefixScale = original.autoSIPrefixScale
            axis.labelUnitPrefix = original.labelUnitPrefix
            axis.setScale(original.scale)
            axis._updateLabel()
            axis.setTicks(display_copy(original._tickLevels))
            if name in ('left', 'right'):
                axis.setWidth(original.fixedWidth)
            else:
                axis.setHeight(original.fixedHeight)
            target.showAxis(name, name in ('left', 'bottom') or original.isVisible())
        title = source.titleLabel.text if source.titleLabel.isVisible() else None
        target.setTitle(title or None, **display_copy(source.titleLabel.opts))
        view = source.getViewBox()
        state = view.getState()
        target.getViewBox().invertX(state['xInverted'])
        target.getViewBox().invertY(state['yInverted'])
        aspect = state['aspectLocked']
        target.setAspectLocked(aspect is not False, ratio=aspect if aspect is not False else 1)
        target.disableAutoRange()
        self.copies = {}
        skipped = []
        items = [item for item in source.items
                 if not getattr(item, '_exostriker_cursor_overlay', False)]
        # Fill bands need the cloned curves to exist first.
        for item in sorted(items, key=lambda obj: isinstance(obj, pg.FillBetweenItem)):
            try:
                clone = clone_item(item, self.copies)
            except (TypeError, ValueError, AttributeError, RuntimeError):
                clone = None
            if clone is None:
                skipped.append(type(item).__name__)
                continue
            self.copies[item] = clone
            target.addItem(clone)
        # PlotItem's trace-management controls may show curves while adding
        # later items. Restore the source visibility after all additions.
        for item, clone in self.copies.items():
            clone.setVisible(item.isVisible())
        legend = source.legend
        if legend is not None and legend.isVisible():
            options = display_copy(legend.opts)
            options.pop('offset', None)
            destination = target.addLegend(offset=legend.offset or (10, 10),
                                           colCount=legend.columnCount, **options)
            destination.clear()
            destination.setPen(QtGui.QPen(legend.pen()))
            destination.setBrush(QtGui.QBrush(legend.brush()))
            destination.setColumnCount(legend.columnCount)
            destination.setOffset(legend.offset or (10, 10))
            destination.opts.update(options)
            destination.setVisible(True)
            for sample, label in legend.items:
                clone = self.copies.get(sample.item)
                if clone is not None:
                    destination.addItem(clone, label.text)
        ranges = view.viewRange()
        if ranges != self._last_range:
            target.setRange(xRange=ranges[0], yRange=ranges[1], padding=0)
            self._last_range = copy.deepcopy(ranges)
        if skipped:
            self.statusChanged.emit('Some annotations could not be copied: '
                                    + ', '.join(sorted(set(skipped))))
        elif not items:
            self.statusChanged.emit('This source plot has no plotted data yet.')
        else:
            self.statusChanged.emit('')
        self.refreshed.emit()


class PlotCrosshair(QtCore.QObject):
    """One local cursor per display, including coordinates on logarithmic axes."""
    def __init__(self, target, parent=None):
        super().__init__(parent)
        self.target = target
        self.enabled = False
        self.has_position = False
        self.vertical = pg.InfiniteLine(angle=90, movable=False)
        self.horizontal = pg.InfiniteLine(angle=0, movable=False)
        # Match Exo_striker.cross_hair: default pg line/text colors and no fill.
        self.label = pg.TextItem(anchor=(0, 1))
        self.items = (self.vertical, self.horizontal, self.label)
        for item in self.items:
            item._exostriker_cursor_overlay = True
            item.setZValue(1e6)
            item.hide()
        self.proxy = pg.SignalProxy(target.scene().sigMouseMoved, rateLimit=60, slot=self.mouse_moved)
        self.proxy.setParent(self)
        target.viewport().installEventFilter(self)

    def eventFilter(self, watched, event):
        if event.type() in (QtCore.QEvent.Type.Leave, QtCore.QEvent.Type.Hide):
            self.reset()
        return False

    def reset(self):
        self.has_position = False
        for item in self.items:
            item.hide()

    def set_enabled(self, enabled):
        self.enabled = bool(enabled)
        self.reset()
        if self.enabled:
            self.attach()
        else:
            for item in self.items:
                if item in self.target.getPlotItem().items:
                    self.target.removeItem(item)

    def attach(self):
        # Source refreshes clear the destination. Reuse these three objects and
        # the existing mouse connection rather than accumulating new cursors.
        if not self.enabled:
            return
        for item in self.items:
            if item not in self.target.getPlotItem().items:
                self.target.addItem(item, ignoreBounds=True)
            item.setVisible(self.has_position)

    def mouse_moved(self, event):
        if not self.enabled:
            return
        position = event[0]
        plot = self.target.getPlotItem()
        view = plot.getViewBox()
        if not view.sceneBoundingRect().contains(position):
            self.reset()
            return
        point = view.mapSceneToView(position)
        x, y = point.x(), point.y()
        with np.errstate(over='ignore'):
            display_x = np.power(10., x) if plot.ctrl.logXCheck.isChecked() else x
            display_y = np.power(10., y) if plot.ctrl.logYCheck.isChecked() else y
        self.vertical.setPos(x)
        self.horizontal.setPos(y)
        self.label.setText(f'x={display_x:.12g},  y={display_y:.12g}')
        ranges = view.viewRange()
        self.label.setAnchor((0 if x < np.mean(ranges[0]) else 1,
                              1 if y < np.mean(ranges[1]) else 0))
        self.label.setPos(x, y)
        self.has_position = True
        self.attach()


class ExtraPlotMirror(QtCore.QObject):
    """One source list and one crosshair control for every extra plot."""
    def __init__(self, window, sources):
        super().__init__(window)
        self.window = window
        self.sources = sources
        self.selector = getattr(window, 'comboBox_plot_mirror', None)
        row = getattr(window, 'horizontalLayout_plot_mirror', None)
        if self.selector is None:
            self.selector = QtWidgets.QComboBox(window.tab_186)
            self.selector.setObjectName('comboBox_plot_mirror')
            row = QtWidgets.QHBoxLayout()
            row.addWidget(QtWidgets.QLabel('Plot source', window.tab_186))
            row.addWidget(self.selector, 1)
            window.gridLayout_410.addLayout(row, 2, 0)
            window.comboBox_plot_mirror = self.selector
            window.horizontalLayout_plot_mirror = row
        self.crosshair_checkbox = window.extra_plot_cross_hair_gls
        self.crosshair_checkbox.setText('Cross hair')
        old_controls = getattr(window, 'gridLayout_409', None)
        if old_controls is not None:
            old_controls.removeWidget(self.crosshair_checkbox)
            window.gridLayout_410.removeItem(old_controls)
        row.addWidget(self.crosshair_checkbox)
        self.original_button = getattr(window, 'button_original_plot', None)
        if self.original_button is None:
            self.original_button = QtWidgets.QPushButton('to the Org. Plot', window.tab_186)
            self.original_button.setObjectName('button_original_plot')
            window.button_original_plot = self.original_button
        else:
            row.removeWidget(self.original_button)
        self.original_button.setToolTip('Show the selected source plot in the left panel')
        row.addWidget(self.original_button)
        self.original_button.clicked.connect(self.show_original)
        window.comboBox_extra_plot_gls.hide()
        window.graphicsView_extra_plot_gls.hide()
        self.status = QtWidgets.QLabel(window.tab_186)
        self.status.setWordWrap(True)
        self.status.hide()
        window.gridLayout_410.addWidget(self.status, 3, 0)
        from pg_hack import CustomAxisItem
        self.target = pg.PlotWidget(parent=window.tab_186,
                                    axisItems={'bottom': CustomAxisItem('bottom')})
        self.target.setObjectName('graphicsView_plot_mirror')
        window.gridLayout_408.addWidget(self.target, 0, 1)
        self.mirror = PlotMirror(self.target, self)
        self.mirror.statusChanged.connect(self.set_status)
        self.crosshair = PlotCrosshair(self.target, self)
        self.mirror.refreshed.connect(self.crosshair.attach)
        self.crosshair_checkbox.toggled.connect(self.crosshair.set_enabled)
        self.crosshair.set_enabled(self.crosshair_checkbox.isChecked())
        for key, (label, getter, bottom_getter) in sources.items():
            self.selector.addItem(label, key)
        self.selector.setMinimumContentsLength(24)
        self.selector.setSizeAdjustPolicy(QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        default = max(0, self.selector.findData('p7'))
        self.selector.setCurrentIndex(default)
        self.selector.currentIndexChanged.connect(self.select)
        self.select(default)

    def set_status(self, text):
        self.status.setText(text)
        self.status.setVisible(bool(text))

    def show_original(self):
        key = self.selector.currentData()
        label, getter, bottom_getter = self.sources[key]
        return reveal_plot(getter())

    def select(self, index):
        key = self.selector.itemData(index)
        self.crosshair.reset()
        self.set_status('')
        label, getter, bottom_getter = self.sources[key]
        original = getter()
        widget, pages = plot_navigation(original)
        self.original_button.setEnabled(widget is not None)
        self.original_button.setToolTip('Show the selected source plot in the left panel'
                                       if widget is not None else
                                       'The original plot page is not available in the current layout')
        fallback = plot_item(bottom_getter()).getAxis('bottom') if bottom_getter else None
        self.mirror.set_source(original, fallback)
        self.crosshair.attach()
