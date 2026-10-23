"""Shared appearance options for pyqtgraph axes (sizes in screen pixels)."""
import math

DEFAULTS = dict(axis_width=1.0, tick_width=1.0, tick_length=5,
                tick_direction=0, x_density=1.0, y_density=1.0,
                minor_levels=2, tick_opacity=-1)
RANGES = dict(axis_width=(0.1, 10.), tick_width=(0.1, 10.), tick_length=(0, 40),
              tick_direction=(0, 1), x_density=(0.2, 5.), y_density=(0.2, 5.),
              minor_levels=(0, 2), tick_opacity=(-1, 100))
INTEGER_FIELDS = {'tick_length', 'tick_direction', 'minor_levels', 'tick_opacity'}


def normalize_options(options):
    result = dict(DEFAULTS)
    for name, default in DEFAULTS.items():
        value = float(options.get(name, default))
        low, high = RANGES[name]
        if not math.isfinite(value) or not low <= value <= high:
            raise ValueError('Invalid plot axis option: ' + name)
        if name in INTEGER_FIELDS:
            if not value.is_integer():
                raise ValueError('Plot axis option must be an integer: ' + name)
            value = int(value)
        result[name] = value
    return result


def apply_axis_options(panel, options):
    """Apply pens/ticks without changing ranges, labels or hidden-axis values."""
    from pyqtgraph.Qt import QtGui
    options = normalize_options(options)
    length = options['tick_length']
    if options['tick_direction'] == 0:  # negative is away from the axis text
        length = -length
    alpha = None if options['tick_opacity'] == -1 else options['tick_opacity']/100.
    for side in ('left', 'bottom', 'top', 'right'):
        axis = panel.getAxis(side)
        pen = QtGui.QPen(axis.pen())
        pen.setWidthF(options['axis_width'])
        pen.setCosmetic(True)
        tick_pen = QtGui.QPen(axis.tickPen())
        tick_pen.setWidthF(options['tick_width'])
        tick_pen.setCosmetic(True)
        axis.setPen(pen)
        axis.setTickPen(tick_pen)
        axis.setStyle(tickLength=length, maxTickLevel=options['minor_levels'],
                      tickAlpha=alpha if length else 0.)
        axis.setTickDensity(options['x_density'] if side in ('bottom','top') else options['y_density'])
