"""Reusable options for cached N-body evolution plots; never runs an integration."""
from PyQt6 import QtCore, QtWidgets
import pyqtgraph as pg

DEFAULTS = dict(line=True, dots=False, line_width=1.1, marker_size=1.,
                line_alpha=255, marker_alpha=255, line_style=0, symbol=0,
                color=None, planet_colors=None, hold=False, auto_range=True, xlabel=None, ylabel=None)
LINE_STYLES = (QtCore.Qt.PenStyle.SolidLine, QtCore.Qt.PenStyle.DashLine,
               QtCore.Qt.PenStyle.DotLine, QtCore.Qt.PenStyle.DashDotLine)
SYMBOLS = ('o','s','t','d','+')


class NbodyPlotOptionsController(QtCore.QObject):
    """Manage existing es.ui controls and plot state; never construct/load UI."""
    def __init__(self, container, controls):
        super().__init__(container)
        self.container = container
        self.controls = controls
        for name, control in controls.items():
            setattr(self, name, control)
        self._ready = False
        self._color = None
        self._planet_colors = {}
        self._curve_colors = {}
        self._labels = dict(xlabel=None, ylabel=None)

    def configure(self, key, panel, session, redraw, font_style):
        self.key, self.panel = key, panel
        self.session, self.redraw, self.font_style = session, redraw, font_style
        self.load_session()
        self._ready = True

    def load_session(self):
        options = dict(DEFAULTS)
        if self.key in ('a','omega','period_ratio','delta_omega','resonant_angles'):
            options.update(line=False,dots=True)
        options.update(getattr(self.session(),'nbody_plot_options',{}).get(self.key,{}))
        blockers = [QtCore.QSignalBlocker(child) for child in self.controls.values()]
        for name in ('line','dots','hold','auto_range'):
            getattr(self,name).setChecked(options[name])
        for name in ('line_width','marker_size','line_alpha','marker_alpha'):
            getattr(self,name).setValue(options[name])
        for name in ('line_style','symbol'):
            getattr(self,name).setCurrentIndex(options[name])
        self._color = options['color']
        self._planet_colors = {str(k): v for k, v in (options.get('planet_colors') or {}).items()}
        self._labels = {name:options[name] for name in ('xlabel','ylabel')}
        del blockers

    def values(self):
        result = {name:getattr(self,name).isChecked() for name in ('line','dots','hold','auto_range')}
        result.update({name:getattr(self,name).value() for name in ('line_width','marker_size','line_alpha','marker_alpha')})
        result.update({name:getattr(self,name).currentIndex() for name in ('line_style','symbol')})
        result.update(self._labels)
        result['color'] = self._color
        result['planet_colors'] = dict(self._planet_colors)
        return result

    def changed(self, *args):
        if not self._ready:
            return
        session = self.session()
        if not hasattr(session,'nbody_plot_options'):
            session.nbody_plot_options = {}
        session.nbody_plot_options[self.key] = self.values()
        self.redraw()

    def edit_label(self, name, side):
        text, accepted = QtWidgets.QInputDialog.getText(
            self.container, 'X-axis label' if side=='bottom' else 'Y-axis label',
            'Label (HTML and Greek symbols accepted; blank restores default):',
            text=self._labels[name] or self.panel.getAxis(side).labelText)
        if accepted:
            self._labels[name] = text if text.strip() else None
            self.changed()

    def choose_color(self):
        color = QtWidgets.QColorDialog.getColor(pg.mkColor(self._color or '#0066ff'),self.container,'Plot color')
        if color.isValid():
            self._color = color.name()
            self._planet_colors.clear()
            self.changed()

    def reset_color(self):
        self._color = None
        self._planet_colors.clear()
        self.changed()

    def choose_planet_color(self):
        if self.key not in ('a', 'e', 'omega', 'inclination'):
            QtWidgets.QMessageBox.information(self.container, 'Planet colors',
                'This plot shows a system or planet-pair quantity. Use Color to change its curve.')
            return
        fit = self.session()
        slots = (list(range(fit.npl_arb)) if fit.pl_arb_test else
                 [i for i, active in enumerate(fit.use_planet) if active])
        if not slots:
            return
        labels = [(f'Planet {chr(98+slot)} (pl {slot+1})' if slot < 25 else f'Planet {slot+1}')
                  for slot in slots]
        labels.append('Reset all colors to defaults')
        selected, accepted = QtWidgets.QInputDialog.getItem(
            self.container, 'Planet colors', 'Choose a planet for this plot:', labels, 0, False)
        if not accepted:
            return
        if selected == labels[-1]:
            self.reset_color()
            return
        slot = slots[labels.index(selected)]
        default = self._curve_colors.get(slot, fit.colors[slot])
        current = self._planet_colors.get(str(slot), self._color or default)
        color = QtWidgets.QColorDialog.getColor(pg.mkColor(current), self.container,
                                                f'Planet {chr(98+slot)} color')
        if color.isValid():
            self._planet_colors[str(slot)] = color.name()
            self.changed()

    def draw(self, curves, xlabel, ylabel):
        """Curves: (x, y, default_color[, planet_slot]); data are already cached."""
        options = self.values()
        panel = self.panel
        if not options['hold']:
            panel.clear()
        panel.enableAutoRange(x=options['auto_range'],y=options['auto_range'])
        self._curve_colors = {curve[3]: curve[2] for curve in curves if len(curve) == 4}
        for curve in curves:
            x, y, color = curve[:3]
            slot = curve[3] if len(curve) == 4 else None
            color = self._planet_colors.get(str(slot), options['color'] or color)
            line_color = pg.mkColor(color)
            line_color.setAlpha(int(options['line_alpha']))
            marker_color = pg.mkColor(color)
            marker_color.setAlpha(int(options['marker_alpha']))
            pen = pg.mkPen(line_color,width=options['line_width'],style=LINE_STYLES[options['line_style']]) if options['line'] else None
            panel.plot(x,y,pen=pen,symbol=SYMBOLS[options['symbol']] if options['dots'] else None,
                       symbolSize=options['marker_size'],symbolPen=marker_color,symbolBrush=marker_color)
        panel.setLabel('bottom',options['xlabel'] if options['xlabel'] is not None else xlabel,units='',**self.font_style())
        panel.setLabel('left',options['ylabel'] if options['ylabel'] is not None else ylabel,units='',**self.font_style())
        if options['auto_range'] and curves:
            panel.autoRange()
