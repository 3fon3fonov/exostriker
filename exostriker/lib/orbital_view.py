"""Interactive orthographic orbital preview, independent of the fitted model."""
from types import SimpleNamespace
import numpy as np
from PyQt6 import QtCore, QtWidgets
import pyqtgraph as pg


def project_xyz(xyz, azimuth=0., elevation=0.):
    """Rotate the camera; zero angles preserve the historical XY projection."""
    a, b = np.radians([azimuth, elevation])
    rz = np.array([[np.cos(a), -np.sin(a), 0],
                   [np.sin(a), np.cos(a), 0], [0, 0, 1]])
    rx = np.array([[1, 0, 0], [0, np.cos(b), -np.sin(b)],
                   [0, np.sin(b), np.cos(b)]])
    return (rx @ rz @ np.asarray(xyz)[:3]).reshape(3, -1)


class OrbitalViewController(QtCore.QObject):
    def __init__(self, parent, panel, session, orbit):
        super().__init__(parent)
        self.panel, self.session, self.orbit = panel, session, orbit
        self.overrides = {}
        self.source = None
        self.dialog = QtWidgets.QDialog(parent)
        self.dialog.setWindowTitle('Orbital view options')
        layout = QtWidgets.QVBoxLayout(self.dialog)
        layout.addWidget(QtWidgets.QLabel(
            'Preview only: orbital edits do not change the fit.\n'
            'Rotate with the camera controls or Shift + left drag; scroll to zoom.'))
        form = QtWidgets.QFormLayout()
        layout.addLayout(form)
        self.angles = []
        for label in ('Camera azimuth [deg]', 'Camera elevation [deg]'):
            box = QtWidgets.QDoubleSpinBox()
            box.setRange(-180, 180)
            box.setSingleStep(5)
            box.valueChanged.connect(self.draw)
            form.addRow(label, box)
            self.angles.append(box)
        self.axes = QtWidgets.QCheckBox('Show XYZ directions')
        self.axes.setChecked(False)
        self.peri = QtWidgets.QCheckBox('Show periastron directions')
        self.peri.setChecked(True)
        for box in (self.axes, self.peri):
            layout.addWidget(box)
            box.toggled.connect(self.draw)
        self.status = QtWidgets.QLabel()
        self.status.setWordWrap(True)
        layout.addWidget(self.status)
        self.table = QtWidgets.QTableWidget(0, 9)
        self.table.setHorizontalHeaderLabels(['Planet', 'a [au]', 'e', 'i [deg]',
                                              'Ω [deg]', 'ω [deg]', 'M₀ [deg]',
                                              'P [days]', 'K [m/s]'])
        self.table.itemChanged.connect(self.edit)
        layout.addWidget(self.table)
        reset = QtWidgets.QPushButton('Reset to fitted orbits and view')
        reset.clicked.connect(self.reset)
        layout.addWidget(reset)
        self.dialog.resize(930, 430)
        menu = panel.getPlotItem().getViewBox().menu
        menu.addAction('Orbital view options…', self.show)
        menu.addAction('Reset orbital view to fit', self.reset)
        button = QtWidgets.QPushButton('Orbital view options…', parent)
        reset_button = QtWidgets.QPushButton('Reset orbital view to fit', parent)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(button)
        row.addWidget(reset_button)
        parent.gridLayout_32.addLayout(row, 1, 0)
        button.clicked.connect(self.show)
        reset_button.clicked.connect(self.reset)
        self.drag = None
        panel.viewport().installEventFilter(self)

    def eventFilter(self, obj, event):
        kind = event.type()
        if kind == QtCore.QEvent.Type.MouseButtonPress and event.button() == QtCore.Qt.MouseButton.LeftButton and event.modifiers() & QtCore.Qt.KeyboardModifier.ShiftModifier:
            self.drag = event.position()
            return True
        if kind == QtCore.QEvent.Type.MouseMove and self.drag is not None:
            delta = event.position() - self.drag
            self.drag = event.position()
            for box, value in zip(self.angles, (delta.x(), delta.y())):
                box.setValue((box.value() + value * .5 + 180) % 360 - 180)
            return True
        if kind == QtCore.QEvent.Type.MouseButtonRelease and self.drag is not None:
            self.drag = None
            return True
        return super().eventFilter(obj, event)

    def show(self):
        self.populate()
        self.dialog.show()
        self.dialog.raise_()

    def populate(self):
        fit = self.session()
        self.slots = [i for i in range(len(fit.use_planet))
                      if fit.use_planet[i] and i < len(fit.fit_results.a)]
        blocker = QtCore.QSignalBlocker(self.table)
        self.table.setRowCount(len(self.slots))
        for row, slot in enumerate(self.slots):
            values = [slot + 1, fit.fit_results.a[slot]] + [getattr(fit, name)[slot] for name in ('e', 'i', 'Node', 'w', 'M0', 'P', 'K')]
            for col, value in enumerate(values):
                value = self.overrides.get((slot, col), value)
                item = QtWidgets.QTableWidgetItem(f'{value:.9g}')
                if col in (0, 7, 8):
                    item.setFlags(item.flags() & ~QtCore.Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(row, col, item)
        del blocker

    def edit(self, item):
        try:
            value = float(item.text())
            if not np.isfinite(value) or (item.column() == 1 and value <= 0) or (item.column() == 2 and not 0 <= value < 1):
                raise ValueError
        except ValueError:
            self.populate()
            return
        self.overrides[self.slots[item.row()], item.column()] = value
        self.draw()

    def reset(self):
        self.overrides.clear()
        for box in self.angles:
            blocker = QtCore.QSignalBlocker(box)
            box.setValue(0)
            del blocker
        self.axes.setChecked(False)
        self.peri.setChecked(True)
        self.populate()
        self.draw()
        self.panel.autoRange()

    def draw(self, *args):
        fit = self.session()
        # New sessions and refits invalidate previews, including same-object refits.
        signature = (id(fit), tuple(fit.fit_results.a),
                     *(tuple(getattr(fit, name)[i] for i in range(len(fit.use_planet))) for name in ('e', 'i', 'Node', 'w', 'M0')))
        if signature != self.source:
            self.source = signature
            self.overrides.clear()
            if self.dialog.isVisible():
                self.populate()
        panel = self.panel
        panel.clear()
        panel.setLabel('bottom', 'x [au]')
        panel.setLabel('left', 'y [au]')
        self.status.clear()
        if not len(fit.fit_results.a):
            self.status.setText('No fitted orbital semi-major axes available.')
            return
        preview = SimpleNamespace(params=fit.params,
            fit_results=SimpleNamespace(a=np.array(fit.fit_results.a, dtype=float, copy=True)))
        for name in ('e', 'i', 'Node', 'w', 'M0'):
            setattr(preview, name, {i: float(getattr(fit, name)[i]) for i in range(len(fit.use_planet))})
        for (slot, col), value in self.overrides.items():
            target = preview.fit_results.a if col == 1 else getattr(preview, {2:'e', 3:'i', 4:'Node', 5:'w', 6:'M0'}[col])
            target[slot] = value
        angles = [box.value() for box in self.angles]
        extent = 0.
        skipped = []
        drawn = 0
        stellar_mass = float(getattr(fit.params, 'stellar_mass', np.nan))
        for slot in range(len(fit.use_planet)):
            if not fit.use_planet[slot]:
                continue
            # Transit-only and unevaluated sessions can contain zero RV axes.
            # The Cartesian elliptic-orbit converter divides by a and 1-e*cos(E).
            values = [getattr(preview, name)[slot] for name in ('e', 'i', 'Node', 'w', 'M0')]
            if (slot >= len(preview.fit_results.a)
                    or not np.isfinite(preview.fit_results.a[slot])
                    or preview.fit_results.a[slot] <= 0
                    or not np.isfinite(stellar_mass) or stellar_mass <= 0
                    or not np.all(np.isfinite(values)) or not 0 <= preview.e[slot] < 1):
                skipped.append(str(slot + 1))
                continue
            xyz, planet, peri, _ = self.orbit(preview, slot)
            drawn += 1
            points = project_xyz(xyz, *angles)
            extent = max(extent, float(np.max(np.linalg.norm(xyz[:3], axis=0))))
            color = fit.colors[slot]
            panel.plot(*points[:2], pen=pg.mkPen(color, width=1.5))
            if self.peri.isChecked():
                p = project_xyz(np.asarray(peri).reshape(6, -1), *angles)
                panel.plot([0, p[0, 0]], [0, p[1, 0]], pen=pg.mkPen(color, style=QtCore.Qt.PenStyle.DashLine))
            p = project_xyz(planet, *angles)
            panel.plot(p[0], p[1], pen=None, symbol='o', symbolBrush=color, symbolSize=7)
        message = ('Cannot draw planet(s) ' + ', '.join(skipped)
                   + ': require a positive semi-major axis, stellar mass, and valid bound orbital elements.') if skipped else ''
        self.status.setText(message)
        if not drawn:
            text = pg.TextItem('No valid bound orbits available', color='#666666', anchor=(.5, .5))
            panel.addItem(text)
        panel.plot([0], [0], pen=None, symbol='o', symbolBrush='r', symbolSize=9)
        if self.axes.isChecked():
            ends = project_xyz(np.eye(3) * max(extent * .35, .01), *angles)
            for j, (name, color) in enumerate(zip('XYZ', ('#e05050', '#40a060', '#5080e0'))):
                x, y = ends[:2, j]
                panel.plot([0, x], [0, y], pen=pg.mkPen(color, width=2))
                text = pg.TextItem('+' + name + (' ⊙' if np.hypot(x, y) < 1.e-10 else ''), color=color, anchor=(0, 1))
                text.setPos(x, y)
                panel.addItem(text)
