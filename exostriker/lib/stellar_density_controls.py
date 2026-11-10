"""Connect the existing Stellar density bounds/priors pages to session fields."""
from PyQt6 import QtCore
import numpy as np
from stellar_density import ensure_density_settings


class StellarDensityControls:
    def __init__(self, gui, session):
        self.gui, self.session = gui, session
        self.groups = (
            ('stellar_rho_bounds', ('St_rho_min_1', 'St_rho_max_1'), None),
            ('stellar_rho_norm_pr', ('St_rho_mean_1', 'St_rho_sigma_1'), 'use_St_rho_norm_pr_1'),
            ('stellar_rho_jeff_pr', ('St_rho_jeff_alpha_1', 'St_rho_jeff_beta_1'), 'use_St_rho_jeff_pr_1'),
        )
        for field, names, flag in self.groups:
            for name in names:
                widget = getattr(gui, name)
                widget.valueChanged.connect(lambda value, f=field: self.read(f))
            if flag:
                getattr(gui, flag).toggled.connect(lambda value, f=field: self.toggle(f, value))
        gui.Stellar_rho_input.valueChanged.connect(self.value_changed)
        for slot in range(9):
            gui.param_gui_tr[slot*3+2].valueChanged.connect(lambda value, i=slot: self.planet_value_changed(i, value))
        self.sync()

    def read(self, field):
        for key, names, flag in self.groups:
            if key == field:
                values = [getattr(self.gui, name).value() for name in names]
                if flag:
                    values.append(getattr(self.gui, flag).isChecked())
                setattr(self.session(), key, values)

    def toggle(self, field, checked):
        self.read(field)
        if checked:
            other = 'stellar_rho_jeff_pr' if field == 'stellar_rho_norm_pr' else 'stellar_rho_norm_pr'
            for key, names, flag in self.groups:
                if key == other:
                    getattr(self.gui, flag).setChecked(False)

    def value_changed(self, value):
        session = self.session()
        session.stellar_rho = value
        slots = np.flatnonzero(session.use_planet)
        if session.rho and len(slots):
            index = int(slots[0])
            session.pl_a[index] = value
            widget = self.gui.param_gui_tr[index*3+2]
            blocker = QtCore.QSignalBlocker(widget)
            widget.setValue(value)
            del blocker

    def planet_value_changed(self, slot, value):
        session = self.session()
        slots = np.flatnonzero(session.use_planet)
        if session.rho and len(slots) and slot == slots[0]:
            session.stellar_rho = value
            session.pl_a[slot] = value
            blocker = QtCore.QSignalBlocker(self.gui.Stellar_rho_input)
            self.gui.Stellar_rho_input.setValue(value)
            del blocker

    def sync(self):
        session = self.session()
        ensure_density_settings(session)
        self.gui.label_ap.setText(self.gui.label_ap.property(
            'stellarDensityText' if session.rho else 'scaledAxisText'))
        for field, names, flag in self.groups:
            for name, value in zip(names, getattr(session, field)):
                widget = getattr(self.gui, name)
                blocker = QtCore.QSignalBlocker(widget)
                widget.setValue(value)
                del blocker
            if flag:
                widget = getattr(self.gui, flag)
                blocker = QtCore.QSignalBlocker(widget)
                widget.setChecked(bool(getattr(session, field)[2]))
                del blocker
        slots = np.flatnonzero(session.use_planet)
        if session.rho and len(slots):
            session.stellar_rho = float(session.pl_a[int(slots[0])])
        blocker = QtCore.QSignalBlocker(self.gui.Stellar_rho_input)
        self.gui.Stellar_rho_input.setValue(session.stellar_rho)
        del blocker
        for row in range(9):
            for group in (self.gui.param_bounds_gui, self.gui.param_nr_priors_gui, self.gui.param_jeff_priors_gui):
                for widget in (group[row*10+9][:2] if group is self.gui.param_bounds_gui else group[row*10+9]):
                    widget.setEnabled(not session.rho)
                    widget.setToolTip('In density mode, use Limits and Priors > Stellar.' if session.rho else '')
