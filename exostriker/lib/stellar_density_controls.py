"""Connect the existing stellar editors, fit flags and prior pages to sessions."""
from PyQt6 import QtCore
import numpy as np
from stellar_parameters import ensure_stellar_settings, sync_stellar, radius_mode


class StellarDensityControls:
    def __init__(self, gui, session):
        self.gui, self.session = gui, session
        self.groups = (
            ('stellar_mass_bounds', ('St_mass_min_1', 'St_mass_max_1'), None),
            ('stellar_mass_norm_pr', ('St_mass_mean_1', 'St_mass_sigma_1'), 'use_St_mass_norm_pr_1'),
            ('stellar_mass_jeff_pr', ('St_mass_jeff_alpha_1', 'St_mass_jeff_beta_1'), 'use_St_mass_jeff_pr_1'),
            ('stellar_radius_bounds', ('St_radius_min_1', 'St_radius_max_1'), None),
            ('stellar_radius_norm_pr', ('St_radius_mean_1', 'St_radius_sigma_1'), 'use_St_radius_norm_pr_1'),
            ('stellar_radius_jeff_pr', ('St_radius_jeff_alpha_1', 'St_radius_jeff_beta_1'), 'use_St_radius_jeff_pr_1'),
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
        gui.use_st_mass.toggled.connect(self.mass_use_changed)
        gui.use_st_radius.toggled.connect(self.radius_use_changed)
        gui.use_rho.toggled.connect(self.density_use_changed)
        gui.Stellar_mass_input.valueChanged.connect(self.mass_changed)
        gui.Stellar_radius_input.valueChanged.connect(self.radius_changed)
        gui.stellar_mass_radius_mode.toggled.connect(self.radius_mode_changed)
        gui.planet_radius_stellar_source.currentIndexChanged.connect(self.radius_source_changed)
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
            if field.endswith('_norm_pr'):
                other = field.replace('_norm_pr', '_jeff_pr')
            else:
                other = field.replace('_jeff_pr', '_norm_pr')
            for key, names, flag in self.groups:
                if key == other:
                    getattr(self.gui, flag).setChecked(False)

    def value_changed(self, value):
        session = self.session()
        if radius_mode(session):
            return
        session.stellar_rho = value
        slots = np.flatnonzero(session.use_planet)
        if session.rho and len(slots):
            index = int(slots[0])
            session.pl_a[index] = value
            widget = self.gui.param_gui_tr[index*3+2]
            blocker = QtCore.QSignalBlocker(widget)
            widget.setValue(value)
            del blocker
        sync_stellar(session)
        self.sync_stellar_controls()

    def planet_value_changed(self, slot, value):
        session = self.session()
        slots = np.flatnonzero(session.use_planet)
        if session.rho and len(slots) and slot == slots[0]:
            session.stellar_rho = value
            session.pl_a[slot] = value
            blocker = QtCore.QSignalBlocker(self.gui.Stellar_rho_input)
            self.gui.Stellar_rho_input.setValue(value)
            del blocker
            sync_stellar(session)
            self.sync_stellar_controls()

    def sync(self):
        session = self.session()
        ensure_stellar_settings(session)
        sync_stellar(session)
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
                    widget.setEnabled(not session.rho and not radius_mode(session))
                    widget.setToolTip('In density mode, use Limits and Priors > Stellar.' if session.rho else '')

        self.sync_stellar_controls()

    def mass_use_changed(self, checked):
        session = self.session()
        session.stellar_mass_use = bool(checked)
        session.use.use_stellar_mass = bool(checked)

    def radius_use_changed(self, checked):
        self.session().stellar_radius_use = bool(checked)

    def density_use_changed(self, checked):
        session = self.session()
        if checked and not session.rho:
            self.gui.radioButton_rho.setChecked(True)
        session.stellar_rho_use = bool(checked)
        blocker = QtCore.QSignalBlocker(self.gui.use_rho)
        self.gui.use_rho.setChecked(bool(checked))
        del blocker
        slots = np.flatnonzero(session.use_planet)
        if session.rho and len(slots):
            slot = int(slots[0])
            session.pl_a_use[slot] = bool(checked)
            widget = self.gui.use_param_gui_tr[slot * 3 + 2]
            blocker = QtCore.QSignalBlocker(widget)
            widget.setChecked(bool(checked))
            del blocker

    def mass_changed(self, value):
        session = self.session()
        session.stellar_mass = float(value)
        session.params.stellar_mass = float(value)
        sync_stellar(session)
        self.sync_stellar_controls()

    def radius_changed(self, value):
        session = self.session()
        if session.rho:
            return
        session.stellar_radius = float(value)
        sync_stellar(session)
        self.sync_stellar_controls()

    def radius_mode_changed(self, checked):
        if checked:
            self.gui.radioButton_stellar_mass_radius.setChecked(True)
        elif self.gui.radioButton_stellar_mass_radius.isChecked():
            self.gui.radioButton_rho.setChecked(True)

    def radius_source_changed(self, index):
        """Keep the corner-plot radius source independent of fitted errors."""
        self.session().planet_radius_stellar_source = ('auto', 'external', 'model')[index]

    def sync_stellar_controls(self):
        """Publish stellar fit flags, derived values and local uncertainties."""
        session = self.session()
        mode = radius_mode(session)
        source = getattr(session, 'planet_radius_stellar_source', 'auto')
        blocker = QtCore.QSignalBlocker(self.gui.planet_radius_stellar_source)
        self.gui.planet_radius_stellar_source.setCurrentIndex(('auto', 'external', 'model').index(source))
        del blocker
        slots = np.flatnonzero(session.use_planet)
        if session.rho and len(slots):
            session.stellar_rho_use = bool(session.pl_a_use[int(slots[0])])
        for name, checked in (
                ('use_st_mass', session.stellar_mass_use or session.use.use_stellar_mass),
                ('use_st_radius', session.stellar_radius_use if mode else False),
                ('use_rho', session.stellar_rho_use if session.rho else False),
                ('stellar_mass_radius_mode', mode)):
            widget = getattr(self.gui, name)
            blocker = QtCore.QSignalBlocker(widget)
            widget.setChecked(checked)
            del blocker
        self.gui.use_st_mass.setEnabled(True)
        self.gui.use_st_radius.setEnabled(mode)
        self.gui.use_rho.setEnabled(not mode)
        self.gui.Stellar_radius_input.setEnabled(not session.rho)
        self.gui.Stellar_rho_input.setEnabled(not mode)
        for root, control, error in (
                ('stellar_mass', 'Stellar_mass_input', 'err_st_mass'),
                ('stellar_radius', 'Stellar_radius_input', 'err_st_radius'),
                ('stellar_rho', 'Stellar_rho_input', 'err_rho')):
            widget = getattr(self.gui, control)
            blocker = QtCore.QSignalBlocker(widget)
            widget.setValue(getattr(session, root))
            del blocker
            uncertainty = np.max(np.abs(getattr(session, root + '_err')))
            getattr(self.gui, error).setText('+/- %.3f' % uncertainty)
        for field, names, flag in self.groups:
            active = not (mode and field.startswith('stellar_rho'))
            if session.rho and field.startswith('stellar_radius'):
                active = False
            for name in names:
                getattr(self.gui, name).setEnabled(active)
            if flag:
                getattr(self.gui, flag).setEnabled(active)
        if mode:
            for slot in range(9):
                self.gui.param_gui_tr[slot * 3 + 2].setEnabled(False)
                self.gui.use_param_gui_tr[slot * 3 + 2].setEnabled(False)
