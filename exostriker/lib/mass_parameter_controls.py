"""Connect the separate Designer mass editors and priors to session state."""
import numpy as np
from PyQt6 import QtCore, QtWidgets
from mass_parameterization import ensure_mass_settings


class MassParameterControls:
    """Keep fitted and derived orbital amplitudes visible in their own rows."""

    def __init__(self, gui, session):
        self.gui = gui
        self.session = session
        gui.mp_bounds_gui = []
        gui.mp_norm_pr_gui = []
        gui.mp_jeff_pr_gui = []
        for slot in range(9):
            index = str(slot + 1)
            gui.mp_bounds_gui.append([
                getattr(gui, 'mp_min_' + index), getattr(gui, 'mp_max_' + index)])
            gui.mp_norm_pr_gui.append([
                getattr(gui, 'mp_mean_' + index), getattr(gui, 'mp_sigma_' + index),
                getattr(gui, 'use_mp_norm_pr_' + index)])
            gui.mp_jeff_pr_gui.append([
                getattr(gui, 'mp_jeff_alpha_' + index), getattr(gui, 'mp_jeff_beta_' + index),
                getattr(gui, 'use_mp_jeff_pr_' + index)])
            for group in (gui.mp_bounds_gui[-1], gui.mp_norm_pr_gui[-1], gui.mp_jeff_pr_gui[-1]):
                for widget in group[:2]:
                    widget.valueChanged.connect(self.read_priors)
                if len(group) == 3:
                    group[2].toggled.connect(self.read_priors)
            for prefix in ('use_K', 'use_mp'):
                getattr(gui, prefix + index).toggled.connect(
                    lambda checked, i=slot: self.read_use(i, checked))
            for name in ('K', 'mp', 'P', 'e', 'om', 'ma', 'incl', 'Omega'):
                getattr(gui, name + index).editingFinished.connect(self.refresh_derived)
        gui.Stellar_mass_input.editingFinished.connect(self.refresh_derived)
        self.sync()

    @staticmethod
    def set_value(widget, value):
        """Update a control without publishing intermediate model values."""
        blocker = QtCore.QSignalBlocker(widget)
        widget.setValue(float(value))
        del blocker

    @staticmethod
    def set_checked(widget, checked):
        blocker = QtCore.QSignalBlocker(widget)
        widget.setChecked(bool(checked))
        del blocker

    def read_use(self, slot, checked):
        """Link the two fit flags while only the chosen coordinate is editable."""
        session = self.session()
        session.K_use[slot] = bool(checked)
        session.mp_use[slot] = bool(checked)
        for prefix in ('use_K', 'use_mp'):
            self.set_checked(getattr(self.gui, prefix + str(slot+1)), checked)

    def read_priors(self, *args):
        """Store mass priors independently of every K prior control."""
        session = self.session()
        ensure_mass_settings(session)
        for slot in range(9):
            session.mp_bound[slot] = np.array([
                widget.value() for widget in self.gui.mp_bounds_gui[slot]])
            for field, groups in (
                    ('mp_norm_pr', self.gui.mp_norm_pr_gui),
                    ('mp_jeff_pr', self.gui.mp_jeff_pr_gui)):
                first, second, enabled = groups[slot]
                getattr(session, field)[slot] = np.array([
                    first.value(), second.value(), enabled.isChecked()])

    def sync_priors(self):
        """Keep mass priors editable and disable inactive K prior selection."""
        session = self.session()
        ensure_mass_settings(session)
        mass_mode = session.mass_parameterization
        for slot in range(9):
            for field, groups in (
                    ('mp_bound', self.gui.mp_bounds_gui),
                    ('mp_norm_pr', self.gui.mp_norm_pr_gui),
                    ('mp_jeff_pr', self.gui.mp_jeff_pr_gui)):
                values = getattr(session, field)[slot]
                for widget, value in zip(groups[slot][:2], values[:2]):
                    self.set_value(widget, value)
                if len(groups[slot]) == 3:
                    self.set_checked(groups[slot][2], values[2])
                for widget in groups[slot]:
                    widget.setEnabled(True)
            for group in (self.gui.param_bounds_gui, self.gui.param_nr_priors_gui,
                          self.gui.param_jeff_priors_gui):
                for widget in group[10*slot]:
                    widget.setEnabled(not mass_mode)
            for prefix in ('use_K_norm_pr_', 'use_K_jeff_pr_'):
                getattr(self.gui, prefix + str(slot+1)).setEnabled(not mass_mode)
            for prefix in ('use_K_bound_', 'use_mp_bound_'):
                widget = getattr(self.gui, prefix + str(slot+1))
                widget.setEnabled(False)

    def sync_values(self):
        """Publish both amplitudes while preserving the fitted coordinate."""
        session = self.session()
        ensure_mass_settings(session)
        for slot in range(9):
            index = str(slot + 1)
            mass = session.mp[slot]
            mass_widget = getattr(self.gui, 'mp' + index)
            self.set_value(mass_widget, mass if np.isfinite(mass) else 0.)
            self.set_value(getattr(self.gui, 'K' + index), session.K[slot])
            selected = session.mp_use[slot] if session.mass_parameterization else session.K_use[slot]
            self.read_use(slot, selected)
            if session.mass_parameterization:
                error = max(np.abs(session.mp_err[slot]))
                getattr(self.gui, 'err_mp' + index).setText('+/- %.3f' % error)
                getattr(self.gui, 'err_K' + index).setText('+/- 0.000')
            else:
                getattr(self.gui, 'err_mp' + index).setText('+/- 0.000')
                error = max(np.abs(session.K_err[slot]))
                getattr(self.gui, 'err_K' + index).setText('+/- %.3f' % error)
            amplitude_available = session.mod_dynamical or any(
                session.type_fit.get(kind, False) for kind in ('RV', 'TTV', 'AST'))
            mass_widget.setEnabled(amplitude_available and session.mass_parameterization)
            getattr(self.gui, 'use_mp' + index).setEnabled(amplitude_available and session.mass_parameterization)
            getattr(self.gui, 'K' + index).setEnabled(amplitude_available and not session.mass_parameterization)
            getattr(self.gui, 'use_K' + index).setEnabled(amplitude_available and not session.mass_parameterization)

    def sync(self):
        """Restore the selected mode, labels, values and independent priors."""
        session = self.session()
        ensure_mass_settings(session)
        for name, checked in (
                ('radioButton_plmass', session.mass_parameterization),
                ('radioButton_K', not session.mass_parameterization)):
            self.set_checked(getattr(self.gui, name), checked)
        property_name = 'massText' if session.mod_dynamical else 'minimumMassText'
        for label in self.gui.findChildren(QtWidgets.QLabel):
            if label.objectName().startswith('label_mp_input'):
                label.setText(label.property(property_name))
        self.sync_values()
        self.sync_priors()

    def refresh_derived(self):
        """Recalculate the greyed-out value after an orbital edit is committed."""
        try:
            self.gui.update_params()
        except ValueError as exc:
            self.gui.statusBar().showMessage(str(exc), 10000)
            return
        self.sync_values()
