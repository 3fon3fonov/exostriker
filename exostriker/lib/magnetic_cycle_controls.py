"""Keep the magnetic-cycle form in sync with the current session."""

from magnetic_cycle import ensure, setting, GLOBAL


class MagneticCycleControls:
    """Connect the Designer controls to the cycle parameters."""

    def __init__(self, view, session):
        self.view = view
        self.session = session
        self.updating = False

        view.cycle_readme.clicked.connect(
            lambda: view.dialog_GP_help.show_help('magnetic_cycle'))
        view.cycle_dataset.currentIndexChanged.connect(self.refresh)
        view.cycle_enabled.toggled.connect(self.store)
        view.cycle_activity_enabled.toggled.connect(self.store)
        view.cycle_reference_epoch.valueChanged.connect(self.store)
        view.tabWidget_14.currentChanged.connect(self.refresh)

        for name in (*GLOBAL, 'amplitude', 'phase_offset', 'offset', 'jitter'):
            for field in ('value', 'lower', 'upper', 'prior1', 'prior2'):
                control = getattr(view, f'cycle_{name}_{field}')
                control.valueChanged.connect(self.store)

            fit_control = getattr(view, f'cycle_{name}_fit')
            fit_control.toggled.connect(self.store)

            prior_control = getattr(view, f'cycle_{name}_prior')
            prior_control.currentIndexChanged.connect(self.prior_changed)

    def settings(self):
        """Return the global settings and the displayed activity settings."""
        fit = self.session()
        ensure(fit)
        slot = self.view.cycle_dataset.currentData()

        parameters = {}
        for name in GLOBAL:
            parameters[name] = setting(fit, None, name)

        if slot in fit.magnetic_cycle['indicators'] and len(fit.act_data_sets.get(slot, [])):
            for name in ('amplitude', 'phase_offset', 'offset', 'jitter'):
                parameters[name] = setting(fit, slot, name)

        return parameters

    def refresh(self, *args):
        """Show session values without writing GUI changes back to the fit."""
        if self.updating:
            return

        self.updating = True

        try:
            view = self.view
            fit = self.session()
            cycle = ensure(fit)

            # Keep the selected file when rebuilding the loaded-file list.
            slot = view.cycle_dataset.currentData()
            view.cycle_dataset.clear()

            for index, data in fit.act_data_sets.items():
                if len(data):
                    view.cycle_dataset.addItem(f'{index+1}: {data[-1]}', index)

            index = view.cycle_dataset.findData(slot)
            if index >= 0:
                view.cycle_dataset.setCurrentIndex(index)

            if not view.cycle_dataset.count():
                view.cycle_dataset.addItem('Load an activity file first', None)

            slot = view.cycle_dataset.currentData()
            loaded = slot is not None

            view.cycle_dataset.setEnabled(loaded)
            view.cycle_activity_enabled.setEnabled(loaded)
            view.cycle_activity_enabled.setChecked(
                loaded and cycle['indicators'][slot]['enabled'])
            view.cycle_enabled.setChecked(cycle['enabled'])
            view.cycle_reference_epoch.setValue(cycle['reference_epoch'])

            parameters = self.settings()

            for name in (*GLOBAL, 'amplitude', 'phase_offset', 'offset', 'jitter'):
                parameter = parameters.get(name)
                error_label = getattr(view, f'cycle_{name}_error')
                error_label.setEnabled(parameter is not None)
                error_label.clear()

                for field in ('value', 'fit', 'lower', 'upper', 'prior', 'prior1', 'prior2'):
                    control = getattr(view, f'cycle_{name}_{field}')
                    control.setEnabled(parameter is not None)

                if parameter is None:
                    continue

                values = [('value', parameter['value']),
                          ('lower', parameter['bounds'][0]),
                          ('upper', parameter['bounds'][1])]

                for field, value in values:
                    getattr(view, f'cycle_{name}_{field}').setValue(value)

                getattr(view, f'cycle_{name}_fit').setChecked(parameter['fit'])

                # The last two boxes belong to the selected prior family.
                if parameter['normal'][2]:
                    prior = 1
                elif parameter['jeffreys'][2]:
                    prior = 2
                else:
                    prior = 0

                getattr(view, f'cycle_{name}_prior').setCurrentIndex(prior)
                prior_values = parameter['jeffreys' if prior == 2 else 'normal']

                for index, field in enumerate(('prior1', 'prior2')):
                    control = getattr(view, f'cycle_{name}_{field}')
                    control.setValue(prior_values[index])
                    control.setEnabled(prior != 0)

                # Match the planet labels, but keep tiny RHK errors visible.
                errors = parameter.get('error', [0., 0.])
                uncertainty = max(abs(value) for value in errors)

                if 0 < uncertainty < .0005:
                    display = f'{uncertainty:.3e}'
                else:
                    display = f'{uncertainty:.3f}'

                error_label.setText('+/- ' + display)
                getattr(view, f'cycle_{name}_value').setToolTip(
                    f'Last uncertainty: -{errors[0]:.12g} +{errors[1]:.12g}')

            names = []
            for index, parameter in cycle['indicators'].items():
                if parameter['enabled'] and len(fit.act_data_sets.get(index, [])):
                    names.append(str(fit.act_data_sets[index][-1]))

            if names:
                view.cycle_selection_summary.setText('Cycle indices: ' + ', '.join(names))
            else:
                view.cycle_selection_summary.setText(
                    'No activity selected: cycle fitted to RVs only.')

        finally:
            self.updating = False

    def prior_changed(self, *args):
        """Restore the saved values for the newly selected prior family."""
        if self.updating:
            return

        self.updating = True

        try:
            for name, parameter in self.settings().items():
                index = getattr(self.view, f'cycle_{name}_prior').currentIndex()
                values = parameter['jeffreys' if index == 2 else 'normal']

                for index, field in enumerate(('prior1', 'prior2')):
                    getattr(self.view, f'cycle_{name}_{field}').setValue(values[index])

        finally:
            self.updating = False

        self.store()

    def store(self, *args):
        """Save edited controls and refresh the shared activity nuisance fields."""
        if self.updating:
            return

        view = self.view
        fit = self.session()
        cycle = ensure(fit)

        cycle['enabled'] = view.cycle_enabled.isChecked()
        cycle['reference_epoch'] = view.cycle_reference_epoch.value()

        slot = view.cycle_dataset.currentData()
        if slot in cycle['indicators']:
            cycle['indicators'][slot]['enabled'] = view.cycle_activity_enabled.isChecked()

        for name, parameter in self.settings().items():
            parameter['value'] = getattr(view, f'cycle_{name}_value').value()
            parameter['fit'] = getattr(view, f'cycle_{name}_fit').isChecked()
            parameter['bounds'] = [getattr(view, f'cycle_{name}_{field}').value()
                                   for field in ('lower', 'upper')]

            prior = getattr(view, f'cycle_{name}_prior').currentIndex()
            parameter['normal'][2] = prior == 1
            parameter['jeffreys'][2] = prior == 2

            if prior:
                family = 'jeffreys' if prior == 2 else 'normal'
                parameter[family][:2] = [getattr(view, f'cycle_{name}_{field}').value()
                                         for field in ('prior1', 'prior2')]

        self.refresh()

        other = getattr(view, 'activity_gp_controls', None)
        if other is not None:
            other.refresh()
