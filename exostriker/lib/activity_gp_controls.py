"""Session synchronization for the joint activity GP form defined in es.ui."""
import copy
from activity_gp import ensure, FIELDS, DIFFERENTIABLE


class ActivityGPControls:
    def __init__(self, view, session):
        self.view, self.session = view, session
        self.updating = False
        view.activity_gp_open.clicked.connect(self.open_selection)
        view.activity_gp_readme.clicked.connect(self.show_readme)
        view.activity_gp_dataset.currentIndexChanged.connect(self.refresh)
        view.activity_gp_enabled.toggled.connect(self.store)
        for name in (*FIELDS, 'rv_derivative'):
            for field in ('value', 'lower', 'upper', 'prior1', 'prior2'):
                getattr(view, f'activity_gp_{name}_{field}').valueChanged.connect(self.store)
            getattr(view, f'activity_gp_{name}_fit').toggled.connect(self.store)
            getattr(view, f'activity_gp_{name}_prior').currentIndexChanged.connect(self.prior_changed)

    def settings(self, slot):
        session = self.session()
        ensure(session)
        data = session.activity_gp.get(slot)
        if data is None:
            return {}
        return dict(**{name: data[name] for name in FIELDS}, rv_derivative=session.activity_gp_rv_derivative)

    def refresh(self, *args):
        if self.updating:
            return
        self.updating = True
        try:
            view, session = self.view, self.session()
            ensure(session)
            derivatives_supported = session.gp_kernel in DIFFERENTIABLE
            # Retain the user's differentiable-kernel setup across kernel
            # switches, but leave a valid zero-derivative model for Rotation/DRW.
            derivatives = [session.activity_gp_rv_derivative]
            derivatives += [p['derivative'] for p in session.activity_gp.values()]
            for parameter in derivatives:
                if derivatives_supported:
                    saved = parameter.pop('_differentiable_settings', None)
                    if saved is not None:
                        parameter.update(saved)
                else:
                    if '_differentiable_settings' not in parameter:
                        parameter['_differentiable_settings'] = copy.deepcopy(parameter)
                    parameter.update(value=0., fit=False)
                    parameter['bounds'] = [min(parameter['bounds'][0], 0.),
                                           max(parameter['bounds'][1], 0.)]
                    parameter['normal'][2] = parameter['jeffreys'][2] = False
            slot = view.activity_gp_dataset.currentData()
            view.activity_gp_dataset.clear()
            loaded = [i for i, row in session.act_data_sets.items() if len(row)]
            for i in loaded:
                row = session.act_data_sets[i]
                view.activity_gp_dataset.addItem(f'{i+1}: {row[-1]}', i)
            if loaded:
                if slot not in loaded:
                    slot = next((i for i in loaded if session.activity_gp[i]['enabled']), loaded[0])
                view.activity_gp_dataset.setCurrentIndex(view.activity_gp_dataset.findData(slot))
            else:
                slot = None
                view.activity_gp_dataset.addItem('Load an activity file under Activity data first', None)
            view.activity_gp_dataset.setEnabled(bool(loaded))
            settings = self.settings(slot)
            view.activity_gp_enabled.setEnabled(bool(settings))
            view.activity_gp_enabled.setChecked(bool(settings and session.activity_gp[slot]['enabled']))
            for name in (*FIELDS, 'rv_derivative'):
                available = bool(settings) and (derivatives_supported or name not in ('derivative', 'rv_derivative'))
                getattr(view, f'activity_gp_{name}_label').setEnabled(available)
                for field in ('value', 'lower', 'upper', 'fit', 'prior', 'prior1', 'prior2'):
                    getattr(view, f'activity_gp_{name}_{field}').setEnabled(available)
                if not settings:
                    continue
                parameter = settings[name]
                getattr(view, f'activity_gp_{name}_value').setValue(parameter['value'])
                getattr(view, f'activity_gp_{name}_fit').setChecked(parameter['fit'])
                getattr(view, f'activity_gp_{name}_lower').setValue(parameter['bounds'][0])
                getattr(view, f'activity_gp_{name}_upper').setValue(parameter['bounds'][1])
                prior = 'normal' if parameter['normal'][2] else 'jeffreys' if parameter['jeffreys'][2] else None
                getattr(view, f'activity_gp_{name}_prior').setCurrentIndex({None: 0, 'normal': 1, 'jeffreys': 2}[prior])
                values = parameter[prior or 'normal']
                getattr(view, f'activity_gp_{name}_prior1').setValue(values[0])
                getattr(view, f'activity_gp_{name}_prior2').setValue(values[1])
                for field in ('prior1', 'prior2'):
                    getattr(view, f'activity_gp_{name}_{field}').setEnabled(available and prior is not None)
                error = parameter.get('error', [0., 0.])
                getattr(view, f'activity_gp_{name}_value').setToolTip(f'Last uncertainty: -{error[0]:.6g} +{error[1]:.6g}')
            enabled = [str(i+1)+': '+str(session.act_data_sets[i][-1]) for i, p in session.activity_gp.items()
                       if p['enabled'] and len(session.act_data_sets[i])]
            view.activity_gp_selection_summary.setText('Joint RV GP: '+', '.join(enabled) if enabled else
                'No activity selected: use the existing RV-only time GP.')
        finally:
            self.updating = False

    def open_selection(self):
        self.refresh()
        self.view.tabWidget_14.setCurrentWidget(self.view.tab_activity_joint_gp)

    def show_readme(self):
        self.view.dialog_GP_help.show_help('joint_activity_gp')

    def prior_changed(self, *args):
        if self.updating:
            return
        # Do not reinterpret Gaussian mean/sigma as Jeffreys lower/upper.
        view = self.view
        slot = view.activity_gp_dataset.currentData()
        settings = self.settings(slot)
        if not settings:
            return
        self.updating = True
        try:
            for name, setting in settings.items():
                index = getattr(view, f'activity_gp_{name}_prior').currentIndex()
                chosen = setting['normal' if index != 2 else 'jeffreys']
                getattr(view, f'activity_gp_{name}_prior1').setValue(chosen[0])
                getattr(view, f'activity_gp_{name}_prior2').setValue(chosen[1])
        finally:
            self.updating = False
        self.store()

    def store(self, *args):
        if self.updating:
            return
        view, session = self.view, self.session()
        slot = view.activity_gp_dataset.currentData()
        settings = self.settings(slot)
        if not settings:
            return
        enabled = view.activity_gp_enabled.isChecked()
        was_enabled = session.activity_gp[slot]['enabled']
        session.activity_gp[slot]['enabled'] = enabled
        if enabled and not was_enabled:
            view.do_RV_GP.setChecked(True)
        for name, setting in settings.items():
            setting['value'] = getattr(view, f'activity_gp_{name}_value').value()
            setting['fit'] = getattr(view, f'activity_gp_{name}_fit').isChecked()
            setting['bounds'] = [getattr(view, f'activity_gp_{name}_{field}').value() for field in ('lower', 'upper')]
            index = getattr(view, f'activity_gp_{name}_prior').currentIndex()
            setting['normal'][2] = index == 1
            setting['jeffreys'][2] = index == 2
            if index:
                chosen = setting['normal' if index == 1 else 'jeffreys']
                chosen[:2] = [getattr(view, f'activity_gp_{name}_{field}').value() for field in ('prior1', 'prior2')]
        self.refresh()
        other = getattr(view, 'magnetic_cycle_controls', None)
        if other is not None:
            other.refresh()
