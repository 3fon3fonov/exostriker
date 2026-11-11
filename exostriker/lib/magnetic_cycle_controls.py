"""Session synchronization for the magnetic-cycle form defined in es.ui."""
from magnetic_cycle import ensure, setting, GLOBAL


class MagneticCycleControls:
    def __init__(self, view, session):
        self.view, self.session = view, session
        self.updating = False
        view.cycle_readme.clicked.connect(lambda: view.dialog_GP_help.show_help('magnetic_cycle'))
        view.cycle_dataset.currentIndexChanged.connect(self.refresh)
        view.cycle_enabled.toggled.connect(self.store)
        view.cycle_activity_enabled.toggled.connect(self.store)
        view.cycle_reference_epoch.valueChanged.connect(self.store)
        view.tabWidget_14.currentChanged.connect(self.refresh)
        for name in (*GLOBAL, 'amplitude', 'phase_offset', 'offset', 'jitter'):
            for field in ('value','lower','upper','prior1','prior2'):
                getattr(view, f'cycle_{name}_{field}').valueChanged.connect(self.store)
            getattr(view, f'cycle_{name}_fit').toggled.connect(self.store)
            getattr(view, f'cycle_{name}_prior').currentIndexChanged.connect(self.prior_changed)

    def settings(self):
        s = self.session(); ensure(s)
        slot = self.view.cycle_dataset.currentData()
        result = {name: setting(s,None,name) for name in GLOBAL}
        if slot in s.magnetic_cycle['indicators'] and len(s.act_data_sets.get(slot, [])):
            result.update({name:setting(s,slot,name) for name in ('amplitude','phase_offset','offset','jitter')})
        return result

    def refresh(self, *args):
        if self.updating:return
        self.updating = True
        try:
            v, s = self.view, self.session(); c = ensure(s)
            slot = v.cycle_dataset.currentData()
            v.cycle_dataset.clear()
            for i,row in s.act_data_sets.items():
                if len(row):v.cycle_dataset.addItem(f'{i+1}: {row[-1]}',i)
            index = v.cycle_dataset.findData(slot)
            if index>=0:v.cycle_dataset.setCurrentIndex(index)
            if not v.cycle_dataset.count():v.cycle_dataset.addItem('Load an activity file first',None)
            slot = v.cycle_dataset.currentData(); loaded = slot is not None
            v.cycle_dataset.setEnabled(loaded)
            v.cycle_activity_enabled.setEnabled(loaded)
            v.cycle_activity_enabled.setChecked(loaded and c['indicators'][slot]['enabled'])
            v.cycle_enabled.setChecked(c['enabled'])
            v.cycle_reference_epoch.setValue(c['reference_epoch'])
            settings = self.settings()
            for name in (*GLOBAL,'amplitude','phase_offset','offset','jitter'):
                p = settings.get(name)
                error_label = getattr(v, f'cycle_{name}_error')
                error_label.setEnabled(p is not None)
                error_label.clear()
                for field in ('value','fit','lower','upper','prior','prior1','prior2'):
                    getattr(v,f'cycle_{name}_{field}').setEnabled(p is not None)
                if p is None:continue
                for field,value in [('value',p['value']),('lower',p['bounds'][0]),('upper',p['bounds'][1])]:
                    getattr(v,f'cycle_{name}_{field}').setValue(value)
                getattr(v,f'cycle_{name}_fit').setChecked(p['fit'])
                prior = 1 if p['normal'][2] else 2 if p['jeffreys'][2] else 0
                getattr(v,f'cycle_{name}_prior').setCurrentIndex(prior)
                values = p['jeffreys' if prior==2 else 'normal']
                for i,field in enumerate(('prior1','prior2')):
                    getattr(v,f'cycle_{name}_{field}').setValue(values[i])
                    getattr(v,f'cycle_{name}_{field}').setEnabled(prior!=0)
                error=p.get('error',[0.,0.])
                uncertainty = max(abs(value) for value in error)
                # Keep ordinary errors consistent with the planet fields, but
                # do not hide small, nonzero activity uncertainties as zero.
                display = f'{uncertainty:.3e}' if 0 < uncertainty < .0005 else f'{uncertainty:.3f}'
                error_label.setText('+/- '+display)
                getattr(v,f'cycle_{name}_value').setToolTip(f'Last uncertainty: -{error[0]:.12g} +{error[1]:.12g}')
            names = [str(s.act_data_sets[i][-1]) for i,p in c['indicators'].items() if p['enabled'] and len(s.act_data_sets.get(i,[]))]
            v.cycle_selection_summary.setText('Cycle indices: '+', '.join(names) if names else 'No activity selected: cycle fitted to RVs only.')
        finally:self.updating=False

    def prior_changed(self, *args):
        if self.updating:return
        self.updating=True
        try:
            for name,p in self.settings().items():
                index=getattr(self.view,f'cycle_{name}_prior').currentIndex()
                values=p['jeffreys' if index==2 else 'normal']
                for i,field in enumerate(('prior1','prior2')):getattr(self.view,f'cycle_{name}_{field}').setValue(values[i])
        finally:self.updating=False
        self.store()

    def store(self, *args):
        if self.updating:return
        v,s=self.view,self.session();c=ensure(s)
        c['enabled']=v.cycle_enabled.isChecked();c['reference_epoch']=v.cycle_reference_epoch.value()
        slot=v.cycle_dataset.currentData()
        if slot in c['indicators']:c['indicators'][slot]['enabled']=v.cycle_activity_enabled.isChecked()
        for name,p in self.settings().items():
            p['value']=getattr(v,f'cycle_{name}_value').value();p['fit']=getattr(v,f'cycle_{name}_fit').isChecked()
            p['bounds']=[getattr(v,f'cycle_{name}_{field}').value() for field in ('lower','upper')]
            prior=getattr(v,f'cycle_{name}_prior').currentIndex();p['normal'][2]=prior==1;p['jeffreys'][2]=prior==2
            if prior:p['jeffreys' if prior==2 else 'normal'][:2]=[getattr(v,f'cycle_{name}_{field}').value() for field in ('prior1','prior2')]
        self.refresh()
        other = getattr(v, 'activity_gp_controls', None)
        if other is not None:
            other.refresh()
