"""Controls in the existing Fortran Simplex and LM pages; values belong to sessions."""
import copy
from PyQt6 import QtCore
from optimizer_spinbox import OptimizerDoubleSpinBox
from fortran_optimizer_options import DEFAULTS, ensure_options



class FortranOptimizerControls:
    """Connect Designer-defined controls to the currently selected session."""
    def __init__(self, gui, session):
        self.session = session
        self.widgets = {}
        for method, keys in (
            ('simplex', ('maxfev','maxiter','ftol','xtol','step_scale','restarts','restart_ftol')),
            ('lm', ('maxfev','maxiter','ftol','xtol','gtol','initial_damping','diff_step_scale')),
        ):
            for key in keys:
                widget = getattr(gui, f'fortran_{method}_{key}')
                widget.valueChanged.connect(lambda value, m=method, k=key: self.changed(m,k,value))
                self.widgets[method,key] = widget
            getattr(gui, f'fortran_{method}_restore_defaults').clicked.connect(
                lambda checked=False, m=method: self.reset(m))
        self.sync()

    def sync(self):
        settings = ensure_options(self.session())
        for (method,key), widget in self.widgets.items():
            blocker = QtCore.QSignalBlocker(widget)
            widget.setValue(settings[method][key])
            del blocker

    def changed(self, method, key, value):
        settings = ensure_options(self.session())
        settings[method][key] = value

    def reset(self, method):
        ensure_options(self.session())[method] = copy.deepcopy(DEFAULTS[method])
        self.sync()
