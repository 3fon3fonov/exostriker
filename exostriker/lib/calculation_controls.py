"""Keep display controls usable while a private scientific job is running."""
import copy
from PyQt6 import QtWidgets

# These actions alter data/parameters or launch another numerical job, even
# though their buttons happen to live alongside a plot.
_LOCKED_BUTTONS = {
    'adopt_best_RV__o_c_GLS_per',
    'adopt_mpl_param',
    'calc_MLP',
    'calc_TLS',
    'calc_TLS_o_c',
    'adopt_tls_o_c_param',
    'button_make_mcmc_cornerplot',
    'button_make_nest_cornerplot',
    'remove_mcmc_samples',
    'remove_ns_samples',
}
_DISPLAY_FIELDS = (
    'nbody_plot_options',
    'ast_colors',
    'gls_colors',
    'rvs_colors',
    'tra_colors',
    'ttv_colors',
    'act_colors',
    'pyqt_color_alpha_act',
    'pyqt_symbols_act',
    'pyqt_symbols_size_act',
    'pyqt_color_alpha_rvs',
    'pyqt_color_alpha_tra',
    'pyqt_symbols_ast',
    'pyqt_symbols_rvs',
    'pyqt_symbols_tra',
    'pyqt_symbols_ttv',
    'pyqt_symbols_size_ast',
    'pyqt_symbols_size_rvs',
    'pyqt_symbols_size_tra',
    'pyqt_symbols_size_ttv',
)


def lockable_controls(gui, controls):
    """Exclude presentation controls; leave their existing enabled state intact."""
    roots = [getattr(gui, name, None) for name in ('plot_tabs', 'plot_options_tabs')]
    roots = [root for root in roots if isinstance(root, QtWidgets.QWidget)]
    locked = []
    for control in controls:
        if isinstance(control, QtWidgets.QMenuBar):
            # Display/export/help menus stay usable; session/data edits do not.
            locked.extend(action.menu() for action in control.actions()
                          if action.menu() is not None
                          and action.menu().objectName() not in ('menuview','menuHelp','menuSettings'))
            continue
        if isinstance(control, QtWidgets.QWidget):
            in_plot = any((root is control or root.isAncestorOf(control) for root in roots))
            if in_plot and control.objectName() not in _LOCKED_BUTTONS:
                # Include editors inside a plot spin box/combo box as well.
                locked_parent = control.parentWidget()
                if locked_parent is None or locked_parent.objectName() not in _LOCKED_BUTTONS:
                    continue
        locked.append(control)
    return locked


def preserve_display_settings(source, result):
    """Adopting fitted parameters must not revert colors edited during the job."""
    for name in _DISPLAY_FIELDS:
        if hasattr(source, name):
            setattr(result, name, copy.deepcopy(getattr(source, name)))
