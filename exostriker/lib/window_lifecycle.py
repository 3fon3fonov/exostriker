"""Close Exo-Striker helper windows after the main quit is accepted."""
from PyQt6 import QtCore, QtWidgets


def close_auxiliary_windows(main):
    app = QtWidgets.QApplication.instance()
    if app is None:
        return
    referenced = {id(value) for value in vars(main).values()
                  if isinstance(value, QtWidgets.QWidget)}
    for window in list(app.topLevelWidgets()):
        if window is main:
            continue
        owner = window
        seen = set()
        belongs = id(window) in referenced
        while owner is not None and id(owner) not in seen:
            if owner is main:
                belongs = True
                break
            seen.add(id(owner))
            # Some older dialogs shadow QObject.parent with a Python attribute.
            legacy_parent = getattr(owner, '__dict__', {}).get('parent')
            owner = legacy_parent if isinstance(legacy_parent, QtWidgets.QWidget) else QtCore.QObject.parent(owner)
        if belongs:
            window.setProperty('exostriker_shutdown', True)
            window.close()
            # A legacy close handler may reject closing; it must not leave a
            # visible orphan after the user has accepted the main quit prompt.
            window.hide()
