"""Consistent color selection for Qt and Matplotlib plot colors."""
from PyQt6 import QtGui, QtWidgets
from matplotlib.colors import to_hex


def choose_color(parent, initial, title):
    """Return a QColor; cancellation returns an invalid color without mutation."""
    try:
        initial = to_hex(initial)
    except (ValueError, TypeError):
        pass
    color = QtGui.QColor(initial)
    if not color.isValid():
        color = QtGui.QColor('#0066ff')
    # Plot colors are stored as RGB, so do not offer an ignored alpha channel.
    return QtWidgets.QColorDialog.getColor(
        color, parent, title,
        QtWidgets.QColorDialog.ColorDialogOption.DontUseNativeDialog,
    )
