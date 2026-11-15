"""Shared file browser with visible detail columns for data imports."""
from PyQt6 import QtCore, QtWidgets

from file_open_ui import Ui_FileOpenDialog


def open_data_file(parent, caption, directory='', file_filter=''):
    """Return the selected file/filter using a screen-sized Designer dialog."""
    return _select_file(parent, caption, directory, file_filter, save=False)


def save_data_file(parent, caption, directory='', file_filter=''):
    """Choose an output filename with the same large, sortable browser."""
    return _select_file(parent, caption, directory, file_filter, save=True)


def _select_file(parent, caption, directory, file_filter, save):
    """Configure open/save behavior without changing the Designer layout."""
    dialog = QtWidgets.QFileDialog(parent)
    Ui_FileOpenDialog().setupUi(dialog)
    dialog.setWindowTitle(caption)
    if save:
        dialog.setAcceptMode(QtWidgets.QFileDialog.AcceptMode.AcceptSave)
        dialog.setFileMode(QtWidgets.QFileDialog.FileMode.AnyFile)
    if directory:
        location = QtCore.QFileInfo(directory)
        if location.isDir():
            dialog.setDirectory(directory)
        else:
            dialog.setDirectory(location.absolutePath())
            dialog.selectFile(location.fileName())
    if file_filter:
        dialog.setNameFilter(file_filter)

    # Keep every detail column visible, even if Qt restored an oversized Name
    # column from an earlier file browser. Sorting remains under user control.
    tree = dialog.findChild(QtWidgets.QTreeView, 'treeView')
    if tree is not None:
        header = tree.header()
        header.setStretchLastSection(False)
        header.setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeMode.Stretch)

    # Respect small screens without replacing the Designer default size.
    screen = parent.screen() if parent is not None else dialog.screen()
    if screen is not None:
        available = screen.availableGeometry()
        dialog.resize(min(dialog.width(), max(1, available.width() - 40)),
                      min(dialog.height(), max(1, available.height() - 80)))
    splitter = dialog.findChild(QtWidgets.QSplitter, 'splitter')
    if splitter is not None:
        sidebar_width = min(180, dialog.width() // 5)
        splitter.setSizes([sidebar_width, dialog.width() - sidebar_width])
    if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
        selected = dialog.selectedFiles()
        if selected:
            return selected[0], dialog.selectedNameFilter()
    return '', ''
