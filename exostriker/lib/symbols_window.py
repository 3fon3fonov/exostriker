"""Choose a plot marker with previews of the actual pyqtgraph paths."""
import sys
from PyQt6 import QtCore, QtGui, QtWidgets
from pyqtgraph.graphicsItems.ScatterPlotItem import Symbols


SYMBOL_OPTIONS = (
    ('o', 'Circle'), ('t', 'Triangle up'), ('t1', 'Triangle down'),
    ('t2', 'Triangle right'), ('t3', 'Triangle left'), ('s', 'Square'),
    ('p', 'Pentagon'), ('h', 'Hexagon'), ('star', 'Star'),
    ('+', 'Plus'), ('d', 'Diamond'),
)


def symbol_icon(symbol, palette):
    pixmap = QtGui.QPixmap(40, 40)
    pixmap.fill(QtCore.Qt.GlobalColor.transparent)
    painter = QtGui.QPainter(pixmap)
    painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
    painter.translate(20, 20)
    painter.scale(24, -24)
    pen = QtGui.QPen(palette.color(QtGui.QPalette.ColorRole.Text))
    pen.setWidthF(0.06)
    painter.setPen(pen)
    painter.setBrush(palette.color(QtGui.QPalette.ColorRole.Highlight))
    painter.drawPath(Symbols[symbol])
    painter.end()
    return QtGui.QIcon(pixmap)


class show_symbols(QtWidgets.QDialog):
    def __init__(self, parent=None, current='o'):
        super().__init__(parent)
        self.setWindowTitle('Plot symbol options')
        self.layout = QtWidgets.QVBoxLayout(self)
        label = QtWidgets.QLabel('Choose a marker for this dataset.', self)
        self.layout.addWidget(label)
        self.radio_group = QtWidgets.QButtonGroup(self)
        grid = QtWidgets.QGridLayout()
        for index, (symbol, name) in enumerate(SYMBOL_OPTIONS, 1):
            button = QtWidgets.QRadioButton(f'{name} ({symbol})', self)
            button.setIcon(symbol_icon(symbol, self.palette()))
            button.setIconSize(QtCore.QSize(40, 40))
            button.setToolTip(f'pyqtgraph marker: {symbol}')
            self.radio_group.addButton(button, index)
            setattr(self, f'button{index}', button)
            grid.addWidget(button, (index - 1) % 6, (index - 1) // 6)
            button.setChecked(symbol == current)
        if self.radio_group.checkedId() == -1:
            self.button1.setChecked(True)
        self.layout.addLayout(grid)
        self.buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok |
            QtWidgets.QDialogButtonBox.StandardButton.Cancel, parent=self)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        self.layout.addWidget(self.buttons)

    def return_but_N(self):
        return self.radio_group.checkedId()

    @staticmethod
    def get_radio(parent=None, current='o'):
        dialog = show_symbols(parent, current)
        if dialog.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            return dialog.return_but_N()
        return None


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    show_symbols.get_radio()
