"""Reusable, searchable information window for help and fit reports."""
import sys
from PyQt6 import QtCore, QtGui, QtWidgets


class print_info(QtWidgets.QMainWindow):
    def __init__(self, parent=None, *, show_tools=True):
        super().__init__(parent)
        self.setWindowTitle('Exo-Striker information')
        self.resize(720, 640)
        self.setMinimumSize(400, 300)
        self.widget = QtWidgets.QWidget(self)
        layout = QtWidgets.QVBoxLayout(self.widget)
        layout.setContentsMargins(16, 12, 16, 12)
        self.text = QtWidgets.QTextBrowser(self.widget)
        self.text.setOpenExternalLinks(True)
        self.text.setReadOnly(True)
        font = self.text.font()
        font.setPointSize(max(10, font.pointSize()))
        self.text.setFont(font)
        self.text.document().setDefaultStyleSheet(
            'h1 { font-size: 20pt; margin-bottom: 12px; } '
            'h2 { font-size: 13pt; margin-top: 18px; } '
            'p, li { line-height: 140%; } '
            'pre { margin: 12px; } th { font-weight: bold; }'
        )

        toolbar = QtWidgets.QHBoxLayout()
        self.find_field = QtWidgets.QLineEdit(self.widget)
        self.find_field.setPlaceholderText('Find in this window…')
        self.find_field.setClearButtonEnabled(True)
        self.find_field.setAccessibleName('Find in this window')
        toolbar.addWidget(self.find_field, 1)
        self.previous_button = QtWidgets.QPushButton('Previous', self.widget)
        self.next_button = QtWidgets.QPushButton('Next', self.widget)
        self.zoom_out_button = QtWidgets.QPushButton('A−', self.widget)
        self.zoom_in_button = QtWidgets.QPushButton('A+', self.widget)
        self.zoom_out_button.setToolTip('Decrease text size')
        self.zoom_in_button.setToolTip('Increase text size')
        self.zoom_out_button.setAccessibleName('Decrease text size')
        self.zoom_in_button.setAccessibleName('Increase text size')
        for button in (self.previous_button, self.next_button,
                       self.zoom_out_button, self.zoom_in_button):
            toolbar.addWidget(button)
        layout.addLayout(toolbar)
        layout.addWidget(self.text, 1)
        footer = QtWidgets.QHBoxLayout()
        self.find_status = QtWidgets.QLabel(self.widget)
        footer.addWidget(self.find_status, 1)
        self.close_button = QtWidgets.QPushButton('Close', self.widget)
        footer.addWidget(self.close_button)
        layout.addLayout(footer)
        self.setCentralWidget(self.widget)

        self.find_field.textChanged.connect(self._find_changed)
        self.find_field.returnPressed.connect(self.find_next)
        self.next_button.clicked.connect(self.find_next)
        self.previous_button.clicked.connect(lambda: self.find_next(backward=True))
        self.zoom_in_button.clicked.connect(lambda: self.text.zoomIn(1))
        self.zoom_out_button.clicked.connect(lambda: self.text.zoomOut(1))
        self.close_button.clicked.connect(self.close)
        self._find_changed('')
        if not show_tools:
            for control in (self.find_field, self.previous_button, self.next_button,
                            self.zoom_out_button, self.zoom_in_button, self.find_status):
                control.hide()
        for key, callback in (('Ctrl+F', self.focus_find), ('F3', self.find_next),
                              ('Shift+F3', lambda: self.find_next(backward=True)),
                              ('Escape', self.close)):
            if not show_tools and key != 'Escape':
                continue
            shortcut = QtGui.QShortcut(QtGui.QKeySequence(key), self)
            shortcut.activated.connect(callback)

    def focus_find(self):
        self.find_field.setFocus()
        self.find_field.selectAll()

    def _find_changed(self, query):
        self.previous_button.setEnabled(bool(query))
        self.next_button.setEnabled(bool(query))
        self.find_status.clear()
        cursor = self.text.textCursor()
        cursor.clearSelection()
        cursor.movePosition(QtGui.QTextCursor.MoveOperation.Start)
        self.text.setTextCursor(cursor)
        if query:
            self.find_next()

    def find_next(self, checked=False, *, backward=False):
        query = self.find_field.text()
        if not query:
            return
        flags = QtGui.QTextDocument.FindFlag.FindBackward if backward else QtGui.QTextDocument.FindFlag(0)
        if not self.text.find(query, flags):
            cursor = self.text.textCursor()
            cursor.movePosition(QtGui.QTextCursor.MoveOperation.End if backward
                                else QtGui.QTextCursor.MoveOperation.Start)
            self.text.setTextCursor(cursor)
            if not self.text.find(query, flags):
                self.find_status.setText('No matches')
                return
        self.find_status.clear()

    def show_help(self, topic):
        from help_content import HELP_PAGES
        title, html = HELP_PAGES[topic]
        self.setWindowTitle(title)
        self.resize(760, 680)
        self.find_field.clear()
        self.text.setHtml(html)
        self.show()

    def showEvent(self, event):
        super().showEvent(event)
        self.find_field.clear()
        self.find_status.clear()
        self.text.moveCursor(QtGui.QTextCursor.MoveOperation.Start)
        # Keep a newly opened helper within the available screen area.
        screen = self.screen()
        if screen:
            available = screen.availableGeometry()
            self.resize(min(self.width(), available.width() - 40),
                        min(self.height(), available.height() - 60))


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    main = print_info()
    main.show_help('ttv')
    sys.exit(app.exec())
