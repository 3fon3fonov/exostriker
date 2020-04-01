import sys

from PyQt5.QtWidgets import *
from PyQt5.QtCore    import *
from PyQt5.QtGui     import *


#font = QtGui.QFont()
#font.setPointSize(8)
#font.setBold(False)

class RVBank_window(QDialog):


    def __init__(self, parent=None):
        super(RVBank_window, self).__init__(parent=parent)
        vLayout = QVBoxLayout(self)
        hLayout = QHBoxLayout()

        self.lineEdit = QLineEdit(self)
        hLayout.addWidget(self.lineEdit)    

        self.filter = QPushButton("Search", self)
        hLayout.addWidget(self.filter)
        self.filter.clicked.connect(self.filterClicked)

        self.list = QListView(self)

        vLayout.addLayout(hLayout)
        vLayout.addWidget(self.list)

        self.model = QStandardItemModel(self.list)

        codes = [
            'Test',
            'HD27894',
            'GJ436',
            'GJ1111',
            'BD-102030'
        ]

        for code in codes:
            item = QStandardItem(code)
            item.setCheckable(True)
            self.model.appendRow(item)
        self.list.setModel(self.model)

    def filterClicked(self):
        filter_text = str(self.lineEdit.text()).lower()
        for row in range(self.model.rowCount()):
            if filter_text in str(self.model.item(row).text()).lower():
                self.list.setRowHidden(row, False)
            else:
                self.list.setRowHidden(row, True)

 