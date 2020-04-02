import sys
from PyQt5 import QtWidgets,QtGui,QtCore

 
#font = QtGui.QFont()
#font.setPointSize(8)
#font.setBold(False)

class RVBank_window(QtWidgets.QDialog):


    def __init__(self, parent=None):
        super(RVBank_window, self).__init__(parent=parent)
        vLayout = QtWidgets.QVBoxLayout(self)
        hLayout = QtWidgets.QHBoxLayout()


        self.lineEdit = QtWidgets.QLineEdit(self)
        hLayout.addWidget(self.lineEdit)    

        self.filter = QtWidgets.QPushButton("Search", self)
        hLayout.addWidget(self.filter)
        self.filter.clicked.connect(self.filterClicked)

        self.list = QtWidgets.QListView(self)

        vLayout.addLayout(hLayout)
        vLayout.addWidget(self.list)

        self.model = QtGui.QStandardItemModel(self.list)

        self.setGeometry(1,1, 495, 325) 

        codes = [
            'Test',
            'HD27894',
            'GJ436',
            'GJ1111',
            'BD-102030'
        ]

        for code in codes:
            item = QtGui.QStandardItem(code)
            item.setCheckable(True)
            self.model.appendRow(item)
        self.list.setModel(self.model)

        self.list.clicked.connect(self.on_clicked)



    def filterClicked(self):
        filter_text = str(self.lineEdit.text()).lower()
        for row in range(self.model.rowCount()):
            if filter_text in str(self.model.item(row).text()).lower():
                self.list.setRowHidden(row, False)
            else:
                self.list.setRowHidden(row, True)


    def on_clicked(self, index):
#        path = self.dirModel.fileInfo(index).absoluteFilePath()
#        self.listview.setRootIndex(self.fileModel.setRootPath(path))
#        return
        
        row = index.row()
        print(row,codes[row],self.model.item(row))
        #print("id = %s" %self.model.record(row).field(0).value().toString())
        #print(self.list.selectedItems())









 