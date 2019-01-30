import sys

from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui  import *

font = QFont()
font.setPointSize(8)
font.setBold(False)
        #font.setWeight(75)
        
 

class Widget_tree(QWidget):
    def __init__(self, *args, **kwargs):
        QWidget.__init__(self, *args, **kwargs)
        hlay = QVBoxLayout(self)
        self.treeview = QTreeView() 
        self.listview = QListView()
        hlay.addWidget(self.treeview)
        hlay.addWidget(self.listview)

        path = QDir.currentPath() #QDir.rootPath()


        self.dirModel = QFileSystemModel()
        self.dirModel.setRootPath(QDir.currentPath())
        self.dirModel.setFilter(QDir.NoDotAndDotDot | QDir.AllDirs)

        self.fileModel = QFileSystemModel()
        self.fileModel.setFilter(QDir.NoDotAndDotDot |  QDir.Files)

        #filter = ['*.vels', '*.act','*.dat']
        filter = ['*.vels']
       
        
        self.fileModel.setNameFilters(filter)


        self.treeview.setModel(self.dirModel)
        self.listview.setModel(self.fileModel)
        self.treeview.setFont(font)
        self.listview.setFont(font)
       
        self.treeview.hideColumn(1)
        self.treeview.hideColumn(2)
        self.treeview.hideColumn(3)
       
        self.treeview.setRootIndex(self.dirModel.index(path))
        self.listview.setRootIndex(self.fileModel.index(path))

        self.treeview.clicked.connect(self.on_clicked)



    def on_clicked(self, index):
        path = self.dirModel.fileInfo(index).absoluteFilePath()
        self.listview.setRootIndex(self.fileModel.setRootPath(path))



if __name__ == '__main__':
    app = QApplication(sys.argv)
    w = Widget_tree()
    w.show()
    sys.exit(app.exec_())