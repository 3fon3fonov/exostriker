import sys

from PyQt6.QtWidgets import *
from PyQt6.QtCore import *
from PyQt6.QtGui  import *

font = QFont()
font.setPointSize(7)
font.setBold(False)
        #font.setWeight(75)
        
 

class Widget_tree(QWidget):
    
    def __init__(self, parent=None, *args, **kwargs):
        
        super(Widget_tree, self).__init__(parent)
    #def __init__(self, *args, **kwargs):
    
        self.title = 'Text Editor'
        self.setFixedSize(550, 800)        
        self.widget =  QWidget(self)
    
       # QWidget.__init__(self, *args, **kwargs)
        hlay = QVBoxLayout(self)
        
        self.setMaximumWidth(250)


        self.treeview = QTreeView() 
        self.listview = QListView()
        hlay.addWidget(self.treeview)
        hlay.addWidget(self.listview)
      #  hlay.addStretch(1)
          

        path = QDir.homePath() #QDir.currentPath() #QDir.rootPath()
#        path = ""
#        self.tree_view_tab.setRootIndex(QtGui.QFileSystemModel().index(fit.cwd))


        self.dirModel = QFileSystemModel()
        self.dirModel.setRootPath(path) #QDir.currentPath())
        self.dirModel.setFilter(QDir.NoDotAndDotDot | QDir.AllDirs)

        self.fileModel = QFileSystemModel()
        self.fileModel.setFilter(QDir.NoDotAndDotDot |  QDir.Files)

        filter = ['*.vels', '*.act','*.tran','*.dat']
        #filter = ['*.vels']
       
        
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


        self.widget.setLayout(hlay)
        self.widget.layout().addWidget(self.widget)
        self.setCentralWidget(self.widget)

    def on_clicked(self, index):
        path = self.dirModel.fileInfo(index).absoluteFilePath()
        self.listview.setRootIndex(self.fileModel.setRootPath(path))

if __name__ == "__main__":
    app = QApplication(sys.argv)
    main = Widget_tree()
    main.show()
    sys.exit(app.exec_())        

 
