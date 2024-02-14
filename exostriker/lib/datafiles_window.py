import sys #,os
from PyQt6 import QtWidgets,QtGui,QtCore


font = QtGui.QFont()
font.setPointSize(8)
font.setBold(False)

class datafiles_window(QtWidgets.QDialog):


    def __init__(self, parent = None):
       # super(show_symbols, self).__init__(parent)
        super(datafiles_window, self).__init__()

        self.layout = QtWidgets.QVBoxLayout(self)
        self.title = 'Select valid data file'
       # self.setFixedSize(550, 800)        
        self.widget=QtWidgets.QWidget(self)  # central widget
        self.setGeometry(1,1, 495, 325) 

        
 
        self.treeview = QtWidgets.QTreeView() 
        self.listview = QtWidgets.QListView()
        self.layout.addWidget(self.treeview)
        self.layout.addWidget(self.listview)
        
        
        path = QtCore.QDir.homePath()
        
        self.dirModel = QtGui.QFileSystemModel()
        self.dirModel.setRootPath(path) #QDir.currentPath())
        self.dirModel.setFilter(QtCore.QDir.NoDotAndDotDot | QtCore.QDir.AllDirs)

        self.fileModel = QtGui.QFileSystemModel()
        self.fileModel.setFilter(QtCore.QDir.NoDotAndDotDot |  QtCore.QDir.Files)

        filter = ['*.vels', '*.act','*.tran','*.dat']

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


 
 


#        self.cancel_button = QtGui.QPushButton('Close', self)
#        self.layout.addWidget(self.cancel_button)

#        self.cancel_button.clicked.connect(self.close)
        #self.Ok_button.clicked.connect(self.get_radio)


    def on_clicked(self, index):
        path = self.dirModel.fileInfo(index).absoluteFilePath()
        self.listview.setRootIndex(self.fileModel.setRootPath(path))
        return
 
 

if __name__ == '__main__':
   # app = QtWidgets.QApplication(sys.argv)
    #w = show_symbols()
   # w.show()
    #sys.exit(app.exec_())    
    app = QtWidgets.QApplication([])
    app.exec_()        
