import sys #,os
from PyQt6 import QtWidgets


class print_info(QtWidgets.QMainWindow):
    
    def __init__(self, parent=None):
        
        super(print_info, self).__init__(parent)
              
        self.title = 'Text Editor'
        self.setFixedSize(550, 800)        
        self.widget = QtWidgets.QWidget(self)
        #self.text = QtWidgets.QTextEdit(self.widget)
        self.text = QtWidgets.QTextBrowser(self.widget)
        self.text.setOpenExternalLinks(True)
        self.widget.setLayout(QtWidgets.QVBoxLayout())
        self.widget.layout().addWidget(self.text)
        self.setCentralWidget(self.widget)
        

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    main = print_info()
    main.show()
    sys.exit(app.exec_())        
        
        
