import sys #,os
from PyQt5 import QtWidgets


import sys #,os
from PyQt5 import QtWidgets



class show_symbols(QtWidgets.QMainWindow):
    
    def __init__(self, parent=None):
        
        super(show_symbols, self).__init__(parent)
              
        self.title = 'This is a test window'
        self.setFixedSize(550, 800)        
        self.widget = QtWidgets.QWidget(self)
 
        self.widget.setLayout(QtWidgets.QGridLayout())
       # self.widget.layout().addWidget(self.text)

        
        #self.btn=QtGui.QPushButton("Send Signal", self)        

        self.button1 = QtWidgets.QRadioButton('this is but1', self)
        self.button1.clicked.connect(self.do_test)      
        
        self.button2 = QtWidgets.QRadioButton('this is but2', self)
        self.button2.clicked.connect(self.do_test)      
        self.setCentralWidget(self.widget)
         
    def do_test(self):        
        print("TESSSST")
        

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    main = show_symbols()
    main.show()
    sys.exit(app.exec_())        
        
    
        