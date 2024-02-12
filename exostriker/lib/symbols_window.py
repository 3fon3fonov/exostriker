import sys #,os
from PyQt6 import QtWidgets,QtGui 

 
class show_symbols(QtWidgets.QDialog):


    def __init__(self, parent = None):
       # super(show_symbols, self).__init__(parent)
        super(show_symbols, self).__init__()

        self.layout = QtWidgets.QVBoxLayout(self)
        
 
        self.title = 'This is a test window'
       # self.setFixedSize(550, 800)        
       # self.widget = QtWidgets.QWidget(self)
       # self.widget.setLayout(QtWidgets.QGridLayout())
       # self.widget.layout().addWidget(self.text)

       #self.layout=QtWidgets.QGridLayout() # layout for the central widget
        self.widget=QtWidgets.QWidget(self)  # central widget
      #  self.widget.setLayout(layout)    

        self.init_buttons()     

    
    def init_buttons(self):


        self.radio_group=QtWidgets.QButtonGroup(self.widget) # Number group
        self.button1  = QtWidgets.QRadioButton('symbol "o"', self)
        self.button2  = QtWidgets.QRadioButton('symbol "t"', self)
        self.button3  = QtWidgets.QRadioButton('symbol "t1"', self)
        self.button4  = QtWidgets.QRadioButton('symbol "t2"', self)
        self.button5  = QtWidgets.QRadioButton('symbol "t3"', self)
        self.button6  = QtWidgets.QRadioButton('symbol "s"', self)
        self.button7  = QtWidgets.QRadioButton('symbol "p"', self)
        self.button8  = QtWidgets.QRadioButton('symbol "h"', self)
        self.button9  = QtWidgets.QRadioButton('symbol "star"', self)
        self.button10 = QtWidgets.QRadioButton('symbol "+"', self)
        self.button11 = QtWidgets.QRadioButton('symbol "d"', self)


        self.radio_group.addButton(self.button1)
        self.radio_group.addButton(self.button2)
        self.radio_group.addButton(self.button3)
        self.radio_group.addButton(self.button4)
        self.radio_group.addButton(self.button5)
        self.radio_group.addButton(self.button6)
        self.radio_group.addButton(self.button7)
        self.radio_group.addButton(self.button8)
        self.radio_group.addButton(self.button9)
        self.radio_group.addButton(self.button10)
        self.radio_group.addButton(self.button11)

        self.radio_group.setId(self.button1,1)
        self.radio_group.setId(self.button2,2)
        self.radio_group.setId(self.button3,3)
        self.radio_group.setId(self.button4,4)
        self.radio_group.setId(self.button5,5)
        self.radio_group.setId(self.button6,6)
        self.radio_group.setId(self.button7,7)
        self.radio_group.setId(self.button8,8)
        self.radio_group.setId(self.button9,9)
        self.radio_group.setId(self.button10,10)
        self.radio_group.setId(self.button11,11)

        self.layout.addWidget(self.button1)
        self.layout.addWidget(self.button2)
        self.layout.addWidget(self.button3)
        self.layout.addWidget(self.button4)
        self.layout.addWidget(self.button5)
        self.layout.addWidget(self.button6)
        self.layout.addWidget(self.button7)
        self.layout.addWidget(self.button8)
        self.layout.addWidget(self.button9)
        self.layout.addWidget(self.button10)
        self.layout.addWidget(self.button11)


        #self.setCentralWidget(self.widget)

        #self.Ok_button = QtGui.QPushButton('OK', self)
        #self.layout.addWidget(self.Ok_button)        
        
        self.cancel_button = QtWidgets.QPushButton('Accept', self)
        self.layout.addWidget(self.cancel_button)

        self.cancel_button.clicked.connect(self.close)
        #self.Ok_button.clicked.connect(self.get_radio)



    def return_but_N(self):
        #   Return list of values. It need map with str (self.lineedit.text() will return QString)
        return self.radio_group.checkedId() 

    # static method to create the dialog and return button pressed
    @staticmethod
    def get_radio(parent = None):
        dialog = show_symbols(parent)
        result = dialog.exec_()
        rad_but = dialog.return_but_N()
        return (rad_but)



if __name__ == '__main__':
   # app = QtWidgets.QApplication(sys.argv)
    #w = show_symbols()
   # w.show()
    #sys.exit(app.exec_())    
    app = QtWidgets.QApplication([])
    but = show_symbols.DateDialog.get_radio()
#    print("{} {} {}".format(date, time, ok))
    app.exec_()        
