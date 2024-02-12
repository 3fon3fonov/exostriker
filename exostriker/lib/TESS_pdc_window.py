import sys #,os
from PyQt6 import QtWidgets,QtGui 

 
class pdc(QtWidgets.QDialog):


    def __init__(self, parent = None):
       # super(show_symbols, self).__init__(parent)
        super(pdc, self).__init__()

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
        self.button1  = QtWidgets.QRadioButton('SAP FLUX', self)
        self.button2  = QtWidgets.QRadioButton('PDCSAP FLUX', self)
 


        self.radio_group.addButton(self.button1)
        self.radio_group.addButton(self.button2)
 

        self.radio_group.setId(self.button1,1)
        self.radio_group.setId(self.button2,2)
 

        self.layout.addWidget(self.button1)
        self.layout.addWidget(self.button2)
 
        self.cancel_button = QtWidgets.QPushButton('Accept', self)
        self.layout.addWidget(self.cancel_button)

        self.button1.setChecked(True)
        self.cancel_button.clicked.connect(self.close)
        #self.Ok_button.clicked.connect(self.get_radio)



    def return_but_N(self):
        #   Return list of values. It need map with str (self.lineedit.text() will return QString)
        return self.radio_group.checkedId() 

    # static method to create the dialog and return button pressed
    @staticmethod
    def get_radio(parent = None):
        dialog = pdc(parent)
        result = dialog.exec_()
        rad_but = dialog.return_but_N()
        return (rad_but)



if __name__ == '__main__':
   # app = QtWidgets.QApplication(sys.argv)
    #w = show_symbols()
   # w.show()
    #sys.exit(app.exec_())    
    app = QtWidgets.QApplication([])
    but = pdc.DateDialog.get_radio()
#    print("{} {} {}".format(date, time, ok))
    app.exec_()        
