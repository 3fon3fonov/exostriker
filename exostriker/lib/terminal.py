import sys
from PyQt6 import QtCore, QtWidgets
import subprocess
import os


class EmbTerminal(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super(EmbTerminal, self).__init__(parent)

        layout = QtWidgets.QGridLayout()
       # layout.setSpacing(0)  # for simplicity
       # layout.setSizeConstraint(QtWidgets.QLayout.SetNoConstraint)

        self.setLayout(layout)
        
        self.process = QtCore.QProcess(self)
        #self.terminal = QtWidgets.QWidget(self)
        
        
        # Works also with urxvt:
        if subprocess.call(["which", 'urxvt2'], stdout=open(os.devnull, 'wb')) == 1:
            self.process.start('xterm',['-into', str(int(self.winId()))])
        else:
            self.process.start('urxvt',['-embed', str(int(self.winId()))])
            
        #self.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)

        self.setGeometry(1,1, 490, 310) 
        #self.setLayout(layout)
        #self.setFixedSize(395, 290)
 

    def close(self):
        self.process.kill()
        
class mainWindow(QtWidgets.QMainWindow):
    
    resized = QtCore.pyqtSignal()
    
    def __init__(self, parent=None):
        super(mainWindow, self).__init__(parent)

        layout = QtWidgets.QGridLayout()
        self.term = EmbTerminal()
        layout.addWidget(self.term)

        container = QtWidgets.QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)
        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)

        #self.term.setFixedSize(395, 290)

        self.resized.connect(self.someFunction)

    def someFunction(self):
        self.term.setGeometry(1, 1, self.width(), self.height())
        
 
    def resizeEvent(self, event):
        self.resized.emit()
        return super(mainWindow, self).resizeEvent(event)



    def close(self):
        self.term.process.kill()



if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    main = mainWindow()
    main.show()
    sys.exit(app.exec_())

 
