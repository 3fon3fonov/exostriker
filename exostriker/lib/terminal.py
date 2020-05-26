import sys
from PyQt5 import QtCore, QtWidgets
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
        self.terminal = QtWidgets.QWidget(self)
        
        
        # Works also with urxvt:
        if subprocess.call(["which", 'urxvt'], stdout=open(os.devnull, 'wb')) == 1:
            self.process.start('xterm',['-into', str(int(self.winId()))])
        else:
            #self.process.start('xterm',['-into', str(int(self.winId()))])
            self.process.start('urxvt',['-embed', str(int(self.winId()))])
            
        #self.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)

        #self.setFixedSize(450, 340)
        self.setGeometry(1,1, 495, 390) 
        #self.terminal.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
   # def sizeHint(self):
   #     size = super(EmbTerminal, self).sizeHint()
   #     size.setHeight(size.height() + 30)
   #     size.setWidth(max(size.width(), size.height()))
   #     return size
    
    def someFunction(self):
        self.terminal.setGeometry(1, 1, self.width(), self.height())
        
        print(self.terminal.width(), self.terminal.height(),self.width(), self.height())
        #print("someFunction")
        #return super(EmbTerminal, self).__init__()
        #self.setFixedSize(self.width(), self.height())

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

 
        #self.resized.connect(self.term.someFunction)

    def resizeEvent(self, event):
        self.resized.emit()
        return super(mainWindow, self).resizeEvent(event)



    def close(self):
        self.term.process.kill()

#if __name__ == "__main__":
    
#    app = QtWidgets.QApplication(sys.argv)
#    main = mainWindow()
#    main.show()
#    sys.exit(app.exec_())
 
