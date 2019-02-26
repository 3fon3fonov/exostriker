import sys
from PyQt5 import QtCore, QtGui

import logging
logger = logging.getLogger(__name__)

class QtHandler(logging.Handler):

    def __init__(self):
        logging.Handler.__init__(self)

    def emit(self, record):
        record = self.format(record)
        XStream.stdout().write("{}\n".format(record))

handler = QtHandler()
handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)

class XStream(QtCore.QObject):
    _stdout = None
    _stderr = None
    messageWritten = QtCore.pyqtSignal(str)
    def flush( self ):
        pass
    def fileno( self ):
        return -1
    def write( self, msg ):
        if ( not self.signalsBlocked() ):
            self.messageWritten.emit(msg)
    @staticmethod
    def stdout():
        if ( not XStream._stdout ):
            XStream._stdout = XStream()
            sys.stdout = XStream._stdout
        return XStream._stdout
    @staticmethod
    def stderr():
        if ( not XStream._stderr ):
            XStream._stderr = XStream()
            sys.stderr = XStream._stderr

        return XStream._stderr

class MyDialog(QtGui.QDialog):
    def __init__( self, parent = None ):
        super(MyDialog, self).__init__(parent)

        self._console = QtGui.QTextBrowser(self)
        
       # horScrollBar = self._console.horizontalScrollBar()
       # verScrollBar = self._console.verticalScrollBar()

        layout = QtGui.QVBoxLayout()
        layout.addWidget(self._console)
        #layout.addWidget(self._button)
        self.setLayout(layout)


        XStream.stdout().messageWritten.connect(self._console.insertPlainText)
        XStream.stderr().messageWritten.connect(self._console.insertPlainText)
 
      #  scrollIsAtEnd = verScrollBar.maximum() - verScrollBar.value() <= 10 
       
      #  if scrollIsAtEnd:
      #      verScrollBar.setValue(verScrollBar.maximum()) # Scrolls to the bottom
      #      horScrollBar.setValue(0) # scroll to the left

        self.pipe_output

    def pipe_output( self ):
        logger.debug('debug message')
        logger.info('info message')
        logger.warning('warning message')
        logger.error('error message')
        #print('Old school hand made print message')

if ( __name__ == '__main__' ):
    #app = None
   # if ( not QtGui.QApplication.instance() ):
    app = QtGui.QApplication([])
    dlg = MyDialog()
    dlg.show()
    #if ( app ):
    app.exec_()
