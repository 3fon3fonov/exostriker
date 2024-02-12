import sys, traceback
from PyQt6 import QtCore, QtGui, QtWidgets

#sys.stdout.isatty = lambda: False


import logging
logger = logging.getLogger(__name__)

#from dynesty import __version__ as dynesty_version

#sys.stdout = sys.__stdout__
#sys.stderr = sys.__stderr__

sys._excepthook = sys.excepthook 
def exception_hook(exctype, value, traceback):
    print(exctype, value, traceback)
    sys._excepthook(exctype, value, traceback) 
    #sys.exit(1) 
sys.excepthook = exception_hook 


#def excepthook(exc_type, exc_value, exc_tb):
#    tb = "".join(traceback.format_exception(exc_type, exc_value, exc_tb))
#    print("error catched!:")
#    print("error message:\n", tb)
#    QtWidgets.QApplication.quit()
    # or QtWidgets.QApplication.exit(0)

#sys.excepthook = excepthook


class QtHandler(logging.Handler):

    def __init__(self):
        logging.Handler.__init__(self)

    def emit(self, record):
        record = self.format(record)
        XStream.stdout().write("{}\n".format(record))
        #XStream.stderr().write("{}\n".format(record))

handler = QtHandler()
handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)

#sys.stdout.isatty = lambda: False
#sys.stdout.encoding = sys.getdefaultencoding()

def ensureUtf(string, encoding='utf8'):
  """Converts input to unicode if necessary."""
  if type(string) == bytes:
    return string.decode(encoding, 'ignore')
  else:
    return string


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
            self.messageWritten.emit(ensureUtf(msg))
    @staticmethod
    def stdout():
        if ( not XStream._stdout ):
            XStream._stdout = XStream()
            #XStream._stdout.isatty = lambda: False
            sys.stdout = XStream._stdout
            sys.stdout.isatty = lambda: False  
            sys.stdout.encoding = sys.getdefaultencoding()
            
        return XStream._stdout
    @staticmethod
    def stderr():
        if ( not XStream._stderr ):
            XStream._stderr = XStream()
            sys.stderr = XStream._stderr
            sys.stderr.isatty = lambda: False
            sys.stderr.encoding = sys.getdefaultencoding()
 
        return XStream._stderr

class LogMessageViewer(QtWidgets.QTextBrowser):

    def __init__(self, parent=None):
        super(LogMessageViewer,self).__init__(parent)
        self.setReadOnly(True)
        #self.setLineWrapMode(QtGui.QTextEdit.NoWrap)


    @QtCore.pyqtSlot(str)
    def appendLogMessage(self, msg):
        #print(msg)
        horScrollBar = self.horizontalScrollBar()
        verScrollBar = self.verticalScrollBar()
        scrollIsAtEnd = verScrollBar.maximum() - verScrollBar.value() <= 10

        self.insertPlainText(ensureUtf(msg))

        if scrollIsAtEnd:
            verScrollBar.setValue(verScrollBar.maximum()) # Scrolls to the bottom
            horScrollBar.setValue(0) # scroll to the left

class MyDialog(QtWidgets.QDialog):
    def __init__( self, parent = None ):
        super(MyDialog, self).__init__(parent)

        self._console = LogMessageViewer(self)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self._console)
        self.setLayout(layout)

        XStream.stdout().messageWritten.connect(self._console.appendLogMessage)
        #if sys.version_info[0] == 2:
        #    XStream.stderr().messageWritten.connect(self._console.appendLogMessage)
        #elif  int(dynesty_version[0]) < 1 and int(dynesty_version[4]) <= 2:   # ignors stderr pipe if dynesty>0.9.2. TBFixed !!!!! 
        XStream.stderr().messageWritten.connect(self._console.appendLogMessage)

    def pipe_output( self ):
        logger.debug('debug message')
        logger.info('info message')
        logger.warning('warning message')
        logger.error('error message')


class DebugDialog(QtWidgets.QDialog):
    def __init__( self, parent = None ):
        super(DebugDialog, self).__init__(parent)

        self._console = LogMessageViewer(self)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self._console)
        self.setLayout(layout)

        self._console.appendLogMessage("The Exo-Striker is in a '-debug' mode. All stdout/stderr output is redirected to the main Terminal.")
 

    def pipe_output( self ):
        logger.debug('debug message')
        logger.info('info message')
        logger.warning('warning message')
        logger.error('error message')
        #print('Old school hand made print message')

#if ( __name__ == '__main__' ):
    
    

    #app = None
    # if ( not QtWidgets.QApplication.instance() ):
#    app = QtWidgets.QApplication(sys.argv)
#    dlg = MyDialog()
#    dlg.show()

    #if ( app ):
    #app.exec_()
#    sys.exit(app.exec_())

 
    
    
    
    
    
    
    
    
    
    
