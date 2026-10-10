import numpy as np
import sys #,os
from PyQt6 import QtCore, QtGui, QtWidgets, uic

from qtconsole.rich_jupyter_widget import RichJupyterWidget
from qtconsole.inprocess import QtInProcessKernelManager
from qtconsole.console_widget import ConsoleWidget


class ConsoleWidget_embed(RichJupyterWidget,ConsoleWidget):
    global fit

    def __init__(self, customBanner=None, *args, **kwargs):
        super(ConsoleWidget_embed, self).__init__(*args, **kwargs)
         
        if customBanner is not None:
            self.banner = customBanner

        #self.font_size = 4
        self.kernel_manager =   QtInProcessKernelManager()
        self.kernel_manager.start_kernel(show_banner=True)
        self.kernel_manager.kernel.gui = 'qt'
        self.kernel = self.kernel_manager.kernel
        self.kernel_client = self._kernel_manager.client()
        self.kernel_client.start_channels()

        def _abort_queues(kernel):
            pass
        self.kernel_manager.kernel._abort_queues = _abort_queues
        #self._execute("kernel = %s"%fit, False) 
     
        def stop():
            self.kernel_client.stop_channels()
            self.kernel_manager.shutdown_kernel()
            self.guisupport.get_app_qt().exit()

        self.exit_requested.connect(stop)


    def push_vars(self, variableDict):
        """
        Given a dictionary containing name / value pairs, push those variables
        to the Jupyter console widget
        """
        self.kernel_manager.kernel.shell.push(variableDict)

    def clear(self):
        """
        Clears the terminal
        """
        self._control.clear()

        # self.kernel_manager

#    @staticmethod
#    def get_fit(parent = None):
#        fit2 = self.kernel_manager.kernel.shell.user_ns.get('fit')
#        return (fit2)


    def print_text(self, text, before_prompt=True):
        """
        Prints some plain text to the console
        """
        self._append_plain_text(text, before_prompt=before_prompt)

    def execute_command(self, command):
        """
        Execute a command in the frame of the console widget
        """
        self._execute(command, False)


 

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    main = mainWindow()
    main.show()
    sys.exit(app.exec_())        
        
        
