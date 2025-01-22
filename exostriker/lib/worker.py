import sys, traceback 
from PyQt6 import QtCore #, QtWidgets
#import subprocess
#import os

import threading
import sys, traceback 
  
class WorkerSignals(QtCore.QObject):
    '''
    Defines the signals available from a running worker thread.
    Supported signals are:

    finished
        No data

    error
        `tuple` (exctype, value, traceback.format_exc() )

    result
        `object` data returned from processing, anything

    progress
        `int` indicating % progress
    '''
    finished = QtCore.pyqtSignal()
    error = QtCore.pyqtSignal(tuple)
    result = QtCore.pyqtSignal(object)
    progress = QtCore.pyqtSignal(int)


class Worker(QtCore.QRunnable):
    '''
    Worker thread

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.

    :param callback: The function callback to run on this worker thread. Supplied args and 
                     kwargs will be passed through to the runner.
    :type callback: function
    :param args: Arguments to pass to the callback function
    :param kwargs: Keywords to pass to the callback function

    '''

    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()

        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

        # Add the callback to our kwargs
        self.kwargs['progress_callback'] = self.signals.progress

        # Add a threading.Event to signal stopping
        self.stop_event = threading.Event()

    @QtCore.pyqtSlot()
    def run(self):
        '''
        Initialise the runner function with passed args, kwargs.
        '''
        try:
            # Ensure the function is aware of the stop_event
            if 'stop_event' in self.kwargs:
                self.kwargs['stop_event'] = self.stop_event
            result = self.fn()
            #result = self.fn(*self.args, **self.kwargs)
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.signals.result.emit(result)  # Return the result of the processing
        finally:
            self.signals.finished.emit()  # Done

    def stop(self):
        '''
        Signal the worker to stop.
        '''
        print('TEST!!!')
        self.stop_event.set()
        #self.signals.finished.emit()  # Done        
        
        
        
        
        
        
        
