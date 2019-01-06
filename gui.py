#!/usr/bin/python
__author__ = 'Trifon Trifonov'


import numpy as np
#import matplotlib as mpl
#mpl.use('Qt5Agg')
import sys,  traceback 
from PyQt5 import QtCore, QtGui, QtWidgets, uic

sys.path.insert(0, './addons')

import RV_mod as rv

import pyqtgraph as pg
import pyqtgraph.console as pg_console

import wordprocessor as ted
import calculator as calc 
import stdout_pipe as stdout_pipe
import gls as gls 

import time

#import BKR as bkr
from doublespinbox import DoubleSpinBox
from Jupyter_emb import ConsoleWidget_embed

from scipy.signal import argrelextrema

import batman as batman

import webbrowser
 
#try:
#    import cPickle as pickle
#except ModuleNotFoundError:
#    import pickle

import dill
import terminal

qtCreatorFile = "rvmod_gui.ui" # Enter file here.
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

pg.setConfigOption('background', '#ffffff')
pg.setConfigOption('foreground', 'k')
pg.setConfigOptions(antialias=True)  
 
 

global fit, colors,ses_list
 

fit=rv.signal_fit(name='init')
ses_list = []

#'#cc0000',

colors = ['#0066ff',  '#ff0000','#66ff66','#00ffff','#cc33ff','#ff9900','#cccc00','#3399ff','#990033','#339933','#666699']

QtGui.QApplication.processEvents()

 

class print_info(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        super(print_info, self).__init__(parent)
        
        
        self.title = 'Text Editor'
        #self.left = 10
        #self.top = 10
        #self.width = 1080
        #self.height = 920
        #self.setGeometry(3,30,450,800)
        self.setFixedSize(550, 800)
        

        self.widget = QtWidgets.QWidget(self)

        #self.text = QtWidgets.QTextEdit(self.widget)
        self.text = QtWidgets.QTextBrowser(self.widget)
        self.text.setOpenExternalLinks(True)
       
        
        self.widget.setLayout(QtWidgets.QVBoxLayout())
        self.widget.layout().addWidget(self.text)
        self.setCentralWidget(self.widget)
        
     
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

    @QtCore.pyqtSlot()
    def run(self):
        '''
        Initialise the runner function with passed args, kwargs.
        '''

        # Retrieve args/kwargs here; and fire processing using them
        try:
            result = self.fn()
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.signals.result.emit(result)  # Return the result of the processing
        finally:
            self.signals.finished.emit()  # Done

class MyApp(QtWidgets.QMainWindow, Ui_MainWindow):

    def update_labels(self):
        global fit

        self.value_stellar_mass.setText("%.2f"%(fit.params.stellar_mass))
        self.value_epoch.setText(str(fit.epoch))
        self.value_rms.setText("%.4f"%(fit.fit_results.rms))
        self.value_chi2.setText("%.4f"%(fit.fit_results.chi2)) 
        self.value_reduced_chi2.setText("%.4f"%(fit.fit_results.reduced_chi2))        
        self.value_loglik.setText("%.4f"%(fit.fit_results.loglik)) 
        self.value_loglik.setText("%.4f"%(fit.fit_results.loglik)) 
        self.value_Ndata.setText("%s"%(len(fit.fit_results.jd))) 
        self.value_DOF.setText("%s"%(len(fit.fit_results.jd) - fit.fit_results.mfit))                 

    def update_gui_params(self):
        global fit

        param_gui = [self.K1, self.P1, self.e1, self.om1, self.ma1, self.incl1, self.Omega1,
                     self.K2, self.P2, self.e2, self.om2, self.ma2, self.incl2, self.Omega2,
                     self.K3, self.P3, self.e3, self.om3, self.ma3, self.incl3, self.Omega3,
                     self.K4, self.P4, self.e4, self.om4, self.ma4, self.incl4, self.Omega4, 
                     self.K5, self.P5, self.e5, self.om5, self.ma5, self.incl5, self.Omega5,
                     self.K6, self.P6, self.e6, self.om6, self.ma6, self.incl6, self.Omega6,
                     self.K7, self.P7, self.e7, self.om7, self.ma7, self.incl7, self.Omega7, 
                     self.K8, self.P8, self.e8, self.om8, self.ma8, self.incl8, self.Omega8,
                     self.K9, self.P9, self.e9, self.om9, self.ma9, self.incl9, self.Omega9,
                     ]
         
        for i in range(fit.npl*7):
            param_gui[i].setValue(fit.params.planet_params[i]) 
            

        data_gui = [self.Data1,self.Data2,self.Data3,self.Data4,self.Data5,
                    self.Data6,self.Data7,self.Data8,self.Data9,self.Data10]
        data_jitter_gui = [self.jitter_Data1,self.jitter_Data2,self.jitter_Data3,self.jitter_Data4,self.jitter_Data5,
                           self.jitter_Data6,self.jitter_Data7,self.jitter_Data8,self.jitter_Data9,self.jitter_Data10]

        for i in range(10): 
            data_gui[i].setValue(fit.params.offsets[i]) 
            data_jitter_gui[i].setValue(fit.params.jitters[i])
            
        self.St_mass_input.setValue(fit.params.stellar_mass)
        
        self.RV_lin_trend.setValue(fit.params.linear_trend)
        
        self.Epoch.setValue(fit.epoch)


    def update_params(self):
        global fit

        param_gui = [self.K1, self.P1, self.e1, self.om1, self.ma1, self.incl1, self.Omega1,
                     self.K2, self.P2, self.e2, self.om2, self.ma2, self.incl2, self.Omega2,
                     self.K3, self.P3, self.e3, self.om3, self.ma3, self.incl3, self.Omega3,
                     self.K4, self.P4, self.e4, self.om4, self.ma4, self.incl4, self.Omega4, 
                     self.K5, self.P5, self.e5, self.om5, self.ma5, self.incl5, self.Omega5,
                     self.K6, self.P6, self.e6, self.om6, self.ma6, self.incl6, self.Omega6,
                     self.K7, self.P7, self.e7, self.om7, self.ma7, self.incl7, self.Omega7, 
                     self.K8, self.P8, self.e8, self.om8, self.ma8, self.incl8, self.Omega8,
                     self.K9, self.P9, self.e9, self.om9, self.ma9, self.incl9, self.Omega9,
                     ]
                     

        for i in range(fit.npl*7):
            fit.params.planet_params[i] = param_gui[i].value() 

        data_gui = [self.Data1,self.Data2,self.Data3,self.Data4,self.Data5,
                    self.Data6,self.Data7,self.Data8,self.Data9,self.Data10]
        data_jitter_gui = [self.jitter_Data1,self.jitter_Data2,self.jitter_Data3,self.jitter_Data4,self.jitter_Data5,
                           self.jitter_Data6,self.jitter_Data7,self.jitter_Data8,self.jitter_Data9,self.jitter_Data10]

        for i in range(10): 
            fit.params.offsets[i] = data_gui[i].value() 
            fit.params.jitters[i] = data_jitter_gui[i].value()

        fit.params.stellar_mass = self.St_mass_input.value() 
 
        fit.params.linear_trend = self.RV_lin_trend.value()
        
        fit.epoch =  self.Epoch.value()
       


    def update_errors(self):
        global fit

        param_errors_gui = [self.err_K1,self.err_P1,self.err_e1,self.err_om1,self.err_ma1, self.err_i1, self.err_Om1,
                            self.err_K2,self.err_P2,self.err_e2,self.err_om2,self.err_ma2, self.err_i2, self.err_Om2,
                            self.err_K3,self.err_P3,self.err_e3,self.err_om3,self.err_ma3, self.err_i3, self.err_Om3,
                            self.err_K4,self.err_P4,self.err_e4,self.err_om4,self.err_ma4, self.err_i4, self.err_Om4,  
                            self.err_K5,self.err_P5,self.err_e5,self.err_om5,self.err_ma5, self.err_i5, self.err_Om5,
                            self.err_K6,self.err_P6,self.err_e6,self.err_om6,self.err_ma6, self.err_i6, self.err_Om6,
                            self.err_K7,self.err_P7,self.err_e7,self.err_om7,self.err_ma7, self.err_i7, self.err_Om7, 
                            self.err_K8,self.err_P8,self.err_e8,self.err_om8,self.err_ma8, self.err_i8, self.err_Om8,
                            self.err_K9,self.err_P9,self.err_e9,self.err_om9,self.err_ma9, self.err_i9, self.err_Om9,                       
                            ]
        for i in range(fit.npl*7):
            param_errors_gui[i].setText("+/- %.3f"%max(np.abs(fit.param_errors.planet_params_errors[i])))


        data_errors_gui        = [self.err_Data1,self.err_Data2,self.err_Data3,self.err_Data4,self.err_Data5,
                                  self.err_Data6,self.err_Data7,self.err_Data8,self.err_Data9,self.err_Data10]
        data_errors_jitter_gui = [self.err_jitter_Data1,self.err_jitter_Data2,self.err_jitter_Data3,self.err_jitter_Data4,
                                  self.err_jitter_Data5,self.err_jitter_Data6,self.err_jitter_Data7,self.err_jitter_Data8,
                                  self.err_jitter_Data9,self.err_jitter_Data10]

        for i in range(10):
            data_errors_gui[i].setText("+/- %.3f"%max(np.abs(fit.param_errors.offset_errors[i])))
            data_errors_jitter_gui[i].setText("+/- %.3f"%max(np.abs(fit.param_errors.jitter_errors[i])))


        self.err_RV_lin_trend.setText("+/- %.8f"%(max(fit.param_errors.linear_trend_error)))




    def update_a_mass(self):
        global fit

        param_a_gui = [self.label_a1, self.label_a2, self.label_a3, self.label_a4, self.label_a5, 
                       self.label_a6, self.label_a7, self.label_a8, self.label_a9]
        param_mass_gui = [self.label_mass1, self.label_mass2, self.label_mass3, self.label_mass4, self.label_mass5, 
                       self.label_mass6, self.label_mass7, self.label_mass8, self.label_mass9]
        param_t_peri_gui = [self.label_t_peri1, self.label_t_peri2, self.label_t_peri3, self.label_t_peri4, self.label_t_peri5, 
                       self.label_t_peri6, self.label_t_peri7, self.label_t_peri8, self.label_t_peri9]


        for i in range(fit.npl):
            param_a_gui[i].setText("%.3f"%(fit.fit_results.a[i])) 
            param_mass_gui[i].setText("%.3f"%(fit.fit_results.mass[i])) 
            param_t_peri_gui[i].setText("%.3f"%( (float(fit.epoch) - (np.radians(fit.params.planet_params[7*i + 4])/(2*np.pi))*fit.params.planet_params[7*i + 1] ))) # epoch  - ((ma/TWOPI)*a[1])




    def update_use_from_input_file(self):
        global fit


        use_param_gui =  [self.use_K1, self.use_P1, self.use_e1, self.use_om1, self.use_ma1, self.use_incl1, self.use_Omega1,
                          self.use_K2, self.use_P2, self.use_e2, self.use_om2, self.use_ma2, self.use_incl2, self.use_Omega2,
                          self.use_K3, self.use_P3, self.use_e3, self.use_om3, self.use_ma3, self.use_incl3, self.use_Omega3,                        
                          self.use_K4, self.use_P4, self.use_e4, self.use_om4, self.use_ma4, self.use_incl4, self.use_Omega4,    
                          self.use_K5, self.use_P5, self.use_e5, self.use_om5, self.use_ma5, self.use_incl5, self.use_Omega5,    
                          self.use_K6, self.use_P6, self.use_e6, self.use_om6, self.use_ma6, self.use_incl6, self.use_Omega6, 
                          self.use_K7, self.use_P7, self.use_e7, self.use_om7, self.use_ma7, self.use_incl7, self.use_Omega7,    
                          self.use_K8, self.use_P8, self.use_e8, self.use_om8, self.use_ma8, self.use_incl8, self.use_Omega8,    
                          self.use_K9, self.use_P9, self.use_e9, self.use_om9, self.use_ma9, self.use_incl9, self.use_Omega9,                       
                          ]
        
        for i in range(fit.npl*7):
            use_param_gui[i].setChecked(bool(fit.use.use_planet_params[i]))
#            print(fit.use.use_planet_params[i])

        #use_data_gui = [self.use_Data1,self.use_Data2,self.use_Data3,self.use_Data4,self.use_Data5,
        #            self.use_Data6,self.use_Data7,self.use_Data8,self.use_Data9,self.use_Data10]

        use_data_offset_gui = [self.use_offset_Data1,self.use_offset_Data2,self.use_offset_Data3,self.use_offset_Data4,
                               self.use_offset_Data5,self.use_offset_Data6,self.use_offset_Data7,self.use_offset_Data8,
                               self.use_offset_Data9,self.use_offset_Data10]
        use_data_jitter_gui = [self.use_jitter_Data1,self.use_jitter_Data2,self.use_jitter_Data3,self.use_jitter_Data4,self.use_jitter_Data5,
                               self.use_jitter_Data6,self.use_jitter_Data7,self.use_jitter_Data8,self.use_jitter_Data9,self.use_jitter_Data10]

       # print(fit.filelist.ndset)

        for i in range(10): 
            #use_data_gui[i].setChecked(bool(fit.use.use_offsets[i])) # attention, TBF
            use_data_jitter_gui[i].setChecked(bool(fit.use.use_jitters[i]))
            use_data_offset_gui[i].setChecked(bool(fit.use.use_offsets[i])) 

        planet_checked_gui = [self.use_Planet1,self.use_Planet2,self.use_Planet3,self.use_Planet4,self.use_Planet5,
                              self.use_Planet6,self.use_Planet7,self.use_Planet8,self.use_Planet9]
        for i in range(9):  
            if i < fit.npl:
                planet_checked_gui[i].setChecked(True)  
            else:
                planet_checked_gui[i].setChecked(False)  

            
        self.use_RV_lin_trend.setChecked(bool(fit.use.use_linear_trend)) 
            

    def update_use(self):
        global fit
        
        use_planet_gui = [self.use_Planet1,self.use_Planet2,self.use_Planet3,self.use_Planet4,self.use_Planet5,
                          self.use_Planet6,self.use_Planet7,self.use_Planet8,self.use_Planet9]
        #for i in range(len(use_planet_gui)):  
        npl_old = fit.npl
        checked = int(np.sum( [use_planet_gui[i].isChecked() for i in range(len(use_planet_gui))] ))
         
         
        #print fit.params.planet_params
         
        #print     fit.npl, nchecked
        if npl_old < checked:
            fit.add_planet()
        elif npl_old >= checked:
            fit.npl = checked
 
        

        use_param_gui2 = [self.use_K1, self.use_P1, self.use_e1, self.use_om1, self.use_ma1, self.use_incl1, self.use_Omega1,
                          self.use_K2, self.use_P2, self.use_e2, self.use_om2, self.use_ma2, self.use_incl2, self.use_Omega2,
                          self.use_K3, self.use_P3, self.use_e3, self.use_om3, self.use_ma3, self.use_incl3, self.use_Omega3,                        
                          self.use_K4, self.use_P4, self.use_e4, self.use_om4, self.use_ma4, self.use_incl4, self.use_Omega4,    
                          self.use_K5, self.use_P5, self.use_e5, self.use_om5, self.use_ma5, self.use_incl5, self.use_Omega5,    
                          self.use_K6, self.use_P6, self.use_e6, self.use_om6, self.use_ma6, self.use_incl6, self.use_Omega6, 
                          self.use_K7, self.use_P7, self.use_e7, self.use_om7, self.use_ma7, self.use_incl7, self.use_Omega7,    
                          self.use_K8, self.use_P8, self.use_e8, self.use_om8, self.use_ma8, self.use_incl8, self.use_Omega8,    
                          self.use_K9, self.use_P9, self.use_e9, self.use_om9, self.use_ma9, self.use_incl9, self.use_Omega9,                       
                          ]

        for i in range(fit.npl*7):
            fit.use.use_planet_params[i] = int(use_param_gui2[i].isChecked())

        use_data_offset_gui = [self.use_offset_Data1,self.use_offset_Data2,self.use_offset_Data3,self.use_offset_Data4,
                               self.use_offset_Data5,self.use_offset_Data6,self.use_offset_Data7,self.use_offset_Data8,
                               self.use_offset_Data9,self.use_offset_Data10]
        use_data_jitter_gui = [self.use_jitter_Data1,self.use_jitter_Data2,self.use_jitter_Data3,self.use_jitter_Data4,self.use_jitter_Data5,
                               self.use_jitter_Data6,self.use_jitter_Data7,self.use_jitter_Data8,self.use_jitter_Data9,self.use_jitter_Data10]

        for i in range(10): 
            fit.use.use_jitters[i] = int(use_data_jitter_gui[i].isChecked())
            fit.use.use_offsets[i] = int(use_data_offset_gui[i].isChecked())   
           # print("test")         
           # print(fit.use.use_jitters[i]) 
           # print(fit.use.use_offsets[i]) 
        fit.use.use_linear_trend = int(self.use_RV_lin_trend.isChecked()) 

 
####################################################        
  
    def initialize_buttons(self):

        # for some reason this does not work!
        #[self.buttonGroup_4.setId(bg4, ii) for ii, bg4 in enumerate(self.buttonGroup_4.buttons())]
        #[self.buttonGroup_remove_RV_data.setId(bg5, jj) for jj, bg5 in enumerate(self.buttonGroup_remove_RV_data.buttons())]   
        
        self.buttonGroup_4.setId(self.Button_RV_data_1,1)
        self.buttonGroup_4.setId(self.Button_RV_data_2,2)
        self.buttonGroup_4.setId(self.Button_RV_data_3,3)
        self.buttonGroup_4.setId(self.Button_RV_data_4,4)
        self.buttonGroup_4.setId(self.Button_RV_data_5,5)
        self.buttonGroup_4.setId(self.Button_RV_data_6,6)
        self.buttonGroup_4.setId(self.Button_RV_data_7,7)
        self.buttonGroup_4.setId(self.Button_RV_data_8,8)
        self.buttonGroup_4.setId(self.Button_RV_data_9,9)
        self.buttonGroup_4.setId(self.Button_RV_data_10,10)

        self.buttonGroup_remove_RV_data.setId(self.remove_rv_data1,1)
        self.buttonGroup_remove_RV_data.setId(self.remove_rv_data2,2)
        self.buttonGroup_remove_RV_data.setId(self.remove_rv_data3,3)
        self.buttonGroup_remove_RV_data.setId(self.remove_rv_data4,4)
        self.buttonGroup_remove_RV_data.setId(self.remove_rv_data5,5)
        self.buttonGroup_remove_RV_data.setId(self.remove_rv_data6,6)
        self.buttonGroup_remove_RV_data.setId(self.remove_rv_data7,7)
        self.buttonGroup_remove_RV_data.setId(self.remove_rv_data8,8)
        self.buttonGroup_remove_RV_data.setId(self.remove_rv_data9,9)
        self.buttonGroup_remove_RV_data.setId(self.remove_rv_data10,10)
        
        
        self.buttonGroup_transit_data.setId(self.Button_transit_data_1,1)
        self.buttonGroup_transit_data.setId(self.Button_transit_data_2,2)
        self.buttonGroup_transit_data.setId(self.Button_transit_data_3,3)
        self.buttonGroup_transit_data.setId(self.Button_transit_data_4,4)
        self.buttonGroup_transit_data.setId(self.Button_transit_data_5,5)
        self.buttonGroup_transit_data.setId(self.Button_transit_data_6,6)
        self.buttonGroup_transit_data.setId(self.Button_transit_data_7,7)
        self.buttonGroup_transit_data.setId(self.Button_transit_data_8,8)
        self.buttonGroup_transit_data.setId(self.Button_transit_data_9,9)
        self.buttonGroup_transit_data.setId(self.Button_transit_data_10,10)

        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data1,1)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data2,2)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data3,3)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data4,4)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data5,5)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data6,6)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data7,7)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data8,8)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data9,9)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data10,10)

        self.buttonGroup_activity_data.setId(self.Button_activity_data_1,1)
        self.buttonGroup_activity_data.setId(self.Button_activity_data_2,2)
        self.buttonGroup_activity_data.setId(self.Button_activity_data_3,3)
        self.buttonGroup_activity_data.setId(self.Button_activity_data_4,4)
        self.buttonGroup_activity_data.setId(self.Button_activity_data_5,5)
        self.buttonGroup_activity_data.setId(self.Button_activity_data_6,6)
        self.buttonGroup_activity_data.setId(self.Button_activity_data_7,7)
        self.buttonGroup_activity_data.setId(self.Button_activity_data_8,8)
        self.buttonGroup_activity_data.setId(self.Button_activity_data_9,9)
        self.buttonGroup_activity_data.setId(self.Button_activity_data_10,10)

        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data1,1)
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data2,2)
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data3,3)
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data4,4)
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data5,5)
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data6,6)
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data7,7)
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data8,8)
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data9,9)
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data10,10)
      
        
    def initialize_plots(self):

        global p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,pe

        p1  = self.graphicsView_timeseries_RV
        p2  = self.graphicsView_timeseries_RV_o_c
        p3  = self.graphicsView_timeseries_phot
        p4  = self.graphicsView_timeseries_phot_o_c
        p5  = self.graphicsView_timeseries_activity
        p6  = self.graphicsView_timeseries_correlations
                
        p7  = self.graphicsView_peridogram_RV 
        p8  = self.graphicsView_periodogram_RV_o_c  
        p9  = self.graphicsView_peridogram_phot
        p10 = self.graphicsView_peridogram_phot_o_c        
        p11 = self.graphicsView_periodogram_activity
        p12 = self.graphicsView_periodogram_window  
        
        p13 = self.graphicsView_orb_evol_elements_a
        p14 = self.graphicsView_orb_evol_elements_e        
        p15 = self.graphicsView_orbital_view
        
        pe  = self.graphicsView_extra_plot

        xaxis = ['BJD','BJD','BJD','BJD','BJD','','days','days','days','days','days','days','yrs','yrs','a','']
        yaxis = ['RV','RV','Relative Flux','Relative Flux','','','power','power','power','power','power','power','a','e','a','']       
        xunit = ['d' ,'d','d','d','d','','','','','','','','','','au','']
        yunit = ['m/s' ,'m/s' , '','','','','','','','','','','','','au','']

        zzz = [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,pe]
        font=QtGui.QFont()
        font.setPixelSize(12) 
        for i in range(len(zzz)):
 

                zzz[i].getAxis("bottom").tickFont = font
                zzz[i].getAxis("bottom").setStyle(tickTextOffset = 12)
                zzz[i].getAxis("top").tickFont = font
                zzz[i].getAxis("top").setStyle(tickTextOffset = 12)
                zzz[i].getAxis("left").tickFont = font
                zzz[i].getAxis("left").setStyle(tickTextOffset = 12)
                zzz[i].getAxis("right").tickFont = font
                zzz[i].getAxis("right").setStyle(tickTextOffset = 12)
                zzz[i].getAxis('left').setWidth(50)
                zzz[i].getAxis('right').setWidth(10)
                zzz[i].getAxis('top').setHeight(10)
                zzz[i].getAxis('bottom').setHeight(50)
                #zzz[i].getAxis('bottom').setScale(1e6)
                            
                zzz[i].setLabel('bottom', '%s'%xaxis[i], units='%s'%xunit[i],  **{'font-size':'12pt'})
                zzz[i].setLabel('left',   '%s'%yaxis[i], units='%s'%yunit[i],  **{'font-size':'12pt'})       
                zzz[i].showAxis('top') 
                zzz[i].showAxis('right') 
                zzz[i].getAxis('bottom').enableAutoSIPrefix(enable=False)


        p15.getViewBox().setAspectLocked(True)

        return   
        
        
        
        
    def identify_power_peaks(self,x,y,sig_level=np.array([]), power_level=np.array([])):
 
        per_ind = argrelextrema(y, np.greater)
        per_x   = x[per_ind]
        per_y   = y[per_ind]     

        peaks_sort = sorted(range(len(per_y)), key=lambda k: per_y[k], reverse=True)

        per_x   = per_x[peaks_sort]   
        per_y   = per_y[peaks_sort]  
        
        ################## text generator #################
        text_peaks = """ 
"""
        if power_level.size != 0 and sig_level.size != 0:
         
            text_peaks = text_peaks +"""FAP levels
-----------------------------------  
"""        
            for ii in range(len(power_level)):     
                text_peaks = text_peaks +"""
%.2f per cent = %.4f"""%(power_level[ii]*100.0,sig_level[ii])       
        
        text_peaks = text_peaks + """
----------------------------------------------
The 10 strongest peaks
----------------------------------------------
"""         
        for j in range(10):
            text_peaks = text_peaks +"""
period = %.2f [d], power = %.4f"""%(per_x[j],per_y[j])  
            if sig_level.size != 0 and per_y[j] > sig_level[-1]:
                text_peaks = text_peaks +"""  significant"""
                
        ################################################        
    
        return text_peaks  
        
        
######################## Correlation plots ######################################        
    def init_correlations_combo(self):
        global fit
        self.comboBox_corr_1.clear()
        self.comboBox_corr_2.clear()
       
        if fit.filelist.ndset > 0:
            for i in range(0,fit.filelist.ndset,1):

                self.comboBox_corr_1.addItem('RV %s'%(i+1),i+1)       
                self.comboBox_corr_2.addItem('RV %s'%(i+1),i+1) 
            
        for i in range(0,10,1):
            
            if len(fit.act_data_sets[i]) != 0: 
                self.comboBox_corr_1.addItem('act. data %s'%(i+1),i+1)       
                self.comboBox_corr_2.addItem('act. data %s'%(i+1),i+1) 
            
          
            
 ######################## Correlation plots END ##################################         
       
        
        
######################## Activity plots ######################################        
    def init_activity_combo(self):
        global fit
        
        for i in range(10):
            self.comboBox_act_data_gls.addItem('act. data %s'%(i+1),i+1)       
            self.comboBox_act_data.addItem('act. data %s'%(i+1),i+1)       
        
        
        #self.comboBox_act_data_gls.activated.connect(self.handleActivated_act_gls) 
        
        
   # def handleActivated_act_gls(self, index):
   #     global fit 
    
   #     ind = self.comboBox_act_data_gls.itemData(index)         
   #     self.update_activity_gls_plots(ind-1)
        
        
        
 
    def update_activity_gls_plots(self,ind):
        global fit, colors,  p11 
 
        omega = 1/ np.logspace(-0.05, 4, num=1000)
        power_levels = np.array([0.1,0.01,0.001])
  
        if len(fit.act_data_sets[ind]) != 0 and len(fit.act_data_sets[ind][0]) > 5:

            p11.plot(clear=True,)        
 
            act_per = gls.Gls((fit.act_data_sets[ind][0], fit.act_data_sets[ind][1],fit.act_data_sets[ind][2]), 
            fast=True,  verbose=False, norm= "ZK",ofac=5, fbeg=omega[999], fend=omega[ 0],)
            
            ######################## GLS ##############################
            if self.radioButton_act_GLS_period.isChecked():
                p11.setLogMode(True,False)        
                p11.plot(1/act_per.freq, act_per.power,pen=colors[ind],symbol=None ) 
                p11.setLabel('bottom', 'days', units='',  **{'font-size':'12pt'}) 
               # print("I ahve been here no")

            else:
                p11.setLogMode(False,False)        
                p11.plot(act_per.freq, act_per.power,pen=colors[ind],symbol=None )                    
                p11.setLabel('bottom', 'frequency', units='',  **{'font-size':'12pt'}) 
                #print("I ahve been here yes")

                                               
            [p11.addLine(x=None, y=fap, pen=pg.mkPen('k', width=0.8, style=QtCore.Qt.DotLine)) for ii,fap in enumerate(act_per.powerLevel(np.array(power_levels)))]
  
            self.act_periodogram_print_info.clicked.connect(lambda: self.print_info_for_object(
            act_per.info(stdout=False) + 
            self.identify_power_peaks(1/act_per.freq, act_per.power, power_level = power_levels, sig_level = act_per.powerLevel(np.array(power_levels)) )))   
    
            return
        else:   
            p11.plot(clear=True,)        

            return



    def update_activity_data_plots(self,ind):
        global fit, colors,  p5 
 
 
  
        if len(fit.act_data_sets[ind]) != 0:

            p5.plot(clear=True,)  
            
            
            err1 = pg.ErrorBarItem(x=fit.act_data_sets[ind][0], y=fit.act_data_sets[ind][1],symbol='o', 
            height=fit.act_data_sets[ind][2], beam=0.0, pen=colors[ind])  

            p5.addItem(err1)      
            p5.addLine(x=None, y=0, pen=pg.mkPen('#ff9933', width=0.8))
 

            p5.plot(fit.act_data_sets[ind][0],fit.act_data_sets[ind][1], pen=None,symbol='o',
            #symbolPen=,
            symbolSize=6,enableAutoRange=True,viewRect=True,
            symbolBrush=colors[ind]
            )      
            return
        else:   
            p5.plot(clear=True,)        

            return  

######################## Activity plots END ######################################        




    def update_RV_GLS_plots(self):
        global fit, colors, RV_per, p7 
 
        p7.plot(clear=True,)        
                          
        omega = 1/ np.logspace(-0.05, 4, num=1000)
        power_levels = np.array([0.1,0.01,0.001])
  
        if len(fit.fit_results.rv_model.jd) > 5:

            RV_per = gls.Gls((fit.fit_results.rv_model.jd, fit.fit_results.rv_model.rvs, fit.fit_results.rv_model.rv_err), 
            fast=True,  verbose=False, norm= "ZK",ofac=5, fbeg=omega[999], fend=omega[ 0],)
            
            ######################## GLS ##############################
            if self.radioButton_RV_GLS_period.isChecked():
                p7.setLogMode(True,False)        
                p7.plot(1/RV_per.freq, RV_per.power,pen='r',symbol=None ) 
                p7.setLabel('bottom', 'days', units='',  **{'font-size':'12pt'})                
            else:
                p7.setLogMode(False,False)        
                p7.plot(RV_per.freq, RV_per.power,pen='r',symbol=None )                    
                p7.setLabel('bottom', 'frequency', units='',  **{'font-size':'12pt'}) 

                                               
            [p7.addLine(x=None, y=fap, pen=pg.mkPen('k', width=0.8, style=QtCore.Qt.DotLine)) for ii,fap in enumerate(RV_per.powerLevel(np.array(power_levels)))]
 
 
            self.RV_periodogram_print_info.clicked.connect(lambda: self.print_info_for_object(
            RV_per.info(stdout=False) + 
            self.identify_power_peaks(1/RV_per.freq, RV_per.power, power_level = power_levels, sig_level = RV_per.powerLevel(np.array(power_levels)) )))   
    
    
    

    def update_RV_o_c_GLS_plots(self):
        global fit, colors, RV_per_res
        global p8  
 
        p8.plot(clear=True,)  
         
        omega = 1/ np.logspace(-0.05, 4, num=1000)
        power_levels = np.array([0.1,0.01,0.001])

        if len(fit.fit_results.rv_model.jd) > 5:
 
            RV_per_res = gls.Gls((fit.fit_results.rv_model.jd, fit.fit_results.rv_model.o_c, fit.fit_results.rv_model.rv_err), 
            fast=True,  verbose=False, norm= "ZK",ofac=5, fbeg=omega[999], fend=omega[ 0],)            
                
 
            ######################## GLS o-c ##############################
            if self.radioButton_RV_o_c_GLS_period.isChecked():
                p8.setLogMode(True,False)        
                p8.plot(1/RV_per_res.freq, RV_per_res.power,pen='r',symbol=None ) 
                p8.setLabel('bottom', 'days', units='',  **{'font-size':'12pt'})
            else:
                p8.setLogMode(False,False)        
                p8.plot(RV_per_res.freq, RV_per_res.power,pen='r',symbol=None )    
                p8.setLabel('bottom', 'frequency', units='',  **{'font-size':'12pt'})                
                
            
            [p8.addLine(x=None, y=fap, pen=pg.mkPen('k', width=0.8, style=QtCore.Qt.DotLine)) for ii,fap in enumerate(RV_per_res.powerLevel(np.array(power_levels)))]            


            self.RV_res_periodogram_print_info.clicked.connect(lambda: self.print_info_for_object(RV_per_res.info(stdout=False)+
            self.identify_power_peaks(1/RV_per_res.freq, RV_per_res.power, power_level = power_levels, sig_level = RV_per.powerLevel(np.array(power_levels)) ) )  )      
 
    def update_WF_plots(self):
        global fit, colors, RV_per_res
        global p12  
 
        p12.plot(clear=True,) 
        p12.setLogMode(True,False)
                        
        omega = 1/ np.logspace(-0.05, 4, num=1000)
        power_levels = np.array([0.1,0.01,0.001])
        
        if len(fit.fit_results.rv_model.jd) > 5:
            ######################## DFT (Window) ##############################
            WF_power = []
            for omi in 2*np.pi*omega: 
                phase = (fit.fit_results.rv_model.jd-fit.fit_results.rv_model.jd[0]) * omi                 
                WC = np.sum(np.cos(phase))
                WS = np.sum(np.sin(phase))
                WF_power.append((WC**2 + WS**2)/len(fit.fit_results.rv_model.jd)**2) 

            WF_power = np.array(WF_power)
            p12.plot(1/np.array(omega), WF_power,pen='k',symbol=None )   
                        
            self.WF_print_info.clicked.connect(lambda: self.print_info_for_object(self.identify_power_peaks(1/np.array(omega), WF_power)))        
         
         
        

    def update_RV_plots(self):
        global fit, colors
        global p1,p2
        

       # brush_list = [pg.mkColor(c) for c in [colors[i] for i in fit.filelist.idset]]
 
        p1.plot(clear=True,)
        p2.plot(clear=True,)
        
       # print(fit.filelist.idset)

        #inf1 = pg.InfiniteLine(movable=False, angle=0, label=None, span=(0, 1), 
        #              labelOpts={'position':0.0, 'color': 'k', 'fill': (200,200,200,50), 'movable': False} )
        #p1.addItem(inf1)    

        if self.jitter_to_plots.isChecked():
            error_list = self.add_jitter(fit.fit_results.rv_model.rv_err, fit.filelist.idset)
        else:
            error_list = fit.fit_results.rv_model.rv_err
          
                      
    
        p1.addLine(x=None, y=0, pen=pg.mkPen('#ff9933', width=0.8))
 

        p1.plot(fit.fit_results.model_jd,fit.fit_results.model, 
        pen={'color': 0.5, 'width': 1.1},enableAutoRange=True,viewRect=True, labels =  {'left':'RV', 'bottom':'JD'}) 
        
        for i in range(max(fit.filelist.idset)+1):
            p1.plot(fit.fit_results.rv_model.jd[fit.filelist.idset==i],fit.fit_results.rv_model.rvs[fit.filelist.idset==i], 
            pen=None, #{'color': colors[i], 'width': 1.1},
            symbol='o',
            symbolPen={'color': colors[i], 'width': 1.1},
            symbolSize=6,enableAutoRange=True,viewRect=True,
            symbolBrush=colors[i]
            )        
            err1 = pg.ErrorBarItem(x=fit.fit_results.rv_model.jd[fit.filelist.idset==i], 
                                   y=fit.fit_results.rv_model.rvs[fit.filelist.idset==i],symbol='o', 
            height=error_list[fit.filelist.idset==i], beam=0.0, pen=colors[i])  

            p1.addItem(err1)  

  

        p2.addLine(x=None, y=0, pen=pg.mkPen('#ff9933', width=0.8))
        
        for i in range(max(fit.filelist.idset)+1):
            p2.plot(fit.fit_results.rv_model.jd[fit.filelist.idset==i],fit.fit_results.rv_model.o_c[fit.filelist.idset==i], 
            pen=None, #{'color': colors[i], 'width': 1.1},
            symbol='o',
            symbolPen={'color': colors[i], 'width': 1.1},
            symbolSize=6,enableAutoRange=True,viewRect=True,
            symbolBrush=colors[i]
            )        
            err2 = pg.ErrorBarItem(x=fit.fit_results.rv_model.jd[fit.filelist.idset==i], 
                                   y=fit.fit_results.rv_model.o_c[fit.filelist.idset==i],symbol='o', 
            height=error_list[fit.filelist.idset==i], beam=0.0, pen=colors[i])  

            p2.addItem(err2)  
        
  #pen={'color': 'r', 'width': 1.1},  
  
      
      
     
        
    def update_plots(self):
        self.update_RV_GLS_plots()
        self.update_RV_o_c_GLS_plots()    
        self.update_WF_plots()                
        self.update_RV_plots()
        self.update_extra_plots()
        self.update_orb_plot()
        #self.change_extra_plot()
        
        

################################ RV files #######################################################
        
    def showDialog_fortran_input_file(self):
        global fit
 
        input_files = QtGui.QFileDialog.getOpenFileName(self, 'Open session', '', 'Data (*.init)')
        #print(input_files[0])
        if str(input_files[0]) != '':
            fit=rv.signal_fit(str(input_files[0]), 'Test',readinputfile=True)
            self.update_use_from_input_file()
            self.init_fit()
            self.update_RV_file_buttons()

    def showDialog_RV_input_file(self):
        global fit

        but_ind = self.buttonGroup_4.checkedId()   
        input_files = QtGui.QFileDialog.getOpenFileName(self, 'Open RV data', '', 'Data (*.vels)')
       # print(input_files[0])        
        if str(input_files[0]) != '':
 
            fit.add_dataset('test', str(input_files[0]),0.0,1.0)
            self.init_fit()            
            self.update_use_from_input_file()            
            self.update_use()
            self.update_params()
            self.update_RV_file_buttons()

    def remove_RV_file(self):
        global fit

        but_ind = self.buttonGroup_remove_RV_data.checkedId()   
        fit.remove_dataset(but_ind -1)
        self.init_fit()         
        self.update_use_from_input_file()   
        self.update_use()
        self.update_gui_params()
        self.update_params()
        self.update_RV_file_buttons()

    def update_RV_file_buttons(self):
        global fit, colors          

        for i in range(10):
            if i < fit.filelist.ndset:
                self.buttonGroup_4.button(i+1).setStyleSheet("color: %s;"%colors[i])
                self.buttonGroup_remove_RV_data.button(i+1).setStyleSheet("color: %s;"%colors[i])
            else:
                self.buttonGroup_4.button(i+1).setStyleSheet("")
                self.buttonGroup_remove_RV_data.button(i+1).setStyleSheet("")
                #"background-color: #333399;""background-color: yellow;" "selection-color: yellow;"  "selection-background-color: blue;")               
        self.init_correlations_combo()


################################ RV files END #######################################################


################################ transit files #######################################################
        

    def showDialog_tra_input_file(self):
        global fit

        but_ind = self.buttonGroup_transit_data.checkedId()   
        input_files = QtGui.QFileDialog.getOpenFileName(self, 'Open Transit data', '', 'Data (*.dat)')

        if str(input_files[0]) != '':
 
            fit.add_transit_dataset('test', str(input_files[0]),tra_idset =but_ind-1)
            #self.init_fit()            
            #self.update_use_from_input_file()            
            #self.update_use()
            #self.update_params()
            self.update_tra_file_buttons()

    def remove_tra_file(self):
        global fit

        but_ind = self.buttonGroup_remove_transit_data.checkedId()   
        fit.remove_transit_dataset(but_ind -1)
       # self.init_fit()         
      #  self.update_use_from_input_file()   
      #  self.update_use()
      #  self.update_gui_params()
     #   self.update_params()
        self.update_tra_file_buttons()

    def update_tra_file_buttons(self):
        global fit, colors          

        for i in range(10):
            if len(fit.tra_data_sets[i]) != 0:
                self.buttonGroup_transit_data.button(i+1).setStyleSheet("color: %s;"%colors[i])
                self.buttonGroup_remove_transit_data.button(i+1).setStyleSheet("color: %s;"%colors[i])
            else:
                self.buttonGroup_transit_data.button(i+1).setStyleSheet("")
                self.buttonGroup_remove_transit_data.button(i+1).setStyleSheet("")
                #"background-color: #333399;""background-color: yellow;" "selection-color: yellow;"  "selection-background-color: blue;")               

 
################################ transit files END #######################################################


################################ activity files #######################################################
        

    def showDialog_act_input_file(self):
        global fit

        but_ind = self.buttonGroup_activity_data.checkedId()   
        input_files = QtGui.QFileDialog.getOpenFileName(self, 'Open Activity data', '', 'Data (*.dat)')

        if str(input_files[0]) != '':
 
            fit.add_act_dataset('test', str(input_files[0]),act_idset =but_ind-1)
            #self.init_fit()            
            #self.update_use_from_input_file()            
            #self.update_use()
            #self.update_params()
            self.update_act_file_buttons()
            self.update_activity_gls_plots(but_ind-1)
            #self.handleActivated_act_gls(but_ind-1)
            
    def remove_act_file(self):
        global fit

        but_ind = self.buttonGroup_remove_activity_data.checkedId()   
        fit.remove_act_dataset(but_ind -1)
       # self.init_fit()         
      #  self.update_use_from_input_file()   
      #  self.update_use()
      #  self.update_gui_params()
     #   self.update_params()
        self.update_act_file_buttons()

    def update_act_file_buttons(self):
        global fit, colors          

        for i in range(10):
            if len(fit.act_data_sets[i]) != 0:
                self.buttonGroup_activity_data.button(i+1).setStyleSheet("color: %s;"%colors[i])
                self.buttonGroup_remove_activity_data.button(i+1).setStyleSheet("color: %s;"%colors[i])
            else:
                self.buttonGroup_activity_data.button(i+1).setStyleSheet("")
                self.buttonGroup_remove_activity_data.button(i+1).setStyleSheet("")
                #"background-color: #333399;""background-color: yellow;" "selection-color: yellow;"  "selection-background-color: blue;")               
        self.init_correlations_combo()

################################ activity files END #######################################################

    def init_fit(self): 
        global fit
        minimize_fortran=True
        fit.fitting(fileinput=False,outputfiles=[1,1,1], minimize_fortran=minimize_fortran,  fortran_kill=30, timeout_sec=300,minimize_loglik=True,amoeba_starts=0, print_stat=False, eps=self.dyn_model_accuracy.value(), dt=self.time_step_model.value(), npoints=self.points_to_draw_model.value(), model_max= self.model_max_range.value())
        
        self.update_labels()
        self.update_gui_params()
        self.update_errors() 
        self.update_a_mass() 
        
        #start_time = time.time()
        self.update_plots() 
        #print("--- %s seconds ---" % (time.time() - start_time))      
     
        self.jupiter_push_vars()       
        
        
    def update_orb_plot(self):
        global fit, p15
        
        p15.plot(clear=True,)    
        
        for i in range(fit.npl):
            orb_xyz, pl_xyz, peri_xyz, apo_xyz = rv.planet_orbit_xyz(fit,i)        
            p15.plot(orb_xyz[0],orb_xyz[1], pen={'color': 0.5, 'width': 1.1},enableAutoRange=True,viewRect=True)   
            p15.plot((0,peri_xyz[0]),(0,peri_xyz[1]), pen={'color': 0.5, 'width': 1.1},enableAutoRange=True,viewRect=True)               
            
            p15.plot((pl_xyz[0],pl_xyz[0]), (pl_xyz[1],pl_xyz[1] ), pen=None,symbol='o', symbolSize=6,enableAutoRange=True,viewRect=True, symbolBrush='b') 
            
                           
        
        p15.plot(np.array([0,0]), np.array([0,0]), pen=None,symbol='o', symbolSize=8,enableAutoRange=True,viewRect=True, symbolBrush='r')                
        
     
        
    def update_extra_plots(self):
        global fit

        self.comboBox_extra_plot.clear()
        self.comboBox_extra_plot.setObjectName("which plot")        

        if fit.npl != 0:
            for i in range(fit.npl):
                self.comboBox_extra_plot.addItem('phase pl %s'%(i+1),i+1)
            
            self.comboBox_extra_plot.addItem('gls',fit.npl+1)
            self.comboBox_extra_plot.addItem('gls o-c',fit.npl+1)
           
            
            self.phase_plots(1)   
            
        self.comboBox_extra_plot.activated.connect(self.handleActivated)        
        
    def handleActivated(self, index):
        global fit, pe
        
        ind = self.comboBox_extra_plot.itemData(index) 
        
        if ind <= fit.npl:
            self.phase_plots(ind)
        elif ind >= fit.npl+1:
            pe.clear()
            #pe = pg.ViewBox()
            #pe.setXLink(zzz[ind])
            #pe.setYLink(zzz[ind])
            #pe.addItem(p8.getPlotItem)
            #pe.setPlotItem(p8.getPlotItem)
            #pe.scene().addItem(p8.plotItem) 
            #pe = dill.copy(p8)
            #pe.addItem(p8.getPlotItem)             
        else:
   
            return
            
            
  
    def add_jitter(self, errors, ind):
        global fit
       
        errors_with_jitt = np.array([np.sqrt(errors[i]**2 + fit.params.jitters[ii]**2)  for i,ii in enumerate(ind)])
        
        return errors_with_jitt




    def phase_plots(self, ind):
        global fit, colors   
        
        pe.plot(clear=True,)    
        
        
        ph_data,ph_model = rv.phase_planet_signal(fit,ind)

             

        if len(ph_data) == 1:
            return
 
       # brush_list2 = [pg.mkColor(c) for c in [colors[i] for i in ph_data[3]]]

        if self.jitter_to_plots.isChecked():
            error_list = self.add_jitter(ph_data[2], ph_data[3])
        else:
            error_list = ph_data[2]
        
        
        pe.addLine(x=None, y=0, pen=pg.mkPen('#ff9933', width=0.8))   
        pe.plot(ph_model[0],ph_model[1], pen={'color': 0.5, 'width': 2.0},
        enableAutoRange=True,viewRect=True, labels =  {'left':'RV', 'bottom':'JD'})   
        
        for i in range(max(ph_data[3])+1):
        
            pe.plot(ph_data[0][ph_data[3]==i],ph_data[1][ph_data[3]==i],             
            pen=None, #{'color': colors[i], 'width': 1.1},
            symbol='o',
            symbolPen={'color': colors[i], 'width': 1.1},
            symbolSize=6,enableAutoRange=True,viewRect=True,
            symbolBrush=colors[i]
            )  
               
            err_ = pg.ErrorBarItem(x=ph_data[0][ph_data[3]==i], y=ph_data[1][ph_data[3]==i],
            symbol='o', height=error_list[ph_data[3]==i], beam=0.0, pen=colors[i])   
         
            pe.addItem(err_)
        
        pe.setLabel('bottom', 'days', units='',  **{'font-size':'12pt'})
        pe.setLabel('left',   'RV', units='m/s',  **{'font-size':'12pt'})  
      
        
 
###################################################### 
    def run_batman_test(self): 
        global fit, p3
    
        p3.plot(clear=True,)    
        # from the example in github
        params = batman.TransitParams()       #object to store transit parameters
        params.t0  = self.t0_1_trans.value()                       #time of inferior conjunction
        params.per = self.P1_trans.value()  #orbital period
        params.ecc = self.e1_trans.value()                     
        params.rp  = self.pl1_radii.value()                       #planet radius (in units of stellar radii)
        params.a   = self.a1_trans.value()                        #semi-major axis (in units of stellar radii)
        params.inc = self.incl1_trans.value()                       #orbital inclination (in degrees)
        params.w   = self.om1_trans.value()                                 #longitude of periastron (in degrees)
        
        
        params.limb_dark = "nonlinear"        #limb darkening model
        params.u = [self.u1_1_trans.value(), self.u2_1_trans.value() , 0.1, -0.1]      #limb darkening coefficients [u1, u2, u3, u4]

        if len(fit.tra_data_sets[0]) != 0:
            t = fit.tra_data_sets[0][0]
        else:    
            t = np.linspace(-0.25, 0.25, 1000)  #times at which to calculate light curve
        
        m = batman.TransitModel(params, t)    #initializes model
 
        flux = m.light_curve(params)          #calculates light curve
 
        p3.plot(t, flux,pen='k',symbol=None )     
       
      
        
############################# N-Body ########################################        

        
    def worker_Nbody_complete(self):
        global fit, colors, p13, p14  

        p13.plot(clear=True,)
        p14.plot(clear=True,)


        for i in range(fit.npl):
            p13.plot(fit.evol_T[i], fit.evol_a[i] ,pen=colors[i],symbol=None )     
            p14.plot(fit.evol_T[i], fit.evol_e[i] ,pen=colors[i],symbol=None )   
            
            
        self.button_orb_evol.setEnabled(True)       
        self.statusBar().showMessage('')           
          
 
    def worker_Nbody(self):
        global fit  
        
        
        self.button_orb_evol.setEnabled(False)         
        
        # check if any fits where performed, and tus planets present
        if fit.fit_results.mass <= 0:
             choice = QtGui.QMessageBox.information(self, 'Warning!',
             "Not possible to integrate a fit that does not exist. First perform an orbital fitting and then test the orbital stability. Okay?", QtGui.QMessageBox.Ok)      
             self.button_orb_evol.setEnabled(True)         
             return        
        
        
        if fit.npl < 2:
            choice = QtGui.QMessageBox.information(self, 'Warning!'," With less than two planets this makes no sense. Okay?",
                                            QtGui.QMessageBox.Ok) 
            return
 

        self.statusBar().showMessage('Running Orbital Evolution......')   
        
        # Pass the function to execute
        worker3 = Worker(lambda:  self.run_orbital_simulations()) # Any other args, kwargs are passed to the run  
        # Execute
        worker3.signals.finished.connect(self.worker_Nbody_complete)
        
        # worker.signals.result.connect(self.print_output)
        #worker.signals.finished.connect(self.thread_complete)
       # worker.signals.progress.connect(self.progress_fn)
        
        
        self.threadpool.start(worker3)  

 
    def run_orbital_simulations(self):
        global fit
        
     

        self.max_time_of_evol

        if self.radioButton_SyMBA.isChecked():
            integrator = 'symba'
        elif self.radioButton_MVS.isChecked():
            integrator = 'mvs'        
        elif self.radioButton_MVS_GR.isChecked():       
             integrator = 'mvs_gr'       
        
        fit.run_stability_last_fit_params(timemax=self.max_time_of_evol.value(), timestep=self.time_step_of_evol.value(), integrator=integrator)      
        
 
        
        
############################# END N-Body #####################################                
               
       
       
        
        
############################# Fortran fitting ###############################        
        
    def worker_RV_fitting_complete(self):
        global fit  
        #fit.print_info(short_errors=False)
        
        self.update_labels()
        self.update_gui_params()
        self.update_errors() 
        self.update_a_mass()                    
        self.update_plots()                   
        self.statusBar().showMessage('')   
        
        self.run_batman_test()   
              
        self.jupiter_push_vars()   
        self.button_fit.setEnabled(True)         
 
    def worker_RV_fitting(self, ff=20,m_ln=True, auto_fit = False ):
        global fit  
        
        self.button_fit.setEnabled(False)         
        
        # check if RV data is present
        if fit.filelist.ndset <= 0:
             choice = QtGui.QMessageBox.information(self, 'Warning!',
             "Not possible to look for planets if there are no RV data loaded. Please add your RV data first. Okay?", QtGui.QMessageBox.Ok)      
             self.button_fit.setEnabled(True)         
             return        
        
 
        
        if self.radioButton_fortran77.isChecked():
            self.statusBar().showMessage('Minimizing parameters....')                 
            # Pass the function to execute
            worker2 = Worker(lambda:  self.optimize_fit(ff=ff,  minimize_fortran=True, m_ln=self.amoeba_radio_button.isChecked(), auto_fit = auto_fit)) # Any other args, kwargs are passed to the run  
        else:
            
            gp_params = [self.GP_rot_kernel_Amp.value(),
                     self.GP_rot_kernel_time_sc.value(),
                     self.GP_rot_kernel_Per.value(),
                     self.GP_rot_kernel_fact.value()]

            use_gp_params = [self.use_GP_rot_kernel_Amp.isChecked(),
                     self.use_GP_rot_kernel_time_sc.isChecked(),
                     self.use_GP_rot_kernel_Per.isChecked(),
                     self.use_GP_rot_kernel_fact.isChecked()]
 
             
            
            
            self.statusBar().showMessage('Minimizing parameters using SciPyOp (might be slow)....')                 
            worker2 = Worker(lambda:  self.optimize_fit(ff=1, doGP=self.goGP.isChecked(), gp_par=np.array(gp_params),use_gp_par=np.array(use_gp_params), gp_kernel_id=-1, minimize_fortran=False, m_ln=self.amoeba_radio_button.isChecked(), auto_fit = auto_fit)) # Any other args, kwargs are passed to the run  
               # Execute
            
            
        worker2.signals.finished.connect(self.worker_RV_fitting_complete)
        
        # worker.signals.result.connect(self.print_output)
        #worker.signals.finished.connect(self.thread_complete)
       # worker.signals.progress.connect(self.progress_fn)
        
        
        self.threadpool.start(worker2)       
     
                     
        
    def optimize_fit(self,ff=20,m_ln=True, doGP=False, gp_par=None, gp_kernel_id=-1, use_gp_par=[True,True,True,True], auto_fit = False, minimize_fortran=True):  
        global fit
        
        if not auto_fit:
        #self.update_use()
            self.update_params()
            
            

            
        if self.radioButton_Dynamical.isChecked():
            fit.mod_dynamical = True
            f_kill = self.dyn_model_to_kill.value()
            ff = 1
        else:
            fit.mod_dynamical = False
            f_kill = self.kep_model_to_kill.value()        

        if m_ln and minimize_fortran:
            if ff > 0:        
                """
                run one time using the L-M method ignorring the jitter (for speed)
                """
                fit.fitting(fileinput=False,outputfiles=[1,0,0], doGP=doGP, gp_par=use_gp_par, use_gp_par=use_gp_par, kernel_id=gp_kernel_id, minimize_fortran=minimize_fortran, fortran_kill=f_kill, timeout_sec=300,minimize_loglik=False,amoeba_starts=ff, print_stat=False, eps=self.dyn_model_accuracy.value(), dt=self.time_step_model.value())
            """
            now run the amoeba code modeling the jitters
            """
            fit.fitting(fileinput=False,outputfiles=[1,0,0], doGP=doGP, gp_par=use_gp_par, use_gp_par=use_gp_par, kernel_id=gp_kernel_id,  minimize_fortran=minimize_fortran,  fortran_kill=f_kill, timeout_sec=300,minimize_loglik=True,amoeba_starts=ff, print_stat=False, eps=self.dyn_model_accuracy.value(), dt=self.time_step_model.value())
            fit.fitting(fileinput=False,outputfiles=[1,1,1], doGP=doGP, gp_par=use_gp_par, use_gp_par=use_gp_par, kernel_id=gp_kernel_id,  minimize_fortran=minimize_fortran,  fortran_kill=f_kill, timeout_sec=300,minimize_loglik=True,amoeba_starts=0, print_stat=False, eps=self.dyn_model_accuracy.value(), dt=self.time_step_model.value(), npoints=self.points_to_draw_model.value(), model_max= self.model_max_range.value())

       # elif m_ln and not minimize_fortran:       
      #      fit.fitting(fileinput=False,outputfiles=[1,1,1], doGP=doGP, gp_par=use_gp_par, use_gp_par=use_gp_par, kernel_id=gp_kernel_id,  minimize_fortran=minimize_fortran,  fortran_kill=f_kill, timeout_sec=300,minimize_loglik=True,amoeba_starts=0, print_stat=False, eps=self.dyn_model_accuracy.value(), dt=self.time_step_model.value(), npoints=self.points_to_draw_model.value(), model_max= self.model_max_range.value())
        
        else:        
                fit.fitting(fileinput=False,outputfiles=[1,1,1], doGP=doGP,gp_par=use_gp_par, use_gp_par=use_gp_par, kernel_id=gp_kernel_id,  minimize_fortran=minimize_fortran, fortran_kill=f_kill, timeout_sec=300,minimize_loglik=m_ln,amoeba_starts=ff, print_stat=False,eps=self.dyn_model_accuracy.value(), dt=self.time_step_model.value(), npoints=self.points_to_draw_model.value(), model_max= self.model_max_range.value())


        if auto_fit:
            self.update_labels()
            self.update_gui_params()
            self.update_errors() 
            self.update_a_mass()                    
            self.update_plots()                   
            self.statusBar().showMessage('')           
            self.run_batman_test()   
              
            self.jupiter_push_vars()



############################# END Fortran fitting ###############################         
       

    def print_info_for_object(self,text):
        #self.dialog.statusBar().showMessage('Ready')
        self.dialog.setGeometry(300, 300, 450, 250)
        self.dialog.setWindowTitle('Detailed Info')  
 
        
        self.dialog.text.setPlainText(text)
        self.dialog.text.setReadOnly(True)       
        #self.dialog.setWindowIcon (QtGui.QIcon('logo.png'))        
        
        self.dialog.show()

    def print_info_credits(self, image=False):
        #self.dialog.statusBar().showMessage('Ready')
        self.dialog_credits.setFixedSize(900, 900)
        self.dialog_credits.setWindowTitle('Credits')  
        #self.dialog.setGeometry(300, 300, 800, 800)
        #self.dialog_credits.acceptRichText(True)
        
        text = ''
        self.dialog_credits.text.setText(text) 
        
        text = "You are using 'The Exo-Striker' (ver. 0.01) \n developed by 3fon3fonov"
        
        self.dialog_credits.text.append(text)

        text = "\n"*15 +"CREDITS:"+"\n"*2 + "This tool uses the publically \n available packages: \n" 
        self.dialog_credits.text.append(text)
        
    

        text = "* " + "<a href='https://github.com/pyqtgraph/pyqtgraph'>pyqtgraph</a>"
        self.dialog_credits.text.append(text)

        text = "* " + "<a href='https://github.com/dfm/emcee'>emcee</a>" 
        self.dialog_credits.text.append(text) 
        
        text = "* " + "<a href='https://github.com/dfm/celerite'>celerite</a>" 
        self.dialog_credits.text.append(text)  
                                
        text = "* " + "<a href='https://github.com/lkreidberg/batman'>batman-package</a>" 
        self.dialog_credits.text.append(text)

        text = "* " + "<a href='https://www.boulder.swri.edu/~hal/swift.html'>swift</a>" 
        self.dialog_credits.text.append(text)        
                        
        text = "* " + "<a href='https://github.com/jupyter/qtconsole'>qtconsole</a>"
        self.dialog_credits.text.append(text)        

        text = "* " + "<a href='https://github.com/mfitzp/15-minute-apps/tree/master/wordprocessor'>megasolid idiom</a>" 
        self.dialog_credits.text.append(text)    
        
        text = "(A few more to be added) \n" 
        self.dialog_credits.text.append(text)   
        
 
        
        #self.dialog_credits.text.setText(text)
        #self.dialog_credits.text.insertHtml(text)
        
        
        text = "\n"*5 + """Note:
Please keep in mind that this software is developed 
mostly for my needs and for fun. I hope, however, 
that you may find it capable to solve your scientific 
problems, too. 

Feedback and help in further developing will be 
highly appreciated!
"""
        self.dialog_credits.text.append(text)   
  
        self.dialog_credits.text.setReadOnly(True)       
        
        self.dialog_credits.setStyleSheet(" QTextEdit{border-image: url(33_striker.png) 0 0 0 0 stretch stretch;} ")

        #self.dialog.setWindowIcon (QtGui.QIcon('logo.png'))        
        
        self.dialog_credits.show()
        
  

    def run_bootstrap(self):
        global fit
        choice = QtGui.QMessageBox.information(self, 'Warning!', "Not available yet. Okay?", QtGui.QMessageBox.Ok) 



    def find_planets(self):
        global fit,RV_per,RV_per_res
 
 
        # check if RV data is present
        if fit.filelist.ndset <= 0:
             choice = QtGui.QMessageBox.information(self, 'Warning!',
             "Not possible to look for planets if there are no RV data loaded. Please add your RV data first. Okay?", QtGui.QMessageBox.Ok)      
             self.button_auto_fit.setEnabled(True)         
             return        
 
 
        # the first one on the data GLS
        if RV_per.power.max() <= RV_per.powerLevel(self.auto_fit_FAP_level.value()):
             choice = QtGui.QMessageBox.information(self, 'Warning!',
             "No significant power on the GLS. Therefore no planets to fit OK?", QtGui.QMessageBox.Ok)      
             self.button_auto_fit.setEnabled(True)                                                           
             return
        
        else:
            if fit.npl !=0:
                for j in range(fit.npl):
                    fit.remove_planet(fit.npl-(j+1))

            mean_anomaly_from_gls = np.degrees((((fit.epoch - float(RV_per.hpstat["T0"]) )% (RV_per.hpstat["P"]) )/ (RV_per.hpstat["P"]) ) * 2*np.pi)
             
            fit.add_planet(RV_per.hpstat["amp"],RV_per.hpstat["P"],0.0,0.0,mean_anomaly_from_gls -90.0,90.0,0.0)
            fit.use.update_use_planet_params_one_planet(0,True,True,True,True,True,False,False)     
            self.update_use_from_input_file()   
            self.update_use()                     
            self.optimize_fit(20,m_ln=self.amoeba_radio_button.isChecked(),auto_fit = True)
            
            #now inspect the residuals
            
            
            for i in range(1,int(self.auto_fit_N_planets.value())):
                
                if RV_per_res.power.max() <= RV_per_res.powerLevel(self.auto_fit_FAP_level.value()):
                    for j in range(fit.npl):
                        fit.use.update_use_planet_params_one_planet(j,True,True,True,True,True,False,False)     
            
                    self.update_use_from_input_file()   
                    self.update_use()                     
                    self.optimize_fit(20,m_ln=self.amoeba_radio_button.isChecked(),auto_fit = True) 
                    self.button_auto_fit.setEnabled(True)     
                    return
                #elif (1/RV_per_res.hpstat["fbest"]) > 1.5:
                else:    
                    mean_anomaly_from_gls = np.degrees((((fit.epoch - float(RV_per_res.hpstat["T0"]) )% (RV_per_res.hpstat["P"]) )/ (RV_per_res.hpstat["P"]) ) * 2*np.pi)
             
                    fit.add_planet(RV_per_res.hpstat["amp"],RV_per_res.hpstat["P"],0.0,0.0,mean_anomaly_from_gls -90.0,90.0,0.0)
                    fit.use.update_use_planet_params_one_planet(i,True,True,True,True,True,False,False)  
                    
                    print(fit.params.planet_params[2 + 7*(i-1)])
                    #fit.use.update_use_planet_params_one_planet(i,True,True,True,True,True,False,False)  
                   
                    self.update_use_from_input_file()   
                    self.update_use()                     
                    self.optimize_fit(20,m_ln=self.amoeba_radio_button.isChecked(),auto_fit = True)  
                    
                #else:
                 #   continue
                    
                    
            for j in range(fit.npl):
                fit.use.update_use_planet_params_one_planet(j,True,True,True,True,True,False,False)     
    
            self.update_use_from_input_file()   
            self.update_use()                     
            self.optimize_fit(20,m_ln=self.amoeba_radio_button.isChecked(),auto_fit = True)                 
                
 
        self.button_auto_fit.setEnabled(True)   
 

    def run_auto_fit(self):
        global fit 
        
        self.radioButton_Keplerian.setChecked(True) # this is to be fixed! Only with keplerian fitting th autofit works fine so far.
        #self.button_fit.setEnabled(False) 
        self.button_auto_fit.setEnabled(False)         
        
        if fit.npl != 0:        
            choice = QtGui.QMessageBox.information(self, 'Warning!',
                                            "Planets already exist. Do you want to overwrite the analysis?",
                                            QtGui.QMessageBox.No | QtGui.QMessageBox.Yes)  
         
            if choice == QtGui.QMessageBox.No:
                return
                self.button_auto_fit.setEnabled(True)         
            elif choice == QtGui.QMessageBox.Yes:
                self.find_planets()
        else:
            self.find_planets()
                



    def minimize_1param(self):
        global fit
        """
        This function must be completely refurbished!!! How to check 
        which QDoubleSpinBox is trigerred? Does everytime one needs to call 
        self.init_fit() ? ? ?
        
        """
        
        self.K1.minimize_signal.connect(lambda: fit.minimize_one_param_K(0)) #TBD!
        self.K1.minimize_signal.connect(self.init_fit) #TBD!       
        self.P1.minimize_signal.connect(lambda: fit.minimize_one_param_P(0)) #TBD!
        self.P1.minimize_signal.connect(self.init_fit) #TBD!     
        self.e1.minimize_signal.connect(lambda: fit.minimize_one_param_e(0)) #TBD!
        self.e1.minimize_signal.connect(self.init_fit) #TBD!  
        self.om1.minimize_signal.connect(lambda: fit.minimize_one_param_w(0)) #TBD!
        self.om1.minimize_signal.connect(self.init_fit) #TBD!  
        self.ma1.minimize_signal.connect(lambda: fit.minimize_one_param_M0(0)) #TBD!
        self.ma1.minimize_signal.connect(self.init_fit) #TBD!  
        
        self.K2.minimize_signal.connect(lambda: fit.minimize_one_param_K(1)) #TBD!
        self.K2.minimize_signal.connect(self.init_fit) #TBD!       
        self.P2.minimize_signal.connect(lambda: fit.minimize_one_param_P(1)) #TBD!
        self.P2.minimize_signal.connect(self.init_fit) #TBD!     
        self.e2.minimize_signal.connect(lambda: fit.minimize_one_param_e(1)) #TBD!
        self.e2.minimize_signal.connect(self.init_fit) #TBD!  
        self.om2.minimize_signal.connect(lambda: fit.minimize_one_param_w(1)) #TBD!
        self.om2.minimize_signal.connect(self.init_fit) #TBD!  
        self.ma2.minimize_signal.connect(lambda: fit.minimize_one_param_M0(1)) #TBD!
        self.ma2.minimize_signal.connect(self.init_fit) #TBD!         
        
        self.K3.minimize_signal.connect(lambda: fit.minimize_one_param_K(2)) #TBD!
        self.K3.minimize_signal.connect(self.init_fit) #TBD!       
        self.P3.minimize_signal.connect(lambda: fit.minimize_one_param_P(2)) #TBD!
        self.P3.minimize_signal.connect(self.init_fit) #TBD!     
        self.e3.minimize_signal.connect(lambda: fit.minimize_one_param_e(2)) #TBD!
        self.e3.minimize_signal.connect(self.init_fit) #TBD!  
        self.om3.minimize_signal.connect(lambda: fit.minimize_one_param_w(2)) #TBD!
        self.om3.minimize_signal.connect(self.init_fit) #TBD!  
        self.ma3.minimize_signal.connect(lambda: fit.minimize_one_param_M0(2)) #TBD!
        self.ma3.minimize_signal.connect(self.init_fit) #TBD!               
        
    def jupiter_push_vars(self):
        global fit        
        self.console_widget.push_vars({'fit':fit})    
        #self.console_widget.push_vars({'pg':pg})    
        #self.console_widget.push_vars({'p1':p1})    

        #self.console_widget.clear()         
        #self.console_widget.print_text(str("Welcome!")) 

########################## work in progress ##################################
 
    def getNewses(self):
        global fit, ses_list  
        
        text, okPressed = QtGui.QInputDialog.getText(self, "New session","Name session: (No space and special characters!)", QtGui.QLineEdit.Normal, "")
        if okPressed and text != '':
            
            if len(ses_list) == 0:
                ses_list.append(fit)
                
                
            file_pi = open('.sessions/empty.ses', 'rb')
            fit_new = dill.load(file_pi)
            file_pi.close()
            fit_new.name=text
            ses_list.append(fit_new)
 
            self.session_list()
            
    def rem_ses(self):
        global fit, ses_list  
        
        ind = self.comboBox_select_ses.currentIndex()
        
        choice = QtGui.QMessageBox.information(self, 'Warning!',
        "Do you really want to remove Session %s"%(ind+1),
                                            QtGui.QMessageBox.No | QtGui.QMessageBox.Yes)  
         
        if choice == QtGui.QMessageBox.No:
            return
        elif choice == QtGui.QMessageBox.Yes and ind <=0:
            self.new_session()
            #ses_list[0] = fit
        elif choice == QtGui.QMessageBox.Yes and ind > 0:
            ses_list.pop(ind)
            self.session_list() 
            self.select_session(ind-1)
            
    def cop_ses(self):
        global fit, ses_list  
        
        ind = self.comboBox_select_ses.currentIndex()
        
        if len(ses_list) == 0:
            ses_list.append(fit)
        if ind <=0:
            fit_new =dill.copy(ses_list[0])
        elif ind > 0:
            fit_new =dill.copy(ses_list[ind])
 
        
        ses_list.append(fit_new)
        self.session_list() 
        self.select_session(ind-1)            
  

    def session_list(self):
        global fit, ses_list
        
        
        if len(ses_list) == 0:
            self.comboBox_select_ses.clear()
            self.comboBox_select_ses.addItem("session 1") 
        

        elif len(ses_list) != 0:
            self.comboBox_select_ses.clear()
            for i in range(len(ses_list)):
                self.comboBox_select_ses.addItem('session %s'%(i+1),i)
                     
        #self.select_session(0)

    def select_session(self, index):
        global fit, ses_list

        ind = self.comboBox_select_ses.itemData(index) 
        #print(ind,index,len(ses_list))
        if ind == None:
            fit = ses_list[0]
        else:
            fit = ses_list[ind]

        #ses_list[ind-1] = fit

        self.init_fit()
       
        
        self.update_use_from_input_file()   
        self.update_use()
        self.update_gui_params()
        self.update_params()
        self.update_RV_file_buttons() 
        
        if not ind == None:    
            ses_list[ind] = fit

#######################################################################            


    def new_session(self):
        global fit
        
        file_pi = open('.sessions/empty.ses', 'rb')
        fit_new = dill.load(file_pi)
        file_pi.close()     
 
        ses_list.append(fit_new)
        self.session_list()
            
    def open_session(self):
        global fit,ses_list
        
        input_file = QtGui.QFileDialog.getOpenFileName(self, 'Open session', '', 'Data (*.ses)')

        file_pi = open(input_file[0], 'rb')
        fit_new = dill.load(file_pi)
        file_pi.close()     
        ses_list.append(fit_new)
        self.session_list()

    def save_session(self):
        global fit
        
        output_file = QtGui.QFileDialog.getSaveFileName(self, 'Save session', '', 'Data (*.ses)')

        file_pi = open(output_file[0], 'wb')
        dill.dump(fit, file_pi)
        file_pi.close()


    def open_sessions(self):
        global fit, ses_list
        
        input_file = QtGui.QFileDialog.getOpenFileName(self, 'Open session', '', 'Data (*.mses)')

        file_pi = open(input_file[0], 'rb')
        fit2 = dill.load(file_pi)
        file_pi.close()   
        if len(ses_list) == 0: 
            ses_list = fit2
        else:  
            ses_list = ses_list + fit2

        self.session_list()
        self.select_session(0)

    def save_sessions(self):
        global fit, ses_list
        
        output_file = QtGui.QFileDialog.getSaveFileName(self, 'Save multi-session', '', 'Data (*.mses)')

        file_pi = open(output_file[0], 'wb')
        dill.dump(ses_list, file_pi)
        file_pi.close()



    def quit(self):
        #os.system("rm temp*.vels")
        choice = QtGui.QMessageBox.information(self, 'Warning!',
                                            "Do you want to save the session before you Quit?",
                                            QtGui.QMessageBox.Cancel | QtGui.QMessageBox.No | QtGui.QMessageBox.Yes)  
         
        if choice == QtGui.QMessageBox.No:
            self.close()
        elif choice == QtGui.QMessageBox.Yes:
            self.save_session()
        elif choice == QtGui.QMessageBox.Cancel:
            return 
 
            
    def keyPressEvent(self, event):
        global fit
        if event.key() in (QtCore.Qt.Key_Enter, QtCore.Qt.Key_Return):
            self.update_use()
            self.update_params() 
            self.init_fit()
            return
       # super(Settings, self).keyPressEvent(event) 

################################## MCMC #######################################

    def worker_mcmc_complete(self):
        global fit  
        fit.print_info(short_errors=False)
        self.statusBar().showMessage('') 

    def worker_mcmc(self):
        global fit  
        
        self.button_MCMC.setEnabled(False)
        self.statusBar().showMessage('MCMC in progress....')        

        # check if RV data is present
        if fit.filelist.ndset <= 0:
             choice = QtGui.QMessageBox.information(self, 'Warning!',
             "Not possible to run MCMC if there are no data loaded. Please add your RV or transit data first. Okay?", QtGui.QMessageBox.Ok)      
             self.button_MCMC.setEnabled(True)         
             return        
        


        choice = QtGui.QMessageBox.information(self, 'Warning!',
                                            "This will run in the background and may take some time. Results are printed in the 'Stdout/Stderr' tab. Okay?",
                                            QtGui.QMessageBox.Ok) 
         
        # Pass the function to execute
        worker = Worker(lambda: self.run_mcmc()) # Any other args, kwargs are passed to the run  
        # Execute
        worker.signals.finished.connect(self.worker_mcmc_complete)
        
        # worker.signals.result.connect(self.print_output)
        #worker.signals.finished.connect(self.thread_complete)
       # worker.signals.progress.connect(self.progress_fn)
        
        
        self.threadpool.start(worker)
        
    def run_mcmc(self):
        global fit

        
        
        gp_params = [self.GP_rot_kernel_Amp.value(),
                     self.GP_rot_kernel_time_sc.value(),
                     self.GP_rot_kernel_Per.value(),
                     self.GP_rot_kernel_fact.value()]

        use_gp_params = [self.use_GP_rot_kernel_Amp.isChecked(),
                     self.use_GP_rot_kernel_time_sc.isChecked(),
                     self.use_GP_rot_kernel_Per.isChecked(),
                     self.use_GP_rot_kernel_fact.isChecked()]
 
        fit.mcmc(doGP=self.goGP.isChecked(), gp_par=np.array(gp_params),use_gp_par=np.array(use_gp_params), 
        burning_ph=self.burning_phase.value(), mcmc_ph=self.mcmc_phase.value(), threads=int(self.N_threads.value()), output=True)  
 
        self.button_MCMC.setEnabled(True)
  

    def grab_screen(self):
        p = QtWidgets.QWidget.grab(self)
        filename = QtGui.QFileDialog.getSaveFileName(self, 'Save image', '', '')
        p.save(filename[0], 'jpg')
        #label.setPixmap(p)        # just for fun :)
        
################################## END MCMC ################################### 

 
    def __init__(self):
        global fit

        QtWidgets.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
       # self.showMaximized()
       # self.setGeometry(3,30,800,800)
       # self.setFixedSize(1024,1024) 
       
        self.setupUi(self)
        
        self.initialize_buttons()
        self.initialize_plots()    
#        self.init_fit()
        self.console_widget = ConsoleWidget_embed(font_size = 10)
        self.terminal_embeded.addTab(self.console_widget, "Jupyter")
        
        if sys.platform[0:5] == "linux":
            self.terminal_embeded.addTab(terminal.EmbTerminal(), "Bash shell")        
        self.terminal_embeded.addTab(pg_console.ConsoleWidget(), "pqg shell")  
 

        self.gridLayout_text_editor.addWidget(ted.MainWindow())       
        self.gridLayout_calculator.addWidget(calc.Calculator())  
        
        
        if sys.version_info[0] == 2:
            self.gridLayout_stdout.addWidget(stdout_pipe.MyDialog())  
       
        self.dialog = print_info(self)
        self.dialog_credits = print_info(self)
       

        self.load_fort_in_file.clicked.connect(self.showDialog_fortran_input_file)

        self.buttonGroup_4.buttonClicked.connect(self.showDialog_RV_input_file)
        self.buttonGroup_remove_RV_data.buttonClicked.connect(self.remove_RV_file)
 
        self.buttonGroup_activity_data.buttonClicked.connect(self.showDialog_act_input_file)
        self.buttonGroup_remove_activity_data.buttonClicked.connect(self.remove_act_file)     
        
        self.buttonGroup_transit_data.buttonClicked.connect(self.showDialog_tra_input_file)
        self.buttonGroup_remove_transit_data.buttonClicked.connect(self.remove_tra_file)         
        
        
        self.buttonGroup_use.buttonClicked.connect(self.update_use)


        self.run_batman_test()

        
        self.button_orb_evol.clicked.connect(self.worker_Nbody) 
        self.button_MCMC.clicked.connect(self.worker_mcmc)
        self.button_Bootstrap.clicked.connect(lambda: self.run_bootstrap())
        
    
        
        ########## RV fitting ########################
        
        self.button_init_fit.clicked.connect(lambda: self.worker_RV_fitting(0))
        self.button_fit.clicked.connect(lambda: self.worker_RV_fitting())        
        self.button_auto_fit.clicked.connect(lambda: self.run_auto_fit())
        self.minimize_1param()


        ############ Sessions #################
        
        self.actionNew_session.triggered.connect(self.new_session)
        self.actionOpen_session.triggered.connect(self.open_session)
        self.actionSave_session.triggered.connect(self.save_session)
        self.actionSave_multi_sesion.triggered.connect(self.save_sessions)
        self.actionOpen_multi_session.triggered.connect(self.open_sessions) 
        #self.comboBox_extra_plot.activated.connect(self.change_extra_plot)      
        self.comboBox_select_ses.activated.connect(self.select_session)
        self.session_list()
        self.new_ses.clicked.connect(self.getNewses)
        self.copy_ses.clicked.connect(self.cop_ses)
        self.remove_ses.clicked.connect(self.rem_ses)
  
    
    
        self.actiongrab_screen.triggered.connect(self.grab_screen) 
        self.actionvisit_TRIFON_on_GitHub.triggered.connect(lambda: webbrowser.open('https://github.com/3fon3fonov/trifon'))    
        self.actionCredits.triggered.connect(lambda: self.print_info_credits())
        
        
        self.init_correlations_combo()
  
  
        self.init_activity_combo()
        
        self.radioButton_act_GLS_period.toggled.connect(lambda: self.update_activity_gls_plots(self.comboBox_act_data_gls.currentIndex()))
       
        self.comboBox_act_data_gls.activated.connect(lambda: self.update_activity_gls_plots(self.comboBox_act_data_gls.currentIndex())) 
        self.comboBox_act_data.activated.connect(lambda: self.update_activity_data_plots(self.comboBox_act_data.currentIndex())) 
       
        self.radioButton_RV_o_c_GLS_period.toggled.connect(self.update_RV_o_c_GLS_plots)
        self.radioButton_RV_GLS_period.toggled.connect(self.update_RV_GLS_plots)

        self.quit_button.clicked.connect(self.quit)
 
        self.jupiter_push_vars()
        
        self.threadpool = QtCore.QThreadPool()

        print("Hi there! Here you can get some more information from the tool's workflow, stdout/strerr, and the mcmc and bootstrap results.")

#Function Main START
def main():
    app = QtWidgets.QApplication(sys.argv)
    window = MyApp()
    window.show()
    sys.exit(app.exec_())
#Function Main END


if __name__ == '__main__':
    main() 











