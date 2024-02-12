# -*- coding: utf-8 -*-
import pyqtgraph as pg
from PyQt6 import QtCore, QtGui, QtWidgets, uic

import numpy as np
import os,sys
from print_info_window import print_info
from worker import Worker
from multiprocessing import cpu_count
import gls as gls
import dill
import RV_mod as rv
from scipy import signal
import pg_hack

#from wotan import flatten


#qtCreatorFile = "./lib/UI/tdt.ui" 
#Ui_DetrendWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)
#import time
#start_time = time.time()

try:
    from rvs import Ui_Rvs_data as Ui_Rvs_dataWindow
except (ImportError, KeyError,ModuleNotFoundError) as e:
    qtCreatorFile = "./lib/UI/rvs.ui" #%lib_path 
    Ui_Rvs_dataWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)
#print("--- %s seconds ---" % (time.time() - start_time))

try:
    import sklearn
    sklearn_found = False
except (ImportError, KeyError,ModuleNotFoundError) as e:
    sklearn_found = False
    pass

try:
    import statsmodels
    statsmodels_found = True
except (ImportError, KeyError,ModuleNotFoundError) as e:
    statsmodels_found = False
    pass

try:
    import pygam
    pygam_found = True
except (ImportError, KeyError,ModuleNotFoundError) as e:
    pygam_found = False
    pass

try:
    import supersmoother
    supersmoother_found = True
except (ImportError, KeyError,ModuleNotFoundError) as e:
    supersmoother_found = False
    pass

try:
    from wotan import flatten
    wotan_found = True
except (ImportError, KeyError,ModuleNotFoundError) as e:
    wotan_found = False
    print("wotan not found!")
    pass

class RvsWindow(QtWidgets.QWidget, Ui_Rvs_dataWindow):
    def __init__(self,parent):
        super(Ui_Rvs_dataWindow, self).__init__()

        QtWidgets.QWidget.__init__(self)
        Ui_Rvs_dataWindow.__init__(self)
       # self.setWindowTitle('Transit detrending options')
        self.font = QtGui.QFont()
        self.font.setPointSize(9)
        self.font.setBold(False)
        self.threadpool_DT = QtCore.QThreadPool()
        self.threadpool_DT.setMaxThreadCount(cpu_count())    
        
        self.parent=parent
        # Create the main window
        self.ui = Ui_Rvs_dataWindow()
        self.ui.setupUi(self)
        self.setWindowIcon(QtGui.QIcon('./lib/UI/33_striker.png'))

        self.sklearn_found = sklearn_found
        self.statsmodels_found = statsmodels_found
        self.pygam_found = pygam_found
        self.supersmoother_found = supersmoother_found
        self.wotan_found = wotan_found
      
        self.ui.radio_GPs.setEnabled(sklearn_found)
        self.ui.comboBox_GP.setEnabled(sklearn_found)
        self.ui.kernel_size.setEnabled(sklearn_found)
        
        if wotan_found == False:
            self.ui.radio_timeW.setEnabled(wotan_found)
            self.ui.radio_Splines.setEnabled(wotan_found)
            self.ui.radio_Polynomials.setEnabled(wotan_found)
            self.ui.radio_Regressions.setEnabled(wotan_found)
            self.ui.radio_GPs.setEnabled(wotan_found)


        self.t_store             = {k: [] for k in range(20)}
        self.flux_store          = {k: [] for k in range(20)}
        self.flux_err_store      = {k: [] for k in range(20)}
        self.flux_o_c_store      = {k: [] for k in range(20)}
        self.flux_err_o_c_store  = {k: [] for k in range(20)}
        self.trend_store         = {k: [] for k in range(20)}
        self.flux_idset_store    = {k: [] for k in range(20)}

        self.t = []
        self.old_t = []
        self.flux_o_c = []
        self.initialize_plots()

        self.init_comboBox_regres()
        self.init_comboBox_sliders()
        self.init_comboBox_poly()
        self.init_comboBox_splines()
        self.init_comboBox_GP()

        self.ui.try_button.clicked.connect(self.worker_options)
        self.ui.saveProduct.clicked.connect(self.save_data_product)
        self.ui.readme_button.clicked.connect(self.info)
        self.ui.print_stat.clicked.connect(self.print_stat_window)
        
        

        self.info_dialog = print_info(self)
        self.stat_dialog = print_info(self)


        self.ui.buttonGroup_plot2.buttonClicked.connect(self.replot)

        self.ui.buttonGroup_trendOptions.buttonClicked.connect(self.update_labels)

        self.ui.click_to_reject.clicked.connect(self.top_plot)
        self.ui.reset_data.clicked.connect(self.reset_data)
        self.ui.add_epoch.clicked.connect(self.add_bjd)

        self.ui.button_bin_data.clicked.connect(self.bin_data)        
        self.ui.button_sigma_clip.clicked.connect(self.sigma_clip_func)

    def replot(self):
        
        if self.ui.GLS_of_data.isChecked():
            self.make_GLS()
            self.plot_GLS()
        elif self.ui.GLS_of_model.isChecked():
            self.make_GLS(model=True)
            self.plot_GLS()
        elif self.ui.GLS_of_detr_data.isChecked():
            self.make_GLS(model=False,o_c=True)
            self.plot_GLS()
        elif self.ui.flatten_data.isChecked():
            self.bottom_plot_lc()


    def init_data(self):

 
        self.t      = self.parent.rvs_data[0]
        self.flux   = self.parent.rvs_data[1]
        self.flux_err = self.parent.rvs_data[2]
        self.flux_idset = self.parent.rvs_data[3]
        self.data_file_name = self.parent.rvs_data[4]
        self.old_t = dill.copy(self.t)
        return


    def smooth(self, y, box_pts):
        box = np.ones(box_pts)/box_pts
        y_smooth = np.convolve(y, box, mode='same')
        return y_smooth 


    def add_bjd(self):

        self.t      = self.t + self.ui.extra_BJD.value()
        self.ui.radio_no_remove.setChecked(True)
        #self.plot()
        self.worker_options()
 
        return

    def bin_data(self):
 
        self.ui.radio_no_remove.setChecked(True)
        
        self.ui.try_button.setText("Working!!!")
        t_, flux_, flux_err_, ind = rv.bin_data(self.t,self.flux,self.flux_err, self.flux_idset, bin_size =self.ui.bin_data.value())
        self.ui.try_button.setText("Try !")

        self.t = t_         
        self.flux = flux_
        self.flux_err = flux_err_
        self.flux_idset = ind
        
        self.worker_options()
 
        return

    def sigma_clip_func(self):
 
        self.ui.radio_no_remove.setChecked(True)
        
        self.ui.try_button.setText("Working!!!")
        t_, flux_, flux_err_, ind = rv.sigma_clip_rvs(self.t,self.flux,self.flux_err, self.flux_idset, sigma_clip =self.ui.rvs_sigma_clip.value())
        self.ui.try_button.setText("Try !")

        self.t = t_         
        self.flux = flux_
        self.flux_err = flux_err_
        self.flux_idset = ind
        
        self.worker_options()
 
        return



    def reset_data(self):
 
        self.ui.radio_no_remove.setChecked(True)
        self.t = []
        self.worker_options()
        
        return


    def calculate(self):

        if self.ui.radio_no_remove.isChecked():

            flatten_lc1 = 0 #np.median(self.flux)
            trend_lc1 = np.ones(len(self.flux)) 
        
        elif self.ui.radio_remove_median.isChecked():

            flatten_lc1 = np.median(self.flux)
            trend_lc1 = np.ones(len(self.flux))*np.median(self.flux)

        elif self.ui.radio_remove_mean.isChecked():

            flatten_lc1 = np.mean(self.flux)
            trend_lc1 = np.ones(len(self.flux))*np.mean(self.flux)

        elif self.ui.radio_timeW.isChecked():

            time_gls = dill.copy(self.t)
            
            #smoothed_  = smooth(self.flux,30)
        
            #fs_ = int(max(time_gls)-min(time_gls))
           # nyq = 0.5 * (1./(fs_*86400.0))
            #2.0*np.pi ,86400.0
            
            if self.ui.comboBox_sliders.currentText() == "lowpass":
                filter_freq = (1 / (float(self.ui.filter_low_freq.value()) )) #/ nyq
            elif  self.ui.comboBox_sliders.currentText() == "highpass":
                filter_freq = (1 / (float(self.ui.filter_high_freq.value()) )) #/ nyq
            else:
                filter_freq = [(1 /  (float(self.ui.filter_low_freq.value()) ))# / nyq, 
                            (1 /  (float(self.ui.filter_high_freq.value()) ))  ]#/ nyq]
                
                
            
            #print(filter_freq)
            
            sos = signal.butter(int(self.ui.filter_order.value()), filter_freq, 
                                str(self.ui.comboBox_sliders.currentText()), analog = True, # fs=fs_, 
                                output='sos')
            
#fs=1000
            filtered = signal.sosfilt(sos, self.flux)
            flatten_lc1 = filtered
            trend_lc1 = filtered
            #flatten_lc1, trend_lc1 = flatten(
           #     self.t,                 # Array of time values
           #     self.flux ,                 # Array of flux values
           #     method=str(self.ui.comboBox_sliders.currentText()),
          #      window_length=self.ui.sliders_wl.value(),    # The length of the filter window in units of ``time``
#                break_tolerance=self.ui.spline_bt.value(),  # Split into segments at breaks longer than that
         #       return_trend=True,    # Return trend and flattened light curve
         #       )

        elif self.ui.radio_Splines.isChecked():
            flatten_lc1, trend_lc1 = flatten(
                self.t,                 # Array of time values
                self.flux ,                 # Array of flux values
                method=str(self.ui.comboBox_splines.currentText()),
                window_length=self.ui.spline_wl.value(),    # The length of the filter window in units of ``time``
                break_tolerance=self.ui.spline_bt.value(),  # Split into segments at breaks longer than that
                return_trend=True,    # Return trend and flattened light curve
                )

        elif self.ui.radio_Polynomials.isChecked():

            flatten_lc1, trend_lc1 = flatten(
                self.t,                 # Array of time values
                self.flux ,                 # Array of flux values
                method=str(self.ui.comboBox_poly.currentText()),
                window_length=self.ui.poly_wl.value(),    # The length of the filter window in units of ``time``
                break_tolerance=self.ui.poly_bt.value(),  # Split into segments at breaks longer than that
                return_trend=True,    # Return trend and flattened light curve
                )

        elif self.ui.radio_Regressions.isChecked():

            flatten_lc1, trend_lc1 = flatten(
                self.t,                 # Array of time values
                self.flux ,                 # Array of flux values
                method=str(self.ui.comboBox_regs.currentText()),
                window_length=self.ui.regres_wl.value(),    # The length of the filter window in units of ``time``
                break_tolerance=self.ui.regres_bt.value(),  # Split into segments at breaks longer than that
                return_trend=True,    # Return trend and flattened light curve
                )
            
        elif self.ui.radio_GPs.isChecked():

            flatten_lc1, trend_lc1 = flatten(
                self.t,                 # Array of time values
                self.flux ,                 # Array of flux values
                method='gp',
                kernel = str(self.ui.comboBox_GP.currentText()),
                kernel_size=self.ui.kernel_size.value(),
                break_tolerance=self.ui.regres_bt.value(),  # Split into segments at breaks longer than that
                kernel_period = self.ui.GP_period.value(),
                robust = self.ui.checkBox_GP_robust.isChecked(),
                return_trend=True    # Return trend and flattened light curve
                )

        else:
            flatten_lc1 = self.flux 
            trend_lc1 = np.ones(len(self.flux))*np.median(self.flux)


        self.flux_o_c = self.flux - flatten_lc1
        self.trend = trend_lc1
        self.flux_err_o_c = dill.copy(self.flux_err)# / trend_lc1
        
        #print(self.flux_o_c,self.flux_err_o_c)
        



    def worker_options_complete(self):

        #self.ui.label_working.setText("")
        self.ui.try_button.setText("Try !")
        self.ui.try_button.setEnabled(True)
        self.ui.flatten_data.setChecked(True)
        self.old_t = dill.copy(self.t)


        self.plot()
        self.show()

        return
        
    def worker_options(self):
        
        if len(self.t) == 0:
            self.init_data()

        self.ui.try_button.setEnabled(False)
        #self.ui.label_working.setText("Working!!!")
        self.ui.try_button.setText("Working!!!")
        worker_options_wk = Worker(self.calculate)# Any other args, kwargs are passed to the run  
        worker_options_wk.signals.finished.connect(self.worker_options_complete)

        # worker.signals.result.connect(self.print_output)
        #worker.signals.finished.connect(self.thread_complete)
       # worker.signals.progress.connect(self.progress_fn)
        self.threadpool_DT.start(worker_options_wk)



    def plot(self):

        self.top_plot()
        self.bottom_plot_lc()


    def top_plot(self):
        global p_1
        
        if self.ui.click_to_reject.isChecked():
            symbolSize=6
        else:
            symbolSize=2

        self.ui.plot.plot(clear=True,)

        ######## Top plot ############

        self.ui.plot.plot(self.t,self.flux, pen=None,
            symbol='o', symbolPen={'color': '#0066ff', 'width': 1.1},
            symbolSize=symbolSize,enableAutoRange=True,viewRect=True,
            symbolBrush='#0066ff')
            
        err_ = pg.ErrorBarItem(x=self.t, y=self.flux, symbol = 'o',
                                   top=self.flux_err, 
                                   bottom=self.flux_err,
                                   beam=0.0, pen='#0066ff')

        self.ui.plot.addItem(err_)

        model_curve = self.ui.plot.plot(self.t, self.trend , pen={'color': '#000000', 'width': 3}, enableAutoRange=True,viewRect=True ) 
        model_curve.setZValue(1)

        self.ui.plot.plotItem.items[1].sigPointsClicked.connect(self.plotClicked)


    def bottom_plot_lc(self):
        #global p_2
        
        self.ui.plot_2.plot(clear=True,)
        ######## Bottom plot ############
        self.ui.plot_2.setLogMode(False,False)

        self.ui.plot_2.plot(self.t,self.flux_o_c, pen=None,
            symbol='o', symbolPen={'color': '#0066ff', 'width': 1.1},
            symbolSize=2,enableAutoRange=True,viewRect=True,
            symbolBrush='#0066ff')
            
        self.ui.plot_2.setLabel('left', 'Flux', units='',  **{'font-size':'9pt'})


        err_ = pg.ErrorBarItem(x=self.t, y=self.flux_o_c, symbol = 'o',
                                   top=self.flux_err_o_c, 
                                   bottom=self.flux_err_o_c,
                                   beam=0.0, pen='#0066ff')
                                   
                                
                                   
        self.ui.plot_2.addItem(err_)
        
        
        

    def plotClicked(self,curve,datas):

        if self.ui.click_to_reject.isChecked() == False:
            return

        rem_x,rem_y = datas[0].pos()
        print("Removed x,y: ",rem_x,rem_y)

        self.old_t = dill.copy(self.t)

        self.t            = dill.copy(self.t[self.old_t != rem_x])
        self.flux         = dill.copy(self.flux[self.old_t != rem_x])
        self.flux_err     = dill.copy(self.flux_err[self.old_t != rem_x])
        self.flux_o_c     = dill.copy(self.flux_o_c[self.old_t != rem_x])
        self.flux_err_o_c = dill.copy(self.flux_err_o_c[self.old_t != rem_x])
        self.trend        = dill.copy(self.trend[self.old_t != rem_x])
        self.flux_idset   = dill.copy(self.flux_idset[self.old_t != rem_x])

        self.ui.plot.plotItem.items[1].setData(x=self.t,y=self.flux)
        self.ui.plot.plotItem.items[2].setData(x=self.t, y=self.flux,  
                                   top=self.flux_err, 
                                   bottom=self.flux_err)

        self.ui.plot_2.plotItem.items[1].setData(x=self.t, y=self.flux_o_c)
        self.ui.plot_2.plotItem.items[2].setData(x=self.t, y=self.flux_o_c,  
                                   top=self.flux_err_o_c, 
                                   bottom=self.flux_err_o_c)
        #self.plot()



    def make_GLS(self, model=False, o_c=False):

        #omega = 1/ np.logspace(np.log10(self.parent.gls_min_period.value()), np.log10(self.parent.gls_max_period.value()), num=int(self.parent.gls_n_omega.value()))
        time_gls = dill.copy(self.t)
        ind_norm = self.parent.gls_norm_combo.currentIndex()        
        #omega = 1/ np.logspace(np.log10(0.9), np.log10((max(time_gls)-min(time_gls))*2.0), num=int(self.parent.gls_n_omega.value()))
        
        
        if model == False and o_c == False:
            data_for_GLS   = dill.copy(self.flux)
            e_data_for_GLS = dill.copy(self.flux_err)
        elif model == False and  o_c == True:
            data_for_GLS   = dill.copy(self.flux_o_c)
            e_data_for_GLS = dill.copy(self.flux_err_o_c)
        else:
            data_for_GLS   = dill.copy(self.trend)
            e_data_for_GLS = dill.copy(self.flux_err)

        self.trend_per = gls.Gls((time_gls, data_for_GLS, e_data_for_GLS), 
           # fast=True,  verbose=False, norm= "ZK",ofac=self.parent.gls_ofac.value(), fbeg=omega[-1], fend=omega[ 0],)
            fast=True,  verbose=False, norm=self.parent.norms[ind_norm],ofac=self.parent.gls_ofac.value(), fbeg=1/self.parent.gls_max_period.value(), fend=1/self.parent.gls_min_period.value())            


    def plot_GLS(self):
        #global p_2

        self.ui.plot_2.plot(clear=True,)

        power_levels = np.array([self.parent.gls_fap1.value(),
                                 self.parent.gls_fap2.value(),
                                 self.parent.gls_fap3.value()])

        ######################## GLS ##############################
        if self.parent.radioButton_act_GLS_period.isChecked():
            self.ui.plot_2.setLogMode(True,False)
            self.ui.plot_2.plot(1.0/self.trend_per.freq, self.trend_per.power,pen='r',symbol=None ) 
            self.ui.plot_2.setLabel('bottom', 'period [d]', units='',  **{'font-size':'9pt'}) 
            self.ui.plot_2.setLabel('left', 'Power', units='',  **{'font-size':'9pt'})

        else:
            self.ui.plot_2.setLogMode(False,False)        
            self.ui.plot_2.plot(self.trend_per.freq, self.trend_per.power,pen='r',symbol=None )
            self.ui.plot_2.setLabel('bottom', 'frequency [1/d]', units='',  **{'font-size':'9pt'}) 

        [self.ui.plot_2.addLine(x=None, y=fap, pen=pg.mkPen('k', width=0.8, style=QtCore.Qt.DotLine)) for ii,fap in enumerate(self.trend_per.powerLevel(np.array(power_levels)))]




    def closeEvent(self, event):
        
        if len(self.old_t)  != len(self.t):
            choice = QtWidgets.QMessageBox.information(self, 'Warning!',
            "It seems that you removed data, but you did not refit! This is not allowed. Please press the 'Try!' button and then close", QtWidgets.QMessageBox.Ok)
            event.ignore()

 
    
        else:
            ret = QtWidgets.QMessageBox.question(None, 'Close request', 'Are you sure you want to quit?',
                                             QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                                             QtWidgets.QMessageBox.Yes)
            if ret == QtWidgets.QMessageBox.Yes:
                
                self.t_store[self.parent.rvs_data_index]             = self.t
                self.flux_store[self.parent.rvs_data_index]          = self.flux
                self.flux_err_store[self.parent.rvs_data_index]      = self.flux_err
                self.flux_o_c_store[self.parent.rvs_data_index]      = self.flux_o_c
                self.flux_err_o_c_store[self.parent.rvs_data_index]  = self.flux_err_o_c
                self.trend_store[self.parent.rvs_data_index]         = self.trend
                self.flux_idset_store[self.parent.rvs_data_index]         = self.flux_idset
                

                self.ui.radio_no_remove.setChecked(True)
                #QtWidgets.QMainWindow.closeEvent(self, event)
                self.close()
            else:
                event.ignore()


    def save_data_product(self):

        output_file = QtWidgets.QFileDialog.getSaveFileName(self, 'Save modified data file', 'rvs_%s'%self.data_file_name, 'All (*.*);;Data (*.tran)', options=QtWidgets.QFileDialog.DontUseNativeDialog)
        
        if str(output_file[0]) != '':
            f = open(output_file[0], 'w')
            f.write("# BJD     RVS data    RVS data errors,   Original data   Original data errors    Model applied \n")
            for i in range(len(self.t)):
                f.write('{0:{width}.{precision}f}  {1:{width}.{precision}f}  {2:{width}.{precision}f} {3:{width}.{precision}f} {4:{width}.{precision}f} {5:{width}.{precision}f}\n'.format(
                        float(self.t[i]), 
                        float(self.flux_o_c[i]), 
                        float(self.flux_err_o_c[i]), 
                        float(self.flux[i]), 
                        float(self.flux_err[i]), 
                        float(self.trend[i]), 
                        width = 14, precision = 7 ))
            f.close()

    def print_stat_window(self):

        self.stat_dialog.setFixedSize(550, 600)
        self.stat_dialog.setWindowTitle('Stat. info')
 

        text_info = """ 
"""
        self.stat_dialog.text.setText(text_info) 

        ################## text generator #################
        text_info = """ 
-----------------------------------  

N data          :  %s

first epoch     :  %.3f
last epoch      :  %.3f
time span       :  %.3f

min. value      :  %.4f
max. value      :  %.4f
end-to-end      :  %.4f
mean            :  %.4f
median          :  %.4f
rms             :  %.4f

min error       :  %.4f
max error       :  %.4f
mean error      :  %.4f
median error    :  %.4f

"""%(len(self.t), 
self.t[0], 
self.t[-1], 
self.t[-1]-self.t[0], 
np.min(self.flux_o_c), 
np.max(self.flux_o_c), 
np.max(self.flux_o_c)-np.min(self.flux_o_c), 
np.mean(self.flux_o_c), 
np.median(self.flux_o_c), 
np.sqrt(np.mean(np.square(self.flux_o_c))),
np.min(self.flux_err_o_c), 
np.max(self.flux_err_o_c),   
np.mean(self.flux_err_o_c),  
np.median(self.flux_err_o_c))


        self.stat_dialog.text.append(text_info)

        self.stat_dialog.text.setReadOnly(True)
        #self.dialog.setWindowIcon (QtGui.QIcon('logo.png'))
        self.stat_dialog.show()



    def info(self):
        
        #self.info_dialog.setGeometry(300, 200, 150, 150)
        self.info_dialog.setFixedSize(550, 600)
        self.info_dialog.setWindowTitle('Options info')
 
    
        text = ''
        self.info_dialog.text.setText(text) 
        
        text = "" 
        self.info_dialog.text.append(text)

        text = """
<br>
<br>
To be done
<br> 
"""
        self.info_dialog.text.append(text)


    
        self.info_dialog.text.setReadOnly(True)
        #self.dialog.setWindowIcon (QtGui.QIcon('logo.png'))
        self.info_dialog.show()


    def update_labels(self):
        
        if self.ui.radio_GPs.isChecked():
            self.ui.label_method.setText("Kernel")
            self.ui.label_wl.setText("Kernel size")
            self.ui.label_high_freq.setText("Kernel period")
        elif self.ui.radio_timeW.isChecked():
            self.ui.label_method.setText("Method")
            self.ui.label_wl.setText("Order")
            self.ui.label_high_freq.setText("High freq")
            self.ui.label_low_freq.setText("Low freq")
        elif self.ui.radio_Splines.isChecked():
            self.ui.label_method.setText("Method")           
            self.ui.label_wl.setText("Window length")
            self.ui.label_high_freq.setText("break tolerance")
        else:
            self.ui.label_method.setText("Method")           
            self.ui.label_wl.setText("Window length")
            self.ui.label_high_freq.setText("break tolerance")
 
    def init_comboBox_sliders(self):

        sliders     = ["lowpass", "highpass", "bandpass", "bandstop"] 
        sliders_use = [True,True,True,True] 
        
        for i in range(len(sliders)):
            if sliders_use[i] == True:
                self.ui.comboBox_sliders.addItem(sliders[i],i) 

    def init_comboBox_poly(self):

        poly = ["cofiam","cosine","savgol"]
        for i in range(len(poly)):
            self.ui.comboBox_poly.addItem(poly[i],i) 

    def init_comboBox_splines(self):

        splines     = ["rspline","hspline","pspline"]
        splines_use = [True,self.sklearn_found,self.pygam_found]

        for i in range(len(splines)):
            if splines_use[i] == True:
                self.ui.comboBox_splines.addItem(splines[i],i) 

    def init_comboBox_regres(self):

        regres = ["lowess","supersmoother","ridge","lasso"]
        regres_use = [True,self.supersmoother_found,self.sklearn_found,self.sklearn_found]
        
        for i in range(len(regres)):
            if regres_use[i] == True:
                self.ui.comboBox_regs.addItem(regres[i],i) 

    def init_comboBox_GP(self):

        gps = ["squared_exp","matern","periodic","periodic_auto"]
        for i in range(len(gps)):
            self.ui.comboBox_GP.addItem(gps[i],i) 


    def initialize_plots(self):

        #global p_1,p_2
        
        xaxis = ['BJD [days]','BJD [days]']
        yaxis = ['Flux','Flux']
        xunit = ['' ,'']
        yunit = ['' ,'' ]

        #p_1 = self.ui.plot
        #p_2 = self.ui.plot_2

        zzz = [self.ui.plot,self.ui.plot_2]


        for i in range(len(zzz)):
            
            zzz[i].setAxisItems({'bottom': pg_hack.CustomAxisItem('bottom')})

            #zzz[i].getAxis("bottom").tickFont = self.font
            zzz[i].getAxis("bottom").setStyle(tickTextOffset = 12, tickFont = self.font)
            #zzz[i].getAxis("top").tickFont = self.font
            zzz[i].getAxis("top").setStyle(tickTextOffset = 12, tickFont = self.font)
            #zzz[i].getAxis("left").tickFont = self.font
            zzz[i].getAxis("left").setStyle(tickTextOffset = 12, tickFont = self.font)
            #zzz[i].getAxis("right").tickFont = self.font
            zzz[i].getAxis("right").setStyle(tickTextOffset = 12, tickFont = self.font)
            zzz[i].getAxis('left').setWidth(50)
            zzz[i].getAxis('right').setWidth(10)
            zzz[i].getAxis('top').setHeight(10)
            zzz[i].getAxis('bottom').setHeight(50)

            zzz[i].setLabel('bottom', '%s'%xaxis[i], units='%s'%xunit[i],  **{'font-size':'9pt'})
            zzz[i].setLabel('left',   '%s'%yaxis[i], units='%s'%yunit[i],  **{'font-size':'9pt'})       
            zzz[i].showAxis('top') 
            zzz[i].showAxis('right') 
            zzz[i].getAxis('bottom').enableAutoSIPrefix(enable=False)

        #zzz[i].getViewBox().setAspectLocked(True)

        return



def main():
    app = QtWidgets.QApplication(sys.argv)
#    main = DetrendWindow()
#    main.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()

