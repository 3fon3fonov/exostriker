# -*- coding: utf-8 -*-
import pyqtgraph as pg
from PyQt6 import QtCore, QtGui, QtWidgets, uic

import numpy as np
import os, sys
from print_info_window import print_info
from worker import Worker
from multiprocessing import cpu_count
import gls as gls
import dill
import RV_mod as rv
import pg_hack


#qtCreatorFile = "./lib/UI/tdt.ui" 
#Ui_DetrendWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)
#import time
#start_time = time.time()

try:
    from tdt import Ui_Detrend as Ui_DetrendWindow
except (ImportError, KeyError,ModuleNotFoundError) as e:
    qtCreatorFile = "./lib/UI/tdt.ui" #%lib_path 
    Ui_DetrendWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)
#print("--- %s seconds ---" % (time.time() - start_time))

try:
    import sklearn
    sklearn_found = True
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




class DetrendWindow(QtWidgets.QWidget, Ui_DetrendWindow):
    def __init__(self,parent):
        super(DetrendWindow, self).__init__()

        QtWidgets.QWidget.__init__(self)
        Ui_DetrendWindow.__init__(self)
       # self.setWindowTitle('Transit detrending options')
        self.font = QtGui.QFont()
        self.font.setPointSize(9)
        self.font.setBold(False)
        self.threadpool_DT = QtCore.QThreadPool()
        self.threadpool_DT.setMaxThreadCount(cpu_count())    
        
        self.parent=parent
        # Create the main window
        self.ui = Ui_DetrendWindow()
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
        self.airmass_store       = {k: [] for k in range(20)}


        self.t = []
        self.old_t = []
        self.flux_o_c = []
        self.airmass = []
        
        self.initialize_plots()

        self.init_comboBox_regres()
        self.init_comboBox_sliders()
        self.init_comboBox_poly()
        self.init_comboBox_splines()
        self.init_comboBox_GP()

        self.ui.try_button.clicked.connect(self.worker_detrend)
        self.ui.saveProduct.clicked.connect(self.save_data_product)
        self.ui.readme_button.clicked.connect(self.info)
        self.ui.print_stat.clicked.connect(self.print_stat_window)

        self.info_dialog = print_info(self)
        self.stat_dialog = print_info(self)


        self.ui.buttonGroup_plot2.buttonClicked.connect(self.replot)

        self.ui.buttonGroup_trendOptions.buttonClicked.connect(self.update_labels)

        self.ui.click_to_reject.clicked.connect(self.top_plot)
        self.ui.addROI.clicked.connect(self.top_plot)

        self.ui.reset_data.clicked.connect(self.reset_data)
        self.ui.add_epoch.clicked.connect(self.add_bjd)

        self.ui.removeROIdata_in.clicked.connect(lambda: self.removeROI_data(inside=True))
        self.ui.removeROIdata_out.clicked.connect(lambda: self.removeROI_data(inside=False))

        self.ui.button_bin_data.clicked.connect(self.bin_data)        
        
        self.ui.apply_dilution.clicked.connect(self.add_dilution)

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

        self.t      = self.parent.tra_data[0]
        self.flux   = self.parent.tra_data[4]
        self.airmass   = self.parent.tra_data[3]
        
        self.flux_err = self.parent.tra_data[2]
        self.data_file_name = self.parent.tra_data[-1]
        self.old_t = dill.copy(self.t)
        return


    def add_dilution(self):

        D_flux = self.flux/(self.ui.Dilution_fact.value())
        self.flux      = D_flux - np.median(D_flux) * (1.0 - self.ui.Dilution_fact.value())

        D_flux_err = self.flux_err/(self.ui.Dilution_fact.value())
        self.flux_err      = D_flux_err

        self.ui.radio_remove_median.setChecked(True)
        #self.plot()
        self.worker_detrend()

 
        return


    def add_bjd(self):

        self.t      = self.t + self.ui.extra_BJD.value()
        self.ui.radio_remove_median.setChecked(True)
        #self.plot()
        self.worker_detrend()
 
        return

    def bin_data(self):
 
        self.ui.radio_remove_median.setChecked(True)


        self.ui.try_button.setEnabled(False)
        #self.ui.label_working.setText("Working!!!")
        self.ui.try_button.setText("Working!!!")        
        #self.ui.try_button.setText("Working!!!")
        #print('test 1')
        t_, flux_, flux_err_, ind = rv.bin_data(self.t,self.flux,self.flux_err, np.zeros(len(self.t)), bin_size =self.ui.bin_data.value())         
        #print('test 2')
        self.ui.try_button.setEnabled(True)
        self.ui.try_button.setText("Try !")

        self.t = t_         
        self.flux = flux_
        self.flux_err = flux_err_
        
        self.worker_detrend()
 
        return



    def reset_data(self):
 
        self.ui.radio_remove_median.setChecked(True)
        self.t = []
        self.worker_detrend()
        
        return


    def calculate(self):

        
        if self.ui.radio_remove_median.isChecked():

            flatten_lc1 = self.flux/np.median(self.flux)
            trend_lc1 = np.ones(len(self.flux))*np.median(self.flux)

        elif self.ui.radio_remove_mean.isChecked():

            flatten_lc1 = self.flux/np.mean(self.flux)
            trend_lc1 = np.ones(len(self.flux))*np.mean(self.flux)

        elif self.ui.radio_timeW.isChecked():
            flatten_lc1, trend_lc1 = flatten(
                self.t,                 # Array of time values
                self.flux ,                 # Array of flux values
                method=str(self.ui.comboBox_sliders.currentText()),
                window_length=self.ui.sliders_wl.value(),    # The length of the filter window in units of ``time``
#                break_tolerance=self.ui.spline_bt.value(),  # Split into segments at breaks longer than that
                return_trend=True,    # Return trend and flattened light curve
                )

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


        self.flux_o_c = flatten_lc1
        self.trend = trend_lc1
        self.flux_err_o_c = self.flux_err/trend_lc1
        



    def worker_detrend_complete(self):

        #self.ui.label_working.setText("")
        self.ui.try_button.setText("Try !")
        self.ui.try_button.setEnabled(True)
        self.ui.flatten_data.setChecked(True)
        self.old_t = dill.copy(self.t)


        self.plot()
        self.show()

        return
        
    def worker_detrend(self):
        
        if len(self.t) == 0:
            self.init_data()

        self.ui.try_button.setEnabled(False)
        #self.ui.label_working.setText("Working!!!")
        self.ui.try_button.setText("Working!!!")
        worker_detrend_wk = Worker(self.calculate)# Any other args, kwargs are passed to the run  
        worker_detrend_wk.signals.finished.connect(self.worker_detrend_complete)

        # worker.signals.result.connect(self.print_output)
        #worker.signals.finished.connect(self.thread_complete)
       # worker.signals.progress.connect(self.progress_fn)
        self.threadpool_DT.start(worker_detrend_wk)



    def plot(self):

        self.top_plot()
        self.bottom_plot_lc()


    def top_plot(self):
        global p_1
        
        
        if self.ui.addROI.isChecked():
            self.ui.removeROIdata_in.setEnabled(True)
            self.ui.removeROIdata_out.setEnabled(True)
        else:
            self.ui.removeROIdata_in.setEnabled(False)
            self.ui.removeROIdata_out.setEnabled(False)                
    
        if self.ui.click_to_reject.isChecked() or self.ui.addROI.isChecked():
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


        if self.ui.addROI.isChecked() == True:

            range_plot = self.ui.plot.viewRange()

            self.roi = pg.ROI([range_plot[0][0],range_plot[1][0]],[(range_plot[0][1]-range_plot[0][0])/10.0,(range_plot[1][1]-range_plot[1][0])/2.0],pen=pg.mkPen(color='r', width=3))
            ## handles scaling horizontally around center
            self.roi.addScaleHandle([1, 0.5], [0.5, 0.5])
            self.roi.addScaleHandle([0, 0.5], [0.5, 0.5])

            ## handles scaling vertically from opposite edge
            self.roi.addScaleHandle([0.5, 0], [0.5, 1])
            self.roi.addScaleHandle([0.5, 1], [0.5, 0])

            ## handles scaling both vertically and horizontally
            self.roi.addScaleHandle([1, 1], [0, 0])
            self.roi.addScaleHandle([0, 0], [1, 1])
            self.ui.plot.addItem(self.roi)
 
            self.ui.plot.plotItem.items[4].sigRegionChanged.connect(self.getROIbox)
            self.getROIbox(self.roi)

        self.ui.plot.plotItem.items[1].sigPointsClicked.connect(self.plotClicked)


    def getROIbox(self,roi):
        self.roi_xwindow = (roi.pos()[0] , roi.pos()[0] + roi.size()[0])
        self.roi_ywindow = (roi.pos()[1] , roi.pos()[1] + roi.size()[1])



    def bottom_plot_lc(self):

        
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
        



    def removeROI_data(self, inside=True):
 
        self.old_t = dill.copy(self.t)
        self.old_f = dill.copy(self.flux)

        if inside==True:
            subsetter = np.where((self.roi_xwindow[0] > self.old_t) | (self.old_t > self.roi_xwindow[1]) | (self.roi_ywindow[0] > self.old_f) | (self.old_f > self.roi_ywindow[1])  )  
        else:
            subsetter = np.where((self.roi_xwindow[0] < self.old_t) & (self.old_t < self.roi_xwindow[1]) & (self.roi_ywindow[0] < self.old_f) & (self.old_f < self.roi_ywindow[1])  )       
 

        self.t            = dill.copy(self.t[subsetter])
        self.flux         = dill.copy(self.flux[subsetter])
        self.flux_err     = dill.copy(self.flux_err[subsetter])
        self.flux_o_c     = dill.copy(self.flux_o_c[subsetter])
        self.flux_err_o_c = dill.copy(self.flux_err_o_c[subsetter])
        self.trend        = dill.copy(self.trend[subsetter])
        self.airmass      = dill.copy(self.airmass[subsetter])

        self.ui.plot.plotItem.items[1].setData(x=self.t, y=self.flux)
        self.ui.plot.plotItem.items[2].setData(x=self.t, y=self.flux,  
                                   top=self.flux_err, 
                                   bottom=self.flux_err)

        
        if self.ui.flatten_data.isChecked():
            self.ui.plot_2.plotItem.items[1].setData(x=self.t, y=self.flux_o_c)
            self.ui.plot_2.plotItem.items[2].setData(x=self.t, y=self.flux_o_c,  
                                       top=self.flux_err_o_c, 
                                       bottom=self.flux_err_o_c)
        else:
            self.replot()
 
   

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
        self.airmass      = dill.copy(self.airmass[self.old_t != rem_x])

        self.ui.plot.plotItem.items[1].setData(x=self.t, y=self.flux)
        self.ui.plot.plotItem.items[2].setData(x=self.t, y=self.flux,  
                                   top=self.flux_err, 
                                   bottom=self.flux_err)

        
        if self.ui.flatten_data.isChecked():
            self.ui.plot_2.plotItem.items[1].setData(x=self.t,y=self.flux_o_c)
            self.ui.plot_2.plotItem.items[2].setData(x=self.t, y=self.flux_o_c,  
                                       top=self.flux_err_o_c, 
                                       bottom=self.flux_err_o_c)
        else:
            self.replot()
        #self.plot()



    def make_GLS(self, model=False, o_c=False):

        #omega = 1/ np.logspace(np.log10(self.parent.gls_min_period.value()), np.log10(self.parent.gls_max_period.value()), num=int(self.parent.gls_n_omega.value()))        
        #omega = 1/ np.logspace(np.log10(0.9), np.log10((max(self.t)-min(self.t))*2.0), num=int(self.parent.gls_n_omega.value()))
        ind_norm = self.parent.gls_norm_combo.currentIndex()

        if model == False and o_c == False:
            data_for_GLS   = self.flux
            e_data_for_GLS = self.flux_err
        elif model == False and  o_c == True:
            data_for_GLS   = self.flux_o_c
            e_data_for_GLS = self.flux_err_o_c
        else:
            data_for_GLS   = self.trend
            e_data_for_GLS = self.flux_err

        self.trend_per = gls.Gls((self.t, data_for_GLS, e_data_for_GLS), 
          #  fast=True,  verbose=False, norm= "ZK",ofac=self.parent.gls_ofac.value(), fbeg=omega[-1], fend=omega[ 0],)
            fast=True,  verbose=False, norm=self.parent.norms[ind_norm],ofac=self.parent.gls_ofac.value(), fbeg=1/self.parent.gls_max_period.value(), fend=1/self.parent.gls_min_period.value())            


    def plot_GLS(self):
        #global p_2

        self.ui.plot_2.plot(clear=True,)

        power_levels = np.array([self.parent.gls_fap1.value(),self.parent.gls_fap2.value(),self.parent.gls_fap3.value()])

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
                
                self.t_store[self.parent.tra_data_index]             = self.t
                self.flux_store[self.parent.tra_data_index]          = self.flux
                self.flux_err_store[self.parent.tra_data_index]      = self.flux_err
                self.flux_o_c_store[self.parent.tra_data_index]      = self.flux_o_c
                self.flux_err_o_c_store[self.parent.tra_data_index]  = self.flux_err_o_c
                self.trend_store[self.parent.tra_data_index]         = self.trend

                self.airmass_store[self.parent.tra_data_index]       = self.airmass
 
                
                self.ui.radio_remove_median.setChecked(True)
                #QtWidgets.QMainWindow.closeEvent(self, event)
                self.close()
            else:
                event.ignore()


    def save_data_product(self):

        output_file = QtWidgets.QFileDialog.getSaveFileName(self, 'Save detrended data file', 'detrended_%s'%self.data_file_name, 'All (*.*);;Data (*.tran)', options=QtWidgets.QFileDialog.DontUseNativeDialog)
        
        if str(output_file[0]) != '':
            f = open(output_file[0], 'w')
            f.write("# BJD     Detrended data    Detrended data errors  Airmass   Original data   Original data errors    Model applied \n")
            for i in range(len(self.t)):
                f.write('{0:{width}.{precision}f}  {1:{width}.{precision}f}  {2:{width}.{precision}f} {3:{width}.{precision}f} {4:{width}.{precision}f} {5:{width}.{precision}f} {6:{width}.{precision}f}\n'.format(
                        float(self.t[i]), 
                        float(self.flux_o_c[i]), 
                        float(self.flux_err_o_c[i]), 
                        float(self.airmass[i]),
                        float(self.flux[i]), 
                        float(self.flux_err[i]), 
                        float(self.trend[i]), 
                        width = 14, precision = 7 ))
            f.close()

    def print_stat_window(self):

        self.stat_dialog.setFixedSize(550, 600)
        self.stat_dialog.setWindowTitle('Detrending stat. info')
 

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
        self.info_dialog.setWindowTitle('Detrending options info')
 
    
        text = ''
        self.info_dialog.text.setText(text) 
        
        text = "For more info on the detrending algorithms see <a href='https://github.com/hippke/wotan'>wotan</a>" 
        self.info_dialog.text.append(text)

        text = """
<br>
<br>
As explained in "wotan", some algorithms request
additional dependencies, which are not included in "wotan", and thus, not included in the Exo-Striker 
dependencies list. For example:
<br> 
<br> "huber", "ramsay", and "hampel" depend on "statsmodels"
<br> "Gaussian processes", "hspline", "ridge", and "lasso"  depend on "sklearn"
<br> "pspline" depends on "pygam"
<br> "supersmoother" depends on "supersmoother"
<br> 
<br> To install all additional dependencies, try to install these python packages:
<br> 
<br> * pip install statsmodels 
<br> * pip install sklearn  
<br> * pip install supersmoother
<br> * pip install pygam
<br>  
<br>  Also, "wotan" depends on "numba" and "llvmlite" so it might be a good idea to update those two if you already have older versions.
<br>  
<br>
<br> If you made the use of the detrending options for your paper, please also cite: <a href='https://ui.adsabs.harvard.edu/abs/2019AJ....158..143H/abstract'> Hippke et al. (2019)</a>
"""
        self.info_dialog.text.append(text)


    
        self.info_dialog.text.setReadOnly(True)
        #self.dialog.setWindowIcon (QtGui.QIcon('logo.png'))
        self.info_dialog.show()


    def update_labels(self):
        
        if self.ui.radio_GPs.isChecked():
            self.ui.label_method.setText("Kernel")
            self.ui.label_wl.setText("Kernel size")
            self.ui.label_tolerance.setText("Kernel period")
        else:
            self.ui.label_method.setText("Method")
            self.ui.label_wl.setText("Window length")
            self.ui.label_tolerance.setText("break tolerance")



    def init_comboBox_sliders(self):

        sliders     = ["biweight","huber","huber_psi","hampel","andrewsinewave","welsch","ramsay","tau","hodges","median",
"medfilt","mean","trim_mean","winsorize","hampelfilt"] 
        sliders_use = [self.wotan_found, self.statsmodels_found, self.wotan_found, self.statsmodels_found,
                       self.wotan_found, self.wotan_found, self.statsmodels_found,self.wotan_found,self.wotan_found,
                       self.wotan_found, self.wotan_found, self.wotan_found, self.wotan_found,self.wotan_found,self.wotan_found] 
        
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
        regres_use = [self.wotan_found,self.supersmoother_found,self.sklearn_found,self.sklearn_found]
        
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

