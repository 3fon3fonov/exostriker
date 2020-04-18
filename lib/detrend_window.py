# -*- coding: utf-8 -*-
import pyqtgraph as pg
from PyQt5 import QtCore, QtGui, QtWidgets, uic

import numpy as np
import os
from print_info_window import print_info
from worker import Worker
from multiprocessing import cpu_count
import gls as gls

from wotan import flatten

qtCreatorFile = "./lib/UI/tdt.ui" 
Ui_DetrendWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)


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

        self.sklearn_found = True
        self.statsmodels_found = True
        self.pygam_found = True
        self.supersmoother_found = True

        try:
            import sklearn
        except (ImportError, KeyError) as e:
            self.ui.radio_GPs.setEnabled(False)
            self.ui.comboBox_GP.setEnabled(False)
            self.ui.kernel_size.setEnabled(False)
            self.sklearn_found = False
            pass
        
        try:
            import statsmodels
        except (ImportError, KeyError) as e:
            self.statsmodels_found = False
            pass
        
        try:
            import pygam
        except (ImportError, KeyError) as e:
            self.pygam_found = False
            pass
        
        try:
            import supersmoother
        except (ImportError, KeyError) as e:
            self.supersmoother_found = False
            pass

        self.initialize_plots()

        self.init_comboBox_regres()
        self.init_comboBox_sliders()
        self.init_comboBox_poly()
        self.init_comboBox_splines()
        self.init_comboBox_GP()

        self.ui.try_button.clicked.connect(self.worker_detrend)
        self.ui.saveProduct.clicked.connect(self.save_data_product)
        self.ui.readme_button.clicked.connect(self.info)

        self.info_dialog = print_info(self)


        self.ui.buttonGroup_plot2.buttonClicked.connect(self.replot)

        self.ui.buttonGroup_trendOptions.buttonClicked.connect(self.update_labels)



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

 
    def calculate(self):

        
        #self.ui.label_working.setText("Working!")
        
        self.t      = self.parent.tra_data[0]
        self.flux   = self.parent.tra_data[4]
        self.flux_err = self.parent.tra_data[2]
        self.data_file_name = self.parent.tra_data[5]


        if self.ui.radio_remove_mean.isChecked():

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
                return_trend=True,    # Return trend and flattened light curve
                )

        else:
            flatten_lc1 = self.flux 
            trend_lc1 = np.ones(len(self.flux))*np.mean(self.flux)


        self.flux_o_c = flatten_lc1
        self.trend = trend_lc1
        self.flux_err_o_c = self.flux_err/trend_lc1
        



    def worker_detrend_complete(self):

        #self.ui.label_working.setText("")
        self.ui.try_button.setText("Try !")
        self.ui.try_button.setEnabled(True)
        self.ui.flatten_data.setChecked(True)

        self.plot()
        self.show()

        return
        
    def worker_detrend(self):

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
        
        p_1.plot(clear=True,)

        ######## Top plot ############

        p_1.plot(self.t,self.flux, pen=None,
            symbol='o', symbolPen={'color': '#0066ff', 'width': 1.1},
            symbolSize=2,enableAutoRange=True,viewRect=True,
            symbolBrush='#0066ff')
            
        err_ = pg.ErrorBarItem(x=self.t, y=self.flux, symbol = 'o',
                                   top=self.flux_err, 
                                   bottom=self.flux_err,
                                   beam=0.0, pen='#0066ff')

        p_1.addItem(err_)

        model_curve = p_1.plot(self.t, self.trend , pen={'color': '#000000', 'width': 3}, enableAutoRange=True,viewRect=True ) 
        model_curve.setZValue(1)



    def bottom_plot_lc(self):
        global p_2
        
        p_2.plot(clear=True,)
        ######## Bottom plot ############
        p_2.setLogMode(False,False)

        p_2.plot(self.t,self.flux_o_c, pen=None,
            symbol='o', symbolPen={'color': '#0066ff', 'width': 1.1},
            symbolSize=2,enableAutoRange=True,viewRect=True,
            symbolBrush='#0066ff')
            
        p_2.setLabel('left', 'Flux', units='',  **{'font-size':'9pt'})


        err_ = pg.ErrorBarItem(x=self.t, y=self.flux_o_c, symbol = 'o',
                                   top=self.flux_err_o_c, 
                                   bottom=self.flux_err_o_c,
                                   beam=0.0, pen='#0066ff')
                                   
        p_2.addItem(err_)
        
        
        

    def make_GLS(self, model=False, o_c=False):

        #omega = 1/ np.logspace(np.log10(self.parent.gls_min_period.value()), np.log10(self.parent.gls_max_period.value()), num=int(self.parent.gls_n_omega.value()))
        
        omega = 1/ np.logspace(np.log10(0.9), np.log10((max(self.t)-min(self.t))*2.0), num=int(self.parent.gls_n_omega.value()))
        
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
            fast=True,  verbose=False, norm= "ZK",ofac=self.parent.gls_ofac.value(), fbeg=omega[-1], fend=omega[ 0],)



    def plot_GLS(self):
        global p_2

        p_2.plot(clear=True,)

        power_levels = np.array([self.parent.gls_fap1.value(),self.parent.gls_fap2.value(),self.parent.gls_fap3.value()])

        ######################## GLS ##############################
        if self.parent.radioButton_act_GLS_period.isChecked():
            p_2.setLogMode(True,False)
            p_2.plot(1.0/self.trend_per.freq, self.trend_per.power,pen='r',symbol=None ) 
            p_2.setLabel('bottom', 'period [d]', units='',  **{'font-size':'9pt'}) 
            p_2.setLabel('left', 'Power', units='',  **{'font-size':'9pt'})

        else:
            p_2.setLogMode(False,False)        
            p_2.plot(self.trend_per.freq, self.trend_per.power,pen='r',symbol=None )
            p_2.setLabel('bottom', 'frequency [1/d]', units='',  **{'font-size':'9pt'}) 

        [p_2.addLine(x=None, y=fap, pen=pg.mkPen('k', width=0.8, style=QtCore.Qt.DotLine)) for ii,fap in enumerate(self.trend_per.powerLevel(np.array(power_levels)))]




    def closeEvent(self, event):
        ret = QtGui.QMessageBox.question(None, 'Close request', 'Are you sure you want to quit?',
                                         QtGui.QMessageBox.Yes | QtGui.QMessageBox.No,
                                         QtGui.QMessageBox.Yes)
        if ret == QtGui.QMessageBox.Yes:
            self.ui.radio_remove_mean.setChecked(True)
            QtGui.QMainWindow.closeEvent(self, event)
        else:
            event.ignore()


    def save_data_product(self):

        output_file = QtGui.QFileDialog.getSaveFileName(self, 'Save detrended data file', 'detrended_%s'%self.data_file_name, 'All (*.*);;Data (*.tran)', options=QtGui.QFileDialog.DontUseNativeDialog)
        
        if str(output_file[0]) != '':
            f = open(output_file[0], 'w')
            for i in range(len(self.t)):
                f.write('{0:{width}.{precision}f}  {1:{width}.{precision}f}  {2:{width}.{precision}f}  \n'.format(float(self.t[i]), float(self.flux_o_c[i]), float(self.flux_err_o_c[i]),  width = 10, precision = 7 )   )
            f.close()



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
            self.ui.label_method.setText("kernel")
            self.ui.label_wl.setText("kernel size")
        else:
            self.ui.label_method.setText("method")
            self.ui.label_wl.setText("window length")



    def init_comboBox_sliders(self):

        sliders     = ["biweight","huber","huber_psi","hampel","andrewsinewave","welsch","ramsay","tau","hodges","median",
"medfilt","mean","trim_mean","winsorize","hampelfilt"] 
        sliders_use = [True, self.statsmodels_found, True, self.statsmodels_found,True, True, self.statsmodels_found,True,True,True,True,True,True,True,True] 
        
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

        global p_1,p_2
        
        xaxis = ['BJD [days]','BJD [days]']
        yaxis = ['Flux','Flux']
        xunit = ['' ,'']
        yunit = ['' ,'' ]

        p_1 = self.ui.plot
        p_2 = self.ui.plot_2

        zzz = [p_1,p_2]


        for i in range(len(zzz)):

            zzz[i].getAxis("bottom").tickFont = self.font
            zzz[i].getAxis("bottom").setStyle(tickTextOffset = 12)
            zzz[i].getAxis("top").tickFont = self.font
            zzz[i].getAxis("top").setStyle(tickTextOffset = 12)
            zzz[i].getAxis("left").tickFont = self.font
            zzz[i].getAxis("left").setStyle(tickTextOffset = 12)
            zzz[i].getAxis("right").tickFont = self.font
            zzz[i].getAxis("right").setStyle(tickTextOffset = 12)
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
