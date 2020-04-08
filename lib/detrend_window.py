# -*- coding: utf-8 -*-
import pyqtgraph as pg
from PyQt5 import QtCore, QtGui, QtWidgets, uic

import numpy as np
import os
from print_info_window import print_info

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
        
        self.parent=parent
        # Create the main window
        self.ui = Ui_DetrendWindow()
        self.ui.setupUi(self)
        self.setWindowIcon(QtGui.QIcon('./lib/UI/33_striker.png'))

        self.initialize_plots()

        self.init_comboBox_regres()
        self.init_comboBox_sliders()
        self.init_comboBox_poly()
        self.init_comboBox_splines()
        self.init_comboBox_GP()

        self.ui.try_button.clicked.connect(self.plot)
        self.ui.readme_button.clicked.connect(self.info)

        self.info_dialog = print_info(self)


        self.ui.buttonGroup_trendOptions.buttonClicked.connect(self.update_labels)



    def calculate(self):

        self.t      = self.parent.tra_data[0]
        self.flux   = self.parent.tra_data[4]
        self.flux_err = self.parent.tra_data[2]


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


 
    def plot(self):

        self.calculate()

        ######## Top plot ############

        self.ui.plot.plot(self.t,self.flux, clear=True, pen=None,
            symbol='o', symbolPen={'color': '#0066ff', 'width': 1.1},
            symbolSize=2,enableAutoRange=True,viewRect=True,
            symbolBrush='#0066ff')
            
        err_ = pg.ErrorBarItem(x=self.t, y=self.flux, symbol = 'o',
                                   top=self.flux_err, 
                                   bottom=self.flux_err,
                                   beam=0.0, pen='#0066ff')
                                   
        self.ui.plot.addItem(err_)

        model_curve = self.ui.plot.plot(self.t, self.trend , pen={'color': '#000000', 'width': 3}, enableAutoRange=True,viewRect=True ) 
        model_curve.setZValue(1)


        ######## Bottom plot ############

        self.ui.plot_2.plot(self.t,self.flux_o_c, clear=True, pen=None,
            symbol='o', symbolPen={'color': '#0066ff', 'width': 1.1},
            symbolSize=2,enableAutoRange=True,viewRect=True,
            symbolBrush='#0066ff')

        err_ = pg.ErrorBarItem(x=self.t, y=self.flux_o_c, symbol = 'o',
                                   top=self.flux_err_o_c, 
                                   bottom=self.flux_err_o_c,
                                   beam=0.0, pen='#0066ff')
                                   
        self.ui.plot_2.addItem(err_)
        self.show()


    def info(self):
        
        self.info_dialog.setGeometry(300, 300, 150, 150)
        self.info_dialog.setWindowTitle('Detrending options info')  
 
    
        text = ''
        self.info_dialog.text.setText(text) 
        
        text = "For more info on the detrending algorithms see <a href='https://github.com/hippke/wotan'>wotan</a>" 
        self.info_dialog.text.append(text)


        text = "\n"*3 + """If you made the use of the detrending options for your paper, please also cite Hippke et al. (2019, AJ, 158, 143)
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

        sliders = ["biweight","huber","huber_psi","hampel","andrewsinewave","welsch","ramsay","tau","hodges","median",
"medfilt","mean","trim_mean","winsorize","hampelfilt"] 
        for i in range(len(sliders)):
            self.ui.comboBox_sliders.addItem(sliders[i],i) 

    def init_comboBox_poly(self):

        poly = ["cofiam","cosine","savgol"]
        for i in range(len(poly)):
            self.ui.comboBox_poly.addItem(poly[i],i) 

    def init_comboBox_splines(self):

        splines = ["rspline","hspline"]#,"pspline"]
        for i in range(len(splines)):
            self.ui.comboBox_splines.addItem(splines[i],i) 

    def init_comboBox_regres(self):

        regres = ["lowess","supersmoother","ridge","lasso"]
        for i in range(len(regres)):
            self.ui.comboBox_regs.addItem(regres[i],i) 

    def init_comboBox_GP(self):

        gps = ["squared_exp","matern","periodic","periodic_auto"]
        for i in range(len(gps)):
            self.ui.comboBox_GP.addItem(gps[i],i) 


    def initialize_plots(self):

        xaxis = ['BJD [days]','BJD [days]']
        yaxis = ['Flux','Flux']
        xunit = ['' ,'']
        yunit = ['' ,'' ]

        zzz = [self.ui.plot,self.ui.plot_2]

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

