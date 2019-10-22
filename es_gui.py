#!/usr/bin/python
__author__ = 'Trifon Trifonov'

from pathos.multiprocessing import freeze_support
freeze_support()

import numpy as np
#import matplotlib as mpl
#mpl.use('Qt5Agg')
import sys, os, traceback 
from PyQt5 import QtCore, QtGui, QtWidgets, uic


sys.path.insert(0, './lib')

import RV_mod as rv

import pyqtgraph as pg
import pyqtgraph.console as pg_console

import word_processor_es as text_editor_es
import calculator as calc 
import gls as gls 
from worker import Worker #, WorkerSignals

#from multiprocessing import cpu_count
#import time

#import BKR as bkr
from doublespinbox import DoubleSpinBox
#from Jupyter_emb import ConsoleWidget_embed
from Jupyter_emb import ConsoleWidget_embed

from stdout_pipe import MyDialog
from print_info_window import print_info
from symbols_window import show_symbols

import terminal
from tree_view import Widget_tree

import ntpath

from scipy.signal import argrelextrema
from scipy.stats.stats import pearsonr   
import scipy.stats as stat

#import batman as batman

try:
    from transitleastsquares import transitleastsquares    
    tls_not_found = False 
except (ImportError, KeyError) as e:
    tls_not_found = True
    pass               
  
try:
    import batman as batman   
    
    try: 
        bat_test = batman.TransitParams()
        batman_not_found = False 
        bat_test = 0     
    except (ImportError, KeyError, AttributeError) as e:     
        batman_not_found = True 
        
except (ImportError, KeyError) as e:
    batman_not_found = True
    pass       

    
import webbrowser
 
#try:
#    import cPickle as pickle
#except ModuleNotFoundError:
#    import pickle
import dill


#if hasattr(QtCore.Qt, 'AA_EnableHighDpiScaling'):
#    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)

if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps,True)



qtCreatorFile = "./lib/UI/es.ui" 
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

pg.setConfigOption('background', '#ffffff')
pg.setConfigOption('foreground', 'k')
 
#pg.setConfigOptions(useOpenGL=True) 
pg.setConfigOptions(antialias=True)  


global fit, colors, ses_list
 

arguments = len(sys.argv) - 1

if arguments==2 and sys.argv[1] == '-ses' and os.path.exists(sys.argv[2]):
    try:
        file_pi = open(sys.argv[2], 'rb')
        fit_ses = dill.load(file_pi)
        file_pi.close()   
        fit = fit_ses 
        ses_list = [fit_ses] 
        fit.init_pl_arb()
#        rv.check_temp_RV_file(fit)
        
        start_arg_ses = True  
    except (ImportError, KeyError, AttributeError) as e:
        print("You have entered non-RVmod session. %s cannot be recognaized"%sys.argv[2])
        fit=rv.signal_fit(name='session')
        ses_list = [fit]            
        start_arg_ses = False    
        
elif arguments==2 and sys.argv[1] == '-mses' and os.path.exists(sys.argv[2]):
    try:
        file_pi = open(sys.argv[2], 'rb')
        fit_ses = dill.load(file_pi)
        file_pi.close()   
        ses_list = fit_ses
        fit = ses_list[0]
        fit.init_pl_arb()
#        rv.check_temp_RV_file(fit)
        start_arg_ses = True  
    except (ImportError, KeyError, TypeError, AttributeError) as e:
        print("You have entered non-RVmod multi-session. %s cannot be recognaized"%sys.argv[2])
        fit=rv.signal_fit(name='session')
        ses_list = [fit]            
        start_arg_ses = False          
        
elif arguments==2 and sys.argv[1] == '-rv_init' and os.path.exists(sys.argv[2]):
    try:
        
        fit=rv.signal_fit(str(sys.argv[2]), 'RVmod session',readinputfile=True)
        fit.init_pl_arb()
        ses_list = [fit]      
        start_arg_ses = True  
    except (ImportError, KeyError, TypeError, AttributeError) as e:
        print("You have entered non-RVmod .init file. %s cannot be recognaized"%sys.argv[2])
        fit=rv.signal_fit(name='session')
        ses_list = [fit]            
        start_arg_ses = False             
        
  
else:    
    fit=rv.signal_fit(name='session')
    ses_list = [fit]            
    start_arg_ses = False         
     
 


colors      = ['#0066ff',  '#ff0000','#66ff66','#00ffff','#cc33ff','#ff9900','#cccc00','#3399ff','#990033','#339933','#666699']
colors_gls  = ['#0066ff',  '#ff0000','#66ff66','#00ffff','#cc33ff','#ff9900','#cccc00','#3399ff','#990033','#339933','#666699']
         
symbols = ['o','t','t1','t2','t3','s','p','h','star','+','d'] 
colors_delta_om = colors
colors_theta = colors


QtGui.QApplication.processEvents()



 
class Exo_striker(QtWidgets.QMainWindow, Ui_MainWindow):

    def update_labels(self):
        global fit

 
        self.value_stellar_mass.setText("%.4f"%(fit.params.stellar_mass))
        self.value_epoch.setText(str(fit.epoch))
        self.value_rms.setText("%.4f"%(fit.fit_results.rms))
        self.value_chi2.setText("%.4f"%(fit.fit_results.chi2)) 
        self.value_reduced_chi2.setText("%.4f"%(fit.fit_results.reduced_chi2))        
        #self.value_loglik.setText("%.4f"%(fit.fit_results.loglik)) 
        self.value_loglik.setText("%.4f"%(fit.loglik)) 
       
        self.value_Ndata.setText("%s"%(len(fit.fit_results.jd))) 
        self.value_DOF.setText("%s"%(len(fit.fit_results.jd) - fit.fit_results.mfit))        

        if fit.mod_dynamical == True:
            self.radioButton_Dynamical.setChecked(True)        
        else:
            self.radioButton_Keplerian.setChecked(True)       
            
        if fit.type_fit["RV"] == True and fit.type_fit["Transit"] == False:
            self.radioButton_RV.setChecked(True)        
        elif fit.type_fit["RV"] == False and fit.type_fit["Transit"] == True:           
            self.radioButton_transit.setChecked(True)                    
        elif fit.type_fit["RV"] == True and fit.type_fit["Transit"] == True:          
            self.radioButton_transit_RV.setChecked(True)    
            
        if fit.hkl == True:
            self.radioButton_hkl.setChecked(True)    
        else:
            self.radioButton_ewm.setChecked(True)    
            
           
            
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
            
        param_gui_wd = [self.om_dot_1, self.om_dot_2, self.om_dot_3, 
                        self.om_dot_4, self.om_dot_5, self.om_dot_6, 
                        self.om_dot_7, self.om_dot_8, self.om_dot_9]

        for i in range(9):
            param_gui_wd[i].setValue(fit.omega_dot[i])        
           
            
            
        param_gui_tr = [
                     self.t0_1, self.pl_rad_1, self.a_sol_1,
                     self.t0_2, self.pl_rad_2, self.a_sol_2,
                     self.t0_3, self.pl_rad_3, self.a_sol_3,
                     self.t0_4, self.pl_rad_4, self.a_sol_4, 
                     self.t0_5, self.pl_rad_5, self.a_sol_5,
                     self.t0_6, self.pl_rad_6, self.a_sol_6,
                     self.t0_7, self.pl_rad_7, self.a_sol_7, 
                     self.t0_8, self.pl_rad_8, self.a_sol_8,
                     self.t0_9, self.pl_rad_9, self.a_sol_9,
                     ]
         
        for i in range(fit.npl):
            param_gui_tr[i*3].setValue(fit.t0[i])           
            param_gui_tr[i*3+1].setValue(fit.pl_rad[i]) 
            param_gui_tr[i*3+2].setValue(fit.pl_a[i]) 
            
        rvs_data_gui = [self.Data1,self.Data2,self.Data3,self.Data4,self.Data5,
                    self.Data6,self.Data7,self.Data8,self.Data9,self.Data10]
        rvs_data_jitter_gui = [self.jitter_Data1,self.jitter_Data2,self.jitter_Data3,self.jitter_Data4,self.jitter_Data5,
                           self.jitter_Data6,self.jitter_Data7,self.jitter_Data8,self.jitter_Data9,self.jitter_Data10]

        for i in range(10): 
            rvs_data_gui[i].setValue(fit.params.offsets[i]) 
            rvs_data_jitter_gui[i].setValue(fit.params.jitters[i])
        
        tra_data_gui = [self.trans_Data1,self.trans_Data2,self.trans_Data3,self.trans_Data4,self.trans_Data5,
                        self.trans_Data6,self.trans_Data7,self.trans_Data8,self.trans_Data9,self.trans_Data10]
        tra_data_jitter_gui = [self.jitter_trans_Data1,self.jitter_trans_Data2,self.jitter_trans_Data3,self.jitter_trans_Data4,self.jitter_trans_Data5,
                               self.jitter_trans_Data6,self.jitter_trans_Data7,self.jitter_trans_Data8,self.jitter_trans_Data9,self.jitter_trans_Data10]
        
            
        for i in range(10): 
            tra_data_gui[i].setValue(fit.tra_off[i]) 
            tra_data_jitter_gui[i].setValue(fit.tra_jitt[i])            
            
        gp_rot_params = [self.GP_rot_kernel_Amp,
                     self.GP_rot_kernel_time_sc,
                     self.GP_rot_kernel_Per,
                     self.GP_rot_kernel_fact]
        
        for i in range(len(gp_rot_params)):
            gp_rot_params[i].setValue(fit.GP_rot_params[i])
 
        gp_sho_params = [self.GP_sho_kernel_S,
                     self.GP_sho_kernel_Q,
                     self.GP_sho_kernel_omega]
        
        for i in range(len(gp_sho_params)):
            gp_sho_params[i].setValue(fit.GP_sho_params[i])    
            
            
        tra_gp_rot_params = [self.tra_GP_rot_kernel_Amp,
                     self.tra_GP_rot_kernel_time_sc,
                     self.tra_GP_rot_kernel_Per,
                     self.tra_GP_rot_kernel_fact]
        
        for i in range(len(tra_gp_rot_params)):
            tra_gp_rot_params[i].setValue(fit.tra_GP_rot_params[i])
 
        tra_gp_sho_params = [self.tra_GP_sho_kernel_S,
                     self.tra_GP_sho_kernel_Q,
                     self.tra_GP_sho_kernel_omega]
        
        for i in range(len(tra_gp_sho_params)):
            tra_gp_sho_params[i].setValue(fit.tra_GP_sho_params[i])                
            

            
        self.St_mass_input.setValue(fit.params.stellar_mass)  
        self.St_radius_input.setValue(fit.stellar_radius)  
        
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
            
            
        param_gui_wd = [self.om_dot_1, self.om_dot_2, self.om_dot_3, 
                        self.om_dot_4, self.om_dot_5, self.om_dot_6, 
                        self.om_dot_7, self.om_dot_8, self.om_dot_9]

        for i in range(9):
            fit.omega_dot[i] = param_gui_wd[i].value()             
        
        fit.hack_around_rv_params() 
         
        param_gui_tr = [self.t0_1, self.pl_rad_1, self.a_sol_1,
             self.t0_2, self.pl_rad_2, self.a_sol_2,
             self.t0_3, self.pl_rad_3, self.a_sol_3,
             self.t0_4, self.pl_rad_4, self.a_sol_4, 
             self.t0_5, self.pl_rad_5, self.a_sol_5,
             self.t0_6, self.pl_rad_6, self.a_sol_6,
             self.t0_7, self.pl_rad_7, self.a_sol_7, 
             self.t0_8, self.pl_rad_8, self.a_sol_8,
             self.t0_9, self.pl_rad_9, self.a_sol_9,
             ]
         
        for i in range(fit.npl):
            fit.t0[i]     = param_gui_tr[i*3].value()   
            fit.pl_rad[i] = param_gui_tr[i*3+1].value() 
            fit.pl_a[i]   = param_gui_tr[i*3+2].value() 

        rvs_data_gui = [self.Data1,self.Data2,self.Data3,self.Data4,self.Data5,
                    self.Data6,self.Data7,self.Data8,self.Data9,self.Data10]
        rvs_data_jitter_gui = [self.jitter_Data1,self.jitter_Data2,self.jitter_Data3,self.jitter_Data4,self.jitter_Data5,
                           self.jitter_Data6,self.jitter_Data7,self.jitter_Data8,self.jitter_Data9,self.jitter_Data10]

        for i in range(10): 
            fit.params.offsets[i] = rvs_data_gui[i].value() 
            fit.params.jitters[i] = rvs_data_jitter_gui[i].value()

        tra_data_gui = [self.trans_Data1,self.trans_Data2,self.trans_Data3,self.trans_Data4,self.trans_Data5,
                        self.trans_Data6,self.trans_Data7,self.trans_Data8,self.trans_Data9,self.trans_Data10]
        tra_data_jitter_gui = [self.jitter_trans_Data1,self.jitter_trans_Data2,self.jitter_trans_Data3,self.jitter_trans_Data4,self.jitter_trans_Data5,
                               self.jitter_trans_Data6,self.jitter_trans_Data7,self.jitter_trans_Data8,self.jitter_trans_Data9,self.jitter_trans_Data10]
 
        for i in range(10): 
            fit.tra_off[i]  = tra_data_gui[i].value() 
            fit.tra_jitt[i] = tra_data_jitter_gui[i].value() 
 
 
        self.read_RV_GP() 
        self.read_tra_GP()

        fit.params.stellar_mass = self.St_mass_input.value() 
        fit.params.linear_trend = self.RV_lin_trend.value()   

        fit.stellar_radius = self.St_radius_input.value()        
        
        if self.checkBox_first_RV_epoch.isChecked() and len(fit.fit_results.rv_model.jd) != 0:
            fit.epoch = min(fit.fit_results.rv_model.jd)
        else:
            fit.epoch =  self.Epoch.value()
       


    def read_tra_GP(self):
        global fit  

        tra_gp_rot_params = [self.tra_GP_rot_kernel_Amp,
                     self.tra_GP_rot_kernel_time_sc,
                     self.tra_GP_rot_kernel_Per,
                     self.tra_GP_rot_kernel_fact]
        
        for i in range(len(tra_gp_rot_params)):
            fit.tra_GP_rot_params[i] = tra_gp_rot_params[i].value()    
            
            
        tra_gp_sho_params = [self.tra_GP_sho_kernel_S,
                     self.tra_GP_sho_kernel_Q,
                     self.tra_GP_sho_kernel_omega]
        
        for i in range(len(tra_gp_sho_params)):
            fit.tra_GP_sho_params[i] = tra_gp_sho_params[i].value()   
            


    def read_RV_GP(self):
        global fit  
            
        gp_rot_params = [self.GP_rot_kernel_Amp,
                     self.GP_rot_kernel_time_sc,
                     self.GP_rot_kernel_Per,
                     self.GP_rot_kernel_fact]
        
        for i in range(len(gp_rot_params)):
            fit.GP_rot_params[i] = gp_rot_params[i].value()    
            
            
        gp_sho_params = [self.GP_sho_kernel_S,
                     self.GP_sho_kernel_Q,
                     self.GP_sho_kernel_omega]
        
        for i in range(len(gp_sho_params)):
            fit.GP_sho_params[i] = gp_sho_params[i].value()             
            
   
    
    def set_hkl(self):
        global fit  
        
    
        param_gui = [self.e1, self.om1, self.ma1,
                     self.e2, self.om2, self.ma2,
                     self.e3, self.om3, self.ma3,
                     self.e4, self.om4, self.ma4, 
                     self.e5, self.om5, self.ma5,
                     self.e6, self.om6, self.ma6,
                     self.e7, self.om7, self.ma7,  
                     self.e8, self.om8, self.ma8, 
                     self.e9, self.om9, self.ma9]    
    
    
        if self.radioButton_ewm.isChecked():
            
            fit.hkl = False
#            fit.hack_around_rv_params()  
            
            self.label_ecc.setText("e")
            self.label_omega.setText("<html><head/><body><p>&omega; [deg]</p></body></html>")
            self.label_Ma.setText("Ma [deg]")   
            self.label_ecc2.setText("e")
            self.label_omega2.setText("<html><head/><body><p>&omega; [deg]</p></body></html>")
            self.label_Ma2.setText("Ma [deg]")               
            self.label_ecc3.setText("e")
            self.label_omega3.setText("<html><head/><body><p>&omega; [deg]</p></body></html>")
            self.label_Ma3.setText("Ma [deg]")        

            
            for i in range(9):
                param_gui[i*3].setRange(0.0,1.0)
               #param_gui[i*3].singleStep(0.001)
                param_gui[i*3+1].setRange(0.0,360.0)
                param_gui[i*3+2].setRange(0.0,360.0)                
                #param_gui[i*3+1].singleStep(0.001)                
                param_gui[i*3].setValue(fit.e[i])             
                param_gui[i*3+1].setValue(fit.w[i])             
                param_gui[i*3+2].setValue(fit.M0[i])             

#                fit.params.update_e(i,fit.e[i]) # update e for a given planet
#                fit.params.update_w(i,fit.w[i]) # update w for a given planet
#                fit.params.update_M0(i,fit.M0[i]) # update w for a given planet
           
         
        elif self.radioButton_hkl.isChecked():
            
            fit.hkl = True
           
            self.label_ecc.setText("<html><head/><body><p>h=esin(&omega;)</p></body></html>")
            self.label_omega.setText("<html><head/><body><p>k=ecos(&omega;)</p></body></html>")
            self.label_Ma.setText("<html><head/><body><p>&lambda; [deg]</p></body></html>")      
            self.label_ecc2.setText("<html><head/><body><p>h=esin(&omega;)</p></body></html>")
            self.label_omega2.setText("<html><head/><body><p>k=ecos(&omega;)</p></body></html>")
            self.label_Ma2.setText("<html><head/><body><p>&lambda; [deg]</p></body></html>")    
            self.label_ecc3.setText("<html><head/><body><p>h=esin(&omega;)</p></body></html>")
            self.label_omega3.setText("<html><head/><body><p>k=ecos(&omega;)</p></body></html>")
            self.label_Ma3.setText("<html><head/><body><p>&lambda; [deg]</p></body></html>")    

  
            for i in range(9):
                param_gui[i*3].setRange(-1.0,1.0)
               # param_gui[i*3].singleStep(0.01)                
                param_gui[i*3+1].setRange(-1.0,1.0) 
                param_gui[i*3+2].setRange(0.0,360.0)  
                param_gui[i*3].setValue(fit.e_sinw[i])              
                param_gui[i*3+1].setValue(fit.e_cosw[i])            
                param_gui[i*3+2].setValue(fit.lamb[i])   

        self.update_params()
        self.update_gui_params() 
    
    
    
    def set_tra_ld(self):
        global fit   


        uni_ld_models = [self.use_uniform_ld_1,self.use_uniform_ld_2,self.use_uniform_ld_3,self.use_uniform_ld_4,self.use_uniform_ld_5,
                         self.use_uniform_ld_6,self.use_uniform_ld_7,self.use_uniform_ld_8,self.use_uniform_ld_9,self.use_uniform_ld_10]
        
        lin_ld_models = [self.use_linear_ld_1,self.use_linear_ld_2,self.use_linear_ld_3,self.use_linear_ld_4,self.use_linear_ld_5,
                         self.use_linear_ld_6,self.use_linear_ld_7,self.use_linear_ld_8,self.use_linear_ld_9,self.use_linear_ld_10]        

        quad_ld_models = [self.use_quadratic_ld_1,self.use_quadratic_ld_2,self.use_quadratic_ld_3,self.use_quadratic_ld_4,self.use_quadratic_ld_5,
                          self.use_quadratic_ld_6,self.use_quadratic_ld_7,self.use_quadratic_ld_8,self.use_quadratic_ld_9,self.use_quadratic_ld_10] 
      
        nonlin_ld_models = [self.use_nonlinear_ld_1,self.use_nonlinear_ld_2,self.use_nonlinear_ld_3,self.use_nonlinear_ld_4,self.use_nonlinear_ld_5,
                          self.use_nonlinear_ld_6,self.use_nonlinear_ld_7,self.use_nonlinear_ld_8,self.use_nonlinear_ld_9,self.use_nonlinear_ld_10] 

        lin_u1 = [self.u1_linear_1, self.u1_linear_2, self.u1_linear_3, self.u1_linear_4, self.u1_linear_5,
                  self.u1_linear_6, self.u1_linear_7, self.u1_linear_8, self.u1_linear_9, self.u1_linear_10]

        quad_u1 = [self.u1_quadratic_1, self.u1_quadratic_2, self.u1_quadratic_3, self.u1_quadratic_4, self.u1_quadratic_5,
                   self.u1_quadratic_6, self.u1_quadratic_7, self.u1_quadratic_8, self.u1_quadratic_9, self.u1_quadratic_10]
        quad_u2 = [self.u2_quadratic_1, self.u2_quadratic_2, self.u2_quadratic_3, self.u2_quadratic_4, self.u2_quadratic_5,
                   self.u2_quadratic_6, self.u2_quadratic_7, self.u2_quadratic_8, self.u2_quadratic_9, self.u2_quadratic_10]
        
        nonlin_u1 = [self.u1_nonlin_1, self.u1_nonlin_2, self.u1_nonlin_3, self.u1_nonlin_4, self.u1_nonlin_5,
                     self.u1_nonlin_6, self.u1_nonlin_7, self.u1_nonlin_8, self.u1_nonlin_9, self.u1_nonlin_10]
        nonlin_u2 = [self.u2_nonlin_1, self.u2_nonlin_2, self.u2_nonlin_3, self.u2_nonlin_4, self.u2_nonlin_5,
                     self.u2_nonlin_6, self.u2_nonlin_7, self.u2_nonlin_8, self.u2_nonlin_9, self.u2_nonlin_10]   
        nonlin_u3 = [self.u3_nonlin_1, self.u3_nonlin_2, self.u3_nonlin_3, self.u3_nonlin_4, self.u3_nonlin_5,
                     self.u3_nonlin_6, self.u3_nonlin_7, self.u3_nonlin_8, self.u3_nonlin_9, self.u3_nonlin_10]
        nonlin_u4 = [self.u4_nonlin_1, self.u4_nonlin_2, self.u4_nonlin_3, self.u4_nonlin_4, self.u4_nonlin_5,
                     self.u4_nonlin_6, self.u4_nonlin_7, self.u4_nonlin_8, self.u4_nonlin_9, self.u4_nonlin_10]
        
        
        for i in range(10):
            if uni_ld_models[i].isChecked():
                fit.ld_m[i] = "uniform"
                fit.ld_u[i] = []               
            elif lin_ld_models[i].isChecked():
                fit.ld_m[i] = "linear" 
                fit.ld_u[i] = [lin_u1[i].value()]                              
            elif quad_ld_models[i].isChecked():
                fit.ld_m[i] = "quadratic"        
                fit.ld_u[i] = [quad_u1[i].value(),quad_u2[i].value()]                                              
            elif nonlin_ld_models[i].isChecked():
                fit.ld_m[i] = "nonlinear"        
                fit.ld_u[i] = [nonlin_u1[i].value(),nonlin_u2[i].value(),nonlin_u3[i].value(),nonlin_u4[i].value()]                  
            else:
                fit.ld_m[i] = "quadratic"              
                fit.ld_u[i] = [quad_u1[i].value(),quad_u2[i].value()]                                              
     

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
            
            
        param_errors_gui_wd = [self.err_om_dot_1,self.err_om_dot_2,self.err_om_dot_3,
                               self.err_om_dot_4,self.err_om_dot_5,self.err_om_dot_6,
                               self.err_om_dot_7,self.err_om_dot_8,self.err_om_dot_9,                      
                            ]
        
        for i in range(9):
            param_errors_gui_wd[i].setText("+/- %.3f"%max(np.abs(fit.omega_dot_err[i])))            
            

        data_errors_gui        = [self.err_Data1,self.err_Data2,self.err_Data3,self.err_Data4,self.err_Data5,
                                  self.err_Data6,self.err_Data7,self.err_Data8,self.err_Data9,self.err_Data10]
        data_errors_jitter_gui = [self.err_jitter_Data1,self.err_jitter_Data2,self.err_jitter_Data3,self.err_jitter_Data4,
                                  self.err_jitter_Data5,self.err_jitter_Data6,self.err_jitter_Data7,self.err_jitter_Data8,
                                  self.err_jitter_Data9,self.err_jitter_Data10]
        

        for i in range(10):
            data_errors_gui[i].setText("+/- %.3f"%max(np.abs(fit.param_errors.offset_errors[i])))
            data_errors_jitter_gui[i].setText("+/- %.3f"%max(np.abs(fit.param_errors.jitter_errors[i])))
            
        tra_data_errors_gui        = [self.err_trans_Data1,self.err_trans_Data2,self.err_trans_Data3,self.err_trans_Data4,self.err_trans_Data5,
                                      self.err_trans_Data6,self.err_trans_Data7,self.err_trans_Data8,self.err_trans_Data9,self.err_trans_Data10]
        tra_data_errors_jitter_gui = [self.err_jitter_trans_Data1,self.err_jitter_trans_Data2,self.err_jitter_trans_Data3,self.err_jitter_trans_Data4,
                                      self.err_jitter_trans_Data5,self.err_jitter_trans_Data6,self.err_jitter_trans_Data7,self.err_jitter_trans_Data8,
                                      self.err_jitter_trans_Data9,self.err_jitter_trans_Data10]

        for i in range(10):
            tra_data_errors_gui[i].setText("+/- %.3f"%max(np.abs(fit.tra_off_err[i])))
            tra_data_errors_jitter_gui[i].setText("+/- %.3f"%max(np.abs(fit.tra_jitt_err[i])))            
                      

        self.err_RV_lin_trend.setText("+/- %.8f"%(max(fit.param_errors.linear_trend_error)))
        
        gp_rot_errors_gui = [self.err_rot_kernel_Amp,
                     self.err_rot_kernel_time_sc,
                     self.err_rot_kernel_Per,
                     self.err_rot_kernel_fact]
        
        for i in range(len(gp_rot_errors_gui)):
            gp_rot_errors_gui[i].setText("+/- %.3f"%max(np.abs(fit.param_errors.GP_params_errors[i])))  
            
        gp_sho_errors_gui = [self.err_sho_kernel_S,
                     self.err_sho_kernel_Q,
                     self.err_sho_kernel_omega]
        
        for i in range(len(gp_sho_errors_gui)):
            gp_sho_errors_gui[i].setText("+/- %.3f"%max(np.abs(fit.param_errors.GP_params_errors[i])))  

    def update_a_mass(self):
        global fit

        param_a_gui = [self.label_a1, self.label_a2, self.label_a3, self.label_a4, self.label_a5, 
                       self.label_a6, self.label_a7, self.label_a8, self.label_a9]
        param_mass_gui = [self.label_mass1, self.label_mass2, self.label_mass3, self.label_mass4, self.label_mass5, 
                       self.label_mass6, self.label_mass7, self.label_mass8, self.label_mass9]
        param_t_peri_gui = [self.label_t_peri1, self.label_t_peri2, self.label_t_peri3, self.label_t_peri4, self.label_t_peri5, 
                       self.label_t_peri6, self.label_t_peri7, self.label_t_peri8, self.label_t_peri9]

        if self.radioButton_RV.isChecked():
            for i in range(fit.npl):
                param_a_gui[i].setText("%.3f"%(fit.fit_results.a[i])) 
                param_mass_gui[i].setText("%.3f"%(fit.fit_results.mass[i])) 
                #param_t_peri_gui[i].setText("%.3f"%( (float(fit.epoch) - (np.radians(fit.params.planet_params[7*i + 4])/(2*np.pi))*fit.params.planet_params[7*i + 1] ))) # epoch  - ((ma/TWOPI)*a[1])
                param_t_peri_gui[i].setText("%.3f"%(fit.t_peri[i]))




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
            
            
        use_param_gui_wd = [self.use_om_dot_1, self.use_om_dot_2, self.use_om_dot_3, 
                            self.use_om_dot_4, self.use_om_dot_5, self.use_om_dot_6, 
                            self.use_om_dot_7, self.use_om_dot_8, self.use_om_dot_9]

        for i in range(fit.npl):
            use_param_gui_wd[i].setChecked(bool(fit.omega_dot_use[i]))           
 
            
        use_param_gui_tr = [self.use_t0_1, self.use_pl_rad_1, self.use_a_sol_1,
             self.use_t0_2, self.use_pl_rad_2, self.use_a_sol_2,
             self.use_t0_3, self.use_pl_rad_3, self.use_a_sol_3,
             self.use_t0_4, self.use_pl_rad_4, self.use_a_sol_4, 
             self.use_t0_5, self.use_pl_rad_5, self.use_a_sol_5,
             self.use_t0_6, self.use_pl_rad_6, self.use_a_sol_6,
             self.use_t0_7, self.use_pl_rad_7, self.use_a_sol_7, 
             self.use_t0_8, self.use_pl_rad_8, self.use_a_sol_8,
             self.use_t0_9, self.use_pl_rad_9, self.use_a_sol_9,
             ]
         
        for i in range(fit.npl):         
            use_param_gui_tr[i*3].setChecked(bool(fit.t0_use[i]) )        
            use_param_gui_tr[i*3+1].setChecked(bool(fit.pl_rad_use[i]) )
            use_param_gui_tr[i*3+2].setChecked(bool(fit.pl_a_use [i]) )
 

        use_data_offset_gui = [self.use_offset_Data1,self.use_offset_Data2,self.use_offset_Data3,self.use_offset_Data4,
                               self.use_offset_Data5,self.use_offset_Data6,self.use_offset_Data7,self.use_offset_Data8,
                               self.use_offset_Data9,self.use_offset_Data10]
        use_data_jitter_gui = [self.use_jitter_Data1,self.use_jitter_Data2,self.use_jitter_Data3,self.use_jitter_Data4,self.use_jitter_Data5,
                               self.use_jitter_Data6,self.use_jitter_Data7,self.use_jitter_Data8,self.use_jitter_Data9,self.use_jitter_Data10]

        for i in range(10): 
            #use_data_gui[i].setChecked(bool(fit.use.use_offsets[i])) # attention, TBF
            use_data_jitter_gui[i].setChecked(bool(fit.use.use_jitters[i]))
            use_data_offset_gui[i].setChecked(bool(fit.use.use_offsets[i])) 
            
            
        use_tra_data_offset_gui = [self.use_offset_trans_Data1,self.use_offset_trans_Data2,self.use_offset_trans_Data3,self.use_offset_trans_Data4,
                                   self.use_offset_trans_Data5,self.use_offset_trans_Data6,self.use_offset_trans_Data7,self.use_offset_trans_Data8,
                                   self.use_offset_trans_Data9,self.use_offset_trans_Data10]
        use_tra_data_jitter_gui = [self.use_jitter_trans_Data1,self.use_jitter_trans_Data2,self.use_jitter_trans_Data3,self.use_jitter_trans_Data4,self.use_jitter_trans_Data5,
                                   self.use_jitter_trans_Data6,self.use_jitter_trans_Data7,self.use_jitter_trans_Data8,self.use_jitter_trans_Data9,self.use_jitter_trans_Data10]

        for i in range(10): 
            #use_data_gui[i].setChecked(bool(fit.use.use_offsets[i])) # attention, TBF
            use_tra_data_jitter_gui[i].setChecked(bool(fit.tra_jitt_use[i]))
            use_tra_data_offset_gui[i].setChecked(bool(fit.tra_off_use[i]))             
            

        planet_checked_gui = [self.use_Planet1,self.use_Planet2,self.use_Planet3,self.use_Planet4,self.use_Planet5,
                              self.use_Planet6,self.use_Planet7,self.use_Planet8,self.use_Planet9]
        for i in range(9):  
            if i < fit.npl:
                planet_checked_gui[i].setChecked(True)  
            else:
                planet_checked_gui[i].setChecked(False)  
            
        self.use_RV_lin_trend.setChecked(bool(fit.use.use_linear_trend)) 
        

        use_gp_rot_params = [self.use_GP_rot_kernel_Amp,
                         self.use_GP_rot_kernel_time_sc,
                         self.use_GP_rot_kernel_Per,
                         self.use_GP_rot_kernel_fact]
                    
        
        for i in range(len(use_gp_rot_params)):
            use_gp_rot_params[i].setChecked(bool(fit.GP_rot_use[i]))
 
    
        use_gp_sho_params = [self.use_GP_sho_kernel_S,
                         self.use_GP_sho_kernel_Q,
                         self.use_GP_sho_kernel_omega]
                    
        
        for i in range(len(use_gp_sho_params)):
            use_gp_sho_params[i].setChecked(bool(fit.GP_sho_use[i]))


    def update_mixed_fitting(self):
        global fit
        
        fit.mixed_fit[0] =  int(self.use_mix_fitting.isChecked())
        fit.mixed_fit[1] = [int(self.mix_pl_1.isChecked()),int(self.mix_pl_2.isChecked()),int(self.mix_pl_3.isChecked()),
                            int(self.mix_pl_4.isChecked()),int(self.mix_pl_5.isChecked()),int(self.mix_pl_6.isChecked()),
                            int(self.mix_pl_7.isChecked()),int(self.mix_pl_8.isChecked()),int(self.mix_pl_9.isChecked()),
                            ]       
        
 
    def update_use(self):
        global fit
        
        use_planet_gui = [self.use_Planet1,self.use_Planet2,self.use_Planet3,self.use_Planet4,self.use_Planet5,
                          self.use_Planet6,self.use_Planet7,self.use_Planet8,self.use_Planet9]
        #for i in range(len(use_planet_gui)):  
        npl_old = fit.npl
        checked = int(np.sum( [use_planet_gui[i].isChecked() for i in range(len(use_planet_gui))] ))
 
        if npl_old < checked:
            fit.add_planet()
        elif npl_old >= checked:
            fit.npl = checked     
            
        #for i in range(len(use_planet_gui)):
        #    if use_planet_gui[i].isChecked() == False:
       #         fit.remove_planet(i)  
       #     else:
       #         fit.add_planet(i)  
       #         

        use_param_gui  = [self.use_K1, self.use_P1, self.use_e1, self.use_om1, self.use_ma1, self.use_incl1, self.use_Omega1,
                          self.use_K2, self.use_P2, self.use_e2, self.use_om2, self.use_ma2, self.use_incl2, self.use_Omega2,
                          self.use_K3, self.use_P3, self.use_e3, self.use_om3, self.use_ma3, self.use_incl3, self.use_Omega3,                        
                          self.use_K4, self.use_P4, self.use_e4, self.use_om4, self.use_ma4, self.use_incl4, self.use_Omega4,    
                          self.use_K5, self.use_P5, self.use_e5, self.use_om5, self.use_ma5, self.use_incl5, self.use_Omega5,    
                          self.use_K6, self.use_P6, self.use_e6, self.use_om6, self.use_ma6, self.use_incl6, self.use_Omega6, 
                          self.use_K7, self.use_P7, self.use_e7, self.use_om7, self.use_ma7, self.use_incl7, self.use_Omega7,    
                          self.use_K8, self.use_P8, self.use_e8, self.use_om8, self.use_ma8, self.use_incl8, self.use_Omega8,    
                          self.use_K9, self.use_P9, self.use_e9, self.use_om9, self.use_ma9, self.use_incl9, self.use_Omega9,                       
                          ]

        for i in range(9*7):
            fit.use.use_planet_params[i] = int(use_param_gui[i].isChecked())         


        use_param_gui_wd = [self.use_om_dot_1, self.use_om_dot_2, self.use_om_dot_3, 
                            self.use_om_dot_4, self.use_om_dot_5, self.use_om_dot_6, 
                            self.use_om_dot_7, self.use_om_dot_8, self.use_om_dot_9]

        for i in range(9):
            fit.omega_dot_use[i] = int(use_param_gui_wd[i].isChecked())
           # print(self.buttonGroup_use_planets.buttons()[i].checkedId())
       # print("###")
              
       # print([i for i, button in enumerate(self.buttonGroup_use_planets.buttons()) if button.isChecked()])
            
        use_param_gui_tr = [
             self.use_t0_1, self.use_pl_rad_1, self.use_a_sol_1,
             self.use_t0_2, self.use_pl_rad_2, self.use_a_sol_2,
             self.use_t0_3, self.use_pl_rad_3, self.use_a_sol_3,
             self.use_t0_4, self.use_pl_rad_4, self.use_a_sol_4, 
             self.use_t0_5, self.use_pl_rad_5, self.use_a_sol_5,
             self.use_t0_6, self.use_pl_rad_6, self.use_a_sol_6,
             self.use_t0_7, self.use_pl_rad_7, self.use_a_sol_7, 
             self.use_t0_8, self.use_pl_rad_8, self.use_a_sol_8,
             self.use_t0_9, self.use_pl_rad_9, self.use_a_sol_9,
             ]
         
        for i in range(9):        

            fit.t0_use[i] =  use_param_gui_tr[i*3].isChecked()  
            fit.pl_rad_use[i] = use_param_gui_tr[i*3+1].isChecked()  
            fit.pl_a_use[i] =  use_param_gui_tr[i*3+2].isChecked()  
       

        use_data_offset_gui = [self.use_offset_Data1,self.use_offset_Data2,self.use_offset_Data3,self.use_offset_Data4,
                               self.use_offset_Data5,self.use_offset_Data6,self.use_offset_Data7,self.use_offset_Data8,
                               self.use_offset_Data9,self.use_offset_Data10]
        use_data_jitter_gui = [self.use_jitter_Data1,self.use_jitter_Data2,self.use_jitter_Data3,self.use_jitter_Data4,self.use_jitter_Data5,
                               self.use_jitter_Data6,self.use_jitter_Data7,self.use_jitter_Data8,self.use_jitter_Data9,self.use_jitter_Data10]

        for i in range(10): 
            fit.use.use_jitters[i] = int(use_data_jitter_gui[i].isChecked())
            fit.use.use_offsets[i] = int(use_data_offset_gui[i].isChecked())   


        use_tra_data_offset_gui = [self.use_offset_trans_Data1,self.use_offset_trans_Data2,self.use_offset_trans_Data3,self.use_offset_trans_Data4,
                                   self.use_offset_trans_Data5,self.use_offset_trans_Data6,self.use_offset_trans_Data7,self.use_offset_trans_Data8,
                                   self.use_offset_trans_Data9,self.use_offset_trans_Data10]
        use_tra_data_jitter_gui = [self.use_jitter_trans_Data1,self.use_jitter_trans_Data2,self.use_jitter_trans_Data3,self.use_jitter_trans_Data4,self.use_jitter_trans_Data5,
                                   self.use_jitter_trans_Data6,self.use_jitter_trans_Data7,self.use_jitter_trans_Data8,self.use_jitter_trans_Data9,self.use_jitter_trans_Data10]

        for i in range(10): 
            fit.tra_jitt_use[i] = int(use_tra_data_jitter_gui[i].isChecked())
            fit.tra_off_use[i]  = int(use_tra_data_offset_gui[i].isChecked())

        fit.use.use_linear_trend = int(self.use_RV_lin_trend.isChecked()) 


        use_gp_rot_params = [self.use_GP_rot_kernel_Amp,
                         self.use_GP_rot_kernel_time_sc,
                         self.use_GP_rot_kernel_Per,
                         self.use_GP_rot_kernel_fact]
                    
        
        for i in range(len(use_gp_rot_params)):
            fit.GP_rot_use[i] = int(use_gp_rot_params[i].isChecked())
            
        use_gp_sho_params = [self.use_GP_sho_kernel_S,
                         self.use_GP_sho_kernel_Q,
                         self.use_GP_sho_kernel_omega]
                    
        
        for i in range(len(use_gp_sho_params)):
            fit.GP_sho_use[i] = int(use_gp_sho_params[i].isChecked())  
            
        self.update_tra_GP_use()
        
        

    def update_tra_GP_use(self):
        global fit
            
        use_tra_gp_rot_params = [self.use_tra_GP_rot_kernel_Amp,
                         self.use_tra_GP_rot_kernel_time_sc,
                         self.use_tra_GP_rot_kernel_Per,
                         self.use_tra_GP_rot_kernel_fact]
                    
        
        for i in range(len(use_tra_gp_rot_params)):
            fit.tra_GP_rot_use[i] = int(use_tra_gp_rot_params[i].isChecked())
            
        use_tra_gp_sho_params = [self.use_tra_GP_sho_kernel_S,
                         self.use_tra_GP_sho_kernel_Q,
                         self.use_tra_GP_sho_kernel_omega]
                    
        
        for i in range(len(use_tra_gp_sho_params)):
            fit.tra_GP_sho_use[i] = int(use_tra_gp_sho_params[i].isChecked())              
            

 
    def check_bounds(self):
        global fit
 
        
        param_bounds_gui = [
        [self.K_min_1.value(),self.K_max_1.value()],[self.P_min_1.value(),self.P_max_1.value()], [self.e_min_1.value(),self.e_max_1.value()],[self.om_min_1.value(),self.om_max_1.value()], [self.ma_min_1.value(),self.ma_max_1.value()],[self.incl_min_1.value(),self.incl_max_1.value()], [self.Omega_min_1.value(),self.Omega_max_1.value()],[self.t0_min_1.value(),self.t0_max_1.value()],[self.pl_rad_min_1.value(),self.pl_rad_max_1.value()],[self.a_sol_min_1.value(),self.a_sol_max_1.value()],
        [self.K_min_2.value(),self.K_max_2.value()],[self.P_min_2.value(),self.P_max_2.value()], [self.e_min_2.value(),self.e_max_2.value()],[self.om_min_2.value(),self.om_max_2.value()], [self.ma_min_2.value(),self.ma_max_2.value()],[self.incl_min_2.value(),self.incl_max_2.value()], [self.Omega_min_2.value(),self.Omega_max_2.value()],[self.t0_min_2.value(),self.t0_max_2.value()],[self.pl_rad_min_2.value(),self.pl_rad_max_2.value()],[self.a_sol_min_2.value(),self.a_sol_max_2.value()],
        [self.K_min_3.value(),self.K_max_3.value()],[self.P_min_3.value(),self.P_max_3.value()], [self.e_min_3.value(),self.e_max_3.value()],[self.om_min_3.value(),self.om_max_3.value()], [self.ma_min_3.value(),self.ma_max_3.value()],[self.incl_min_3.value(),self.incl_max_3.value()], [self.Omega_min_3.value(),self.Omega_max_3.value()],[self.t0_min_3.value(),self.t0_max_3.value()],[self.pl_rad_min_3.value(),self.pl_rad_max_3.value()],[self.a_sol_min_3.value(),self.a_sol_max_3.value()],
        [self.K_min_4.value(),self.K_max_4.value()],[self.P_min_4.value(),self.P_max_4.value()], [self.e_min_4.value(),self.e_max_4.value()],[self.om_min_4.value(),self.om_max_4.value()], [self.ma_min_4.value(),self.ma_max_4.value()],[self.incl_min_4.value(),self.incl_max_4.value()], [self.Omega_min_4.value(),self.Omega_max_4.value()],[self.t0_min_4.value(),self.t0_max_4.value()],[self.pl_rad_min_4.value(),self.pl_rad_max_4.value()],[self.a_sol_min_4.value(),self.a_sol_max_4.value()],
        [self.K_min_5.value(),self.K_max_5.value()],[self.P_min_5.value(),self.P_max_5.value()], [self.e_min_5.value(),self.e_max_5.value()],[self.om_min_5.value(),self.om_max_5.value()], [self.ma_min_5.value(),self.ma_max_5.value()],[self.incl_min_5.value(),self.incl_max_5.value()], [self.Omega_min_5.value(),self.Omega_max_5.value()],[self.t0_min_5.value(),self.t0_max_5.value()],[self.pl_rad_min_5.value(),self.pl_rad_max_5.value()],[self.a_sol_min_5.value(),self.a_sol_max_5.value()],
        [self.K_min_6.value(),self.K_max_6.value()],[self.P_min_6.value(),self.P_max_6.value()], [self.e_min_6.value(),self.e_max_6.value()],[self.om_min_6.value(),self.om_max_6.value()], [self.ma_min_6.value(),self.ma_max_6.value()],[self.incl_min_6.value(),self.incl_max_6.value()], [self.Omega_min_6.value(),self.Omega_max_6.value()],[self.t0_min_6.value(),self.t0_max_6.value()],[self.pl_rad_min_6.value(),self.pl_rad_max_6.value()],[self.a_sol_min_6.value(),self.a_sol_max_6.value()],
        [self.K_min_7.value(),self.K_max_7.value()],[self.P_min_7.value(),self.P_max_7.value()], [self.e_min_7.value(),self.e_max_7.value()],[self.om_min_7.value(),self.om_max_7.value()], [self.ma_min_7.value(),self.ma_max_7.value()],[self.incl_min_7.value(),self.incl_max_7.value()], [self.Omega_min_7.value(),self.Omega_max_7.value()],[self.t0_min_7.value(),self.t0_max_7.value()],[self.pl_rad_min_7.value(),self.pl_rad_max_7.value()],[self.a_sol_min_7.value(),self.a_sol_max_7.value()],
        [self.K_min_8.value(),self.K_max_8.value()],[self.P_min_8.value(),self.P_max_8.value()], [self.e_min_8.value(),self.e_max_8.value()],[self.om_min_8.value(),self.om_max_8.value()], [self.ma_min_8.value(),self.ma_max_8.value()],[self.incl_min_8.value(),self.incl_max_8.value()], [self.Omega_min_8.value(),self.Omega_max_8.value()],[self.t0_min_8.value(),self.t0_max_8.value()],[self.pl_rad_min_8.value(),self.pl_rad_max_8.value()],[self.a_sol_min_8.value(),self.a_sol_max_8.value()],
        [self.K_min_9.value(),self.K_max_9.value()],[self.P_min_9.value(),self.P_max_9.value()], [self.e_min_9.value(),self.e_max_9.value()],[self.om_min_9.value(),self.om_max_9.value()], [self.ma_min_9.value(),self.ma_max_9.value()],[self.incl_min_9.value(),self.incl_max_9.value()], [self.Omega_min_9.value(),self.Omega_max_9.value()],[self.t0_min_9.value(),self.t0_max_9.value()],[self.pl_rad_min_9.value(),self.pl_rad_max_9.value()],[self.a_sol_min_9.value(),self.a_sol_max_9.value()]               
        ]
 
        for i in range(fit.npl):
            fit.K_bound[i] = param_bounds_gui[10*i + 0]    
            fit.P_bound[i] = param_bounds_gui[10*i + 1]    
            fit.e_bound[i] = param_bounds_gui[10*i + 2]    
            fit.w_bound[i] = param_bounds_gui[10*i + 3]    
            fit.M0_bound[i] = param_bounds_gui[10*i + 4]    
            fit.i_bound[i] = param_bounds_gui[10*i + 5]    
            fit.Node_bound[i] = param_bounds_gui[10*i + 6]    
            fit.t0_bound[i]  =  param_bounds_gui[10*i + 7]
            fit.pl_rad_bound[i]  =   param_bounds_gui[10*i + 8]
            fit.pl_a_bound[i]   =   param_bounds_gui[10*i + 9]

        offset_bounds_gui = [
        [self.Data1_min.value(),self.Data1_max.value()], [self.Data2_min.value(),self.Data2_max.value()], [self.Data3_min.value(),self.Data3_max.value()], [self.Data4_min.value(),self.Data4_max.value()], [self.Data5_min.value(),self.Data5_max.value()],   
        [self.Data6_min.value(),self.Data6_max.value()], [self.Data7_min.value(),self.Data7_max.value()], [self.Data8_min.value(),self.Data8_max.value()], [self.Data9_min.value(),self.Data9_max.value()], [self.Data10_min.value(),self.Data10_max.value()]
        ]
        
        jitter_bounds_gui = [
        [self.jitter1_min.value(),self.jitter1_max.value()], [self.jitter2_min.value(),self.jitter2_max.value()], [self.jitter3_min.value(),self.jitter3_max.value()], [self.jitter4_min.value(),self.jitter4_max.value()], [self.jitter5_min.value(),self.jitter5_max.value()],   
        [self.jitter6_min.value(),self.jitter6_max.value()], [self.jitter7_min.value(),self.jitter7_max.value()], [self.jitter8_min.value(),self.jitter8_max.value()], [self.jitter9_min.value(),self.jitter9_max.value()], [self.jitter10_min.value(),self.Data10_max.value()]   
        ]  
    
    
        for i in range(10): 
            fit.rvoff_bounds[i] = offset_bounds_gui[i]
            fit.jitt_bounds[i]  = jitter_bounds_gui[i] 
            
        om_dot_bounds_gui = [
        [self.omega_dot_min_1.value(),self.omega_dot_max_1.value()], [self.omega_dot_min_2.value(),self.omega_dot_max_2.value()], 
        [self.omega_dot_min_3.value(),self.omega_dot_max_3.value()], [self.omega_dot_min_4.value(),self.omega_dot_max_4.value()], 
        [self.omega_dot_min_5.value(),self.omega_dot_max_5.value()], [self.omega_dot_min_6.value(),self.omega_dot_max_6.value()], 
        [self.omega_dot_min_7.value(),self.omega_dot_max_7.value()], [self.omega_dot_min_8.value(),self.omega_dot_max_8.value()], 
        [self.omega_dot_min_9.value(),self.omega_dot_max_9.value()] 
        ]  
    
    
        for i in range(9): 
            fit.omega_dot_bounds[i] = om_dot_bounds_gui[i]
             
 
        fit.rv_lintr_bounds[0]  = [self.lin_trend_min.value(),self.lin_trend_max.value()]

        offset_bounds_gui_tra = [
        [self.tra_Data_min_1.value(),self.tra_Data_max_1.value()], [self.tra_Data_min_2.value(),self.tra_Data_max_2.value()], [self.tra_Data_min_3.value(),self.tra_Data_max_3.value()], [self.tra_Data_min_4.value(),self.tra_Data_max_4.value()], [self.tra_Data_min_5.value(),self.tra_Data_max_5.value()],   
        [self.tra_Data_min_6.value(),self.tra_Data_max_6.value()], [self.tra_Data_min_7.value(),self.tra_Data_max_7.value()], [self.tra_Data_min_8.value(),self.tra_Data_max_8.value()], [self.tra_Data_min_9.value(),self.tra_Data_max_9.value()], [self.tra_Data_min_10.value(),self.tra_Data_max_10.value()]
        ]
        
        jitter_bounds_gui_tra = [
        [self.tra_jitter_min_1.value(),self.tra_jitter_max_1.value()], [self.tra_jitter_min_2.value(),self.tra_jitter_max_2.value()], [self.tra_jitter_min_3.value(),self.tra_jitter_max_3.value()], [self.tra_jitter_min_4.value(),self.tra_jitter_max_4.value()], [self.tra_jitter_min_5.value(),self.tra_jitter_max_5.value()],   
        [self.tra_jitter_min_6.value(),self.tra_jitter_max_6.value()], [self.tra_jitter_min_7.value(),self.tra_jitter_max_7.value()], [self.tra_jitter_min_8.value(),self.tra_jitter_max_8.value()], [self.tra_jitter_min_9.value(),self.tra_jitter_max_9.value()], [self.tra_jitter_min_10.value(),self.tra_jitter_max_10.value()]
        ]
        
    
        for i in range(10): 
            fit.tra_off_bounds[i] = offset_bounds_gui_tra[i]
            fit.tra_jitt_bounds[i]  = jitter_bounds_gui_tra[i] 
            

        self.check_RV_GP_bounds()            
        self.check_tra_GP_bounds()
        
        
        
        
    def check_RV_GP_bounds(self):
        global fit
        
        GP_rot_bounds_gui = [
        [self.GP_rot_kernel_Amp_min.value(),self.GP_rot_kernel_Amp_max.value()],  
        [self.GP_rot_kernel_time_sc_min.value(),self.GP_rot_kernel_time_sc_max.value()],  
        [self.GP_rot_kernel_Per_min.value(),self.GP_rot_kernel_Per_max.value()],  
        [self.GP_rot_kernel_fact_min.value(),self.GP_rot_kernel_fact_max.value()],  
        ]
 
        for i in range(4): 
            fit.GP_rot_bounds[i] = GP_rot_bounds_gui[i]
            
        GP_sho_bounds_gui = [
        [self.GP_sho_kernel_S_min.value(),self.GP_sho_kernel_S_max.value()],  
        [self.GP_sho_kernel_Q_min.value(),self.GP_sho_kernel_Q_max.value()],  
        [self.GP_sho_kernel_omega_min.value(),self.GP_sho_kernel_omega_max.value()],  
        ]
 
        for i in range(3): 
            fit.GP_sho_bounds[i] = GP_sho_bounds_gui[i]        
                               
    def check_tra_GP_bounds(self):
        global fit

        tra_GP_rot_bounds_gui = [
        [self.tra_GP_rot_kernel_Amp_min.value(),self.tra_GP_rot_kernel_Amp_max.value()],  
        [self.tra_GP_rot_kernel_time_sc_min.value(),self.tra_GP_rot_kernel_time_sc_max.value()],  
        [self.tra_GP_rot_kernel_Per_min.value(),self.tra_GP_rot_kernel_Per_max.value()],  
        [self.tra_GP_rot_kernel_fact_min.value(),self.tra_GP_rot_kernel_fact_max.value()],  
        ]
 
        for i in range(4): 
            fit.tra_GP_rot_bounds[i] = tra_GP_rot_bounds_gui[i]
            
        tra_GP_sho_bounds_gui = [
        [self.tra_GP_sho_kernel_S_min.value(),self.tra_GP_sho_kernel_S_max.value()],  
        [self.tra_GP_sho_kernel_Q_min.value(),self.tra_GP_sho_kernel_Q_max.value()],  
        [self.tra_GP_sho_kernel_omega_min.value(),self.tra_GP_sho_kernel_omega_max.value()],  
        ]
 
        for i in range(3): 
            fit.tra_GP_sho_bounds[i] = tra_GP_sho_bounds_gui[i]                  
            
            
    def check_priors_nr(self):
        global fit

        param_nr_priors_gui = [
        [self.K_mean_1.value(),self.K_sigma_1.value(),self.use_K_norm_pr_1.isChecked()],[self.P_mean_1.value(),self.P_sigma_1.value(),self.use_P_norm_pr_1.isChecked()], [self.e_mean_1.value(),self.e_sigma_1.value(),self.use_e_norm_pr_1.isChecked()],[self.om_mean_1.value(),self.om_sigma_1.value(),self.use_om_norm_pr_1.isChecked()], [self.ma_mean_1.value(),self.ma_sigma_1.value(),self.use_ma_norm_pr_1.isChecked()],[self.incl_mean_1.value(),self.incl_sigma_1.value(),self.use_incl_norm_pr_1.isChecked()], [self.Omega_mean_1.value(),self.Omega_sigma_1.value(), self.use_Omega_norm_pr_1.isChecked()],[self.t0_mean_1.value(),self.t0_sigma_1.value(), self.use_t0_norm_pr_1.isChecked()],[self.pl_rad_mean_1.value(),self.pl_rad_sigma_1.value(),self.use_pl_rad_norm_pr_1.isChecked()],[self.a_sol_mean_1.value(),self.a_sol_sigma_1.value(),self.use_a_sol_norm_pr_1.isChecked()],
        [self.K_mean_2.value(),self.K_sigma_2.value(),self.use_K_norm_pr_2.isChecked()],[self.P_mean_2.value(),self.P_sigma_2.value(),self.use_P_norm_pr_2.isChecked()], [self.e_mean_2.value(),self.e_sigma_2.value(),self.use_e_norm_pr_2.isChecked()],[self.om_mean_2.value(),self.om_sigma_2.value(),self.use_om_norm_pr_2.isChecked()], [self.ma_mean_2.value(),self.ma_sigma_2.value(),self.use_ma_norm_pr_2.isChecked()],[self.incl_mean_2.value(),self.incl_sigma_2.value(),self.use_incl_norm_pr_2.isChecked()], [self.Omega_mean_2.value(),self.Omega_sigma_2.value(), self.use_Omega_norm_pr_2.isChecked()],[self.t0_mean_2.value(),self.t0_sigma_2.value(), self.use_t0_norm_pr_2.isChecked()],[self.pl_rad_mean_2.value(),self.pl_rad_sigma_2.value(),self.use_pl_rad_norm_pr_2.isChecked()],[self.a_sol_mean_2.value(),self.a_sol_sigma_2.value(),self.use_a_sol_norm_pr_2.isChecked()],
        [self.K_mean_3.value(),self.K_sigma_3.value(),self.use_K_norm_pr_3.isChecked()],[self.P_mean_3.value(),self.P_sigma_3.value(),self.use_P_norm_pr_3.isChecked()], [self.e_mean_3.value(),self.e_sigma_3.value(),self.use_e_norm_pr_3.isChecked()],[self.om_mean_3.value(),self.om_sigma_3.value(),self.use_om_norm_pr_3.isChecked()], [self.ma_mean_3.value(),self.ma_sigma_3.value(),self.use_ma_norm_pr_3.isChecked()],[self.incl_mean_3.value(),self.incl_sigma_3.value(),self.use_incl_norm_pr_3.isChecked()], [self.Omega_mean_3.value(),self.Omega_sigma_3.value(), self.use_Omega_norm_pr_3.isChecked()],[self.t0_mean_3.value(),self.t0_sigma_3.value(), self.use_t0_norm_pr_3.isChecked()],[self.pl_rad_mean_3.value(),self.pl_rad_sigma_3.value(),self.use_pl_rad_norm_pr_3.isChecked()],[self.a_sol_mean_3.value(),self.a_sol_sigma_3.value(),self.use_a_sol_norm_pr_3.isChecked()],
        [self.K_mean_4.value(),self.K_sigma_4.value(),self.use_K_norm_pr_4.isChecked()],[self.P_mean_4.value(),self.P_sigma_4.value(),self.use_P_norm_pr_4.isChecked()], [self.e_mean_4.value(),self.e_sigma_4.value(),self.use_e_norm_pr_4.isChecked()],[self.om_mean_4.value(),self.om_sigma_4.value(),self.use_om_norm_pr_4.isChecked()], [self.ma_mean_4.value(),self.ma_sigma_4.value(),self.use_ma_norm_pr_4.isChecked()],[self.incl_mean_4.value(),self.incl_sigma_4.value(),self.use_incl_norm_pr_4.isChecked()], [self.Omega_mean_4.value(),self.Omega_sigma_4.value(), self.use_Omega_norm_pr_4.isChecked()],[self.t0_mean_4.value(),self.t0_sigma_4.value(), self.use_t0_norm_pr_4.isChecked()],[self.pl_rad_mean_4.value(),self.pl_rad_sigma_4.value(),self.use_pl_rad_norm_pr_4.isChecked()],[self.a_sol_mean_4.value(),self.a_sol_sigma_4.value(),self.use_a_sol_norm_pr_4.isChecked()],
        [self.K_mean_5.value(),self.K_sigma_5.value(),self.use_K_norm_pr_5.isChecked()],[self.P_mean_5.value(),self.P_sigma_5.value(),self.use_P_norm_pr_5.isChecked()], [self.e_mean_5.value(),self.e_sigma_5.value(),self.use_e_norm_pr_5.isChecked()],[self.om_mean_5.value(),self.om_sigma_5.value(),self.use_om_norm_pr_5.isChecked()], [self.ma_mean_5.value(),self.ma_sigma_5.value(),self.use_ma_norm_pr_5.isChecked()],[self.incl_mean_5.value(),self.incl_sigma_5.value(),self.use_incl_norm_pr_5.isChecked()], [self.Omega_mean_5.value(),self.Omega_sigma_5.value(), self.use_Omega_norm_pr_5.isChecked()],[self.t0_mean_5.value(),self.t0_sigma_5.value(), self.use_t0_norm_pr_5.isChecked()],[self.pl_rad_mean_5.value(),self.pl_rad_sigma_5.value(),self.use_pl_rad_norm_pr_5.isChecked()],[self.a_sol_mean_5.value(),self.a_sol_sigma_5.value(),self.use_a_sol_norm_pr_5.isChecked()],
        [self.K_mean_6.value(),self.K_sigma_6.value(),self.use_K_norm_pr_6.isChecked()],[self.P_mean_6.value(),self.P_sigma_6.value(),self.use_P_norm_pr_6.isChecked()], [self.e_mean_6.value(),self.e_sigma_6.value(),self.use_e_norm_pr_6.isChecked()],[self.om_mean_6.value(),self.om_sigma_6.value(),self.use_om_norm_pr_6.isChecked()], [self.ma_mean_6.value(),self.ma_sigma_6.value(),self.use_ma_norm_pr_6.isChecked()],[self.incl_mean_6.value(),self.incl_sigma_6.value(),self.use_incl_norm_pr_6.isChecked()], [self.Omega_mean_6.value(),self.Omega_sigma_6.value(), self.use_Omega_norm_pr_6.isChecked()],[self.t0_mean_6.value(),self.t0_sigma_6.value(), self.use_t0_norm_pr_6.isChecked()],[self.pl_rad_mean_6.value(),self.pl_rad_sigma_6.value(),self.use_pl_rad_norm_pr_6.isChecked()],[self.a_sol_mean_6.value(),self.a_sol_sigma_6.value(),self.use_a_sol_norm_pr_6.isChecked()],
        [self.K_mean_7.value(),self.K_sigma_7.value(),self.use_K_norm_pr_7.isChecked()],[self.P_mean_7.value(),self.P_sigma_7.value(),self.use_P_norm_pr_7.isChecked()], [self.e_mean_7.value(),self.e_sigma_7.value(),self.use_e_norm_pr_7.isChecked()],[self.om_mean_7.value(),self.om_sigma_7.value(),self.use_om_norm_pr_7.isChecked()], [self.ma_mean_7.value(),self.ma_sigma_7.value(),self.use_ma_norm_pr_7.isChecked()],[self.incl_mean_7.value(),self.incl_sigma_7.value(),self.use_incl_norm_pr_7.isChecked()], [self.Omega_mean_7.value(),self.Omega_sigma_7.value(), self.use_Omega_norm_pr_7.isChecked()],[self.t0_mean_7.value(),self.t0_sigma_7.value(), self.use_t0_norm_pr_7.isChecked()],[self.pl_rad_mean_7.value(),self.pl_rad_sigma_7.value(),self.use_pl_rad_norm_pr_7.isChecked()],[self.a_sol_mean_7.value(),self.a_sol_sigma_7.value(),self.use_a_sol_norm_pr_7.isChecked()],
        [self.K_mean_8.value(),self.K_sigma_8.value(),self.use_K_norm_pr_8.isChecked()],[self.P_mean_8.value(),self.P_sigma_8.value(),self.use_P_norm_pr_8.isChecked()], [self.e_mean_8.value(),self.e_sigma_8.value(),self.use_e_norm_pr_8.isChecked()],[self.om_mean_8.value(),self.om_sigma_8.value(),self.use_om_norm_pr_8.isChecked()], [self.ma_mean_8.value(),self.ma_sigma_8.value(),self.use_ma_norm_pr_8.isChecked()],[self.incl_mean_8.value(),self.incl_sigma_8.value(),self.use_incl_norm_pr_8.isChecked()], [self.Omega_mean_8.value(),self.Omega_sigma_8.value(), self.use_Omega_norm_pr_8.isChecked()],[self.t0_mean_8.value(),self.t0_sigma_8.value(), self.use_t0_norm_pr_8.isChecked()],[self.pl_rad_mean_8.value(),self.pl_rad_sigma_8.value(),self.use_pl_rad_norm_pr_8.isChecked()],[self.a_sol_mean_8.value(),self.a_sol_sigma_8.value(),self.use_a_sol_norm_pr_8.isChecked()],
        [self.K_mean_9.value(),self.K_sigma_9.value(),self.use_K_norm_pr_9.isChecked()],[self.P_mean_9.value(),self.P_sigma_9.value(),self.use_P_norm_pr_9.isChecked()], [self.e_mean_9.value(),self.e_sigma_9.value(),self.use_e_norm_pr_9.isChecked()],[self.om_mean_9.value(),self.om_sigma_9.value(),self.use_om_norm_pr_9.isChecked()], [self.ma_mean_9.value(),self.ma_sigma_9.value(),self.use_ma_norm_pr_9.isChecked()],[self.incl_mean_9.value(),self.incl_sigma_9.value(),self.use_incl_norm_pr_9.isChecked()], [self.Omega_mean_9.value(),self.Omega_sigma_9.value(), self.use_Omega_norm_pr_9.isChecked()],[self.t0_mean_9.value(),self.t0_sigma_9.value(), self.use_t0_norm_pr_9.isChecked()],[self.pl_rad_mean_9.value(),self.pl_rad_sigma_9.value(),self.use_pl_rad_norm_pr_9.isChecked()],[self.a_sol_mean_9.value(),self.a_sol_sigma_9.value(),self.use_a_sol_norm_pr_9.isChecked()],
        ]
 
        for i in range(fit.npl):
            fit.K_norm_pr[i]  = param_nr_priors_gui[10*i + 0]    
            fit.P_norm_pr[i]  = param_nr_priors_gui[10*i + 1]    
            fit.e_norm_pr[i]  = param_nr_priors_gui[10*i + 2]    
            fit.w_norm_pr[i]  = param_nr_priors_gui[10*i + 3]    
            fit.M0_norm_pr[i] = param_nr_priors_gui[10*i + 4]    
            fit.i_norm_pr[i]  = param_nr_priors_gui[10*i + 5]    
            fit.Node_norm_pr[i] = param_nr_priors_gui[10*i + 6]    
            fit.t0_norm_pr[i]   = param_nr_priors_gui[10*i + 7]
            fit.pl_rad_norm_pr[i]  = param_nr_priors_gui[10*i + 8]
            fit.pl_a_norm_pr[i]    = param_nr_priors_gui[10*i + 9]

        offset_nr_priors_gui = [
        [self.RV_Data_mean_1.value(),self.RV_Data_sigma_1.value(),self.use_rvoff_nr_1.isChecked()], 
        [self.RV_Data_mean_2.value(),self.RV_Data_sigma_2.value(),self.use_rvoff_nr_2.isChecked()], 
        [self.RV_Data_mean_3.value(),self.RV_Data_sigma_3.value(),self.use_rvoff_nr_3.isChecked()], 
        [self.RV_Data_mean_4.value(),self.RV_Data_sigma_4.value(),self.use_rvoff_nr_4.isChecked()], 
        [self.RV_Data_mean_5.value(),self.RV_Data_sigma_5.value(),self.use_rvoff_nr_5.isChecked()], 
        [self.RV_Data_mean_6.value(),self.RV_Data_sigma_6.value(),self.use_rvoff_nr_6.isChecked()], 
        [self.RV_Data_mean_7.value(),self.RV_Data_sigma_7.value(),self.use_rvoff_nr_7.isChecked()], 
        [self.RV_Data_mean_8.value(),self.RV_Data_sigma_8.value(),self.use_rvoff_nr_8.isChecked()], 
        [self.RV_Data_mean_9.value(),self.RV_Data_sigma_9.value(),self.use_rvoff_nr_9.isChecked()], 
        [self.RV_Data_mean_10.value(),self.RV_Data_sigma_10.value(),self.use_rvoff_nr_10.isChecked()]
        ]
        
        jitter_nr_priors_gui = [
        [self.RV_jitter_mean_1.value(),self.RV_jitter_sigma_1.value(),self.use_rvjitt_nr_1.isChecked()], 
        [self.RV_jitter_mean_2.value(),self.RV_jitter_sigma_2.value(),self.use_rvjitt_nr_2.isChecked()], 
        [self.RV_jitter_mean_3.value(),self.RV_jitter_sigma_3.value(),self.use_rvjitt_nr_3.isChecked()], 
        [self.RV_jitter_mean_4.value(),self.RV_jitter_sigma_4.value(),self.use_rvjitt_nr_4.isChecked()], 
        [self.RV_jitter_mean_5.value(),self.RV_jitter_sigma_5.value(),self.use_rvjitt_nr_5.isChecked()], 
        [self.RV_jitter_mean_6.value(),self.RV_jitter_sigma_6.value(),self.use_rvjitt_nr_6.isChecked()],
        [self.RV_jitter_mean_7.value(),self.RV_jitter_sigma_7.value(),self.use_rvjitt_nr_7.isChecked()], 
        [self.RV_jitter_mean_8.value(),self.RV_jitter_sigma_8.value(),self.use_rvjitt_nr_8.isChecked()], 
        [self.RV_jitter_mean_9.value(),self.RV_jitter_sigma_9.value(),self.use_rvjitt_nr_9.isChecked()], 
        [self.RV_jitter_mean_10.value(),self.RV_jitter_sigma_10.value(),self.use_rvjitt_nr_10.isChecked()]   
        ]  
    
    
        for i in range(10): 
            fit.rvoff_norm_pr[i] = offset_nr_priors_gui[i]
            fit.jitt_norm_pr[i]  = jitter_nr_priors_gui[i] 
            
            
        om_nr_priors_gui = [
        [self.omega_dot_mean_1.value(),self.omega_dot_sigma_1.value(),self.use_omega_dot_norm_pr_1.isChecked()], [self.omega_dot_mean_2.value(),self.omega_dot_sigma_2.value(),self.use_omega_dot_norm_pr_2.isChecked()], 
        [self.omega_dot_mean_3.value(),self.omega_dot_sigma_3.value(),self.use_omega_dot_norm_pr_3.isChecked()], [self.omega_dot_mean_4.value(),self.omega_dot_sigma_4.value(),self.use_omega_dot_norm_pr_4.isChecked()], 
        [self.omega_dot_mean_5.value(),self.omega_dot_sigma_5.value(),self.use_omega_dot_norm_pr_5.isChecked()], [self.omega_dot_mean_6.value(),self.omega_dot_sigma_6.value(),self.use_omega_dot_norm_pr_6.isChecked()], 
        [self.omega_dot_mean_7.value(),self.omega_dot_sigma_7.value(),self.use_omega_dot_norm_pr_7.isChecked()], [self.omega_dot_mean_8.value(),self.omega_dot_sigma_8.value(),self.use_omega_dot_norm_pr_8.isChecked()], 
        [self.omega_dot_mean_9.value(),self.omega_dot_sigma_9.value(),self.use_omega_dot_norm_pr_9.isChecked()]   
        ]  
    
    
        for i in range(9): 
            fit.omega_dot_norm_pr[i] = om_nr_priors_gui[i]            
    
 
        fit.rv_lintr_norm_pr[0]  = [self.lin_trend_mean.value(),self.lin_trend_sigma.value(),self.use_lin_tr_nr_pr.isChecked()]

        offset_nr_priors_gui_tra = [
        [self.tra_Data_mean_1.value(),self.tra_Data_sigma_1.value(),self.use_traoff_nr_1.isChecked()], 
        [self.tra_Data_mean_2.value(),self.tra_Data_sigma_2.value(),self.use_traoff_nr_2.isChecked()], 
        [self.tra_Data_mean_3.value(),self.tra_Data_sigma_3.value(),self.use_traoff_nr_3.isChecked()],
        [self.tra_Data_mean_4.value(),self.tra_Data_sigma_4.value(),self.use_traoff_nr_4.isChecked()], 
        [self.tra_Data_mean_5.value(),self.tra_Data_sigma_5.value(),self.use_traoff_nr_5.isChecked()],   
        [self.tra_Data_mean_6.value(),self.tra_Data_sigma_6.value(),self.use_traoff_nr_6.isChecked()],
        [self.tra_Data_mean_7.value(),self.tra_Data_sigma_7.value(),self.use_traoff_nr_7.isChecked()],
        [self.tra_Data_mean_8.value(),self.tra_Data_sigma_8.value(),self.use_traoff_nr_8.isChecked()], 
        [self.tra_Data_mean_9.value(),self.tra_Data_sigma_9.value(),self.use_traoff_nr_9.isChecked()], 
        [self.tra_Data_mean_10.value(),self.tra_Data_sigma_10.value(),self.use_traoff_nr_10.isChecked()]
        ]
        
        jitter_nr_priors_gui_tra = [
        [self.tra_jitter_mean_1.value(),self.tra_jitter_sigma_1.value(),self.use_trajitt_nr_1.isChecked()], 
        [self.tra_jitter_mean_2.value(),self.tra_jitter_sigma_2.value(),self.use_trajitt_nr_2.isChecked()], 
        [self.tra_jitter_mean_3.value(),self.tra_jitter_sigma_3.value(),self.use_trajitt_nr_3.isChecked()], 
        [self.tra_jitter_mean_4.value(),self.tra_jitter_sigma_4.value(),self.use_trajitt_nr_4.isChecked()], 
        [self.tra_jitter_mean_5.value(),self.tra_jitter_sigma_5.value(),self.use_trajitt_nr_5.isChecked()],   
        [self.tra_jitter_mean_6.value(),self.tra_jitter_sigma_6.value(),self.use_trajitt_nr_6.isChecked()], 
        [self.tra_jitter_mean_7.value(),self.tra_jitter_sigma_7.value(),self.use_trajitt_nr_7.isChecked()], 
        [self.tra_jitter_mean_8.value(),self.tra_jitter_sigma_8.value(),self.use_trajitt_nr_8.isChecked()], 
        [self.tra_jitter_mean_9.value(),self.tra_jitter_sigma_9.value(),self.use_trajitt_nr_9.isChecked()], 
        [self.tra_jitter_mean_10.value(),self.tra_jitter_sigma_10.value(),self.use_trajitt_nr_10.isChecked()]
        ]
        
    
        for i in range(10): 
            fit.tra_off_norm_pr[i] = offset_nr_priors_gui_tra[i]
            fit.tra_jitt_norm_pr[i] = jitter_nr_priors_gui_tra[i] 

        self.check_RV_GP_priors_nr()
        self.check_tra_GP_priors_nr()
    
    

    def check_RV_GP_priors_nr(self):
        global fit

 
        GP_rot_nr_priors_gui = [
        [self.GP_rot_kernel_Amp_mean.value(),self.GP_rot_kernel_Amp_sigma.value(),self.use_GP_rot_kernel_Amp_nr_pr.isChecked()],  
        [self.GP_rot_kernel_time_sc_mean.value(),self.GP_rot_kernel_time_sc_sigma.value(),self.use_GP_rot_kernel_time_sc_nr_pr.isChecked()],  
        [self.GP_rot_kernel_Per_mean.value(),self.GP_rot_kernel_Per_sigma.value(),self.use_GP_rot_kernel_Per_sigma_nr_pr.isChecked()],  
        [self.GP_rot_kernel_fact_mean.value(),self.GP_rot_kernel_fact_sigma.value(),self.use_GP_rot_kernel_fact_nr_pr.isChecked()],  
        ]
 
        for i in range(4): 
            fit.GP_rot_norm_pr[i] = GP_rot_nr_priors_gui[i]            
    

        GP_sho_nr_priors_gui = [
        [self.GP_sho_kernel_S_mean.value(),self.GP_sho_kernel_S_sigma.value(), self.use_GP_sho_kernel_S_nr_pr.isChecked()],  
        [self.GP_sho_kernel_Q_mean.value(),self.GP_sho_kernel_Q_sigma.value(), self.use_GP_sho_kernel_Q_nr_pr.isChecked()],  
        [self.GP_sho_kernel_omega_mean.value(),self.GP_sho_kernel_omega_sigma.value(), self.use_GP_sho_kernel_omega_nr_pr.isChecked()],  
        ]
 
        for i in range(3): 
            fit.GP_sho_norm_pr[i] = GP_sho_nr_priors_gui[i]   
            

    def check_tra_GP_priors_nr(self):
        global fit

 
        tra_GP_rot_nr_priors_gui = [
        [self.tra_GP_rot_kernel_Amp_mean.value(),self.tra_GP_rot_kernel_Amp_sigma.value(),self.use_tra_GP_rot_kernel_Amp_nr_pr.isChecked()],  
        [self.tra_GP_rot_kernel_time_sc_mean.value(),self.tra_GP_rot_kernel_time_sc_sigma.value(),self.use_tra_GP_rot_kernel_time_sc_nr_pr.isChecked()],  
        [self.tra_GP_rot_kernel_Per_mean.value(),self.tra_GP_rot_kernel_Per_sigma.value(),self.use_tra_GP_rot_kernel_Per_sigma_nr_pr.isChecked()],  
        [self.tra_GP_rot_kernel_fact_mean.value(),self.tra_GP_rot_kernel_fact_sigma.value(),self.use_tra_GP_rot_kernel_fact_nr_pr.isChecked()],  
        ]
 
        for i in range(4): 
            fit.tra_GP_rot_norm_pr[i] = tra_GP_rot_nr_priors_gui[i]            
    

        tra_GP_sho_nr_priors_gui = [
        [self.tra_GP_sho_kernel_S_mean.value(),self.tra_GP_sho_kernel_S_sigma.value(), self.use_tra_GP_sho_kernel_S_nr_pr.isChecked()],  
        [self.tra_GP_sho_kernel_Q_mean.value(),self.tra_GP_sho_kernel_Q_sigma.value(), self.use_tra_GP_sho_kernel_Q_nr_pr.isChecked()],  
        [self.tra_GP_sho_kernel_omega_mean.value(),self.tra_GP_sho_kernel_omega_sigma.value(), self.use_tra_GP_sho_kernel_omega_nr_pr.isChecked()],  
        ]
 
        for i in range(3): 
            fit.tra_GP_sho_norm_pr[i] = tra_GP_sho_nr_priors_gui[i]               
            

            
    def check_priors_jeff(self):
        global fit

        param_jeff_priors_gui = [
        [self.K_jeff_alpha_1.value(),self.K_jeff_beta_1.value(),self.use_K_jeff_pr_1.isChecked()],[self.P_jeff_alpha_1.value(),self.P_jeff_beta_1.value(),self.use_P_jeff_pr_1.isChecked()], [self.e_jeff_alpha_1.value(),self.e_jeff_beta_1.value(),self.use_e_jeff_pr_1.isChecked()],[self.om_jeff_alpha_1.value(),self.om_jeff_beta_1.value(),self.use_om_jeff_pr_1.isChecked()], [self.ma_jeff_alpha_1.value(),self.ma_jeff_beta_1.value(),self.use_ma_jeff_pr_1.isChecked()],[self.incl_jeff_alpha_1.value(),self.incl_jeff_beta_1.value(),self.use_incl_jeff_pr_1.isChecked()], [self.Omega_jeff_alpha_1.value(),self.Omega_jeff_beta_1.value(), self.use_Omega_jeff_pr_1.isChecked()],[self.t0_jeff_alpha_1.value(),self.t0_jeff_beta_1.value(), self.use_t0_jeff_pr_1.isChecked()],[self.pl_rad_jeff_alpha_1.value(),self.pl_rad_jeff_beta_1.value(),self.use_pl_rad_jeff_pr_1.isChecked()],[self.a_sol_jeff_alpha_1.value(),self.a_sol_jeff_beta_1.value(),self.use_a_sol_jeff_pr_1.isChecked()],
        [self.K_jeff_alpha_2.value(),self.K_jeff_beta_2.value(),self.use_K_jeff_pr_2.isChecked()],[self.P_jeff_alpha_2.value(),self.P_jeff_beta_2.value(),self.use_P_jeff_pr_2.isChecked()], [self.e_jeff_alpha_2.value(),self.e_jeff_beta_2.value(),self.use_e_jeff_pr_2.isChecked()],[self.om_jeff_alpha_2.value(),self.om_jeff_beta_2.value(),self.use_om_jeff_pr_2.isChecked()], [self.ma_jeff_alpha_2.value(),self.ma_jeff_beta_2.value(),self.use_ma_jeff_pr_2.isChecked()],[self.incl_jeff_alpha_2.value(),self.incl_jeff_beta_2.value(),self.use_incl_jeff_pr_2.isChecked()], [self.Omega_jeff_alpha_2.value(),self.Omega_jeff_beta_2.value(), self.use_Omega_jeff_pr_2.isChecked()],[self.t0_jeff_alpha_2.value(),self.t0_jeff_beta_2.value(), self.use_t0_jeff_pr_2.isChecked()],[self.pl_rad_jeff_alpha_2.value(),self.pl_rad_jeff_beta_2.value(),self.use_pl_rad_jeff_pr_2.isChecked()],[self.a_sol_jeff_alpha_2.value(),self.a_sol_jeff_beta_2.value(),self.use_a_sol_jeff_pr_2.isChecked()],
        [self.K_jeff_alpha_3.value(),self.K_jeff_beta_3.value(),self.use_K_jeff_pr_3.isChecked()],[self.P_jeff_alpha_3.value(),self.P_jeff_beta_3.value(),self.use_P_jeff_pr_3.isChecked()], [self.e_jeff_alpha_3.value(),self.e_jeff_beta_3.value(),self.use_e_jeff_pr_3.isChecked()],[self.om_jeff_alpha_3.value(),self.om_jeff_beta_3.value(),self.use_om_jeff_pr_3.isChecked()], [self.ma_jeff_alpha_3.value(),self.ma_jeff_beta_3.value(),self.use_ma_jeff_pr_3.isChecked()],[self.incl_jeff_alpha_3.value(),self.incl_jeff_beta_3.value(),self.use_incl_jeff_pr_3.isChecked()], [self.Omega_jeff_alpha_3.value(),self.Omega_jeff_beta_3.value(), self.use_Omega_jeff_pr_3.isChecked()],[self.t0_jeff_alpha_3.value(),self.t0_jeff_beta_3.value(), self.use_t0_jeff_pr_3.isChecked()],[self.pl_rad_jeff_alpha_3.value(),self.pl_rad_jeff_beta_3.value(),self.use_pl_rad_jeff_pr_3.isChecked()],[self.a_sol_jeff_alpha_3.value(),self.a_sol_jeff_beta_3.value(),self.use_a_sol_jeff_pr_3.isChecked()],
        [self.K_jeff_alpha_4.value(),self.K_jeff_beta_4.value(),self.use_K_jeff_pr_4.isChecked()],[self.P_jeff_alpha_4.value(),self.P_jeff_beta_4.value(),self.use_P_jeff_pr_4.isChecked()], [self.e_jeff_alpha_4.value(),self.e_jeff_beta_4.value(),self.use_e_jeff_pr_4.isChecked()],[self.om_jeff_alpha_4.value(),self.om_jeff_beta_4.value(),self.use_om_jeff_pr_4.isChecked()], [self.ma_jeff_alpha_4.value(),self.ma_jeff_beta_4.value(),self.use_ma_jeff_pr_4.isChecked()],[self.incl_jeff_alpha_4.value(),self.incl_jeff_beta_4.value(),self.use_incl_jeff_pr_4.isChecked()], [self.Omega_jeff_alpha_4.value(),self.Omega_jeff_beta_4.value(), self.use_Omega_jeff_pr_4.isChecked()],[self.t0_jeff_alpha_4.value(),self.t0_jeff_beta_4.value(), self.use_t0_jeff_pr_4.isChecked()],[self.pl_rad_jeff_alpha_4.value(),self.pl_rad_jeff_beta_4.value(),self.use_pl_rad_jeff_pr_4.isChecked()],[self.a_sol_jeff_alpha_4.value(),self.a_sol_jeff_beta_4.value(),self.use_a_sol_jeff_pr_4.isChecked()],
        [self.K_jeff_alpha_5.value(),self.K_jeff_beta_5.value(),self.use_K_jeff_pr_5.isChecked()],[self.P_jeff_alpha_5.value(),self.P_jeff_beta_5.value(),self.use_P_jeff_pr_5.isChecked()], [self.e_jeff_alpha_5.value(),self.e_jeff_beta_5.value(),self.use_e_jeff_pr_5.isChecked()],[self.om_jeff_alpha_5.value(),self.om_jeff_beta_5.value(),self.use_om_jeff_pr_5.isChecked()], [self.ma_jeff_alpha_5.value(),self.ma_jeff_beta_5.value(),self.use_ma_jeff_pr_5.isChecked()],[self.incl_jeff_alpha_5.value(),self.incl_jeff_beta_5.value(),self.use_incl_jeff_pr_5.isChecked()], [self.Omega_jeff_alpha_5.value(),self.Omega_jeff_beta_5.value(), self.use_Omega_jeff_pr_5.isChecked()],[self.t0_jeff_alpha_5.value(),self.t0_jeff_beta_5.value(), self.use_t0_jeff_pr_5.isChecked()],[self.pl_rad_jeff_alpha_5.value(),self.pl_rad_jeff_beta_5.value(),self.use_pl_rad_jeff_pr_5.isChecked()],[self.a_sol_jeff_alpha_5.value(),self.a_sol_jeff_beta_5.value(),self.use_a_sol_jeff_pr_5.isChecked()],
        [self.K_jeff_alpha_6.value(),self.K_jeff_beta_6.value(),self.use_K_jeff_pr_6.isChecked()],[self.P_jeff_alpha_6.value(),self.P_jeff_beta_6.value(),self.use_P_jeff_pr_6.isChecked()], [self.e_jeff_alpha_6.value(),self.e_jeff_beta_6.value(),self.use_e_jeff_pr_6.isChecked()],[self.om_jeff_alpha_6.value(),self.om_jeff_beta_6.value(),self.use_om_jeff_pr_6.isChecked()], [self.ma_jeff_alpha_6.value(),self.ma_jeff_beta_6.value(),self.use_ma_jeff_pr_6.isChecked()],[self.incl_jeff_alpha_6.value(),self.incl_jeff_beta_6.value(),self.use_incl_jeff_pr_6.isChecked()], [self.Omega_jeff_alpha_6.value(),self.Omega_jeff_beta_6.value(), self.use_Omega_jeff_pr_6.isChecked()],[self.t0_jeff_alpha_6.value(),self.t0_jeff_beta_6.value(), self.use_t0_jeff_pr_6.isChecked()],[self.pl_rad_jeff_alpha_6.value(),self.pl_rad_jeff_beta_6.value(),self.use_pl_rad_jeff_pr_6.isChecked()],[self.a_sol_jeff_alpha_6.value(),self.a_sol_jeff_beta_6.value(),self.use_a_sol_jeff_pr_6.isChecked()],
        [self.K_jeff_alpha_7.value(),self.K_jeff_beta_7.value(),self.use_K_jeff_pr_7.isChecked()],[self.P_jeff_alpha_7.value(),self.P_jeff_beta_7.value(),self.use_P_jeff_pr_7.isChecked()], [self.e_jeff_alpha_7.value(),self.e_jeff_beta_7.value(),self.use_e_jeff_pr_7.isChecked()],[self.om_jeff_alpha_7.value(),self.om_jeff_beta_7.value(),self.use_om_jeff_pr_7.isChecked()], [self.ma_jeff_alpha_7.value(),self.ma_jeff_beta_7.value(),self.use_ma_jeff_pr_7.isChecked()],[self.incl_jeff_alpha_7.value(),self.incl_jeff_beta_7.value(),self.use_incl_jeff_pr_7.isChecked()], [self.Omega_jeff_alpha_7.value(),self.Omega_jeff_beta_7.value(), self.use_Omega_jeff_pr_7.isChecked()],[self.t0_jeff_alpha_7.value(),self.t0_jeff_beta_7.value(), self.use_t0_jeff_pr_7.isChecked()],[self.pl_rad_jeff_alpha_7.value(),self.pl_rad_jeff_beta_7.value(),self.use_pl_rad_jeff_pr_7.isChecked()],[self.a_sol_jeff_alpha_7.value(),self.a_sol_jeff_beta_7.value(),self.use_a_sol_jeff_pr_7.isChecked()],
        [self.K_jeff_alpha_8.value(),self.K_jeff_beta_8.value(),self.use_K_jeff_pr_8.isChecked()],[self.P_jeff_alpha_8.value(),self.P_jeff_beta_8.value(),self.use_P_jeff_pr_8.isChecked()], [self.e_jeff_alpha_8.value(),self.e_jeff_beta_8.value(),self.use_e_jeff_pr_8.isChecked()],[self.om_jeff_alpha_8.value(),self.om_jeff_beta_8.value(),self.use_om_jeff_pr_8.isChecked()], [self.ma_jeff_alpha_8.value(),self.ma_jeff_beta_8.value(),self.use_ma_jeff_pr_8.isChecked()],[self.incl_jeff_alpha_8.value(),self.incl_jeff_beta_8.value(),self.use_incl_jeff_pr_8.isChecked()], [self.Omega_jeff_alpha_8.value(),self.Omega_jeff_beta_8.value(), self.use_Omega_jeff_pr_8.isChecked()],[self.t0_jeff_alpha_8.value(),self.t0_jeff_beta_8.value(), self.use_t0_jeff_pr_8.isChecked()],[self.pl_rad_jeff_alpha_8.value(),self.pl_rad_jeff_beta_8.value(),self.use_pl_rad_jeff_pr_8.isChecked()],[self.a_sol_jeff_alpha_8.value(),self.a_sol_jeff_beta_8.value(),self.use_a_sol_jeff_pr_8.isChecked()],
        [self.K_jeff_alpha_9.value(),self.K_jeff_beta_9.value(),self.use_K_jeff_pr_9.isChecked()],[self.P_jeff_alpha_9.value(),self.P_jeff_beta_9.value(),self.use_P_jeff_pr_9.isChecked()], [self.e_jeff_alpha_9.value(),self.e_jeff_beta_9.value(),self.use_e_jeff_pr_9.isChecked()],[self.om_jeff_alpha_9.value(),self.om_jeff_beta_9.value(),self.use_om_jeff_pr_9.isChecked()], [self.ma_jeff_alpha_9.value(),self.ma_jeff_beta_9.value(),self.use_ma_jeff_pr_9.isChecked()],[self.incl_jeff_alpha_9.value(),self.incl_jeff_beta_9.value(),self.use_incl_jeff_pr_9.isChecked()], [self.Omega_jeff_alpha_9.value(),self.Omega_jeff_beta_9.value(), self.use_Omega_jeff_pr_9.isChecked()],[self.t0_jeff_alpha_9.value(),self.t0_jeff_beta_9.value(), self.use_t0_jeff_pr_9.isChecked()],[self.pl_rad_jeff_alpha_9.value(),self.pl_rad_jeff_beta_9.value(),self.use_pl_rad_jeff_pr_9.isChecked()],[self.a_sol_jeff_alpha_9.value(),self.a_sol_jeff_beta_9.value(),self.use_a_sol_jeff_pr_9.isChecked()],
        ]
 
        for i in range(fit.npl):
            fit.K_jeff_pr[i]  = param_jeff_priors_gui[10*i + 0]    
            fit.P_jeff_pr[i]  = param_jeff_priors_gui[10*i + 1]    
            fit.e_jeff_pr[i]  = param_jeff_priors_gui[10*i + 2]    
            fit.w_jeff_pr[i]  = param_jeff_priors_gui[10*i + 3]    
            fit.M0_jeff_pr[i] = param_jeff_priors_gui[10*i + 4]    
            fit.i_jeff_pr[i]  = param_jeff_priors_gui[10*i + 5]    
            fit.Node_jeff_pr[i] = param_jeff_priors_gui[10*i + 6]    
            fit.t0_jeff_pr[i]   = param_jeff_priors_gui[10*i + 7]
            fit.pl_rad_jeff_pr[i]  = param_jeff_priors_gui[10*i + 8]
            fit.pl_a_jeff_pr[i]    = param_jeff_priors_gui[10*i + 9]

        offset_jeff_priors_gui = [
        [self.RV_Data_jeff_alpha_1.value(),self.RV_Data_jeff_beta_1.value(),self.use_rvoff_jeff_1.isChecked()], 
        [self.RV_Data_jeff_alpha_2.value(),self.RV_Data_jeff_beta_2.value(),self.use_rvoff_jeff_2.isChecked()], 
        [self.RV_Data_jeff_alpha_3.value(),self.RV_Data_jeff_beta_3.value(),self.use_rvoff_jeff_3.isChecked()], 
        [self.RV_Data_jeff_alpha_4.value(),self.RV_Data_jeff_beta_4.value(),self.use_rvoff_jeff_4.isChecked()], 
        [self.RV_Data_jeff_alpha_5.value(),self.RV_Data_jeff_beta_5.value(),self.use_rvoff_jeff_5.isChecked()], 
        [self.RV_Data_jeff_alpha_6.value(),self.RV_Data_jeff_beta_6.value(),self.use_rvoff_jeff_6.isChecked()], 
        [self.RV_Data_jeff_alpha_7.value(),self.RV_Data_jeff_beta_7.value(),self.use_rvoff_jeff_7.isChecked()], 
        [self.RV_Data_jeff_alpha_8.value(),self.RV_Data_jeff_beta_8.value(),self.use_rvoff_jeff_8.isChecked()], 
        [self.RV_Data_jeff_alpha_9.value(),self.RV_Data_jeff_beta_9.value(),self.use_rvoff_jeff_9.isChecked()], 
        [self.RV_Data_jeff_alpha_10.value(),self.RV_Data_jeff_beta_10.value(),self.use_rvoff_jeff_10.isChecked()]
        ]
        
        jitter_jeff_priors_gui = [
        [self.RV_jitter_jeff_alpha_1.value(),self.RV_jitter_jeff_beta_1.value(),self.use_rvjitt_jeff_1.isChecked()], 
        [self.RV_jitter_jeff_alpha_2.value(),self.RV_jitter_jeff_beta_2.value(),self.use_rvjitt_jeff_2.isChecked()], 
        [self.RV_jitter_jeff_alpha_3.value(),self.RV_jitter_jeff_beta_3.value(),self.use_rvjitt_jeff_3.isChecked()], 
        [self.RV_jitter_jeff_alpha_4.value(),self.RV_jitter_jeff_beta_4.value(),self.use_rvjitt_jeff_4.isChecked()], 
        [self.RV_jitter_jeff_alpha_5.value(),self.RV_jitter_jeff_beta_5.value(),self.use_rvjitt_jeff_5.isChecked()], 
        [self.RV_jitter_jeff_alpha_6.value(),self.RV_jitter_jeff_beta_6.value(),self.use_rvjitt_jeff_6.isChecked()],
        [self.RV_jitter_jeff_alpha_7.value(),self.RV_jitter_jeff_beta_7.value(),self.use_rvjitt_jeff_7.isChecked()], 
        [self.RV_jitter_jeff_alpha_8.value(),self.RV_jitter_jeff_beta_8.value(),self.use_rvjitt_jeff_8.isChecked()], 
        [self.RV_jitter_jeff_alpha_9.value(),self.RV_jitter_jeff_beta_9.value(),self.use_rvjitt_jeff_9.isChecked()], 
        [self.RV_jitter_jeff_alpha_10.value(),self.RV_jitter_jeff_beta_10.value(),self.use_rvjitt_jeff_10.isChecked()]   
        ]  
    
    
        for i in range(10): 
            fit.rvoff_jeff_pr[i] = offset_jeff_priors_gui[i]
            fit.jitt_jeff_pr[i]  = jitter_jeff_priors_gui[i] 


        om_dot_jeff_priors_gui = [
        [self.omega_dot_alpha_1.value(),self.omega_dot_beta_1.value(),self.use_omega_dot_jeff_pr_1.isChecked()], [self.omega_dot_alpha_2.value(),self.omega_dot_beta_2.value(),self.use_omega_dot_jeff_pr_2.isChecked()], 
        [self.omega_dot_alpha_3.value(),self.omega_dot_beta_3.value(),self.use_omega_dot_jeff_pr_3.isChecked()], [self.omega_dot_alpha_4.value(),self.omega_dot_beta_4.value(),self.use_omega_dot_jeff_pr_4.isChecked()], 
        [self.omega_dot_alpha_5.value(),self.omega_dot_beta_5.value(),self.use_omega_dot_jeff_pr_5.isChecked()], [self.omega_dot_alpha_6.value(),self.omega_dot_beta_6.value(),self.use_omega_dot_jeff_pr_6.isChecked()], 
        [self.omega_dot_alpha_7.value(),self.omega_dot_beta_7.value(),self.use_omega_dot_jeff_pr_7.isChecked()], [self.omega_dot_alpha_8.value(),self.omega_dot_beta_8.value(),self.use_omega_dot_jeff_pr_8.isChecked()], 
        [self.omega_dot_alpha_9.value(),self.omega_dot_beta_9.value(),self.use_omega_dot_jeff_pr_9.isChecked()]    
        ]  
    
    
        for i in range(9): 
            fit.omega_dot_jeff_pr[i] = om_dot_jeff_priors_gui[i]   
    
 
        fit.rv_lintr_jeff_pr[0]  = [self.lin_trend_jeff_alpha.value(),self.lin_trend_jeff_beta.value(),self.use_lin_tr_jeff_pr.isChecked()]


        offset_jeff_priors_gui_tra = [
        [self.tra_Data_alpha_1.value(),self.tra_Data_beta_1.value(),self.use_traoff_jeff_1.isChecked()], 
        [self.tra_Data_alpha_2.value(),self.tra_Data_beta_2.value(),self.use_traoff_jeff_2.isChecked()], 
        [self.tra_Data_alpha_3.value(),self.tra_Data_beta_3.value(),self.use_traoff_jeff_3.isChecked()], 
        [self.tra_Data_alpha_4.value(),self.tra_Data_beta_4.value(),self.use_traoff_jeff_4.isChecked()], 
        [self.tra_Data_alpha_5.value(),self.tra_Data_beta_5.value(),self.use_traoff_jeff_5.isChecked()],   
        [self.tra_Data_alpha_6.value(),self.tra_Data_beta_6.value(),self.use_traoff_jeff_6.isChecked()], 
        [self.tra_Data_alpha_7.value(),self.tra_Data_beta_7.value(),self.use_traoff_jeff_7.isChecked()], 
        [self.tra_Data_alpha_8.value(),self.tra_Data_beta_8.value(),self.use_traoff_jeff_8.isChecked()], 
        [self.tra_Data_alpha_9.value(),self.tra_Data_beta_9.value(),self.use_traoff_jeff_9.isChecked()], 
        [self.tra_Data_alpha_10.value(),self.tra_Data_beta_10.value(),self.use_traoff_jeff_10.isChecked()]
        ]
        
        jitter_jeff_priors_gui_tra = [
        [self.tra_jitter_alpha_1.value(),self.tra_jitter_beta_1.value(),self.use_trajitt_jeff_1.isChecked()],
        [self.tra_jitter_alpha_2.value(),self.tra_jitter_beta_2.value(),self.use_trajitt_jeff_2.isChecked()], 
        [self.tra_jitter_alpha_3.value(),self.tra_jitter_beta_3.value(),self.use_trajitt_jeff_3.isChecked()], 
        [self.tra_jitter_alpha_4.value(),self.tra_jitter_beta_4.value(),self.use_trajitt_jeff_4.isChecked()], 
        [self.tra_jitter_alpha_5.value(),self.tra_jitter_beta_5.value(),self.use_trajitt_jeff_5.isChecked()],   
        [self.tra_jitter_alpha_6.value(),self.tra_jitter_beta_6.value(),self.use_trajitt_jeff_6.isChecked()], 
        [self.tra_jitter_alpha_7.value(),self.tra_jitter_beta_7.value(),self.use_trajitt_jeff_7.isChecked()], 
        [self.tra_jitter_alpha_8.value(),self.tra_jitter_beta_8.value(),self.use_trajitt_jeff_8.isChecked()], 
        [self.tra_jitter_alpha_9.value(),self.tra_jitter_beta_9.value(),self.use_trajitt_jeff_9.isChecked()], 
        [self.tra_jitter_alpha_10.value(),self.tra_jitter_beta_10.value(),self.use_trajitt_jeff_10.isChecked()]
        ]
        
    
        for i in range(10): 
            fit.tra_off_jeff_pr[i] = offset_jeff_priors_gui_tra[i]
            fit.tra_jitt_jeff_pr[i] = jitter_jeff_priors_gui_tra[i] 

        self.check_RV_GP_priors_jeff()
        self.check_tra_GP_priors_jeff()
    
    

    def check_RV_GP_priors_jeff(self):
        global fit

        GP_rot_jeff_priors_gui = [
        [self.GP_rot_kernel_Amp_jeff_alpha.value(),self.GP_rot_kernel_Amp_jeff_beta.value(),self.use_GP_rot_kernel_Amp_jeff_pr.isChecked()],  
        [self.GP_rot_kernel_time_sc_jeff_alpha.value(),self.GP_rot_kernel_time_sc_jeff_beta.value(),self.use_GP_rot_kernel_time_sc_jeff_pr.isChecked()],  
        [self.GP_rot_kernel_Per_jeff_alpha.value(),self.GP_rot_kernel_Per_jeff_beta.value(),self.use_GP_rot_kernel_Per_jeff_pr.isChecked()],  
        [self.GP_rot_kernel_fact_jeff_alpha.value(),self.GP_rot_kernel_fact_jeff_beta.value(),self.use_GP_rot_kernel_fact_jeff_pr.isChecked()],  
        ]
 
        for i in range(4): 
            fit.GP_rot_jeff_pr[i] = GP_rot_jeff_priors_gui[i]               

        GP_sho_jeff_priors_gui = [
        [self.GP_sho_kernel_S_jeff_alpha.value(),self.GP_sho_kernel_S_jeff_beta.value(), self.use_GP_sho_kernel_S_jeff_pr.isChecked()],  
        [self.GP_sho_kernel_Q_jeff_alpha.value(),self.GP_sho_kernel_Q_jeff_beta.value(), self.use_GP_sho_kernel_Q_jeff_pr.isChecked()],  
        [self.GP_sho_kernel_omega_jeff_alpha.value(),self.GP_sho_kernel_omega_jeff_beta.value(), self.use_GP_sho_kernel_omega_jeff_pr.isChecked()],  
        ]
 
        for i in range(3): 
            fit.GP_sho_jeff_pr[i] = GP_sho_jeff_priors_gui[i]   
            
            
    def check_tra_GP_priors_jeff(self):
        global fit

        tra_GP_rot_jeff_priors_gui = [
        [self.tra_GP_rot_kernel_Amp_jeff_alpha.value(),self.tra_GP_rot_kernel_Amp_jeff_beta.value(),self.use_tra_GP_rot_kernel_Amp_jeff_pr.isChecked()],  
        [self.tra_GP_rot_kernel_time_sc_jeff_alpha.value(),self.tra_GP_rot_kernel_time_sc_jeff_beta.value(),self.use_tra_GP_rot_kernel_time_sc_jeff_pr.isChecked()],  
        [self.tra_GP_rot_kernel_Per_jeff_alpha.value(),self.tra_GP_rot_kernel_Per_jeff_beta.value(),self.use_tra_GP_rot_kernel_Per_jeff_pr.isChecked()],  
        [self.tra_GP_rot_kernel_fact_jeff_alpha.value(),self.tra_GP_rot_kernel_fact_jeff_beta.value(),self.use_tra_GP_rot_kernel_fact_jeff_pr.isChecked()],  
        ]
 
        for i in range(4): 
            fit.tra_GP_rot_jeff_pr[i] = tra_GP_rot_jeff_priors_gui[i]            
    
        tra_GP_sho_jeff_priors_gui = [
        [self.tra_GP_sho_kernel_S_jeff_alpha.value(),self.tra_GP_sho_kernel_S_jeff_beta.value(), self.use_tra_GP_sho_kernel_S_jeff_pr.isChecked()],  
        [self.tra_GP_sho_kernel_Q_jeff_alpha.value(),self.tra_GP_sho_kernel_Q_jeff_beta.value(), self.use_tra_GP_sho_kernel_Q_jeff_pr.isChecked()],  
        [self.tra_GP_sho_kernel_omega_jeff_alpha.value(),self.tra_GP_sho_kernel_omega_jeff_beta.value(), self.use_tra_GP_sho_kernel_omega_jeff_pr.isChecked()],  
        ]
 
        for i in range(3): 
            fit.tra_GP_sho_jeff_pr[i] = tra_GP_sho_jeff_priors_gui[i]   
            
 
 
    def check_arb_pl(self):
        global fit


        fit.arb_st_mass = self.arb_st_mass.value()
        
        arb_param_gui_use = [self.use_arb_Planet_1,self.use_arb_Planet_2,self.use_arb_Planet_3,
                             self.use_arb_Planet_4,self.use_arb_Planet_5,self.use_arb_Planet_6,
                             self.use_arb_Planet_7,self.use_arb_Planet_8,self.use_arb_Planet_9]
 
        
        arb_param_gui = [
                     self.arb_K_1, self.arb_P_1, self.arb_e_1, self.arb_om_1, self.arb_ma_1, self.arb_incl_1, self.arb_Om_1,
                     self.arb_K_2, self.arb_P_2, self.arb_e_2, self.arb_om_2, self.arb_ma_2, self.arb_incl_2, self.arb_Om_2,
                     self.arb_K_3, self.arb_P_3, self.arb_e_3, self.arb_om_3, self.arb_ma_3, self.arb_incl_3, self.arb_Om_3,
                     self.arb_K_4, self.arb_P_4, self.arb_e_4, self.arb_om_4, self.arb_ma_4, self.arb_incl_4, self.arb_Om_4, 
                     self.arb_K_5, self.arb_P_5, self.arb_e_5, self.arb_om_5, self.arb_ma_5, self.arb_incl_5, self.arb_Om_5,
                     self.arb_K_6, self.arb_P_6, self.arb_e_6, self.arb_om_6, self.arb_ma_6, self.arb_incl_6, self.arb_Om_6,
                     self.arb_K_7, self.arb_P_7, self.arb_e_7, self.arb_om_7, self.arb_ma_7, self.arb_incl_7, self.arb_Om_7, 
                     self.arb_K_8, self.arb_P_8, self.arb_e_8, self.arb_om_8, self.arb_ma_8, self.arb_incl_8, self.arb_Om_8,
                     self.arb_K_9, self.arb_P_9, self.arb_e_9, self.arb_om_9, self.arb_ma_9, self.arb_incl_9, self.arb_Om_9,
                     ]
        
        j = 0
        for i in range(9):
            fit.pl_arb_use[i] = arb_param_gui_use[i].isChecked()
            
            if fit.pl_arb_use[i] == True:
                j += 1
            
            fit.e_arb[i]    = arb_param_gui[7*i + 2].value()    
            fit.w_arb[i]    = arb_param_gui[7*i + 3].value()    
            fit.M0_arb[i]   = arb_param_gui[7*i + 4].value()    
            fit.i_arb[i]    = arb_param_gui[7*i + 5].value()    
            fit.Node_arb[i] = arb_param_gui[7*i + 6].value()    
 
            if self.radioButton_KP.isChecked():
                fit.K_arb[i]    = arb_param_gui[7*i + 0].value()    
                fit.P_arb[i]    = arb_param_gui[7*i + 1].value()                 
               # mass_,a_ = rv.mass_a_from_Kepler_fit([fit.K_arb[i],  fit.P_arb[i], fit.e_arb[i],  fit.w_arb[i], fit.M0_arb[i]],1,fit.arb_st_mass)
                mass_,a_ = rv.mass_a_from_Kepler_fit([fit.K_arb[i]/np.sin(np.radians(fit.i_arb[i])),  fit.P_arb[i], fit.e_arb[i],  fit.w_arb[i], fit.M0_arb[i]],1,fit.arb_st_mass)

                fit.mass_arb[i] = float(mass_[0])
                fit.a_arb[i] = float(a_[0])
            else:                
                fit.mass_arb[i] = arb_param_gui[7*i + 0].value()  
                fit.a_arb[i]    = arb_param_gui[7*i + 1].value()  
            
        fit.npl_arb = j# np.count_nonzero(fit.pl_arb_use.values())
 
    

    
####################################################        

    def initialize_color_dialog(self):

        self.colorDialog = QtGui.QColorDialog()
        self.colorDialog.setOption(QtGui.QColorDialog.ShowAlphaChannel, True)
        self.colorDialog.setOption(QtGui.QColorDialog.DontUseNativeDialog, True)


    def initialize_buttons(self):

        # for some reason this does not work!
        #[self.buttonGroup_add_RV_data.setId(bg4, ii) for ii, bg4 in enumerate(self.buttonGroup_add_RV_data.buttons())]
        #[self.buttonGroup_remove_RV_data.setId(bg5, jj) for jj, bg5 in enumerate(self.buttonGroup_remove_RV_data.buttons())]   
        
        self.buttonGroup_add_RV_data.setId(self.Button_RV_data_1,1)
        self.buttonGroup_add_RV_data.setId(self.Button_RV_data_2,2)
        self.buttonGroup_add_RV_data.setId(self.Button_RV_data_3,3)
        self.buttonGroup_add_RV_data.setId(self.Button_RV_data_4,4)
        self.buttonGroup_add_RV_data.setId(self.Button_RV_data_5,5)
        self.buttonGroup_add_RV_data.setId(self.Button_RV_data_6,6)
        self.buttonGroup_add_RV_data.setId(self.Button_RV_data_7,7)
        self.buttonGroup_add_RV_data.setId(self.Button_RV_data_8,8)
        self.buttonGroup_add_RV_data.setId(self.Button_RV_data_9,9)
        self.buttonGroup_add_RV_data.setId(self.Button_RV_data_10,10)

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
        
        self.buttonGroup_color_picker.setId(self.rv_pushButton_color_1,1)
        self.buttonGroup_color_picker.setId(self.rv_pushButton_color_2,2)
        self.buttonGroup_color_picker.setId(self.rv_pushButton_color_3,3)
        self.buttonGroup_color_picker.setId(self.rv_pushButton_color_4,4)
        self.buttonGroup_color_picker.setId(self.rv_pushButton_color_5,5)
        self.buttonGroup_color_picker.setId(self.rv_pushButton_color_6,6)
        self.buttonGroup_color_picker.setId(self.rv_pushButton_color_7,7)
        self.buttonGroup_color_picker.setId(self.rv_pushButton_color_8,8)
        self.buttonGroup_color_picker.setId(self.rv_pushButton_color_9,9)
        self.buttonGroup_color_picker.setId(self.rv_pushButton_color_10,10)
        self.buttonGroup_color_picker.setId(self.rv_pushButton_color_11,11)
        
        self.buttonGroup_color_picker_tra.setId(self.trans_pushButton_color_1,1)
        self.buttonGroup_color_picker_tra.setId(self.trans_pushButton_color_2,2)
        self.buttonGroup_color_picker_tra.setId(self.trans_pushButton_color_3,3)
        self.buttonGroup_color_picker_tra.setId(self.trans_pushButton_color_4,4)
        self.buttonGroup_color_picker_tra.setId(self.trans_pushButton_color_5,5)
        self.buttonGroup_color_picker_tra.setId(self.trans_pushButton_color_6,6)
        self.buttonGroup_color_picker_tra.setId(self.trans_pushButton_color_7,7)
        self.buttonGroup_color_picker_tra.setId(self.trans_pushButton_color_8,8)
        self.buttonGroup_color_picker_tra.setId(self.trans_pushButton_color_9,9)
        self.buttonGroup_color_picker_tra.setId(self.trans_pushButton_color_10,10)
        self.buttonGroup_color_picker_tra.setId(self.trans_pushButton_color_11,11)
     
        self.buttonGroup_symbol_picker.setId(self.rv_pushButton_symbol_1,1)
        self.buttonGroup_symbol_picker.setId(self.rv_pushButton_symbol_2,2)
        self.buttonGroup_symbol_picker.setId(self.rv_pushButton_symbol_3,3)
        self.buttonGroup_symbol_picker.setId(self.rv_pushButton_symbol_4,4)
        self.buttonGroup_symbol_picker.setId(self.rv_pushButton_symbol_5,5)
        self.buttonGroup_symbol_picker.setId(self.rv_pushButton_symbol_6,6)
        self.buttonGroup_symbol_picker.setId(self.rv_pushButton_symbol_7,7)
        self.buttonGroup_symbol_picker.setId(self.rv_pushButton_symbol_8,8)
        self.buttonGroup_symbol_picker.setId(self.rv_pushButton_symbol_9,9)
        self.buttonGroup_symbol_picker.setId(self.rv_pushButton_symbol_10,10)        
        
        self.buttonGroup_symbol_picker_tra.setId(self.trans_pushButton_symbol_1,1)
        self.buttonGroup_symbol_picker_tra.setId(self.trans_pushButton_symbol_2,2)
        self.buttonGroup_symbol_picker_tra.setId(self.trans_pushButton_symbol_3,3)
        self.buttonGroup_symbol_picker_tra.setId(self.trans_pushButton_symbol_4,4)
        self.buttonGroup_symbol_picker_tra.setId(self.trans_pushButton_symbol_5,5)
        self.buttonGroup_symbol_picker_tra.setId(self.trans_pushButton_symbol_6,6)
        self.buttonGroup_symbol_picker_tra.setId(self.trans_pushButton_symbol_7,7)
        self.buttonGroup_symbol_picker_tra.setId(self.trans_pushButton_symbol_8,8)
        self.buttonGroup_symbol_picker_tra.setId(self.trans_pushButton_symbol_9,9)
        self.buttonGroup_symbol_picker_tra.setId(self.trans_pushButton_symbol_10,10)      
        
        
        
        self.buttonGroup_use_planets = QtWidgets.QButtonGroup()
        self.buttonGroup_use_planets.setExclusive(False)
        
        self.buttonGroup_use_planets.addButton(self.use_Planet1,1)
        self.buttonGroup_use_planets.addButton(self.use_Planet2,2)
        self.buttonGroup_use_planets.addButton(self.use_Planet3,3)
        self.buttonGroup_use_planets.addButton(self.use_Planet4,4)
        self.buttonGroup_use_planets.addButton(self.use_Planet5,5)
        self.buttonGroup_use_planets.addButton(self.use_Planet6,6)
        self.buttonGroup_use_planets.addButton(self.use_Planet7,7)
        self.buttonGroup_use_planets.addButton(self.use_Planet8,8)
        self.buttonGroup_use_planets.addButton(self.use_Planet9,9)
        
       # self.buttonGroup_use_planets.setId(self.use_Planet1,1)
       # self.buttonGroup_use_planets.setId(self.use_Planet2,2)
       # self.buttonGroup_use_planets.setId(self.use_Planet3,3)      
       # self.buttonGroup_use_planets.setId(self.use_Planet4,4)
       # self.buttonGroup_use_planets.setId(self.use_Planet5,5)
       # self.buttonGroup_use_planets.setId(self.use_Planet6,6)
       # self.buttonGroup_use_planets.setId(self.use_Planet7,7)
       # self.buttonGroup_use_planets.setId(self.use_Planet8,8)
       # self.buttonGroup_use_planets.setId(self.use_Planet9,9)
        
        
        self.colors_gls.setFont(self.font) 
        self.colors_gls_o_c.setFont(self.font) 

      
        
    def initialize_plots(self):

        global p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,pe,pdi,pcor

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
        p15 = self.graphicsView_orb_evol_elements_om   
       
        p16 = self.graphicsView_orbital_view
        
        p17 = self.graphicsView_orb_evol_res_dom
        p18 = self.graphicsView_orb_evol_res_theta
        
        p19 = self.graphicsView_orb_evol_elements_i
        p20 = self.graphicsView_orb_evol_energy
        
        pe  = self.graphicsView_extra_plot
        
        pdi = self.load_data_plot

        pcor = self.graphicsView_corner

        xaxis = ['BJD [days]','BJD [days]','BJD [days]','BJD [days]','BJD [days]','x','period [d]','period [d]','period [d]','period [d]','period [d]','period [d]','t [yr]','t [yr]','t [yr]','a [au]','t [yr]','t [yr]','t [yr]','t [yr]','','x','x']
        yaxis = ['RV [m/s]','RV [m/s]','Rel. Flux','Rel. Flux','y','y','power','power','SDE','SDE','power','power','a [au]','e','omega [deg]','a [au]','delta omega [deg]','theta [deg]','inc [deg]','energy','','y','y']       
        xunit = ['' ,'','','','','','','','','','','','','','','','','','','','','','']
        yunit = ['' ,'' , '','','','','','','','','','','','','','','','','','','','','']

        zzz = [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,pe,pdi,pcor]
 
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

        p16.getViewBox().setAspectLocked(True)

        return   
        

        
    def identify_power_peaks(self,x,y,sig_level=np.array([]), power_level=np.array([]) ):
 
        per_ind = argrelextrema(y, np.greater)
        per_x   = x[per_ind]
        per_y   = y[per_ind]     

        peaks_sort = sorted(range(len(per_y)), key=lambda k: per_y[k], reverse=True)

        per_x   = per_x[peaks_sort]   
        per_y   = per_y[peaks_sort]  
        
        peaks_pos = [per_x,per_y]
        
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
 
        if len(per_x)  < 10:
            max_peaks = len(per_x) 
        else:
            max_peaks = 10

        for j in range(max_peaks):
            text_peaks = text_peaks +"""
period = %.2f [d], power = %.4f"""%(per_x[j],per_y[j])  
            if sig_level.size != 0 and per_y[j] > sig_level[-1]:
                text_peaks = text_peaks +"""  significant"""
                
        ################################################        
    
        return text_peaks , peaks_pos 
        
  
        
######################## Correlation plots ###################################### 
        
    def init_correlations_combo(self):
        global fit
        self.comboBox_corr_1.clear()
        self.comboBox_corr_2.clear()
        
        self.initialize_corr_y = {k: [] for k in range(20)}
        z = 0 

        if fit.filelist.ndset != 0:

            for i in range(max(fit.filelist.idset)+1):
                self.comboBox_corr_1.addItem('RV %s'%(i+1),i+1) 
                self.comboBox_corr_2.addItem('RV %s'%(i+1),i+1) 
            
                self.initialize_corr_y[z] = np.array([fit.fit_results.rv_model.jd[fit.filelist.idset==i],
                                                      fit.fit_results.rv_model.rvs[fit.filelist.idset==i], 
                                                      fit.fit_results.rv_model.rv_err[fit.filelist.idset==i]])  
                z +=1
                
            for i in range(max(fit.filelist.idset)+1):
                self.comboBox_corr_1.addItem('RV o-c %s'%(i+1),i+1)         
                self.comboBox_corr_2.addItem('RV o-c %s'%(i+1),i+1)  
    
                self.initialize_corr_y[z] = np.array([fit.fit_results.rv_model.jd[fit.filelist.idset==i],
                                                      fit.fit_results.rv_model.o_c[fit.filelist.idset==i], 
                                                      fit.fit_results.rv_model.rv_err[fit.filelist.idset==i]]) 
                z +=1                          
         

        for i in range(0,10,1):         
            if len(fit.act_data_sets[i]) != 0: 
                self.comboBox_corr_1.addItem('act. data %s'%(i+1),i+1)       
                self.comboBox_corr_2.addItem('act. data %s'%(i+1),i+1) 
                
                self.initialize_corr_y[z] = fit.act_data_sets[i] 
                z +=1                                         
                
        #return                
                
                
    def update_correlations_data_plots(self):
        global fit, colors,  p6 
        
        ind1 = self.comboBox_corr_1.currentIndex()
        ind2 = self.comboBox_corr_2.currentIndex()
 
        p6.plot(clear=True,)  
        
        self.color_corr.setStyleSheet("color: %s;"%colors[0]) 

        #p6.autoRange()     
        
        if not ind1 == None and not ind2 == None:

            
            if len(self.initialize_corr_y[ind1][0]) == len(self.initialize_corr_y[ind2][0]):
                p6.plot(self.initialize_corr_y[ind1][1],self.initialize_corr_y[ind2][1], pen=None,symbol='o',
                #symbolPen=,
                symbolSize=self.act_data_size.value(),enableAutoRange=True,viewRect=True,
                symbolBrush=colors[0]
                )    
                
                
                pears = pearsonr(self.initialize_corr_y[ind1][1],self.initialize_corr_y[ind2][1] )
                
                if pears[0] < 0:
                    pos_neg = "negative"
                else:
                    pos_neg = "positive"
                    
                if abs(pears[0]) < 0.3:                      
                    strong_mod_weak = "very weak"
                elif 0.3 <= abs(pears[0]) <= 0.5: 
                     strong_mod_weak = "weak"
                elif 0.5 <= abs(pears[0]) <= 0.7: 
                     strong_mod_weak = "moderate"                       
                elif 0.7 <= abs(pears[0]) <= 1: 
                     strong_mod_weak = "strong"  
                else:
                     strong_mod_weak = "n/a"  
                     
                m, c = np.polyfit(self.initialize_corr_y[ind1][1],self.initialize_corr_y[ind2][1], 1,
                                  w=1/np.sqrt(self.initialize_corr_y[ind1][2]**2 + self.initialize_corr_y[ind2][2]**2),
                                  full=False,cov=True)  
                
                e = np.sqrt(np.diag(c))
                
                text = '''Pearson's correlation coefficient 2-tailed p-value: 
%s, %s
(A %s %s correlation)

Polyfit coefficients: 
%s +/- %s, 
%s +/- %s   

'''%(pears[0],pears[1], pos_neg, strong_mod_weak, m[0],e[0], m[1],e[1])

                self.corr_print_info.clicked.connect(lambda: self.print_info_for_object(text))

                
                if self.plot_corr_err.isChecked():
                    err1 = pg.ErrorBarItem(x=self.initialize_corr_y[ind1][1], y=self.initialize_corr_y[ind2][1],symbol='o', 
                    top=self.initialize_corr_y[ind2][2],bottom=self.initialize_corr_y[ind2][2], 
                    left=self.initialize_corr_y[ind1][2],right=self.initialize_corr_y[ind1][2],                     
                    beam=0.0, pen=colors[0])  

                    p6.addItem(err1)   
                    
                if self.plot_corr_coef.isChecked():

                    p6.plot(self.initialize_corr_y[ind1][1], self.initialize_corr_y[ind1][1]*m[0] +m[1] , pen='k')                            

                
                p6.autoRange()
            

                return   
            
            else:               
                text_err = pg.TextItem('Not the same time series!',color=(0,0,0))#, anchor=(0,0), border='w',color) #, fill=(0, 0, 255, 100))
                p6.addItem(text_err, ignoreBounds=True)   
        else:   
             
                return

    def get_corr_color(self):
        global fit
        
        colorz = self.colorDialog.getColor()
        colors[0]=colorz.name()   

        self.update_correlations_data_plots()

    def corr_plot_x_labels(self):
        global fit
        
        text, okPressed = QtGui.QInputDialog.getText(self, "x-axis label","(No special characters!)", QtGui.QLineEdit.Normal, "")
        
        if okPressed and text != '':
            p6.setLabel('bottom', '%s'%text, units='',  **{'font-size':'9pt'})
 
        else:
            return
    
        self.update_correlations_data_plots()
 

    def corr_plot_y_labels(self):
        global fit
        
        text, okPressed = QtGui.QInputDialog.getText(self, "y-axis label","(No special characters!)", QtGui.QLineEdit.Normal, "")
        
        if okPressed and text != '':
            p6.setLabel('left', '%s'%text, units='',  **{'font-size':'9pt'})
 
        else:
            return
    
        self.update_correlations_data_plots()
              
        
        
######################## Activity plots ######################################  
                
    def init_activity_combo(self):
        global fit
        
        for i in range(10):
            self.comboBox_act_data_gls.addItem('act. data %s'%(i+1),i+1)       
            self.comboBox_act_data.addItem('act. data %s'%(i+1),i+1)       
                
        
 
    def update_activity_gls_plots(self,ind):
        global fit, colors,  p11 
 
        omega = 1/ np.logspace(np.log(self.gls_min_period.value()), np.log(self.gls_max_period.value()), num=self.gls_n_omega.value())
        power_levels = np.array([self.gls_fap1.value(),self.gls_fap2.value(),self.gls_fap3.value()])
  
        if len(fit.act_data_sets[ind]) != 0 and len(fit.act_data_sets[ind][0]) > 5:

            p11.plot(clear=True,)        
 
            act_per = gls.Gls((fit.act_data_sets[ind][0], fit.act_data_sets[ind][1],fit.act_data_sets[ind][2]), 
            fast=True,  verbose=False, norm= "ZK",ofac=self.gls_ofac.value(), fbeg=omega[-1], fend=omega[ 0],)
            
            ######################## GLS ##############################
            if self.radioButton_act_GLS_period.isChecked():
                p11.setLogMode(True,False)        
                p11.plot(1/act_per.freq, act_per.power,pen=fit.colors[ind],symbol=None ) 
                p11.setLabel('bottom', 'period [d]', units='',  **{'font-size':'9pt'}) 

            else:
                p11.setLogMode(False,False)        
                p11.plot(act_per.freq, act_per.power,pen=fit.colors[ind],symbol=None )                    
                p11.setLabel('bottom', 'frequency [1/d]', units='',  **{'font-size':'9pt'}) 

                                               
            [p11.addLine(x=None, y=fap, pen=pg.mkPen('k', width=0.8, style=QtCore.Qt.DotLine)) for ii,fap in enumerate(act_per.powerLevel(np.array(power_levels)))]
  
            text_peaks, pos_peaks = self.identify_power_peaks(1/act_per.freq, act_per.power, power_level = power_levels, sig_level = act_per.powerLevel(np.array(power_levels)) )    
    
            self.act_periodogram_print_info.clicked.connect(lambda: self.print_info_for_object(
            act_per.info(stdout=False) + text_peaks ))   
    
            return
        else:   
            p11.plot(clear=True,)        

            return



    def update_activity_data_plots(self,ind):
        global fit, colors,  p5 

        if len(fit.act_data_sets[ind]) != 0:

            p5.plot(clear=True,)  

            err1 = pg.ErrorBarItem(x=fit.act_data_sets[ind][0], y=fit.act_data_sets[ind][1],symbol='o', 
           # height=fit.act_data_sets[ind][2], beam=0.0, pen=fit.colors[ind])  
            top=fit.act_data_sets[ind][2],
            bottom=fit.act_data_sets[ind][2],           
            beam=0.0, pen=fit.colors[ind])


            p5.addItem(err1)      
            p5.addLine(x=None, y=0, pen=pg.mkPen('#ff9933', width=0.8))

            p5.plot(fit.act_data_sets[ind][0],fit.act_data_sets[ind][1], pen=None,symbol='o',
            symbolSize=self.act_data_size.value(),enableAutoRange=True,viewRect=True,
            symbolBrush=fit.colors[ind]            
            )      

            p5.setLabel('left', 'y', units='',  **{'font-size':'9pt'})     

            return
        else:   
            p5.plot(clear=True,)        

            return  



######################## SciPy setup ######################################        

    def init_scipy_combo(self):    
        global fit 

        for i in range(len(fit.SciPy_min)):
            self.comboBox_scipy_minimizer_1.addItem('%s'%(fit.SciPy_min[i]),i) 
            self.comboBox_scipy_minimizer_2.addItem('%s'%(fit.SciPy_min[i]),i) 
           
        self.comboBox_scipy_minimizer_1.setCurrentIndex(6)
        self.comboBox_scipy_minimizer_2.setCurrentIndex(0)
            
    def check_scipy_min(self):    
        global fit             
            
        ind_min_1 = self.comboBox_scipy_minimizer_1.currentIndex()
        ind_min_2 = self.comboBox_scipy_minimizer_2.currentIndex()
       

        fit.SciPy_min_use_1 = fit.SciPy_min[ind_min_1]
        fit.SciPy_min_use_2 = fit.SciPy_min[ind_min_2]
        fit.SciPy_min_N_use_1 = int(self.scipy_N_consecutive_iter_1.value())
        fit.SciPy_min_N_use_2 = int(self.scipy_N_consecutive_iter_2.value())
        
        
        fit.Simplex_opt    = {'disp': True, 'maxiter': int(self.simplex_maxiter.value()), 'return_all': False, 'maxfev': int(self.simplex_maxfev.value()), 'xtol':self.simplex_xtol.value() , 'ftol': self.simplex_ftol.value() ,'adaptive':True }
        fit.Powell_opt     = {'disp': True, 'return_all': False, 'maxiter': int(self.powell_maxiter.value()), 'direc': None, 'func': None, 'maxfev': int(self.powell_maxfev.value()), 'xtol': self.powell_xtol.value(), 'ftol': self.powell_ftol.value()}
        fit.CG_opt         = {'disp': True, 'gtol': self.cg_gtol.value(), 'eps': 1.4901161193847656e-08, 'return_all': False, 'maxiter': int(self.cg_maxiter.value()), 'norm': np.inf}
        fit.BFGS_opt       = {'disp': True, 'gtol': self.bfgs_gtol.value(), 'eps': 1.4901161193847656e-08, 'return_all': False, 'maxiter': int(self.bfgs_maxiter.value()), 'norm': np.inf}
        fit.Newton_cg_opt  = {'disp': True, 'xtol': self.Newton_cg_xtol.value(), 'eps': 1.4901161193847656e-08, 'return_all': False, 'maxiter': int(self.Newton_cg_maxiter.value())} 
        fit.L_BFGS_B_opt   = {'disp': True, 'maxcor': int(self.LBFGSB_maxcor.value()), 'ftol': 2.220446049250313e-09, 'gtol': self.LBFGSB_gtol.value(), 'eps': 1e-08, 'maxfun': int(self.LBFGSB_maxiter.value()), 'maxiter': int(self.LBFGSB_maxiter.value()), 'iprint': -1, 'maxls': 20}    
        fit.TNC_opt        = {'disp': True, 'eps': self.TNC_eps.value(), 'scale': None, 'offset': None, 'mesg_num': None, 'maxCGit': int(self.TNC_maxcgit.value()), 'maxiter': int(self.TNC_maxiter.value()), 'eta': self.TNC_eta.value(), 'stepmx':self.TNC_stepmx.value(), 'accuracy': self.TNC_accuracy.value(), 'minfev': self.TNC_minfev.value(), 'ftol': self.TNC_ftol.value(), 'xtol':self.TNC_ftol.value(), 'gtol': self.TNC_gtol.value(), 'rescale': -1 }  
       # fit.COBYLA_opt     = {'disp': True, 'rhobeg': self.cobyla_rhobeg.value(), 'maxiter':  int(self.cobyla_maxiter.value()), 'catol': self.cobyla_catol.value() }
        fit.SLSQP_opt      = {'disp': True, 'maxiter': int(self.slsqp_maxiter.value()),  'eps': 1.4901161193847656e-08, 'ftol': self.slsqp_ftol.value(), 'iprint': 1}
      
        
    def cross_hair(self, plot_wg, log=False ):
        global fit 

        vLine = pg.InfiniteLine(angle=90, movable=False)#, pos=0)
        hLine = pg.InfiniteLine(angle=0,  movable=False)#, pos=2450000)
        plot_wg.addItem(vLine, ignoreBounds=True)
        plot_wg.addItem(hLine, ignoreBounds=True)
        label = pg.TextItem(anchor=(1, 1))
        plot_wg.addItem(label, ignoreBounds=True) 
        
        vb = plot_wg.getViewBox()   
       # viewrange = vb.viewRange()

        def mouseMoved(evt):
            pos = evt[0]  ## using signal proxy turns original arguments into a tuple
            if plot_wg.sceneBoundingRect().contains(pos):             

                mousePoint = vb.mapSceneToView(pos)

                if log == True:
                    label.setText("x=%0.3f,  y=%0.3f"%(10**mousePoint.x(), mousePoint.y()))
                else:
                    label.setText("x=%0.3f,  y=%0.3f"%(mousePoint.x(), mousePoint.y()))
                    
                vLine.setPos(mousePoint.x())
                hLine.setPos(mousePoint.y())
                #print(mousePoint.x(),mousePoint.y())
                label.setPos(mousePoint.x(), mousePoint.y())
        plot_wg.getViewBox().setAutoVisible(y=True)

        proxy = pg.SignalProxy(plot_wg.scene().sigMouseMoved, rateLimit=60, slot=mouseMoved)    
        plot_wg.proxy = proxy
  
   
            
    def label_peaks(self, plot_wg2, pos_peaks, GLS = True, o_c = False):
    
        if GLS == True and self.avoid_GLS_RV_alias.isChecked():
            x_peaks = pos_peaks[0][pos_peaks[0]>1.2]
            y_peaks = pos_peaks[1][pos_peaks[0]>1.2]
        else:            
            x_peaks = pos_peaks[0] 
            y_peaks = pos_peaks[1]          
            
            
        if GLS == True:
            N_peaks = int(self.N_GLS_peak_to_point.value())
            if o_c == True:
                log = self.radioButton_RV_o_c_GLS_period.isChecked()
            else:
                log = self.radioButton_RV_GLS_period.isChecked()
                
            type_per = "GLS"
        else:
            N_peaks = int(self.N_TLS_peak_to_point.value())
            log = False
            type_per = "TLS"

        if len(x_peaks) <  N_peaks:
            N_peaks = len(x_peaks)              
            print("You have reached the maximum number of %s peaks."%type_per)

        for i in range(N_peaks):

            text_arrow = pg.TextItem("test", anchor=(0.5,1.9))
            
            if log == True:
                arrow = pg.ArrowItem(pos=(np.log10(x_peaks[i]), y_peaks[i]), angle=270)
                text_arrow.setText('%0.2f d' % (x_peaks[i]))               
                text_arrow.setPos(np.log10(x_peaks[i]),y_peaks[i])
                
            elif log == False and GLS == True:   
                arrow = pg.ArrowItem(pos=(1/x_peaks[i], y_peaks[i]), angle=270)
                text_arrow.setText('%0.2f d' % (x_peaks[i]))               
                text_arrow.setPos(1/x_peaks[i],y_peaks[i])     
                
            elif log == False and GLS == False: 
                arrow = pg.ArrowItem(pos=(x_peaks[i], y_peaks[i]), angle=270)
                text_arrow.setText('%0.2f d' % (x_peaks[i]))               
                text_arrow.setPos(x_peaks[i],y_peaks[i])                    
  
            plot_wg2.addItem(arrow) 
            plot_wg2.addItem(text_arrow)        
        
        
        
######################## RV plots ######################################        
      
 
    def init_gls_norm_combo(self):    
        global fit
        
        self.norms = ['ZK',  'HorneBaliunas', 'Cumming', 'wrms', 'chisq', 'lnL', 'dlnL']
        #'Scargle',
        for i in range(len(self.norms)):
            self.gls_norm_combo.addItem('%s'%(self.norms[i]),i+1)   
            
         
    def run_gls(self):
        global fit
                
        omega = 1/ np.logspace(np.log10(self.gls_min_period.value()), np.log10(self.gls_max_period.value()), num=int(self.gls_n_omega.value()))
        ind_norm = self.gls_norm_combo.currentIndex()

        if len(fit.fit_results.rv_model.jd) > 5:      
            RV_per = gls.Gls((fit.fit_results.rv_model.jd, fit.fit_results.rv_model.rvs, fit.fit_results.rv_model.rv_err), 
            fast=True,  verbose=False, norm=self.norms[ind_norm],ofac=self.gls_ofac.value(), fbeg=omega[-1], fend=omega[0],)
            
            fit.gls = RV_per
        else:
            return
        
        self.update_RV_GLS_plots()
        
        
    def run_gls_o_c(self):
        global fit
                        
        omega = 1/ np.logspace(np.log10(self.gls_min_period.value()), np.log10(self.gls_max_period.value()), num=int(self.gls_n_omega.value()))
        ind_norm = self.gls_norm_combo.currentIndex()
 
        if len(fit.fit_results.rv_model.jd) > 5:
            RV_per_res = gls.Gls((fit.fit_results.rv_model.jd, fit.fit_results.rv_model.o_c, fit.fit_results.rv_model.rv_err), 
            fast=True,  verbose=False, norm= self.norms[ind_norm],ofac=self.gls_ofac.value(), fbeg=omega[-1], fend=omega[ 0],)            
    
            fit.gls_o_c = RV_per_res        
        else:
            return
        
        self.update_RV_o_c_GLS_plots()  
        
        

    def get_RV_GLS_plot_color(self):
        global fit
        
        colorz = self.colorDialog.getColor()
        fit.gls_colors[0]=colorz.name()   

        self.update_RV_GLS_plots() 
        
    def get_RV_o_c_GLS_plot_color(self):
        global fit
        
        colorz = self.colorDialog.getColor()
        fit.gls_colors[1]=colorz.name()   
         
        self.update_RV_o_c_GLS_plots()  
        

    def update_RV_GLS_plots(self):
        global fit, p7 
 
        p7.plot(clear=True,)   
        
        self.colors_gls.setStyleSheet("color: %s;"%fit.gls_colors[0])               
                          
        power_levels = np.array([self.gls_fap1.value(),self.gls_fap2.value(),self.gls_fap3.value()])
        gls_model_width = float(self.gls_model_width.value())
    
        if len(fit.fit_results.rv_model.jd) > 5:

            ######################## GLS ##############################
            if self.radioButton_RV_GLS_period.isChecked():
                p7.setLogMode(True,False)        
                p7.plot(1/fit.gls.freq, fit.gls.power,pen={'color': fit.gls_colors[0], 'width': gls_model_width},symbol=None ) 

                p7.setLabel('bottom', 'period [d]', units='',  **{'font-size':'9pt'})    
                
            else:
                p7.setLogMode(False,False)        
                p7.plot(fit.gls.freq, fit.gls.power,pen={'color': fit.gls_colors[0], 'width': self.gls_model_width.value()},symbol=None )                
                p7.setLabel('bottom', 'frequency [1/d]', units='',  **{'font-size':'9pt'}) 
                
                
            if fit.gls.norm == 'ZK':
                [p7.addLine(x=None, y=fap, pen=pg.mkPen('k', width=0.8, style=QtCore.Qt.DotLine)) for ii,fap in enumerate(fit.gls.powerLevel(np.array(power_levels)))]
 
            text_peaks, pos_peaks = self.identify_power_peaks(1/fit.gls.freq, fit.gls.power, power_level = power_levels, sig_level = fit.gls.powerLevel(np.array(power_levels)) )   

            self.label_peaks(p7, pos_peaks, GLS = True, o_c = False)

            self.RV_periodogram_print_info.clicked.connect(lambda: self.print_info_for_object(
            fit.gls.info(stdout=False) + text_peaks   ))   


        if self.gls_cross_hair.isChecked():
            self.cross_hair(p7,log=self.radioButton_RV_GLS_period.isChecked())    
 
 
    def update_RV_o_c_GLS_plots(self):
        global fit,  p8  
 
        p8.plot(clear=True,)  

        self.colors_gls_o_c.setStyleSheet("color: %s;"%fit.gls_colors[1]) 
        gls_o_c_model_width = float(self.gls_o_c_model_width.value())
        
        power_levels = np.array([self.gls_fap1.value(),self.gls_fap2.value(),self.gls_fap3.value()])

        if len(fit.fit_results.rv_model.jd) > 5:
 

            ######################## GLS o-c ##############################
            if self.radioButton_RV_o_c_GLS_period.isChecked():
                p8.setLogMode(True,False)        
                p8.plot(1/fit.gls_o_c.freq, fit.gls_o_c.power, pen={'color': fit.gls_colors[1], 'width': gls_o_c_model_width},symbol=None ) 
                p8.setLabel('bottom', 'period [d]', units='',  **{'font-size':'9pt'})
 
            else:
                p8.setLogMode(False,False)        
                p8.plot(fit.gls_o_c.freq, fit.gls_o_c.power, pen={'color': fit.gls_colors[1], 'width': gls_o_c_model_width},symbol=None )   
                p8.setLabel('bottom', 'frequency [1/d]', units='',  **{'font-size':'9pt'})                
                
            if fit.gls_o_c.norm == 'ZK':
                [p8.addLine(x=None, y=fap, pen=pg.mkPen('k', width=0.8, style=QtCore.Qt.DotLine)) for ii,fap in enumerate(fit.gls_o_c.powerLevel(np.array(power_levels)))]            

            text_peaks, pos_peaks = self.identify_power_peaks(1/fit.gls_o_c.freq, fit.gls_o_c.power, power_level = power_levels, sig_level = fit.gls_o_c.powerLevel(np.array(power_levels)) )

            self.label_peaks(p8, pos_peaks, GLS=True, o_c = True)
 
            self.RV_res_periodogram_print_info.clicked.connect(lambda: self.print_info_for_object(fit.gls_o_c.info(stdout=False)+ text_peaks  )  )      
 
        if self.gls_o_c_cross_hair.isChecked():
            self.cross_hair(p8,log=self.radioButton_RV_o_c_GLS_period.isChecked())    
    
    def update_WF_plots(self):
        global fit, p12  
 
        p12.plot(clear=True,) 
        p12.setLogMode(True,False)
                        
        omega = 1/ np.logspace(np.log10(self.gls_min_period.value()), np.log10(self.gls_max_period.value()),  num=int(self.gls_n_omega.value()))
        #power_levels = np.array([0.1,0.01,0.001])
        
        if len(fit.fit_results.rv_model.jd) > 5:
            ######################## DFT (Window) ##############################
            WF_power = []
            for omi in 2*np.pi*omega: 
                phase = (fit.fit_results.rv_model.jd-fit.fit_results.rv_model.jd[0]) * omi                 
                WC = np.sum(np.cos(phase))
                WS = np.sum(np.sin(phase))
                WF_power.append((WC**2 + WS**2)/len(fit.fit_results.rv_model.jd)**2) 

            WF_power = np.array(WF_power)
            ######################## GLS o-c ##############################
            if self.radioButton_RV_WF_period.isChecked():
                p12.setLogMode(True,False)        
                p12.plot(1/np.array(omega), WF_power,pen='k',symbol=None )   
                p12.setLabel('bottom', 'period [d]', units='',  **{'font-size':'9pt'})
            else:
                p12.setLogMode(False,False)        
                p12.plot(np.array(omega), WF_power,pen='k',symbol=None )   
                p12.setLabel('bottom', 'frequency [1/d]', units='',  **{'font-size':'9pt'})      

            text_peaks, pos_peaks = self.identify_power_peaks(1/np.array(omega), WF_power)

                        
            self.WF_print_info.clicked.connect(lambda: self.print_info_for_object(text_peaks))        


    def rv_GP_set_use(self):

        if self.do_RV_GP.isChecked():
            fit.doGP = True
        else:
            fit.doGP = False
            
    def tra_GP_set_use(self):

        if self.do_tra_GP.isChecked():
            fit.tra_doGP = True
        else:
            fit.tra_doGP = False            
            
            

    def update_RV_plots(self):
        global fit, p1,p2
 
        p1.plot(clear=True,)
        p2.plot(clear=True,)
 
    
        self.check_RV_symbol_sizes()

        if len(fit.filelist.idset)==0:
            return
 
        p1.addLine(x=None, y=0,   pen=pg.mkPen('#ff9933', width=0.8))
  
        if fit.doGP == True:
            #rv.get_gps_model(self) 
            y_model = fit.fit_results.model + fit.gp_model_curve[0]
            y_model_o_c = fit.gp_model_curve[0]
        else:
            y_model = fit.fit_results.model 
            y_model_o_c = np.zeros(len(y_model))
            
            
        model_curve = p1.plot(fit.fit_results.model_jd,y_model, 
        pen={'color': fit.colors[-1], 'width': self.rv_model_width.value()},enableAutoRange=True, #symbolPen={'color': 0.5, 'width': 0.1}, symbolSize=1,symbol='o',
        viewRect=True, labels =  {'left':'RV', 'bottom':'JD'}) 
        
        model_curve.setZValue(self.RV_model_z.value()) 
        
        
        if  fit.doGP == True:
            pfill = pg.FillBetweenItem(p1.plot(fit.fit_results.model_jd, fit.fit_results.model + fit.gp_model_curve[0]+fit.gp_model_curve[2]), 
                                       p1.plot(fit.fit_results.model_jd, fit.fit_results.model + fit.gp_model_curve[0]-fit.gp_model_curve[2]), 
                                       brush = pg.mkColor(244,140,66,128))
            p1.addItem(pfill)  



        if self.jitter_to_plots.isChecked() and not self.split_jitter.isChecked():
            error_list = self.add_jitter(fit.fit_results.rv_model.rv_err, fit.filelist.idset)
        elif self.jitter_to_plots.isChecked() and self.split_jitter.isChecked():
            error_list = fit.fit_results.rv_model.rv_err
            error_list2 = self.add_jitter(fit.fit_results.rv_model.rv_err, fit.filelist.idset)            
        else:
            error_list = fit.fit_results.rv_model.rv_err
            
            
        for i in range(max(fit.filelist.idset)+1):
            p1.plot(fit.fit_results.rv_model.jd[fit.filelist.idset==i],fit.fit_results.rv_model.rvs[fit.filelist.idset==i], 
            pen=None, #{'color': colors[i], 'width': 1.1},
            symbol=fit.pyqt_symbols_rvs[i],
            symbolPen={'color': fit.colors[i], 'width': 1.1},
            symbolSize=fit.pyqt_symbols_size_rvs[i],enableAutoRange=True,viewRect=True,
            symbolBrush=fit.colors[i]
            )    
              
            err1 = pg.ErrorBarItem(x=fit.fit_results.rv_model.jd[fit.filelist.idset==i], 
                                   y=fit.fit_results.rv_model.rvs[fit.filelist.idset==i],symbol='o', 
            #height=error_list[fit.filelist.idset==i],
            top=error_list[fit.filelist.idset==i],
            bottom=error_list[fit.filelist.idset==i],           
            beam=0.0, pen=fit.colors[i])  
 
            p1.addItem(err1)
            
            
            if self.jitter_to_plots.isChecked() and self.split_jitter.isChecked():
            
                err1a = pg.ErrorBarItem(x=fit.fit_results.rv_model.jd[fit.filelist.idset==i], 
                                       y=fit.fit_results.rv_model.rvs[fit.filelist.idset==i],symbol='o', 
                #height=error_list[fit.filelist.idset==i],
                top=error_list2[fit.filelist.idset==i],
                bottom=error_list2[fit.filelist.idset==i],           
                beam=0.0, pen='#000000')  
                err1a.setZValue(-10)
                p1.addItem(err1a)            
            
            
            
            
            
        if self.RV_plot_cross_hair.isChecked():
            self.cross_hair(p1,log=False)  
            
            
            
        p2.addLine(x=None, y=0, pen=pg.mkPen('#ff9933', width=0.8))

        p2.plot(fit.fit_results.model_jd,y_model_o_c, 
        pen={'color':  fit.colors[-1], 'width': self.rv_model_width.value()},enableAutoRange=True, #symbolPen={'color': 0.5, 'width': 0.1}, symbolSize=1,symbol='o',
        viewRect=True, labels =  {'left':'RV', 'bottom':'JD'}) 
        
        if fit.doGP == True:
            pfill_o_c = pg.FillBetweenItem(p2.plot(fit.fit_results.model_jd, fit.gp_model_curve[0]+fit.gp_model_curve[2]), 
                                           p2.plot(fit.fit_results.model_jd, fit.gp_model_curve[0]-fit.gp_model_curve[2]), 
                                           brush = pg.mkColor(244,140,66,128))
            p2.addItem(pfill_o_c)  
        
        
        for i in range(max(fit.filelist.idset)+1):
            p2.plot(fit.fit_results.rv_model.jd[fit.filelist.idset==i],fit.fit_results.rv_model.o_c[fit.filelist.idset==i], 
            pen=None, #{'color': colors[i], 'width': 1.1},
            symbol=fit.pyqt_symbols_rvs[i],
            symbolPen={'color': fit.colors[i], 'width': 1.1},
            symbolSize=fit.pyqt_symbols_size_rvs[i],enableAutoRange=True,viewRect=True,
            symbolBrush=fit.colors[i]
            )        
            err2 = pg.ErrorBarItem(x=fit.fit_results.rv_model.jd[fit.filelist.idset==i], 
                                   y=fit.fit_results.rv_model.o_c[fit.filelist.idset==i],symbol='o', 
            #height=error_list[fit.filelist.idset==i],
            top=error_list[fit.filelist.idset==i],
            bottom=error_list[fit.filelist.idset==i],           
            beam=0.0, pen=fit.colors[i])  
            
            p2.addItem(err2)  
            
            
            if self.jitter_to_plots.isChecked() and self.split_jitter.isChecked():
            
                err2a = pg.ErrorBarItem(x=fit.fit_results.rv_model.jd[fit.filelist.idset==i], 
                                       y=fit.fit_results.rv_model.o_c[fit.filelist.idset==i],symbol='o', 
                #height=error_list[fit.filelist.idset==i],
                top=error_list2[fit.filelist.idset==i],
                bottom=error_list2[fit.filelist.idset==i],           
                beam=0.0, pen='#000000')  
                err2a.setZValue(-10)
                p2.addItem(err2a)            
                        
            
            
 
        if self.RV_o_c_plot_cross_hair.isChecked():
            self.cross_hair(p2,log=False)       
            
            
        
    def update_plots(self):
        self.update_RV_GLS_plots()
        self.update_RV_o_c_GLS_plots()     
        
        self.update_WF_plots()                
        self.update_RV_plots()
        self.update_extra_plots()
        self.update_orb_plot()
        #self.change_extra_plot()
        self.update_transit_plots()    
        

################################ RV files #######################################################
        
    def showDialog_fortran_input_file(self):
        global fit, ses_list
 
        input_files = QtGui.QFileDialog.getOpenFileName(self, 'Open session', '', 'Data (*.init)')
        
        if str(input_files[0]) != '':
            fit_new=rv.signal_fit(str(input_files[0]), 'RVmod session',readinputfile=True)
            
            if len(ses_list) == 1:
                ses_list[0] = fit_new
                fit = fit_new
            else:
                ses_list.append(fit_new)
                
            self.session_list()
            self.update_use_from_input_file()
            self.init_fit()
            self.update_RV_file_buttons()
            
            
    def showDialog_RVbank_input_file(self):
        global fit, ses_list
 
        input_files = QtGui.QFileDialog.getOpenFileName(self, 'Open RVBank data', '', 'All (*.*);;Data (*.csv)')

        if str(input_files[0]) != '':
 
            fit.add_RVbank_dataset(self.file_from_path(input_files[0]), str(input_files[0]))

            ##################
            self.init_fit()            
            self.update_use_from_input_file()            
            self.update_use()
            self.update_params()
            self.update_RV_file_buttons()
            self.update_act_file_buttons()
            #self.update_activity_gls_plots(0)
            #self.buttonGroup_activity_data.button(but_ind).setText(self.file_from_path(input_files[0]))          
            

    def showDialog_RV_input_file(self):
        global fit

        but_ind = self.buttonGroup_add_RV_data.checkedId()   
        input_files = QtGui.QFileDialog.getOpenFileName(self, 'Open RV data', '', 'All (*.*);;Data (*.vels)')

        if str(input_files[0]) != '':
 
            fit.add_dataset(self.file_from_path(input_files[0]), str(input_files[0]),0.0,1.0)
            #### new stuf ####
            fit.add_rv_dataset('test', str(input_files[0]),rv_idset =but_ind-1)
            ##################
            self.init_fit()            
            self.update_use_from_input_file()            
            self.update_use()
            self.update_params()
            self.update_RV_file_buttons()
            


    def remove_RV_file(self):
        global fit

        but_ind = self.buttonGroup_remove_RV_data.checkedId()   
        
       # try:
       #     dirname, basename = os.path.split(fit.filelist.files[but_ind-1].path)
       #     os.system('rm -r %s'%dirname)
       # except:
       #     return
        
        fit.remove_dataset(but_ind -1)
        
        #### new stuf ####
        fit.remove_rv_dataset(but_ind -1)
        #### new stuf ####
           
        self.init_fit()         
        self.update_use_from_input_file()   
        self.update_use()
        self.update_gui_params()
        self.update_params()
        self.update_RV_file_buttons()

    def update_RV_file_buttons(self):
        global fit, colors    
        
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        #font.setWeight(75)
        
        for i in range(10):
            if i < fit.filelist.ndset:
                self.buttonGroup_add_RV_data.button(i+1).setStyleSheet("color: %s;"%fit.colors[i])
                self.buttonGroup_remove_RV_data.button(i+1).setStyleSheet("color: %s;"%fit.colors[i])
                font.setPointSize(9)
                self.buttonGroup_add_RV_data.button(i+1).setText(fit.filelist.files[i].name) 
                self.buttonGroup_add_RV_data.button(i+1).setFont(font)

            else:
                self.buttonGroup_add_RV_data.button(i+1).setStyleSheet("")
                self.buttonGroup_remove_RV_data.button(i+1).setStyleSheet("")
                self.buttonGroup_add_RV_data.button(i+1).setText("data %s"%(i+1))

                #"background-color: #333399;""background-color: yellow;" "selection-color: yellow;"  "selection-background-color: blue;")               
        self.init_correlations_combo()


################################ transit files #######################################################      

    def showDialog_tra_input_file(self):
        global fit

        but_ind = self.buttonGroup_transit_data.checkedId()   
        input_files = QtGui.QFileDialog.getOpenFileName(self, 'Open Transit data', '', 'All (*.*);;Data (*.tran)')

        if str(input_files[0]) != '':
 
            fit.add_transit_dataset('test', str(input_files[0]),tra_idset =but_ind-1)
              
            self.update_use_from_input_file()            
            self.update_use()
            self.update_gui_params()
          
            self.init_fit()            

            
            self.update_params()
            self.update_tra_file_buttons()
            self.buttonGroup_transit_data.button(but_ind).setText(self.file_from_path(input_files[0]))
            
            self.tab_timeseries_RV.setCurrentWidget(self.Phot_timeseries_plot)

            
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
                self.buttonGroup_transit_data.button(i+1).setStyleSheet("color: %s;"%fit.colors[i])
                self.buttonGroup_remove_transit_data.button(i+1).setStyleSheet("color: %s;"%fit.colors[i])
            else:
                self.buttonGroup_transit_data.button(i+1).setStyleSheet("")
                self.buttonGroup_remove_transit_data.button(i+1).setStyleSheet("")
                self.buttonGroup_transit_data.button(i+1).setText("data %s"%(i+1))

                #"background-color: #333399;""background-color: yellow;" "selection-color: yellow;"  "selection-background-color: blue;")               
        self.update_transit_plots()
 

################################ activity files #######################################################
        
    def showDialog_act_input_file(self):
        global fit

        but_ind = self.buttonGroup_activity_data.checkedId()   
        input_files = QtGui.QFileDialog.getOpenFileName(self, 'Open Activity data', '', 'Data (*.act);;All (*.*)')

        if str(input_files[0]) != '':
 
            fit.add_act_dataset('test', str(input_files[0]),act_idset =but_ind-1)
            #self.init_fit()            
            #self.update_use_from_input_file()            
            #self.update_use()
            #self.update_params()
            self.update_act_file_buttons()
            self.update_activity_gls_plots(but_ind-1)
            self.buttonGroup_activity_data.button(but_ind).setText(self.file_from_path(input_files[0]))

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
                self.buttonGroup_activity_data.button(i+1).setStyleSheet("color: %s;"%fit.colors[i])
                self.buttonGroup_remove_activity_data.button(i+1).setStyleSheet("color: %s;"%fit.colors[i])
            else:
                self.buttonGroup_activity_data.button(i+1).setStyleSheet("")
                self.buttonGroup_remove_activity_data.button(i+1).setStyleSheet("")
                self.buttonGroup_activity_data.button(i+1).setText("data %s"%(i+1))

                #"background-color: #333399;""background-color: yellow;" "selection-color: yellow;"  "selection-background-color: blue;")               
        self.init_correlations_combo()


##################################### Various ################################# 


    def init_fit(self): 
        global fit
        
        # A hack when .ses are imported from the Example directory... TBFixed
        fit.cwd = os.getcwd()

        minimize_fortran=True
        if fit.model_saved == False or len(fit.fit_results.rv_model.jd) != len(fit.filelist.idset):
            fit.fitting(fileinput=False,outputfiles=[1,1,1], minimize_fortran=minimize_fortran,  fortran_kill=self.dyn_model_to_kill.value(), timeout_sec=self.master_timeout.value(), minimize_loglik=True,amoeba_starts=0, print_stat=False, eps=self.dyn_model_accuracy.value(), dt=self.time_step_model.value(), npoints=self.points_to_draw_model.value(), model_max= self.model_max_range.value(), model_min= self.model_min_range.value())
            #self.worker_RV_fitting(, ff=0, m_ln=True, auto_fit = False , init = True ):
            #self.fit_dispatcher(init=False)
            for i in range(fit.npl):
                 rv.phase_RV_planet_signal(fit,i+1)       
                 
        self.update_labels()
        self.update_gui_params()
        self.update_errors() 
        self.update_a_mass() 
        
        self.run_gls()
        self.run_gls_o_c()
        
        self.update_plots() 
        self.update_transit_plots() 
        self.plot_evol_all()

        self.jupiter_push_vars() 

        
 
    def add_jitter(self, errors, ind):
        global fit
        
        errors_with_jitt = np.array([np.sqrt(errors[i]**2 + fit.params.jitters[ii]**2)  for i,ii in enumerate(ind)])
        return errors_with_jitt




################ Extra Plots (work in progress) ######################
         
    
    def update_extra_plots(self):
        global fit

        self.comboBox_extra_plot.clear()

        if fit.npl != 0:
            for i in range(fit.npl):
                self.comboBox_extra_plot.addItem('phase pl %s'%(i+1),i+1)
            
            self.comboBox_extra_plot.addItem('RV GLS',fit.npl+1)
            self.comboBox_extra_plot.addItem('RV GLS o-c',fit.npl+2)

            self.phase_plots(1)  
            
        elif fit.filelist.ndset != 0:
            self.comboBox_extra_plot.addItem('RV GLS',fit.npl+1)
            self.comboBox_extra_plot.addItem('RV GLS o-c',fit.npl+2)     
            
            self.extra_RV_GLS_plots()

        self.comboBox_extra_plot.activated.connect(self.handleActivated)        


 
    def handleActivated(self, index):
        global fit, pe 
        
        ind = self.comboBox_extra_plot.itemData(index) 

 
        if ind <= fit.npl:
            self.phase_plots(ind)
        elif ind == fit.npl+1: 
            self.extra_RV_GLS_plots()
        elif ind == fit.npl+2: 
            self.extra_RV_GLS_o_c_plots()            
            #pe.setYLink(p2)
           # pe.setXLink(p2)

            #gg = p2.getPlotItem().getViewBox()
         #   hh = gg.getViewBox()
         #   pe.scene().addItem(gg)   
         #   pe.scene().addItem(gg)
        else:
            return

    def phase_plots(self, ind, offset = 0):
        global fit, colors   
        
        pe.plot(clear=True,)    
        pe.setLogMode(False,False)        
        
        ######## TBF #############
        if self.radioButton_transit.isChecked():
            return
        ########################
    
        ph_data = fit.ph_data[ind-1]
        ph_model = fit.ph_model[ind-1]


        offset = (self.RV_phase_slider.value()/100.0)* fit.params.planet_params[7*(ind-1)+1] 

        if len(ph_data) == 1:
            return

        if self.jitter_to_plots.isChecked() and len(ph_data) != 0 and not self.split_jitter.isChecked() :
            error_list = self.add_jitter(ph_data[2], ph_data[3])
        elif self.jitter_to_plots.isChecked() and len(ph_data) != 0 and self.split_jitter.isChecked() :
            error_list = ph_data[2]
            error_list2 = self.add_jitter(ph_data[2], ph_data[3])
        else:
            if len(ph_data) != 0:
                error_list = ph_data[2]
            else:
                return
 
            
        #rv_data = ph_data[0]   
       # if fit.doGP == True:
        #    rv_data = fit.gp_model_data[0]
        #else:
        #    rv_data = ph_data[1]
        rv_data = ph_data[1]

        model_time_phase = np.array((ph_model[0]-offset)%fit.params.planet_params[7*(ind-1)+1] )
                             
        sort = sorted(range(len(model_time_phase)), key=lambda k: model_time_phase[k])                        
        model_time_phase  = model_time_phase[sort] 
        ph_model =  ph_model[1][sort] 
         
        pe.addLine(x=None, y=0, pen=pg.mkPen('#ff9933', width=0.8))  
                        
        model_curve = pe.plot(model_time_phase,ph_model, pen={'color':  fit.colors[-1], 'width': self.rv_model_width.value()+1},
        enableAutoRange=True,viewRect=True, labels =  {'left':'RV', 'bottom':'JD'})   
 
        model_curve.setZValue(self.RV_model_z.value())        
        
        
        for i in range(max(ph_data[3])+1):
        
            pe.plot((ph_data[0][ph_data[3]==i]-offset)%fit.params.planet_params[7*(ind-1)+1],rv_data[ph_data[3]==i],             
            pen=None, #{'color': colors[i], 'width': 1.1},
            symbol=fit.pyqt_symbols_rvs[i],
            symbolPen={'color': fit.colors[i], 'width': 1.1},
            symbolSize=fit.pyqt_symbols_size_rvs[i],enableAutoRange=True,viewRect=True,
            symbolBrush=fit.colors[i]
            )  
               
            err_ = pg.ErrorBarItem(x=(ph_data[0][ph_data[3]==i]-offset)%fit.params.planet_params[7*(ind-1)+1], y=rv_data[ph_data[3]==i],
            symbol=fit.pyqt_symbols_rvs[i], 
            top=error_list[ph_data[3]==i],
            bottom=error_list[ph_data[3]==i],           
            beam=0.0, pen=fit.colors[i]) 
            
            pe.addItem(err_)
            
            if self.jitter_to_plots.isChecked() and self.split_jitter.isChecked():
 
                err_2 = pg.ErrorBarItem(x=(ph_data[0][ph_data[3]==i]-offset)%fit.params.planet_params[7*(ind-1)+1], y=rv_data[ph_data[3]==i],
                symbol=fit.pyqt_symbols_rvs[i], 
                top=error_list2[ph_data[3]==i],
                bottom=error_list2[ph_data[3]==i],           
                beam=0.0, pen='#000000')  
                err_2.setZValue(-10)
                pe.addItem(err_2) 
                
            
            
        
        pe.setLabel('bottom', 'phase [days]', units='',  **{'font-size':'9pt'})
        pe.setLabel('left',   'RV [m/s]', units='',  **{'font-size':'9pt'})  

 
        if self.extra_plot_cross_hair.isChecked():
            self.cross_hair(pe,log=False)   



    ############### VERY VERY VERY Ugly fix !!!! 
    
    def extra_RV_GLS_plots(self):
        global fit,  pe 
 
        pe.plot(clear=True,)   
        power_levels = np.array([self.gls_fap1.value(),self.gls_fap2.value(),self.gls_fap3.value()])

 
        ######################## GLS ##############################
        if self.radioButton_RV_GLS_period.isChecked():
            pe.setLogMode(True,False)        
            pe.plot(1/fit.gls.freq, fit.gls.power, pen='r',symbol=None ) 
            pe.setLabel('bottom', 'period [d]', units='',  **{'font-size':'9pt'})    
            pe.setLabel('left', 'Power', units='',  **{'font-size':'9pt'})    
           
        else:
            pe.setLogMode(False,False)        
            pe.plot(fit.gls.freq, fit.gls.power, pen='r',symbol=None )                    
            pe.setLabel('bottom', 'frequency [1/d]', units='',  **{'font-size':'9pt'}) 
            pe.setLabel('left', 'Power', units='',  **{'font-size':'9pt'})    

        if fit.gls.norm == 'ZK':
            [pe.addLine(x=None, y=fap, pen=pg.mkPen('k', width=0.8, style=QtCore.Qt.DotLine)) for ii,fap in enumerate(fit.gls.powerLevel(np.array(power_levels)))]


        if self.extra_plot_cross_hair.isChecked():
            self.cross_hair(pe,log=self.radioButton_RV_GLS_period.isChecked())   

    def extra_RV_GLS_o_c_plots(self):
        global fit,  pe 
 
        pe.plot(clear=True,)           
        power_levels = np.array([self.gls_fap1.value(),self.gls_fap2.value(),self.gls_fap3.value()])
        
        
        ######################## GLS o-c ##############################
        if self.radioButton_RV_o_c_GLS_period.isChecked():
            pe.setLogMode(True,False)        
            pe.plot(1/fit.gls_o_c.freq, fit.gls_o_c.power, pen='r',symbol=None ) 
            pe.setLabel('bottom', 'period [d]', units='',  **{'font-size':'9pt'})    
            pe.setLabel('left', 'Power', units='',  **{'font-size':'9pt'})    
           
        else:
            pe.setLogMode(False,False)        
            pe.plot(fit.gls_o_c.freq, fit.gls_o_c.power, pen='r',symbol=None )                    
            pe.setLabel('bottom', 'frequency [1/d]', units='',  **{'font-size':'9pt'}) 
            pe.setLabel('left', 'Power', units='',  **{'font-size':'9pt'})    


        if fit.gls.norm == 'ZK':
            [pe.addLine(x=None, y=fap, pen=pg.mkPen('k', width=0.8, style=QtCore.Qt.DotLine)) for ii,fap in enumerate(fit.gls_o_c.powerLevel(np.array(power_levels)))]

        if self.extra_plot_cross_hair.isChecked():
            self.cross_hair(pe,log=self.radioButton_RV_o_c_GLS_period.isChecked())   



############ TLS ##############################      
       
    def worker_tls_complete(self, resid = False):
        global fit  
 
        if resid == False:     
            self.update_tls_plots() 
        else:    
            self.update_tls_o_c_plots() 
                 
        self.statusBar().showMessage('')   
 
        self.jupiter_push_vars()   
        self.calc_TLS.setEnabled(True)         
        self.calc_TLS_o_c.setEnabled(True)  
 
    def worker_tls(self, resid = False):
        global fit  
        
        if tls_not_found==True:
            print("TLS Not found, try to install with 'pip install transitleastsquares'") 
            return


        self.calc_TLS.setEnabled(False)         
        self.calc_TLS_o_c.setEnabled(False)  
        
        z=0
        for i in range(10):
            if len(fit.tra_data_sets[i]) != 0:
                z=z+1
        
        if z <= 0:
            choice = QtGui.QMessageBox.information(self, 'Warning!',
            "Not possible to look for planets if there are no transit data loaded. Please add your transit data first. Okay?", QtGui.QMessageBox.Ok)      
            self.calc_TLS.setEnabled(True)         
            return   

        self.statusBar().showMessage('Looking for Transit events (TLS).... ')                 
        worker_tls_ = Worker(lambda:  self.tls_search(resid = resid) )# Any other args, kwargs are passed to the run  
 
        worker_tls_.signals.finished.connect(lambda:  self.worker_tls_complete(resid = resid))
        
        self.tabWidget_helper.setCurrentWidget(self.tab_info)

        
        # worker.signals.result.connect(self.print_output)
        #worker.signals.finished.connect(self.thread_complete)
       # worker.signals.progress.connect(self.progress_fn)
        self.threadpool.start(worker_tls_)       
     

    def tls_search(self, resid = False):
        global fit
        
        if resid == True:
            lc_data = fit.tra_data_sets[0][3]
        else:
            lc_data = fit.tra_data_sets[0][1]
            
        
        tls_model = transitleastsquares(fit.tra_data_sets[0][0], lc_data)
        tls_results = tls_model.power(oversampling_factor=int(self.tls_ofac.value()), duration_grid_step=self.tls_grid_step.value())
    
        if resid == True:
            fit.tls_o_c = tls_results  # TB Fixed with an rvmod object (i.e. fit.tls_obj)
        else:
            fit.tls = tls_results  # TB Fixed with an rvmod object (i.e. fit.tls_obj)



    def update_tls_plots(self): 
        global fit, p9, colors

        if len(fit.tls) == 0:
            return
    
        p9.plot(clear=True,) 
        
        if self.tls_cross_hair.isChecked():
            self.cross_hair(p9,log=False)      
            
        if len(fit.tra_data_sets[0]) != 0:

            text = '''
Best results from TLS:
 
Period: %s d   
Transit depth: %s 
Transit duration: %s d
'''%(fit.tls.period,fit.tls.depth,fit.tls.duration)
           
            p9.plot(fit.tls.periods, fit.tls.power,        
            pen='r',  enableAutoRange=True,viewRect=True)
#0.9      5.7
#0.95     6.1
#0.99     7.0
#0.999    8.3
#0.9999   9.1
            [p9.addLine(x=None, y=fap, pen=pg.mkPen('k', width=0.8, style=QtCore.Qt.DotLine)) for ii,fap in enumerate(np.array([5.7,7.0,8.3]))]

            text_peaks, pos_peaks = self.identify_power_peaks(fit.tls.periods,fit.tls.power,  sig_level = np.array([5.7,7.0,8.3] )   )

            self.label_peaks(p9, pos_peaks, GLS = False)
 
            self.tls_print_info.clicked.connect(lambda: self.print_info_for_object(text + text_peaks))   
            
            return

        else:    
            text_err = pg.TextItem('Nothing to plot',color=(0,0,0))#, anchor=(0,0), border='w',color) #, fill=(0, 0, 255, 100))
            p9.addItem(text_err, ignoreBounds=True)   
            self.tls_print_info.clicked.connect(lambda: self.print_info_for_object(""))            
            return
     
        
    def update_tls_o_c_plots(self): 
        global fit, p10, colors
    
        if len(fit.tls_o_c) == 0:
            return
    
        p10.plot(clear=True,) 

            
        if self.tls_o_c_cross_hair.isChecked():
            self.cross_hair(p10,log=False) 
            
        if len(fit.tra_data_sets[0]) != 0:
            #t = fit.tra_data_sets[0][0]
            #flux = fit.tra_data_sets[0][1]
           # flux_err = fit.tra_data_sets[0][2]
           
            text = '''
Best results from TLS:
          
Period: %s d   
Transit depth: %s 
Transit duration: %s d
'''%(fit.tls_o_c.period,fit.tls_o_c.depth,fit.tls_o_c.duration)
           
            p10.plot(fit.tls_o_c.periods, fit.tls_o_c.power,        
            pen='r',  enableAutoRange=True,viewRect=True)
#0.9      5.7
#0.95     6.1
#0.99     7.0
#0.999    8.3
#0.9999   9.1
            [p10.addLine(x=None, y=fap, pen=pg.mkPen('k', width=0.8, style=QtCore.Qt.DotLine)) for ii,fap in enumerate(np.array([5.7,7.0,8.3]))]
   
            text_peaks, pos_peaks = self.identify_power_peaks(fit.tls_o_c.periods,fit.tls_o_c.power,  sig_level = np.array([5.7,7.0,8.3] )   )
 
            self.label_peaks(p10, pos_peaks, GLS = False)
            
            self.tls_o_c_print_info.clicked.connect(lambda: self.print_info_for_object(text + text_peaks))          
            
            return

        else:    
            text_err = pg.TextItem('Nothing to plot',color=(0,0,0))#, anchor=(0,0), border='w',color) #, fill=(0, 0, 255, 100))
            p10.addItem(text_err, ignoreBounds=True)   
            self.tls_o_c_print_info.clicked.connect(lambda: self.print_info_for_object(""))            
            return


        
 ############ transit fitting (Work in progress here) ##############################      
       
    def worker_transit_fitting_complete(self):
        global fit  
 
   
       
        self.update_labels()
        self.update_gui_params()
        self.update_errors() 
        self.update_a_mass()                 
        self.update_transit_plots()  

        fit=rv.get_xyz(fit)
                         
        self.statusBar().showMessage('')  
        
        if fit.type_fit["RV"] == True:
            for i in range(fit.npl):
                rv.phase_RV_planet_signal(fit,i+1) 
            self.run_gls()
            self.run_gls_o_c()                
            self.update_plots()  
            
            
        self.jupiter_push_vars()   
        self.button_fit.setEnabled(True)         
 
    def worker_transit_fitting(self, ff=1, auto_fit = False ):
        global fit  
        
        self.button_fit.setEnabled(False)         
        self.update_params() 
        self.update_use()   
        
        # check if transit data is present
        z=0
        for i in range(10):
            if len(fit.tra_data_sets[i]) != 0:
                z=z+1
        
        if z <= 0:
            choice = QtGui.QMessageBox.information(self, 'Warning!',
            "Not possible to look for planets if there are no transit data loaded. Please add your transit data first. Okay?", QtGui.QMessageBox.Ok)      
            self.button_fit.setEnabled(True)         
            return 
        
        if fit.type_fit["RV"] == True:
             if fit.filelist.ndset <= 0:
                 choice = QtGui.QMessageBox.information(self, 'Warning!',
                 "Not possible to look for planets if there are no RV data loaded. Please add your RV data first. Okay?", QtGui.QMessageBox.Ok)      
                 self.button_fit.setEnabled(True)         
                 return   

        if fit.type_fit["RV"] == True :        
            self.statusBar().showMessage('Minimizing Transit + RV parameters.... SciPy in action, please be patient.  ')       
        else:
            self.statusBar().showMessage('Minimizing Transit parameters.... SciPy in action, please be patient. ')       
           

        self.set_tra_ld()            
        self.check_bounds()
        self.check_priors_nr()   
        self.check_priors_jeff()   

        self.check_scipy_min()
        fit.model_npoints = self.points_to_draw_model.value()

          
        worker4 = Worker(lambda:  self.transit_fit(ff=ff ) )# Any other args, kwargs are passed to the run  
 
        worker4.signals.finished.connect(self.worker_transit_fitting_complete)
        
        # worker.signals.result.connect(self.print_output)
        #worker.signals.finished.connect(self.thread_complete)
       # worker.signals.progress.connect(self.progress_fn)
        self.threadpool.start(worker4)       
     

    def transit_fit(self, ff=0 ):
        global fit
        
        if ff ==0:        
            fit.init_fit = True
        else:
            fit.init_fit = False

        # this is only a simple hack.. This junk will be removed later on
        if ff ==3034:
            old_t0_use = dill.copy(fit.t0_use)
            old_pl_a_use = dill.copy(fit.pl_a_use)
            old_pl_rad_use = dill.copy(fit.pl_rad_use)
            old_rv_use = dill.copy(fit.use.use_planet_params)
            old_rvoff_use = dill.copy(fit.use.use_offsets)
            old_rvjitt_use = dill.copy(fit.use.use_jitters)
            old_tra_off_use = dill.copy(fit.tra_off)
            old_tra_jitt_use = dill.copy(fit.tra_jitt)
            
            for i in range(fit.npl):
                fit.t0_use[i] = False
                fit.pl_a_use[i] = False
                fit.pl_rad_use[i] = False 
                fit.use.use_planet_params[i*7] = False 
                fit.use.use_planet_params[i*7+1] = False 
                fit.use.use_planet_params[i*7+2] = False 
                fit.use.use_planet_params[i*7+3] = False 
                fit.use.use_planet_params[i*7+4] = False 
                fit.use.use_planet_params[i*7+5] = False 
                fit.use.use_planet_params[i*7+6] = False 
            for i in range(10): 
                fit.use.use_jitters[i] = False
                fit.use.use_offsets[i] = False   
                fit.tra_off_use[i] = False
                fit.tra_jitt_use[i] = False
            #old_tra_use = fit.tr_params_use 
            #fit.tr_params_use = [False, False,False,False,False,False,False]
            #rv.run_SciPyOp_transit(fit)
            rv.run_SciPyOp(fit)
            
            for i in range(fit.npl):
                fit.t0_use[i] = old_t0_use[i]
                fit.pl_a_use[i] = old_pl_a_use[i]
                fit.pl_rad_use[i] = old_pl_rad_use[i]             
                fit.use.use_planet_params[i*7] = dill.copy(old_rv_use[i*7])  
                fit.use.use_planet_params[i*7+1] = dill.copy(old_rv_use[i*7+1])  
                fit.use.use_planet_params[i*7+2] = dill.copy(old_rv_use[i*7+2]) 
                fit.use.use_planet_params[i*7+3] = dill.copy(old_rv_use[i*7+3]) 
                fit.use.use_planet_params[i*7+4] = dill.copy(old_rv_use[i*7+4])  
                fit.use.use_planet_params[i*7+5] = dill.copy(old_rv_use[i*7+5])  
                fit.use.use_planet_params[i*7+6] = dill.copy(old_rv_use[i*7+6]) 
            for i in range(10): 
                fit.use.use_jitters[i] =  dill.copy(old_rvjitt_use[i]) 
                fit.use.use_offsets[i] =  dill.copy(old_rvoff_use[i])  
                fit.tra_off_use[i] = dill.copy(old_tra_off_use[i])
                fit.tra_jitt_use[i] = dill.copy(old_tra_jitt_use[i])   

                
        else:

            rv.run_SciPyOp(fit)

#### Transit plots ################ 
    def update_transit_plots(self): 
        global fit, p3, colors
    
        p3.plot(clear=True,) 
        p4.plot(clear=True,)         
           
        self.check_tra_symbol_sizes()

        tr_files = []
        
        for i in range(10):
            if len(fit.tra_data_sets[i]) != 0:
                tr_files.append(fit.tra_data_sets[i])
        
        for j in range(len(tr_files)):        
        
        #if len(fit.tra_data_sets[0]) != 0:
            t = np.array(tr_files[j][0])
            flux = np.array(tr_files[j][1] + fit.tra_off[j])
            flux_err = np.sqrt(tr_files[j][2]**2 + fit.tra_jitt[j]**2)
            
            
            
            fit.prepare_for_mcmc(rtg = fit.rtg)    
            par = np.array(fit.parameters)  

            flux_model = np.ones(len(flux))
            m =  {k: [] for k in range(9)}
             
            
            #### a quick fix, TBD! ########
            if fit.rtg[1]:
                if fit.gp_kernel == 'RotKernel':
                    rv_gp_npar = 4
                if fit.gp_kernel == 'SHOKernel':
                    rv_gp_npar = 3
                #fit.gps = []
            else:
                rv_gp_npar = 0   
                
            fit.tr_params.limb_dark = str(fit.ld_m[j])      #limb darkening model       
            fit.tr_params.u = fit.ld_u[j]
            
            for i in range(fit.npl):


                if fit.hkl == True:
                    fit.tr_params.ecc = np.sqrt(par[fit.filelist.ndset*2 +7*i+2]**2 + par[fit.filelist.ndset*2 +7*i+3]**2)
                    fit.tr_params.w  = np.degrees(np.arctan2(par[fit.filelist.ndset*2 +7*i+2],par[fit.filelist.ndset*2 +7*i+3]))%360
                else:
                    fit.tr_params.ecc = par[fit.filelist.ndset*2 +7*i+2] #0.0  
                    fit.tr_params.w   = par[fit.filelist.ndset*2 +7*i+3] #90
                
                fit.tr_params.per = par[fit.filelist.ndset*2 +7*i+1] #1.0    #orbital period
                fit.tr_params.inc = par[fit.filelist.ndset*2 +7*i+5]#90. #orbital inclination (in degrees)
                    
                fit.tr_params.t0  = par[fit.filelist.ndset*2  +7*fit.npl +1+rv_gp_npar + 3*i]                
                fit.tr_params.a   = par[fit.filelist.ndset*2  +7*fit.npl +1+rv_gp_npar + 3*i+1] #15  #semi-major axis (in units of stellar radii)
                fit.tr_params.rp  = par[fit.filelist.ndset*2  +7*fit.npl +1+rv_gp_npar + 3*i+2] #0.15   #planet radius (in units of stellar radii)
        
                m[i] = batman.TransitModel(fit.tr_params, t)    #initializes model
     
                flux_model = flux_model * m[i].light_curve(fit.tr_params)     
                
                ############### Phase signal TBD this should not be here! ####################################
                
                if self.plot_phase_pholded_tran.isChecked() and fit.tra_doGP != True:
                    data_time_phase = np.array( (t  - t[0]- fit.tr_params.per/2.0)% fit.tr_params.per  )  
                 
                    sort = np.array(sorted(range(len(data_time_phase)), key=lambda k: data_time_phase[k])    )                    
                     
                    t      = data_time_phase[sort] 
                    flux          = flux[sort] 
                    flux_err      = flux_err[sort]  
                    flux_model    = flux_model[sort] 
                                       
                    fit.ph_data_tra[i] = [data_time_phase[sort] ,flux[sort], flux_err[sort]]
                    fit.ph_model_tra[i] = [data_time_phase[sort] ,flux_model[sort]]
                
                    p3.setLabel('bottom', 'phase [days]', units='',  **{'font-size':'9pt'})
                else:
                    p3.setLabel('bottom', 'BJD [days]', units='',  **{'font-size':'9pt'})
                    
               
                
            tr_o_c = flux -flux_model     
            ######## TBD this should not be here!
            fit.tra_data_sets[j][3] = tr_o_c + 1
            fit.tra_data_sets[j][4] = tr_o_c 

                
            if fit.tra_doGP == True:
                y_model = flux_model + fit.tra_gp_model_curve[0]
                y_model_o_c = fit.tra_gp_model_curve[0]
            else:
                y_model = flux_model 
                y_model_o_c = np.zeros(len(flux_model))
                

            
            p3.plot(t, flux,        
            pen=None,  
            symbol=fit.pyqt_symbols_tra[i],
            symbolPen={'color': fit.tra_colors[j], 'width': 1.1},
            symbolSize=fit.pyqt_symbols_size_tra[i],enableAutoRange=True,viewRect=True,
            symbolBrush=fit.tra_colors[j] ) 
            
            err_ = pg.ErrorBarItem(x=t, y=flux, symbol='o',
                                  # height=flux_err, 
                                   top=flux_err, 
                                   bottom=flux_err,                                    
                                   beam=0.0, pen=fit.tra_colors[j])   
     
            p3.addItem(err_)            
            
           # m = batman.TransitModel(fit.tr_params, t)    #initializes model
 
            #flux_model = m.light_curve(fit.tr_params)          #calculates light curve           
            #p3.plot(t, flux_model,pen=fit.tra_colors[-],symbol=None )   
            
 
            model_curve = p3.plot(t,y_model,  pen={'color':  fit.tra_colors[-1], 'width': self.tra_model_width.value()+1},
            enableAutoRange=True,viewRect=True ) 
            
            model_curve.setZValue(self.tra_model_z.value())            
            
            if self.trans_plot_cross_hair.isChecked():
                self.cross_hair(p3,log=False)     
            
            
            p4.plot(t, tr_o_c,        
            pen=None,  
            symbol=fit.pyqt_symbols_tra[i],
            symbolPen={'color': fit.tra_colors[j], 'width': 1.1},
            symbolSize=fit.pyqt_symbols_size_tra[i],enableAutoRange=True,viewRect=True,
            symbolBrush=fit.tra_colors[j] )             

            err_ = pg.ErrorBarItem(x=t, y=flux-flux_model, symbol='o', 
           # height=flux_err,
            top=flux_err,
            bottom=flux_err,            
            beam=0.0, pen=fit.tra_colors[j])               
            p4.addItem(err_)   
            
            #model_curve_o_c = p4.plot(t,y_model_o_c,  pen={'color':  fit.tra_colors[-1], 'width': self.tra_model_width.value()+1}, enableAutoRange=True,viewRect=True ) 
            
            #model_curve_o_c.setZValue(self.tra_model_z.value())                
            
  
            if self.trans_o_c_plot_cross_hair.isChecked():
                self.cross_hair(p4,log=False)  
          
            #model_curve = p4.plot(t, flux_model, pen={'color':  fit.tra_colors[-1], 'width': self.tra_model_width.value()+1},
            #enableAutoRange=True,viewRect=True )               
 
           # model_curve.setZValue(self.tra_model_z.value())   

        #else:    
        #    t = np.linspace(-0.25, 0.25, 1000)  #times at which to calculate light curve   
        #    m = batman.TransitModel(fit.tr_params, t)    #initializes model
        #    flux_model = m.light_curve(fit.tr_params)          #calculates light curve
 
        #    p3.plot(t, flux_model,pen='k',symbol=None )     
  

 
        
############################# N-Body ########################################     


    def update_orb_plot(self):
        global fit, p16
        
        p16.plot(clear=True,)    

        if fit.pl_arb_test == True:
            npl = fit.npl_arb
        else:
            npl = fit.npl     
        
        for i in range(npl):
            orb_xyz, pl_xyz, peri_xyz, apo_xyz = rv.planet_orbit_xyz(fit,i)        
            p16.plot(orb_xyz[0],orb_xyz[1], pen={'color': 0.5, 'width': 1.1},enableAutoRange=True,viewRect=True)   
            p16.plot((0,peri_xyz[0]),(0,peri_xyz[1]), pen={'color': 0.5, 'width': 1.1},enableAutoRange=True,viewRect=True)               
            
            p16.plot((pl_xyz[0],pl_xyz[0]), (pl_xyz[1],pl_xyz[1] ), pen=None,symbol='o', symbolSize=6,enableAutoRange=True,viewRect=True, symbolBrush='b') 
            
        p16.plot(np.array([0,0]), np.array([0,0]), pen=None,symbol='o', symbolSize=8,enableAutoRange=True,viewRect=True, symbolBrush='r')                


        
    def delta_omega_combo(self):

        if fit.pl_arb_test == True:
            npl = fit.npl_arb
        else:
            npl = fit.npl           
 
        self.comboBox_pl_1.clear()
        self.comboBox_pl_2.clear()
        
        for i in range(npl):
            self.comboBox_pl_1.addItem('omega %s'%str(i+1),i+1) 
            self.comboBox_pl_2.addItem('omega %s'%str(i+1),i+1) 
            
        self.comboBox_pl_1.setCurrentIndex(0)
        self.comboBox_pl_2.setCurrentIndex(1)


    def plot_delta_omega(self):
        global fit, colors_delta_om, p17 
        

        self.color_delta_om.setStyleSheet("color: %s;"%colors_delta_om[0]) 


        pl1_ind = self.comboBox_pl_1.currentIndex()
        pl2_ind = self.comboBox_pl_2.currentIndex()
 
        if pl1_ind ==-1 or pl2_ind ==-1 or len(fit.evol_p[pl1_ind]) ==0 or len(fit.evol_p[pl2_ind]) ==0:
            return
        else:
            last_stable = min(len(fit.evol_p[pl1_ind]),len(fit.evol_p[pl2_ind]))
        
        dom = (fit.evol_p[pl1_ind][0:last_stable] - fit.evol_p[pl2_ind][0:last_stable])%360

        if self.radioButton_dom_180_fold.isChecked():
            dom[dom>=180.0] -= 360.0
        
        
        
        p17.plot(clear=True,)
        p17.plot(fit.evol_T[0][0:last_stable], dom ,pen=None, #{'color': colors[i], 'width': 1.1},
        symbol='o',
        symbolPen={'color': colors_delta_om[0], 'width': 1.1},
        symbolSize=1,enableAutoRange=True,viewRect=True,
        symbolBrush=fit.colors[0]
        )  
 

        
    def get_delta_omega_color(self):
        global fit, colors_delta_om
        
        #colorz = QtGui.QColorDialog.getColor()
        colorz = self.colorDialog.getColor()
        colors_delta_om[0]=colorz.name()   
 
        self.plot_delta_omega()

    def delta_omega_plot_x_labels(self):
        global fit, p17
        
        text, okPressed = QtGui.QInputDialog.getText(self, "x-axis label","(No special characters!)", QtGui.QLineEdit.Normal, "")
        
        if okPressed and text != '':
            p17.setLabel('bottom', '%s'%text, units='',  **{'font-size':'9pt'})
 
        else:
            return
    
        self.plot_delta_omega()
 

    def delta_omega_plot_y_labels(self):
        global fit, p17
        
        text, okPressed = QtGui.QInputDialog.getText(self, "y-axis label","(No special characters!)", QtGui.QLineEdit.Normal, "")
        
        if okPressed and text != '':
            p17.setLabel('left', '%s'%text, units='',  **{'font-size':'9pt'})
 
        else:
            return
    
        self.plot_delta_omega()        
        
    ######### theta ############################## 
 

    def theta_which_pl_combo(self):
        
        if fit.pl_arb_test == True:
            npl = fit.npl_arb
        else:
            npl = fit.npl    

        self.comboBox_MMR_which_pl_1.clear()
        self.comboBox_MMR_which_pl_2.clear()
        
        for i in range(npl):
            self.comboBox_MMR_which_pl_1.addItem('pl. %s'%str(i+1),i+1) 
            self.comboBox_MMR_which_pl_2.addItem('pl. %s'%str(i+1),i+1) 
            
        self.comboBox_MMR_which_pl_1.setCurrentIndex(0)
        self.comboBox_MMR_which_pl_2.setCurrentIndex(1)

  
    def theta_combo(self):
 
        self.comboBox_MMR_theta.clear()
        
        ind1 = self.comboBox_MMR_pl_1.currentIndex()
        ind2 = self.comboBox_MMR_pl_2.currentIndex()              
        
        if ind1 == ind2:
            return
        
        for i in range(abs(ind2-ind1)+1):
            self.comboBox_MMR_theta.addItem('%s'%str(i+1),i+1) 
            
        self.comboBox_MMR_theta.setCurrentIndex(0)
  
        
    def MMR_combo(self):
 
        self.comboBox_MMR_pl_1.clear()
        self.comboBox_MMR_pl_2.clear()
        
        for i in range(10):
            self.comboBox_MMR_pl_1.addItem('%s'%str(i+1),i+1) 
            self.comboBox_MMR_pl_2.addItem('%s'%str(i+1),i+1) 
            
        self.comboBox_MMR_pl_1.setCurrentIndex(1)
        self.comboBox_MMR_pl_2.setCurrentIndex(0)        
        
        self.theta_combo()


    def plot_theta(self):
        global fit, colors_delta_om, p18

        self.color_theta.setStyleSheet("color: %s;"%colors_theta[0]) 

        pl1_ind = self.comboBox_MMR_which_pl_1.currentIndex()
        pl2_ind = self.comboBox_MMR_which_pl_2.currentIndex()
        
        Per_1 = self.comboBox_MMR_pl_1.currentIndex()
        Per_2 = self.comboBox_MMR_pl_2.currentIndex()  
        
        tet_n = self.comboBox_MMR_theta.currentIndex()
        
        if Per_1 == Per_2:
            return
 
        if pl1_ind ==-1 or pl2_ind ==-1 or len(fit.evol_p[pl1_ind]) ==0 or len(fit.evol_p[pl2_ind]) ==0:
            return
        else:
            last_stable = min(len(fit.evol_p[pl1_ind]),len(fit.evol_p[pl2_ind]))
      
        lambda1  = (fit.evol_M[pl1_ind][0:last_stable]    + fit.evol_p[pl1_ind][0:last_stable]   + 0)%360
        lambda2  = (fit.evol_M[pl2_ind][0:last_stable]    + fit.evol_p[pl2_ind][0:last_stable]   + 0)%360
 
        theta = {k: [ ] for k in range(10)}    
        coef1 = Per_2 +1
        coef2 = Per_1 +1
        order = abs(coef2 - coef1)
 
        for i in range(order+1):
 
            theta[i] = (coef1*lambda1%360 - coef2*lambda2%360 )%360 + (((coef2 -coef1) -i)*fit.evol_p[pl1_ind][0:last_stable] + i*fit.evol_p[pl2_ind][0:last_stable])%360 
            theta[i] = theta[i]%360

            if self.radioButton_theta_180_fold.isChecked():
                theta[i][theta[i]>=180.0] -= 360.0
 
        
        p18.plot(clear=True,)
        p18.plot(fit.evol_T[0][0:last_stable], theta[tet_n] ,pen=None, #{'color': colors[i], 'width': 1.1},
        symbol='o',
        symbolPen={'color': colors_theta[0], 'width': 1.1},
        symbolSize=1,enableAutoRange=True,viewRect=True,
        symbolBrush=fit.colors[0]
        )  
               
    def get_theta_color(self):
        global fit, colors_delta_om
        
        colorz = self.colorDialog.getColor()
        colors_theta[0]=colorz.name()   
 
        self.plot_theta()

    def theta_plot_x_labels(self):
        global fit, p18
        
        text, okPressed = QtGui.QInputDialog.getText(self, "x-axis label","", QtGui.QLineEdit.Normal, "")
        
        if okPressed and text != '':
            p18.setLabel('bottom', '%s'%text, units='',  **{'font-size':'9pt'})
 
        else:
            return
    
        self.plot_theta()
 

    def theta_plot_y_labels(self):
        global fit, p18
        
        text, okPressed = QtGui.QInputDialog.getText(self, "y-axis label","", QtGui.QLineEdit.Normal, "")
        
        if okPressed and text != '':
            p18.setLabel('left', '%s'%text, units='',  **{'font-size':'9pt'})
 
        else:
            return
    
        self.plot_theta()             
        
 
    def plot_i_Om(self):
        global fit, colors, p19

        if fit.pl_arb_test == True:
            npl = fit.npl_arb
        else:
            npl = fit.npl  
            
        p19.plot(clear=True,)   

        if self.plot_i.isChecked():
            for i in range(npl):
                p19.plot(fit.evol_T[i], fit.evol_i[i] ,pen=fit.colors[i],symbol=None )    
            p19.setLabel('left', 'i [deg]', units='',  **{'font-size':'9pt'})    
    
        elif self.plot_Om.isChecked():
            for i in range(npl):
                
                Om_evol = np.array(fit.evol_Om[i])
                
                if self.radioButton_Omega_180_fold.isChecked():
                    Om_evol[Om_evol>=180.0] -= 360.0

                p19.plot(fit.evol_T[i], Om_evol ,pen=None, #{'color': colors[i], 'width': 1.1},
                symbol='o',
                symbolPen={'color': fit.colors[i], 'width': 1.1},
                symbolSize=1,enableAutoRange=True,viewRect=True,
                symbolBrush=fit.colors[i]
                )
                
                
            p19.setLabel('left', 'Omega [deg]', units='',  **{'font-size':'9pt'})    
            Om_evol = 0
 
    def plot_energy(self):
        global fit, colors, p20

        p20.plot(clear=True,)  
        
        if len(np.atleast_1d(fit.evol_T_energy)) < 3:
            return
        

        if self.radioButton_energy.isChecked():
 
            p20.plot(fit.evol_T_energy, fit.evol_energy ,pen=fit.colors[0],symbol=None )    
            p20.setLabel('left', 'Energy', units='',  **{'font-size':'9pt'})    
    
        elif self.radioButton_lx.isChecked():
            p20.plot(fit.evol_T_energy, fit.evol_momentum['lx'] ,pen=fit.colors[0],symbol=None )    
            p20.setLabel('left', 'Momentum lx', units='',  **{'font-size':'9pt'})    

        elif self.radioButton_ly.isChecked():
            p20.plot(fit.evol_T_energy, fit.evol_momentum['ly'] ,pen=fit.colors[0],symbol=None )    
            p20.setLabel('left', 'Momentum ly', units='',  **{'font-size':'9pt'}) 
            
        elif self.radioButton_lz.isChecked():
            p20.plot(fit.evol_T_energy, fit.evol_momentum['lz'] ,pen=fit.colors[0],symbol=None )    
            p20.setLabel('left', 'Momentum lz', units='',  **{'font-size':'9pt'})                    
 
    
    def plot_evol_a(self):
        global fit, colors, p13
        
        if fit.pl_arb_test == True:
            npl = fit.npl_arb
        else:
            npl = fit.npl   
            
        if not self.hold_old_plot_a.isChecked():    
            p13.plot(clear=True,)
        else:
            p13.plot(clear=False,)
            
 
        for i in range(npl):
            p13.plot(fit.evol_T[i], fit.evol_a[i] ,pen=fit.colors[i],symbol=None )     
 
            
 
    def plot_evol_e(self):
        global fit, colors,   p14 
        
        if fit.pl_arb_test == True:
            npl = fit.npl_arb
        else:
            npl = fit.npl      
    
 
        if not self.hold_old_plot_e.isChecked():    
            p14.plot(clear=True,)
        else:
            p14.plot(clear=False,)
 

        for i in range(npl):
            p14.plot(fit.evol_T[i], fit.evol_e[i] ,pen=fit.colors[i],symbol=None )  
 
            
    def plot_evol_p(self):
        global fit, colors, p15  
        
        if fit.pl_arb_test == True:
            npl = fit.npl_arb
        else:
            npl = fit.npl      
    
 
        if not self.hold_old_plot_om.isChecked():    
            p15.plot(clear=True,)
        else:
            p15.plot(clear=False,)
 
        for i in range(npl):
 
            p15.plot(fit.evol_T[i], fit.evol_p[i] ,pen=None, #{'color': colors[i], 'width': 1.1},
        symbol='o',
        symbolPen={'color': fit.colors[i], 'width': 1.1},
        symbolSize=1,enableAutoRange=True,viewRect=True,
        symbolBrush=fit.colors[i]
        )    
            
      
            
            
    def plot_evol_all(self):
        global fit        
        
        self.plot_evol_a()
        self.plot_evol_e()
        self.plot_evol_p()
        

        self.delta_omega_combo()
        self.theta_which_pl_combo()
        self.plot_delta_omega()
        self.plot_theta()
        self.plot_i_Om()    
        self.plot_energy()              
        
    def worker_Nbody_complete(self):
        global fit, colors, p13, p14, p15  
          


        self.plot_evol_all()
#        self.update_orb_plot()

        self.plot_tabs.setCurrentWidget(self.tab_Orbital_evol)
            
             
        self.button_orb_evol.setEnabled(True)       
        self.statusBar().showMessage('')      
          
 
    def worker_Nbody(self):
        global fit  

        self.button_orb_evol.setEnabled(False)         
        npl_to_fit = np.atleast_1d(fit.fit_results.mass)
 
        # check if any fits where performed, and tus planets present
        if len(np.atleast_1d(npl_to_fit)) == 1 and npl_to_fit[0] <= 0:
            choice = QtGui.QMessageBox.information(self, 'Warning!',
            "Not possible to integrate a fit that does not exist. First perform an orbital fitting and then test the orbital stability. Okay?", QtGui.QMessageBox.Ok)      
            self.button_orb_evol.setEnabled(True)         
            return        

        if fit.npl < 2:
            choice = QtGui.QMessageBox.information(self, 'Warning!'," With less than two planets this makes no sense. Okay?",
                                            QtGui.QMessageBox.Ok) 
            self.button_orb_evol.setEnabled(True)                    
            return
        
        ######## this is a fix in case one adds another planet without initialize it 
        if fit.npl != len(fit.fit_results.mass):  
            fit.model_saved = False 
            self.init_fit()

        fit.pl_arb_test = False        
 
        self.statusBar().showMessage('Running Orbital Evolution......')   


        # Pass the function to execute
        worker3 = Worker(lambda: self.run_orbital_simulations()) # Any other args, kwargs are passed to the run  
        # Execute
        worker3.signals.finished.connect(self.worker_Nbody_complete)
        
        # worker.signals.result.connect(self.print_output)
        #worker.signals.finished.connect(self.thread_complete)
       # worker.signals.progress.connect(self.progress_fn)
        self.threadpool.start(worker3)  



    def worker_Nbody_arb(self):
        global fit  

        self.run_orb_evol_arbitary.setEnabled(False)           
        self.check_arb_pl()

        print(fit.npl_arb)

        if fit.npl_arb < 2:
            choice = QtGui.QMessageBox.information(self, 'Warning!'," With less than two planets this makes no sense. Okay?",
                                            QtGui.QMessageBox.Ok) 
            self.run_orb_evol_arbitary.setEnabled(True)                    
            return
 
        fit.pl_arb_test = True
        
        self.statusBar().showMessage('Running Orbital Evolution......')   
        
        # Pass the function to execute
        worker_arb = Worker(lambda: self.run_orbital_simulations(arbitary=True)) # Any other args, kwargs are passed to the run  
        # Execute
        worker_arb.signals.finished.connect(self.worker_Nbody_complete)
        
        # worker.signals.result.connect(self.print_output)
        #worker.signals.finished.connect(self.thread_complete)
       # worker.signals.progress.connect(self.progress_fn)
        self.threadpool.start(worker_arb)  

        self.run_orb_evol_arbitary.setEnabled(True)   
        
        
    def run_orbital_simulations(self, arbitary=False):
        global fit

        #self.max_time_of_evol

        if self.radioButton_SyMBA.isChecked():
            integrator = 'symba'
        elif self.radioButton_MVS.isChecked():
            integrator = 'mvs'        
        elif self.radioButton_MVS_GR.isChecked():       
             integrator = 'mvs_gr'
             fit.GR_step = self.mvs_GR_steps.value()
             
        import time
        start_time = time.time()        
       # fit.run_stability_last_fit_params(timemax=self.max_time_of_evol.value(), timestep=self.time_step_of_evol.value(), integrator=integrator)      
       
        if arbitary == True:
            fit = rv.run_stability_arb(fit, timemax=self.max_time_of_evol.value(), timestep=self.time_step_of_evol.value(), integrator=integrator)      
        else:         
            fit = rv.run_stability(fit, timemax=self.max_time_of_evol.value(), timestep=self.time_step_of_evol.value(), integrator=integrator)      
        
        print("--- %s seconds ---" % (time.time() - start_time))          


         
       
############################# Fortran fitting ###############################        
        
    def worker_RV_fitting_complete(self):
        global fit  
        
        fit=rv.get_xyz(fit)
        
        self.update_labels()
        self.update_gui_params()
        self.update_errors() 
        self.update_a_mass()    
                
                 
        self.statusBar().showMessage('')   
        #self.console_widget.print_text(str(fit.print_info(short_errors=False))) 

        self.jupiter_push_vars()   
        self.button_fit.setEnabled(True)  
        self.run_gls()
        self.run_gls_o_c()        
        self.update_plots()  
 
    def worker_RV_fitting(self, ff=20, m_ln=True, auto_fit = False , init = False ):
        global fit  
        
        self.button_fit.setEnabled(False)         
        
        # check if RV data is present
        if fit.filelist.ndset <= 0:
             choice = QtGui.QMessageBox.information(self, 'Warning!',
             "Not possible to look for planets if there are no RV data loaded. Please add your RV data first. Okay?", QtGui.QMessageBox.Ok)      
             self.button_fit.setEnabled(True)         
             return   
         
       # self.set_tra_ld()   

        self.check_model_params()
         
        self.check_bounds()
        self.check_priors_nr()   
        self.check_priors_jeff()   
        
        
        fit.model_npoints = self.points_to_draw_model.value()
        #self.tabWidget_helper.setCurrentWidget(self.tab_info)

        if init == True:
            fit.init_fit= True
            ff = 0
            doGP=False
        else:
            doGP=self.do_RV_GP.isChecked()            
            fit.init_fit= False  
       
        
        if self.radioButton_fortran77.isChecked() and not self.do_RV_GP.isChecked() or init == True:
            self.statusBar().showMessage('Minimizing parameters....')    
             # Pass the function to execute
            worker2 = Worker(lambda:  self.optimize_fit(ff=ff, doGP=doGP, minimize_fortran=True, m_ln=m_ln, auto_fit = auto_fit)) # Any other args, kwargs are passed to the run  
 
        else:    
              
            self.check_scipy_min()

            self.statusBar().showMessage('Minimizing parameters using SciPyOp (might be slow)....')                 
            worker2 = Worker(lambda:  self.optimize_fit(ff=0, doGP=self.do_RV_GP.isChecked(),  gp_kernel_id=-1, minimize_fortran=False, m_ln=m_ln, auto_fit = auto_fit)) # Any other args, kwargs are passed to the run  
            self.tabWidget_helper.setCurrentWidget(self.tab_info)
            
            
        worker2.signals.finished.connect(self.worker_RV_fitting_complete)
        
        # worker.signals.result.connect(self.print_output)
        #worker.signals.finished.connect(self.thread_complete)
       # worker.signals.progress.connect(self.progress_fn)
        self.threadpool.start(worker2) 
        

     
    def update_RV_jitter_flag(self):
        global fit

        use_data_jitter_gui = [self.use_jitter_Data1,self.use_jitter_Data2,self.use_jitter_Data3,self.use_jitter_Data4,self.use_jitter_Data5,
                               self.use_jitter_Data6,self.use_jitter_Data7,self.use_jitter_Data8,self.use_jitter_Data9,self.use_jitter_Data10]

        if self.amoeba_radio_button.isChecked():
            for i in range(10):     
                use_data_jitter_gui[i].setEnabled(True)
        else:
            for i in range(10):     
                use_data_jitter_gui[i].setEnabled(False)
                
  
    def update_dyn_kep_flag(self):
        global fit

        if self.radioButton_Dynamical.isChecked():
            fit.mod_dynamical = True
        else:
            fit.mod_dynamical = False
            
          
    def optimize_fit(self,ff=20,m_ln=True, doGP=False, gp_kernel_id=-1, auto_fit = False, minimize_fortran=True):  
        global fit
        
        if not auto_fit:
            self.update_params()
 
            
        if self.radioButton_Dynamical.isChecked():
            #fit.mod_dynamical = True
            f_kill = self.dyn_model_to_kill.value()
            if ff > 1:
                ff = 1 #TBF
           # ff = 1
        else:
            #fit.mod_dynamical = False
            f_kill = self.kep_model_to_kill.value()    
        
        if minimize_fortran==False:
            ff = 0 
            
      #  print(ff)   

        if m_ln == True and doGP == False:
           # if ff > 0:        
           #     """
           #     run one time using the L-M method ignorring the jitter (for speed)
          #      """
               # fit.fitting(fileinput=False,outputfiles=[1,0,0], doGP=doGP, kernel_id=gp_kernel_id, minimize_fortran=minimize_fortran, fortran_kill=f_kill, timeout_sec=self.master_timeout.value(),minimize_loglik=False,amoeba_starts=ff, print_stat=False, eps=self.dyn_model_accuracy.value(), dt=self.time_step_model.value())
            """
            now run the amoeba code modeling the jitters
            """
#            fit.fitting(fileinput=False,outputfiles=[1,0,0], doGP=doGP,  kernel_id=gp_kernel_id,  minimize_fortran=minimize_fortran,  fortran_kill=f_kill, timeout_sec=self.master_timeout.value(),minimize_loglik=True,amoeba_starts=ff, print_stat=False, eps=self.dyn_model_accuracy.value(), dt=self.time_step_model.value())
            fit.fitting(fileinput=False,outputfiles=[1,1,1], doGP=doGP,  kernel_id=gp_kernel_id,  minimize_fortran=minimize_fortran,  fortran_kill=f_kill, timeout_sec=self.master_timeout.value(),minimize_loglik=True,amoeba_starts=ff, print_stat=False, eps=self.dyn_model_accuracy.value(), dt=self.time_step_model.value(), npoints=self.points_to_draw_model.value(), model_max= self.model_max_range.value(), model_min= self.model_min_range.value())

        elif m_ln == True and doGP == True:       
            fit.fitting(fileinput=False,outputfiles=[1,1,1], doGP=doGP,  kernel_id=gp_kernel_id,  minimize_fortran=minimize_fortran,  fortran_kill=f_kill, timeout_sec=self.master_timeout.value(),minimize_loglik=True,amoeba_starts=ff,  print_stat=False, eps=self.dyn_model_accuracy.value(), dt=self.time_step_model.value(), npoints=self.points_to_draw_model.value(), model_max= self.model_max_range.value(), model_min= self.model_min_range.value())
        
        elif m_ln == False and   minimize_fortran==False:      
            fit.fitting(fileinput=False,outputfiles=[1,1,1], doGP=doGP,  kernel_id=gp_kernel_id,  minimize_fortran=minimize_fortran, fortran_kill=f_kill, timeout_sec=self.master_timeout.value(),minimize_loglik=True,amoeba_starts=0, print_stat=False,eps=self.dyn_model_accuracy.value(), dt=self.time_step_model.value(), npoints=self.points_to_draw_model.value(), model_max= self.model_max_range.value(), model_min= self.model_min_range.value())
        
        else:      
            fit.fitting(fileinput=False,outputfiles=[1,1,1], doGP=doGP,  kernel_id=gp_kernel_id,  minimize_fortran=minimize_fortran, fortran_kill=f_kill, timeout_sec=self.master_timeout.value(),minimize_loglik=m_ln,amoeba_starts=ff, print_stat=False,eps=self.dyn_model_accuracy.value(), dt=self.time_step_model.value(), npoints=self.points_to_draw_model.value(), model_max= self.model_max_range.value(), model_min= self.model_min_range.value())

        if fit.doGP == True:
            rv.get_gps_model(fit)

        for i in range(fit.npl):
             rv.phase_RV_planet_signal(fit,i+1)  

        if auto_fit:
                                          
            self.update_labels()
            self.update_gui_params()
            self.update_errors() 
            self.update_a_mass() 

            #self.run_gls()
            self.run_gls_o_c()                   
            self.update_plots()                   
            self.statusBar().showMessage('')           
            self.jupiter_push_vars()




       
    def print_info_for_object(self,text):
        #self.dialog.statusBar().showMessage('Ready')
        self.dialog.setGeometry(300, 300, 450, 250)
        self.dialog.setWindowTitle('Detailed Info')  
 
        self.dialog.text.setPlainText(text)
        self.dialog.text.setReadOnly(True)       
        #self.dialog.setWindowIcon (QtGui.QIcon('logo.png'))        
        self.dialog.show()




    def print_chi_table(self):
        #self.dialog.statusBar().showMessage('Ready')
        self.dialog_chi_table.setFixedSize(900,350)
        self.dialog_chi_table.setWindowTitle('Confidence intervals table')  
        #self.dialog.setGeometry(300, 300, 800, 800)
        #self.dialog_credits.acceptRichText(True)
        
        # standard deviations to calculate 
        # (You can ofcourse simplify here for only 1,2 and 3 sigma)
        sigma = [1.0,  np.sqrt(stat.chi2.ppf(0.8,1)),
            np.sqrt(stat.chi2.ppf(0.9,1)),
            np.sqrt(stat.chi2.ppf(0.95,1)),          2.0,
            np.sqrt(stat.chi2.ppf(0.99,1)),        3.0,
            np.sqrt(stat.chi2.ppf(0.999,1)),       4.0,5.0  ]
        
        #confidence intervals these sigmas represent:
        conf_int = [ stat.chi2.cdf( s**2,1) for s in sigma ]
        
        #degrees of freedom to calculate instead of 10 you can put k = 20, 30 or even 100
        dof = range(1,12)        
        
        text = ''
        self.dialog_chi_table.text.setText(text) 
        
        text = "sigma     \t" + "\t".join(["%1.2f"%(s) for s in sigma])
        self.dialog_chi_table.text.append(text)        
        
        text = "conf_int  \t" + "\t".join(["%1.2f%%"%(100*ci) for ci in conf_int])
        self.dialog_chi_table.text.append(text)        

        text = "p-value   \t" + "\t".join(["%1.5f"%(1-ci) for ci in conf_int])        
        self.dialog_chi_table.text.append(text)        
 
    
        for d in dof:
            chi_squared = [stat.chi2.ppf( ci, d) for ci in conf_int ]
            text = "chi2(k=%d)\t"%d + "\t".join(["%1.2f" % c for c in chi_squared])     
            self.dialog_chi_table.text.append(text)  


        self.dialog_chi_table.text.setReadOnly(True)       
   
        self.dialog_chi_table.show()



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
        
        text = "* " + "<a href='https://github.com/hippke/tls'>Transit Least Squares</a>" 
        self.dialog_credits.text.append(text)             

        text = "* " + "<a href='https://www.boulder.swri.edu/~hal/swift.html'>swift</a>" 
        self.dialog_credits.text.append(text)        
                        
        text = "* " + "<a href='https://github.com/jupyter/qtconsole'>qtconsole</a>"
        self.dialog_credits.text.append(text)        

        text = "* " + "<a href='https://github.com/mfitzp/15-minute-apps/tree/master/wordprocessor'>megasolid idiom</a>" 
        self.dialog_credits.text.append(text)  
        
        text = "* " + "<a href='https://dynesty.readthedocs.io/en/latest/'>dynesty</a>" 
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
        
        self.dialog_credits.setStyleSheet(" QTextEdit{border-image: url(./lib/33_striker.png) 0 0 0 0 stretch stretch;} ")

        #self.dialog.setWindowIcon (QtGui.QIcon('logo.png'))        
        
        self.dialog_credits.show()
        
  

    def run_nest_samp(self):
        global fit
        choice = QtGui.QMessageBox.information(self, 'Warning!', "Not available yet. Okay?", QtGui.QMessageBox.Ok) 


    def find_planets(self):
        global fit 

        # check if RV data is present
        if fit.filelist.ndset <= 0:
             choice = QtGui.QMessageBox.information(self, 'Warning!',
             "Not possible to look for planets if there are no RV data loaded. Please add your RV data first. Okay?", QtGui.QMessageBox.Ok)      
             self.button_auto_fit.setEnabled(True)         
             return        

        # the first one on the data GLS
        if fit.gls.power.max() <= fit.gls.powerLevel(self.auto_fit_FAP_level.value()):
             choice = QtGui.QMessageBox.information(self, 'Warning!',
             "No significant power on the GLS. Therefore no planets to fit OK?", QtGui.QMessageBox.Ok)      
             self.button_auto_fit.setEnabled(True)                                                           
             return
        
        else:
            if fit.npl !=0:
                for j in range(fit.npl):
                    fit.remove_planet(fit.npl-(j+1))

            mean_anomaly_from_gls = np.degrees((((fit.epoch - float(fit.gls.hpstat["T0"]) )% (fit.gls.hpstat["P"]) )/ (fit.gls.hpstat["P"]) ) * 2*np.pi)
             
            fit.add_planet(fit.gls.hpstat["amp"],fit.gls.hpstat["P"],0.0,0.0,mean_anomaly_from_gls -90.0,90.0,0.0)
            fit.use.update_use_planet_params_one_planet(0,True,True,True,True,True,False,False)     
            self.update_use_from_input_file()   
            self.update_use()                     
            self.optimize_fit(20,m_ln=self.amoeba_radio_button.isChecked(),auto_fit = True)
            
            #now inspect the residuals
            
            for i in range(1,int(self.auto_fit_N_planets.value())):
                
                if fit.gls_o_c.power.max() <= fit.gls_o_c.powerLevel(self.auto_fit_FAP_level.value()):
                    for j in range(fit.npl):
                        fit.use.update_use_planet_params_one_planet(j,True,True,True,True,True,False,False)     
            
                    self.update_use_from_input_file()   
                    self.update_use()                     
                    self.optimize_fit(20,m_ln=self.amoeba_radio_button.isChecked(),auto_fit = True) 
                    self.button_auto_fit.setEnabled(True)     
                    return
                #elif (1/RV_per_res.hpstat["fbest"]) > 1.5:
                else:    
                    mean_anomaly_from_gls = np.degrees((((fit.epoch - float(fit.gls_o_c.hpstat["T0"]) )% (fit.gls_o_c.hpstat["P"]) )/ (fit.gls_o_c.hpstat["P"]) ) * 2*np.pi)
             
                    fit.add_planet(fit.gls_o_c.hpstat["amp"],fit.gls_o_c.hpstat["P"],0.0,0.0,mean_anomaly_from_gls -90.0,90.0,0.0)
                    fit.use.update_use_planet_params_one_planet(i,True,True,True,True,True,False,False)  

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
        self.button_auto_fit.setEnabled(False)         
        
        if fit.npl != 0:        
            choice = QtGui.QMessageBox.information(self, 'Warning!',
                                            "Planets already exist. Do you want to overwrite the analysis?",
                                            QtGui.QMessageBox.No | QtGui.QMessageBox.Yes)  
         
            if choice == QtGui.QMessageBox.No:
                self.button_auto_fit.setEnabled(True)         
                return
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


        self.K4.minimize_signal.connect(lambda: fit.minimize_one_param_K(3)) #TBD!
        self.K4.minimize_signal.connect(self.init_fit) #TBD!       
        self.P4.minimize_signal.connect(lambda: fit.minimize_one_param_P(3)) #TBD!
        self.P4.minimize_signal.connect(self.init_fit) #TBD!     
        self.e4.minimize_signal.connect(lambda: fit.minimize_one_param_e(3)) #TBD!
        self.e4.minimize_signal.connect(self.init_fit) #TBD!  
        self.om4.minimize_signal.connect(lambda: fit.minimize_one_param_w(3)) #TBD!
        self.om4.minimize_signal.connect(self.init_fit) #TBD!  
        self.ma4.minimize_signal.connect(lambda: fit.minimize_one_param_M0(3)) #TBD!
        self.ma4.minimize_signal.connect(self.init_fit) #TBD!  
        
        self.K5.minimize_signal.connect(lambda: fit.minimize_one_param_K(4)) #TBD!
        self.K5.minimize_signal.connect(self.init_fit) #TBD!       
        self.P5.minimize_signal.connect(lambda: fit.minimize_one_param_P(4)) #TBD!
        self.P5.minimize_signal.connect(self.init_fit) #TBD!     
        self.e5.minimize_signal.connect(lambda: fit.minimize_one_param_e(4)) #TBD!
        self.e5.minimize_signal.connect(self.init_fit) #TBD!  
        self.om5.minimize_signal.connect(lambda: fit.minimize_one_param_w(4)) #TBD!
        self.om5.minimize_signal.connect(self.init_fit) #TBD!  
        self.ma5.minimize_signal.connect(lambda: fit.minimize_one_param_M0(4)) #TBD!
        self.ma5.minimize_signal.connect(self.init_fit) #TBD!         
        
        self.K6.minimize_signal.connect(lambda: fit.minimize_one_param_K(5)) #TBD!
        self.K6.minimize_signal.connect(self.init_fit) #TBD!       
        self.P6.minimize_signal.connect(lambda: fit.minimize_one_param_P(5)) #TBD!
        self.P6.minimize_signal.connect(self.init_fit) #TBD!     
        self.e6.minimize_signal.connect(lambda: fit.minimize_one_param_e(5)) #TBD!
        self.e6.minimize_signal.connect(self.init_fit) #TBD!  
        self.om6.minimize_signal.connect(lambda: fit.minimize_one_param_w(5)) #TBD!
        self.om6.minimize_signal.connect(self.init_fit) #TBD!  
        self.ma6.minimize_signal.connect(lambda: fit.minimize_one_param_M0(5)) #TBD!
        self.ma6.minimize_signal.connect(self.init_fit) #TBD!   

 
        
    def jupiter_push_vars(self):
        global fit        
        self.console_widget.push_vars({'fit':fit})    
        #self.console_widget.push_vars({'pg':pg})    

        #self.console_widget.clear()         
        #self.console_widget.print_text(str("Welcome!"+"\n")) 


########################## Sessions ##################################
 
    def getNewses(self):
        global fit, ses_list  
        
        text, okPressed = QtGui.QInputDialog.getText(self, "New session","Name session: (No special characters!)", QtGui.QLineEdit.Normal, "")
        if okPressed and text != '':
            
            if len(ses_list) == 0:
                ses_list.append(fit)
                      
            #file_pi = open('.sessions/empty.ses', 'rb')
            #fit_new = dill.load(file_pi)
            #file_pi.close()
            fit_new=rv.signal_fit(name=text)

            fit_new.name=text
            ses_list.append(fit_new)
 
            self.session_list()
            
    def rem_ses(self):
        global fit, ses_list  
        
        ind = self.comboBox_select_ses.currentIndex()
        
        choice = QtGui.QMessageBox.information(self, 'Warning!',
        "Do you really want to remove Session %s"%(ses_list[ind].name),
                                            QtGui.QMessageBox.No | QtGui.QMessageBox.Yes)  
         
        if choice == QtGui.QMessageBox.No:
            return
        elif choice == QtGui.QMessageBox.Yes: # and ind <=0:
            if len(ses_list)==1:
                ses_list.pop(0)         
                self.new_session()
                self.session_list() 
                self.select_session(0)
            else:
                ses_list.pop(ind)
                self.session_list() 
                self.select_session(ind-1)
            
    def cop_ses(self):
        global fit, ses_list  
        
        ind = self.comboBox_select_ses.currentIndex()
        
        if len(ses_list) == 0:
            ses_list.append(fit)
        if ind < 0:
            fit_new =dill.copy(ses_list[0])
        elif ind >= 0:
            fit_new =dill.copy(ses_list[ind])
      
        ses_list.append(fit_new)
           
        self.session_list() 
        self.select_session(ind-1)            
  

    def new_session(self):
        global fit, ses_list
        
        #file_pi = open('.sessions/empty.ses', 'rb')
        #fit_new = dill.load(file_pi)
        #file_pi.close()     
        fit_new=rv.signal_fit(name='new Session')

        ses_list.append(fit_new)
        self.session_list()
              
            
    def open_session(self):
        global fit,ses_list
        
        input_file = QtGui.QFileDialog.getOpenFileName(self, 'Open session', '', 'Data (*.ses)')

        if str(input_file[0]) != '':

            file_pi = open(input_file[0], 'rb')
            fit_new = dill.load(file_pi)
            file_pi.close()     
            
            ses_list.append(fit_new)
            
            #self.check_settings()
            rv.check_temp_RV_file(fit_new)
            self.session_list()
            self.select_session(-1)
            
 
        

    def save_session(self):
        global fit
        
        output_file = QtGui.QFileDialog.getSaveFileName(self, 'Save session', '%s.ses'%fit.name, 'Data (*.ses)')
        
        if str(output_file[0]) != '':
            file_pi = open(output_file[0], 'wb')
            dill.dump(fit, file_pi)
            file_pi.close()


    def open_sessions(self):
        global fit, ses_list
        
        input_file = QtGui.QFileDialog.getOpenFileName(self, 'Open session', '', 'Data (*.mses)')

        if str(input_file[0]) != '':

            file_pi = open(input_file[0], 'rb')
            fit2 = dill.load(file_pi)
            file_pi.close()           
        
            choice = QtGui.QMessageBox.information(self, 'Warning!',
                                            "Do you want to overwrite the current sessions? If you choose 'No' will add the session, 'Cancel' will exit",
                                            QtGui.QMessageBox.Cancel | QtGui.QMessageBox.No | QtGui.QMessageBox.Yes)  
         
            if choice == QtGui.QMessageBox.No:
                ses_list = ses_list + fit2
            elif choice == QtGui.QMessageBox.Yes:
                ses_list = fit2
            elif choice == QtGui.QMessageBox.Cancel:        
                return         

            #fit_new.check_temp_RV_file()
    
            self.check_settings()  
            self.session_list()
            self.select_session(-1)

    def save_sessions(self):
        global fit, ses_list
      
        output_file = QtGui.QFileDialog.getSaveFileName(self, 'Save multi-session', 'my_sessions.mses', 'Data (*.mses)')

        if str(output_file[0]) != '':
            file_pi = open(output_file[0], 'wb')
            dill.dump(ses_list, file_pi)
            file_pi.close()

        
    def session_list(self):
        global fit, ses_list
              
        if len(ses_list) == 0:
            self.comboBox_select_ses.clear()
            self.comboBox_select_ses.addItem('1 %s'%fit.name) 

        elif len(ses_list) != 0:
            self.comboBox_select_ses.clear()
            for i in range(len(ses_list)):
                self.comboBox_select_ses.addItem('%s'%(ses_list[i].name),i) 
                self.comboBox_select_ses.setItemText(i, '%s'%(ses_list[i].name))
 

    def select_session(self, index):
        global fit, ses_list
        
        if index == -1:
            ind = -1
            self.comboBox_select_ses.setCurrentIndex(len(ses_list)-1)
        else:
            ind = self.comboBox_select_ses.itemData(index) 

        if ind == None:
            return
        else:
            fit = ses_list[ind]

        
        self.check_settings()

        self.init_fit()

        self.update_use_from_input_file()   
        self.update_use()
        self.update_gui_params()
        self.update_params()
        self.update_RV_file_buttons() 
        self.update_color_picker()
        
        if not ind == None:    
            ses_list[ind] = fit 
 
    def change_session_label(self):
        global fit, ses_list 
        
        ind = self.comboBox_select_ses.currentIndex()
        if ind == None:
            return
        elif ind >= 0:
            ses_list[ind].name = self.comboBox_select_ses.currentText() 
            self.comboBox_select_ses.setItemText(ind, '%s'%(ses_list[ind].name))
            

            
################################## Nest Samp. #######################################

    def worker_nest_complete(self):
        global fit  
        #fit.print_info(short_errors=False)
        
        self.update_labels()
        self.update_gui_params()
        self.update_errors() 
        self.update_a_mass() 
        
        self.statusBar().showMessage('') 
        #self.console_widget.print_text(str(fit.print_info(short_errors=False))) 
        
        if self.adopt_nest_means_as_par.isChecked() or self.adopt_nest_best_lnL_as_pars.isChecked() or self.adopt_nest_mode_as_par.isChecked():
            self.init_fit()
 
 

    def worker_nest(self):
        global fit  
        
        
        if self.radioButton_RV.isChecked():
            fit.rtg = [True,self.do_RV_GP.isChecked(), False, self.do_tra_GP.isChecked()]
        elif self.radioButton_transit.isChecked():
            fit.rtg = [False, self.do_RV_GP.isChecked(), True, self.do_tra_GP.isChecked()]
        elif self.radioButton_transit_RV.isChecked():
            fit.rtg = [True,self.do_RV_GP.isChecked(), True, self.do_tra_GP.isChecked()]
        
        self.button_nest_samp.setEnabled(False)
        self.statusBar().showMessage('Nested Sampling in progress....')        
        # check if RV data is present
        if fit.type_fit["RV"] == True and fit.filelist.ndset <= 0:
             choice = QtGui.QMessageBox.information(self, 'Warning!',
             "Not possible to run MCMC if there are no RV data loaded. Please add your RV data first. Okay?", QtGui.QMessageBox.Ok)      
             self.button_nest_samp.setEnabled(True)  
             self.statusBar().showMessage('') 

             return   
        
        ntran_data = 0
        for i in range(0,10,1):         
            ntran_data += len(fit.tra_data_sets[i]) 
            
        if fit.type_fit["Transit"] == True  and ntran_data == 0:
             choice = QtGui.QMessageBox.information(self, 'Warning!',
             "Not possible to run MCMC if there are no transit data loaded. Please add your transit data first. Okay?", QtGui.QMessageBox.Ok)      
             self.button_nest_samp.setEnabled(True)  
             self.statusBar().showMessage('') 

             return             
            
            

        choice = QtGui.QMessageBox.information(self, 'Warning!',
                                            "This will run in the background and may take some time. Results are printed in the 'Stdout/Stderr' tab. Okay?",
                                            QtGui.QMessageBox.Cancel | QtGui.QMessageBox.Ok)       
 
         
        if choice == QtGui.QMessageBox.Cancel:
            self.statusBar().showMessage('') 
            self.button_nest_samp.setEnabled(True)
            return        
        
        self.set_tra_ld()
        self.check_bounds()
        self.check_priors_nr() 
        fit.model_npoints = self.points_to_draw_model.value()
        
        self.tabWidget_helper.setCurrentWidget(self.tab_info)
        
        
        if self.use_nest_percentile_level.isChecked():
            fit.nest_percentile_level = self.nest_percentile_level.value()
        else:
            fit.nest_percentile_level = 68.3
           
        
        # Pass the function to execute
        worker_n = Worker(lambda: self.run_nest()) # Any other args, kwargs are passed to the run  
        # Execute
        worker_n.signals.finished.connect(self.worker_nest_complete)
        
        # worker.signals.result.connect(self.print_output)
        #worker.signals.finished.connect(self.thread_complete)
       # worker.signals.progress.connect(self.progress_fn)
        self.threadpool.start(worker_n)
        
    def run_nest(self):
        global fit
        
        self.check_model_params() 
        self.check_nested_params()
      
        fit = rv.run_nestsamp(fit, threads=int(self.nest_N_threads.value()), std_output=False, stop_crit = self.stop_crit.value(), 
        Dynamic_nest = self.radioButton_dyn_nest_samp.isChecked(), live_points = int(self.live_points.value()),fileoutput=self.save_samples.isChecked(),
        save_means=self.adopt_nest_means_as_par.isChecked(), save_mode=self.adopt_nest_mode_as_par.isChecked(), save_maxlnL=self.adopt_nest_best_lnL_as_pars.isChecked())
     
        self.button_nest_samp.setEnabled(True)            
 
    def change_nest_samples_file_name(self):
        global fit
        
        output_file = QtGui.QFileDialog.getSaveFileName(self, 'path and name of the nested samples', '', '')
        
        if output_file[0] != '':
            fit.nest_sample_file = output_file[0] 
            self.nest_samples_change_name.setText(output_file[0])
        else:
            return

    def check_nested_params(self):
        global fit
        #fit.gaussian_ball = self.init_gauss_ball.value() 
        fit.live_points_fact = int(self.live_points.value())
        


    def force_nest_check_box(self):
        if self.button_make_cornerplot_nested.isChecked():
            self.save_samples_nested.setChecked(True)
 



################################## MCMC #######################################

    def worker_mcmc_complete(self):
        global fit  
        #fit.print_info(short_errors=False)
        
        self.update_labels()
        self.update_gui_params()
        self.update_errors() 
        self.update_a_mass() 
        
        self.init_plot_corr()
        
        self.statusBar().showMessage('') 
        #self.console_widget.print_text(str(fit.print_info(short_errors=False))) 
        
        if self.adopt_mcmc_means_as_par.isChecked() or self.adopt_best_lnL_as_pars.isChecked() or self.adopt_mcmc_mode_as_par.isChecked():
            self.init_fit()
 
       # if sys.version_info[0] == 3:
       #     self.print_py3_warning()

    def worker_mcmc(self):
        global fit  
        
        
        if self.radioButton_RV.isChecked():
            fit.rtg = [True,self.do_RV_GP.isChecked(), False, self.do_tra_GP.isChecked()]
        elif self.radioButton_transit.isChecked():
            fit.rtg = [False, self.do_RV_GP.isChecked(), True, self.do_tra_GP.isChecked()]
        elif self.radioButton_transit_RV.isChecked():
            fit.rtg = [True,self.do_RV_GP.isChecked(), True, self.do_tra_GP.isChecked()]
        
        self.button_MCMC.setEnabled(False)
        self.statusBar().showMessage('MCMC in progress....')        
        # check if RV data is present
        if fit.type_fit["RV"] == True and fit.filelist.ndset <= 0:
             choice = QtGui.QMessageBox.information(self, 'Warning!',
             "Not possible to run MCMC if there are no RV data loaded. Please add your RV data first. Okay?", QtGui.QMessageBox.Ok)      
             self.button_MCMC.setEnabled(True)  
             self.statusBar().showMessage('') 

             return   
        
        ntran_data = 0
        for i in range(0,10,1):         
            ntran_data += len(fit.tra_data_sets[i]) 
            
        if fit.type_fit["Transit"] == True  and ntran_data == 0:
             choice = QtGui.QMessageBox.information(self, 'Warning!',
             "Not possible to run MCMC if there are no transit data loaded. Please add your transit data first. Okay?", QtGui.QMessageBox.Ok)      
             self.button_MCMC.setEnabled(True)  
             self.statusBar().showMessage('') 

             return             
    

        choice = QtGui.QMessageBox.information(self, 'Warning!',
                                            "This will run in the background and may take some time. Results are printed in the 'Stdout/Stderr' tab. Okay?",
                                            QtGui.QMessageBox.Cancel | QtGui.QMessageBox.Ok)       
         
        if choice == QtGui.QMessageBox.Cancel:
            self.statusBar().showMessage('') 
            self.button_MCMC.setEnabled(True)
            return        

        self.set_tra_ld()        
        self.check_bounds()
        self.check_priors_nr() 
        fit.model_npoints = self.points_to_draw_model.value()
        
        self.tabWidget_helper.setCurrentWidget(self.tab_info)
        
        
        if self.use_percentile_level.isChecked():
            fit.percentile_level = self.percentile_level.value()
        else:
            fit.percentile_level = 68.3
           
        
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
        
 
        self.check_model_params()
        self.check_mcmc_params()
      
        fit = rv.run_mcmc(fit, burning_ph=self.burning_phase.value(), mcmc_ph=self.mcmc_phase.value(), threads=int(self.N_threads.value()), output=False,
        fileoutput=self.save_samples.isChecked(),save_means=self.adopt_mcmc_means_as_par.isChecked(), save_mode=self.adopt_mcmc_mode_as_par.isChecked(),
        save_maxlnL=self.adopt_best_lnL_as_pars.isChecked(),save_sampler=True)
        
    
        self.button_MCMC.setEnabled(True)            
 
    def change_mcmc_samples_file_name(self):
        global fit
        
        output_file = QtGui.QFileDialog.getSaveFileName(self, 'path and name of the mcmc samples', '', '')
        
        if output_file[0] != '':
            fit.mcmc_sample_file = output_file[0] 
            self.mcmc_samples_change_name.setText(output_file[0])
        else:
            return

    def check_mcmc_params(self):
        global fit
        fit.gaussian_ball = self.init_gauss_ball.value() 
        fit.nwalkers_fact = int(self.nwalkers_fact.value()) 

    def check_model_params(self):
        global fit
        fit.time_step_model = self.time_step_model.value()
        fit.dyn_model_accuracy = self.dyn_model_accuracy.value()
        fit.master_timeout = self.master_timeout.value()   


    def force_mcmc_check_box(self):
        if self.make_corner_plot.isChecked():
            self.save_samples.setChecked(True)
 
                       
################################## Cornerplot #######################################

    def worker_cornerplot_complete(self):
        global fit  
        self.statusBar().showMessage('') 
        self.button_make_mcmc_cornerplot.setEnabled(True)
        self.button_make_nest_cornerplot.setEnabled(True)
       


    def worker_cornerplot(self, type_plot = "mcmc"):
        global fit  
        
        self.button_make_mcmc_cornerplot.setEnabled(False)
        self.button_make_nest_cornerplot.setEnabled(False)
       
        self.statusBar().showMessage('Cornerplot in progress....')        
        # check if RV data is present
        if type_plot == "mcmc":
            samp_file = fit.mcmc_sample_file
            type_samp = "MCMC"
        elif type_plot == "nest":
            samp_file = fit.nest_sample_file
            type_samp = "Nest. Samp."            
        
        if not os.path.exists(samp_file):
             choice = QtGui.QMessageBox.information(self, 'Warning!',
             "%s file not found. Generate one and try again?"%type_samp, QtGui.QMessageBox.Ok)      
             self.button_make_mcmc_cornerplot.setEnabled(True)
             self.button_make_nest_cornerplot.setEnabled(True)
             
             return  
 
 
        # Pass the function to execute
        worker_cor = Worker(lambda: self.make_cornerplot(type_plot = type_plot)) # Any other args, kwargs are passed to the run  
        # Execute
        worker_cor.signals.finished.connect(self.worker_cornerplot_complete)
        
        # worker.signals.result.connect(self.print_output)
        #worker.signals.finished.connect(self.thread_complete)
       # worker.signals.progress.connect(self.progress_fn)
        self.threadpool.start(worker_cor)
 
            
    def make_cornerplot(self,type_plot = 'mcmc'):
        global fit
        rv.cornerplot(fit, fileinput=True, type_plot = type_plot )
            
            
      
    def change_corner_plot_file_name(self, type_plot = "mcmc"):
        global fit
        
        output_file = QtGui.QFileDialog.getSaveFileName(self, 'path and name of the corener plot', '', 'Data (*.png)')
        if output_file[0] != '':
            if type_plot == "mcmc":
                fit.mcmc_corner_plot_file = output_file[0] 
                self.mcmc_corner_plot_change_name.setText(output_file[0])   
            elif type_plot == "nest":
                fit.nest_corner_plot_file = output_file[0] 
                self.nest_corner_plot_change_name.setText(output_file[0])               
            
       
################################# data inspector ###################################  
        
        
        
    def load_data_inspect(self): 
        global fit
         
        #path = self.tree_view_tab.listview.model().filePath(index)
        path = self.inspector_file
        if path == '':
            return 
        
        filename, file_extension = os.path.splitext(path)  
            
        if file_extension == '.vels':
            fit.add_dataset(self.file_from_path(path), str(path),0.0,1.0)
            self.init_fit()            
            self.update_use_from_input_file()            
            self.update_use()
            self.update_params()
            self.update_RV_file_buttons()        
 
        elif file_extension == '.act':
            
            for i in range(10):
                if len(fit.act_data_sets[i]) == 0:    
                    but_ind = i +1
            
                    fit.add_act_dataset('test', str(path),act_idset =but_ind-1)
         
        
                    self.update_act_file_buttons()
                    self.update_activity_gls_plots(but_ind-1)
                    self.buttonGroup_activity_data.button(but_ind).setText(self.file_from_path(path))
                    return
            
        elif  file_extension == '.tran':
            for i in range(10):
                if len(fit.tra_data_sets[i]) == 0:                 
                    but_ind = i +1
 
                    fit.add_transit_dataset('test', str(path),tra_idset =but_ind-1)
                    self.update_use_from_input_file()            
                    self.update_use()
                    self.update_gui_params()           
                    self.update_params()
                    self.update_tra_file_buttons()
                    self.buttonGroup_transit_data.button(but_ind).setText(self.file_from_path(path))
                    return
           
        else: 
            return
        
    def cross_data_inspect(self):        
        if self.inpector_plot_cross_hair.isChecked():
            self.cross_hair(pdi,log=False)   
            
            
    def plot_data_inspect(self, index):
        global fit, colors, pdi 
        # self.sender() == self.treeView
        # self.sender().model() == self.fileSystemModel

        path = self.sender().model().filePath(index)
 
   
        pdi.plot(clear=True,)  
        
        try:    
            x     = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [0])
            y     = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [1])
            y_err = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [2]) 
        except:
            #pdi.addLine(x=[0,1], y=0, pen=pg.mkPen('#ff9933', width=0.8)) 
            pdi.setLabel('bottom', 'x', units='',  **{'font-size':'9pt'})
            pdi.setLabel('left',   'y', units='',  **{'font-size':'9pt'})
            return
 
        pdi.addLine(x=None, y=np.mean(y), pen=pg.mkPen('#ff9933', width=0.8))   
 
        if self.insp_data_size.value() > 2:
            symbolsize = self.insp_data_size.value() -2
        else:
            symbolsize = self.insp_data_size.value() 

    
        pdi.plot(x,y,             
        pen=None, #{'color': colors[i], 'width': 1.1},
        symbol='o',
        symbolPen={'color': fit.colors[0], 'width': 1.1},
        symbolSize=symbolsize,enableAutoRange=True,viewRect=True,
        symbolBrush=fit.colors[0]
        )  
               
        err_ = pg.ErrorBarItem(x=x, y=y, symbol='o', 
        top = y_err, bottom = y_err,
        #height=y_err, 
        beam=0.0, pen=fit.colors[0])   
     
        pdi.addItem(err_)
        
        pdi.autoRange()
 
        
        filename, file_extension = os.path.splitext(path)  
            
        if file_extension == '.vels':
            pdi.setLabel('bottom', 'BJD', units='d',  **{'font-size':'9pt'})
            pdi.setLabel('left',   'RV', units='m/s',  **{'font-size':'9pt'})         
 
        elif file_extension == '.act':
            pdi.setLabel('bottom', 'BJD', units='d',  **{'font-size':'9pt'})
            pdi.setLabel('left',   'y', units='',  **{'font-size':'9pt'})    
            
        elif file_extension == '.tran':      
            pdi.setLabel('bottom', 'BJD', units='d',  **{'font-size':'9pt'})
            pdi.setLabel('left',   'flux', units='',  **{'font-size':'9pt'})

        else:      
            pdi.setLabel('bottom', 'x', units='',  **{'font-size':'9pt'})
            pdi.setLabel('left',   'y', units='',  **{'font-size':'9pt'})
        
 
        self.inspector_file = path     
        self.data_insp_print_info.clicked.connect(lambda: self.print_info_for_object(self.stat_info(x,y,y_err,path)))   

        #self.data_insp_load_data.clicked.connect(lambda: self.load_data_inspect(path))
 
        self.cross_data_inspect()


        
        
    def stat_info(self,x,y,y_err,path):
 
 
        ################## text generator #################
        text_info = """ 
"""
        text_info = text_info +"""%s
-----------------------------------  
"""%path    
        text_info = text_info +"""
N data          :  %s        
        
first epoch     :  %.3f
last epoch      :  %.3f
time span       :  %.3f

min. value      :  %.4f
max. value      :  %.4f
end-to-end      :  %.4f
mean            :  %.4f
median          :  %.4f
r.m.s.          :  %.4f

min error       :  %.4f
max error       :  %.4f
mean error      :  %.4f
median error    :  %.4f

"""%(len(x), x[0], x[-1], x[-1]-x[0], np.min(y), np.max(y), np.max(y)-np.min(y), np.mean(y),  np.median(y), np.sqrt(np.mean(np.square(y))),
np.min(y_err), np.max(y_err),   np.mean(y_err),  np.median(y_err))       
      
     
        return text_info  
              
        
             
            
       
############################# Dispatcher #####################################  

    def fit_dispatcher(self, init=False):
        global fit
   
        if self.radioButton_RV.isChecked():
            fit.rtg = [True,self.do_RV_GP.isChecked(),False, self.do_tra_GP.isChecked()]            
            if(init):
                self.worker_RV_fitting(ff=0,m_ln=True, init = init )  
                #print('test')
            else:
                self.worker_RV_fitting(m_ln=self.amoeba_radio_button.isChecked())  
                               
        elif self.radioButton_transit.isChecked(): 
            fit.rtg = [False,False,True, self.do_tra_GP.isChecked()]
            if(init):             
                self.worker_transit_fitting(ff=0)  
            else:
                self.worker_transit_fitting()
                                               
        elif self.radioButton_transit_RV.isChecked():
            
            fit.rtg=[True,self.do_RV_GP.isChecked(),True, self.do_tra_GP.isChecked()]
            if(init):
                self.worker_transit_fitting(ff=0 )  
            else:
                self.worker_transit_fitting()

 
    
 
###########################  GUI events #############################            
 
    def mute_boxes(self):
        
        if batman_not_found == True:
            self.radioButton_transit.setEnabled(False)          
            self.radioButton_transit_RV.setEnabled(False)       
            print("""
You dont have batman installed or you have the wrong 'batman'! 
Therefore, you cannnot use the transit modelling. 
Please install via ' pip install batman-package' not 'pip install batman' and try again. 
For more info on the used 'batman' in the 'Exo-Striker', please check 'Help --> Credits' 
""")
            self.tabWidget_helper.setCurrentWidget(self.tab_info)

        ######### TESTS!!!!!!!!!!!###########
        
        self.mute_boxes_dyn()
        
        if self.radioButton_transit_RV.isChecked():
            
            ma_flag = False
            t0_flag = True
            K_flag = True
            pl_rad_flag = True
            a_sol_flag = True
            fit.type_fit["RV"] = True           
            fit.type_fit["Transit"] = True 
            
        elif self.radioButton_transit.isChecked():            
            
            ma_flag = False
            t0_flag = True
            K_flag = False
            pl_rad_flag = True
            a_sol_flag = True
            fit.type_fit["RV"] = False           
            fit.type_fit["Transit"] = True 
            
        elif self.radioButton_RV.isChecked():            
            
            ma_flag = True
            t0_flag = False
            K_flag = True
            pl_rad_flag = False
            a_sol_flag = False                        
            fit.type_fit["RV"] = True           
            fit.type_fit["Transit"] = False 

            
        self.ma1.setEnabled(ma_flag)
        self.use_ma1.setEnabled(ma_flag)
        self.ma2.setEnabled(ma_flag)
        self.use_ma2.setEnabled(ma_flag)           
        self.ma3.setEnabled(ma_flag)
        self.use_ma3.setEnabled(ma_flag)  
        self.ma4.setEnabled(ma_flag)
        self.use_ma4.setEnabled(ma_flag)
        self.ma5.setEnabled(ma_flag)
        self.use_ma5.setEnabled(ma_flag)           
        self.ma6.setEnabled(ma_flag)
        self.use_ma6.setEnabled(ma_flag)  
        self.ma7.setEnabled(ma_flag)
        self.use_ma7.setEnabled(ma_flag)
        self.ma8.setEnabled(ma_flag)
        self.use_ma8.setEnabled(ma_flag)           
        self.ma9.setEnabled(ma_flag)
        self.use_ma9.setEnabled(ma_flag)  
        
        self.t0_1.setEnabled(t0_flag)
        self.use_t0_1.setEnabled(t0_flag)
        self.t0_2.setEnabled(t0_flag)
        self.use_t0_2.setEnabled(t0_flag)           
        self.t0_3.setEnabled(t0_flag)
        self.use_t0_3.setEnabled(t0_flag)  
        self.t0_4.setEnabled(t0_flag)
        self.use_t0_4.setEnabled(t0_flag)
        self.t0_5.setEnabled(t0_flag)
        self.use_t0_5.setEnabled(t0_flag)           
        self.t0_6.setEnabled(t0_flag)
        self.use_t0_6.setEnabled(t0_flag)              
        self.t0_7.setEnabled(t0_flag)
        self.use_t0_7.setEnabled(t0_flag)
        self.t0_8.setEnabled(t0_flag)
        self.use_t0_8.setEnabled(t0_flag)           
        self.t0_9.setEnabled(t0_flag)
        self.use_t0_9.setEnabled(t0_flag)              
        
        self.K1.setEnabled(K_flag)
        self.use_K1.setEnabled(K_flag)
        self.K2.setEnabled(K_flag)
        self.use_K2.setEnabled(K_flag)           
        self.K3.setEnabled(K_flag)
        self.use_K3.setEnabled(K_flag)        
        self.K4.setEnabled(K_flag)
        self.use_K4.setEnabled(K_flag)
        self.K5.setEnabled(K_flag)
        self.use_K5.setEnabled(K_flag)           
        self.K6.setEnabled(K_flag)
        self.use_K6.setEnabled(K_flag)                    
        self.K7.setEnabled(K_flag)
        self.use_K7.setEnabled(K_flag)
        self.K8.setEnabled(K_flag)
        self.use_K8.setEnabled(K_flag)           
        self.K9.setEnabled(K_flag)
        self.use_K9.setEnabled(K_flag)                    

        self.pl_rad_1.setEnabled(pl_rad_flag)
        self.use_pl_rad_1.setEnabled(pl_rad_flag)
        self.pl_rad_2.setEnabled(pl_rad_flag)
        self.use_pl_rad_2.setEnabled(pl_rad_flag)           
        self.pl_rad_3.setEnabled(pl_rad_flag)
        self.use_pl_rad_3.setEnabled(pl_rad_flag) 
        self.pl_rad_4.setEnabled(pl_rad_flag)
        self.use_pl_rad_4.setEnabled(pl_rad_flag)
        self.pl_rad_5.setEnabled(pl_rad_flag)
        self.use_pl_rad_5.setEnabled(pl_rad_flag)           
        self.pl_rad_6.setEnabled(pl_rad_flag)
        self.use_pl_rad_6.setEnabled(pl_rad_flag)             
        self.pl_rad_7.setEnabled(pl_rad_flag)
        self.use_pl_rad_7.setEnabled(pl_rad_flag)
        self.pl_rad_8.setEnabled(pl_rad_flag)
        self.use_pl_rad_8.setEnabled(pl_rad_flag)           
        self.pl_rad_9.setEnabled(pl_rad_flag)
        self.use_pl_rad_9.setEnabled(pl_rad_flag)                   
  
        self.a_sol_1.setEnabled(a_sol_flag)
        self.use_a_sol_1.setEnabled(a_sol_flag)
        self.a_sol_2.setEnabled(a_sol_flag)
        self.use_a_sol_2.setEnabled(a_sol_flag)           
        self.a_sol_3.setEnabled(a_sol_flag)
        self.use_a_sol_3.setEnabled(a_sol_flag)                        
        self.a_sol_4.setEnabled(a_sol_flag)
        self.use_a_sol_4.setEnabled(a_sol_flag)
        self.a_sol_5.setEnabled(a_sol_flag)
        self.use_a_sol_5.setEnabled(a_sol_flag)           
        self.a_sol_6.setEnabled(a_sol_flag)
        self.use_a_sol_6.setEnabled(a_sol_flag)                 
        self.a_sol_7.setEnabled(a_sol_flag)
        self.use_a_sol_7.setEnabled(a_sol_flag)
        self.a_sol_8.setEnabled(a_sol_flag)
        self.use_a_sol_8.setEnabled(a_sol_flag)           
        self.a_sol_9.setEnabled(a_sol_flag)
        self.use_a_sol_9.setEnabled(a_sol_flag)                
 
    
    def set_gr_flag(self):
        global fit
        
        if self.radioButton_omega_dot_free.isChecked():
            fit.gr_flag = False
        elif self.deattach_omega_dot.isChecked() and self.radioButton_omega_dot_GR.isChecked():
            fit.gr_flag = True
        else:
            fit.gr_flag = False

    def mute_boxes_dyn(self):
                
        if self.radioButton_Keplerian.isChecked() and self.radioButton_RV.isChecked()==True:
            
            om_flag = False
            incl_flag = False
            Dom_flag = False
            if self.deattach_omega_dot.isChecked():
                self.set_gr_flag()
                if not self.radioButton_omega_dot_GR.isChecked(): 
                    om_flag = True                                

        elif self.radioButton_Keplerian.isChecked() and self.radioButton_RV.isChecked()==False:

            om_flag = False
            incl_flag = True
            Dom_flag = False
                                    
        elif self.radioButton_Dynamical.isChecked():
            
            om_flag = False
            incl_flag = True
            Dom_flag = True 
            
        self.om_dot_1.setEnabled(om_flag)
        self.use_om_dot_1.setEnabled(om_flag)
        self.om_dot_2.setEnabled(om_flag)
        self.use_om_dot_2.setEnabled(om_flag)           
        self.om_dot_3.setEnabled(om_flag)
        self.use_om_dot_3.setEnabled(om_flag) 
        self.om_dot_4.setEnabled(om_flag)
        self.use_om_dot_4.setEnabled(om_flag)
        self.om_dot_5.setEnabled(om_flag)
        self.use_om_dot_5.setEnabled(om_flag)           
        self.om_dot_6.setEnabled(om_flag)
        self.use_om_dot_6.setEnabled(om_flag)             
        self.om_dot_7.setEnabled(om_flag)
        self.use_om_dot_7.setEnabled(om_flag)
        self.om_dot_8.setEnabled(om_flag)
        self.use_om_dot_8.setEnabled(om_flag)           
        self.om_dot_9.setEnabled(om_flag)
        self.use_om_dot_9.setEnabled(om_flag)             
                                
        self.incl1.setEnabled(incl_flag)  
        self.use_incl1.setEnabled(incl_flag)  
        self.incl2.setEnabled(incl_flag)  
        self.use_incl2.setEnabled(incl_flag)              
        self.incl3.setEnabled(incl_flag)  
        self.use_incl3.setEnabled(incl_flag)  
        self.incl4.setEnabled(incl_flag)  
        self.use_incl4.setEnabled(incl_flag)  
        self.incl5.setEnabled(incl_flag)  
        self.use_incl5.setEnabled(incl_flag)  
        self.incl6.setEnabled(incl_flag)  
        self.use_incl6.setEnabled(incl_flag)  
        self.incl7.setEnabled(incl_flag)  
        self.use_incl7.setEnabled(incl_flag)  
        self.incl8.setEnabled(incl_flag)  
        self.use_incl8.setEnabled(incl_flag)  
        self.incl9.setEnabled(incl_flag)  
        self.use_incl9.setEnabled(incl_flag)    
        
        self.Omega1.setEnabled(Dom_flag)  
        self.use_Omega1.setEnabled(Dom_flag)  
        self.Omega2.setEnabled(Dom_flag)  
        self.use_Omega2.setEnabled(Dom_flag)             
        self.Omega3.setEnabled(Dom_flag)  
        self.use_Omega3.setEnabled(Dom_flag)  
        self.Omega4.setEnabled(Dom_flag)  
        self.use_Omega4.setEnabled(Dom_flag)
        self.Omega5.setEnabled(Dom_flag)  
        self.use_Omega5.setEnabled(Dom_flag)  
        self.Omega6.setEnabled(Dom_flag)  
        self.use_Omega6.setEnabled(Dom_flag)             
        self.Omega7.setEnabled(Dom_flag)  
        self.use_Omega7.setEnabled(Dom_flag)  
        self.Omega8.setEnabled(Dom_flag)  
        self.use_Omega8.setEnabled(Dom_flag) 
        self.Omega9.setEnabled(Dom_flag)  
        self.use_Omega9.setEnabled(Dom_flag)             
            
       
    def keyPressEvent(self, event):
        global fit
        if event.key() in (QtCore.Qt.Key_Enter, QtCore.Qt.Key_Return):
            self.update_veiw()
            #self.update_use()
            #self.update_params() 
            #self.init_fit()
            #self.fit_dispatcher(init=True)
            return

    def update_veiw(self):
        global fit
 
        self.update_use()
        self.update_params() 
        #self.init_fit()
        self.fit_dispatcher(init=True)
        return


############################# Tab selector (not ready) ################################  

    def tab_selected(self,ind):

        if ind == 4:
            self.update_activity_data_plots(self.comboBox_act_data.currentIndex())        
        if ind == 5:
            self.update_correlations_data_plots()
 
    
  
 
#############################  Color control ################################  
    def update_color_picker_tra(self):
        global fit
        
        font = QtGui.QFont()
        font.setPointSize(11)
        font.setBold(False)
        #font.setWeight(75)
        
 
        for i in range(11):
            self.buttonGroup_color_picker_tra.button(i+1).setStyleSheet("color: %s;"%fit.tra_colors[i])
            self.buttonGroup_color_picker_tra.button(i+1).setFont(font)              
        for i in range(10):    
            self.buttonGroup_symbol_picker_tra.button(i+1).setStyleSheet("color: %s;"%fit.tra_colors[i])  
            self.buttonGroup_symbol_picker_tra.button(i+1).setText(fit.pyqt_symbols_tra[i]) 
            self.buttonGroup_symbol_picker_tra.button(i+1).setFont(font)         

          
    def get_color_tra(self):
        global fit
        
        but_ind = self.buttonGroup_color_picker_tra.checkedId()       
        colorz = self.colorDialog.getColor()       
        
        #QtGui.QColorDialog.setOption(QtGui.QColorDialog.ShowAlphaChannel,True)
        #colorz = QtGui.QColorDialog.getColor()
        #print(but_ind-1)
        #print(fit.colors[but_ind-1])       
        if colorz.isValid():
            fit.tra_colors[but_ind-1]=colorz.name()   
            self.update_color_picker_tra()
            self.update_act_file_buttons()      
            self.update_RV_file_buttons() 
            self.update_tra_file_buttons() 
            self.update_RV_plots() 
            self.update_extra_plots()            
            self.update_transit_plots() 
            #self.update_activity_data_plots() 
            #self.update_activity_gls_plots()     
        else:
            return



    def update_color_picker(self):
        global fit
        
        font = QtGui.QFont()
        font.setPointSize(11)
        font.setBold(False)
        #font.setWeight(75)
        
 
        for i in range(11):
            self.buttonGroup_color_picker.button(i+1).setStyleSheet("color: %s;"%fit.colors[i])
            self.buttonGroup_color_picker.button(i+1).setFont(font)              
        for i in range(10):    
            self.buttonGroup_symbol_picker.button(i+1).setStyleSheet("color: %s;"%fit.colors[i])  
            self.buttonGroup_symbol_picker.button(i+1).setText(fit.pyqt_symbols_rvs[i]) 
            self.buttonGroup_symbol_picker.button(i+1).setFont(font)              
            
          
    def get_color(self):
        global fit
        
        but_ind = self.buttonGroup_color_picker.checkedId()       
        colorz = self.colorDialog.getColor()       
        
        #QtGui.QColorDialog.setOption(QtGui.QColorDialog.ShowAlphaChannel,True)
        #colorz = QtGui.QColorDialog.getColor()
        print(but_ind-1)
        print(fit.colors[but_ind-1])       
        if colorz.isValid():
            fit.colors[but_ind-1]=colorz.name()   
            self.update_color_picker()
            self.update_act_file_buttons()      
            self.update_RV_file_buttons() 
            self.update_tra_file_buttons() 
            self.update_RV_plots() 
            self.update_extra_plots()            
            self.update_transit_plots() 
            #self.update_activity_data_plots() 
            #self.update_activity_gls_plots()     
        else:
            return
        
        



############################# Symbol controls ################################  


    def get_symbol_tra(self):
        global fit
 
        but_ind = self.buttonGroup_symbol_picker_tra.checkedId()   
        but_n = self.dialog_symbols.get_radio()
            
        if but_n != None:
            fit.pyqt_symbols_tra[but_ind-1] = symbols[but_n-1]
            self.update_color_picker_tra()
            self.update_act_file_buttons()      
            self.update_RV_file_buttons() 
            self.update_tra_file_buttons() 
            self.update_RV_plots() 
            self.update_extra_plots()
            self.update_transit_plots()     
        else:
            return    
        
    def check_tra_symbol_sizes(self):
        global fit
       
       # for i in range(10):
        fit.pyqt_symbols_size_tra[0] = self.trans_data_size_1.value()
        fit.pyqt_symbols_size_tra[1] = self.trans_data_size_2.value()
        fit.pyqt_symbols_size_tra[2] = self.trans_data_size_3.value()
        fit.pyqt_symbols_size_tra[3] = self.trans_data_size_4.value()
        fit.pyqt_symbols_size_tra[4] = self.trans_data_size_5.value()
        fit.pyqt_symbols_size_tra[5] = self.trans_data_size_6.value()
        fit.pyqt_symbols_size_tra[6] = self.trans_data_size_7.value()
        fit.pyqt_symbols_size_tra[7] = self.trans_data_size_8.value()
        fit.pyqt_symbols_size_tra[8] = self.trans_data_size_9.value()
        fit.pyqt_symbols_size_tra[9] = self.trans_data_size_10.value()




            
    def get_symbol(self):
        global fit
 
        but_ind = self.buttonGroup_symbol_picker.checkedId()   
        but_n = self.dialog_symbols.get_radio()
            
        if but_n != None:
            fit.pyqt_symbols_rvs[but_ind-1] = symbols[but_n-1]
            self.update_color_picker()
            self.update_act_file_buttons()      
            self.update_RV_file_buttons() 
            self.update_tra_file_buttons() 
            self.update_RV_plots() 
            self.update_extra_plots()
            self.update_transit_plots()     
        else:
            return    
        

        
      
    def check_RV_symbol_sizes(self):
        global fit
       
       # for i in range(10):
        fit.pyqt_symbols_size_rvs[0] = self.rv_data_size_1.value()
        fit.pyqt_symbols_size_rvs[1] = self.rv_data_size_2.value()
        fit.pyqt_symbols_size_rvs[2] = self.rv_data_size_3.value()
        fit.pyqt_symbols_size_rvs[3] = self.rv_data_size_4.value()
        fit.pyqt_symbols_size_rvs[4] = self.rv_data_size_5.value()
        fit.pyqt_symbols_size_rvs[5] = self.rv_data_size_6.value()
        fit.pyqt_symbols_size_rvs[6] = self.rv_data_size_7.value()
        fit.pyqt_symbols_size_rvs[7] = self.rv_data_size_8.value()
        fit.pyqt_symbols_size_rvs[8] = self.rv_data_size_9.value()
        fit.pyqt_symbols_size_rvs[9] = self.rv_data_size_10.value()


################################## View Actions #######################################

    def grab_screen(self):
        p = QtWidgets.QWidget.grab(self)
       # p.scaled(40, 40)
        #painter = QtGui.QPainter(p)
       # painter.setRenderHint(QtGui.QPainter.TextAntialiasing, True)
        
        #painter.setRenderHint(QtGui.QPainter.Antialiasing, True)       
        #p.scaled(8000, 8000, transformMode=QtCore.Qt.SmoothTransformation)
       # p.setDevicePixelRatio(2.0)
        #screen = QtWidgets.QApplication.primaryScreen()
        #p =  screen.grabWindow(0)
        #painter.setRenderHint(QtGui.QPainter.Antialiasing, self.params['antialias'])
       # filename = QtGui.QFileDialog.getSaveFileName(self, 'Save image', filter="PNG(*.png);; JPEG(*.jpg)")
        #p.save(filename[0], 'jpg')        
        #label.setPixmap(p)        # just for fun :)
        img, _ = QtGui.QFileDialog.getSaveFileName(self,"Save image",
                                            filter="PNG(*.png);; JPEG(*.jpg)")
        if img[-3:] == "png":
            p.save(img, "png")
        elif img[-3:] == "jpg":
            p.save(img, "jpg")       
            
            
            
            
        
    def print_f_test_stat(self):
        global fit

        rv.f_test(fit, alpha = 0.01)
        self.tabWidget_helper.setCurrentWidget(self.tab_info)

    def get_latex_table(self):
        global fit       

        self.console_widget.print_text("rv.latex_pl_param_table(fit, width = 10, precision = 2, asymmetric = False, file_name='test.tex',path='./')", before_prompt=False)  
#        self.console_widget.execute_command("rv.latex_pl_param_table(fit, width = 10, precision = 2, asymmetric = False, file_name='test.tex',path='./')")  
        self.tabWidget_helper.setCurrentWidget(self.tab_shells) 
        self.terminal_embeded.setCurrentWidget(self.console_widget)
        
    def get_RV_model(self):
        global fit       

        self.console_widget.print_text("rv.export_RV_model(fit, file='RV_model.txt', width = 10, precision = 4)", before_prompt=False)  
        self.tabWidget_helper.setCurrentWidget(self.tab_shells) 
        self.terminal_embeded.setCurrentWidget(self.console_widget)        
        
    def get_RV_data(self):
        global fit       

        self.console_widget.print_text("rv.export_RV_data(fit, [0], file='RV_data.txt',  jitter=False, o_c=False, print_data=False, width = 10, precision = 3)", before_prompt=False)  
        self.tabWidget_helper.setCurrentWidget(self.tab_shells) 
        self.terminal_embeded.setCurrentWidget(self.console_widget)   
        
        
    def get_orb_evol(self):
        global fit       

        self.console_widget.print_text("rv.export_orbital_evol(fit, file='planet_N.txt', planet = 1, width = 10, precision = 6)", before_prompt=False)  
        self.tabWidget_helper.setCurrentWidget(self.tab_shells) 
        self.terminal_embeded.setCurrentWidget(self.console_widget)           
        


################################## System #######################################
            
    def quit(self):
        global fit            
        #os.system("rm temp*.vels")
        choice = QtGui.QMessageBox.information(self, 'Warning!',
                                            "Do you want to save the session before you Quit?",
                                            QtGui.QMessageBox.Cancel | QtGui.QMessageBox.No | QtGui.QMessageBox.Yes)  
         
        if choice == QtGui.QMessageBox.No:
            
            # to remove the /tmp created files. Unfortunatly this leads to problems later (TBF)
           # for i in range(fit.filelist.ndset): 
           #     dirname, basename = os.path.split(fit.filelist.files[i].path)
           #     os.system('rm -r %s'%dirname) 
            
            self.term_emb.close()            
            self.close()
        elif choice == QtGui.QMessageBox.Yes:
            self.save_session()
            self.term_emb.close()            
            self.close()            
        elif choice == QtGui.QMessageBox.Cancel:
            return    
        
    def file_from_path(self, path):
        head, tail = ntpath.split(path)
        return tail or ntpath.basename(head)
        
       
    def check_settings(self):
        global fit

        mixed_fit_use = [self.mix_pl_1,self.mix_pl_2,self.mix_pl_3,
                         self.mix_pl_4,self.mix_pl_5,self.mix_pl_6,
                         self.mix_pl_7,self.mix_pl_8,self.mix_pl_9 ]
        
        self.use_mix_fitting.setChecked(bool(fit.mixed_fit[0][0]))
        
        for i in range(9):
            mixed_fit_use[i].setChecked(bool(fit.mixed_fit[1][i]))
            
        self.time_step_model.setValue(fit.time_step_model)
        self.dyn_model_accuracy.setValue(fit.dyn_model_accuracy)
        self.dyn_model_to_kill.setValue(fit.dyn_model_to_kill)
        self.kep_model_to_kill.setValue(fit.kep_model_to_kill)
        self.master_timeout.setValue(fit.master_timeout)    
        
        
    def adopt_RV_GLS_param(self):
        global fit   
        
        mean_anomaly_from_gls = np.degrees((((fit.epoch - float(fit.gls_o_c.hpstat["T0"]) )% (fit.gls_o_c.hpstat["P"]) )/ (fit.gls_o_c.hpstat["P"]) ) * 2*np.pi)
 
        fit.add_planet(fit.gls_o_c.hpstat["amp"],fit.gls_o_c.hpstat["P"],0.0,0.0,mean_anomaly_from_gls -90.0,90.0,0.0)
        fit.use.update_use_planet_params_one_planet(fit.npl+1,True,True,True,True,True,False,False)   
       
        self.update_use_from_input_file()   
        self.update_use() 
        #self.update_params()  
       # self.update_gui_params()                  
        self.optimize_fit(0,m_ln=self.amoeba_radio_button.isChecked(),auto_fit = True)          
 
    def adopt_trans_TLS_param(self):
        global fit   
 
        
        mean_anomaly_from_tls = np.degrees((((fit.epoch - fit.tls_o_c.transit_times[0] )% (fit.tls_o_c.period) )/ (fit.tls_o_c.period) ) * 2*np.pi)
       
        fit.t0[fit.npl]     = fit.tls_o_c.transit_times[0]
        fit.pl_rad[fit.npl] = fit.stellar_radius*np.sqrt(1.0 - fit.tls_o_c.depth_mean[0])  # alternativly fit.tls_o_c.rp_rs ?
        fit.pl_a[fit.npl]   = 11.44
 
        fit.add_planet(10.0,fit.tls_o_c.period,0.0,0.0,mean_anomaly_from_tls,90.0,0.0)
        fit.use.update_use_planet_params_one_planet(fit.npl+1,True,True,True,True,True,False,False)   
       
 
            
        self.update_use_from_input_file()   
        self.update_use() 
       # self.update_params()  
        self.update_gui_params()                  
        self.radioButton_transit.setChecked(True)                    

        self.worker_transit_fitting(ff=0 ) 
        
#############################  TEST ZONE ################################  
        
#     def update_stellar_params(self):
#        global fit         
        
#        fit.stellar_radius = self.St_radius_input.value()
        
        
    def init_plot_corr(self):
        global fit  
  
        self.comboBox_samp_corr_1.clear()
        self.comboBox_samp_corr_2.clear()
      
        for i in range(len(fit.e_for_mcmc)):
            self.comboBox_samp_corr_1.addItem(fit.e_for_mcmc[i],i)
            self.comboBox_samp_corr_2.addItem(fit.e_for_mcmc[i],i)

            #self.comboBox_pl_2.setItemText(i, '%s'%(ses_list[i].name))
        self.comboBox_samp_corr_1.setCurrentIndex(0)
        self.comboBox_samp_corr_2.setCurrentIndex(0)



    def update_plot_corr(self):
        global fit,pcor
 
        corr1_ind = self.comboBox_samp_corr_1.currentIndex()
        corr2_ind = self.comboBox_samp_corr_2.currentIndex()
 
        if corr1_ind ==-1 or corr2_ind ==-1 or len(fit.e_for_mcmc) ==0 or fit.sampler==None:
            return
        #else:
       #     last_stable = min(len(fit.evol_p[pl1_ind]),len(fit.evol_p[pl2_ind]))
        
 
        pcor.plot(clear=True,)
        pcor.plot(fit.sampler.samples[:,corr1_ind], fit.sampler.samples[:,corr2_ind] ,pen=None, #{'color': colors[i], 'width': 1.1},
        symbol='o',
        symbolPen={'color': 'b', 'width': 1},
        symbolSize=1,enableAutoRange=True,viewRect=True,
        symbolBrush='b'
        )  

        pcor.setLabel('bottom', '%s'%fit.e_for_mcmc[corr1_ind], units='',  **{'font-size':'9pt'}) 
        pcor.setLabel('left', '%s'%fit.e_for_mcmc[corr2_ind], units='',  **{'font-size':'9pt'}) 
        

        ## generate empty curves
        #curves = []
        #levels = np.linspace(0, 1, 3)
        #for i in range(len(levels)):
        #    v = levels[i]
            ## generate isocurve with automatic color selection
        #    c = pg.IsocurveItem(level=v, pen=(i, len(levels)*1.5))
           # c.setParentItem(pcor)  ## make sure isocurve is always correctly displayed over image
        #    c.setZValue(10)
        #    curves.append(c)        
        
        
        
        #c = pg.IsocurveItem(level=0.6, pen=('r'))#i, len(levels)*1.5))
        #c.setParentItem(pcor)  ## make sure isocurve is always correctly displayed over image
       # c.setZValue(10)
        #curves.append(c)  
        

    def set_kp_ma(self):
        global fit  
 
        if self.radioButton_KP.isChecked():
            
            self.label_K_mass_arb.setText("K [m/s]")
            self.label_P_a_arb.setText("P [day]")
        else:
            
            self.label_K_mass_arb.setText("mass [Mj]")
            self.label_P_a_arb.setText("a [au]")
     
 
   # def layout_widgets(layout):
    #   return (layout.itemAt(i) for i in range(layout.count()))

    def rv_plot_phase_change(self):
        global fit        
        
        #RVphase = self.RV_phase_slider.value()
        #print(RVphase)
        #self.phase_plots(1, offset = RVphase)
        ind = self.comboBox_extra_plot.currentIndex()
        if ind+1 <= fit.npl:
            self.phase_plots(ind+1)
        else:
            return
        

    def set_RV_GP(self):
        global fit
        
        if self.use_GP_sho_kernel.isChecked():
            fit.gp_kernel = 'SHOKernel'  
        elif self.use_GP_rot_kernel.isChecked():
            fit.gp_kernel = 'RotKernel'
            
            
    def set_tra_GP(self):
        global fit
        
        if self.use_tra_GP_sho_kernel.isChecked():
            fit.tra_gp_kernel = 'SHOKernel'  
        elif self.use_tra_GP_rot_kernel.isChecked():
            fit.tra_gp_kernel = 'RotKernel'                 
            
            
    def set_use_GP(self):
        global fit            
            
        if  self.do_RV_GP.isChecked():
            fit.doGP = True
        else:
            fit.doGP = False
            
        if  self.do_tra_GP.isChecked():
            fit.tra_doGP = True
        else:
            fit.tra_doGP = False   

      
    def set_force_copl_incl(self):
        global fit   
        fit.copl_incl = self.force_copl_incl.isChecked()

    #def update_inspector(self):
   #     self.tree_view_tab.listview.clicked.connect(self.plot_data_inspect)
   
   
    def update_St_params(self):
        global fit  

        fit.stellar_mass     = self.St_mass_input.value()
        fit.stellar_mass_err = self.err_St_mass_input.value()         
        
        fit.stellar_radius     = self.St_radius_input.value()
        fit.stellar_radius_err = self.err_St_radius_input.value()    
        
        fit.stellar_luminosity     = self.St_lumin_input.value()
        fit.stellar_luminosity_err = self.err_St_lumin_input .value()           
        
        fit.stellar_Teff       = self.St_teff_input.value()
        fit.stellar_Teff_err   = self.err_St_teff_input.value()             
 
        fit.stellar_vsini       = self.St_vsini_input.value()
        fit.stellar_vsini_err   = self.err_St_vsini_input.value()       
        
        st_rot = rv.get_stellar_rotation(fit)
        kb1995 = rv.get_rv_scatter(fit)
        kb2011 = rv.get_rv_scatter(fit,use_kb2011=True)
        
        self.label_St_rot_value.setText("%.4f +/- %.4f [days]"%(st_rot[0],st_rot[1]))
        self.label_kb1995.setText("%.4f +/- %.4f [m/sec]"%(kb1995[0],kb1995[1]))
        self.label_kb2011.setText("%.4f +/- %.4f [m/sec]"%(kb2011[0],kb2011[1]))
       
        
        #fit.stellar_rotation = 25.0
       # fit.stellar_rotation_err = 0.0     
   
   
   
    def initialize_font(self):    
            
        self.font = QtGui.QFont()
        self.font.setPointSize(9)
        self.font.setBold(False)
  
################################################################################################
    

    def __init__(self):
        global fit

        QtWidgets.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
       # self.showMaximized()

        self.setupUi(self)
        self.initialize_font()

        
        self.initialize_buttons()
        self.initialize_plots()   
 
        self.initialize_color_dialog()               

        
        ###################### Console #############################
        self.console_widget = ConsoleWidget_embed(font_size = 9)
        # add the console widget to the user interface
        # push some variables to the console
        self.console_widget.push_vars({"rv": rv,
                                "np": np,
                                "fit": fit,
                                #"plt": plt,
                                #"clc": self.clc,
                                #'app': self
                                })        
        
        
        #self.console_widget = ConsoleWidget_embed(font_size = 10)
        
        self.terminal_embeded.addTab(self.console_widget, "Jupyter")
        
        ###################### Console #############################
        
       
       # self.terminal_embeded.addTab(self.tree_view_tab, "tree")
      
        
        if sys.platform[0:5] == "linux":
            self.term_emb = terminal.EmbTerminal()
            self.terminal_embeded.addTab(self.term_emb, "Bash shell")        
        self.terminal_embeded.addTab(pg_console.ConsoleWidget(), "pqg shell")  
        
        
        #self.gridLayout_116.addWidget(terminal.EmbTerminal())
 
        self.gridLayout_text_editor.addWidget(text_editor_es.MainWindow())       
        self.gridLayout_calculator.addWidget(calc.Calculator())  
        
        #################### data inspector ########################
                        
        self.tree_view_tab = Widget_tree()        
       # self.gridLayout_file_tree.setRowStretch(0, 6)
        self.gridLayout_file_tree.setRowStretch(1, 4)
        self.gridLayout_file_tree.addWidget(self.tree_view_tab)

        
        self.tree_view_tab.listview.clicked.connect(self.plot_data_inspect)
        self.data_insp_load_data.clicked.connect(self.load_data_inspect)  
        
    
        
        #################### stdout pipe ########################
       
        #if sys.version_info[0] == 3:
        self.pipe_text = MyDialog()
        self.gridLayout_stdout.addWidget(self.pipe_text)  
   
        #################### credits  ########################
    
        self.dialog = print_info(self)
        self.dialog_credits = print_info(self)
        self.dialog_chi_table = print_info(self)      

        self.buttonGroup_add_RV_data.buttonClicked.connect(self.showDialog_RV_input_file)
        self.buttonGroup_remove_RV_data.buttonClicked.connect(self.remove_RV_file)
 
        self.buttonGroup_activity_data.buttonClicked.connect(self.showDialog_act_input_file)
        self.buttonGroup_remove_activity_data.buttonClicked.connect(self.remove_act_file)     
        
        self.buttonGroup_transit_data.buttonClicked.connect(self.showDialog_tra_input_file)
        self.buttonGroup_remove_transit_data.buttonClicked.connect(self.remove_tra_file)         
        
        self.buttonGroup_use.buttonClicked.connect(self.update_use)
        self.buttonGroup_mixed_fitting.buttonClicked.connect(self.update_mixed_fitting)
        
        self.button_orb_evol.clicked.connect(self.worker_Nbody) 
        self.button_MCMC.clicked.connect(self.worker_mcmc)
       # self.button_nest_samp.clicked.connect(lambda: self.run_nest_samp())
        self.button_nest_samp.clicked.connect(self.worker_nest)
        
        self.run_orb_evol_arbitary.clicked.connect(self.worker_Nbody_arb) 
 
        
        self.button_make_mcmc_cornerplot.clicked.connect(lambda: self.worker_cornerplot(type_plot = "mcmc"))
        self.button_make_nest_cornerplot.clicked.connect(lambda: self.worker_cornerplot(type_plot = "nest"))
        
        self.mcmc_corner_plot_change_name.clicked.connect(lambda: self.change_corner_plot_file_name(type_plot = "mcmc"))
        self.nest_corner_plot_change_name.clicked.connect(lambda: self.change_corner_plot_file_name(type_plot = "nest"))
        
        self.mcmc_samples_change_name.clicked.connect(self.change_mcmc_samples_file_name)
        self.nest_samples_change_name.clicked.connect(self.change_nest_samples_file_name)
        
        self.comboBox_samp_corr_1.activated.connect(self.update_plot_corr)
        self.comboBox_samp_corr_2.activated.connect(self.update_plot_corr)        
        
        
        ########## RV fitting ########################
        
        self.button_init_fit.clicked.connect(lambda: self.fit_dispatcher(init=True))
        self.button_fit.clicked.connect(lambda: self.fit_dispatcher())        
        self.button_auto_fit.clicked.connect(lambda: self.run_auto_fit())
        self.minimize_1param()

        self.radioButton_Dynamical.toggled.connect(self.update_dyn_kep_flag)
        self.radioButton_Dynamical.toggled.connect(self.mute_boxes_dyn)

        self.radioButton_omega_dot_free.toggled.connect(self.mute_boxes_dyn)
#        self.radioButton_omega_dot_GR.toggled.connect(self.set_gr_flag)

 

        self.radioButton_Keplerian.toggled.connect(self.update_dyn_kep_flag)
        
        self.deattach_omega_dot.stateChanged.connect(self.mute_boxes_dyn)
                
        self.amoeba_radio_button.toggled.connect(self.update_RV_jitter_flag)
        self.lm_radio_button.toggled.connect(self.update_RV_jitter_flag)       
        self.radioButton_Keplerian.toggled.connect(self.mute_boxes_dyn)
        
        
        self.adopt_best_RV__o_c_GLS_per.clicked.connect(self.adopt_RV_GLS_param)
        
        self.adopt_tls_o_c_param.clicked.connect(self.adopt_trans_TLS_param)

        ############ View #################

        self.actiongrab_screen.triggered.connect(self.grab_screen) 
        self.actionprint_f_test_FAP.triggered.connect(self.print_f_test_stat)
        self.actionGet_LaTeX_table_with_parameters.triggered.connect(self.get_latex_table)
        self.actionGet_RV_model.triggered.connect(self.get_RV_model)
        self.actionGet_RV_data.triggered.connect(self.get_RV_data)     
        
        self.actionGet_Orb_evol.triggered.connect(self.get_orb_evol)   
        
        
        ############ Sessions #################
        
        self.actionNew_session.triggered.connect(self.new_session)
        self.actionOpen_session.triggered.connect(self.open_session)
        self.actionSave_session.triggered.connect(self.save_session)
        self.actionSave_multi_sesion.triggered.connect(self.save_sessions)
        self.actionOpen_multi_session.triggered.connect(self.open_sessions) 
        #self.comboBox_extra_plot.activated.connect(self.change_extra_plot)      
        self.comboBox_select_ses.activated.connect(self.select_session)
        #self.comboBox_select_ses.lineEdit.editingFinished.connect(self.change_session_label)
        self.comboBox_select_ses.currentTextChanged.connect(self.change_session_label)
        
        self.session_list()
        self.new_ses.clicked.connect(self.getNewses)
        self.copy_ses.clicked.connect(self.cop_ses)
        self.remove_ses.clicked.connect(self.rem_ses)
  
        self.actionvisit_TRIFON_on_GitHub.triggered.connect(lambda: webbrowser.open('https://github.com/3fon3fonov/trifon'))    
        self.actionCredits.triggered.connect(lambda: self.print_info_credits())
        self.actionConfidence_Intervals_Table.triggered.connect(lambda: self.print_chi_table())

        ############### Orb. Evol. plotting ####################
        self.comboBox_pl_1.activated.connect(self.plot_delta_omega)
        self.comboBox_pl_2.activated.connect(self.plot_delta_omega)
        self.radioButton_dom_180_fold.toggled.connect(self.plot_delta_omega)
        self.radioButton_theta_180_fold.toggled.connect(self.plot_theta)        
        
        self.plot_i.toggled.connect(self.plot_i_Om)
        self.plot_Om.toggled.connect(self.plot_i_Om)
        #self.radioButton_Omega_no_fold.toggled.connect(self.plot_i_Om)
        self.radioButton_Omega_180_fold.toggled.connect(self.plot_i_Om)
       
        self.radioButton_energy.clicked.connect(self.plot_energy)       
        self.radioButton_lx.clicked.connect(self.plot_energy)
        self.radioButton_ly.clicked.connect(self.plot_energy)
        self.radioButton_lz.clicked.connect(self.plot_energy)



        
        self.MMR_combo()
        self.comboBox_MMR_pl_1.activated.connect(self.theta_combo)
        self.comboBox_MMR_pl_1.activated.connect(self.plot_theta)
        self.comboBox_MMR_pl_2.activated.connect(self.theta_combo)
        self.comboBox_MMR_pl_2.activated.connect(self.plot_theta)
        self.comboBox_MMR_theta.activated.connect(self.plot_theta)
        #self.comboBox_MMR_which_pl_1.activated.connect(self.theta_combo)
        self.comboBox_MMR_which_pl_1.activated.connect(self.plot_theta)
        #self.comboBox_MMR_which_pl_2.activated.connect(self.theta_combo)
        self.comboBox_MMR_which_pl_2.activated.connect(self.plot_theta)        
        
        
      #  self.insp_data_size.valueChanged.connect(self.update_inspector)
        
        ############### transit plotting controll ####################      
        self.plot_phase_pholded_tran.stateChanged.connect(self.update_transit_plots)
        self.tra_model_width.valueChanged.connect(self.update_transit_plots)
#        self.tra_model_width.valueChanged.connect(self.update_transit_plots)    
        self.tra_model_z.valueChanged.connect(self.update_transit_plots)
#        self.tra_model_z.valueChanged.connect(self.update_transit_plots)    

        ############### RV GLS plotting controll ####################      
        self.rv_model_width.valueChanged.connect(self.update_RV_plots)
        self.rv_model_width.valueChanged.connect(self.update_extra_plots) 

        ############### RV plotting controll ####################      
        self.rv_model_width.valueChanged.connect(self.update_RV_plots)
        self.rv_model_width.valueChanged.connect(self.update_extra_plots)    
        self.RV_model_z.valueChanged.connect(self.update_RV_plots)
        self.RV_model_z.valueChanged.connect(self.update_extra_plots)    
  
        ############### RV GLS plotting controll ####################      
        self.gls_model_width.valueChanged.connect(self.update_RV_GLS_plots)
        self.gls_o_c_model_width.valueChanged.connect(self.update_RV_o_c_GLS_plots) 
      
        self.N_GLS_peak_to_point.valueChanged.connect(self.update_RV_GLS_plots)
        self.N_GLS_peak_to_point.valueChanged.connect(self.update_RV_o_c_GLS_plots)
        self.avoid_GLS_RV_alias.stateChanged.connect(self.update_RV_GLS_plots)
        self.avoid_GLS_RV_alias.stateChanged.connect(self.update_RV_o_c_GLS_plots)
        
        self.N_TLS_peak_to_point.valueChanged.connect(self.update_tls_plots)
        self.N_TLS_peak_to_point.valueChanged.connect(self.update_tls_o_c_plots)        
        
        
        self.jitter_to_plots.stateChanged.connect(self.update_plots)
        self.split_jitter.stateChanged.connect(self.update_plots)
                
        self.buttonGroup_use_planets.buttonClicked.connect(self.update_veiw)               
#        self.use_Planet1.stateChanged.connect(self.update_veiw)        
        
                
        self.init_correlations_combo()
        self.init_activity_combo()
        self.init_scipy_combo()
        self.comboBox_scipy_minimizer_1.activated.connect(self.check_scipy_min)
        self.comboBox_scipy_minimizer_2.activated.connect(self.check_scipy_min)
        
        
        self.init_gls_norm_combo()
        self.gls_norm_combo.activated.connect(self.update_plots) 

        
        self.setWindowIcon(QtGui.QIcon('./lib/33_striker.png'))
        
        self.radioButton_act_GLS_period.toggled.connect(lambda: self.update_activity_gls_plots(self.comboBox_act_data_gls.currentIndex()))
       
        self.comboBox_act_data_gls.activated.connect(lambda: self.update_activity_gls_plots(self.comboBox_act_data_gls.currentIndex())) 
        self.comboBox_act_data.activated.connect(lambda: self.update_activity_data_plots(self.comboBox_act_data.currentIndex())) 
       
        self.comboBox_corr_1.activated.connect(self.update_correlations_data_plots) 
        self.comboBox_corr_2.activated.connect(self.update_correlations_data_plots) 
        self.plot_corr_err.stateChanged.connect(self.update_correlations_data_plots)
        self.plot_corr_coef.stateChanged.connect(self.update_correlations_data_plots)        

        self.do_RV_GP.stateChanged.connect(self.rv_GP_set_use)
        self.do_tra_GP.stateChanged.connect(self.tra_GP_set_use)


        ############### Cross hair ####################      

        self.gls_cross_hair.stateChanged.connect(self.update_RV_GLS_plots)
        self.gls_o_c_cross_hair.stateChanged.connect(self.update_RV_o_c_GLS_plots)
        self.RV_plot_cross_hair.stateChanged.connect(self.update_RV_plots)
        self.RV_o_c_plot_cross_hair.stateChanged.connect(self.update_RV_plots)
        self.trans_plot_cross_hair.stateChanged.connect(self.update_transit_plots)
        self.trans_o_c_plot_cross_hair.stateChanged.connect(self.update_transit_plots)
        self.tls_cross_hair.stateChanged.connect(self.update_tls_plots)
        self.tls_o_c_cross_hair.stateChanged.connect(self.update_tls_o_c_plots)

        self.extra_plot_cross_hair.stateChanged.connect(self.update_extra_plots)
        self.inpector_plot_cross_hair.stateChanged.connect(lambda: self.plot_data_inspect(self.tree_view_tab.listview))

 

                
        self.color_corr.clicked.connect(self.get_corr_color)
        self.corr_x_label.clicked.connect(self.corr_plot_x_labels)
        self.corr_y_label.clicked.connect(self.corr_plot_y_labels)
        
        self.colors_gls.clicked.connect(self.get_RV_GLS_plot_color)
        self.colors_gls_o_c.clicked.connect(self.get_RV_o_c_GLS_plot_color)


        self.color_delta_om.clicked.connect(self.get_delta_omega_color)
        self.delta_om_x_label.clicked.connect(self.delta_omega_plot_x_labels)
        self.delta_om_y_label.clicked.connect(self.delta_omega_plot_y_labels)

        self.color_theta.clicked.connect(self.get_theta_color)
        self.theta_x_label.clicked.connect(self.theta_plot_x_labels)
        self.theta_y_label.clicked.connect(self.theta_plot_y_labels)

    
        self.tab_timeseries_RV.currentChanged.connect(self.tab_selected)


        self.radioButton_RV_o_c_GLS_period.toggled.connect(self.update_RV_o_c_GLS_plots)
        self.radioButton_RV_GLS_period.toggled.connect(self.update_RV_GLS_plots)
        
        self.mute_boxes()
        self.radioButton_transit_RV.toggled.connect(self.mute_boxes)
        self.radioButton_transit.toggled.connect(self.mute_boxes)
        self.radioButton_RV.toggled.connect(self.mute_boxes)
        
        self.radioButton_ewm.toggled.connect(self.set_hkl)
#        self.radioButton_hkl.toggled.connect(self.set_hkl)
        self.radioButton_KP.toggled.connect(self.set_kp_ma)
 


        self.radioButton_RV_WF_period.toggled.connect(self.update_WF_plots)

        self.calc_TLS.clicked.connect(self.worker_tls)
        self.calc_TLS_o_c.clicked.connect(lambda: self.worker_tls(resid =True))


        self.quit_button.clicked.connect(self.quit)
        self.actionQuit.triggered.connect(self.quit) 

        self.actionopen_RVmod_init_file.triggered.connect(self.showDialog_fortran_input_file)
        self.actionOpen_RVbank_file.triggered.connect(self.showDialog_RVbank_input_file)
        
 
        self.jupiter_push_vars()
        
        self.update_color_picker()
        self.buttonGroup_color_picker.buttonClicked.connect(self.get_color)
        self.update_color_picker_tra()        
        self.buttonGroup_color_picker_tra.buttonClicked.connect(self.get_color_tra)    
       
  
        self.dialog_symbols = show_symbols(self)
        self.buttonGroup_symbol_picker.buttonClicked.connect(self.get_symbol) 
        self.buttonGroup_symbol_picker_tra.buttonClicked.connect(self.get_symbol_tra) 
        
        
        self.buttonGroup_use_RV_GP_kernel.buttonClicked.connect(self.set_RV_GP)   
        self.buttonGroup_use_tra_GP_kernel.buttonClicked.connect(self.set_tra_GP)
        
        self.buttonGroup_use_GP.buttonClicked.connect(self.set_use_GP)
        
        
        ####### LD models #############
 
        #self.buttonGroup_use_ld_1.buttonClicked.connect(self.set_tra_ld)   
        #self.buttonGroup_use_ld_2.buttonClicked.connect(self.set_tra_ld)   
 

       # self.RV_phase_slider.sliderReleased.connect(self.rv_plot_phase_change)       
        self.RV_phase_slider.valueChanged.connect(self.rv_plot_phase_change)       
        
        self.check_settings()       
        self.mute_boxes_dyn()
        self.update_RV_jitter_flag()
        
        self.force_copl_incl.stateChanged.connect(self.set_force_copl_incl)       

        self.threadpool = QtCore.QThreadPool()
        #self.threadpool.setMaxThreadCount(cpu_count())    


        self.update_St_params() 

        ############### Stellar params ####################      
        self.St_mass_input.valueChanged.connect(self.update_St_params)
        self.St_radius_input.valueChanged.connect(self.update_St_params)
        self.St_lumin_input.valueChanged.connect(self.update_St_params)
        self.St_teff_input.valueChanged.connect(self.update_St_params)
        self.St_vsini_input.valueChanged.connect(self.update_St_params)
 
        self.err_St_mass_input.valueChanged.connect(self.update_St_params)
        self.err_St_radius_input.valueChanged.connect(self.update_St_params)
        self.err_St_lumin_input.valueChanged.connect(self.update_St_params)
        self.err_St_teff_input.valueChanged.connect(self.update_St_params)
        self.err_St_vsini_input.valueChanged.connect(self.update_St_params)    
    
        if start_arg_ses == True:
            self.init_fit()         
            self.update_use_from_input_file()   
            self.update_use()
            self.update_gui_params()
            self.update_params()
            self.update_RV_file_buttons()
            self.fit_dispatcher(init=True)  
            self.init_plot_corr()
            self.update_plot_corr()    
    
        print("""Hi there! You are running a demo version of the Exo-Striker (ver. 0.01). 
              
This version is almost full, but there are still some parts of the tool, which are in a 'Work in progress' state. Please, 'git clone' regularly to be up to date with the newest version.
""")

        if sys.version_info[0] == 2:
            print("""
It seems that you started the 'Exo-Striker' with Python 2. Please consider Python 3 for your work with the 'Exo-Striker'.
""") 
            
         
        print("""Here you can get some more information from the tool's workflow, stdout/strerr, and piped results.""")
        #self.use_K1.setStyleSheet("color: red")

def main():
    app = QtWidgets.QApplication(sys.argv)
    app.setStyle('Fusion') #The available styles depend on your platform but are usually 'Fusion', 'Windows', 'WindowsVista' (Windows only) and 'Macintosh' (Mac only). 
    window = Exo_striker()
    screen_resolution = app.desktop().screenGeometry()
    width, height = screen_resolution.width(), screen_resolution.height()
    #print(width, height)
    if height < 920:
        window.setMinimumWidth(width*0.6)
        window.setMinimumHeight(height*0.6)
        window.resize(width*0.8, height*0.8)
    else:
        pass
    window.show()

    sys.exit(app.exec_())

if __name__ == '__main__':
    main() 


