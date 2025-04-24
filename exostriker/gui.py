#!/usr/bin/python3

from exostriker import __version__

import time 
from pathos.multiprocessing import freeze_support
freeze_support()

import numpy as np
#import matplotlib as mpl
#mpl.use('Qt5Agg')

import sys, os, traceback 
from PyQt6 import QtCore, QtGui, QtWidgets, uic
import json

#import collections 
#if sys.version_info.major == 3 and sys.version_info.minor >= 10:
#    from collections.abc import MutableMapping
#else:
#    from collections import MutableMapping

#sys.path.insert(0, './lib') 
es_path = os.path.dirname(os.path.abspath(__file__))
lib_path = os.path.join(es_path, 'lib')

sys.path.insert(0,lib_path)
 
os.chdir(os.path.dirname(os.path.abspath(__file__)))

#sys._excepthook = sys.excepthook
#def exception_hook(exctype, value, traceback):
#    print(exctype, value, traceback)
#    sys._excepthook(exctype, value, traceback)
    #sys.exit(1)
#sys.excepthook = exception_hook


#if QtCore.QT_VERSION >= 0x50501:
#    def excepthook(type_, value, traceback_):
#        traceback.print_exception(type_, value, traceback_)
#        QtCore.qFatal('')
#sys.excepthook = excepthook



if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

import RV_mod as rv

import pyqtgraph as pg
import pyqtgraph.console as pg_console

import word_processor_es as text_editor_es
import calculator as calc
import gls as gls
import mlp as mlp
from worker import Worker #, WorkerSignals
#start_time = time.time()
import gui_groups
import qso

from multiprocessing import cpu_count

#import BKR as bkr
#from doublespinbox import DoubleSpinBox
from Jupyter_emb import ConsoleWidget_embed



from stdout_pipe import MyDialog, DebugDialog
from print_info_window import print_info
from symbols_window import show_symbols
from TESS_pdc_window import pdc

from cornerplot_parameter_window import show_param_boxes



from datafiles_window import datafiles_window
from RVBank_window import RVBank_window

from detrend_window import DetrendWindow
from activity_window import ActivityWindow
from rvs_window import RvsWindow
#from RVBank_window import RVBank_window as DetrendWindow


import terminal
import webbrowser
import ntpath
import pg_hack
import tls_hack

from scipy.signal import argrelextrema
from scipy.stats.stats import pearsonr   
import scipy.stats as stat
from scipy.interpolate import interp1d

#import batman as batman

try:
    import ttvfast as ttvfast
    ttvfast_not_found = False 
except (ImportError, KeyError) as e:
    ttvfast_not_found = True
    pass

try:
    import openai as openai
    openai_not_found = False 
    from exochat import ChatWidget
except (ImportError, KeyError) as e:
    openai_not_found = True
    pass

try:
    from transitleastsquares_ES import transitleastsquares_ES
#    transitleastsquares.helpers.running_median = tls_hack.running_median
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


#try:
#    import cPickle as pickle
#except ModuleNotFoundError:
#    import pickle
import dill
dill._dill._reverse_typemap['ObjectType'] = object
dill.settings['fmode']

#os.system("taskset -p %s" %os.getpid())
os.environ["OPENBLAS_MAIN_FREE"] = "1"
#os.environ["QT_QPA_PLATFORM"] = "offscreen"
os.environ["QT_SCREEN_SCALE_FACTORS"] = "1"
os.environ["QT_SCALE_FACTOR"] = "1"
os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "0"
os.environ["QT_DEVICE_PIXEL_RATIO"] = "1"
os.environ["QT_NO_FT_CACHE"] = "1"
os.environ["PYDEVD_DISABLE_FILE_VALIDATION"] = "1"

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

#os.environ["QT_LOGGING_RULES"] = '*.debug=false'
#os.system("export QT_LOGGING_RULES")

#if hasattr(QtCore.Qt, 'AA_EnableHighDpiScaling'):
#    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)

#QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_Use96Dpi)

#if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
#    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps,True)



#start_time = time.time()


try:
    from es import Ui_MainWindow 
except (ImportError, KeyError) as e:
    qtCreatorFile = "%s/UI/es.ui"%lib_path 
    Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)


#from es import Ui_MainWindow
#print("--- %s seconds ---" % (time.time() - start_time))

pg.setConfigOption('background', '#ffffff')
pg.setConfigOption('foreground', 'k')
pg.setConfigOptions(antialias=True)
#pg.setConfigOptions(useOpenGL=True) 




global fit, colors, ses_list
arguments = len(sys.argv) - 1

if '-debug' in sys.argv:
    debug = True
else:
    debug = False

fit=rv.signal_fit(name='session')


if '-QSO' in sys.argv:
    qso_mode = True
else:
    qso_mode = False
    

if '-last' in sys.argv:
    try:
        file_pi = open("autosave/auto_save.ses", 'rb')
        fit_ses = dill.load(file_pi)
        file_pi.close()   

        fit_sesN = rv.check_for_missing_instances(fit,fit_ses)
 
        fit_sesN.init_pl_arb()
        ses_list = [fit_sesN] 
        
        start_arg_ses = True
        
    except (ImportError, KeyError, AttributeError) as e:
        print("No last session found.")
        ##fit=rv.signal_fit(name='session')
        ses_list = [fit]
        start_arg_ses = False

elif arguments != 0 and sys.argv[1] == '-ses' and os.path.exists(sys.argv[2]):
    try:
        file_pi = open(sys.argv[2], 'rb')
        fit_ses = dill.load(file_pi)
        file_pi.close()   
        #fit = fit_ses 
        fit_sesN = rv.check_for_missing_instances(fit,fit_ses)

        fit_sesN.init_pl_arb()
        ses_list = [fit_sesN] 
        
        start_arg_ses = True  
    except (ImportError, KeyError, AttributeError) as e:
        print("You have entered non-RVmod session. %s cannot be recognaized"%sys.argv[2])
        #fit=rv.signal_fit(name='session')
        ses_list = [fit]
        start_arg_ses = False

elif arguments != 0 and sys.argv[1] == '-mses' and os.path.exists(sys.argv[2]):
    try:
        file_pi = open(sys.argv[2], 'rb')
        fit_ses = dill.load(file_pi)
        file_pi.close()
        ses_list = fit_ses
        fit_sesN = rv.check_for_missing_instances(fit,ses_list[0])
        #fit = ses_list[0]
        fit_sesN.init_pl_arb()
        

        start_arg_ses = True  
    except (ImportError, KeyError, TypeError, AttributeError) as e:
        print("You have entered non-RVmod multi-session. %s cannot be recognaized"%sys.argv[2])
        #fit=rv.signal_fit(name='session')
        ses_list = [fit]
        start_arg_ses = False

#elif  arguments != 0 and sys.argv[1] == '-rv_init' and os.path.exists(sys.argv[2]):
#    try:
        
#        fit=rv.signal_fit(str(sys.argv[2]), 'RVmod session',readinputfile=True)
#        fit.init_pl_arb()
#        ses_list = [fit]
#        start_arg_ses = True
#    except (ImportError, KeyError, TypeError, AttributeError) as e:
#        print("You have entered non-RVmod .init file. %s cannot be recognaized"%sys.argv[2])
#        fit=rv.signal_fit(name='session')
#        ses_list = [fit]
#        start_arg_ses = False

elif arguments != 0 and sys.argv[1] == '-rvbank' and os.path.exists(sys.argv[2]):
    try:
        fit=rv.signal_fit(name='session')
        fit.init_pl_arb()
        fit.add_RVbank_dataset(rv.file_from_path(str(sys.argv[2])), str(sys.argv[2]), split = False)
        fit.fitting(fileinput=False,outputfiles=[1,1,1], minimize_fortran=True, fortran_kill=3, timeout_sec=3, minimize_loglik=True,amoeba_starts=0)
        ses_list = [fit]
        start_arg_ses = True
    except (ImportError, KeyError, TypeError, AttributeError) as e:
        print("You have entered non-RVBank file. %s cannot be recognaized"%sys.argv[2])
        fit=rv.signal_fit(name='session')
        ses_list = [fit]
        start_arg_ses = False

else:
   # fit=rv.signal_fit(name='session')
    ses_list = [fit]
    start_arg_ses = False


colors           = ['#0066ff','#ff0000','#66ff66','#00ffff','#cc33ff','#ff9900','#cccc00','#3399ff','#990033','#339933','#666699']
colors_gls       = ['#0066ff','#ff0000']
colors_delta_om  = ['#0066ff','#666699']
colors_theta     = ['#0066ff','#666699']
colors_per_rat   = ['#0066ff','#666699']
colors_GLS_alias = ['#666699']
colors_MLP_alias = ['#666699']
colors_RV_jitter = ['#000000']
                    
               
symbols = ['o','t','t1','t2','t3','s','p','h','star','+','d','o','t','t1','t2','t3','s','p','h','star','+','d'] 
ttv_interpol_opt = ['linear', 'nearest', 'nearest-up', 'zero', 'slinear', 'quadratic', 'cubic'] #, 'previous', ‘next’.


QtWidgets.QApplication.processEvents()




class Exo_striker(QtWidgets.QMainWindow, Ui_MainWindow):
    



    def update_labels(self):
        global fit

        self.value_stellar_mass.setText("%.4f"%(fit.params.stellar_mass))
        self.value_epoch.setText(str(fit.epoch))
        self.value_rms.setText("%.4f"%(fit.fit_results.rms))
        self.value_wrms.setText("%.4f"%(fit.fit_results.wrms))

        #self.value_wrms.setText("%.4f"%(fit.wrms()))
        
        self.value_chi2.setText("%.4f"%(fit.fit_results.chi2)) 
        self.value_reduced_chi2.setText("%.4f"%(fit.fit_results.reduced_chi2))
        #self.value_loglik.setText("%.4f"%(fit.fit_results.loglik)) 
        self.value_loglik.setText("%.4f"%(fit.loglik)) 
        self.value_BIC.setText("%.2f"%(fit.BIC()))
        self.value_AIC.setText("%.2f"%(fit.AIC()))
        
        #self.value_Ndata.setText("%s"%(len(fit.fit_results.jd))) 
        self.value_Ndata.setText("%s"%(fit.fit_results.Ndata)) 
        
        self.value_DOF.setText("%s"%(int(fit.fit_results.stat.dof)))

        amd = rv.get_AMD_stab(fit)

        if fit.npl >1 and amd == True:
            AMD_LED = './lib/UI/green_led.png'
        elif fit.npl >1 and amd == False:
            AMD_LED = './lib/UI/red_led.png' 
        else:
            AMD_LED = './lib/UI/grey_led.png'

        self.AMD_led.setPixmap(QtGui.QPixmap(AMD_LED))

        if fit.mod_dynamical == True:
            self.radioButton_Dynamical.setChecked(True)
        else:
            self.radioButton_Keplerian.setChecked(True)

        if fit.hkl == True:
            self.radioButton_hkl.setChecked(True)
        else:
            self.radioButton_ewm.setChecked(True)
            
        if fit.rho == True:
            self.radioButton_rho.setChecked(True)
        else:
            self.radioButton_aR.setChecked(True)            
            

    def update_gui_params(self):
        global fit

 
        for i in range(9):
            #if not self.buttonGroup_use_planets.buttons()[i].isChecked():
            if not bool(fit.use_planet[i]):
                continue

            self.param_gui_P[i].setValue(fit.P[i])
            self.param_gui_K[i].setValue(fit.K[i])
            self.param_gui_incl[i].setValue(fit.i[i]%360)
            self.param_gui_Omega[i].setValue(fit.Node[i]%360)

            if fit.hkl == True:

                self.param_gui_e[i].setValue(fit.e_sinw[i])
                self.param_gui_om[i].setValue(fit.e_cosw[i])
                self.param_gui_ma[i].setValue(fit.lamb[i]%360)
            else:
                self.param_gui_e[i].setValue(fit.e[i])
                self.param_gui_om[i].setValue(fit.w[i]%360)
                self.param_gui_ma[i].setValue(fit.M0[i]%360)
 
 
            self.param_gui_wd[i].setValue(fit.omega_dot[i])
 
            self.param_gui_tr[i*3].setValue(fit.t0[i])
            self.param_gui_tr[i*3+1].setValue(fit.pl_rad[i])
            self.param_gui_tr[i*3+2].setValue(fit.pl_a[i])
            

 
        for i in range(20): 
            self.rvs_data_gui[i].setValue(fit.rvoff[i]) 
            self.rvs_data_jitter_gui[i].setValue(fit.jitt[i])

 
            self.tra_data_gui[i].setValue(fit.tra_off[i]) 
            self.tra_data_jitter_gui[i].setValue(fit.tra_jitt[i])
            self.tra_dilution[i][0].setValue(fit.tra_dil[i])

            self.tra_data_lin_trend_gui[i].setValue(fit.tra_lintr[i])
            self.tra_data_quad_trend_gui[i].setValue(fit.tra_quadtr[i])
            
            
 
        for i in range(len(self.gp_rot_params)):
            self.gp_rot_params[i].setValue(fit.GP_rot_params[i])

        for i in range(len(self.gp_sho_params)):
            self.gp_sho_params[i].setValue(fit.GP_sho_params[i])
           
        for i in range(len(self.gp_mat_params)):
            self.gp_mat_params[i].setValue(fit.GP_mat_params[i])   

        for i in range(len(self.gp_drw_params)):
            self.gp_drw_params[i].setValue(fit.GP_drw_params[i])   

        for i in range(len(self.gp_double_sho_params)):
            self.gp_double_sho_params[i].setValue(fit.GP_double_sho_params[i])
         

        for i in range(len(self.tra_gp_rot_params)):
            self.tra_gp_rot_params[i].setValue(fit.tra_GP_rot_params[i])

        for i in range(len(self.tra_gp_sho_params)):
            self.tra_gp_sho_params[i].setValue(fit.tra_GP_sho_params[i])

        for i in range(len(self.tra_gp_mat_params)):
            self.tra_gp_mat_params[i].setValue(fit.tra_GP_mat_params[i])

        for i in range(len(self.tra_gp_drw_params)):
            self.tra_gp_drw_params[i].setValue(fit.tra_GP_drw_params[i])

        for i in range(len(self.tra_gp_double_sho_params)):
            self.tra_gp_double_sho_params[i].setValue(fit.tra_GP_double_sho_params[i])


        for i in range(20):
            self.lin_u[i].setValue(fit.ld_u_lin[i][0])
            self.quad_u1[i].setValue(fit.ld_u_quad[i][0])
            self.quad_u2[i].setValue(fit.ld_u_quad[i][1])
            self.nonlin_u1[i].setValue(fit.ld_u_nonlin[i][0])
            self.nonlin_u2[i].setValue(fit.ld_u_nonlin[i][1])
            self.nonlin_u3[i].setValue(fit.ld_u_nonlin[i][2])
            self.nonlin_u4[i].setValue(fit.ld_u_nonlin[i][3])

 
        
        self.param_gui_ast[0].setValue(fit.ast_alpha[0])
        self.param_gui_ast[1].setValue(fit.ast_delta[0])
        self.param_gui_ast[2].setValue(fit.ast_pi[0])        
        self.param_gui_ast[3].setValue(fit.ast_mu_alpha[0])        
        self.param_gui_ast[4].setValue(fit.ast_mu_delta[0])                
 

        #self.St_mass_input.setValue(fit.params.stellar_mass)  
        #self.St_radius_input.setValue(fit.stellar_radius)  


        self.RV_lin_trend.setValue(fit.rv_lintr)   
        self.RV_quad_trend.setValue(fit.rv_quadtr)   

        self.Epoch.setValue(fit.epoch)
        self.Epoch_ttv.setValue(fit.epoch_ttv)
        self.Epoch_ast.setValue(fit.epoch_ast)
        
        self.update_GUI_St_params()

    def update_params(self):
        global fit
 
        for i in range(9):

            fit.P[i] = self.param_gui_P[i].value()
            fit.K[i] = self.param_gui_K[i].value()
            fit.i[i] = self.param_gui_incl[i].value()
            fit.Node[i] = self.param_gui_Omega[i].value()

            if fit.hkl == True:
                fit.e_sinw[i]  = self.param_gui_e[i].value()
                fit.e_cosw[i]  = self.param_gui_om[i].value()
                fit.lamb[i]    = self.param_gui_ma[i].value()

                fit.e[i]   = np.sqrt(fit.e_sinw[i]**2 + fit.e_cosw[i]**2)
                fit.w[i]   = np.degrees(np.arctan2(np.radians(fit.e_sinw[i]),np.radians(fit.e_cosw[i])))
                fit.M0[i]  = (fit.lamb[i] - fit.w[i])%360.0
            else:

                fit.e[i] = self.param_gui_e[i].value()
                fit.w[i] = self.param_gui_om[i].value()
                fit.M0[i] = self.param_gui_ma[i].value()

                fit.e_sinw[i] = fit.e[i]*np.sin(np.radians(fit.w[i]))
                fit.e_cosw[i] = fit.e[i]*np.cos(np.radians(fit.w[i]))
                fit.lamb[i]   = (fit.w[i] + fit.M0[i])%360.0

 

            fit.omega_dot[i] = self.param_gui_wd[i].value()

            fit.t0[i]     = self.param_gui_tr[i*3].value()
            fit.pl_rad[i] = self.param_gui_tr[i*3+1].value()
            fit.pl_a[i]   = self.param_gui_tr[i*3+2].value()

 
        for i in range(20):
            fit.rvoff[i] = self.rvs_data_gui[i].value()
            fit.jitt[i] = self.rvs_data_jitter_gui[i].value()
 
        for i in range(20):
            fit.tra_off[i]  = self.tra_data_gui[i].value()
            fit.tra_jitt[i] = self.tra_data_jitter_gui[i].value()
            fit.tra_dil[i]  = self.tra_dilution[i][0].value()

            fit.tra_lintr[i]   = self.tra_data_lin_trend_gui[i].value()
            fit.tra_quadtr[i]  = self.tra_data_quad_trend_gui[i].value()
            

        self.read_ast()
        self.read_RV_GP()
        self.read_tra_GP()
        self.read_ld()

        fit.params.stellar_mass = self.St_mass_input.value()
        fit.rv_lintr = self.RV_lin_trend.value()
        fit.rv_quadtr = self.RV_quad_trend.value()

       # fit.stellar_radius = self.St_radius_input.value()

        if self.checkBox_first_RV_epoch.isChecked() and len(fit.fit_results.rv_model.jd) != 0:
            fit.epoch = min(fit.fit_results.rv_model.jd)
        else:
            fit.epoch = self.Epoch.value()

        #### TESTS TTV #####
        if self.checkBox_first_TTV_epoch.isChecked() and len(fit.ttv_data_sets[0]) != 0 and  self.radioButton_ttv.isChecked() ==True:
            fit.epoch_ttv = min(fit.ttv_data_sets[0][1])
            self.Epoch_ttv.setValue(fit.epoch_ttv)
            #self.checkBox_first_TTV_epoch.setEnabled(True)
        elif self.radioButton_ttv.isChecked() == False:
            fit.epoch_ttv = float(fit.epoch)
            self.Epoch_ttv.setValue(float(fit.epoch))
            #self.checkBox_first_TTV_epoch.setEnabled(False)
        else:
            fit.epoch_ttv = self.Epoch_ttv.value()
           # self.checkBox_first_TTV_epoch.setEnabled(True)

        if self.Epoch_ttv_end_plus_1000.isChecked()==True:
            if len(fit.ttv_data_sets[0]) != 0:
                ttv_TS = np.concatenate([fit.ttv_data_sets[x][1] for x in range(10) if len(np.atleast_1d(fit.ttv_data_sets[x])) != 0])
                fit.epoch_ttv_end = max(ttv_TS) + max(fit.P[x] for x in range(9))*2
            else:
                fit.epoch_ttv_end = float(fit.epoch_ttv) + 1000.0 

        elif self.Epoch_ttv_end_plus_1000.isChecked()==False:
            
            fit.epoch_ttv_end = self.Epoch_ttv_end.value()
 
        self.Epoch_ttv_end.setValue(fit.epoch_ttv_end)
 
        #### TESTS AST #####

        if self.checkBox_first_ast_epoch.isChecked() and len(fit.ast_data_sets[0]) != 0 and  self.radioButton_ast.isChecked() ==True:
            fit.epoch_ast = min(fit.ast_data_sets[0][0])
            self.Epoch_ast.setValue(fit.epoch_ast)
 
        elif self.radioButton_ast.isChecked() == False:
            fit.epoch_ast = float(fit.epoch)
            self.Epoch_ast.setValue(float(fit.epoch))
        else:
            fit.epoch_ast = self.Epoch_ast.value()

        if self.Epoch_ast_end_plus_1000.isChecked()==True:
            if len(fit.ast_data_sets[0]) != 0:
                ast_TS = np.concatenate([fit.ast_data_sets[x][0] for x in range(10) if len(np.atleast_1d(fit.ast_data_sets[x])) != 0])
                fit.epoch_ast_end = max(ast_TS) + max(fit.P[x] for x in range(9))*2
            else:
                fit.epoch_ast_end = float(fit.epoch_ast) + 1000.0 

        elif self.Epoch_ast_end_plus_1000.isChecked()==False:
            
            fit.epoch_ast_end = self.Epoch_ast_end.value()
 
        self.Epoch_ast_end.setValue(fit.epoch_ast_end)


    def read_tra_GP(self):
        global fit  

        for i in range(len(self.tra_gp_rot_params)):
            fit.tra_GP_rot_params[i] = self.tra_gp_rot_params[i].value()

        for i in range(len(self.tra_gp_sho_params)):
            fit.tra_GP_sho_params[i] = self.tra_gp_sho_params[i].value()
            
        for i in range(len(self.tra_gp_mat_params)):
            fit.tra_GP_mat_params[i] = self.tra_gp_mat_params[i].value()

        for i in range(len(self.tra_gp_drw_params)):
            fit.tra_GP_drw_params[i] = self.tra_gp_drw_params[i].value()

        for i in range(len(self.tra_gp_double_sho_params)):
            fit.tra_GP_double_sho_params[i] = self.tra_gp_double_sho_params[i].value()

    def read_RV_GP(self):
        global fit

        for i in range(len(self.gp_rot_params)):
            fit.GP_rot_params[i] = self.gp_rot_params[i].value()

        for i in range(len(self.gp_sho_params)):
            fit.GP_sho_params[i] = self.gp_sho_params[i].value()

        for i in range(len(self.gp_mat_params)):
            fit.GP_mat_params[i] = self.gp_mat_params[i].value()

        for i in range(len(self.gp_drw_params)):
            fit.GP_drw_params[i] = self.gp_drw_params[i].value()

        for i in range(len(self.gp_double_sho_params)):
            fit.GP_double_sho_params[i] = self.gp_double_sho_params[i].value()

            
    def read_ld(self):
        global fit

        for i in range(20):
            fit.ld_u_lin[i]    = [self.lin_u[i].value()]
            fit.ld_u_quad[i]   = [self.quad_u1[i].value(),self.quad_u2[i].value()]
            fit.ld_u_nonlin[i] = [self.nonlin_u1[i].value(),self.nonlin_u2[i].value(),
                                      self.nonlin_u3[i].value(),self.nonlin_u4[i].value()]

    def read_ast(self):
        global fit
            
        fit.ast_alpha[0]    = self.param_gui_ast[0].value()
        fit.ast_delta[0]    = self.param_gui_ast[1].value()
        fit.ast_pi[0]       = self.param_gui_ast[2].value()
        fit.ast_mu_alpha[0] = self.param_gui_ast[3].value()
        fit.ast_mu_delta[0] = self.param_gui_ast[4].value()

            

    def set_hkl(self):
        global fit  

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

                self.param_gui_e[i].setRange(0.0,1.0)
               #param_gui[i*3].singleStep(0.001)
                self.param_gui_om[i].setRange(0.0,360.0)
                self.param_gui_ma[i].setRange(0.0,360.0)
                #param_gui[i*3+1].singleStep(0.001)
                self.param_gui_e[i].setValue(fit.e[i])
                self.param_gui_om[i].setValue(fit.w[i]%360.0)
                self.param_gui_ma[i].setValue(fit.M0[i]%360.0)
 

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
                self.param_gui_e[i].setRange(-1.0,1.0)
               # param_gui[i*3].singleStep(0.01)
                self.param_gui_om[i].setRange(-1.0,1.0)
                self.param_gui_ma[i].setRange(0.0,360.0)
                self.param_gui_e[i].setValue(fit.e_sinw[i])
                self.param_gui_om[i].setValue(fit.e_cosw[i])
                self.param_gui_ma[i].setValue(fit.lamb[i])

        self.update_use()
        self.update_params()
        self.update_gui_params() 


            

    def set_aR(self):
        global fit  

        if self.radioButton_aR.isChecked():
            
            fit.rho = False
#            fit.hack_around_rv_params()  
            
 
            self.label_ap.setText("<html><head/><body><p>a<span style='vertical-align:sub;'>pl</span>/R<span style='vertical-align:sub;'>★</span></p></body></html>")
 

        elif self.radioButton_rho.isChecked():

            fit.rho = True
            self.label_ap.setText("<html><head/><body><p>&rho;<span style='vertical-align:sub;'>★</span> [kg/m^3]</p></body></html>")


        self.mute_boxes()

        self.update_use()
        self.update_params()
        self.update_gui_params() 



    def set_tra_reg(self, ind = 0):
        global fit

        #for i in range(10):
        if len(fit.tra_data_sets[ind]) != 0:
       #     continue
       # else:
            fit.tra_data_sets[ind][10] = self.data_tra_reg_group[ind][1].isChecked()

    def update_GUI_tra_reg(self):
        global fit

        for i in range(20):
            if len(fit.tra_data_sets[i]) == 0:
                self.data_tra_reg_group[i][0].setChecked(True)
                #print(self.data_tra_reg_group[i][1].isChecked())
              #  continue
            else:
                #print(fit.tra_data_sets[i][10])
                self.data_tra_reg_group[i][1].setChecked(bool(fit.tra_data_sets[i][10]))



    def set_tra_ld(self):
        global fit

        for i in range(20):
            
            if self.use_uni_ld_models[i].isChecked():
                fit.ld_m[i] = "uniform"
                fit.ld_u[i] = []
            elif self.use_lin_ld_models[i].isChecked():
                fit.ld_m[i] = "linear"
                fit.ld_u[i] = [self.lin_u[i].value()]
            elif self.use_quad_ld_models[i].isChecked():
                fit.ld_m[i] = "quadratic"
                fit.ld_u[i] = [self.quad_u1[i].value(),self.quad_u2[i].value()]
            elif self.use_nonlin_ld_models[i].isChecked():
                fit.ld_m[i] = "nonlinear"
                fit.ld_u[i] = [self.nonlin_u1[i].value(),self.nonlin_u2[i].value(),self.nonlin_u3[i].value(),self.nonlin_u4[i].value()]
            else:
                fit.ld_m[i] = "quadratic"
                fit.ld_u[i] = [self.quad_u1[i].value(),self.quad_u2[i].value()]

        self.set_tra_gr_index()


    def set_tra_gr_index(self):
        global fit

        k =0
 
        for i in range(20):
            if len(fit.tra_data_sets[i]) == 0:
                fit.ld_gr_ind[i] = 0
                continue
            elif fit.ld_gr[i] == i:
                if fit.ld_m[i] == "linear":
                    fit.ld_gr_ind[i] = k
                    k +=1
                elif fit.ld_m[i] == "quadratic":
                    fit.ld_gr_ind[i] = k
                    k +=2
                elif fit.ld_m[i] == "nonlinear":
                    fit.ld_gr_ind[i] = k
                    k +=4
            elif fit.ld_gr[i] != i:
                fit.ld_gr_ind[i] = fit.ld_gr_ind[fit.ld_gr[i]]
                fit.ld_m[i] = fit.ld_m[fit.ld_gr[i]]
                fit.ld_u[i] = fit.ld_u[fit.ld_gr[i]]


##############################
#    k =0
#    
#    z = [0]*10

#    for i in range(10):
#        if len(tr_files[i]) == 0:
#            z[i] = 0
#            continue
#        elif tr_model[2][i] == i:
#            if tr_model[0][i] == "linear":
#                z[i] = k
#                k +=1
#            elif tr_model[0][i] == "quadratic":
#                z[i] = k
#                k +=2
#            elif tr_model[0][i] == "nonlinear":
#                z[i] = k
#                k +=4
#        elif tr_model[2][i] != i:
#            z[i] = z[tr_model[2][i]]
#            tr_model[1][i] = tr_model[1][tr_model[2][i]]
#            tr_model[0][i] = tr_model[0][tr_model[2][i]]
##############################
    


    def set_tra_gr(self):
        global fit
        
        for i in range(20):
            
            if len(fit.tra_data_sets[i]) == 0:
                print("Data set # %s is not present: The request is ignored"%str(i+1))
                self.data_ld_group[i].setValue(i+1)
                fit.ld_gr[i] = i
            else:
                if self.data_ld_group[i].value()-1 > i:
                    print("""
Data set # %s is present, but you cannot tie it to a Data set with a larger index; in this case, # %s. You should do the opposite; (see README): The request is ignored.
"""%(str(i+1), self.data_ld_group[i].value()))
                    self.data_ld_group[i].setValue(i+1)
                    fit.ld_gr[i] = i
                elif len(fit.tra_data_sets[self.data_ld_group[i].value()-1]) == 0:
                    print("Data set # %s is present, but Data set # %s is not: The request is ignored"%(str(i+1), self.data_ld_group[i].value()))
                    self.data_ld_group[i].setValue(i+1)
                    fit.ld_gr[i] = i
                else:
                    fit.ld_gr[i] = self.data_ld_group[i].value() -1
        print("\n")
        
        self.set_tra_gr_index()
        self.tabWidget_helper.setCurrentWidget(self.tab_info)

#        print(fit.ld_gr)

    def update_errors(self):
        global fit
 
        for i in range(9):
            if not self.buttonGroup_use_planets.buttons()[i].isChecked():
                continue

            self.param_errors_gui_P[i].setText("+/- %.3f"%max(np.abs(fit.P_err[i])))
            self.param_errors_gui_K[i].setText("+/- %.3f"%max(np.abs(fit.K_err[i])))

            if fit.hkl == True:
                self.param_errors_gui_e[i].setText("+/- %.3f"%max(np.abs(fit.e_sinw_err[i])))
                self.param_errors_gui_om[i].setText("+/- %.3f"%max(np.abs(fit.e_cosw_err[i])))
                self.param_errors_gui_ma[i].setText("+/- %.3f"%max(np.abs(fit.lamb_err[i])))
            else:

                self.param_errors_gui_e[i].setText("+/- %.3f"%max(np.abs(fit.e_err[i])))
                self.param_errors_gui_om[i].setText("+/- %.3f"%max(np.abs(fit.w_err[i])))
                self.param_errors_gui_ma[i].setText("+/- %.3f"%max(np.abs(fit.M0_err[i])))
 
            self.param_errors_gui_incl[i].setText("+/- %.3f"%max(np.abs(fit.i_err[i])))
            self.param_errors_gui_Omega[i].setText("+/- %.3f"%max(np.abs(fit.Node_err[i])))
 
            self.param_errors_gui_wd[i].setText("+/- %.3f"%max(np.abs(fit.omega_dot_err[i])))

            self.err_t0[i].setText("+/- %.3f"%max(np.abs(fit.t0_err[i])))
            self.err_pl_rad[i].setText("+/- %.3f"%max(np.abs(fit.pl_rad_err[i])))
            self.err_a_sol[i].setText("+/- %.3f"%max(np.abs(fit.pl_a_err[i])))

 

        for i in range(20):

            self.data_errors_gui[i].setText("+/- %.3f"%max(np.abs(fit.rvoff_err[i])))
            self.data_errors_jitter_gui[i].setText("+/- %.3f"%max(np.abs(fit.jitt_err[i])))
 
            self.tra_data_errors_gui[i].setText("+/- %.3f"%max(np.abs(fit.tra_off_err[i])))
            self.tra_data_errors_jitter_gui[i].setText("+/- %.3f"%max(np.abs(fit.tra_jitt_err[i])))

            self.err_tra_data_lin_trend_gui[i].setText("+/- %.3f"%max(np.abs(fit.tra_lintr_err[i])))
            self.err_tra_data_quad_trend_gui[i].setText("+/- %.3f"%max(np.abs(fit.tra_quadtr_err[i])))
 

        self.err_RV_lin_trend.setText("+/- %.8f"%(max(fit.rv_lintr_err)))
        self.err_RV_quad_trend.setText("+/- %.8f"%(max(fit.rv_quadtr_err)))
 

        if fit.gp_kernel == 'RotKernel':
            for i in range(len(self.gp_rot_errors_gui)):
                self.gp_rot_errors_gui[i].setText("+/- %.3f"%max(np.abs(fit.param_errors.GP_params_errors[i])))

        elif fit.gp_kernel == 'SHOKernel':
            for i in range(len(self.gp_sho_errors_gui)):
                self.gp_sho_errors_gui[i].setText("+/- %.3f"%max(np.abs(fit.param_errors.GP_params_errors[i])))
            
        elif fit.gp_kernel == 'Matern32':
            for i in range(len(self.gp_mat_errors_gui)):
                self.gp_mat_errors_gui[i].setText("+/- %.3f"%max(np.abs(fit.param_errors.GP_params_errors[i])))

        elif fit.gp_kernel == 'RealTerm':
            for i in range(len(self.gp_drw_errors_gui)):
                self.gp_drw_errors_gui[i].setText("+/- %.3f"%max(np.abs(fit.param_errors.GP_params_errors[i])))

        elif fit.gp_kernel == 'dSHOKernel':
            for i in range(len(self.gp_double_sho_errors_gui)-1):
                self.gp_double_sho_errors_gui[i].setText("+/- %.3f"%max(np.abs(fit.param_errors.GP_params_errors[i])))

 
        for i in range(20):
            self.err_lin_u[i].setText("+/- %.3f"%max(np.abs(fit.ld_u_lin_err[i][0])))
            self.err_quad_u1[i].setText("+/- %.3f"%max(np.abs(fit.ld_u_quad_err[i][0])))
            self.err_quad_u2[i].setText("+/- %.3f"%max(np.abs(fit.ld_u_quad_err[i][1])))
            self.err_nonlin_u1[i].setText("+/- %.3f"%max(np.abs(fit.ld_u_nonlin_err[i][0])))
            self.err_nonlin_u2[i].setText("+/- %.3f"%max(np.abs(fit.ld_u_nonlin_err[i][1])))
            self.err_nonlin_u3[i].setText("+/- %.3f"%max(np.abs(fit.ld_u_nonlin_err[i][2])))
            self.err_nonlin_u4[i].setText("+/- %.3f"%max(np.abs(fit.ld_u_nonlin_err[i][3])))


        self.param_errors_gui_ast[0].setText("+/- %.3f"%max(np.abs(fit.ast_alpha_err[0])))
        self.param_errors_gui_ast[1].setText("+/- %.3f"%max(np.abs(fit.ast_delta_err[0])))
        self.param_errors_gui_ast[2].setText("+/- %.3f"%max(np.abs(fit.ast_pi_err[0])))
        self.param_errors_gui_ast[3].setText("+/- %.3f"%max(np.abs(fit.ast_mu_alpha_err[0])))
        self.param_errors_gui_ast[4].setText("+/- %.3f"%max(np.abs(fit.ast_mu_delta_err[0])))    
    
    
    
    def update_a_mass(self):
        global fit

        if fit.type_fit["RV"] == True and len(fit.fit_results.a) != 0:
            for i in range(9):
                self.param_a_gui[i].setText("%.5f"%(fit.fit_results.a[i])) 
                self.param_mass_gui[i].setText("%.5f"%(fit.fit_results.mass[i])) 
                self.param_t_peri_gui[i].setText("%.4f"%(fit.t_peri[i]))
 



    def update_use_from_session(self):
        global fit
 
        
        for i in range(9):
            #self.use_param_gui[i].setChecked(bool(fit.use.use_planet_params[i]))
            self.use_param_gui_K[i].setChecked(bool(fit.K_use[i]))
            self.use_param_gui_P[i].setChecked(bool(fit.P_use[i]))
            self.use_param_gui_incl[i].setChecked(bool(fit.i_use[i]))
            self.use_param_gui_Omega[i].setChecked(bool(fit.Node_use[i]))

            if fit.hkl == True:
                self.use_param_gui_e[i].setChecked(bool(fit.e_sinw_use[i]))
                self.use_param_gui_om[i].setChecked(bool(fit.e_cosw_use[i]))
                self.use_param_gui_ma[i].setChecked(bool(fit.lamb_use[i]))
            else:
                self.use_param_gui_e[i].setChecked(bool(fit.e_use[i]))
                self.use_param_gui_om[i].setChecked(bool(fit.w_use[i]))
                self.use_param_gui_ma[i].setChecked(bool(fit.M0_use[i]))


            self.use_param_gui_wd[i].setChecked(bool(fit.omega_dot_use[i]))           
 
        for i in range(fit.npl):         
            self.use_param_gui_tr[i*3].setChecked(bool(fit.t0_use[i]) )        
            self.use_param_gui_tr[i*3+1].setChecked(bool(fit.pl_rad_use[i]) )
            self.use_param_gui_tr[i*3+2].setChecked(bool(fit.pl_a_use [i]) )
 

        for i in range(20): 
            #use_data_gui[i].setChecked(bool(fit.use.use_offsets[i])) # attention, TBF
            self.use_data_jitter_gui[i].setChecked(bool(fit.jitt_use[i]))
            self.use_data_offset_gui[i].setChecked(bool(fit.rvoff_use[i])) 

            
        for i in range(20): 
            self.use_tra_data_jitter_gui[i].setChecked(bool(fit.tra_jitt_use[i]))
            self.use_tra_data_offset_gui[i].setChecked(bool(fit.tra_off_use[i]))
            self.tra_dilution[i][1].setChecked(bool(fit.tra_dil_use[i]))
            self.use_tra_data_lin_trend_gui[i].setChecked(bool(fit.tra_lintr_use[i]))
            self.use_tra_data_quad_trend_gui[i].setChecked(bool(fit.tra_quadtr_use[i]))
              
 
        for i in range(9):  
            if bool(fit.use_planet[i]):
                self.planet_checked_gui[i].setChecked(True)  
            else:
                self.planet_checked_gui[i].setChecked(False)  
        #for i in range(9):  
        #    self.planet_checked_gui[i].setChecked(bool(fit.use_planet[i]))
            
       
    
            
        self.use_RV_lin_trend.setChecked(bool(fit.rv_lintr_use)) 
        self.use_RV_quad_trend.setChecked(bool(fit.rv_quadtr_use)) 
        
        
        for i in range(20):
            self.use_lin_u[i].setChecked(bool(fit.ld_u_lin_use[i][0]))
            self.use_quad_u1[i].setChecked(bool(fit.ld_u_quad_use[i][0]))
            self.use_quad_u2[i].setChecked(bool(fit.ld_u_quad_use[i][1]))
            self.use_nonlin_u1[i].setChecked(bool(fit.ld_u_nonlin_use[i][0]))
            self.use_nonlin_u2[i].setChecked(bool(fit.ld_u_nonlin_use[i][1]))
            self.use_nonlin_u3[i].setChecked(bool(fit.ld_u_nonlin_use[i][2]))
            self.use_nonlin_u4[i].setChecked(bool(fit.ld_u_nonlin_use[i][3]))
 
            self.data_ld_group[i].setValue(int(fit.ld_gr[i])+1)


        for i in range(len(self.use_gp_rot_params)):
            self.use_gp_rot_params[i].setChecked(bool(fit.GP_rot_use[i]))
 
        for i in range(len(self.use_gp_sho_params)):
            self.use_gp_sho_params[i].setChecked(bool(fit.GP_sho_use[i]))

        for i in range(len(self.use_gp_mat_params)):
            self.use_gp_mat_params[i].setChecked(bool(fit.GP_mat_use[i]))

        for i in range(len(self.use_gp_drw_params)):
            self.use_gp_drw_params[i].setChecked(bool(fit.GP_drw_use[i]))            

        for i in range(len(self.use_gp_double_sho_params)):
            self.use_gp_double_sho_params[i].setChecked(bool(fit.GP_double_sho_use[i]))


        for i in range(len(self.use_tra_gp_rot_params)):
            self.use_tra_gp_rot_params[i].setChecked(bool(fit.tra_GP_rot_use[i]))
 
        for i in range(len(self.use_tra_gp_sho_params)):
            self.use_tra_gp_sho_params[i].setChecked(bool(fit.tra_GP_sho_use[i]))

        for i in range(len(self.use_tra_gp_mat_params)):
            self.use_tra_gp_mat_params[i].setChecked(bool(fit.tra_GP_mat_use[i]))

        for i in range(len(self.use_tra_gp_drw_params)):
            self.use_tra_gp_drw_params[i].setChecked(bool(fit.tra_GP_drw_use[i]))

        for i in range(len(self.use_tra_gp_double_sho_params)):
            self.use_tra_gp_double_sho_params[i].setChecked(bool(fit.tra_GP_double_sho_use[i]))


        self.get_use_astrometry()

        self.use_param_gui_ast[0].setChecked(bool(fit.ast_alpha_use[0]))
        self.use_param_gui_ast[1].setChecked(bool(fit.ast_delta_use[0]))
        self.use_param_gui_ast[2].setChecked(bool(fit.ast_pi_use[0]))
        self.use_param_gui_ast[3].setChecked(bool(fit.ast_mu_alpha_use[0]))
        self.use_param_gui_ast[4].setChecked(bool(fit.ast_mu_delta_use[0]))        
 


    def update_mixed_fitting(self):
        global fit
        
        fit.mixed_fit[0][0] =  int(self.use_mix_fitting.isChecked())
        fit.mixed_fit[1] = [int(self.mix_pl_1.isChecked()),int(self.mix_pl_2.isChecked()),int(self.mix_pl_3.isChecked()),
                            int(self.mix_pl_4.isChecked()),int(self.mix_pl_5.isChecked()),int(self.mix_pl_6.isChecked()),
                            int(self.mix_pl_7.isChecked()),int(self.mix_pl_8.isChecked()),int(self.mix_pl_9.isChecked()),
                            ]


    def update_use(self):
        global fit
        

        for i in range(9):  
            fit.use_planet[i] = int(self.planet_checked_gui[i].isChecked()) 
            
        npl_old = fit.npl
        checked = int(np.sum( [self.planet_checked_gui[i].isChecked() for i in range(9)]))

        if npl_old < checked:
            fit.add_planet(index = max(np.nonzero(fit.use_planet)[0], default=0))
        elif npl_old >= checked:
            fit.npl = checked

        self.ttv_pl_combo()
        self.ast_pl_combo()

        for i in range(9):
            #fit.use.use_planet_params[i] = int(self.use_param_gui[i].isChecked())

            fit.P_use[i] = int(self.use_param_gui_P[i].isChecked())
            fit.K_use[i] = int(self.use_param_gui_K[i].isChecked())
            fit.i_use[i] = int(self.use_param_gui_incl[i].isChecked())
            fit.Node_use[i] = int(self.use_param_gui_Omega[i].isChecked()) 

            if fit.hkl == True:
                fit.e_sinw_use[i] = int(self.use_param_gui_e[i].isChecked())
                fit.e_cosw_use[i] = int(self.use_param_gui_om[i].isChecked())
                fit.lamb_use[i] = int(self.use_param_gui_ma[i].isChecked())
            else:
                fit.e_use[i] = int(self.use_param_gui_e[i].isChecked())
                fit.w_use[i] = int(self.use_param_gui_om[i].isChecked())
                fit.M0_use[i] = int(self.use_param_gui_ma[i].isChecked())

            fit.omega_dot_use[i] = int(self.use_param_gui_wd[i].isChecked())

        #for i in range(9):
            fit.t0_use[i] = self.use_param_gui_tr[i*3].isChecked()  
            fit.pl_rad_use[i] = self.use_param_gui_tr[i*3+1].isChecked()  
            fit.pl_a_use[i] = self.use_param_gui_tr[i*3+2].isChecked()  

        for i in range(20): 
            fit.jitt_use[i] = int(self.use_data_jitter_gui[i].isChecked())
            fit.rvoff_use[i] = int(self.use_data_offset_gui[i].isChecked())  

        for i in range(20): 
            fit.tra_jitt_use[i] = int(self.use_tra_data_jitter_gui[i].isChecked())
            fit.tra_off_use[i]  = int(self.use_tra_data_offset_gui[i].isChecked())
            fit.tra_dil_use[i]  = int(self.tra_dilution[i][1].isChecked())

            fit.tra_lintr_use[i]  = int(self.use_tra_data_lin_trend_gui[i].isChecked())
            fit.tra_quadtr_use[i] = int(self.use_tra_data_quad_trend_gui[i].isChecked())
            
       # fit.use.use_linear_trend = int(self.use_RV_lin_trend.isChecked()) 

        fit.rv_lintr_use = int(self.use_RV_lin_trend.isChecked())
        fit.rv_quadtr_use = int(self.use_RV_quad_trend.isChecked())
 
        self.set_use_astrometry() 
        
        self.update_ast_use()
        self.update_RV_GP_use()        
        self.update_tra_GP_use()
        self.update_ld_use()

        self.set_ttv_dataset_to_planet()
        self.set_ast_dataset_to_planet()
        self.set_ast_dataset_to_planet_2()
        
        
    def update_RV_GP_use(self):
        global fit

        for i in range(len(self.use_gp_rot_params)):
            fit.GP_rot_use[i] = int(self.use_gp_rot_params[i].isChecked())

        for i in range(len(self.use_gp_sho_params)):
            fit.GP_sho_use[i] = int(self.use_gp_sho_params[i].isChecked())  

        for i in range(len(self.use_gp_mat_params)):
            fit.GP_mat_use[i] = int(self.use_gp_mat_params[i].isChecked())  

        for i in range(len(self.use_gp_drw_params)):
            fit.GP_drw_use[i] = int(self.use_gp_drw_params[i].isChecked()) 

        for i in range(len(self.use_gp_double_sho_params)):
            fit.GP_double_sho_use[i] = int(self.use_gp_double_sho_params[i].isChecked())  

    def update_tra_GP_use(self):
        global fit

        for i in range(len(self.use_tra_gp_rot_params)):
            fit.tra_GP_rot_use[i] = int(self.use_tra_gp_rot_params[i].isChecked())

        for i in range(len(self.use_tra_gp_sho_params)):
            fit.tra_GP_sho_use[i] = int(self.use_tra_gp_sho_params[i].isChecked())

        for i in range(len(self.use_tra_gp_mat_params)):
            fit.tra_GP_mat_use[i] = int(self.use_tra_gp_mat_params[i].isChecked())

        for i in range(len(self.use_tra_gp_drw_params)):
            fit.tra_GP_drw_use[i] = int(self.use_tra_gp_drw_params[i].isChecked())

        for i in range(len(self.use_tra_gp_double_sho_params)):
            fit.tra_GP_double_sho_use[i] = int(self.use_tra_gp_double_sho_params[i].isChecked())
            
    def update_ld_use(self):
        global fit

        for i in range(20):
            fit.ld_u_lin_use[i]    = [self.use_lin_u[i].isChecked()]
            fit.ld_u_quad_use[i]   = [self.use_quad_u1[i].isChecked(),self.use_quad_u2[i].isChecked()]
            fit.ld_u_nonlin_use[i] = [self.use_nonlin_u1[i].isChecked(),self.use_nonlin_u2[i].isChecked(),
                                      self.use_nonlin_u3[i].isChecked(),self.use_nonlin_u4[i].isChecked()]

    def update_ast_use(self):
        global fit
        
        fit.ast_alpha_use[0]     = self.use_param_gui_ast[0].isChecked()
        fit.ast_delta_use[0]     = self.use_param_gui_ast[1].isChecked()
        fit.ast_pi_use[0]        = self.use_param_gui_ast[2].isChecked()
        fit.ast_mu_alpha_use[0]  = self.use_param_gui_ast[3].isChecked()        
        fit.ast_mu_delta_use[0]  = self.use_param_gui_ast[4].isChecked()          
 

    def update_bounds(self):
        global fit
       

        for i in range(9):
            for z in range(2):
                self.param_bounds_gui[10*i + 0][z].setValue(fit.K_bound[i][z])
                self.param_bounds_gui[10*i + 1][z].setValue(fit.P_bound[i][z])
                self.param_bounds_gui[10*i + 2][z].setValue(fit.e_bound[i][z])
                self.param_bounds_gui[10*i + 3][z].setValue(fit.w_bound[i][z])
                self.param_bounds_gui[10*i + 4][z].setValue(fit.M0_bound[i][z])
                self.param_bounds_gui[10*i + 5][z].setValue(fit.i_bound[i][z])
                self.param_bounds_gui[10*i + 6][z].setValue(fit.Node_bound[i][z])
                self.param_bounds_gui[10*i + 7][z].setValue(fit.t0_bound[i][z])
                self.param_bounds_gui[10*i + 8][z].setValue(fit.pl_rad_bound[i][z])
                self.param_bounds_gui[10*i + 9][z].setValue(fit.pl_a_bound[i][z])


        for i in range(9): 
            for z in range(2):
                self.om_dot_bounds_gui[i][z].setValue(fit.omega_dot_bounds[i][z])
                self.h_bounds_gui[i][z].setValue(fit.e_sinw_bound[i][z])
                self.k_bounds_gui[i][z].setValue(fit.e_cosw_bound[i][z])
                self.lambda_bounds_gui[i][z].setValue(fit.lamb_bound[i][z])

        for i in range(20): 
            for z in range(2):
                self.offset_bounds_gui[i][z].setValue(fit.rvoff_bounds[i][z])
                self.jitter_bounds_gui[i][z].setValue(fit.jitt_bounds[i][z])

        for i in range(20): 
            for z in range(2):
                self.offset_bounds_gui_tra[i][z].setValue(fit.tra_off_bounds[i][z])
                self.jitter_bounds_gui_tra[i][z].setValue(fit.tra_jitt_bounds[i][z])

        self.lin_trend_min.setValue(fit.rv_lintr_bounds[0][0])
        self.lin_trend_max.setValue(fit.rv_lintr_bounds[0][1])    

        self.quad_trend_min.setValue(fit.rv_quadtr_bounds[0][0])
        self.quad_trend_max.setValue(fit.rv_quadtr_bounds[0][1])             
                
        self.update_RV_GP_bounds()
        self.update_tra_GP_bounds()
                          
        self.update_ld_bounds()
        self.update_ast_bounds()       
        
        
        self.update_nr_prior()
        self.update_jeff_prior()



    def update_RV_GP_bounds(self):
        global fit
        
  
        for i in range(3): 
            for z in range(2):
                self.GP_sho_bounds_gui[i][z].setValue(fit.GP_sho_bounds[i][z])
           # if len(fit.GP_sho_bounds[i]) ==2:
           #     fit.GP_sho_bounds[i].append(True)
           #     self.GP_sho_bounds_gui[i][2].setChecked(fit.GP_sho_bounds[i][2])
           # else:
           #     self.GP_sho_bounds_gui[i][2].setChecked(fit.GP_sho_bounds[i][2])
 

        for i in range(4):
            for z in range(2):
                self.GP_rot_bounds_gui[i][z].setValue(fit.GP_rot_bounds[i][z])
            #if len(fit.GP_rot_bounds[i]) ==2:
            #    fit.GP_rot_bounds[i].append(True)
            #    self.GP_rot_bounds_gui[i][2].setChecked(fit.GP_rot_bounds[i][2])
            #else:
            #    self.GP_rot_bounds_gui[i][2].setChecked(fit.GP_rot_bounds[i][2])
            
 
        for i in range(3): 
            for z in range(2):
                self.GP_mat_bounds_gui[i][z].setValue(fit.GP_mat_bounds[i][z])   
            #if len(fit.GP_mat_bounds[i]) ==2:
            #    fit.GP_mat_bounds[i].append(True)
            #    self.GP_mat_bounds_gui[i][2].setChecked(fit.GP_mat_bounds[i][2])
            #else:
            #    self.GP_mat_bounds_gui[i][2].setChecked(fit.GP_mat_bounds[i][2])
        for i in range(2): 
            for z in range(2):
                self.GP_drw_bounds_gui[i][z].setValue(fit.GP_drw_bounds[i][z])   
    
        for i in range(5): 
            for z in range(2):
                self.GP_double_sho_bounds_gui[i][z].setValue(fit.GP_double_sho_bounds[i][z])
               
                
    def update_tra_GP_bounds(self):
        global fit
        
 
        for i in range(3): 
            for z in range(2):
                self.tra_GP_sho_bounds_gui[i][z].setValue(fit.tra_GP_sho_bounds[i][z])
            #if len(fit.tra_GP_sho_bounds[i]) ==2:
            #    fit.tra_GP_sho_bounds[i].append(True)
            #    self.tra_GP_sho_bounds_gui[i][2].setChecked(fit.tra_GP_sho_bounds[i][2])
            #else:
            #    self.tra_GP_sho_bounds_gui[i][2].setChecked(fit.tra_GP_sho_bounds[i][2])

        for i in range(4):
            for z in range(2):
                self.tra_GP_rot_bounds_gui[i][z].setValue(fit.tra_GP_rot_bounds[i][z])
            #if len(fit.tra_GP_rot_bounds[i]) ==2:
            #    fit.tra_GP_rot_bounds[i].append(True)
            #    self.tra_GP_rot_bounds_gui[i][2].setChecked(fit.tra_GP_rot_bounds[i][2])
            #else:
            #    self.tra_GP_rot_bounds_gui[i][2].setChecked(fit.tra_GP_rot_bounds[i][2])
                   
        for i in range(3): 
            for z in range(2):
                self.tra_GP_mat_bounds_gui[i][z].setValue(fit.tra_GP_mat_bounds[i][z])  
            #if len(fit.tra_GP_mat_bounds[i]) ==2:
            #    fit.tra_GP_mat_bounds[i].append(True)
            #    self.tra_GP_mat_bounds_gui[i][2].setChecked(fit.tra_GP_mat_bounds[i][2])
            #else:
            #    self.tra_GP_mat_bounds_gui[i][2].setChecked(fit.tra_GP_mat_bounds[i][2])   
        for i in range(2): 
            for z in range(2):
                self.tra_GP_drw_bounds_gui[i][z].setValue(fit.tra_GP_drw_bounds[i][z])  
             
        for i in range(5): 
            for z in range(2):
                self.tra_GP_double_sho_bounds_gui[i][z].setValue(fit.tra_GP_double_sho_bounds[i][z])
                

    def update_nr_prior(self):
        global fit


        for i in range(9):
            for z in range(2):

                self.param_nr_priors_gui[10*i + 0][z].setValue(fit.K_norm_pr[i][z])
                self.param_nr_priors_gui[10*i + 1][z].setValue(fit.P_norm_pr[i][z])
                self.param_nr_priors_gui[10*i + 2][z].setValue(fit.e_norm_pr[i][z])
                self.param_nr_priors_gui[10*i + 3][z].setValue(fit.w_norm_pr[i][z])
                self.param_nr_priors_gui[10*i + 4][z].setValue(fit.M0_norm_pr[i][z])
                self.param_nr_priors_gui[10*i + 5][z].setValue(fit.i_norm_pr[i][z])
                self.param_nr_priors_gui[10*i + 6][z].setValue(fit.Node_norm_pr[i][z])
                self.param_nr_priors_gui[10*i + 7][z].setValue(fit.t0_norm_pr[i][z])
                self.param_nr_priors_gui[10*i + 8][z].setValue(fit.pl_rad_norm_pr[i][z])
                self.param_nr_priors_gui[10*i + 9][z].setValue(fit.pl_a_norm_pr[i][z])
                       
            self.param_nr_priors_gui[10*i + 0][2].setChecked(int(fit.K_norm_pr[i][2]))
            self.param_nr_priors_gui[10*i + 1][2].setChecked(int(fit.P_norm_pr[i][2]))
            self.param_nr_priors_gui[10*i + 2][2].setChecked(int(fit.e_norm_pr[i][2]))
            self.param_nr_priors_gui[10*i + 3][2].setChecked(int(fit.w_norm_pr[i][2]))
            self.param_nr_priors_gui[10*i + 4][2].setChecked(int(fit.M0_norm_pr[i][2]))
            self.param_nr_priors_gui[10*i + 5][2].setChecked(int(fit.i_norm_pr[i][2]))
            self.param_nr_priors_gui[10*i + 6][2].setChecked(int(fit.Node_norm_pr[i][2]))
            self.param_nr_priors_gui[10*i + 7][2].setChecked(int(fit.t0_norm_pr[i][2]))
            self.param_nr_priors_gui[10*i + 8][2].setChecked(int(fit.pl_rad_norm_pr[i][2]))
            self.param_nr_priors_gui[10*i + 9][2].setChecked(int(fit.pl_a_norm_pr[i][2]))
 

        for i in range(9): 
            for z in range(2):
                self.om_dot_norm_pr_gui[i][z].setValue(fit.omega_dot_norm_pr[i][z])
                self.h_norm_pr_gui[i][z].setValue(fit.e_sinw_norm_pr[i][z])
                self.k_norm_pr_gui[i][z].setValue(fit.e_cosw_norm_pr[i][z])
                self.lambda_norm_pr_gui[i][z].setValue(fit.lamb_norm_pr[i][z])

            self.om_dot_norm_pr_gui[i][2].setChecked(int(fit.omega_dot_norm_pr[i][2]))
            self.h_norm_pr_gui[i][2].setChecked(int(fit.e_sinw_norm_pr[i][2]))
            self.k_norm_pr_gui[i][2].setChecked(int(fit.e_cosw_norm_pr[i][2]))
            self.lambda_norm_pr_gui[i][2].setChecked(int(fit.lamb_norm_pr[i][2]))

        self.lin_trend_mean.setValue(fit.rv_lintr_norm_pr[0][0])
        self.lin_trend_sigma.setValue(fit.rv_lintr_norm_pr[0][1])
        self.use_lin_tr_nr_pr.setChecked(int(fit.rv_lintr_norm_pr[0][2]))
 
        self.quad_trend_mean.setValue(fit.rv_quadtr_norm_pr[0][0])
        self.quad_trend_sigma.setValue(fit.rv_quadtr_norm_pr[0][1])
        self.use_quad_tr_nr_pr.setChecked(int(fit.rv_quadtr_norm_pr[0][2]))
 
 
        for i in range(20): 
            for z in range(2):

                 self.offset_nr_priors_gui[i][z].setValue(fit.rvoff_norm_pr[i][z])
                 self.jitter_nr_priors_gui[i][z].setValue(fit.jitt_norm_pr[i][z])
                 
                 
            self.offset_nr_priors_gui[i][2].setChecked(int(fit.rvoff_norm_pr[i][2]))
            self.jitter_nr_priors_gui[i][2].setChecked(int(fit.jitt_norm_pr[i][2]))

        for i in range(20): 
            for z in range(2):

                 self.offset_nr_priors_gui_tra[i][z].setValue(fit.tra_off_norm_pr[i][z])
                 self.jitter_nr_priors_gui_tra[i][z].setValue(fit.tra_jitt_norm_pr[i][z])
                 self.tra_lin_trend_nr_priors_gui[i][z].setValue(fit.tra_lintr_norm_pr[i][z])
                 self.tra_quad_trend_nr_priors_gui[i][z].setValue(fit.tra_quadtr_norm_pr[i][z])  
 
            self.offset_nr_priors_gui_tra[i][2].setChecked(int(fit.tra_off_norm_pr[i][2]))
            self.jitter_nr_priors_gui_tra[i][2].setChecked(int(fit.tra_jitt_norm_pr[i][2]))
            self.tra_lin_trend_nr_priors_gui[i][2].setChecked(int(fit.tra_lintr_norm_pr[i][2]))
            self.tra_quad_trend_nr_priors_gui[i][2].setChecked(int(fit.tra_quadtr_norm_pr[i][2]))


        self.update_RV_GP_priors_nr()
        self.update_tra_GP_priors_nr()
        self.update_ld_priors_nr()
        self.update_ast_priors_nr()
        
    def update_RV_GP_priors_nr(self):
        global fit
 
        for i in range(4): 
            for z in range(2):
                self.GP_rot_nr_priors_gui[i][z].setValue(fit.GP_rot_norm_pr[i][z])
            self.GP_rot_nr_priors_gui[i][2].setChecked(int(fit.GP_rot_norm_pr[i][2]))
    
        for i in range(3): 
            for z in range(2):
                self.GP_sho_nr_priors_gui[i][z].setValue(fit.GP_sho_norm_pr[i][z])
            self.GP_sho_nr_priors_gui[i][2].setChecked(int(fit.GP_sho_norm_pr[i][2]))
 
        for i in range(3): 
            for z in range(2):
                self.GP_mat_nr_priors_gui[i][z].setValue(fit.GP_mat_norm_pr[i][z])
            self.GP_mat_nr_priors_gui[i][2].setChecked(int(fit.GP_mat_norm_pr[i][2]))   

        for i in range(2): 
            for z in range(2):
                self.GP_drw_nr_priors_gui[i][z].setValue(fit.GP_drw_norm_pr[i][z])
            self.GP_drw_nr_priors_gui[i][2].setChecked(int(fit.GP_drw_norm_pr[i][2]))  


        for i in range(5): 
            for z in range(2):
                self.GP_double_sho_nr_priors_gui[i][z].setValue(fit.GP_double_sho_norm_pr[i][z])
            self.GP_double_sho_nr_priors_gui[i][2].setChecked(int(fit.GP_double_sho_norm_pr[i][2])     )    
            
    def update_tra_GP_priors_nr(self):
        global fit
 
        for i in range(4): 
            for z in range(2):
                self.tra_GP_rot_nr_priors_gui[i][z].setValue(fit.tra_GP_rot_norm_pr[i][z])
            self.tra_GP_rot_nr_priors_gui[i][2].setChecked(int(fit.tra_GP_rot_norm_pr[i][2]))
    
        for i in range(3): 
            for z in range(2):
                self.tra_GP_sho_nr_priors_gui[i][z].setValue(fit.tra_GP_sho_norm_pr[i][z])
            self.tra_GP_sho_nr_priors_gui[i][2].setChecked(int(fit.tra_GP_sho_norm_pr[i][2]))
 
        for i in range(3): 
            for z in range(2):
                self.tra_GP_mat_nr_priors_gui[i][z].setValue(fit.tra_GP_mat_norm_pr[i][z])
            self.tra_GP_mat_nr_priors_gui[i][2].setChecked(int(fit.tra_GP_mat_norm_pr[i][2]))       

        for i in range(2): 
            for z in range(2):
                self.tra_GP_drw_nr_priors_gui[i][z].setValue(fit.tra_GP_drw_norm_pr[i][z])
            self.tra_GP_drw_nr_priors_gui[i][2].setChecked(int(fit.tra_GP_drw_norm_pr[i][2]))  
            
        for i in range(5): 
            for z in range(2):
                self.tra_GP_double_sho_nr_priors_gui[i][z].setValue(fit.tra_GP_double_sho_norm_pr[i][z])
            self.tra_GP_double_sho_nr_priors_gui[i][2].setChecked(int(fit.tra_GP_double_sho_norm_pr[i][2])  )          






    def update_jeff_prior(self):
        global fit

        for i in range(9):
            for z in range(2):

                self.param_jeff_priors_gui[10*i + 0][z].setValue(fit.K_jeff_pr[i][z])
                self.param_jeff_priors_gui[10*i + 1][z].setValue(fit.P_jeff_pr[i][z])
                self.param_jeff_priors_gui[10*i + 2][z].setValue(fit.e_jeff_pr[i][z])
                self.param_jeff_priors_gui[10*i + 3][z].setValue(fit.w_jeff_pr[i][z])
                self.param_jeff_priors_gui[10*i + 4][z].setValue(fit.M0_jeff_pr[i][z])
                self.param_jeff_priors_gui[10*i + 5][z].setValue(fit.i_jeff_pr[i][z])
                self.param_jeff_priors_gui[10*i + 6][z].setValue(fit.Node_jeff_pr[i][z])
                self.param_jeff_priors_gui[10*i + 7][z].setValue(fit.t0_jeff_pr[i][z])
                self.param_jeff_priors_gui[10*i + 8][z].setValue(fit.pl_rad_jeff_pr[i][z])
                self.param_jeff_priors_gui[10*i + 9][z].setValue(fit.pl_a_jeff_pr[i][z])
                
            self.param_jeff_priors_gui[10*i + 0][2].setChecked(int(fit.K_jeff_pr[i][2]))
            self.param_jeff_priors_gui[10*i + 1][2].setChecked(int(fit.P_jeff_pr[i][2]))
            self.param_jeff_priors_gui[10*i + 2][2].setChecked(int(fit.e_jeff_pr[i][2]))
            self.param_jeff_priors_gui[10*i + 3][2].setChecked(int(fit.w_jeff_pr[i][2]))
            self.param_jeff_priors_gui[10*i + 4][2].setChecked(int(fit.M0_jeff_pr[i][2]))
            self.param_jeff_priors_gui[10*i + 5][2].setChecked(int(fit.i_jeff_pr[i][2]))
            self.param_jeff_priors_gui[10*i + 6][2].setChecked(int(fit.Node_jeff_pr[i][2]))
            self.param_jeff_priors_gui[10*i + 7][2].setChecked(int(fit.t0_jeff_pr[i][2]))
            self.param_jeff_priors_gui[10*i + 8][2].setChecked(int(fit.pl_rad_jeff_pr[i][2]))
            self.param_jeff_priors_gui[10*i + 9][2].setChecked(int(fit.pl_a_jeff_pr[i][2]))
 

        for i in range(9): 
            for z in range(2):
                self.om_dot_jeff_pr_gui[i][z].setValue(fit.omega_dot_jeff_pr[i][z])
                self.h_jeff_pr_gui[i][z].setValue(fit.e_sinw_jeff_pr[i][z])
                self.k_jeff_pr_gui[i][z].setValue(fit.e_cosw_jeff_pr[i][z])
                self.lambda_jeff_pr_gui[i][z].setValue(fit.lamb_jeff_pr[i][z])

            self.om_dot_jeff_pr_gui[i][2].setChecked(int(fit.omega_dot_jeff_pr[i][2]))
            self.h_jeff_pr_gui[i][2].setChecked(int(fit.e_sinw_jeff_pr[i][2]))
            self.k_jeff_pr_gui[i][2].setChecked(int(fit.e_cosw_jeff_pr[i][2]))
            self.lambda_jeff_pr_gui[i][2].setChecked(int(fit.lamb_jeff_pr[i][2]))

        self.lin_trend_jeff_alpha.setValue(fit.rv_lintr_jeff_pr[0][0])
        self.lin_trend_jeff_beta.setValue(fit.rv_lintr_jeff_pr[0][1])
        self.use_lin_tr_jeff_pr.setChecked(int(fit.rv_lintr_jeff_pr[0][2]))

        self.quad_trend_jeff_alpha.setValue(fit.rv_quadtr_jeff_pr[0][0])
        self.quad_trend_jeff_beta.setValue(fit.rv_quadtr_jeff_pr[0][1])
        self.use_quad_tr_jeff_pr.setChecked(int(fit.rv_quadtr_jeff_pr[0][2]) )
 

        for i in range(20): 
            for z in range(2):

                 self.offset_jeff_priors_gui[i][z].setValue(fit.rvoff_jeff_pr[i][z])
                 self.jitter_jeff_priors_gui[i][z].setValue(fit.jitt_jeff_pr[i][z])      
                 
            self.offset_jeff_priors_gui[i][2].setChecked(int(fit.rvoff_jeff_pr[i][2]))
            self.jitter_jeff_priors_gui[i][2].setChecked(int(fit.jitt_jeff_pr[i][2]))


        for i in range(20): 
            for z in range(2):
 
                 self.offset_jeff_priors_gui_tra[i][z].setValue(fit.tra_off_jeff_pr[i][z])
                 self.jitter_jeff_priors_gui_tra[i][z].setValue(fit.tra_jitt_jeff_pr[i][z])
                 self.tra_lin_trend_jeff_priors_gui[i][z].setValue(fit.tra_lintr_jeff_pr[i][z])
                 self.tra_quad_trend_jeff_priors_gui[i][z].setValue(fit.tra_quadtr_jeff_pr[i][z])                 
                 
 
            self.offset_jeff_priors_gui_tra[i][2].setChecked(int(fit.tra_off_jeff_pr[i][2]))
            self.jitter_jeff_priors_gui_tra[i][2].setChecked(int(fit.tra_jitt_jeff_pr[i][2]))
            self.tra_lin_trend_jeff_priors_gui[i][2].setChecked(int(fit.tra_lintr_jeff_pr[i][2]))
            self.tra_quad_trend_jeff_priors_gui[i][2].setChecked(int(fit.tra_quadtr_jeff_pr[i][2]))

            

        self.update_RV_GP_priors_jeff()
        self.update_tra_GP_priors_jeff()
        self.update_ld_priors_jeff()
        self.update_ast_priors_jeff()


    def update_RV_GP_priors_jeff(self):
        global fit
 
        for i in range(4): 
            for z in range(2):
                self.GP_rot_jeff_priors_gui[i][z].setValue(fit.GP_rot_jeff_pr[i][z])
            self.GP_rot_jeff_priors_gui[i][2].setChecked(int(fit.GP_rot_jeff_pr[i][2]))
    
        for i in range(3): 
            for z in range(2):
                self.GP_sho_jeff_priors_gui[i][z].setValue(fit.GP_sho_jeff_pr[i][z])
            self.GP_sho_jeff_priors_gui[i][2].setChecked(int(fit.GP_sho_jeff_pr[i][2]))
 
        for i in range(3): 
            for z in range(2):
                self.GP_mat_jeff_priors_gui[i][z].setValue(fit.GP_mat_jeff_pr[i][z])
            self.GP_mat_jeff_priors_gui[i][2].setChecked(int(fit.GP_mat_jeff_pr[i][2]))            

        for i in range(2): 
            for z in range(2):
                self.GP_drw_jeff_priors_gui[i][z].setValue(fit.GP_drw_jeff_pr[i][z])
            self.GP_drw_jeff_priors_gui[i][2].setChecked(int(fit.GP_drw_jeff_pr[i][2]))  

        for i in range(5): 
            for z in range(2):
                self.GP_double_sho_jeff_priors_gui[i][z].setValue(fit.GP_double_sho_jeff_pr[i][z])
            self.GP_double_sho_jeff_priors_gui[i][2].setChecked(int(fit.GP_double_sho_jeff_pr[i][2]))
 
            
    def update_tra_GP_priors_jeff(self):
        global fit
 
        for i in range(4): 
            for z in range(2):
                self.tra_GP_rot_jeff_priors_gui[i][z].setValue(fit.tra_GP_rot_jeff_pr[i][z])
            self.tra_GP_rot_jeff_priors_gui[i][2].setChecked(int(fit.tra_GP_rot_jeff_pr[i][2]))
    
        for i in range(3): 
            for z in range(2):
                self.tra_GP_sho_jeff_priors_gui[i][z].setValue(fit.tra_GP_sho_jeff_pr[i][z])
            self.tra_GP_sho_jeff_priors_gui[i][2].setChecked(int(fit.tra_GP_sho_jeff_pr[i][2]))
 
        for i in range(3): 
            for z in range(2):
                self.tra_GP_mat_jeff_priors_gui[i][z].setValue(fit.tra_GP_mat_jeff_pr[i][z])
            self.tra_GP_mat_jeff_priors_gui[i][2].setChecked(int(fit.tra_GP_mat_jeff_pr[i][2]))

        for i in range(2): 
            for z in range(2):
                self.tra_GP_drw_jeff_priors_gui[i][z].setValue(fit.tra_GP_drw_jeff_pr[i][z])
            self.tra_GP_drw_jeff_priors_gui[i][2].setChecked(int(fit.tra_GP_drw_jeff_pr[i][2]))

        for i in range(5): 
            for z in range(2):
                self.tra_GP_double_sho_jeff_priors_gui[i][z].setValue(fit.tra_GP_double_sho_jeff_pr[i][z])
            self.tra_GP_double_sho_jeff_priors_gui[i][2].setChecked(int(fit.tra_GP_double_sho_jeff_pr[i][2]))

    def check_bounds(self):
        global fit

        for i in range(9):
            for z in range(2):

                fit.K_bound[i][z] = self.param_bounds_gui[10*i + 0][z].value()
                fit.P_bound[i][z] = self.param_bounds_gui[10*i + 1][z].value()
                fit.e_bound[i][z] = self.param_bounds_gui[10*i + 2][z].value()
                fit.w_bound[i][z] = self.param_bounds_gui[10*i + 3][z].value()
                fit.M0_bound[i][z] = self.param_bounds_gui[10*i + 4][z].value()
                fit.i_bound[i][z] = self.param_bounds_gui[10*i + 5][z].value()
                fit.Node_bound[i][z] =self.param_bounds_gui[10*i + 6][z].value()
                fit.t0_bound[i][z]  =  self.param_bounds_gui[10*i + 7][z].value()
                fit.pl_rad_bound[i][z]  =   self.param_bounds_gui[10*i + 8][z].value()
                fit.pl_a_bound[i][z]   =   self.param_bounds_gui[10*i + 9][z].value()
                fit.K_bound[i][z] = self.param_bounds_gui[10*i + 0][z].value()

        for i in range(20): 
            for z in range(2):
                fit.rvoff_bounds[i][z] = self.offset_bounds_gui[i][z].value()
                fit.jitt_bounds[i][z]  = self.jitter_bounds_gui[i][z].value()

        for i in range(9): 
            for z in range(2):
                fit.omega_dot_bounds[i][z] = self.om_dot_bounds_gui[i][z].value()
                fit.e_sinw_bound[i][z] = self.h_bounds_gui[i][z].value()
                fit.e_cosw_bound[i][z] = self.k_bounds_gui[i][z].value()
                fit.lamb_bound[i][z] = self.lambda_bounds_gui[i][z].value()
 

        fit.rv_lintr_bounds[0]  = [self.lin_trend_min.value(),self.lin_trend_max.value()]
        fit.rv_quadtr_bounds[0] = [self.quad_trend_min.value(),self.quad_trend_max.value()]


        for i in range(20): 
            for z in range(2):
                fit.tra_off_bounds[i][z]    = self.offset_bounds_gui_tra[i][z].value()
                fit.tra_jitt_bounds[i][z]   = self.jitter_bounds_gui_tra[i][z].value()

                fit.tra_lintr_bounds[i][z]  = self.tra_lintr_bounds_gui[i][z].value()
                fit.tra_quadtr_bounds[i][z] = self.tra_quadtr_bounds_gui[i][z].value()
 


        self.check_RV_GP_bounds()
        self.check_tra_GP_bounds()
        self.check_ld_bounds()
        self.check_ast_bounds()


    def check_RV_GP_bounds(self):
        global fit
        
 
        for i in range(4): 
            for z in range(2):
                fit.GP_rot_bounds[i][z] = self.GP_rot_bounds_gui[i][z].value()
            #fit.GP_rot_bounds[i][2] = self.GP_rot_bounds_gui[i][2].isChecked()
            
        for i in range(3): 
            for z in range(2):
                fit.GP_sho_bounds[i][z]  = self.GP_sho_bounds_gui[i][z].value()
           # fit.GP_sho_bounds[i][2]  = self.GP_sho_bounds_gui[i][2].isChecked()
 
        for i in range(3): 
            for z in range(2):
                fit.GP_mat_bounds[i][z]  = self.GP_mat_bounds_gui[i][z].value()           
            #fit.GP_mat_bounds[i][2]  = self.GP_mat_bounds_gui[i][2].isChecked()   
            
        for i in range(2): 
            for z in range(2):
                fit.GP_drw_bounds[i][z]  = self.GP_drw_bounds_gui[i][z].value()   
                
        for i in range(5): 
            for z in range(2):
                fit.GP_double_sho_bounds[i][z]  = self.GP_double_sho_bounds_gui[i][z].value()
           # fit.GP_sho_bounds[i][2]  = self.GP_sho_bounds_gui[i][2].isChecked()            

    def check_tra_GP_bounds(self):
        global fit

        for i in range(4): 
            for z in range(2):
                fit.tra_GP_rot_bounds[i][z] = self.tra_GP_rot_bounds_gui[i][z].value()
            #fit.tra_GP_rot_bounds[i][2] = self.tra_GP_rot_bounds_gui[i][2].isChecked()
            
        for i in range(3): 
            for z in range(2):
                fit.tra_GP_sho_bounds[i][z]  = self.tra_GP_sho_bounds_gui[i][z].value()
           # fit.tra_GP_sho_bounds[i][2]  = self.tra_GP_sho_bounds_gui[i][2].isChecked()
 
        for i in range(3): 
            for z in range(2):
                fit.tra_GP_mat_bounds[i][z]  = self.tra_GP_mat_bounds_gui[i][z].value()           
            #fit.tra_GP_mat_bounds[i][2]  = self.tra_GP_mat_bounds_gui[i][2].isChecked() 
            
        for i in range(2): 
            for z in range(2):
                fit.tra_GP_drw_bounds[i][z]  = self.tra_GP_drw_bounds_gui[i][z].value()           
             
        for i in range(5): 
            for z in range(2):
                fit.tra_GP_double_sho_bounds[i][z]  = self.tra_GP_double_sho_bounds_gui[i][z].value()
           # fit.tra_GP_sho_bounds[i][2]  = self.tra_GP_sho_bounds_gui[i][2].isChecked()
            

    def check_ld_bounds(self):
        global fit

        for i in range(20): 
            for z in range(2):    
                fit.ld_u_lin_bound[i][0][z] = self.ld_u1_bounds_gui[i][z].value()

                fit.ld_u_quad_bound[i][0][z] = self.ld_u1_bounds_gui[i][z].value()
                fit.ld_u_quad_bound[i][1][z] = self.ld_u2_bounds_gui[i][z].value()
                
                fit.ld_u_nonlin_bound[i][0][z] = self.ld_u1_bounds_gui[i][z].value()
                fit.ld_u_nonlin_bound[i][1][z] = self.ld_u2_bounds_gui[i][z].value()
                fit.ld_u_nonlin_bound[i][2][z] = self.ld_u3_bounds_gui[i][z].value()
                fit.ld_u_nonlin_bound[i][3][z] = self.ld_u4_bounds_gui[i][z].value()


    def update_ld_bounds(self):
        global fit

        for i in range(20): 
            for z in range(2):    
                self.ld_u1_bounds_gui[i][z].setValue(fit.ld_u_lin_bound[i][0][z])

                self.ld_u1_bounds_gui[i][z].setValue(fit.ld_u_quad_bound[i][0][z])
                self.ld_u2_bounds_gui[i][z].setValue(fit.ld_u_quad_bound[i][1][z])
                
                self.ld_u1_bounds_gui[i][z].setValue(fit.ld_u_nonlin_bound[i][0][z])
                self.ld_u2_bounds_gui[i][z].setValue(fit.ld_u_nonlin_bound[i][1][z])
                self.ld_u3_bounds_gui[i][z].setValue(fit.ld_u_nonlin_bound[i][2][z])
                self.ld_u4_bounds_gui[i][z].setValue(fit.ld_u_nonlin_bound[i][3][z])


    def check_ld_priors_nr(self):
        global fit

        for i in range(20): 
            for z in range(2):    
                fit.ld_u_lin_norm_pr[i][0][z] = self.ld_u1_norm_pr_gui[i][z].value()

                fit.ld_u_quad_norm_pr[i][0][z] = self.ld_u1_norm_pr_gui[i][z].value()
                fit.ld_u_quad_norm_pr[i][1][z] = self.ld_u2_norm_pr_gui[i][z].value()
                
                fit.ld_u_nonlin_norm_pr[i][0][z] = self.ld_u1_norm_pr_gui[i][z].value()
                fit.ld_u_nonlin_norm_pr[i][1][z] = self.ld_u2_norm_pr_gui[i][z].value()
                fit.ld_u_nonlin_norm_pr[i][2][z] = self.ld_u3_norm_pr_gui[i][z].value()
                fit.ld_u_nonlin_norm_pr[i][3][z] = self.ld_u4_norm_pr_gui[i][z].value()

            fit.ld_u_lin_norm_pr[i][0][2] = self.ld_u1_norm_pr_gui[i][2].isChecked()

            fit.ld_u_quad_norm_pr[i][0][2] = self.ld_u1_norm_pr_gui[i][2].isChecked()
            fit.ld_u_quad_norm_pr[i][1][2] = self.ld_u2_norm_pr_gui[i][2].isChecked()
            
            fit.ld_u_nonlin_norm_pr[i][0][2] = self.ld_u1_norm_pr_gui[i][2].isChecked()
            fit.ld_u_nonlin_norm_pr[i][1][2] = self.ld_u2_norm_pr_gui[i][2].isChecked()
            fit.ld_u_nonlin_norm_pr[i][2][2] = self.ld_u3_norm_pr_gui[i][2].isChecked()
            fit.ld_u_nonlin_norm_pr[i][3][2] = self.ld_u4_norm_pr_gui[i][2].isChecked()


    def update_ld_priors_nr(self):
        global fit



        for i in range(20): 
            for z in range(2):    
                self.ld_u1_norm_pr_gui[i][z].setValue(fit.ld_u_lin_norm_pr[i][0][z])

                self.ld_u1_norm_pr_gui[i][z].setValue(fit.ld_u_quad_norm_pr[i][0][z])
                self.ld_u2_norm_pr_gui[i][z].setValue(fit.ld_u_quad_norm_pr[i][1][z])

                        
                        
                self.ld_u1_norm_pr_gui[i][z].setValue(fit.ld_u_nonlin_norm_pr[i][0][z])
                self.ld_u2_norm_pr_gui[i][z].setValue(fit.ld_u_nonlin_norm_pr[i][1][z])
                self.ld_u3_norm_pr_gui[i][z].setValue(fit.ld_u_nonlin_norm_pr[i][2][z])
                self.ld_u4_norm_pr_gui[i][z].setValue(fit.ld_u_nonlin_norm_pr[i][3][z])

            self.ld_u1_norm_pr_gui[i][2].setChecked(int(fit.ld_u_lin_norm_pr[i][0][2]))

            self.ld_u1_norm_pr_gui[i][2].setChecked(int(fit.ld_u_quad_norm_pr[i][0][2]))
            self.ld_u2_norm_pr_gui[i][2].setChecked(int(fit.ld_u_quad_norm_pr[i][1][2]))
            
            self.ld_u1_norm_pr_gui[i][2].setChecked(int(fit.ld_u_nonlin_norm_pr[i][0][2]))
            self.ld_u2_norm_pr_gui[i][2].setChecked(int(fit.ld_u_nonlin_norm_pr[i][1][2]))
            self.ld_u3_norm_pr_gui[i][2].setChecked(int(fit.ld_u_nonlin_norm_pr[i][2][2]))
            self.ld_u4_norm_pr_gui[i][2].setChecked(int(fit.ld_u_nonlin_norm_pr[i][3][2]))



    def check_ld_priors_jeff(self):
        global fit

        for i in range(20): 
            for z in range(2):    
                fit.ld_u_lin_jeff_pr[i][0][z] = self.ld_u1_jeff_pr_gui[i][z].value()

                fit.ld_u_quad_jeff_pr[i][0][z] = self.ld_u1_jeff_pr_gui[i][z].value()
                fit.ld_u_quad_jeff_pr[i][1][z] = self.ld_u2_jeff_pr_gui[i][z].value()
                
                fit.ld_u_nonlin_jeff_pr[i][0][z] = self.ld_u1_jeff_pr_gui[i][z].value()
                fit.ld_u_nonlin_jeff_pr[i][1][z] = self.ld_u2_jeff_pr_gui[i][z].value()
                fit.ld_u_nonlin_jeff_pr[i][2][z] = self.ld_u3_jeff_pr_gui[i][z].value()
                fit.ld_u_nonlin_jeff_pr[i][3][z] = self.ld_u4_jeff_pr_gui[i][z].value()

            fit.ld_u_lin_jeff_pr[i][0][2] = self.ld_u1_jeff_pr_gui[i][2].isChecked()

            fit.ld_u_quad_jeff_pr[i][0][2] = self.ld_u1_jeff_pr_gui[i][2].isChecked()
            fit.ld_u_quad_jeff_pr[i][1][2] = self.ld_u2_jeff_pr_gui[i][2].isChecked()
            
            fit.ld_u_nonlin_jeff_pr[i][0][2] = self.ld_u1_jeff_pr_gui[i][2].isChecked()
            fit.ld_u_nonlin_jeff_pr[i][1][2] = self.ld_u2_jeff_pr_gui[i][2].isChecked()
            fit.ld_u_nonlin_jeff_pr[i][2][2] = self.ld_u3_jeff_pr_gui[i][2].isChecked()
            fit.ld_u_nonlin_jeff_pr[i][3][2] = self.ld_u4_jeff_pr_gui[i][2].isChecked()



   
    def update_ld_priors_jeff(self):
        global fit

        for i in range(20): 
            for z in range(2):    
                self.ld_u1_jeff_pr_gui[i][z].setValue(fit.ld_u_lin_jeff_pr[i][0][z])

                self.ld_u1_jeff_pr_gui[i][z].setValue(fit.ld_u_quad_jeff_pr[i][0][z])
                self.ld_u2_jeff_pr_gui[i][z].setValue(fit.ld_u_quad_jeff_pr[i][1][z])
                
                self.ld_u1_jeff_pr_gui[i][z].setValue(fit.ld_u_nonlin_jeff_pr[i][0][z])
                self.ld_u2_jeff_pr_gui[i][z].setValue(fit.ld_u_nonlin_jeff_pr[i][1][z])
                self.ld_u3_jeff_pr_gui[i][z].setValue(fit.ld_u_nonlin_jeff_pr[i][2][z])
                self.ld_u4_jeff_pr_gui[i][z].setValue(fit.ld_u_nonlin_jeff_pr[i][3][z])

            self.ld_u1_jeff_pr_gui[i][2].setChecked(int(fit.ld_u_lin_jeff_pr[i][0][2]))

            self.ld_u1_jeff_pr_gui[i][2].setChecked(int(fit.ld_u_quad_jeff_pr[i][0][2]))
            self.ld_u2_jeff_pr_gui[i][2].setChecked(int(fit.ld_u_quad_jeff_pr[i][1][2]))
            
            self.ld_u1_jeff_pr_gui[i][2].setChecked(int(fit.ld_u_nonlin_jeff_pr[i][0][2]))
            self.ld_u2_jeff_pr_gui[i][2].setChecked(int(fit.ld_u_nonlin_jeff_pr[i][1][2]))
            self.ld_u3_jeff_pr_gui[i][2].setChecked(int(fit.ld_u_nonlin_jeff_pr[i][2][2]))
            self.ld_u4_jeff_pr_gui[i][2].setChecked(int(fit.ld_u_nonlin_jeff_pr[i][3][2]))
                        
    ########### Astrometry ##########
    def update_ast_bounds(self):
        global fit

        for z in range(2):    
            self.ast_bounds_gui[0][z].setValue(fit.ast_alpha_bound[0][z])
            self.ast_bounds_gui[1][z].setValue(fit.ast_delta_bound[0][z])
            self.ast_bounds_gui[2][z].setValue(fit.ast_pi_bound[0][z])
            self.ast_bounds_gui[3][z].setValue(fit.ast_mu_alpha_bound[0][z])
            self.ast_bounds_gui[4][z].setValue(fit.ast_mu_delta_bound[0][z])                        


    def update_ast_priors_nr(self):
        global fit

        for z in range(2):    
            self.ast_norm_pr_gui[0][z].setValue(fit.ast_alpha_norm_pr[0][z])
            self.ast_norm_pr_gui[1][z].setValue(fit.ast_delta_norm_pr[0][z])
            self.ast_norm_pr_gui[2][z].setValue(fit.ast_pi_norm_pr[0][z])
            self.ast_norm_pr_gui[3][z].setValue(fit.ast_mu_alpha_norm_pr[0][z])
            self.ast_norm_pr_gui[4][z].setValue(fit.ast_mu_delta_norm_pr[0][z])  
            
        self.ast_norm_pr_gui[0][2].setChecked(int(fit.ast_alpha_norm_pr[0][2]))
        self.ast_norm_pr_gui[1][2].setChecked(int(fit.ast_delta_norm_pr[0][2]))
        self.ast_norm_pr_gui[2][2].setChecked(int(fit.ast_pi_norm_pr[0][2]))
        self.ast_norm_pr_gui[3][2].setChecked(int(fit.ast_mu_alpha_norm_pr[0][2]))
        self.ast_norm_pr_gui[4][2].setChecked(int(fit.ast_mu_delta_norm_pr[0][2]))   
        
    def update_ast_priors_jeff(self):
        global fit

        for z in range(2):    
            self.ast_jeff_pr_gui[0][z].setValue(fit.ast_alpha_jeff_pr[0][z])
            self.ast_jeff_pr_gui[1][z].setValue(fit.ast_delta_jeff_pr[0][z])
            self.ast_jeff_pr_gui[2][z].setValue(fit.ast_pi_jeff_pr[0][z])
            self.ast_jeff_pr_gui[3][z].setValue(fit.ast_mu_alpha_jeff_pr[0][z])
            self.ast_jeff_pr_gui[4][z].setValue(fit.ast_mu_delta_jeff_pr[0][z])  
            
        self.ast_jeff_pr_gui[0][2].setChecked(int(fit.ast_alpha_jeff_pr[0][2]))
        self.ast_jeff_pr_gui[1][2].setChecked(int(fit.ast_delta_jeff_pr[0][2]))
        self.ast_jeff_pr_gui[2][2].setChecked(int(fit.ast_pi_jeff_pr[0][2]))
        self.ast_jeff_pr_gui[3][2].setChecked(int(fit.ast_mu_alpha_jeff_pr[0][2]))
        self.ast_jeff_pr_gui[4][2].setChecked(int(fit.ast_mu_delta_jeff_pr[0][2]))          

    def check_ast_bounds(self):
        global fit
 
        for z in range(2):    
            fit.ast_alpha_bound[0][z] = self.ast_bounds_gui[0][z].value()
            fit.ast_delta_bound[0][z] = self.ast_bounds_gui[1][z].value() 
            fit.ast_pi_bound[0][z] = self.ast_bounds_gui[2][z].value()
            fit.ast_mu_alpha_bound[0][z] = self.ast_bounds_gui[3][z].value()
            fit.ast_mu_delta_bound[0][z] = self.ast_bounds_gui[4][z].value()
            
    def check_ast_priors_nr(self):
        global fit
 
        for z in range(2):    
            fit.ast_alpha_norm_pr[0][z]    = self.ast_norm_pr_gui[0][z].value()
            fit.ast_delta_norm_pr[0][z]    = self.ast_norm_pr_gui[1][z].value() 
            fit.ast_pi_norm_pr[0][z]       = self.ast_norm_pr_gui[2][z].value()
            fit.ast_mu_alpha_norm_pr[0][z] = self.ast_norm_pr_gui[3][z].value()
            fit.ast_mu_delta_norm_pr[0][z] = self.ast_norm_pr_gui[4][z].value()
                      
        fit.ast_alpha_norm_pr[0][2]    = self.ast_norm_pr_gui[0][2].isChecked()
        fit.ast_delta_norm_pr[0][2]    = self.ast_norm_pr_gui[1][2].isChecked() 
        fit.ast_pi_norm_pr[0][2]       = self.ast_norm_pr_gui[2][2].isChecked()
        fit.ast_mu_alpha_norm_pr[0][2] = self.ast_norm_pr_gui[3][2].isChecked()
        fit.ast_mu_delta_norm_pr[0][2] = self.ast_norm_pr_gui[4][2].isChecked()
                            
    def check_ast_priors_jeff(self):
        global fit
 
        for z in range(2):    
            fit.ast_alpha_jeff_pr[0][z]    = self.ast_jeff_pr_gui[0][z].value()
            fit.ast_delta_jeff_pr[0][z]    = self.ast_jeff_pr_gui[1][z].value() 
            fit.ast_pi_jeff_pr[0][z]       = self.ast_jeff_pr_gui[2][z].value()
            fit.ast_mu_alpha_jeff_pr[0][z] = self.ast_jeff_pr_gui[3][z].value()
            fit.ast_mu_delta_jeff_pr[0][z] = self.ast_jeff_pr_gui[4][z].value()
                      
        fit.ast_alpha_jeff_pr[0][2]    = self.ast_jeff_pr_gui[0][2].isChecked()
        fit.ast_delta_jeff_pr[0][2]    = self.ast_jeff_pr_gui[1][2].isChecked() 
        fit.ast_pi_jeff_pr[0][2]       = self.ast_jeff_pr_gui[2][2].isChecked()
        fit.ast_mu_alpha_jeff_pr[0][2] = self.ast_jeff_pr_gui[3][2].isChecked()
        fit.ast_mu_delta_jeff_pr[0][2] = self.ast_jeff_pr_gui[4][2].isChecked()
                                      
                                    
    ################################# 

    def check_priors_nr(self):
        global fit

        for i in range(9):
            for z in range(2):

                fit.K_norm_pr[i][z] = self.param_nr_priors_gui[10*i + 0][z].value()
                fit.P_norm_pr[i][z] = self.param_nr_priors_gui[10*i + 1][z].value()
                fit.e_norm_pr[i][z] = self.param_nr_priors_gui[10*i + 2][z].value()
                fit.w_norm_pr[i][z] = self.param_nr_priors_gui[10*i + 3][z].value()
                fit.M0_norm_pr[i][z] = self.param_nr_priors_gui[10*i + 4][z].value()
                fit.i_norm_pr[i][z] = self.param_nr_priors_gui[10*i + 5][z].value()
                fit.Node_norm_pr[i][z] =self.param_nr_priors_gui[10*i + 6][z].value()
                fit.t0_norm_pr[i][z]  =  self.param_nr_priors_gui[10*i + 7][z].value()
                fit.pl_rad_norm_pr[i][z]  =   self.param_nr_priors_gui[10*i + 8][z].value()
                fit.pl_a_norm_pr[i][z]   =   self.param_nr_priors_gui[10*i + 9][z].value()

            fit.K_norm_pr[i][2] = self.param_nr_priors_gui[10*i + 0][2].isChecked()
            fit.P_norm_pr[i][2] = self.param_nr_priors_gui[10*i + 1][2].isChecked()
            fit.e_norm_pr[i][2] = self.param_nr_priors_gui[10*i + 2][2].isChecked()
            fit.w_norm_pr[i][2] = self.param_nr_priors_gui[10*i + 3][2].isChecked()
            fit.M0_norm_pr[i][2] = self.param_nr_priors_gui[10*i + 4][2].isChecked()
            fit.i_norm_pr[i][2] = self.param_nr_priors_gui[10*i + 5][2].isChecked()
            fit.Node_norm_pr[i][2] =self.param_nr_priors_gui[10*i + 6][2].isChecked()
            fit.t0_norm_pr[i][2]  =  self.param_nr_priors_gui[10*i + 7][2].isChecked()
            fit.pl_rad_norm_pr[i][2]  =   self.param_nr_priors_gui[10*i + 8][2].isChecked()
            fit.pl_a_norm_pr[i][2]   =   self.param_nr_priors_gui[10*i + 9][2].isChecked()
 

        for i in range(20): 
            for z in range(2):

                fit.rvoff_norm_pr[i][z] = self.offset_nr_priors_gui[i][z].value()
                fit.jitt_norm_pr[i][z]  = self.jitter_nr_priors_gui[i][z].value()
            fit.rvoff_norm_pr[i][2] = self.offset_nr_priors_gui[i][2].isChecked()
            fit.jitt_norm_pr[i][2]  = self.jitter_nr_priors_gui[i][2].isChecked()

        for i in range(9): 
            for z in range(2):
                fit.omega_dot_norm_pr[i][z] = self.om_dot_norm_pr_gui[i][z].value()
                fit.e_sinw_norm_pr[i][z] = self.h_norm_pr_gui[i][z].value()
                fit.e_cosw_norm_pr[i][z] = self.k_norm_pr_gui[i][z].value()
                fit.lamb_norm_pr[i][z] = self.lambda_norm_pr_gui[i][z].value()

            fit.omega_dot_norm_pr[i][2] = self.om_dot_norm_pr_gui[i][2].isChecked()
            fit.e_sinw_norm_pr[i][2] = self.h_norm_pr_gui[i][2].isChecked()
            fit.e_cosw_norm_pr[i][2] = self.k_norm_pr_gui[i][2].isChecked()
            fit.lamb_norm_pr[i][2] = self.lambda_norm_pr_gui[i][2].isChecked()


        fit.rv_lintr_norm_pr[0]  = [self.lin_trend_mean.value(),self.lin_trend_sigma.value(),self.use_lin_tr_nr_pr.isChecked()]
        fit.rv_quadtr_norm_pr[0]  = [self.quad_trend_mean.value(),self.quad_trend_sigma.value(),self.use_quad_tr_nr_pr.isChecked()]


        for i in range(20): 
            for z in range(2):

                fit.tra_off_norm_pr[i][z]    = self.offset_nr_priors_gui_tra[i][z].value()
                fit.tra_jitt_norm_pr[i][z]   = self.jitter_nr_priors_gui_tra[i][z].value()
                fit.tra_lintr_norm_pr[i][z]  = self.tra_lin_trend_nr_priors_gui[i][z].value()
                fit.tra_quadtr_norm_pr[i][z] = self.tra_quad_trend_nr_priors_gui[i][z].value()
                


            fit.tra_off_norm_pr[i][2]  = self.offset_nr_priors_gui_tra[i][2].isChecked()
            fit.tra_jitt_norm_pr[i][2] = self.jitter_nr_priors_gui_tra[i][2].isChecked()
            fit.tra_lintr_norm_pr[i][2]  = self.tra_lin_trend_nr_priors_gui[i][2].isChecked()
            fit.tra_quadtr_norm_pr[i][2] = self.tra_quad_trend_nr_priors_gui[i][2].isChecked()
            
            
        self.check_RV_GP_priors_nr()
        self.check_tra_GP_priors_nr()
        self.check_ld_priors_nr()
        self.check_ast_priors_nr()

    def check_RV_GP_priors_nr(self):
        global fit


        for i in range(4): 
            for z in range(2):
                fit.GP_rot_norm_pr[i][z] = self.GP_rot_nr_priors_gui[i][z].value()
            fit.GP_rot_norm_pr[i][2] = self.GP_rot_nr_priors_gui[i][2].isChecked()

        for i in range(3): 
            for z in range(2):
                fit.GP_sho_norm_pr[i][z] = self.GP_sho_nr_priors_gui[i][z].value()   
            fit.GP_sho_norm_pr[i][2] = self.GP_sho_nr_priors_gui[i][2].isChecked()   
 
        for i in range(3): 
            for z in range(2):
                fit.GP_mat_norm_pr[i][z] = self.GP_mat_nr_priors_gui[i][z].value()               
            fit.GP_mat_norm_pr[i][2] = self.GP_mat_nr_priors_gui[i][2].isChecked()               

        for i in range(2): 
            for z in range(2):
                fit.GP_drw_norm_pr[i][z] = self.GP_drw_nr_priors_gui[i][z].value()               
            fit.GP_drw_norm_pr[i][2] = self.GP_drw_nr_priors_gui[i][2].isChecked()   

        for i in range(5): 
            for z in range(2):
                fit.GP_double_sho_norm_pr[i][z] = self.GP_double_sho_nr_priors_gui[i][z].value()   
            fit.GP_double_sho_norm_pr[i][2] = self.GP_double_sho_nr_priors_gui[i][2].isChecked()   

    def check_tra_GP_priors_nr(self):
        global fit
 
        for i in range(4): 
            for z in range(2):
                fit.tra_GP_rot_norm_pr[i][z] = self.tra_GP_rot_nr_priors_gui[i][z].value()
            fit.tra_GP_rot_norm_pr[i][2] = self.tra_GP_rot_nr_priors_gui[i][2].isChecked()
 
        for i in range(3): 
            for z in range(2):
                fit.tra_GP_sho_norm_pr[i][z] = self.tra_GP_sho_nr_priors_gui[i][z].value()   
            fit.tra_GP_sho_norm_pr[i][2] = self.tra_GP_sho_nr_priors_gui[i][2].isChecked()   
 
        for i in range(3): 
            for z in range(2):
                fit.tra_GP_mat_norm_pr[i][z] = self.tra_GP_mat_nr_priors_gui[i][z].value()               
            fit.tra_GP_mat_norm_pr[i][2] = self.tra_GP_mat_nr_priors_gui[i][2].isChecked()               

        for i in range(2): 
            for z in range(2):
                fit.tra_GP_drw_norm_pr[i][z] = self.tra_GP_drw_nr_priors_gui[i][z].value()               
            fit.tra_GP_drw_norm_pr[i][2] = self.tra_GP_drw_nr_priors_gui[i][2].isChecked()  

        for i in range(5): 
            for z in range(2):
                fit.tra_GP_double_sho_norm_pr[i][z] = self.tra_GP_double_sho_nr_priors_gui[i][z].value()   
            fit.tra_GP_double_sho_norm_pr[i][2] = self.tra_GP_double_sho_nr_priors_gui[i][2].isChecked()       
    
    
            
    def check_priors_jeff(self):
        global fit

        for i in range(9):
            for z in range(2):
                
                fit.K_jeff_pr[i][z] = self.param_jeff_priors_gui[10*i + 0][z].value()
                fit.P_jeff_pr[i][z] = self.param_jeff_priors_gui[10*i + 1][z].value()
                fit.e_jeff_pr[i][z] = self.param_jeff_priors_gui[10*i + 2][z].value()
                fit.w_jeff_pr[i][z] = self.param_jeff_priors_gui[10*i + 3][z].value()
                fit.M0_jeff_pr[i][z] = self.param_jeff_priors_gui[10*i + 4][z].value()
                fit.i_jeff_pr[i][z] = self.param_jeff_priors_gui[10*i + 5][z].value()
                fit.Node_jeff_pr[i][z] =self.param_jeff_priors_gui[10*i + 6][z].value()
                fit.t0_jeff_pr[i][z]  =  self.param_jeff_priors_gui[10*i + 7][z].value()
                fit.pl_rad_jeff_pr[i][z]  =   self.param_jeff_priors_gui[10*i + 8][z].value()
                fit.pl_a_jeff_pr[i][z]   =   self.param_jeff_priors_gui[10*i + 9][z].value()
            

            fit.K_jeff_pr[i][2] = self.param_jeff_priors_gui[10*i + 0][2].isChecked()
            fit.P_jeff_pr[i][2] = self.param_jeff_priors_gui[10*i + 1][2].isChecked()
            fit.e_jeff_pr[i][2] = self.param_jeff_priors_gui[10*i + 2][2].isChecked()
            fit.w_jeff_pr[i][2] = self.param_jeff_priors_gui[10*i + 3][2].isChecked()
            fit.M0_jeff_pr[i][2] = self.param_jeff_priors_gui[10*i + 4][2].isChecked()
            fit.i_jeff_pr[i][2] = self.param_jeff_priors_gui[10*i + 5][2].isChecked()
            fit.Node_jeff_pr[i][2] =self.param_jeff_priors_gui[10*i + 6][2].isChecked()
            fit.t0_jeff_pr[i][2]  =  self.param_jeff_priors_gui[10*i + 7][2].isChecked()
            fit.pl_rad_jeff_pr[i][2]  =   self.param_jeff_priors_gui[10*i + 8][2].isChecked()
            fit.pl_a_jeff_pr[i][2]   =   self.param_jeff_priors_gui[10*i + 9][2].isChecked()


        for i in range(20): 
            for z in range(2):

                fit.rvoff_jeff_pr[i][z] = self.offset_jeff_priors_gui[i][z].value()
                fit.jitt_jeff_pr[i][z]  = self.jitter_jeff_priors_gui[i][z].value()
            fit.rvoff_jeff_pr[i][2] = self.offset_jeff_priors_gui[i][2].isChecked()
            fit.jitt_jeff_pr[i][2]  = self.jitter_jeff_priors_gui[i][2].isChecked()


        for i in range(9): 
            for z in range(2):
                fit.omega_dot_jeff_pr[i][z] = self.om_dot_jeff_pr_gui[i][z].value()
                fit.e_sinw_jeff_pr[i][z] = self.h_jeff_pr_gui[i][z].value()
                fit.e_cosw_jeff_pr[i][z] = self.k_jeff_pr_gui[i][z].value()
                fit.lamb_jeff_pr[i][z] = self.lambda_jeff_pr_gui[i][z].value()

            fit.omega_dot_jeff_pr[i][2] = self.om_dot_jeff_pr_gui[i][2].isChecked()
            fit.e_sinw_jeff_pr[i][2] = self.h_jeff_pr_gui[i][2].isChecked()
            fit.e_cosw_jeff_pr[i][2] = self.k_jeff_pr_gui[i][2].isChecked()
            fit.lamb_jeff_pr[i][2] = self.lambda_jeff_pr_gui[i][2].isChecked()

    
 
        fit.rv_lintr_jeff_pr[0]   = [self.lin_trend_jeff_alpha.value(),self.lin_trend_jeff_beta.value(),self.use_lin_tr_jeff_pr.isChecked()]
        fit.rv_quadtr_jeff_pr[0]  = [self.quad_trend_jeff_alpha.value(),self.quad_trend_jeff_beta.value(),self.use_quad_tr_jeff_pr.isChecked()]


        for i in range(20): 
            for z in range(2):
                fit.tra_off_jeff_pr[i][z]    = self.offset_jeff_priors_gui_tra[i][z].value()
                fit.tra_jitt_jeff_pr[i][z]   = self.jitter_jeff_priors_gui_tra[i][z].value() 
                fit.tra_lintr_jeff_pr[i][z]  = self.tra_lin_trend_jeff_priors_gui[i][z].value()
                fit.tra_quadtr_jeff_pr[i][z] = self.tra_quad_trend_jeff_priors_gui[i][z].value()                
                
            fit.tra_off_jeff_pr[i][2]    = self.offset_jeff_priors_gui_tra[i][2].isChecked()
            fit.tra_jitt_jeff_pr[i][2]   = self.jitter_jeff_priors_gui_tra[i][2].isChecked()
            fit.tra_lintr_jeff_pr[i][2]  = self.tra_lin_trend_jeff_priors_gui[i][2].isChecked()
            fit.tra_quadtr_jeff_pr[i][2] = self.tra_quad_trend_jeff_priors_gui[i][2].isChecked()            
            
        self.check_RV_GP_priors_jeff()
        self.check_tra_GP_priors_jeff()
        self.check_ld_priors_jeff()    
        self.check_ast_priors_jeff()        

    def check_RV_GP_priors_jeff(self):
        global fit


        for i in range(4): 
            for z in range(2):
                fit.GP_rot_jeff_pr[i][z] = self.GP_rot_jeff_priors_gui[i][z].value()
            fit.GP_rot_jeff_pr[i][2] = self.GP_rot_jeff_priors_gui[i][2].isChecked()

        for i in range(3): 
            for z in range(2):
                fit.GP_sho_jeff_pr[i][z] = self.GP_sho_jeff_priors_gui[i][z].value()   
            fit.GP_sho_jeff_pr[i][2] = self.GP_sho_jeff_priors_gui[i][2].isChecked()   
 
        for i in range(3): 
            for z in range(2):
                fit.GP_mat_jeff_pr[i][z] = self.GP_mat_jeff_priors_gui[i][z].value()               
            fit.GP_mat_jeff_pr[i][2] = self.GP_mat_jeff_priors_gui[i][2].isChecked()               

        for i in range(2): 
            for z in range(2):
                fit.GP_drw_jeff_pr[i][z] = self.GP_drw_jeff_priors_gui[i][z].value()               
            fit.GP_drw_jeff_pr[i][2] = self.GP_drw_jeff_priors_gui[i][2].isChecked()  

        for i in range(5): 
            for z in range(2):
                fit.GP_double_sho_jeff_pr[i][z] = self.GP_double_sho_jeff_priors_gui[i][z].value()   
            fit.GP_double_sho_jeff_pr[i][2] = self.GP_double_sho_jeff_priors_gui[i][2].isChecked()   
 

    def check_tra_GP_priors_jeff(self):
        global fit
 
        for i in range(4): 
            for z in range(2):
                fit.tra_GP_rot_jeff_pr[i][z] = self.tra_GP_rot_jeff_priors_gui[i][z].value()
            fit.tra_GP_rot_jeff_pr[i][2] = self.tra_GP_rot_jeff_priors_gui[i][2].isChecked()
 
        for i in range(3): 
            for z in range(2):
                fit.tra_GP_sho_jeff_pr[i][z] = self.tra_GP_sho_jeff_priors_gui[i][z].value()   
            fit.tra_GP_sho_jeff_pr[i][2] = self.tra_GP_sho_jeff_priors_gui[i][2].isChecked()   
 
        for i in range(3): 
            for z in range(2):
                fit.tra_GP_mat_jeff_pr[i][z] = self.tra_GP_mat_jeff_priors_gui[i][z].value()               
            fit.tra_GP_mat_jeff_pr[i][2] = self.tra_GP_mat_jeff_priors_gui[i][2].isChecked()  
             
        for i in range(2): 
            for z in range(2):
                fit.tra_GP_drw_jeff_pr[i][z] = self.tra_GP_drw_jeff_priors_gui[i][z].value()               
            fit.tra_GP_drw_jeff_pr[i][2] = self.tra_GP_drw_jeff_priors_gui[i][2].isChecked()   

        for i in range(5): 
            for z in range(2):
                fit.tra_GP_double_sho_jeff_pr[i][z] = self.tra_GP_double_sho_jeff_priors_gui[i][z].value()   
            fit.tra_GP_double_sho_jeff_pr[i][2] = self.tra_GP_double_sho_jeff_priors_gui[i][2].isChecked()               
 

    def check_arb_pl(self):
        global fit


        fit.arb_st_mass = self.arb_st_mass.value()

        j = 0
        for i in range(9):
            fit.pl_arb_use[i] = self.arb_param_gui_use[i].isChecked()
            
            if fit.pl_arb_use[i] == True:
                j += 1
            
            fit.e_arb[i]    = self.arb_param_gui[7*i + 2].value()    
            fit.w_arb[i]    = self.arb_param_gui[7*i + 3].value()    
            fit.M0_arb[i]   = self.arb_param_gui[7*i + 4].value()    
            fit.i_arb[i]    = self.arb_param_gui[7*i + 5].value()    
            fit.Node_arb[i] = self.arb_param_gui[7*i + 6].value()    
 
            if self.radioButton_KP.isChecked():
                fit.K_arb[i]    = self.arb_param_gui[7*i + 0].value()    
                fit.P_arb[i]    = self.arb_param_gui[7*i + 1].value()                 
               # mass_,a_ = rv.mass_a_from_Kepler_fit([fit.K_arb[i],  fit.P_arb[i], fit.e_arb[i],  fit.w_arb[i], fit.M0_arb[i]],1,fit.arb_st_mass)
                mass_,a_ = rv.mass_a_from_Kepler_fit([fit.K_arb[i]/np.sin(np.radians(fit.i_arb[i])),  fit.P_arb[i], fit.e_arb[i],  fit.w_arb[i], fit.M0_arb[i]],1,fit.arb_st_mass)

                fit.mass_arb[i] = float(mass_[0])
                fit.a_arb[i] = float(a_[0])
            else:                
                fit.mass_arb[i] = self.arb_param_gui[7*i + 0].value()  
                fit.a_arb[i]    = self.arb_param_gui[7*i + 1].value()  
            
        fit.npl_arb = j# np.count_nonzero(fit.pl_arb_use.values())
 
    

    
####################################################        

    def initialize_color_dialog(self):

        self.colorDialog = QtWidgets.QColorDialog()
        self.colorDialog.setOption(QtWidgets.QColorDialog.ColorDialogOption.ShowAlphaChannel, True)
        self.colorDialog.setOption(QtWidgets.QColorDialog.ColorDialogOption.DontUseNativeDialog, True)
        #elf.colorDialog.setOptions(QtWidgets.QColorDialog.DontUseNativeDialog |QtWidgets.QColorDialog.NoButtons |QtWidgets.QColorDialog.ShowAlphaChannel)

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
        self.buttonGroup_add_RV_data.setId(self.Button_RV_data_11,11)
        self.buttonGroup_add_RV_data.setId(self.Button_RV_data_12,12)
        self.buttonGroup_add_RV_data.setId(self.Button_RV_data_13,13)
        self.buttonGroup_add_RV_data.setId(self.Button_RV_data_14,14)
        self.buttonGroup_add_RV_data.setId(self.Button_RV_data_15,15)
        self.buttonGroup_add_RV_data.setId(self.Button_RV_data_16,16)
        self.buttonGroup_add_RV_data.setId(self.Button_RV_data_17,17)
        self.buttonGroup_add_RV_data.setId(self.Button_RV_data_18,18)
        self.buttonGroup_add_RV_data.setId(self.Button_RV_data_19,19)
        self.buttonGroup_add_RV_data.setId(self.Button_RV_data_20,20)

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
        self.buttonGroup_remove_RV_data.setId(self.remove_rv_data11,11)
        self.buttonGroup_remove_RV_data.setId(self.remove_rv_data12,12)
        self.buttonGroup_remove_RV_data.setId(self.remove_rv_data13,13)
        self.buttonGroup_remove_RV_data.setId(self.remove_rv_data14,14)
        self.buttonGroup_remove_RV_data.setId(self.remove_rv_data15,15)
        self.buttonGroup_remove_RV_data.setId(self.remove_rv_data16,16)
        self.buttonGroup_remove_RV_data.setId(self.remove_rv_data17,17)
        self.buttonGroup_remove_RV_data.setId(self.remove_rv_data18,18)
        self.buttonGroup_remove_RV_data.setId(self.remove_rv_data19,19)
        self.buttonGroup_remove_RV_data.setId(self.remove_rv_data20,20)
        
        self.buttonGroup_apply_act_data_options.setId(self.Button_apply_act_options_1,1)
        self.buttonGroup_apply_act_data_options.setId(self.Button_apply_act_options_2,2)
        self.buttonGroup_apply_act_data_options.setId(self.Button_apply_act_options_3,3)
        self.buttonGroup_apply_act_data_options.setId(self.Button_apply_act_options_4,4)
        self.buttonGroup_apply_act_data_options.setId(self.Button_apply_act_options_5,5)
        self.buttonGroup_apply_act_data_options.setId(self.Button_apply_act_options_6,6)
        self.buttonGroup_apply_act_data_options.setId(self.Button_apply_act_options_7,7)
        self.buttonGroup_apply_act_data_options.setId(self.Button_apply_act_options_8,8)
        self.buttonGroup_apply_act_data_options.setId(self.Button_apply_act_options_9,9)
        self.buttonGroup_apply_act_data_options.setId(self.Button_apply_act_options_10,10)
        self.buttonGroup_apply_act_data_options.setId(self.Button_apply_act_options_11,11)
        self.buttonGroup_apply_act_data_options.setId(self.Button_apply_act_options_12,12)
        self.buttonGroup_apply_act_data_options.setId(self.Button_apply_act_options_13,13)
        self.buttonGroup_apply_act_data_options.setId(self.Button_apply_act_options_14,14)
        self.buttonGroup_apply_act_data_options.setId(self.Button_apply_act_options_15,15)
        self.buttonGroup_apply_act_data_options.setId(self.Button_apply_act_options_16,16)
        self.buttonGroup_apply_act_data_options.setId(self.Button_apply_act_options_17,17)
        self.buttonGroup_apply_act_data_options.setId(self.Button_apply_act_options_18,18)
        self.buttonGroup_apply_act_data_options.setId(self.Button_apply_act_options_19,19)
        self.buttonGroup_apply_act_data_options.setId(self.Button_apply_act_options_20,20)

        self.buttonGroup_apply_rvs_data_options.setId(self.Button_apply_rvs_options_1,1)
        self.buttonGroup_apply_rvs_data_options.setId(self.Button_apply_rvs_options_2,2)
        self.buttonGroup_apply_rvs_data_options.setId(self.Button_apply_rvs_options_3,3)
        self.buttonGroup_apply_rvs_data_options.setId(self.Button_apply_rvs_options_4,4)
        self.buttonGroup_apply_rvs_data_options.setId(self.Button_apply_rvs_options_5,5)
        self.buttonGroup_apply_rvs_data_options.setId(self.Button_apply_rvs_options_6,6)
        self.buttonGroup_apply_rvs_data_options.setId(self.Button_apply_rvs_options_7,7)
        self.buttonGroup_apply_rvs_data_options.setId(self.Button_apply_rvs_options_8,8)
        self.buttonGroup_apply_rvs_data_options.setId(self.Button_apply_rvs_options_9,9)
        self.buttonGroup_apply_rvs_data_options.setId(self.Button_apply_rvs_options_10,10)
        self.buttonGroup_apply_rvs_data_options.setId(self.Button_apply_rvs_options_11,11)
        self.buttonGroup_apply_rvs_data_options.setId(self.Button_apply_rvs_options_12,12)
        self.buttonGroup_apply_rvs_data_options.setId(self.Button_apply_rvs_options_13,13)
        self.buttonGroup_apply_rvs_data_options.setId(self.Button_apply_rvs_options_14,14)
        self.buttonGroup_apply_rvs_data_options.setId(self.Button_apply_rvs_options_15,15)
        self.buttonGroup_apply_rvs_data_options.setId(self.Button_apply_rvs_options_16,16)
        self.buttonGroup_apply_rvs_data_options.setId(self.Button_apply_rvs_options_17,17)
        self.buttonGroup_apply_rvs_data_options.setId(self.Button_apply_rvs_options_18,18)
        self.buttonGroup_apply_rvs_data_options.setId(self.Button_apply_rvs_options_19,19)
        self.buttonGroup_apply_rvs_data_options.setId(self.Button_apply_rvs_options_20,20)


        self.buttonGroup_apply_tra_data_options.setId(self.Button_apply_tra_options_1,1)
        self.buttonGroup_apply_tra_data_options.setId(self.Button_apply_tra_options_2,2)
        self.buttonGroup_apply_tra_data_options.setId(self.Button_apply_tra_options_3,3)
        self.buttonGroup_apply_tra_data_options.setId(self.Button_apply_tra_options_4,4)
        self.buttonGroup_apply_tra_data_options.setId(self.Button_apply_tra_options_5,5)
        self.buttonGroup_apply_tra_data_options.setId(self.Button_apply_tra_options_6,6)
        self.buttonGroup_apply_tra_data_options.setId(self.Button_apply_tra_options_7,7)
        self.buttonGroup_apply_tra_data_options.setId(self.Button_apply_tra_options_8,8)
        self.buttonGroup_apply_tra_data_options.setId(self.Button_apply_tra_options_9,9)
        self.buttonGroup_apply_tra_data_options.setId(self.Button_apply_tra_options_10,10)
        self.buttonGroup_apply_tra_data_options.setId(self.Button_apply_tra_options_11,11)
        self.buttonGroup_apply_tra_data_options.setId(self.Button_apply_tra_options_12,12)
        self.buttonGroup_apply_tra_data_options.setId(self.Button_apply_tra_options_13,13)
        self.buttonGroup_apply_tra_data_options.setId(self.Button_apply_tra_options_14,14)
        self.buttonGroup_apply_tra_data_options.setId(self.Button_apply_tra_options_15,15)
        self.buttonGroup_apply_tra_data_options.setId(self.Button_apply_tra_options_16,16)
        self.buttonGroup_apply_tra_data_options.setId(self.Button_apply_tra_options_17,17)
        self.buttonGroup_apply_tra_data_options.setId(self.Button_apply_tra_options_18,18)
        self.buttonGroup_apply_tra_data_options.setId(self.Button_apply_tra_options_19,19)
        self.buttonGroup_apply_tra_data_options.setId(self.Button_apply_tra_options_20,20)


        self.buttonGroup_apply_tra_dilution.setId(self.Button_apply_dilution_1,1)
        self.buttonGroup_apply_tra_dilution.setId(self.Button_apply_dilution_2,2)
        self.buttonGroup_apply_tra_dilution.setId(self.Button_apply_dilution_3,3)
        self.buttonGroup_apply_tra_dilution.setId(self.Button_apply_dilution_4,4)
        self.buttonGroup_apply_tra_dilution.setId(self.Button_apply_dilution_5,5)
        self.buttonGroup_apply_tra_dilution.setId(self.Button_apply_dilution_6,6)
        self.buttonGroup_apply_tra_dilution.setId(self.Button_apply_dilution_7,7)
        self.buttonGroup_apply_tra_dilution.setId(self.Button_apply_dilution_8,8)
        self.buttonGroup_apply_tra_dilution.setId(self.Button_apply_dilution_9,9)
        self.buttonGroup_apply_tra_dilution.setId(self.Button_apply_dilution_10,10)
        self.buttonGroup_apply_tra_dilution.setId(self.Button_apply_dilution_11,11)
        self.buttonGroup_apply_tra_dilution.setId(self.Button_apply_dilution_12,12)
        self.buttonGroup_apply_tra_dilution.setId(self.Button_apply_dilution_13,13)
        self.buttonGroup_apply_tra_dilution.setId(self.Button_apply_dilution_14,14)
        self.buttonGroup_apply_tra_dilution.setId(self.Button_apply_dilution_15,15)
        self.buttonGroup_apply_tra_dilution.setId(self.Button_apply_dilution_16,16)
        self.buttonGroup_apply_tra_dilution.setId(self.Button_apply_dilution_17,17)
        self.buttonGroup_apply_tra_dilution.setId(self.Button_apply_dilution_18,18)
        self.buttonGroup_apply_tra_dilution.setId(self.Button_apply_dilution_19,19)
        self.buttonGroup_apply_tra_dilution.setId(self.Button_apply_dilution_20,20)




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
        self.buttonGroup_transit_data.setId(self.Button_transit_data_11,11)
        self.buttonGroup_transit_data.setId(self.Button_transit_data_12,12)
        self.buttonGroup_transit_data.setId(self.Button_transit_data_13,13)
        self.buttonGroup_transit_data.setId(self.Button_transit_data_14,14)
        self.buttonGroup_transit_data.setId(self.Button_transit_data_15,15)
        self.buttonGroup_transit_data.setId(self.Button_transit_data_16,16)
        self.buttonGroup_transit_data.setId(self.Button_transit_data_17,17)
        self.buttonGroup_transit_data.setId(self.Button_transit_data_18,18)
        self.buttonGroup_transit_data.setId(self.Button_transit_data_19,19)
        self.buttonGroup_transit_data.setId(self.Button_transit_data_20,20)

        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data_1,1)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data_2,2)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data_3,3)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data_4,4)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data_5,5)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data_6,6)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data_7,7)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data_8,8)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data_9,9)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data_10,10)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data_11,11)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data_12,12)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data_13,13)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data_14,14)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data_15,15)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data_16,16)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data_17,17)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data_18,18)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data_19,19)
        self.buttonGroup_remove_transit_data.setId(self.remove_transit_data_20,20)

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
        self.buttonGroup_activity_data.setId(self.Button_activity_data_11,11)
        self.buttonGroup_activity_data.setId(self.Button_activity_data_12,12)
        self.buttonGroup_activity_data.setId(self.Button_activity_data_13,13)
        self.buttonGroup_activity_data.setId(self.Button_activity_data_14,14)
        self.buttonGroup_activity_data.setId(self.Button_activity_data_15,15)
        self.buttonGroup_activity_data.setId(self.Button_activity_data_16,16)
        self.buttonGroup_activity_data.setId(self.Button_activity_data_17,17)
        self.buttonGroup_activity_data.setId(self.Button_activity_data_18,18)
        self.buttonGroup_activity_data.setId(self.Button_activity_data_19,19)
        self.buttonGroup_activity_data.setId(self.Button_activity_data_20,20)

        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data_1,1)
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data_2,2)
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data_3,3)
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data_4,4)
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data_5,5)
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data_6,6)
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data_7,7)
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data_8,8)
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data_9,9)
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data_10,10)       
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data_11,11)       
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data_12,12)       
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data_13,13)       
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data_14,14)       
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data_15,15)       
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data_16,16)       
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data_17,17)       
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data_18,18)       
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data_19,19)       
        self.buttonGroup_remove_activity_data.setId(self.remove_activity_data_20,20)       
        

        self.buttonGroup_ttv_data.setId(self.Button_ttv_data_1,1)
        self.buttonGroup_ttv_data.setId(self.Button_ttv_data_2,2)
        self.buttonGroup_ttv_data.setId(self.Button_ttv_data_3,3)
        self.buttonGroup_ttv_data.setId(self.Button_ttv_data_4,4)
        self.buttonGroup_ttv_data.setId(self.Button_ttv_data_5,5)
        self.buttonGroup_ttv_data.setId(self.Button_ttv_data_6,6)
        self.buttonGroup_ttv_data.setId(self.Button_ttv_data_7,7)
        self.buttonGroup_ttv_data.setId(self.Button_ttv_data_8,8)
        self.buttonGroup_ttv_data.setId(self.Button_ttv_data_9,9)
        self.buttonGroup_ttv_data.setId(self.Button_ttv_data_10,10)
        

        self.buttonGroup_remove_ttv_data.setId(self.remove_ttv_data_1,1)
        self.buttonGroup_remove_ttv_data.setId(self.remove_ttv_data_2,2)
        self.buttonGroup_remove_ttv_data.setId(self.remove_ttv_data_3,3)
        self.buttonGroup_remove_ttv_data.setId(self.remove_ttv_data_4,4)
        self.buttonGroup_remove_ttv_data.setId(self.remove_ttv_data_5,5)
        self.buttonGroup_remove_ttv_data.setId(self.remove_ttv_data_6,6)
        self.buttonGroup_remove_ttv_data.setId(self.remove_ttv_data_7,7)
        self.buttonGroup_remove_ttv_data.setId(self.remove_ttv_data_8,8)
        self.buttonGroup_remove_ttv_data.setId(self.remove_ttv_data_9,9)
        self.buttonGroup_remove_ttv_data.setId(self.remove_ttv_data_10,10)
        
        self.buttonGroup_use_ttv_data_to_planet.setId(self.use_ttv_data_1,1)
        self.buttonGroup_use_ttv_data_to_planet.setId(self.use_ttv_data_2,2)
        self.buttonGroup_use_ttv_data_to_planet.setId(self.use_ttv_data_3,3)
        self.buttonGroup_use_ttv_data_to_planet.setId(self.use_ttv_data_4,4)
        self.buttonGroup_use_ttv_data_to_planet.setId(self.use_ttv_data_5,5)
        self.buttonGroup_use_ttv_data_to_planet.setId(self.use_ttv_data_6,6)
        self.buttonGroup_use_ttv_data_to_planet.setId(self.use_ttv_data_7,7)
        self.buttonGroup_use_ttv_data_to_planet.setId(self.use_ttv_data_8,8)
        self.buttonGroup_use_ttv_data_to_planet.setId(self.use_ttv_data_9,9)
        self.buttonGroup_use_ttv_data_to_planet.setId(self.use_ttv_data_10,10)


        self.buttonGroup_ast_data.setId(self.Button_ast_data_1,1)
        self.buttonGroup_ast_data.setId(self.Button_ast_data_2,2)
        self.buttonGroup_ast_data.setId(self.Button_ast_data_3,3)
        self.buttonGroup_ast_data.setId(self.Button_ast_data_4,4)
        self.buttonGroup_ast_data.setId(self.Button_ast_data_5,5)
        self.buttonGroup_ast_data.setId(self.Button_ast_data_6,6)
        self.buttonGroup_ast_data.setId(self.Button_ast_data_7,7)
        self.buttonGroup_ast_data.setId(self.Button_ast_data_8,8)
        self.buttonGroup_ast_data.setId(self.Button_ast_data_9,9)
        self.buttonGroup_ast_data.setId(self.Button_ast_data_10,10)
        
        self.buttonGroup_remove_ast_data.setId(self.remove_ast_data_1,1)
        self.buttonGroup_remove_ast_data.setId(self.remove_ast_data_2,2)
        self.buttonGroup_remove_ast_data.setId(self.remove_ast_data_3,3)
        self.buttonGroup_remove_ast_data.setId(self.remove_ast_data_4,4)
        self.buttonGroup_remove_ast_data.setId(self.remove_ast_data_5,5)
        self.buttonGroup_remove_ast_data.setId(self.remove_ast_data_6,6)
        self.buttonGroup_remove_ast_data.setId(self.remove_ast_data_7,7)
        self.buttonGroup_remove_ast_data.setId(self.remove_ast_data_8,8)
        self.buttonGroup_remove_ast_data.setId(self.remove_ast_data_9,9)
        self.buttonGroup_remove_ast_data.setId(self.remove_ast_data_10,10)
        
        self.buttonGroup_use_ast_data_to_planet.setId(self.use_ast_data_1,1)
        self.buttonGroup_use_ast_data_to_planet.setId(self.use_ast_data_2,2)
        self.buttonGroup_use_ast_data_to_planet.setId(self.use_ast_data_3,3)
        self.buttonGroup_use_ast_data_to_planet.setId(self.use_ast_data_4,4)
        self.buttonGroup_use_ast_data_to_planet.setId(self.use_ast_data_5,5)
        self.buttonGroup_use_ast_data_to_planet.setId(self.use_ast_data_6,6)
        self.buttonGroup_use_ast_data_to_planet.setId(self.use_ast_data_7,7)
        self.buttonGroup_use_ast_data_to_planet.setId(self.use_ast_data_8,8)
        self.buttonGroup_use_ast_data_to_planet.setId(self.use_ast_data_9,9)
        self.buttonGroup_use_ast_data_to_planet.setId(self.use_ast_data_10,10)
     
        self.buttonGroup_ast_data_2.setId(self.Button_ast_data_hipp_1,1)
        self.buttonGroup_ast_data_2.setId(self.Button_ast_data_hipp_2,2)       
        
        self.buttonGroup_remove_ast_data_2.setId(self.remove_ast_data_hipp_1,1)
        self.buttonGroup_remove_ast_data_2.setId(self.remove_ast_data_hipp_2,2)        
        
        self.buttonGroup_use_ast_data_to_planet_2.setId(self.use_ast_data_hipp_1,1)
        self.buttonGroup_use_ast_data_to_planet_2.setId(self.use_ast_data_hipp_2,2)    
        
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
        self.buttonGroup_color_picker.setId(self.rv_pushButton_color_12,12)
        self.buttonGroup_color_picker.setId(self.rv_pushButton_color_13,13)
        self.buttonGroup_color_picker.setId(self.rv_pushButton_color_14,14)
        self.buttonGroup_color_picker.setId(self.rv_pushButton_color_15,15)
        self.buttonGroup_color_picker.setId(self.rv_pushButton_color_16,16)
        self.buttonGroup_color_picker.setId(self.rv_pushButton_color_17,17)
        self.buttonGroup_color_picker.setId(self.rv_pushButton_color_18,18)
        self.buttonGroup_color_picker.setId(self.rv_pushButton_color_19,19)
        self.buttonGroup_color_picker.setId(self.rv_pushButton_color_20,20)
        self.buttonGroup_color_picker.setId(self.rv_pushButton_color_21,21)
        
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
        self.buttonGroup_color_picker_tra.setId(self.trans_pushButton_color_12,12)
        self.buttonGroup_color_picker_tra.setId(self.trans_pushButton_color_13,13)
        self.buttonGroup_color_picker_tra.setId(self.trans_pushButton_color_14,14)
        self.buttonGroup_color_picker_tra.setId(self.trans_pushButton_color_15,15)
        self.buttonGroup_color_picker_tra.setId(self.trans_pushButton_color_16,16)
        self.buttonGroup_color_picker_tra.setId(self.trans_pushButton_color_17,17)
        self.buttonGroup_color_picker_tra.setId(self.trans_pushButton_color_18,18)
        self.buttonGroup_color_picker_tra.setId(self.trans_pushButton_color_19,19)
        self.buttonGroup_color_picker_tra.setId(self.trans_pushButton_color_20,20)
        self.buttonGroup_color_picker_tra.setId(self.trans_pushButton_color_21,21)
        
        self.buttonGroup_color_picker_ttv.setId(self.ttv_color_1,1)
        self.buttonGroup_color_picker_ttv.setId(self.ttv_color_2,2)
        self.buttonGroup_color_picker_ttv.setId(self.ttv_color_3,3)
        self.buttonGroup_color_picker_ttv.setId(self.ttv_color_4,4)
        self.buttonGroup_color_picker_ttv.setId(self.ttv_color_5,5)
        self.buttonGroup_color_picker_ttv.setId(self.ttv_color_6,6)
        self.buttonGroup_color_picker_ttv.setId(self.ttv_color_7,7)
        self.buttonGroup_color_picker_ttv.setId(self.ttv_color_8,8)
        self.buttonGroup_color_picker_ttv.setId(self.ttv_color_9,9)
        self.buttonGroup_color_picker_ttv.setId(self.ttv_color_10,10)
        self.buttonGroup_color_picker_ttv.setId(self.ttv_model_color,11)

        self.buttonGroup_color_picker_ast.setId(self.ast_color_1,1)
        self.buttonGroup_color_picker_ast.setId(self.ast_color_2,2)
        self.buttonGroup_color_picker_ast.setId(self.ast_color_3,3)
        self.buttonGroup_color_picker_ast.setId(self.ast_color_4,4)
        self.buttonGroup_color_picker_ast.setId(self.ast_color_5,5)
        self.buttonGroup_color_picker_ast.setId(self.ast_color_6,6)
        self.buttonGroup_color_picker_ast.setId(self.ast_color_7,7)
        self.buttonGroup_color_picker_ast.setId(self.ast_color_8,8)
        self.buttonGroup_color_picker_ast.setId(self.ast_color_9,9)
        self.buttonGroup_color_picker_ast.setId(self.ast_color_10,10)
        self.buttonGroup_color_picker_ast.setId(self.ast_model_color,11)        
     
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
        self.buttonGroup_symbol_picker.setId(self.rv_pushButton_symbol_11,11)
        self.buttonGroup_symbol_picker.setId(self.rv_pushButton_symbol_12a,12)
        self.buttonGroup_symbol_picker.setId(self.rv_pushButton_symbol_13a,13)
        self.buttonGroup_symbol_picker.setId(self.rv_pushButton_symbol_14a,14)
        self.buttonGroup_symbol_picker.setId(self.rv_pushButton_symbol_15a,15)
        self.buttonGroup_symbol_picker.setId(self.rv_pushButton_symbol_16,16)
        self.buttonGroup_symbol_picker.setId(self.rv_pushButton_symbol_17,17)
        self.buttonGroup_symbol_picker.setId(self.rv_pushButton_symbol_18,18)
        self.buttonGroup_symbol_picker.setId(self.rv_pushButton_symbol_19,19)
        self.buttonGroup_symbol_picker.setId(self.rv_pushButton_symbol_20,20)


        
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
        self.buttonGroup_symbol_picker_tra.setId(self.trans_pushButton_symbol_11,11)
        self.buttonGroup_symbol_picker_tra.setId(self.trans_pushButton_symbol_12,12)
        self.buttonGroup_symbol_picker_tra.setId(self.trans_pushButton_symbol_13,13)
        self.buttonGroup_symbol_picker_tra.setId(self.trans_pushButton_symbol_14,14)
        self.buttonGroup_symbol_picker_tra.setId(self.trans_pushButton_symbol_15,15)
        self.buttonGroup_symbol_picker_tra.setId(self.trans_pushButton_symbol_16,16)
        self.buttonGroup_symbol_picker_tra.setId(self.trans_pushButton_symbol_17,17)
        self.buttonGroup_symbol_picker_tra.setId(self.trans_pushButton_symbol_18,18)
        self.buttonGroup_symbol_picker_tra.setId(self.trans_pushButton_symbol_19,19)
        self.buttonGroup_symbol_picker_tra.setId(self.trans_pushButton_symbol_20,20)


        self.buttonGroup_symbol_picker_ttv.setId(self.ttv_symbol_1,1)
        self.buttonGroup_symbol_picker_ttv.setId(self.ttv_symbol_2,2)
        self.buttonGroup_symbol_picker_ttv.setId(self.ttv_symbol_3,3)
        self.buttonGroup_symbol_picker_ttv.setId(self.ttv_symbol_4,4)
        self.buttonGroup_symbol_picker_ttv.setId(self.ttv_symbol_5,5)
        self.buttonGroup_symbol_picker_ttv.setId(self.ttv_symbol_6,6)
        self.buttonGroup_symbol_picker_ttv.setId(self.ttv_symbol_7,7)
        self.buttonGroup_symbol_picker_ttv.setId(self.ttv_symbol_8,8)
        self.buttonGroup_symbol_picker_ttv.setId(self.ttv_symbol_9,9)
        self.buttonGroup_symbol_picker_ttv.setId(self.ttv_symbol_10,10)

        self.buttonGroup_symbol_picker_ast.setId(self.ast_symbol_1,1)
        self.buttonGroup_symbol_picker_ast.setId(self.ast_symbol_2,2)
        self.buttonGroup_symbol_picker_ast.setId(self.ast_symbol_3,3)
        self.buttonGroup_symbol_picker_ast.setId(self.ast_symbol_4,4)
        self.buttonGroup_symbol_picker_ast.setId(self.ast_symbol_5,5)
        self.buttonGroup_symbol_picker_ast.setId(self.ast_symbol_6,6)
        self.buttonGroup_symbol_picker_ast.setId(self.ast_symbol_7,7)
        self.buttonGroup_symbol_picker_ast.setId(self.ast_symbol_8,8)
        self.buttonGroup_symbol_picker_ast.setId(self.ast_symbol_9,9)
        self.buttonGroup_symbol_picker_ast.setId(self.ast_symbol_10,10)

        self.buttonGroup_detrend_tra.setId(self.tra_detrend_1,1)
        self.buttonGroup_detrend_tra.setId(self.tra_detrend_2,2)
        self.buttonGroup_detrend_tra.setId(self.tra_detrend_3,3)
        self.buttonGroup_detrend_tra.setId(self.tra_detrend_4,4)
        self.buttonGroup_detrend_tra.setId(self.tra_detrend_5,5)
        self.buttonGroup_detrend_tra.setId(self.tra_detrend_6,6)
        self.buttonGroup_detrend_tra.setId(self.tra_detrend_7,7)
        self.buttonGroup_detrend_tra.setId(self.tra_detrend_8,8)
        self.buttonGroup_detrend_tra.setId(self.tra_detrend_9,9)
        self.buttonGroup_detrend_tra.setId(self.tra_detrend_10,10)
        self.buttonGroup_detrend_tra.setId(self.tra_detrend_11,11)
        self.buttonGroup_detrend_tra.setId(self.tra_detrend_12,12)
        self.buttonGroup_detrend_tra.setId(self.tra_detrend_13,13)
        self.buttonGroup_detrend_tra.setId(self.tra_detrend_14,14)
        self.buttonGroup_detrend_tra.setId(self.tra_detrend_15,15)
        self.buttonGroup_detrend_tra.setId(self.tra_detrend_16,16)
        self.buttonGroup_detrend_tra.setId(self.tra_detrend_17,17)
        self.buttonGroup_detrend_tra.setId(self.tra_detrend_18,18)
        self.buttonGroup_detrend_tra.setId(self.tra_detrend_19,19)
        self.buttonGroup_detrend_tra.setId(self.tra_detrend_20,20)

        self.buttonGroup_options_rvs.setId(self.rvs_dataOption_1,1)
        self.buttonGroup_options_rvs.setId(self.rvs_dataOption_2,2)
        self.buttonGroup_options_rvs.setId(self.rvs_dataOption_3,3)
        self.buttonGroup_options_rvs.setId(self.rvs_dataOption_4,4)
        self.buttonGroup_options_rvs.setId(self.rvs_dataOption_5,5)
        self.buttonGroup_options_rvs.setId(self.rvs_dataOption_6,6)
        self.buttonGroup_options_rvs.setId(self.rvs_dataOption_7,7)
        self.buttonGroup_options_rvs.setId(self.rvs_dataOption_8,8)
        self.buttonGroup_options_rvs.setId(self.rvs_dataOption_9,9)
        self.buttonGroup_options_rvs.setId(self.rvs_dataOption_10,10)
        self.buttonGroup_options_rvs.setId(self.rvs_dataOption_11,11)
        self.buttonGroup_options_rvs.setId(self.rvs_dataOption_12,12)
        self.buttonGroup_options_rvs.setId(self.rvs_dataOption_13,13)
        self.buttonGroup_options_rvs.setId(self.rvs_dataOption_14,14)
        self.buttonGroup_options_rvs.setId(self.rvs_dataOption_15,15)
        self.buttonGroup_options_rvs.setId(self.rvs_dataOption_16,16)
        self.buttonGroup_options_rvs.setId(self.rvs_dataOption_17,17)
        self.buttonGroup_options_rvs.setId(self.rvs_dataOption_18,18)
        self.buttonGroup_options_rvs.setId(self.rvs_dataOption_19,19)
        self.buttonGroup_options_rvs.setId(self.rvs_dataOption_20,20)

        self.buttonGroup_options_act.setId(self.act_dataOption_1,1)
        self.buttonGroup_options_act.setId(self.act_dataOption_2,2)
        self.buttonGroup_options_act.setId(self.act_dataOption_3,3)
        self.buttonGroup_options_act.setId(self.act_dataOption_4,4)
        self.buttonGroup_options_act.setId(self.act_dataOption_5,5)
        self.buttonGroup_options_act.setId(self.act_dataOption_6,6)
        self.buttonGroup_options_act.setId(self.act_dataOption_7,7)
        self.buttonGroup_options_act.setId(self.act_dataOption_8,8)
        self.buttonGroup_options_act.setId(self.act_dataOption_9,9)
        self.buttonGroup_options_act.setId(self.act_dataOption_10,10)
        self.buttonGroup_options_act.setId(self.act_dataOption_11,11)
        self.buttonGroup_options_act.setId(self.act_dataOption_12,12)
        self.buttonGroup_options_act.setId(self.act_dataOption_13,13)
        self.buttonGroup_options_act.setId(self.act_dataOption_14,14)
        self.buttonGroup_options_act.setId(self.act_dataOption_15,15)
        self.buttonGroup_options_act.setId(self.act_dataOption_16,16)
        self.buttonGroup_options_act.setId(self.act_dataOption_17,17)
        self.buttonGroup_options_act.setId(self.act_dataOption_18,18)
        self.buttonGroup_options_act.setId(self.act_dataOption_19,19)
        self.buttonGroup_options_act.setId(self.act_dataOption_20,20)


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

        self.colors_gls.setFont(self.plot_font) 
        self.colors_gls_o_c.setFont(self.plot_font) 

        #self.colors_ttv.setFont(self.font) 


    def get_pyqtgraph_plots(self):

        import pyqtgraph.exporters

        global p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,pe,pdi,pcor,p_mlp,p_ttv,p_ttv_oc,p_per_ev,p00,p01,p30,p31,pe2,pe0,pe1,p_ttv_00,p_ttv_01,p_ast, p_ast_00,p_ast_01,p_ast_hipp


        zzz = [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,pe,pdi,pcor,p_mlp,p_ttv,p_ttv_oc,p_per_ev,p_ttv_00,p_ttv_01,p_ast, p_ast_00,p_ast_01,p_ast_hipp]
        zzz_str = ['p1','p2','p3','p4','p5','p6','p7','p8','p9','p10','p11','p12','p13','p14','p15','p16','p17','p18','p19','p20','pe','pdi','pcor',
                  'p_mlp','p_ttv','p_ttv_oc','p_per_ev','pe2',',p_ttv_00','p_ttv_01', 'p_ast', 'p_ast_00', 'p_ast_01','p_ast_hipp']
        for i in range(len(zzz)):

            # create an exporter instance, as an argument give it
            # the item you wish to export

            #zzz[i].scene().getViewWidget()
 
            plots_dir = 'exported_plots'
            if not os.path.exists(plots_dir):
                os.system("mkdir %s"%plots_dir)

            exporter = pg.exporters.ImageExporter(zzz[i].scene())

            # set export parameters if needed
            #exporter.parameters()['width'] = 100   # (note this also affects height parameter)

            # save to file
            exporter.export('%s/%s.png'%(plots_dir,zzz_str[i]))    
        return

    def update_font_plots(self):

        global p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,pe,pdi,pcor,p_mlp,p_ttv,p_ttv_oc,p_per_ev,p00,p01,p30,p31,pe2,pe0,pe1,p_ttv_00,p_ttv_01, p_ast,p_ast_oc, p_ast_00,p_ast_01,p_ast_hipp

        zzz = [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,pe,pdi,pcor,p_mlp,p_ttv,p_ttv_oc,p_per_ev,pe2,p_ast,p_ast_oc,p_ast_hipp]


        for i in range(len(zzz)):

            zzz[i].getAxis('left').setWidth(np.rint(50.0*(float(self.plot_font.pointSize())/11.0)))
            zzz[i].getAxis("left").tickFont = self.plot_font
            zzz[i].getAxis('bottom').setHeight(np.rint(50.0*(float(self.plot_font.pointSize())/11.0)))

            zzz[i].getAxis("bottom").tickFont = self.plot_font
            # zzz[i].getAxis("left").setStyle(tickTextOffset=20)        
            #  zzz[i].getAxis("bottom").setStyle(tickTextOffset=20)        

            zzz[i].setLabel('bottom', '%s'%zzz[i].getAxis("bottom").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
            zzz[i].setLabel('left', '%s'%zzz[i].getAxis("left").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})

        zzz2 = [p00,p30,pe0,p_ttv_00,p_ast_00]
        for i in range(len(zzz2)):

            zzz2[i].getAxis('left').setWidth(np.rint(50.0*(float(self.plot_font.pointSize())/11.0)))
            zzz2[i].getAxis("left").tickFont = self.plot_font
          #  zzz2[i].getAxis('bottom').setHeight(np.rint(50.0*(float(self.plot_font.pointSize())/11.0)))
          #  zzz2[i].getAxis("bottom").tickFont = self.plot_font
          #  zzz[i].getAxis("left").setStyle(tickTextOffset=20)        
          #  zzz[i].getAxis("bottom").setStyle(tickTextOffset=20)        
            
          #  zzz2[i].setLabel('bottom', '%s'%zzz2[i].getAxis("bottom").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
            zzz2[i].setLabel('left', '%s'%zzz2[i].getAxis("left").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})

        zzz3 = [p01,p31,pe1,p_ttv_01,p_ast_01]
        for i in range(len(zzz3)):

            zzz3[i].getAxis('left').setWidth(np.rint(50.0*(float(self.plot_font.pointSize())/11.0)))
            zzz3[i].getAxis("left").tickFont = self.plot_font
            zzz3[i].getAxis('bottom').setHeight(np.rint(50.0*(float(self.plot_font.pointSize())/11.0)))
            zzz3[i].getAxis("bottom").tickFont = self.plot_font
           # zzz3[i].getAxis("left").setStyle(tickTextOffset=20)        
           # zzz3[i].getAxis("bottom").setStyle(tickTextOffset=20)        
            
            zzz3[i].setLabel('bottom', '%s'%zzz3[i].getAxis("bottom").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
            zzz3[i].setLabel('left', '%s'%zzz3[i].getAxis("left").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
 
        return
 

    def initialize_plots(self):

        global p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,pe,pdi,pcor,p_mlp,p_ttv,p_ttv_oc,p_per_ev,pe2,p_ast,p_ast_oc,p_ast_hipp

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
        pe2  = self.graphicsView_extra_plot_gls        
        pdi = self.load_data_plot

        pcor = self.graphicsView_corner
        p_mlp = self.graphicsView_peridogram_RV_mlp
        
        p_ttv = self.graphicsView_timeseries_ttv
        p_ttv_oc = self.graphicsView_timeseries_ttv_o_c
        p_ast = self.graphicsView_timeseries_ast
        p_ast_oc = self.graphicsView_timeseries_ast_o_c
        p_per_ev = self.graphicsView_orb_evol_periods
        p_ast_hipp = self.graphicsView_timeseries_ast_hipp_gaia


        xaxis = ['BJD [days]','BJD [days]','BJD [days]','BJD [days]','BJD [days]','x','Period [d]','Period [d]','Period [d]',
                 'Period [d]','Period [d]','Period [d]','t [yr]','t [yr]','t [yr]','a [au]','t [yr]',
                 't [yr]','t [yr]','t [yr]','','x','x','Period [d]','N transit','N transit','t [yr]','Period [d]','RA [mas]','RA [mas]','RA [mas]']

        if qso_mode:
            yaxis = ['flux [mag.]','flux [mag.]','Rel. Flux','Rel. Flux','y','y','Power','Power','SDE','SDE','Power','Power','a [au]','e',
                 '<html><head/><body><p>&omega; [deg] </p></body></html>','a [au]',
                 '<html><head/><body><p>&Delta;&omega; [deg] </p></body></html>',
                 '<html><head/><body><p>&theta; [deg] </p></body></html>',
                 'i [deg]','energy','','y','y','dlnL','BJD [days]','BJD [days]','Period rat.','Power','au','au','au']
        else:

            yaxis = ['RV [m/s]','RV [m/s]','Rel. Flux','Rel. Flux','y','y','Power','Power','SDE','SDE','Power','Power','a [au]','e',
                 '<html><head/><body><p>&omega; [deg] </p></body></html>','a [au]',
                 '<html><head/><body><p>&Delta;&omega; [deg] </p></body></html>',
                 '<html><head/><body><p>&theta; [deg] </p></body></html>',
                 'i [deg]','energy','','y','y','dlnL','BJD [days]','BJD [days]','Period rat.','Power','DEC [mas]','DEC [mas]','DEC [mas]']
        xunit = ['','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','']
        yunit = ['','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','']
 
        zzz = [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,pe,pdi,pcor,p_mlp,p_ttv,p_ttv_oc,p_per_ev,pe2,p_ast,p_ast_oc,p_ast_hipp]
 

        for i in range(len(zzz)):

            zzz[i].setAxisItems({'bottom': pg_hack.CustomAxisItem('bottom')})

            #zzz[i].getAxis("bottom").tickFont = self.plot_font
            zzz[i].getAxis("bottom").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
            zzz[i].getAxis("top").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
            zzz[i].getAxis("left").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
            zzz[i].getAxis("right").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
           # zzz[i].getAxis('left').setWidth(50)
            zzz[i].getAxis('right').setWidth(10)
            zzz[i].getAxis('top').setHeight(10)
           # zzz[i].getAxis('bottom').setHeight(50)
                        
            zzz[i].setLabel('bottom', '%s'%xaxis[i], units='%s'%xunit[i],  **{'font-size':'%dpt'%self.plot_font.pointSize()})
            zzz[i].setLabel('left',   '%s'%yaxis[i], units='%s'%yunit[i],  **{'font-size':'%dpt'%self.plot_font.pointSize()})       
            zzz[i].showAxis('top') 
            zzz[i].showAxis('right') 
            zzz[i].getAxis('bottom').enableAutoSIPrefix(enable=False)
            zzz[i].getAxis('left').enableAutoSIPrefix(enable=False)            
            #zzz[i].autoRange()

            zzz[i].getViewBox().parentItem().ctrlMenu.actions()[-4].setVisible(False) #removes the "Avarage" junk

        p16.getViewBox().setAspectLocked(True)

     
        self.initialize_RV_subplots()
        self.initialize_tra_subplots()
        self.initialize_TTV_subplots()
        self.initialize_ast_subplots()
        self.initialize_RV_phase_subplots()               


        return


    def initialize_RV_phase_subplots(self):

        global pe,pe0,pe1,legend_RV_phased
        
        l = pg.GraphicsLayout()                                                             
        pe.setCentralItem(l)            

        if qso_mode:
            o_c_label = "o-c [mag.]"
        else:
            o_c_label = "o-c [m/s]"   
            
        pe0 = l.addPlot(0, 0, colspan=3)                                                                
        pe0.hideAxis('bottom')                                                              
        pe1 = l.addPlot(1, 0, colspan=1)                                                                
        pe1.setXLink(pe0)                 

        legend_RV_phased = pe0.addLegend()

        #for i in (1, 2):
        l.layout.setRowMinimumHeight(0, 220)                                                    
        l.layout.setRowMinimumHeight(1, 30)         
        l.layout.setRowMaximumHeight(1, 150)        
        
                                                                                                       
        pe0.showAxis('top') 
        pe0.showAxis('right') 
        pe1.showAxis('top') 
        pe1.showAxis('right') 
 
        pe0.ctrlMenu.actions()[-4].setVisible(False) #removes the submenu "Avarage" junk       
        pe1.ctrlMenu.actions()[-4].setVisible(False) #removes the submenu "Avarage" junk       

        #pe0.getAxis('left').setStyle(autoExpandTextSpace=True)
        pe0.getAxis('left').setWidth(np.rint(60.0*(float(self.plot_font.pointSize())/11.0)))
        pe0.getAxis("left").tickFont = self.plot_font
        #pe0.getAxis('bottom').setHeight(np.rint(50.0*(float(self.plot_font.pointSize())/11.0)))
        pe0.getAxis("bottom").tickFont = self.plot_font
        #pe0.setLabel('bottom', '%s'%p1.getAxis("bottom").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        pe0.setLabel('left', '%s'%p1.getAxis("left").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        pe0.getAxis('right').setWidth(0)
        pe0.getAxis('top').setHeight(0)

        pe0.setAxisItems({'bottom': pg_hack.CustomAxisItem('bottom')})

        #pe0.getViewBox().setAspectLocked(lock=False, ratio=2)

        #pe0.getAxis("bottom").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        pe0.getAxis("top").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        pe0.getAxis("left").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        pe0.getAxis("right").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        #pe0.getAxis('bottom').enableAutoSIPrefix(enable=False)

        #pe1.getAxis('left').setStyle(autoExpandTextSpace=True)

        pe1.getAxis('left').setWidth(np.rint(60.0*(float(self.plot_font.pointSize())/11.0)))
        pe1.getAxis("left").tickFont = self.plot_font
        pe1.getAxis('bottom').setHeight(np.rint(50.0*(float(self.plot_font.pointSize())/11.0)))
        pe1.getAxis("bottom").tickFont = self.plot_font

        #pe1.setLabel('left', '%s'%p2.getAxis("left").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        pe1.setLabel('left', o_c_label, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        pe1.getAxis('right').setWidth(0)
        pe1.getAxis('top').setHeight(0)

        pe1.setAxisItems({'bottom': pg_hack.CustomAxisItem('bottom')})
        pe1.setLabel('bottom', '%s'%p2.getAxis("bottom").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})

 
        pe1.getAxis("bottom").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        pe1.getAxis("top").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        pe1.getAxis("left").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        pe1.getAxis("right").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        pe1.getAxis('bottom').enableAutoSIPrefix(enable=False)
        pe1.getAxis('left').enableAutoSIPrefix(enable=False)

        #pe0.getAxis('left').label.setStyle(offset=20)  # Add this line
        #pe1.getAxis('left').label.setStyle(offset=20)  # Add this line        
        ax0 = pe0.getAxis('bottom')      #get handle to x-axis 0
        ax0.setStyle(showValues=False)

        l.layout.setSpacing(0.)                                                             
        l.setContentsMargins(0., 0., 0., 0.)                                                

 
        return

    def initialize_RV_subplots(self):

        global p1,p00,p01,legend_RV

        if qso_mode:
            o_c_label = "o-c [mag.]"
        else:
            o_c_label = "o-c [m/s]"    
       
        l = pg.GraphicsLayout()                                                             
        p1.setCentralItem(l)            
                                                                                                      

        p00 = l.addPlot(0, 0, colspan=3)                                                                
        p00.hideAxis('bottom')                                                              
        p01 = l.addPlot(1, 0, colspan=1)                                                                
        p01.setXLink(p00)                 

        legend_RV = p00.addLegend()

        #for i in (1, 2):
        l.layout.setRowMinimumHeight(0, 220)                                                    
        l.layout.setRowMinimumHeight(1, 30)         
        l.layout.setRowMaximumHeight(1, 150)                                                                                               
        p00.showAxis('top') 
        p00.showAxis('right') 
        p01.showAxis('top') 
        p01.showAxis('right') 
 
        p00.ctrlMenu.actions()[-4].setVisible(False) #removes the submenu "Avarage" junk       
        p01.ctrlMenu.actions()[-4].setVisible(False) #removes the submenu "Avarage" junk       


        p00.getAxis('left').setWidth(np.rint(60.0*(float(self.plot_font.pointSize())/11.0)))
        p00.getAxis("left").tickFont = self.plot_font
        #p00.getAxis('bottom').setHeight(np.rint(50.0*(float(self.plot_font.pointSize())/11.0)))
        p00.getAxis("bottom").tickFont = self.plot_font
        #p00.setLabel('bottom', '%s'%p1.getAxis("bottom").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        p00.setLabel('left', '%s'%p1.getAxis("left").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        p00.getAxis('right').setWidth(0)
        p00.getAxis('top').setHeight(0)

        p00.setAxisItems({'bottom': pg_hack.CustomAxisItem('bottom')})

        #p00.getViewBox().setAspectLocked(lock=False, ratio=2)

        p00.getAxis("bottom").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p00.getAxis("top").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p00.getAxis("left").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p00.getAxis("right").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        #p00.getAxis('bottom').enableAutoSIPrefix(enable=False)


        p01.getAxis('left').setWidth(np.rint(60.0*(float(self.plot_font.pointSize())/11.0)))
        p01.getAxis("left").tickFont = self.plot_font
        p01.getAxis('bottom').setHeight(np.rint(50.0*(float(self.plot_font.pointSize())/11.0)))
        p01.getAxis("bottom").tickFont = self.plot_font

        #p01.setLabel('left', '%s'%p2.getAxis("left").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        p01.setLabel('left', o_c_label, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        p01.getAxis('right').setWidth(0)
        p01.getAxis('top').setHeight(0)

        p01.setAxisItems({'bottom': pg_hack.CustomAxisItem('bottom')})
        p01.setLabel('bottom', '%s'%p2.getAxis("bottom").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})


        p01.getAxis("bottom").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p01.getAxis("top").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p01.getAxis("left").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p01.getAxis("right").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p01.getAxis('bottom').enableAutoSIPrefix(enable=False)
        p01.getAxis('left').enableAutoSIPrefix(enable=False)

        ax0 = p00.getAxis('bottom')      #get handle to x-axis 0
        ax0.setStyle(showValues=False)

        l.layout.setSpacing(0.)                                                             
        l.setContentsMargins(0., 0., 0., 0.)                                                


        return


    def initialize_tra_subplots(self):

        global p3,p30,p31,legend_tra
        
        ll = pg.GraphicsLayout()                                                             
        p3.setCentralItem(ll)                                                                                                                        
        
        p30 = ll.addPlot(0, 0, colspan=3)                                                               
        p30.hideAxis('bottom')                                                              
        p31 = ll.addPlot(1, 0, colspan=1)                                                               
        p31.setXLink(p30)           

        legend_tra = p30.addLegend()                                                          
         
        ll.layout.setRowMinimumHeight(0, 220)                                                    
        ll.layout.setRowMinimumHeight(1, 30)         
        ll.layout.setRowMaximumHeight(1, 150)                                                                                               
 

        p30.showAxis('top') 
        p30.showAxis('right') 
        p31.showAxis('top') 
        p31.showAxis('right') 
  
        p30.ctrlMenu.actions()[-4].setVisible(False) #removes the submenu "Avarage" junk       
        p31.ctrlMenu.actions()[-4].setVisible(False) #removes the submenu "Avarage" junk       


        #p30.autoRange(False)
        #p31.autoRange(False)
        
        p30.getAxis('left').setWidth(np.rint(60.0*(float(self.plot_font.pointSize())/11.0)))
        p30.getAxis("left").tickFont = self.plot_font
        #p00.getAxis('bottom').setHeight(np.rint(50.0*(float(self.plot_font.pointSize())/11.0)))
        p30.getAxis("bottom").tickFont = self.plot_font
        #p00.setLabel('bottom', '%s'%p1.getAxis("bottom").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        
        p30.getAxis('right').setWidth(0)
        p30.getAxis('top').setHeight(0)
        p30.setAxisItems({'bottom': pg_hack.CustomAxisItem('bottom')})
        
        p30.setLabel('left', '%s'%p3.getAxis("left").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        
        p30.getAxis("bottom").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p30.getAxis("top").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p30.getAxis("left").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p30.getAxis("right").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p30.getAxis('bottom').enableAutoSIPrefix(enable=False)
        p30.getAxis('left').enableAutoSIPrefix(enable=False)        
         
        
        p31.getAxis('left').setWidth(np.rint(60.0*(float(self.plot_font.pointSize())/11.0)))
        p31.getAxis("left").tickFont = self.plot_font
        p31.getAxis('bottom').setHeight(np.rint(50.0*(float(self.plot_font.pointSize())/11.0)))
        p31.getAxis("bottom").tickFont = self.plot_font
        
        #p31.setLabel('left', '%s'%p4.getAxis("left").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        p31.setLabel('left', 'o-c', units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        p31.getAxis('right').setWidth(0)
        p31.getAxis('top').setHeight(0)
        p31.setAxisItems({'bottom': pg_hack.CustomAxisItem('bottom')})
        
        p31.getAxis("bottom").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p31.getAxis("top").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p31.getAxis("left").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p31.getAxis("right").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p31.getAxis('bottom').enableAutoSIPrefix(enable=False)
        p31.getAxis('left').enableAutoSIPrefix(enable=False)

        
        p31.setLabel('bottom', '%s'%p4.getAxis("bottom").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        
        ax30 = p30.getAxis('bottom')      #get handle to x-axis 0
        ax30.setStyle(showValues=False)
        
        ll.layout.setSpacing(0.)                                                             
        ll.setContentsMargins(0., 0., 0., 0.)                                                
        return


    def initialize_TTV_subplots(self):

        global p_ttv_oc, p_ttv_00,p_ttv_01,legend_ttv

        o_c_label = "o-c [d]"    
       
        l = pg.GraphicsLayout()                                                             
        p_ttv_oc.setCentralItem(l)                                                                                                                        

        p_ttv_00 = l.addPlot(0, 0, colspan=3)                                                                
        p_ttv_00.hideAxis('bottom')                                                              
        p_ttv_01 = l.addPlot(1, 0, colspan=1)                                                                
        p_ttv_01.setXLink(p_ttv_00)      

        legend_ttv = p_ttv_00.addLegend()                                                                     


        #for i in (1, 2):
        l.layout.setRowMinimumHeight(0, 220)                                                    
        l.layout.setRowMinimumHeight(1, 30)         
        l.layout.setRowMaximumHeight(1, 150)                                                                                               
        p_ttv_00.showAxis('top') 
        p_ttv_00.showAxis('right') 
        p_ttv_01.showAxis('top') 
        p_ttv_01.showAxis('right') 
 
        p_ttv_00.ctrlMenu.actions()[-4].setVisible(False) #removes the submenu "Avarage" junk       
        p_ttv_01.ctrlMenu.actions()[-4].setVisible(False) #removes the submenu "Avarage" junk       


        p_ttv_00.getAxis('left').setWidth(np.rint(60.0*(float(self.plot_font.pointSize())/11.0)))
        p_ttv_00.getAxis("left").tickFont = self.plot_font
        #p00.getAxis('bottom').setHeight(np.rint(50.0*(float(self.plot_font.pointSize())/11.0)))
        p_ttv_00.getAxis("bottom").tickFont = self.plot_font
        #p00.setLabel('bottom', '%s'%p1.getAxis("bottom").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        p_ttv_00.setLabel('left', '%s'%p_ttv_oc.getAxis("left").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        p_ttv_00.getAxis('right').setWidth(0)
        p_ttv_00.getAxis('top').setHeight(0)

        p_ttv_00.setAxisItems({'bottom': pg_hack.CustomAxisItem('bottom')})

        #p00.getViewBox().setAspectLocked(lock=False, ratio=2)

        p_ttv_00.getAxis("bottom").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p_ttv_00.getAxis("top").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p_ttv_00.getAxis("left").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p_ttv_00.getAxis("right").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p_ttv_00.getAxis('left').enableAutoSIPrefix(enable=False)


        p_ttv_01.getAxis('left').setWidth(np.rint(60.0*(float(self.plot_font.pointSize())/11.0)))
        p_ttv_01.getAxis("left").tickFont = self.plot_font
        p_ttv_01.getAxis('bottom').setHeight(np.rint(50.0*(float(self.plot_font.pointSize())/11.0)))
        p_ttv_01.getAxis("bottom").tickFont = self.plot_font

        #p_ttv_01.setLabel('left', '%s'%p2.getAxis("left").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        p_ttv_01.setLabel('left', o_c_label, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        p_ttv_01.getAxis('right').setWidth(0)
        p_ttv_01.getAxis('top').setHeight(0)

        p_ttv_01.setAxisItems({'bottom': pg_hack.CustomAxisItem('bottom')})
        p_ttv_01.setLabel('bottom', '%s'%p_ttv_oc.getAxis("bottom").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})


        p_ttv_01.getAxis("bottom").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p_ttv_01.getAxis("top").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p_ttv_01.getAxis("left").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p_ttv_01.getAxis("right").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p_ttv_01.getAxis('bottom').enableAutoSIPrefix(enable=False)
        p_ttv_01.getAxis('left').enableAutoSIPrefix(enable=False)

        ax0 = p_ttv_00.getAxis('bottom')      #get handle to x-axis 0
        ax0.setStyle(showValues=False)

        l.layout.setSpacing(0.)                                                             
        l.setContentsMargins(0., 0., 0., 0.)                                                
 
        return



    def initialize_ast_subplots(self):

        global p_ast_oc, p_ast_00,p_ast_01,legend_ast,legend_ast_hipp

        o_c_label = "o-c [d]"    
       
        l = pg.GraphicsLayout()                                                             
        p_ast_oc.setCentralItem(l)                                                                                                                        

        p_ast_00 = l.addPlot(0, 0, colspan=3)                                                                
        p_ast_00.hideAxis('bottom')                                                              
        p_ast_01 = l.addPlot(1, 0, colspan=1)                                                                
        p_ast_01.setXLink(p_ast_00)      

        legend_ast = p_ast_00.addLegend()                                                                     
        legend_ast_hipp = p_ast_hipp.addLegend()                                                                     

        #for i in (1, 2):
        l.layout.setRowMinimumHeight(0, 220)                                                    
        l.layout.setRowMinimumHeight(1, 30)         
        l.layout.setRowMaximumHeight(1, 150)                                                                                               
        p_ast_00.showAxis('top') 
        p_ast_00.showAxis('right') 
        p_ast_01.showAxis('top') 
        p_ast_01.showAxis('right') 
 
        p_ast_00.ctrlMenu.actions()[-4].setVisible(False) #removes the submenu "Avarage" junk       
        p_ast_01.ctrlMenu.actions()[-4].setVisible(False) #removes the submenu "Avarage" junk       


        p_ast_00.getAxis('left').setWidth(np.rint(60.0*(float(self.plot_font.pointSize())/11.0)))
        p_ast_00.getAxis("left").tickFont = self.plot_font
        #p00.getAxis('bottom').setHeight(np.rint(50.0*(float(self.plot_font.pointSize())/11.0)))
        p_ast_00.getAxis("bottom").tickFont = self.plot_font
        #p00.setLabel('bottom', '%s'%p1.getAxis("bottom").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        p_ast_00.setLabel('left', '%s'%p_ast_oc.getAxis("left").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        p_ast_00.getAxis('right').setWidth(0)
        p_ast_00.getAxis('top').setHeight(0)

        p_ast_00.setAxisItems({'bottom': pg_hack.CustomAxisItem('bottom')})

        #p00.getViewBox().setAspectLocked(lock=False, ratio=2)

        p_ast_00.getAxis("bottom").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p_ast_00.getAxis("top").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p_ast_00.getAxis("left").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p_ast_00.getAxis("right").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        #p00.getAxis('bottom').enableAutoSIPrefix(enable=False)


        p_ast_01.getAxis('left').setWidth(np.rint(60.0*(float(self.plot_font.pointSize())/11.0)))
        p_ast_01.getAxis("left").tickFont = self.plot_font
        p_ast_01.getAxis('bottom').setHeight(np.rint(50.0*(float(self.plot_font.pointSize())/11.0)))
        p_ast_01.getAxis("bottom").tickFont = self.plot_font

        #p_ast_01.setLabel('left', '%s'%p2.getAxis("left").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        p_ast_01.setLabel('left', o_c_label, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        p_ast_01.getAxis('right').setWidth(0)
        p_ast_01.getAxis('top').setHeight(0)

        p_ast_01.setAxisItems({'bottom': pg_hack.CustomAxisItem('bottom')})
        p_ast_01.setLabel('bottom', '%s'%p_ast_oc.getAxis("bottom").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})


        p_ast_01.getAxis("bottom").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p_ast_01.getAxis("top").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p_ast_01.getAxis("left").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p_ast_01.getAxis("right").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        p_ast_01.getAxis('bottom').enableAutoSIPrefix(enable=False)
        p_ast_01.getAxis('left').enableAutoSIPrefix(enable=False)

        ax0 = p_ast_00.getAxis('bottom')      #get handle to x-axis 0
        ax0.setStyle(showValues=False)

        l.layout.setSpacing(0.)                                                             
        l.setContentsMargins(0., 0., 0., 0.)                                                
 
        return

    def initialize_phase_subplots(self):

        global p3,pe0,pe1
        
        ll = pg.GraphicsLayout()                                                             
        pe.setCentralItem(ll)                                                                                                                        
        
        pe0 = ll.addPlot(0, 0, colspan=3)                                                               
        pe0.hideAxis('bottom')                                                              
        pe1 = ll.addPlot(1, 0, colspan=1)                                                               
        pe1.setXLink(pe0)                                                                     
         
        ll.layout.setRowMinimumHeight(0, 220)                                                    
        ll.layout.setRowMinimumHeight(1, 30)         
        #ll.layout.setRowMaximumHeight(1, 150)                                                                                               
          
        pe0.showAxis('top') 
        pe0.showAxis('right') 
        pe1.showAxis('top') 
        pe1.showAxis('right') 
  
        pe0.ctrlMenu.actions()[-4].setVisible(False) #removes the submenu "Avarage" junk       
        pe1.ctrlMenu.actions()[-4].setVisible(False) #removes the submenu "Avarage" junk       


        #pe0.autoRange(False)
        #pe1.autoRange(False)
        
        pe0.getAxis('left').setWidth(np.rint(60.0*(float(self.plot_font.pointSize())/11.0)))
        pe0.getAxis("left").tickFont = self.plot_font
        #p00.getAxis('bottom').setHeight(np.rint(50.0*(float(self.plot_font.pointSize())/11.0)))
        pe0.getAxis("bottom").tickFont = self.plot_font
        #p00.setLabel('bottom', '%s'%p1.getAxis("bottom").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        
        pe0.getAxis('right').setWidth(0)
        pe0.getAxis('top').setHeight(0)
        pe0.setAxisItems({'bottom': pg_hack.CustomAxisItem('bottom')})
        
        pe0.setLabel('left', '%s'%p3.getAxis("left").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
       
        pe0.getAxis("bottom").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        pe0.getAxis("top").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        pe0.getAxis("left").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        pe0.getAxis("right").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        pe0.getAxis('bottom').enableAutoSIPrefix(enable=False)
        pe0.getAxis('left').enableAutoSIPrefix(enable=False)        
         
        
        pe1.getAxis('left').setWidth(np.rint(60.0*(float(self.plot_font.pointSize())/11.0)))
        pe1.getAxis("left").tickFont = self.plot_font
        pe1.getAxis('bottom').setHeight(np.rint(50.0*(float(self.plot_font.pointSize())/11.0)))
        pe1.getAxis("bottom").tickFont = self.plot_font
        
        #pe1.setLabel('left', '%s'%p4.getAxis("left").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        pe1.setLabel('left', 'o-c', units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        pe1.getAxis('right').setWidth(0)
        pe1.getAxis('top').setHeight(0)
        pe1.setAxisItems({'bottom': pg_hack.CustomAxisItem('bottom')})
        
        pe1.getAxis("bottom").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        pe1.getAxis("top").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        pe1.getAxis("left").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        pe1.getAxis("right").setStyle(tickTextOffset = 12, tickFont = self.plot_font)
        pe1.getAxis('bottom').enableAutoSIPrefix(enable=False)
        pe1.getAxis('left').enableAutoSIPrefix(enable=False)
        
        pe1.setLabel('bottom', '%s'%pe.getAxis("bottom").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        
        ax30 = pe0.getAxis('bottom')      #get handle to x-axis 0
        ax30.setStyle(showValues=False)
        
        ll.layout.setSpacing(0.)                                                             
        ll.setContentsMargins(0., 0., 0., 0.)                                                
        return



    def identify_power_peaks(self,x,y,sig_level=np.array([]), power_level=np.array([]) ):
 
        per_ind = argrelextrema(y, np.greater)
        per_x   = x[per_ind]
        per_y   = y[per_ind]

        peaks_sort = sorted(range(len(per_y)), key=lambda k: per_y[k], reverse=True)

        per_x   = per_x[peaks_sort]   
        per_y   = per_y[peaks_sort]

        peaks_pos = [per_x,per_y]

        if power_level is None or sig_level is None:
            power_level = np.array([])
            sig_level = np.array([])


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

        return text_peaks , peaks_pos 




######################## SciPy setup ######################################

    def init_scipy_combo(self):
        global fit 

        for i in range(len(fit.SciPy_min)):
            self.comboBox_scipy_minimizer_1.addItem('%s'%(fit.SciPy_min[i]),i) 
            self.comboBox_scipy_minimizer_2.addItem('%s'%(fit.SciPy_min[i]),i) 
           
        self.comboBox_scipy_minimizer_1.setCurrentIndex(6)
        self.comboBox_scipy_minimizer_2.setCurrentIndex(0)



    def init_scipy_ttv_interpolation(self):
        global fit 


        for i in range(len(ttv_interpol_opt)):
            self.comboBox_ttv_model_interpolate.addItem('%s'%(ttv_interpol_opt[i]),i) 
           
        self.comboBox_ttv_model_interpolate.setCurrentIndex(0)


    def set_AMD_consttraints(self):
        global fit
        
        if self.use_optim_stab_constraints.isChecked():
            fit.optim_AMD_stab   = self.optim_AMD.isChecked()
            fit.optim_Nbody_stab = self.optim_Nbody.isChecked()
        else:
            fit.optim_AMD_stab   = False
            fit.optim_Nbody_stab = False

    def check_scipy_min(self):
        global fit

        ind_min_1 = self.comboBox_scipy_minimizer_1.currentIndex()
        ind_min_2 = self.comboBox_scipy_minimizer_2.currentIndex()

        fit.SciPy_min_use_1 = fit.SciPy_min[ind_min_1]
        fit.SciPy_min_use_2 = fit.SciPy_min[ind_min_2]
        fit.SciPy_min_N_use_1 = int(self.scipy_N_consecutive_iter_1.value())
        fit.SciPy_min_N_use_2 = int(self.scipy_N_consecutive_iter_2.value())

        self.set_AMD_consttraints()

        fit.Simplex_opt    = {'disp': True, 'maxiter': int(self.simplex_maxiter.value()), 'return_all': False, 'maxfev': int(self.simplex_maxfev.value()), 'xtol':self.simplex_xtol.value() , 'ftol': self.simplex_ftol.value() ,'adaptive':True }
        fit.Powell_opt     = {'disp': True, 'return_all': False, 'maxiter': int(self.powell_maxiter.value()), 'direc': None, 'func': None, 'maxfev': int(self.powell_maxfev.value()), 'xtol': self.powell_xtol.value(), 'ftol': self.powell_ftol.value()}
        fit.CG_opt         = {'disp': True, 'gtol': self.cg_gtol.value(), 'eps': 1.4901161193847656e-08, 'return_all': False, 'maxiter': int(self.cg_maxiter.value()), 'norm': np.inf}
        fit.BFGS_opt       = {'disp': True, 'gtol': self.bfgs_gtol.value(), 'eps': 1.4901161193847656e-08, 'return_all': False, 'maxiter': int(self.bfgs_maxiter.value()), 'norm': np.inf}
        fit.Newton_cg_opt  = {'disp': True, 'xtol': self.Newton_cg_xtol.value(), 'eps': 1.4901161193847656e-08, 'return_all': False, 'maxiter': int(self.Newton_cg_maxiter.value())} 
        fit.L_BFGS_B_opt   = {'disp': True, 'maxcor': int(self.LBFGSB_maxcor.value()), 'ftol': 2.220446049250313e-09, 'gtol': self.LBFGSB_gtol.value(), 'eps': 1e-08, 'maxfun': int(self.LBFGSB_maxiter.value()), 'maxiter': int(self.LBFGSB_maxiter.value()), 'iprint': -1, 'maxls': 20}    
        fit.TNC_opt        = {'disp': True, 'eps': self.TNC_eps.value(), 'scale': None, 'offset': None, 'mesg_num': None, 'maxCGit': int(self.TNC_maxcgit.value()), 'maxiter': int(self.TNC_maxiter.value()), 'eta': self.TNC_eta.value(), 'stepmx':self.TNC_stepmx.value(), 'accuracy': self.TNC_accuracy.value(), 'minfev': self.TNC_minfev.value(), 'ftol': self.TNC_ftol.value(), 'xtol':self.TNC_ftol.value(), 'gtol': self.TNC_gtol.value(), 'rescale': -1 }  
       # fit.COBYLA_opt     = {'disp': True, 'rhobeg': self.cobyla_rhobeg.value(), 'maxiter':  int(self.cobyla_maxiter.value()), 'catol': self.cobyla_catol.value() }
        fit.SLSQP_opt      = {'disp': True, 'maxiter': int(self.slsqp_maxiter.value()),  'eps': 1.4901161193847656e-08, 'ftol': self.slsqp_ftol.value(), 'iprint': 1}



######################## Cross hair and label picks ######################################


    def cross_hair(self, plot_wg, log=False, alias = [False, 365.250,'#666699']):
        global fit 

        vLine = pg.InfiniteLine(angle=90, movable=False)#, pos=0)
        hLine = pg.InfiniteLine(angle=0,  movable=False)#, pos=2450000)
        plot_wg.addItem(vLine, ignoreBounds=True)
        plot_wg.addItem(hLine, ignoreBounds=True)
        label = pg.TextItem()

        plot_wg.addItem(label, ignoreBounds=True)  
         
        vb = plot_wg.getViewBox()   
        viewrange = vb.viewRange()
        
        if alias[0] == True:
            v2aLine = pg.InfiniteLine(angle=90, movable=False, pen=alias[2])#, pos=0)
            v2bLine = pg.InfiniteLine(angle=90, movable=False, pen=alias[2])#, pos=0)

            plot_wg.addItem(v2aLine, ignoreBounds=True)
            plot_wg.addItem(v2bLine, ignoreBounds=True)

            v3aLine = pg.InfiniteLine(angle=90, movable=False, pen=alias[2])#, pos=0)
            v3bLine = pg.InfiniteLine(angle=90, movable=False, pen=alias[2])#, pos=0)

            plot_wg.addItem(v3aLine, ignoreBounds=True)
            plot_wg.addItem(v3bLine, ignoreBounds=True)

        def mouseMoved(evt):

            pos = evt[0]  ## using signal proxy turns original arguments into a tuple
            if plot_wg.sceneBoundingRect().contains(pos):

                mousePoint = vb.mapSceneToView(pos)

                if log == True:
                    label.setText("x=%0.3f,  y=%0.3f"%(10**mousePoint.x(), mousePoint.y()))
                else:
                    label.setText("x=%0.3f,  y=%0.3f"%(mousePoint.x(), mousePoint.y()))
                    #label.rotateAxis=(1, 0)

                vLine.setPos(mousePoint.x())
                hLine.setPos(mousePoint.y())
                
                if alias[0] == True:
                    if log == True:
                        v2aLine.setPos(np.log10((1.0 / ( (1.0/(10**mousePoint.x()) ) + 1.0/alias[1] )) ))
                        v2bLine.setPos(np.log10((1.0 / ( (1.0/(10**mousePoint.x()) ) - 1.0/alias[1] )) ))
                        v3aLine.setPos(np.log10((1.0 / ( 1.0/alias[1] + (1.0/(10**mousePoint.x()) ) )) ))
                        v3bLine.setPos(np.log10((1.0 / ( 1.0/alias[1] - (1.0/(10**mousePoint.x()) ) )) ))
                    else:
                        v2aLine.setPos(  (mousePoint.x()  + 1.0/alias[1]) )  
                        v2bLine.setPos(  (mousePoint.x()  - 1.0/alias[1]) )  
                        v3aLine.setPos(  (1.0/alias[1] + mousePoint.x() ))   
                        v3bLine.setPos(  (1.0/alias[1] - mousePoint.x() ))                         
                        

                if mousePoint.x() < (viewrange[0][1]+viewrange[0][0])/2.0:
                    label.setAnchor((0,1))
                else:
                    label.setAnchor((1,1))
                label.setPos(mousePoint.x(), mousePoint.y())
                #fit.label = label

        plot_wg.getViewBox().setAutoVisible(y=True)

        proxy = pg.SignalProxy(plot_wg.scene().sigMouseMoved, rateLimit=60, slot=mouseMoved)
        plot_wg.proxy = proxy

    def cross_hair_remove(self, plot_wg):
        global fit 
 
        for kk in plot_wg.items():
            if kk.__class__.__name__ == "InfiniteLine":
                if kk._name != "zero":
                    plot_wg.removeItem(kk)
            elif kk.__class__.__name__ == "TextItem":
                plot_wg.removeItem(kk)
                


    def label_peaks(self, plot_wg2, pos_peaks, GLS = True, o_c = False, activity = False, MLP = False, DFT=False):

        if GLS == True and DFT == False and self.avoid_GLS_RV_alias.isChecked():
            x_peaks = pos_peaks[0][pos_peaks[0]>1.2]
            y_peaks = pos_peaks[1][pos_peaks[0]>1.2]
        else:
            x_peaks = pos_peaks[0]
            y_peaks = pos_peaks[1]


        if GLS == True:
            N_peaks = int(self.N_GLS_peak_to_point.value())
            if o_c == True:
                log = self.radioButton_RV_o_c_GLS_period.isChecked()
            elif activity == True:
                log = self.radioButton_act_GLS_period.isChecked()
            elif MLP == True:
                log = self.radioButton_RV_MLP_period.isChecked()
                N_peaks = int(self.N_MLP_peak_to_point.value())
            elif DFT == True:
                log = self.radioButton_RV_WF_period.isChecked()
                N_peaks = int(self.N_window_peak_to_point.value())

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

            text_arrow = pg.TextItem("", anchor=(0.5,1.9))

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
 

            plot_wg2.addItem(arrow)#, ignoreBounds=True)   
            plot_wg2.addItem(text_arrow, ignoreBounds=True)  

        #  plot_wg2.autoRange()
      #  plot_wg2.enableAutoRange(enable=True)
       # plot_wg2.getViewBox().setAutoVisible(y=True)

        #vb = plot_wg2.getViewBox()   
      #  viewrange = vb.viewRange()
        #plot_wg2.autoRange( padding=0.02)
        #plot_wg2.getViewBox().setAutoVisible(y=True)
        #plot_wg2.getViewBox().enableAutoRange(axis='x',enable=0.49)
        #plot_wg2.setYRange(-0.001, max(y_peaks)+0.3*max(y_peaks), padding=0.001)
        #plot_wg2.setLimits(yMax=max(y_peaks)+0.3*max(y_peaks))
       # plot_wg2.getViewBox().setXRange(-0.001, max(y_peaks)+0.3*max(y_peaks), padding=0.01)
        #plot_wg2.setXRange(min(model_time_phase), max(model_time_phase), padding=0.002)
        #plot_wg2.setYRange(-0.001, max(y_peaks)+0.3*max(y_peaks), padding=0.01)
        if len(y_peaks) != 0:
            plot_wg2.setRange(yRange=(-0.001, max(y_peaks)+0.3*max(y_peaks)), padding=0.01, update=True, disableAutoRange=True)
        #plot_wg2.getViewBox().updateAutoRange()
       
        #plot_wg2.enableAutoRange(plot_wg2.getViewBox().YAxis, y=True)
        #plot_wg2.ViewBox().updateAutoRange()
        #plot_wg2.setAutoVisible(y=True,x=True)   
        #plot_wg2.enableAutoRange()
        return
        #plot_wg2.scene()
        #plot_wg2.setAutoPan()
       # plot_wg2.getViewBox().setAutoVisible(y=True)

######################## RV plots ######################################

    def run_gls(self):
        global fit

      #  omega = 1/ np.logspace(np.log10(self.gls_min_period.value()), np.log10(self.gls_max_period.value()), num=int(self.gls_n_omega.value()))
        ind_norm = self.gls_norm_combo.currentIndex()

        if self.gls_incl_jitter.isChecked():
            error_list = self.add_jitter(fit.fit_results.rv_model.rv_err, fit.fit_results.idset)
        else:
            error_list = fit.fit_results.rv_model.rv_err


        if len(fit.fit_results.rv_model.jd) > 5:      
            RV_per = gls.Gls((fit.fit_results.rv_model.jd, fit.fit_results.rv_model.rvs, error_list), 
            #fast=True,  verbose=False, norm=self.norms[ind_norm],ofac=self.gls_ofac.value(), fbeg=omega[-1], fend=omega[0],)
            fast=True,  verbose=False, norm=self.norms[ind_norm],ofac=self.gls_ofac.value(), fbeg=1/self.gls_max_period.value(), fend=1/self.gls_min_period.value())            
            fit.gls = RV_per
        else:
            return

        self.update_RV_GLS_plots()
        self.update_WF_plots()


    def run_gls_o_c(self):
        global fit
                        
 
        if fit.doGP == True and self.gls_o_c_GP.isChecked() == True:
            data_o_c = fit.fit_results.rv_model.o_c - fit.gp_model_data[0] 
        else:
            data_o_c = fit.fit_results.rv_model.o_c  

        
        #omega = 1/ np.logspace(np.log10(self.gls_min_period.value()), np.log10(self.gls_max_period.value()), num=int(self.gls_n_omega.value()))
        ind_norm = self.gls_norm_combo.currentIndex()

        if self.gls_o_c_incl_jitter.isChecked():
            error_list = self.add_jitter(fit.fit_results.rv_model.rv_err, fit.fit_results.idset)
        else:
            error_list = fit.fit_results.rv_model.rv_err
 
        if len(fit.fit_results.rv_model.jd) > 5:
            RV_per_res = gls.Gls((fit.fit_results.rv_model.jd, data_o_c, error_list), 
            fast=True,  verbose=False, norm= self.norms[ind_norm],ofac=self.gls_ofac.value(), fbeg=1/self.gls_max_period.value(), fend=1/self.gls_min_period.value())            

            fit.gls_o_c = RV_per_res
        else:
            return

        self.update_RV_o_c_GLS_plots()  


#    def run_WF(self, RV = True):
#        global fit
#
#
#        omega = 1/ np.logspace(np.log10(self.gls_min_period.value()), np.log10(self.gls_max_period.value()),  num=int(self.gls_n_omega.value()))
# 
#
#        if len(fit.fit_results.rv_model.jd) > 5:
#            ######################## DFT (Window) ##############################
#            WF_power = []
#            for omi in 2*np.pi*omega: 
#                phase = (fit.fit_results.rv_model.jd-fit.fit_results.rv_model.jd[0]) * omi
#                WC = np.sum(np.cos(phase))
#                WS = np.sum(np.sin(phase))
#                WF_power.append((WC**2 + WS**2)/len(fit.fit_results.rv_model.jd)**2) 
#
#            fit.RV_WF_power = np.array(WF_power)
#            fit.RV_WF_omega = np.array(omega)
# 
#        else:
#            return
#
#        self.update_WF_plots()


    def init_gls_norm_combo(self):    
        global fit

        self.norms = ['ZK',  'HorneBaliunas', 'Cumming', 'wrms', 'chisq', 'lnL', 'dlnL']
        #'Scargle',
        for i in range(len(self.norms)):
            self.gls_norm_combo.addItem('%s'%(self.norms[i]),i+1)

    def get_RV_GLS_plot_color(self):
        global fit

        colorz = self.colorDialog.getColor(options=QtWidgets.QColorDialog.DontUseNativeDialog)
        fit.gls_colors[0]=dill.copy(colorz.name())   

        self.update_RV_GLS_plots() 


    def get_RV_o_c_GLS_plot_color(self):
        global fit
        
        colorz = self.colorDialog.getColor(options=QtWidgets.QColorDialog.DontUseNativeDialog)
        fit.gls_colors[1]=dill.copy(colorz.name())   
         
        self.update_RV_o_c_GLS_plots()  
        
    def get_RV_GLS_alias_color(self):
        global fit

        colorz = self.colorDialog.getColor(options=QtWidgets.QColorDialog.DontUseNativeDialog)
        colors_GLS_alias[0]=dill.copy(colorz.name()) 

        self.update_RV_GLS_plots() 
        self.update_RV_o_c_GLS_plots() 

    def get_RV_MLP_alias_color(self):
        global fit

        colorz = self.colorDialog.getColor(options=QtWidgets.QColorDialog.DontUseNativeDialog)
        colors_MLP_alias[0]=dill.copy(colorz.name()) 

        self.update_RV_MLP_plots() 

    def get_RV_jitter_color(self):
        global fit

        colorz = self.colorDialog.getColor(options=QtWidgets.QColorDialog.DontUseNativeDialog)
        colors_RV_jitter[0]=dill.copy(colorz.name())   

        self.update_RV_plots() 
        self.update_extra_plots()


    def update_RV_GLS_plots(self):
        global fit, p7 
 
        p7.plot(clear=True,)

        self.colors_gls.setStyleSheet("color: %s;"%dill.copy(fit.gls_colors[0]))
        self.colors_alias_gls.setStyleSheet("color: %s;"%colors_GLS_alias[0])

                          
        power_levels = np.array([self.gls_fap1.value(),self.gls_fap2.value(),self.gls_fap3.value()])
        gls_model_width = float(self.gls_model_width.value())
    
        if len(fit.fit_results.rv_model.jd) > 5:

            ######################## GLS ##############################
            if self.radioButton_RV_GLS_period.isChecked():
                p7.setLogMode(True,False)        
                p7.plot(1/fit.gls.freq, fit.gls.power,pen={'color': dill.copy(fit.gls_colors[0]), 'width': gls_model_width},symbol=None ) 

                p7.setLabel('bottom', 'period [d]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})    
                
            else:
                p7.setLogMode(False,False)        
                p7.plot(fit.gls.freq, fit.gls.power,pen={'color': dill.copy(fit.gls_colors[0]), 'width': self.gls_model_width.value()},symbol=None )                
                p7.setLabel('bottom', 'frequency [1/d]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
                
                
            if fit.gls.norm == 'ZK':
                [p7.addLine(x=None, y=fap, pen=pg.mkPen('k', width=0.8, style=QtCore.Qt.PenStyle.DotLine)) for ii,fap in enumerate(fit.gls.powerLevel(np.array(power_levels)))]
 
            text_peaks, pos_peaks = self.identify_power_peaks(1/fit.gls.freq, fit.gls.power, power_level = power_levels, sig_level = fit.gls.powerLevel(np.array(power_levels)) )   

            self.label_peaks(p7, pos_peaks, GLS = True, o_c = False)

            self.RV_periodogram_print_info.clicked.connect(lambda: self.print_info_for_object(
            fit.gls.info(stdout=False) + text_peaks   ))   

        if self.gls_cross_hair.isChecked():
            self.cross_hair(p7,log=self.radioButton_RV_GLS_period.isChecked(), alias=[self.show_alias_GLS.isChecked(), self.alias_days_gls.value(), colors_GLS_alias[0]])    
 
    
       # fit.pgg = p7.getPlotItem()
 
 
    def update_RV_o_c_GLS_plots(self):
        global fit,  p8  
 
        p8.plot(clear=True,)  

        self.colors_gls_o_c.setStyleSheet("color: %s;"%dill.copy(fit.gls_colors[1])) 
        gls_o_c_model_width = float(self.gls_o_c_model_width.value())
        
        power_levels = np.array([self.gls_fap1.value(),self.gls_fap2.value(),self.gls_fap3.value()])

        if len(fit.fit_results.rv_model.jd) > 5:
 

            ######################## GLS o-c ##############################
            if self.radioButton_RV_o_c_GLS_period.isChecked():
                p8.setLogMode(True,False)        
                p8.plot(1/fit.gls_o_c.freq, fit.gls_o_c.power, pen={'color': dill.copy(fit.gls_colors[1]), 'width': gls_o_c_model_width},symbol=None ) 
                p8.setLabel('bottom', 'period [d]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})
 
            else:
                p8.setLogMode(False,False)        
                p8.plot(fit.gls_o_c.freq, fit.gls_o_c.power, pen={'color': dill.copy(fit.gls_colors[1]), 'width': gls_o_c_model_width},symbol=None )   
                p8.setLabel('bottom', 'frequency [1/d]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})                
                
            if fit.gls_o_c.norm == 'ZK':
                [p8.addLine(x=None, y=fap, pen=pg.mkPen('k', width=0.8, style=QtCore.Qt.PenStyle.DotLine)) for ii,fap in enumerate(fit.gls_o_c.powerLevel(np.array(power_levels)))]            

            text_peaks, pos_peaks = self.identify_power_peaks(1/fit.gls_o_c.freq, fit.gls_o_c.power, power_level = power_levels, sig_level = fit.gls_o_c.powerLevel(np.array(power_levels)) )

            self.label_peaks(p8, pos_peaks, GLS=True, o_c = True)
 
            self.RV_res_periodogram_print_info.clicked.connect(lambda: self.print_info_for_object(fit.gls_o_c.info(stdout=False)+ text_peaks  )  )      
 
        if self.gls_o_c_cross_hair.isChecked():
            self.cross_hair(p8,log=self.radioButton_RV_o_c_GLS_period.isChecked(), alias=[self.show_alias_GLS.isChecked(), self.alias_days_gls.value(), colors_GLS_alias[0]])    
    
    def update_WF_plots(self):
        global fit, p12  
 
        p12.plot(clear=True,) 
        p12.setLogMode(True,False)
                        

        ######################## GLS o-c ##############################
        if self.radioButton_RV_WF_period.isChecked():
            p12.setLogMode(True,False)        
            #p12.plot(1/fit.RV_WF_omega, fit.RV_WF_power,pen='k',symbol=None , viewRect=True, enableAutoRange=True)   
            p12.plot(1/fit.gls.WF_omega, fit.gls.WF_power,pen='k',symbol=None , viewRect=True, enableAutoRange=True)   
            p12.setLabel('bottom', 'period [d]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})
        else:
            p12.setLogMode(False,False)        
           # p12.plot(fit.RV_WF_omega, fit.RV_WF_power,pen='k',symbol=None,  viewRect=True, enableAutoRange=True)  
            p12.plot(fit.gls.WF_omega, fit.gls.WF_power,pen='k',symbol=None,  viewRect=True, enableAutoRange=True)   
            p12.setLabel('bottom', 'frequency [1/d]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})

        text_peaks, pos_peaks = self.identify_power_peaks(1/fit.gls.WF_omega, fit.gls.WF_power)

        self.label_peaks(p12, pos_peaks, GLS = True, DFT = True)

        self.WF_print_info.clicked.connect(lambda: self.print_info_for_object(text_peaks))


    def update_RV_plot(self):
        global fit, p1

        p1.plot(clear=True,)
 
        p1.addLine(x=None, y=0,   pen=pg.mkPen('#ff9933', width=0.8))
  
        if fit.doGP == True:
            y_model = fit.fit_results.model + fit.gp_model_curve[0]
        else:
            y_model = fit.fit_results.model 
            
            
        model_curve = p1.plot(fit.fit_results.model_jd,y_model, 
        pen={'color': dill.copy(fit.rvs_colors[-1]), 'width': self.rv_model_width.value()},enableAutoRange=True, #symbolPen={'color': 0.5, 'width': 0.1}, symbolSize=1,symbol='o',
        viewRect=True, labels =  {'left':'RV', 'bottom':'JD'}) 
        
        model_curve.setZValue(self.RV_model_z.value()) 
        
        
        if fit.doGP == True:
            pfill = pg.FillBetweenItem(p1.plot(fit.fit_results.model_jd, fit.fit_results.model + fit.gp_model_curve[0]+fit.gp_model_curve[2]), 
                                       p1.plot(fit.fit_results.model_jd, fit.fit_results.model + fit.gp_model_curve[0]-fit.gp_model_curve[2]), 
                                       brush = pg.mkColor(244,140,66,128))
            p1.addItem(pfill) 
 

        if self.jitter_to_plots.isChecked() and not self.split_jitter.isChecked():
            error_list = self.add_jitter(fit.fit_results.rv_model.rv_err, fit.fit_results.idset)
        elif self.jitter_to_plots.isChecked() and self.split_jitter.isChecked():
            error_list = fit.fit_results.rv_model.rv_err
            error_list2 = self.add_jitter(fit.fit_results.rv_model.rv_err, fit.fit_results.idset)
        else:
            error_list = fit.fit_results.rv_model.rv_err


        #for i in range(max(fit.fit_results.idset)+1):

        for i in range(len(fit.rv_data_sets)):
            if len(fit.rv_data_sets[i]) == 0:
                continue         

            p1.plot(fit.fit_results.rv_model.jd[fit.fit_results.idset==i],fit.fit_results.rv_model.rvs[fit.fit_results.idset==i], 
            pen=None, #{'color': colors[i], 'width': 1.1},
            symbol=dill.copy(fit.pyqt_symbols_rvs[i]),
            symbolPen={'color': dill.copy(fit.rvs_colors[i])+"%02x"%int(dill.copy(fit.pyqt_color_alpha_rvs[i]))
, 'width': 1.1},
            symbolSize=dill.copy(fit.pyqt_symbols_size_rvs[i]),enableAutoRange=True,viewRect=True,
            symbolBrush=dill.copy(fit.rvs_colors[i])+"%02x"%int(dill.copy(fit.pyqt_color_alpha_rvs[i]))
            )

            err1 = pg.ErrorBarItem(x=fit.fit_results.rv_model.jd[fit.fit_results.idset==i], 
                                   y=fit.fit_results.rv_model.rvs[fit.fit_results.idset==i],symbol='o', 
            #height=error_list[fit.filelist.idset==i],
            top=error_list[fit.fit_results.idset==i],
            bottom=error_list[fit.fit_results.idset==i],
            beam=0.0, pen=dill.copy(fit.rvs_colors[i])+"%02x"%int(dill.copy(fit.pyqt_color_alpha_rvs[i])))  

            p1.addItem(err1)

            if self.jitter_to_plots.isChecked() and self.split_jitter.isChecked():

                err1a = pg.ErrorBarItem(x=fit.fit_results.rv_model.jd[fit.fit_results.idset==i], 
                                       y=fit.fit_results.rv_model.rvs[fit.fit_results.idset==i],symbol='o', 
                #height=error_list[fit.filelist.idset==i],
                top=error_list2[fit.fit_results.idset==i],
                bottom=error_list2[fit.fit_results.idset==i],
                beam=0.0, pen=colors_RV_jitter[0])  
                err1a.setZValue(-10)
                p1.addItem(err1a)

        if self.RV_plot_autorange.isChecked():
            p1.autoRange() #padding=0

        if self.RV_plot_cross_hair.isChecked():
            self.cross_hair(p1,log=False)  
            
            
    def update_RV_plot_o_c(self):
        global fit, p2

 
        p2.addLine(x=None, y=0, pen=pg.mkPen('#ff9933', width=0.8))

        if len(fit.fit_results.idset)==0:
            return

        if self.jitter_to_plots.isChecked() and not self.split_jitter.isChecked():
            error_list = self.add_jitter(fit.fit_results.rv_model.rv_err, fit.fit_results.idset)
        elif self.jitter_to_plots.isChecked() and self.split_jitter.isChecked():
            error_list = fit.fit_results.rv_model.rv_err
            error_list2 = self.add_jitter(fit.fit_results.rv_model.rv_err, fit.fit_results.idset)
        else:
            error_list = fit.fit_results.rv_model.rv_err
        
        if fit.doGP == True and self.plot_RV_GP_model.isChecked() == False:
            pfill_o_c = pg.FillBetweenItem(p2.plot(fit.fit_results.model_jd, fit.gp_model_curve[0]+fit.gp_model_curve[2]), 
                                           p2.plot(fit.fit_results.model_jd, fit.gp_model_curve[0]-fit.gp_model_curve[2]), 
                                           brush = pg.mkColor(244,140,66,128))
            p2.addItem(pfill_o_c)
            y_model_o_c = fit.gp_model_curve[0]
            data_o_c = fit.fit_results.rv_model.o_c 
        elif fit.doGP == True and self.plot_RV_GP_model.isChecked() ==  True and len(fit.fit_results.rv_model.o_c)==len(fit.gp_model_data[0]):
            data_o_c = fit.fit_results.rv_model.o_c - fit.gp_model_data[0]
            y_model_o_c = np.zeros(len(fit.fit_results.model))
        else:
            data_o_c = fit.fit_results.rv_model.o_c
            y_model_o_c = np.zeros(len(fit.fit_results.model))



        model_curve_o_c = p2.plot(fit.fit_results.model_jd,y_model_o_c, 
        pen={'color':  dill.copy(fit.rvs_colors[-1]), 'width': self.rv_model_width.value()},enableAutoRange=True, #symbolPen={'color': 0.5, 'width': 0.1}, symbolSize=1,symbol='o',
        viewRect=True, labels =  {'left':'RV', 'bottom':'JD'}) 
        
        
        model_curve_o_c.setZValue(self.RV_model_z.value()) 

 
            
       # for i in range(max(fit.fit_results.idset)+1):

        for i in range(len(fit.rv_data_sets)):
            if len(fit.rv_data_sets[i]) == 0:
                continue         

            p2.plot(fit.fit_results.rv_model.jd[fit.fit_results.idset==i],data_o_c[fit.fit_results.idset==i], 
            pen=None, #{'color': colors[i], 'width': 1.1},
            symbol=dill.copy(fit.pyqt_symbols_rvs[i]),
            symbolPen={'color': dill.copy(fit.rvs_colors[i])+"%02x"%int(dill.copy(fit.pyqt_color_alpha_rvs[i])), 'width': 1.1},
            symbolSize=dill.copy(fit.pyqt_symbols_size_rvs[i]),enableAutoRange=True,viewRect=True,
            symbolBrush=dill.copy(fit.rvs_colors[i])+"%02x"%int(dill.copy(fit.pyqt_color_alpha_rvs[i]))
            )
            err2 = pg.ErrorBarItem(x=fit.fit_results.rv_model.jd[fit.fit_results.idset==i], 
                                   y=data_o_c[fit.fit_results.idset==i],symbol='o', 
            #height=error_list[fit.filelist.idset==i],
            top=error_list[fit.fit_results.idset==i],
            bottom=error_list[fit.fit_results.idset==i],
            beam=0.0, pen=dill.copy(fit.rvs_colors[i])+"%02x"%int(dill.copy(fit.pyqt_color_alpha_rvs[i])))  

            p2.addItem(err2)

            if self.jitter_to_plots.isChecked() and self.split_jitter.isChecked():

                err2a = pg.ErrorBarItem(x=fit.fit_results.rv_model.jd[fit.fit_results.idset==i], 
                                       y=data_o_c[fit.fit_results.idset==i],symbol='o',
                #height=error_list[fit.filelist.idset==i],
                top=error_list2[fit.fit_results.idset==i],
                bottom=error_list2[fit.fit_results.idset==i],
                beam=0.0, pen='#000000')
                err2a.setZValue(-10)
                p2.addItem(err2a)

        if self.RV_o_c_plot_cross_hair.isChecked():
            self.cross_hair(p2,log=False)       
            
        if self.RV_plot_autorange.isChecked():
            p2.autoRange()


#    def control_RV_plot_add_o_c(self):

#        if self.RV_plot_add_o_c.isChecked():
#            self.update_RV_plot_with_o_c()
#        else:
#            p1.removeItem(p00)
#            p1.removeItem(p01) 
#            self.update_RV_plot()

    def update_RV_plots(self):
        global fit, p1,p2,p00,p01

        p1.plot(clear=True,)
        p2.plot(clear=True,)

        self.check_RV_symbol_sizes()
        self.jitter_color_button.setStyleSheet("color: %s;"%colors_RV_jitter[0])


        #print(self.RV_plot_add_o_c.isChecked())

        if self.RV_plot_add_o_c.isChecked():
            self.update_RV_plot_with_o_c()            
        else:
            #if hasattr(p1, 'p00'):
           # p1.removeItem(-1)
            #if hasattr(p1, 'p01'):
            try:
                p1.scene().removeItem(p00)
                p1.scene().removeItem(p01)
            #    p1.centralWidget = None
            except:
                pass
            #p1  = self.graphicsView_timeseries_RV
           # p1 = pz
            self.update_RV_plot()

        self.update_RV_plot_o_c()






    def update_RV_plot_with_o_c(self):
        global fit, p1,p2,p00,p01,legend_RV
 

        #p1.scene().removeItem(-1)


        if not self.hold_old_plot_RV.isChecked():    
            p00.plot(clear=True,)
            p01.plot(clear=True,)
        else:
            p00.plot(clear=False,)
            p01.plot(clear=False,)

        if self.RV_legend.isChecked()==True:
            legend_RV.clear()
            legend_RV.setVisible(True)
        else:
            legend_RV.setVisible(False)

 
        #fit.p1 = p1
        self.check_RV_symbol_sizes()
        self.jitter_color_button.setStyleSheet("color: %s;"%colors_RV_jitter[0])

        if len(fit.fit_results.idset)==0:
            return

       # fit.p00 = p00

        if fit.doGP == True:
            #rv.get_RV_gps_model(self) 
            y_model = fit.fit_results.model + fit.gp_model_curve[0]
        else:
            y_model = fit.fit_results.model 
  

        p00.addLine(x=None, y=0,   pen=pg.mkPen('#ff9933', width=0.8))
          
            
        model_curve = p00.plot(fit.fit_results.model_jd,y_model, 
        pen={'color': dill.copy(fit.rvs_colors[-1]), 'width': self.rv_model_width.value()},enableAutoRange=True, #symbolPen={'color': 0.5, 'width': 0.1}, symbolSize=1,symbol='o',
        viewRect=True, labels =  {'left':'RV', 'bottom':'JD'}) 
        
        model_curve.setZValue(self.RV_model_z.value()) 
        
        
        if fit.doGP == True:
            pfill = pg.FillBetweenItem(p00.plot(fit.fit_results.model_jd, fit.fit_results.model + fit.gp_model_curve[0]+fit.gp_model_curve[2]), 
                                       p00.plot(fit.fit_results.model_jd, fit.fit_results.model + fit.gp_model_curve[0]-fit.gp_model_curve[2]), 
                                       brush = pg.mkColor(244,140,66,128))
            p00.addItem(pfill) 
            
        if self.jitter_to_plots.isChecked() and not self.split_jitter.isChecked():
            error_list = self.add_jitter(fit.fit_results.rv_model.rv_err, fit.fit_results.idset)
        elif self.jitter_to_plots.isChecked() and self.split_jitter.isChecked():
            error_list = fit.fit_results.rv_model.rv_err
            error_list2 = self.add_jitter(fit.fit_results.rv_model.rv_err, fit.fit_results.idset)
        else:
            error_list = fit.fit_results.rv_model.rv_err
 

#        for i in range(max(fit.fit_results.idset)+1):
        for i in range(len(fit.rv_data_sets)):
            if len(fit.rv_data_sets[i]) == 0:
                continue
            p00.plot(fit.fit_results.rv_model.jd[fit.fit_results.idset==i],fit.fit_results.rv_model.rvs[fit.fit_results.idset==i], 
            pen=None, #{'color': colors[i], 'width': 1.1},
            symbol=dill.copy(fit.pyqt_symbols_rvs[i]),
            symbolPen={'color': dill.copy(fit.rvs_colors[i])+"%02x"%int(dill.copy(fit.pyqt_color_alpha_rvs[i]))
, 'width': 1.1},
            symbolSize=dill.copy(fit.pyqt_symbols_size_rvs[i]),enableAutoRange=True,viewRect=True,
            symbolBrush=dill.copy(fit.rvs_colors[i])+"%02x"%int(dill.copy(fit.pyqt_color_alpha_rvs[i])),
#            name=fit.filelist.files[i].name)
            name=fit.rv_data_sets[i][4])

            err1 = pg.ErrorBarItem(x=fit.fit_results.rv_model.jd[fit.fit_results.idset==i], 
                                   y=fit.fit_results.rv_model.rvs[fit.fit_results.idset==i],symbol='o', 
            #height=error_list[fit.filelist.idset==i],
            top=error_list[fit.fit_results.idset==i],
            bottom=error_list[fit.fit_results.idset==i],
            beam=0.0, pen=dill.copy(fit.rvs_colors[i])+"%02x"%int(dill.copy(fit.pyqt_color_alpha_rvs[i])))  

            p00.addItem(err1)

            if self.jitter_to_plots.isChecked() and self.split_jitter.isChecked():

                err1a = pg.ErrorBarItem(x=fit.fit_results.rv_model.jd[fit.fit_results.idset==i], 
                                       y=fit.fit_results.rv_model.rvs[fit.fit_results.idset==i],symbol='o', 
                #height=error_list[fit.filelist.idset==i],
                top=error_list2[fit.fit_results.idset==i],
                bottom=error_list2[fit.fit_results.idset==i],
                beam=0.0, pen=colors_RV_jitter[0])  
                err1a.setZValue(-10)
                p00.addItem(err1a)

        if self.RV_plot_autorange.isChecked():
            p00.autoRange() #padding=0

        if self.RV_plot_cross_hair.isChecked():
            self.cross_hair(p00,log=False)  
            

        p01.addLine(x=None, y=0, pen=pg.mkPen('#ff9933', width=0.8))

        
        if fit.doGP == True and self.plot_RV_GP_model.isChecked() == False:
            pfill_o_c = pg.FillBetweenItem(p01.plot(fit.fit_results.model_jd, fit.gp_model_curve[0]+fit.gp_model_curve[2]), 
                                           p01.plot(fit.fit_results.model_jd, fit.gp_model_curve[0]-fit.gp_model_curve[2]), 
                                           brush = pg.mkColor(244,140,66,128))
            p01.addItem(pfill_o_c)
            y_model_o_c = fit.gp_model_curve[0]
            data_o_c = fit.fit_results.rv_model.o_c 
        elif fit.doGP == True and self.plot_RV_GP_model.isChecked() ==  True and len(fit.fit_results.rv_model.o_c)==len(fit.gp_model_data[0]):
            data_o_c = fit.fit_results.rv_model.o_c - fit.gp_model_data[0]
            y_model_o_c = np.zeros(len(y_model))
        else:
            data_o_c = fit.fit_results.rv_model.o_c
            y_model_o_c = np.zeros(len(y_model))


        if qso_mode:
            label = "flux [mag.]"
        else:
            label = "RV [m/s]" 

        model_curve_o_c = p01.plot(fit.fit_results.model_jd,y_model_o_c, 
        pen={'color':  dill.copy(fit.rvs_colors[-1]), 'width': self.rv_model_width.value()},enableAutoRange=True, #symbolPen={'color': 0.5, 'width': 0.1}, symbolSize=1,symbol='o',
        viewRect=True, labels =  {'left':'%s'%label, 'bottom':'JD'}) 
        
        
        model_curve_o_c.setZValue(self.RV_model_z.value()) 

 
        for i in range(len(fit.rv_data_sets)):
            if len(fit.rv_data_sets[i]) == 0:
                continue            
        #for i in range(max(fit.fit_results.idset)+1):
            p01.plot(fit.fit_results.rv_model.jd[fit.fit_results.idset==i],data_o_c[fit.fit_results.idset==i], 
            pen=None, #{'color': colors[i], 'width': 1.1},
            symbol=dill.copy(fit.pyqt_symbols_rvs[i]),
            symbolPen={'color': dill.copy(fit.rvs_colors[i])+"%02x"%int(dill.copy(fit.pyqt_color_alpha_rvs[i])), 'width': 1.1},
            symbolSize=dill.copy(fit.pyqt_symbols_size_rvs[i]),enableAutoRange=True,viewRect=True,
            symbolBrush=dill.copy(fit.rvs_colors[i])+"%02x"%int(dill.copy(fit.pyqt_color_alpha_rvs[i]))
            )
            err2 = pg.ErrorBarItem(x=fit.fit_results.rv_model.jd[fit.fit_results.idset==i], 
                                   y=data_o_c[fit.fit_results.idset==i],symbol='o', 
            #height=error_list[fit.filelist.idset==i],
            top=error_list[fit.fit_results.idset==i],
            bottom=error_list[fit.fit_results.idset==i],
            beam=0.0, pen=dill.copy(fit.rvs_colors[i])+"%02x"%int(dill.copy(fit.pyqt_color_alpha_rvs[i])))  

            p01.addItem(err2)

            if self.jitter_to_plots.isChecked() and self.split_jitter.isChecked():

                err2a = pg.ErrorBarItem(x=fit.fit_results.rv_model.jd[fit.fit_results.idset==i], 
                                       y=data_o_c[fit.fit_results.idset==i],symbol='o',
                #height=error_list[fit.filelist.idset==i],
                top=error_list2[fit.fit_results.idset==i],
                bottom=error_list2[fit.fit_results.idset==i],
                beam=0.0, pen='#000000')
                err2a.setZValue(-10)
                p01.addItem(err2a)

        p01.setLabel('bottom', '%s'%p01.getAxis("bottom").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
        p01.setLabel('left', '%s'%p01.getAxis("left").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})


        if self.RV_o_c_plot_cross_hair.isChecked():
            self.cross_hair(p01,log=False)       
            
        if self.RV_plot_autorange.isChecked():
            p01.autoRange()
            


    def update_RV_plots_old(self):
        global fit, p1,p2

        p1.plot(clear=True,)
        p2.plot(clear=True,)
 
        self.check_RV_symbol_sizes()
        self.jitter_color_button.setStyleSheet("color: %s;"%colors_RV_jitter[0])

        if len(fit.fit_results.idset)==0:
            return
 

        p1.addLine(x=None, y=0,   pen=pg.mkPen('#ff9933', width=0.8))
  
        if fit.doGP == True:
            #rv.get_RV_gps_model(self) 
            y_model = fit.fit_results.model + fit.gp_model_curve[0]
        else:
            y_model = fit.fit_results.model 
            
            
        model_curve = p1.plot(fit.fit_results.model_jd,y_model, 
        pen={'color': dill.copy(fit.rvs_colors[-1]), 'width': self.rv_model_width.value()},enableAutoRange=True, #symbolPen={'color': 0.5, 'width': 0.1}, symbolSize=1,symbol='o',
        viewRect=True, labels =  {'left':'RV', 'bottom':'JD'}) 
        
        model_curve.setZValue(self.RV_model_z.value()) 
        
        
        if fit.doGP == True:
            pfill = pg.FillBetweenItem(p1.plot(fit.fit_results.model_jd, fit.fit_results.model + fit.gp_model_curve[0]+fit.gp_model_curve[2]), 
                                       p1.plot(fit.fit_results.model_jd, fit.fit_results.model + fit.gp_model_curve[0]-fit.gp_model_curve[2]), 
                                       brush = pg.mkColor(244,140,66,128))
            p1.addItem(pfill) 
            
            
#        if len(fit.fit_results.rv_model.rv_err) != len(fit.filelist.idset) or len(fit.fit_results.rv_model.rv_err) ==0 :
#            return

        if self.jitter_to_plots.isChecked() and not self.split_jitter.isChecked():
            error_list = self.add_jitter(fit.fit_results.rv_model.rv_err, fit.fit_results.idset)
        elif self.jitter_to_plots.isChecked() and self.split_jitter.isChecked():
            error_list = fit.fit_results.rv_model.rv_err
            error_list2 = self.add_jitter(fit.fit_results.rv_model.rv_err, fit.fit_results.idset)
        else:
            error_list = fit.fit_results.rv_model.rv_err


        for i in range(max(fit.fit_results.idset)+1):
            p1.plot(fit.fit_results.rv_model.jd[fit.fit_results.idset==i],fit.fit_results.rv_model.rvs[fit.fit_results.idset==i], 
            pen=None, #{'color': colors[i], 'width': 1.1},
            symbol=dill.copy(fit.pyqt_symbols_rvs[i]),
            symbolPen={'color': dill.copy(fit.rvs_colors[i])+"%02x"%int(dill.copy(fit.pyqt_color_alpha_rvs[i]))
, 'width': 1.1},
            symbolSize=dill.copy(fit.pyqt_symbols_size_rvs[i]),enableAutoRange=True,viewRect=True,
            symbolBrush=dill.copy(fit.rvs_colors[i])+"%02x"%int(dill.copy(fit.pyqt_color_alpha_rvs[i]))
            )

            err1 = pg.ErrorBarItem(x=fit.fit_results.rv_model.jd[fit.fit_results.idset==i], 
                                   y=fit.fit_results.rv_model.rvs[fit.fit_results.idset==i],symbol='o', 
            #height=error_list[fit.filelist.idset==i],
            top=error_list[fit.fit_results.idset==i],
            bottom=error_list[fit.fit_results.idset==i],
            beam=0.0, pen=dill.copy(fit.rvs_colors[i])+"%02x"%int(dill.copy(fit.pyqt_color_alpha_rvs[i])))  

            p1.addItem(err1)

            if self.jitter_to_plots.isChecked() and self.split_jitter.isChecked():

                err1a = pg.ErrorBarItem(x=fit.fit_results.rv_model.jd[fit.fit_results.idset==i], 
                                       y=fit.fit_results.rv_model.rvs[fit.fit_results.idset==i],symbol='o', 
                #height=error_list[fit.filelist.idset==i],
                top=error_list2[fit.fit_results.idset==i],
                bottom=error_list2[fit.fit_results.idset==i],
                beam=0.0, pen=colors_RV_jitter[0])  
                err1a.setZValue(-10)
                p0.addItem(err1a)

        if self.RV_plot_autorange.isChecked():
            p1.autoRange() #padding=0

        if self.RV_plot_cross_hair.isChecked():
            self.cross_hair(p1,log=False)  
            
            
            
        #if fit.doGP == True:
        #    y_model_o_c = fit.gp_model_curve[0]
       # else:
       #     y_model_o_c = np.zeros(len(y_model))
        p2.addLine(x=None, y=0, pen=pg.mkPen('#ff9933', width=0.8))

        
        if fit.doGP == True and self.plot_RV_GP_model.isChecked() == False:
            pfill_o_c = pg.FillBetweenItem(p2.plot(fit.fit_results.model_jd, fit.gp_model_curve[0]+fit.gp_model_curve[2]), 
                                           p2.plot(fit.fit_results.model_jd, fit.gp_model_curve[0]-fit.gp_model_curve[2]), 
                                           brush = pg.mkColor(244,140,66,128))
            p00.addItem(pfill_o_c)
            y_model_o_c = fit.gp_model_curve[0]
            data_o_c = fit.fit_results.rv_model.o_c 
        elif fit.doGP == True and self.plot_RV_GP_model.isChecked() ==  True and len(fit.fit_results.rv_model.o_c)==len(fit.gp_model_data[0]):
            data_o_c = fit.fit_results.rv_model.o_c - fit.gp_model_data[0]
            y_model_o_c = np.zeros(len(y_model))
        else:
            data_o_c = fit.fit_results.rv_model.o_c
            y_model_o_c = np.zeros(len(y_model))



        model_curve_o_c = p00.plot(fit.fit_results.model_jd,y_model_o_c, 
        pen={'color':  dill.copy(fit.rvs_colors[-1]), 'width': self.rv_model_width.value()},enableAutoRange=True, #symbolPen={'color': 0.5, 'width': 0.1}, symbolSize=1,symbol='o',
        viewRect=True, labels =  {'left':'RV', 'bottom':'JD'}) 
        
        
        model_curve_o_c.setZValue(self.RV_model_z.value()) 

 
            
        for i in range(max(fit.fit_results.idset)+1):
            p2.plot(fit.fit_results.rv_model.jd[fit.fit_results.idset==i],data_o_c[fit.fit_results.idset==i], 
            pen=None, #{'color': colors[i], 'width': 1.1},
            symbol=dill.copy(fit.pyqt_symbols_rvs[i]),
            symbolPen={'color': dill.copy(fit.rvs_colors[i])+"%02x"%int(dill.copy(fit.pyqt_color_alpha_rvs[i])), 'width': 1.1},
            symbolSize=dill.copy(fit.pyqt_symbols_size_rvs[i]),enableAutoRange=True,viewRect=True,
            symbolBrush=dill.copy(fit.rvs_colors[i])+"%02x"%int(dill.copy(fit.pyqt_color_alpha_rvs[i]))
            )
            err2 = pg.ErrorBarItem(x=fit.fit_results.rv_model.jd[fit.fit_results.idset==i], 
                                   y=data_o_c[fit.fit_results.idset==i],symbol='o', 
            #height=error_list[fit.filelist.idset==i],
            top=error_list[fit.fit_results.idset==i],
            bottom=error_list[fit.fit_results.idset==i],
            beam=0.0, pen=dill.copy(fit.rvs_colors[i])+"%02x"%int(dill.copy(fit.pyqt_color_alpha_rvs[i])))  

            p00.addItem(err2)

            if self.jitter_to_plots.isChecked() and self.split_jitter.isChecked():

                err2a = pg.ErrorBarItem(x=fit.fit_results.rv_model.jd[fit.fit_results.idset==i], 
                                       y=data_o_c[fit.fit_results.idset==i],symbol='o',
                #height=error_list[fit.filelist.idset==i],
                top=error_list2[fit.fit_results.idset==i],
                bottom=error_list2[fit.fit_results.idset==i],
                beam=0.0, pen='#000000')
                err2a.setZValue(-10)
                p00.addItem(err2a)

        if self.RV_o_c_plot_cross_hair.isChecked():
            self.cross_hair(p00,log=False)       
            
        if self.RV_plot_autorange.isChecked():
            p2.autoRange()
            


#### Transit plots ################ 
    def update_transit_plots(self): 
        global fit, p3, p4, colors,legend_tra

        p3.plot(clear=True,)
        p4.plot(clear=True,)


        #if self.tra_plot_add_o_c.isChecked():
 

        if not self.hold_old_plot_tra.isChecked():    
            p30.plot(clear=True,)
            p31.plot(clear=True,)
        else:
            p30.plot(clear=False,)
            p31.plot(clear=False,)            
            
        if self.tra_legend.isChecked()==True:
            legend_tra.clear()
            legend_tra.setVisible(True)
        else:
            legend_tra.setVisible(False)

        #p30.addLine(x=None, y=0,   pen=pg.mkPen('#ff9933', width=0.8))


        self.check_tra_symbol_sizes()

        if len([x for x in range(len(fit.tra_data_sets)) if len(fit.tra_data_sets[x]) != 0]) == 0:
            return

        transit_results_sep = fit.transit_results[1]
        transit_results_all = fit.transit_results[2]

        if self.use_rich_tra_model.isChecked() and fit.tra_doGP == False:
            transit_model_rich  = fit.transit_results[3]
            t_model             = np.array(transit_model_rich[0], dtype=float)
            flux_model_ex       = np.array(transit_model_rich[1], dtype=float)
        elif self.use_rich_tra_model.isChecked() and fit.tra_doGP == True:
            print("Not yet possible to plot 'rich transit model' with a GP model")
            t_model        = np.concatenate([np.array(transit_results_sep[0][x], dtype=float) for x in range(len(fit.tra_data_sets)) if len(transit_results_sep[0][x]) != 0])
            flux_model_ex  = np.concatenate([np.array(transit_results_sep[3][x], dtype=float) for x in range(len(fit.tra_data_sets)) if len(transit_results_sep[3][x]) != 0])
            flux_model_ex_gp  = np.concatenate([np.array(transit_results_sep[6][x], dtype=float) for x in range(len(fit.tra_data_sets)) if len(transit_results_sep[6][x]) != 0])            
            

        else:
            t_model        = np.concatenate([np.array(transit_results_sep[0][x], dtype=float) for x in range(len(fit.tra_data_sets)) if len(transit_results_sep[0][x]) != 0])
            flux_model_ex  = np.concatenate([np.array(transit_results_sep[3][x], dtype=float) for x in range(len(fit.tra_data_sets)) if len(transit_results_sep[3][x]) != 0])
            flux_model_ex_gp  = np.concatenate([np.array(transit_results_sep[6][x], dtype=float) for x in range(len(fit.tra_data_sets)) if len(transit_results_sep[6][x]) != 0])   



        for j in range(len(fit.tra_data_sets)):

            if len(transit_results_sep[0][j]) == 0:
                continue

            t            = np.array(transit_results_sep[0][j], dtype=float)
            flux         = np.array(transit_results_sep[1][j], dtype=float)
            flux_err     = np.array(transit_results_sep[2][j], dtype=float)
            flux_model   = np.array(transit_results_sep[3][j], dtype=float)

            flux_model_gp   = np.array(transit_results_sep[6][j], dtype=float)
            
            if fit.tra_doGP == True:
                if self.plot_transit_GP_model.isChecked():
                    tr_o_c       = np.array(transit_results_sep[5][j], dtype=float)                
                else:
                    tr_o_c       = np.array(transit_results_sep[4][j], dtype=float)
            else:
                tr_o_c       = np.array(transit_results_sep[4][j], dtype=float)
    
                

#        elif fit.doGP == True and self.plot_RV_GP_model.isChecked() ==  True and len(fit.fit_results.rv_model.o_c)==len(fit.gp_model_data[0]):
#            data_o_c = fit.fit_results.rv_model.o_c - fit.gp_model_data[0]
#            y_model_o_c = np.zeros(len(fit.fit_results.model))
#        else:
#            data_o_c = fit.fit_results.rv_model.o_c
#            y_model_o_c = np.zeros(len(fit.fit_results.model))

            #print(t)
           # print(flux)
            ############### Phase signal TBD this should not be here! ####################################

            if self.plot_phase_pholded_tran.isChecked() and fit.npl > 0:
                
                self.tra_xaxis_offset.setEnabled(True) 
                self.trans_phase_slider.setEnabled(True) 
             
                ph_pl_ind = self.comboBox_phase_pl_tran.currentIndex()
                offset = (self.trans_phase_slider.value()/100.0)*fit.P[ph_pl_ind]
                data_time_phase = np.array( (t  - offset)%fit.P[ph_pl_ind] )
                tra_offset_xaxis = self.tra_xaxis_offset.value()

                if self.extra_plot_Tra_phase_norm.isChecked():
                    data_time_phase = (data_time_phase/max(data_time_phase)) - 0.5

                sort = np.array(sorted(range(len(data_time_phase)), key=lambda k: data_time_phase[k])    )

                t      = data_time_phase[sort] + tra_offset_xaxis
                flux          = flux[sort] 
                flux_err      = flux_err[sort]
                flux_model    = flux_model[sort] 
                tr_o_c        = tr_o_c[sort]
                flux_model_gp = flux_model_gp[sort]
                
                #if self.extra_plot_Tra_phase_norm.isChecked():
               #     t = (t/max(t)) - 0.5   
                             
                if fit.tra_doGP == True: 
                    flux = flux - flux_model_gp
                    #print("GP model was subtracted!!!")
                
                if self.extra_plot_Tra_phase_norm.isChecked():
                    p31.setLabel('bottom', 'phase', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})
                    p3.setLabel('bottom', 'phase', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})            
                    p4.setLabel('bottom', 'phase', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})  
                else:
                    p31.setLabel('bottom', 'phase [days]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})
                    p3.setLabel('bottom', 'phase [days]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})                    
                    p4.setLabel('bottom', 'phase [days]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})               
 
                #if self.tra_plot_add_o_c.isChecked():
               #     p31.setLabel('bottom', 'phase [days]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})
                    
                    
               # else:
               #     p3.setLabel('bottom', 'phase [days]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})
                
                
            else:

                if self.tra_plot_add_o_c.isChecked():
                    p31.setLabel('bottom', 'BJD [days]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})
                   # p31.getAxis("bottom").setTickFont= self.plot_font
                    p31.getAxis("bottom").setStyle(tickTextOffset = 2)
                    p31.getAxis("left").setStyle(tickTextOffset = 2)
                    p30.getAxis("left").setStyle(tickTextOffset = 2)
                    p4.setLabel('bottom', 'BJD [days]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
                  #  p3.getAxis("bottom").setTickFont= self.plot_font
                    p4.getAxis("bottom").setStyle(tickTextOffset = 2)                
                    p4.getAxis("left").setStyle(tickTextOffset = 2) 
                else:
                    p3.setLabel('bottom', 'BJD [days]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
                  #  p3.getAxis("bottom").setTickFont= self.plot_font
                    p3.getAxis("bottom").setStyle(tickTextOffset = 2)                
                    p3.getAxis("left").setStyle(tickTextOffset = 2) 
                    p4.setLabel('bottom', 'BJD [days]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
                  #  p3.getAxis("bottom").setTickFont= self.plot_font
                    p4.getAxis("bottom").setStyle(tickTextOffset = 2)                
                    p4.getAxis("left").setStyle(tickTextOffset = 2) 

                self.tra_xaxis_offset.setEnabled(False) 
                self.trans_phase_slider.setEnabled(False)

            p31.setLabel('bottom', '%s'%p31.getAxis("bottom").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})
            p31.setLabel('left', '%s'%p31.getAxis("left").labelText, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})



            if self.tra_plot_add_o_c.isChecked():

                p30.plot(t, flux,
                pen=None,
                symbol=dill.copy(fit.pyqt_symbols_tra[j]),
                symbolPen={'color': fit.tra_colors[j]+"%02x"%int(dill.copy(fit.pyqt_color_alpha_tra[j])), 'width': 1.1},
                symbolSize=dill.copy(fit.pyqt_symbols_size_tra[j]),enableAutoRange=True,viewRect=True,
                symbolBrush=fit.tra_colors[j]+"%02x"%int(dill.copy(fit.pyqt_color_alpha_tra[j])),name=fit.tra_data_sets[j][-1] ) 
                
                err_ = pg.ErrorBarItem(x=t, y=flux, symbol = dill.copy(fit.pyqt_symbols_tra[j]),
                                      # height=flux_err, 
                                       top=flux_err, 
                                       bottom=flux_err,
                                       beam=0.0, pen=fit.tra_colors[j]+"%02x"%int(dill.copy(fit.pyqt_color_alpha_tra[j])))

                p30.addItem(err_)

                p31.plot(t, tr_o_c,
                pen=None,
                symbol=dill.copy(fit.pyqt_symbols_tra[j]),
                symbolPen={'color': fit.tra_colors[j]+"%02x"%int(dill.copy(fit.pyqt_color_alpha_tra[j])), 'width': 1.1},
                symbolSize=dill.copy(fit.pyqt_symbols_size_tra[j]),enableAutoRange=True,viewRect=True,
                symbolBrush=fit.tra_colors[j]+"%02x"%int(dill.copy(fit.pyqt_color_alpha_tra[j])) )

                err_ = pg.ErrorBarItem(x=t, y=tr_o_c, symbol=dill.copy(fit.pyqt_symbols_tra[j]),
               # height=flux_err,
                top=flux_err,
                bottom=flux_err,
                beam=0.0, pen=fit.tra_colors[j]+"%02x"%int(dill.copy(fit.pyqt_color_alpha_tra[j])))
                p31.addItem(err_)

            else:

                p3.plot(t, flux,
                pen=None,
                symbol=dill.copy(fit.pyqt_symbols_tra[j]),
                symbolPen={'color': fit.tra_colors[j]+"%02x"%int(dill.copy(fit.pyqt_color_alpha_tra[j])), 'width': 1.1},
                symbolSize=dill.copy(fit.pyqt_symbols_size_tra[j]),enableAutoRange=True,viewRect=True,
                symbolBrush=fit.tra_colors[j]+"%02x"%int(dill.copy(fit.pyqt_color_alpha_tra[j])) ) 
                
                err_ = pg.ErrorBarItem(x=t, y=flux, symbol = dill.copy(fit.pyqt_symbols_tra[j]),
                                      # height=flux_err, 
                                       top=flux_err, 
                                       bottom=flux_err,
                                       beam=0.0, pen=fit.tra_colors[j]+"%02x"%int(dill.copy(fit.pyqt_color_alpha_tra[j])))

                p3.addItem(err_)




            p4.plot(t, tr_o_c,
            pen=None,
            symbol=dill.copy(fit.pyqt_symbols_tra[j]),
            symbolPen={'color': fit.tra_colors[j]+"%02x"%int(dill.copy(fit.pyqt_color_alpha_tra[j])), 'width': 1.1},
            symbolSize=dill.copy(fit.pyqt_symbols_size_tra[j]),enableAutoRange=True,viewRect=True,
            symbolBrush=fit.tra_colors[j]+"%02x"%int(dill.copy(fit.pyqt_color_alpha_tra[j])) )

            err_ = pg.ErrorBarItem(x=t, y=tr_o_c, symbol=dill.copy(fit.pyqt_symbols_tra[j]),
           # height=flux_err,
            top=flux_err,
            bottom=flux_err,
            beam=0.0, pen=fit.tra_colors[j]+"%02x"%int(dill.copy(fit.pyqt_color_alpha_tra[j])))
            p4.addItem(err_)


        if fit.tra_doGP == True:
            y_model =  np.array(transit_results_all[6], dtype=float)   #fit.tra_gp_model_curve[0]
            if self.plot_transit_GP_model.isChecked():
                y_model_o_c = np.zeros(len(flux_model_ex))       
            else:
                y_model_o_c = np.array(transit_results_all[6], dtype=float)   - np.array(transit_results_all[3], dtype=float)    #fit.tra_gp_model_curve[0]
        else:
            y_model = flux_model_ex 
            y_model_o_c = np.zeros(len(flux_model_ex))
            
            
        if self.plot_phase_pholded_tran.isChecked() and fit.npl > 0:
 
            model_time_phase = np.array( (t_model - offset)%fit.P[ph_pl_ind]  )
            sort2 = np.array(sorted(range(len(model_time_phase)), key=lambda k: model_time_phase[k])    )

            if self.extra_plot_Tra_phase_norm.isChecked():
                model_time_phase = (model_time_phase/max(model_time_phase)) - 0.5
                #time_phase = ((time_phase)/max(time_phase)) -0.5  

            t_model = model_time_phase[sort2] + tra_offset_xaxis
            flux_model_ex    = flux_model_ex[sort2] 
            flux_model_ex_gp    = flux_model_ex_gp[sort2]
            
   
                    
            y_model = flux_model_ex             
            
            if fit.tra_doGP == True: 
                y_model = flux_model_ex #- flux_model_ex_gp
                print("GP model was subtracted!!!")
 

            

        if len(t_model) != 0:

            if self.tra_plot_add_o_c.isChecked():
                model_curve = p30.plot(t_model,y_model, pen={'color':  fit.tra_colors[-1], 'width': self.tra_model_width.value()+1},
            enableAutoRange=True,viewRect=True )
                model_curve_o_c = p31.plot(t_model,y_model_o_c, pen={'color':  fit.tra_colors[-1], 'width': self.tra_model_width.value()+1},
            enableAutoRange=True,viewRect=True ) 
            else:
                model_curve = p3.plot(t_model,y_model, pen={'color':  fit.tra_colors[-1], 'width': self.tra_model_width.value()+1},
            enableAutoRange=True,viewRect=True ) 

            model_curve_o_c = p4.plot(t_model,y_model_o_c, pen={'color':  fit.tra_colors[-1], 'width': self.tra_model_width.value()+1},
            enableAutoRange=True,viewRect=True ) 

            model_curve.setZValue(self.tra_model_z.value())
            model_curve_o_c.setZValue(self.tra_model_z.value())



        if self.trans_plot_cross_hair.isChecked():

            if self.tra_plot_add_o_c.isChecked():
                self.cross_hair(p30,log=False)
                self.cross_hair(p31,log=False)
            else:            
                self.cross_hair(p3,log=False)

        if self.trans_o_c_plot_cross_hair.isChecked():
            self.cross_hair(p4,log=False)  

        if self.tra_plot_autorange.isChecked():
            if self.tra_plot_add_o_c.isChecked():
                p30.autoRange()                
                p31.autoRange()
            else:
                p3.autoRange()
            p4.autoRange()

#plot_transit_GP_model

######################## Correlation plots ###################################### 

    def init_correlations_combo(self):
        global fit
        self.comboBox_corr_1.clear()
        self.comboBox_corr_2.clear()
        
        self.initialize_corr_y = {k: [] for k in range(20)}
        z = 0 
 
        if fit.ndset != 0:

           # for i in range(max(fit.fit_results.idset)+1):
            for i in range(len(fit.rv_data_sets)):
                if len(fit.rv_data_sets[i]) == 0:
                    continue         


                self.comboBox_corr_1.addItem('RV %s'%(i+1),i+1) 
                self.comboBox_corr_2.addItem('RV %s'%(i+1),i+1) 

                self.initialize_corr_y[z] = np.array([fit.fit_results.rv_model.jd[fit.fit_results.idset==i],
                                                      fit.fit_results.rv_model.rvs[fit.fit_results.idset==i], 
                                                      fit.fit_results.rv_model.rv_err[fit.fit_results.idset==i]])  
                z +=1


            for i in range(len(fit.rv_data_sets)):
                if len(fit.rv_data_sets[i]) == 0:
                    continue         

                self.comboBox_corr_1.addItem('RV o-c %s'%(i+1),i+1)         
                self.comboBox_corr_2.addItem('RV o-c %s'%(i+1),i+1)  
    
                self.initialize_corr_y[z] = np.array([fit.fit_results.rv_model.jd[fit.fit_results.idset==i],
                                                      fit.fit_results.rv_model.o_c[fit.fit_results.idset==i], 
                                                      fit.fit_results.rv_model.rv_err[fit.fit_results.idset==i]]) 
                z +=1                          
 
            
        for i in range(0,20,1):         
            if len(fit.act_data_sets[i]) != 0: 
                self.comboBox_corr_1.addItem('act. data %s'%(i+1),i+1)
                self.comboBox_corr_2.addItem('act. data %s'%(i+1),i+1)
                
                self.initialize_corr_y[z] = fit.act_data_sets[i] 
                z +=1



    def update_correlations_data_plots(self):
        global fit, colors,  p6 
        
        ind1 = self.comboBox_corr_1.currentIndex()
        ind2 = self.comboBox_corr_2.currentIndex()
 
        p6.plot(clear=True,)  
                
        
        self.color_corr.setStyleSheet("color: %s;"%colors[0]) 

        #p6.autoRange()     
 
        
        if ind1 == -1 and ind2 == -1:
        
            text_err = pg.TextItem('No Data!',color=(0,0,0))#, anchor=(0,0), border='w',color) #, fill=(0, 0, 255, 100))
            p6.addItem(text_err, ignoreBounds=True)   
            return        
        
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
                    err1 = pg.ErrorBarItem(x=self.initialize_corr_y[ind1][1], y=self.initialize_corr_y[ind2][1], symbol='o', 
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
        
        colorz = self.colorDialog.getColor(options=QtWidgets.QColorDialog.DontUseNativeDialog)
        colors[0]=colorz.name()

        self.update_correlations_data_plots()

    def corr_plot_x_labels(self):
        global fit

        text, okPressed = QtWidgets.QInputDialog.getText(self, "x-axis label","(No special characters!)", QtWidgets.QLineEdit.Normal, "")
 
        if okPressed and text != '':
            p6.setLabel('bottom', '%s'%text, units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})
        else:
            return
    
        self.update_correlations_data_plots()
 

    def corr_plot_y_labels(self):
        global fit
        
        text, okPressed = QtWidgets.QInputDialog.getText(self, "y-axis label","(No special characters!)", QtWidgets.QLineEdit.Normal, "")
        
        if okPressed and text != '':
            p6.setLabel('left', '%s'%text, units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})

        else:
            return

        self.update_correlations_data_plots()


######################## Activity plots ######################################  

    def init_activity_combo(self):
        global fit

        for i in range(20):
            self.comboBox_act_data_gls.addItem('act. data %s'%(i+1),i+1)       
            self.comboBox_act_data.addItem('act. data %s'%(i+1),i+1)       
                
        
 
    def update_activity_gls_plots(self,ind):
        global fit, colors,  p11 

        p11.plot(clear=True,) 

        ind_norm = self.gls_norm_combo.currentIndex()
      #  omega = 1/ np.logspace(np.log10(self.gls_min_period.value()), np.log10(self.gls_max_period.value()), num=int(self.gls_n_omega.value()))
        power_levels = np.array([self.gls_fap1.value(),self.gls_fap2.value(),self.gls_fap3.value()])

        if len(fit.act_data_sets[ind]) != 0 and len(fit.act_data_sets[ind][0]) > 5:

            act_per = gls.Gls((fit.act_data_sets[ind][0], fit.act_data_sets[ind][1],fit.act_data_sets[ind][2]), 
            #fast=True,  verbose=False, norm= "ZK",ofac=self.gls_ofac.value(), fbeg=omega[-1], fend=omega[ 0],)
            fast=True,  verbose=False, norm=self.norms[ind_norm], ofac=self.gls_ofac.value(), fbeg=1/self.gls_max_period.value(), fend=1/self.gls_min_period.value())            
            ######################## GLS ##############################
            if self.radioButton_act_GLS_period.isChecked():
                p11.setLogMode(True,False)        
                p11.plot(1/act_per.freq, act_per.power,pen=dill.copy(fit.act_colors[ind]),symbol=None ) 
                p11.setLabel('bottom', 'period [d]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 

            else:
                p11.setLogMode(False,False)        
                p11.plot(act_per.freq, act_per.power,pen=dill.copy(fit.act_colors[ind]),symbol=None )
                p11.setLabel('bottom', 'frequency [1/d]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 

            [p11.addLine(x=None, y=fap, pen=pg.mkPen('k', width=0.8, style=QtCore.Qt.PenStyle.DotLine)) for ii,fap in enumerate(act_per.powerLevel(np.array(power_levels)))]
  
            text_peaks, pos_peaks = self.identify_power_peaks(1/act_per.freq, act_per.power, power_level = power_levels, sig_level = act_per.powerLevel(np.array(power_levels)) )    

            self.label_peaks(p11, pos_peaks, GLS = True, activity = True)
            self.act_periodogram_print_info.clicked.connect(lambda: self.print_info_for_object(act_per.info(stdout=False) + text_peaks ))   

        if self.gls_act_cross_hair.isChecked():
            #self.cross_hair(p11,log=self.radioButton_act_GLS_period.isChecked()) 

            self.cross_hair(p11,log=self.radioButton_act_GLS_period.isChecked(), alias=[self.show_alias_GLS.isChecked(), self.alias_days_gls.value(), colors_GLS_alias[0]])    


    def update_activity_data_plots(self,ind):
        global fit, colors,  p5 

        if len(fit.act_data_sets[ind]) != 0:

            p5.plot(clear=True,)  

            err1 = pg.ErrorBarItem(x=fit.act_data_sets[ind][0], y=fit.act_data_sets[ind][1],symbol='o', 
           # height=fit.act_data_sets[ind][2], beam=0.0, pen=fit.colors[ind])  
            top=fit.act_data_sets[ind][2],
            bottom=fit.act_data_sets[ind][2],
            beam=0.0, pen=dill.copy(fit.act_colors[ind]))


            p5.addItem(err1)      
            p5.addLine(x=None, y=0, pen=pg.mkPen('#ff9933', width=0.8))

            p5.plot(fit.act_data_sets[ind][0],fit.act_data_sets[ind][1], pen=None,symbol='o',
            symbolSize=self.act_data_size.value(),enableAutoRange=True,viewRect=True,
            symbolBrush=dill.copy(fit.act_colors[ind])
            )

            p5.setLabel('left', 'y', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})

            return
        else:
            p5.plot(clear=True,)

            return



#### TTV plots ################
 
    def update_ttv_pl_index(self, o_c = False):

        if o_c == True:
            pl_ind = self.ttv_o_c_comboBox_pl.currentIndex()
            self.ttv_comboBox_pl.setCurrentIndex(pl_ind)
        else:
            pl_ind = self.ttv_comboBox_pl.currentIndex()
            #print(pl_ind)
            self.ttv_o_c_comboBox_pl.setCurrentIndex(pl_ind)

        self.update_ttv_plots()


    def update_ttv_plots(self): 
        global fit, p_ttv, p_ttv_oc,p_ttv_00,p_ttv_01, colors,legend_ttv
        
        self.check_ttv_symbol_sizes()

        pl_ind     = self.ttv_comboBox_pl.currentIndex()


        p_ttv.plot(clear=True,) 
        p_ttv_oc.plot(clear=True,)
        p_ttv_00.plot(clear=True,)
        p_ttv_01.plot(clear=True,)
 
        ttv_files = fit.ttv_data_sets

        if self.ttv_legend.isChecked()==True:
            legend_ttv.clear()
            legend_ttv.setVisible(True)
        else:
            legend_ttv.setVisible(False)

        fit.prepare_for_mcmc()
        #times = [float(fit.epoch),fit.time_step_model,float(fit.epoch)+400.0]
        times = fit.ttv_times

        vel_files = []
        for i in range(len(fit.rv_data_sets)):
            if len(fit.rv_data_sets[i]) == 0:
                continue
            vel_files.append(fit.rv_data_sets[i][5]) 
        
 
        for j in range(len(ttv_files)):
            
            if len(ttv_files[j]) == 0 or ttv_files[j][4] == False or ttv_files[j][3] != pl_ind+1:
                continue
            
            first_transit = min([min(ttv_files[x][1]) for x in range(10) if len(ttv_files[x]) != 0 and ttv_files[x][3] == pl_ind+1])
            

            t = np.array(ttv_files[j][0])
            flux = np.array(ttv_files[j][1])
            flux_err = np.array(ttv_files[j][2])
 
            
                #ttv_loglik = ttvs_loglik(par,vel_files,ttv_files,npl,stmass,times,fit_results, return_model = False)
            if fit.npl > 0:
                
                if fit.rtg[0] == False:
                    ttv_loglik = rv.ttvs_loglik(fit.parameters,vel_files, ttv_files,fit.npl,fit.params.stellar_mass,times,fit.hkl,fit_results = False, return_model = True)
                else:
                    ttv_loglik = rv.ttvs_loglik(fit.parameters,vel_files, ttv_files,fit.npl,fit.params.stellar_mass,times,fit.hkl,fit_results =fit.fit_results, return_model = True)

                fit.ttv_results = dill.copy(ttv_loglik)

                if ttv_loglik == None:
                    print("""
Failed to plot! Perhaps the number of computed transits is smaler than the number of the observed transits.
There is no good fix for that at the moment.... Maybe adjust the epoch and try again.
                          """)
                    continue
                
                if isinstance(ttv_loglik, float):
                    return

                if self.ttv_show_Ntransit.isChecked():
                    ttv_model = dill.copy(ttv_loglik[4][j][1]) #- flux[0] #ttv_loglik[4][j][1][0]                  
                    ttv_model_transits = []
                    model_N_transits = dill.copy(ttv_loglik[4][j][0])


                    if self.ttv_apply_mean_period.isChecked():
                        periods_t0 = [ttv_model[k+1] - ttv_model[k] for k in range(len(ttv_model)-1)]
                        mean_P = np.mean(periods_t0)                    
                    else:
                        mean_P = fit.P[int(ttv_files[j][3]-1)]

                    #model_N_transits =  ttv_loglik[4][j][0] + np.rint((ttv_loglik[4][j][1][0] - first_transit)/mean_P)
        
                    for k in range(len(ttv_loglik[4][j][1])):
                       # ttv_model[k] = ttv_model[k] - mean_P*(ttv_loglik[4][j][0][k]-1)

                        ttv_model[k] = ttv_model[k] - (mean_P*(ttv_loglik[4][j][0][k]-1) + first_transit) 
                    for k in range(len(ttv_loglik[3][j][1])):
                        flux[k] = flux[k] - (mean_P*(ttv_loglik[3][j][0][k]-1) + first_transit) #ttv_files[j][1][0])
                        ttv_model_transits.append(ttv_model[ttv_loglik[3][j][0][k]-1])
 

                elif self.ttv_show_BJD.isChecked():

                    ttv_model = dill.copy(ttv_loglik[4][j][1]) #-ttv_loglik[4][j][1][0]                  
                    ttv_model_transits = []
                    model_N_transits = dill.copy(ttv_loglik[4][j][0])

                    if self.ttv_apply_mean_period.isChecked():
                        periods_t0 = [ttv_model[k+1] - ttv_model[k] for k in range(len(ttv_model)-1)]
                        mean_P = np.mean(periods_t0)                    
                    else:
                        mean_P = fit.P[int(ttv_files[j][3]-1)]

                   # model_N_transits =  ttv_loglik[4][j][0] + np.rint((ttv_loglik[4][j][1][0] - first_transit)/mean_P)
        
                    for k in range(len(ttv_loglik[4][j][1])):
                        ttv_model[k] = ttv_model[k] #- mean_P*(ttv_loglik[4][j][0][k]-1)
                        model_N_transits[k] = (mean_P*(ttv_loglik[4][j][0][k]-1) + first_transit) 
                    for k in range(len(ttv_loglik[3][j][1])):
                        flux[k] = flux[k] #- (mean_P*(ttv_loglik[3][j][0][k]-1) + first_transit) #ttv_files[j][1][0])
                        t[k] = (mean_P*(ttv_loglik[3][j][0][k]-1) + first_transit) 
                        ttv_model_transits.append(ttv_model[ttv_loglik[3][j][0][k]-1])
                        #t[k] = t[k]*fit.P[int(ttv_files[j][3]-1)] + first_transit



                flux_o_c = flux-ttv_model_transits
 
                if self.ttv_subtract_mean.isChecked():
                    ttv_model = ttv_model - np.mean(flux)
                    flux = flux - np.mean(flux)
                                        

            else:
                ttv_model = np.zeros(len(flux))+ np.mean(flux)
                ttv_model_transits = np.zeros(len(flux))+ np.mean(flux)
                model_N_transits = t

                flux_o_c = flux-ttv_model_transits


            if self.ttv_o_c.isChecked():   

                ttv_data = flux_o_c         
                ttv_model_ = np.zeros(len(ttv_model)) 

  
            else:
                ttv_data = flux         
                ttv_model_ = ttv_model         

            if self.ttv_model_interpolate.isChecked():

                ind_ttv_inter = self.comboBox_ttv_model_interpolate.currentIndex()

                func_interpol = interp1d(model_N_transits, ttv_model_, kind=ttv_interpol_opt[ind_ttv_inter])
                model_N_transits = np.linspace(min(model_N_transits),max(model_N_transits), num=int(self.ttv_model_interpolate_points.value()), endpoint=True)
                ttv_model_ = func_interpol(model_N_transits) 


            p_ttv.plot(t, ttv_data,
            pen=None,
            symbol=dill.copy(fit.pyqt_symbols_ttv[j]),
            symbolPen={'color': dill.copy(fit.ttv_colors[j]), 'width': 1.1},
            symbolSize=dill.copy(fit.pyqt_symbols_size_ttv[j]),enableAutoRange=True,viewRect=True,
            symbolBrush=dill.copy(fit.ttv_colors[j]) ) 
            
            err_ = pg.ErrorBarItem(x=t, y=ttv_data, symbol=dill.copy(fit.pyqt_symbols_ttv[j]),
                                  # height=flux_err, 
                                   top=flux_err, 
                                   bottom=flux_err,
                                   beam=0.0, pen=dill.copy(fit.ttv_colors[j]))

            p_ttv.addItem(err_)

            model_curve = p_ttv.plot(model_N_transits,ttv_model_,  pen={'color':  fit.ttv_colors[-1], 'width': self.ttv_model_width.value()+1}, enableAutoRange=True, viewRect=True )
            model_curve.setZValue(self.ttv_model_z.value())

            if self.ttv_plot_cross_hair.isChecked():
                self.cross_hair(p_ttv,log=False)

            #######################################
            #p_ttv_00.addLine(x=None, y=0,   pen=pg.mkPen('#ff9933', width=0.8))

            model_curve_oc = p_ttv_00.plot(model_N_transits,ttv_model_,  pen={'color':  fit.ttv_colors[-1], 'width': self.ttv_model_width.value()+1}, enableAutoRange=True, viewRect=True )
            model_curve_oc.setZValue(self.ttv_model_z.value())

            p_ttv_00.plot(t, ttv_data,
            pen=None,  
            symbol=dill.copy(fit.pyqt_symbols_ttv[j]),
            symbolPen={'color': dill.copy(fit.ttv_colors[j]), 'width': 1.1},
            symbolSize=dill.copy(fit.pyqt_symbols_size_ttv[j]),enableAutoRange=True,viewRect=True,
            symbolBrush=dill.copy(fit.ttv_colors[j]),name=fit.ttv_data_sets[j][-1] )

            err_ = pg.ErrorBarItem(x=t, y=ttv_data, symbol=dill.copy(fit.pyqt_symbols_ttv[j]),
           # height=flux_err,
            top=flux_err,
            bottom=flux_err,
            beam=0.0, pen=dill.copy(fit.ttv_colors[j]))

            p_ttv_00.addItem(err_)

            #######################################
            p_ttv_01.addLine(x=None, y=0,   pen=pg.mkPen('#ff9933', width=0.8))

            model_curve_oc = p_ttv_01.plot(model_N_transits,np.zeros(len(ttv_model_)) ,  pen={'color':  fit.ttv_colors[-1], 'width': self.ttv_model_width.value()+1}, enableAutoRange=True, viewRect=True )
            model_curve_oc.setZValue(self.ttv_model_z.value())


            p_ttv_01.plot(t, flux_o_c,
            pen=None,  
            symbol=dill.copy(fit.pyqt_symbols_ttv[j]),
            symbolPen={'color': dill.copy(fit.ttv_colors[j]), 'width': 1.1},
            symbolSize=dill.copy(fit.pyqt_symbols_size_ttv[j]),enableAutoRange=True,viewRect=True,
            symbolBrush=dill.copy(fit.ttv_colors[j]) )

            err_oc = pg.ErrorBarItem(x=t, y=flux_o_c, symbol=dill.copy(fit.pyqt_symbols_ttv[j]),
            #height=flux_err,
            top=flux_err,
            bottom=flux_err,
            beam=0.0, pen=dill.copy(fit.ttv_colors[j]))

            p_ttv_01.addItem(err_oc)


            if self.ttv_o_c_plot_cross_hair.isChecked():
                #self.cross_hair(p_ttv_oc,log=False)  
                self.cross_hair(p_ttv_00,log=False)  
                self.cross_hair(p_ttv_01,log=False)  

        if self.ttv_plot_autorange.isChecked():
            p_ttv.autoRange()           
            p_ttv_oc.autoRange()  




#    def get_ttv_plot_color(self):
#        global fit
        
#        colorz = self.colorDialog.getColor(options=QtWidgets.QColorDialog.DontUseNativeDialog)
#        fit.ttv_colors[0]=colorz.name()   
        
#        self.update_ttv_plots() 
        
#    def get_ttv_model_color(self):
#        global fit
        
#        colorz = self.colorDialog.getColor(options=QtWidgets.QColorDialog.DontUseNativeDialog)
#        fit.ttv_colors[-1]=colorz.name()   
        
#        self.update_ttv_plots() 


    def ttv_pl_combo(self):
        global fit

        self.ttv_comboBox_pl.clear()
        self.ttv_o_c_comboBox_pl.clear()
        
        for i in range(9):

            if not bool(fit.use_planet[i]):
                continue

            self.ttv_comboBox_pl.addItem('Planet %s'%str(i+1),i+1) 
            self.ttv_o_c_comboBox_pl.addItem('Planet %s'%str(i+1),i+1) 
            
        self.ttv_comboBox_pl.setCurrentIndex(0)
        self.ttv_o_c_comboBox_pl.setCurrentIndex(0)


#### Ast plots ################

    def update_ast_pl_index(self, o_c = False):

        if o_c == True:
            pl_ind = self.ast_o_c_comboBox_pl.currentIndex()
            self.ast_comboBox_pl.setCurrentIndex(pl_ind)
        else:
            pl_ind = self.ast_comboBox_pl.currentIndex()
            print(pl_ind)
            self.ast_o_c_comboBox_pl.setCurrentIndex(pl_ind)

        self.update_ast_plots()

    def update_ast_plots(self):
        
        self.update_ast_plots_img()
        
        if self.use_hipp_gaia_ast.isChecked():
        
            self.update_ast_plots_hipp()        
        
        
    def update_ast_plots_img(self):
        global fit, p_ast, p_ast_oc,p_ast_00,p_ast_01, colors,legend_ast
        
        self.check_ast_symbol_sizes()

        pl_ind     = self.ast_comboBox_pl.currentIndex()
 
        p_ast.plot(clear=True,) 
        p_ast_oc.plot(clear=True,)
        p_ast_00.plot(clear=True,)
        p_ast_01.plot(clear=True,)
 
        ast_files = fit.ast_data_sets

#        print(ast_files, pl_ind+1)

        if self.ast_legend.isChecked()==True:
            legend_ast.clear()
            legend_ast.setVisible(True)
        else:
            legend_ast.setVisible(False)

        fit.prepare_for_mcmc()
        #times = [float(fit.epoch),fit.time_step_model,float(fit.epoch)+400.0]
#        fit.ast_times[0] = fit.epoch
        times = fit.ast_times
 
        vel_files = []
        for i in range(len(fit.rv_data_sets)):
            if len(fit.rv_data_sets[i]) == 0:
                continue
            vel_files.append(fit.rv_data_sets[i][5]) 

       
        for j in range(len(ast_files)):
            
            if len(ast_files[j]) == 0 or ast_files[j][6] == False or ast_files[j][5] != pl_ind+1:
                continue
            
            #first_transit = min([min(ast_files[x][1]) for x in range(10) if len(ast_files[x]) != 0 and ast_files[x][5] == pl_ind+1])
            

            t = np.array(ast_files[j][0])
            x_axis = np.array(ast_files[j][1])
            x_axis_err = np.array(ast_files[j][2])
            y_axis = np.array(ast_files[j][3])
            y_axis_err = np.array(ast_files[j][4]) 
            
            if fit.npl > 0:
                
                if fit.rtg[0] == False:
                    ast_loglik = rv.ast_loglik(fit.parameters,vel_files, ast_files,fit.npl,fit.params.stellar_mass,times,fit.hkl,fit_results = False, return_model = True)
                else:
                    ast_loglik = rv.ast_loglik(fit.parameters,vel_files, ast_files,fit.npl,fit.params.stellar_mass,times,fit.hkl,fit_results =fit.fit_results, return_model = True)

               # print("TEst", rv.ast_loglik(fit.parameters,vel_files, ast_files,fit.npl,fit.params.stellar_mass,times,fit.hkl,fit_results = False, return_model = False))   

                if ast_loglik == None:
                    print("Something went wrong when calculating astr. lnL....")
                    continue

                if isinstance(ast_loglik, float):
                    return

                fit.ast_results = dill.copy(ast_loglik)

                ast_model_x = ast_loglik[2][pl_ind][0] 
                ast_model_y = ast_loglik[2][pl_ind][1]
            else:
                ast_model_x = np.zeros(5)
                ast_model_y = np.zeros(5)
                
 
            p_ast.plot(x_axis, y_axis,
            pen=None,
            symbol=dill.copy(fit.pyqt_symbols_ast[j]),
            symbolPen={'color': dill.copy(fit.ast_colors[j]), 'width': 1.1},
            symbolSize=dill.copy(fit.pyqt_symbols_size_ast[j]),enableAutoRange=True,viewRect=True,
            symbolBrush=dill.copy(fit.ast_colors[j]) ) 
            

            err_ = pg.ErrorBarItem(x=x_axis, y=y_axis, symbol=dill.copy(fit.pyqt_symbols_ast[j]),
                                  # height=flux_err, 
                                   top=y_axis_err, 
                                   bottom=y_axis_err,
                                   left=y_axis_err,
                                   right=y_axis_err,
                                   beam=0.0, pen=dill.copy(fit.ast_colors[j]))

            p_ast.addItem(err_)

            model_curve = p_ast.plot(ast_model_x, ast_model_y, pen={'color':  fit.ast_colors[-1], 'width': self.ast_model_width.value()}, enableAutoRange=True, viewRect=True )
            model_curve.setZValue(self.ast_model_z.value())

            if self.ast_plot_cross_hair.isChecked():
                self.cross_hair(p_ast,log=False)

            #######################################
 

            model_curve_oc = p_ast_00.plot(ast_model_x, ast_model_y, pen={'color':  fit.ast_colors[-1], 'width': self.ast_model_width.value()+1}, enableAutoRange=True, viewRect=True )
            model_curve_oc.setZValue(self.ast_model_z.value())

            p_ast_00.plot(x_axis, y_axis,
            pen=None,  
            symbol=dill.copy(fit.pyqt_symbols_ast[j]),
            symbolPen={'color': dill.copy(fit.ast_colors[j]), 'width': 1.1},
            symbolSize=dill.copy(fit.pyqt_symbols_size_ast[j]),enableAutoRange=True,viewRect=True,
            symbolBrush=dill.copy(fit.ast_colors[j]),name=fit.ast_data_sets[j][-1] )

            err_ = pg.ErrorBarItem(x=x_axis, y=y_axis, symbol=dill.copy(fit.pyqt_symbols_ast[j]),
                                  # height=flux_err, 
                                   top=y_axis_err, 
                                   bottom=y_axis_err,
                                   left=y_axis_err,
                                   right=y_axis_err,
                                   beam=0.0, pen=dill.copy(fit.ast_colors[j]))

            p_ast_00.addItem(err_)

            #######################################
#            p_ast_01.addLine(x=None, y=0,   pen=pg.mkPen('#ff9933', width=0.8))

#            model_curve_oc = p_ast_01.plot(model_N_transits,np.zeros(len(ast_model_)) ,  pen={'color':  fit.ast_colors[-1], 'width': self.ast_model_width.value()+1}, enableAutoRange=True, viewRect=True )
#            model_curve_oc.setZValue(self.ast_model_z.value())


#            p_ast_01.plot(t, flux_o_c,
#            pen=None,  
#            symbol=dill.copy(fit.pyqt_symbols_ast[j]),
#            symbolPen={'color': dill.copy(fit.ast_colors[j]), 'width': 1.1},
#            symbolSize=dill.copy(fit.pyqt_symbols_size_ast[j]),enableAutoRange=True,viewRect=True,
#            symbolBrush=dill.copy(fit.ast_colors[j]) )

#            err_oc = pg.ErrorBarItem(x=t, y=flux_o_c, symbol=dill.copy(fit.pyqt_symbols_ast[j]),
            #height=flux_err,
#            top=flux_err,
#            bottom=flux_err,
#            beam=0.0, pen=dill.copy(fit.ast_colors[j]))

#            p_ast_01.addItem(err_oc)


            if self.ast_o_c_plot_cross_hair.isChecked():
#                #self.cross_hair(p_ast_oc,log=False)  
                self.cross_hair(p_ast_00,log=False)  
                self.cross_hair(p_ast_01,log=False)  

        if self.ast_plot_autorange.isChecked():
            p_ast.autoRange()           
            p_ast_oc.autoRange()  

        
    def update_ast_plots_hipp(self):
        global fit, p_ast_hipp, p_ast_oc,p_ast_00,p_ast_01, colors,legend_ast_hipp
        
        self.check_ast_symbol_sizes()
 
        pl_ind     = self.ast_comboBox_pl.currentIndex()
 
        p_ast_hipp.plot(clear=True,) 

            
        ast_files = fit.ast_data_sets_hipp_gaia
 
        if self.ast_legend.isChecked()==True:
            legend_ast_hipp.clear()
            legend_ast_hipp.setVisible(True)
        else:
            legend_ast_hipp.setVisible(False)

        fit.prepare_for_mcmc()
        times = fit.ast_times 
 
        vel_files = []
        for i in range(len(fit.rv_data_sets)):
            if len(fit.rv_data_sets[i]) == 0:
                continue
            vel_files.append(fit.rv_data_sets[i][5]) 

       
        for j in range(len(ast_files)):
            
            if len(ast_files[j]) == 0 or ast_files[j][8] == False or ast_files[j][7] != pl_ind+1:
                continue
 

           # t = np.array(ast_files[j][0])
           # x_axis = np.array(ast_files[j][1])
           # x_axis_err = np.array(ast_files[j][2])
          #  y_axis = np.array(ast_files[j][3])
          #  y_axis_err = np.array(ast_files[j][4]) 
            
            if fit.npl > 0:
                
                if fit.rtg[0] == False:
                    ast_loglik = rv.ast_loglik_hipp(fit.parameters,vel_files, ast_files,fit.npl,fit.params.stellar_mass,times,fit.hkl,fit_results = False, return_model = True)
                else:
                    ast_loglik = rv.ast_loglik_hipp(fit.parameters,vel_files, ast_files,fit.npl,fit.params.stellar_mass,times,fit.hkl,fit_results =fit.fit_results, return_model = True)

 
                if ast_loglik == None:
                    print("Something went wrong when calculating astr. lnL....")
                    continue

                if isinstance(ast_loglik, float):
                    return

                fit.ast_results = dill.copy(ast_loglik)


                if self.ast_o_c_hipp_gaia.isChecked():
                    x_axis = np.array(ast_loglik[1][pl_ind][1][0] )
                    x_axis_err = np.array(ast_loglik[1][pl_ind][1][1])
                    y_axis = np.array(ast_loglik[1][pl_ind][1][2])
                    y_axis_err = np.array(ast_loglik[1][pl_ind][1][3])
                else:
                    x_axis = np.array(ast_loglik[1][pl_ind][0][0] )
                    x_axis_err = np.array(ast_loglik[1][pl_ind][0][1])
                    y_axis = np.array(ast_loglik[1][pl_ind][0][2])
                    y_axis_err = np.array(ast_loglik[1][pl_ind][0][3])                
               
                 
                bjd_dates = np.array(ast_loglik[1][pl_ind][2]) 
                 
                ast_model_x   = ast_loglik[2][pl_ind][0] 
                ast_model_y   = ast_loglik[2][pl_ind][1]
                ast_model_bjd = ast_loglik[2][pl_ind][2]                
            else:
                ast_model_x = np.zeros(5)
                ast_model_y = np.zeros(5)
                
            
            if self.ast_use_plot_ra_dec.isChecked():
               
                p_ast_hipp.invertX(False)
                p_ast_hipp.getViewBox().setAspectLocked(False)
                
                if self.ast_plot_ra.isChecked():
                    p_ast_hipp.setLabel('left', '<html><head/><body><p>&#916; &#945; cos(&#948;) [mas]</p></body></html>', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
                    p_ast_hipp.setLabel('bottom', '<html><head/><body><p>  BJD [days] </p></body></html>', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})                 
                    yaxis = x_axis
                    yaxis_err = x_axis_err
                    yaxis_model = ast_model_x
                    
                elif self.ast_plot_dec.isChecked():
                    p_ast_hipp.setLabel('left', '<html><head/><body><p>&#948; [mas]</p></body></html>', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
                    p_ast_hipp.setLabel('bottom', '<html><head/><body><p>  BJD [days] </p></body></html>', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})                 
                    yaxis = y_axis
                    yaxis_err = y_axis_err                    
                    yaxis_model = ast_model_y                    
                
                p_ast_hipp.plot(bjd_dates, yaxis,
                pen=None,
                symbol=dill.copy(fit.pyqt_symbols_ast[j]),
                symbolPen={'color': dill.copy(fit.ast_colors[j]), 'width': 1.1},
                symbolSize=dill.copy(fit.pyqt_symbols_size_ast[j]),enableAutoRange=True,viewRect=True,
                symbolBrush=dill.copy(fit.ast_colors[j]),name=ast_files[j][-1]  )             
                
                err_ = pg.ErrorBarItem(x=bjd_dates, y=yaxis, symbol=dill.copy(fit.pyqt_symbols_ast[j]),
                                      # height=flux_err, 
                                       top=yaxis_err, 
                                       bottom=yaxis_err,
                                       #left=y_axis_err,
                                      # right=y_axis_err,
                                       beam=0.0, pen=dill.copy(fit.ast_colors[j]))

                p_ast_hipp.addItem(err_)               
                
                
                model_curve = p_ast_hipp.plot(ast_model_bjd, yaxis_model, pen={'color':  fit.ast_colors[-1], 'width': self.ast_model_width.value()}, enableAutoRange=True, viewRect=True )
                 

            else:

                p_ast_hipp.setLabel('bottom', '<html><head/><body><p>&#916; &#945; cos(&#948;) [mas]</p></body></html>', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
                p_ast_hipp.setLabel('left',   '<html><head/><body><p>&#948; [mas]</p></body></html>', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
            
                p_ast_hipp.invertX(True)
                p_ast_hipp.getViewBox().setAspectLocked(True)

                p_ast_hipp.plot(np.array([0,0]), np.array([0,0]), pen=None,symbol='o', symbolSize=8,enableAutoRange=True,viewRect=True, symbolBrush='r')    
                            
                p_ast_hipp.plot(x_axis, y_axis,
                pen=None,
                symbol=dill.copy(fit.pyqt_symbols_ast[j]),
                symbolPen={'color': dill.copy(fit.ast_colors[j]), 'width': 1.1},
                symbolSize=dill.copy(fit.pyqt_symbols_size_ast[j]),enableAutoRange=True,viewRect=True,
                symbolBrush=dill.copy(fit.ast_colors[j]),name=ast_files[j][-1]  ) 
                
                
                #x1=x_axis-x_axis_err
                #y1=y_axis-y_axis_err
                #x2=x_axis+x_axis_err
                #y2=y_axis+y_axis_err
     
                size=0.05*((max(x_axis)-min(x_axis))/2)
                norm=1/((x_axis_err**2 +y_axis_err**2)**0.5)
                x1=x_axis-y_axis_err*norm*size
                y1=y_axis+x_axis_err*norm*size
                x2=x_axis+y_axis_err*norm*size
                y2=y_axis-x_axis_err*norm*size
                #then plots a line for each measurement in size and angle of the error    
     
                # Normalize BJD dates to the range [0, 1]
                #bjd_min, bjd_max = np.min(bjd_dates), np.max(bjd_dates)
                #bjd_normalized = (bjd_dates - bjd_min) / (bjd_max - bjd_min)
                bjd_phased = ((bjd_dates - bjd_dates[0]) % fit.P[ast_files[j][7]]) / fit.P[ast_files[j][7]]  # Normalize to the range [0, 1]
                ast_model_bjd_phased = (ast_model_bjd % fit.P[ast_files[j][7]]) / fit.P[ast_files[j][7]]  # Normalize to the range [0, 1]            
                
                #start_time = time.time()   
                
                test1 = True          
                if test1:
                    # Create a color map
                    colors = pg.colormap.get('viridis').getLookupTable()  # Get colormap as a lookup table
                    positions = np.linspace(0, 1, len(colors))  # Positions for the colormap
                    cmap = pg.ColorMap(pos=positions, color=colors)  # Create a pyqtgraph ColorMap

                    x_coords = np.vstack((x1, x2)).T  # Shape (N, 2)
                    y_coords = np.vstack((y1, y2)).T  # Shape (N, 2)
                    # Create the color gradient for the phased BJD
                    gradient = pg.ColorMap(pos=positions, color=colors)
                    pen_colors = [gradient.map(b, mode='qcolor') for b in bjd_phased]  # List of QColor for each line

                    for i in range(len(x_coords)):
                        p_ast_hipp.plot(
                            x_coords[i],
                            y_coords[i],
                            pen={'color': pen_colors[i], 'width': self.ast_model_width.value()},
                        )
                else:
                    # Plot each line with color based on BJD
                    for i in range(len(x1)):
                        # Map normalized BJD to a color
                        color = cmap.map(bjd_phased[i], mode='qcolor')
                        
                        # Plot the line segment with the corresponding color
                        p_ast_hipp.plot(
                            np.array([x1[i], x2[i]]),
                            np.array([y1[i], y2[i]]),
                            pen={'color': color, 'width': self.ast_model_width.value()},
                        )

                    # Add a colorbar (gradient legend)
                    gradient = pg.GradientLegend((10, 200), (0, 0))  # Size and position
                    gradient.setGradient(pg.colormap.get('viridis').getGradient())  # Match the viridis colormap
                    gradient.setLabels({
                        0: f"0 (Start of Phase)",  # Start of the phase
                        1: f"1 (End of Phase)"     # End of the phase
                    })
                #p_ast_hipp.addItem(gradient)
     
                #print(" plot1 takes --- %s seconds ---" % (time.time() - start_time))     
                # Plot the model curve with gradient
                #start_time = time.time() 
     
                if self.ast_model_color_gradient.isChecked():

                
                    gradient = QtGui.QLinearGradient(QtCore.QPointF(np.min(ast_model_x), 0), QtCore.QPointF(np.max(ast_model_x), 0))
                    gradient.setColorAt(0, pg.mkColor('blue'))  # Start color
                    gradient.setColorAt(1, pg.mkColor('yellow'))  # End color

                    # Create a pen with the gradient
                    pen = QtGui.QPen()
                    pen.setBrush(gradient)  # Use the gradient as the brush for the pen
                    pen.setWidthF(self.ast_model_width.value())  # Adjust the line width with a float
                    pen.setCosmetic(True)  # Match simple pen behavior
                    # Plot the model curve with the gradient pen
                    model_curve = p_ast_hipp.plot(ast_model_x, ast_model_y)
                    model_curve.setPen(pen)
                     
                else:
                                    
                     model_curve = p_ast_hipp.plot(ast_model_x, ast_model_y, pen={'color':  fit.ast_colors[-1], 'width': self.ast_model_width.value()}, enableAutoRange=True, viewRect=True )



                #print(" plot2 takes --- %s seconds ---" % (time.time() - start_time))     


            model_curve.setZValue(self.ast_model_z.value())

            if self.ast_plot_cross_hair_hipp_gaia.isChecked():
                self.cross_hair(p_ast_hipp,log=False)



            #######################################
 
 

    def ast_pl_combo(self):
        global fit

        self.ast_comboBox_pl.clear()
        self.ast_o_c_comboBox_pl.clear()
        
        for i in range(9):

            if not bool(fit.use_planet[i]):
                continue
            self.ast_comboBox_pl.addItem('Planet %s'%str(i+1),i+1) 
            self.ast_o_c_comboBox_pl.addItem('Planet %s'%str(i+1),i+1) 
            
        self.ast_comboBox_pl.setCurrentIndex(0)
        self.ast_o_c_comboBox_pl.setCurrentIndex(0)

############################# N-Body plots ########################################     


    def update_orb_plot(self):
        global fit, p16
        
        p16.plot(clear=True,)    


        if len(fit.fit_results.a) ==0:
            return
        
        if fit.pl_arb_test == True:
            npl = fit.npl_arb
        else:
            npl = fit.npl     
        
        
        for i in range(9):
            if not bool(fit.use_planet[i]):
                continue
            orb_xyz, pl_xyz, peri_xyz, apo_xyz = rv.planet_orbit_xyz(fit,i)        
            p16.plot(orb_xyz[0],orb_xyz[1], pen={'color': 0.5, 'width': 1.1},enableAutoRange=True,viewRect=True)   
            p16.plot((0,peri_xyz[0]),(0,peri_xyz[1]), pen={'color': 0.5, 'width': 1.1},enableAutoRange=True,viewRect=True)               
            
            p16.plot((pl_xyz[0],pl_xyz[0]), (pl_xyz[1],pl_xyz[1] ), pen=None,symbol='o', symbolSize=6,enableAutoRange=True,viewRect=True, symbolBrush='b') 
            
        p16.plot(np.array([0,0]), np.array([0,0]), pen=None,symbol='o', symbolSize=8,enableAutoRange=True,viewRect=True, symbolBrush='r')                


    #### Period evolution plot ##################

        
    def per_rat_combo(self):
        global fit

        if fit.pl_arb_test == True:
            npl = fit.npl_arb
        else:
            npl = fit.npl
 
        self.per_evol_comboBox_pl_1.clear()
        self.per_evol_comboBox_pl_2.clear()
        
        for i in range(npl):
            self.per_evol_comboBox_pl_1.addItem('Planet %s'%str(i+1),i+1) 
            self.per_evol_comboBox_pl_2.addItem('Planet %s'%str(i+1),i+1) 
            
        self.per_evol_comboBox_pl_1.setCurrentIndex(1)
        self.per_evol_comboBox_pl_2.setCurrentIndex(0)


    def plot_per_rat(self):
        global fit, colors_per_rat, p_per_ev 
        

        if not self.hold_old_plot_per_evol.isChecked():    
            p_per_ev.plot(clear=True,)
        else:
            p_per_ev.plot(clear=False,)


        if self.Prat_plot_line.isChecked():
            Prat_plot_pen={'color':colors_per_rat[0], 'width':1.1}
        else:
            Prat_plot_pen=None


        if self.Prat_plot_dot.isChecked():
            Prat_plot_dot='o'
        else:
            Prat_plot_dot=None

        self.color_per_evol.setStyleSheet("color: %s;"%colors_per_rat[0]) 


        pl1_ind = self.per_evol_comboBox_pl_1.currentIndex()
        pl2_ind = self.per_evol_comboBox_pl_2.currentIndex()
 
        if pl1_ind ==-1 or pl2_ind ==-1 or len(fit.evol_Per[pl1_ind]) ==0 or len(fit.evol_Per[pl2_ind]) ==0:
            return
        else:
            last_stable = min(len(fit.evol_Per[pl1_ind]),len(fit.evol_Per[pl2_ind]))
        
        Prat = fit.evol_Per[pl1_ind][0:last_stable] / fit.evol_Per[pl2_ind][0:last_stable]


        p_per_ev.plot(fit.evol_T[0][0:last_stable], Prat ,pen=Prat_plot_pen, #{'color': colors[i], 'width': 1.1},
        symbol=Prat_plot_dot,
        symbolPen={'color': colors_per_rat[0], 'width': 1.1},
        symbolSize=1,enableAutoRange=True,viewRect=True,
        symbolBrush=dill.copy(fit.colors[0])
        )  

        if self.orb_evol_auto_range_per_evol.isChecked():
            p_per_ev.autoRange()   



    def get_per_rat_color(self):
        global fit, colors_per_rat
        
        #colorz = QtWidgets.QColorDialog.getColor()
        colorz = self.colorDialog.getColor(options=QtWidgets.QColorDialog.DontUseNativeDialog)
        colors_per_rat[0]=colorz.name()   
 
        self.plot_per_rat()

    def per_rat_plot_x_labels(self):
        global fit, p_per_ev
        
        text, okPressed = QtWidgets.QInputDialog.getText(self, "x-axis label","(No special characters!)", QtWidgets.QLineEdit.Normal, "")
        
        if okPressed and text != '':
            p_per_ev.setLabel('bottom', '%s'%text, units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})
 
        else:
            return
    
        self.plot_per_rat()
 

    def per_rat_plot_y_labels(self):
        global fit, p_per_ev
        
        text, okPressed = QtWidgets.QInputDialog.getText(self, "y-axis label","(No special characters!)", QtWidgets.QLineEdit.Normal, "")
        
        if okPressed and text != '':
            p_per_ev.setLabel('left', '%s'%text, units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})
 
        else:
            return
    
        self.plot_per_rat()




    ################### Delta \omega plots ################

    def delta_omega_combo(self):
        global fit

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


        if self.Dom_plot_line.isChecked():
            Dom_plot_pen={'color':colors_delta_om[0], 'width':1.1}
        else:
            Dom_plot_pen=None


        if self.Dom_plot_dot.isChecked():
            Dom_plot_dot='o'
        else:
            Dom_plot_dot=None

    
        if self.domega_esin_ecos.isChecked():
            
            if self.domega_esin_ecos_e_in.isChecked():
                e_ind = pl1_ind
            else:
                e_ind = pl2_ind
                     
            x = fit.evol_e[e_ind][0:last_stable]*np.cos(np.radians(dom))
            y = fit.evol_e[e_ind][0:last_stable]*np.sin(np.radians(dom)) 
            p17.setLabel('bottom', '<html><head/><body><p>e<sub>%s</sub> cos(&Delta;&omega;) </p></body></html>'%(e_ind+1), units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
            p17.setLabel('left', '<html><head/><body><p>e<sub>%s</sub> sin(&Delta;&omega;) </p></body></html>'%(e_ind+1), units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 


            if self.domega_esin_ecos_exe_out.isChecked():

                x = fit.evol_e[pl1_ind][0:last_stable]*fit.evol_e[pl2_ind][0:last_stable]*np.cos(np.radians(dom))
                y = fit.evol_e[pl1_ind][0:last_stable]*fit.evol_e[pl2_ind][0:last_stable]*np.sin(np.radians(dom)) 
                p17.setLabel('bottom', '<html><head/><body><p>e<sub>%s</sub>e<sub>%s</sub> cos(&Delta;&omega;) </p></body></html>'%(pl1_ind+1,pl2_ind+1), units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
                p17.setLabel('left',   '<html><head/><body><p>e<sub>%s</sub>e<sub>%s</sub> sin(&Delta;&omega;) </p></body></html>'%(pl1_ind+1,pl2_ind+1), units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 


            p17.getViewBox().setAspectLocked(True)
            p17.plot(np.array([0,0]), np.array([0,0]), pen=None,symbol='o', symbolSize=8,enableAutoRange=True,viewRect=True, symbolBrush='r')                

        else:
            x = fit.evol_T[0][0:last_stable]
            y = dom
            p17.setLabel('bottom', 't [yr]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
            p17.setLabel('left', '<html><head/><body><p>&Delta;&omega; [deg] </p></body></html>', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
            p17.getViewBox().setAspectLocked(False)
           
 
 
        p17.plot(x, y ,pen=Dom_plot_pen, #{'color': colors[i], 'width': 1.1},
        symbol=Dom_plot_dot,
        symbolPen={'color': colors_delta_om[0], 'width': 1.1},
        symbolSize=1,enableAutoRange=True,viewRect=True,
        symbolBrush=dill.copy(fit.colors[0])
        )  



        
    def get_delta_omega_color(self):
        global fit, colors_delta_om
        
        #colorz = QtWidgets.QColorDialog.getColor()
        colorz = self.colorDialog.getColor(options=QtWidgets.QColorDialog.DontUseNativeDialog)
        colors_delta_om[0]=colorz.name()   
 
        self.plot_delta_omega()

    def delta_omega_plot_x_labels(self):
        global fit, p17
        
        text, okPressed = QtWidgets.QInputDialog.getText(self, "x-axis label","(HTML special characters accepted!)", QtWidgets.QLineEdit.Normal, "")
        
        if okPressed and text != '':
           # p17.setLabel('bottom', '%s'%text, units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})
            p17.setLabel('bottom', '<html><head/><body><p>%s </p></body></html>'%text, units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})  
        else:
            return
    
        #self.plot_delta_omega()
 

    def delta_omega_plot_y_labels(self):
        global fit, p17
        
        text, okPressed = QtWidgets.QInputDialog.getText(self, "y-axis label","(HTML special characters accepted!)", QtWidgets.QLineEdit.Normal, "")
        
        if okPressed and text != '':
            #p17.setLabel('left', '%s'%text, units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})
            p17.setLabel('left', '<html><head/><body><p>%s </p></body></html>'%text, units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})   
        else:
            return
    
        #self.plot_delta_omega()        
  


      
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
        
        for i in range(20):
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

        if self.theta_plot_line.isChecked():
            theta_plot_pen={'color':colors_theta[0], 'width':1.1}
        else:
            theta_plot_pen=None


        if self.theta_plot_dot.isChecked():
            theta_plot_dot='o'
        else:
            theta_plot_dot=None

        
        if self.theta_esin_ecos.isChecked():
            
            if self.theta_esin_ecos_e_in.isChecked():
                e_ind = pl1_ind
            else:
                e_ind = pl2_ind
                     
            x = fit.evol_e[e_ind][0:last_stable]*np.cos(np.radians(theta[tet_n]))
            y = fit.evol_e[e_ind][0:last_stable]*np.sin(np.radians(theta[tet_n])) 
            p18.setLabel('bottom', '<html><head/><body><p>e<sub>%s</sub> cos(&theta;%s) </p></body></html>'%(e_ind+1,tet_n+1), units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
            p18.setLabel('left',   '<html><head/><body><p>e<sub>%s</sub> sin(&theta;%s) </p></body></html>'%(e_ind+1,tet_n+1), units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
            p18.getViewBox().setAspectLocked(True)

            p18.plot(np.array([0,0]), np.array([0,0]), pen=None,symbol='o', symbolSize=8,enableAutoRange=True,viewRect=True, symbolBrush='r')                

        else:
            x = fit.evol_T[0][0:last_stable]
            y = theta[tet_n]
            p18.setLabel('bottom', 't [yr]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
            p18.setLabel('left', '<html><head/><body><p>&theta;%s [deg] </p></body></html>'%(tet_n+1), units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
            p18.getViewBox().setAspectLocked(False)
           

        p18.plot(x, y ,pen=theta_plot_pen, #{'color': colors[i], 'width': 1.1},
        symbol=theta_plot_dot,
        symbolPen={'color': colors_theta[0], 'width': 1.1},
        symbolSize=1,enableAutoRange=True,viewRect=True,
        symbolBrush=dill.copy(fit.colors[0])
        )  
               
    def get_theta_color(self):
        global fit, colors_delta_om
        
        colorz = self.colorDialog.getColor(options=QtWidgets.QColorDialog.DontUseNativeDialog)
        colors_theta[0]=colorz.name()   
 
        self.plot_theta()

    def theta_plot_x_labels(self):
        global fit, p18
        
        text, okPressed = QtWidgets.QInputDialog.getText(self, "x-axis label","(HTML special characters accepted!)", QtWidgets.QLineEdit.Normal, "")
        
        if okPressed and text != '':

            p18.setLabel('bottom', '<html><head/><body><p>%s </p></body></html>'%text, units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})  
        else:
            return

    def theta_plot_y_labels(self):
        global fit, p18

        text, okPressed = QtWidgets.QInputDialog.getText(self, "x-axis label","(HTML special characters accepted!)", QtWidgets.QLineEdit.Normal, "")

        if okPressed and text != '':
            p18.setLabel('left', '<html><head/><body><p>%s </p></body></html>'%text, units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})  
        else:
            return

        #self.plot_theta()


    def plot_i_Om(self):
        global fit, colors, p19

        if fit.pl_arb_test == True:
            npl = fit.npl_arb
        else:
            npl = fit.npl  
            
        p19.plot(clear=True,)   

        if self.plot_i.isChecked():
            for i in range(npl):
                p19.plot(fit.evol_T[i], fit.evol_i[i] ,pen=dill.copy(fit.colors[i]),symbol=None )    
            p19.setLabel('left', 'i [deg]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})    
    
        elif self.plot_Om.isChecked():
            for i in range(npl):
                
                Om_evol = np.array(fit.evol_Om[i])
                
                if self.radioButton_Omega_180_fold.isChecked():
                    Om_evol[Om_evol>=180.0] -= 360.0

                p19.plot(fit.evol_T[i], Om_evol ,pen=None, #{'color': colors[i], 'width': 1.1},
                symbol='o',
                symbolPen={'color': dill.copy(fit.colors[i]), 'width': 1.1},
                symbolSize=1,enableAutoRange=True,viewRect=True,
                symbolBrush=dill.copy(fit.colors[i])
                )
                
            p19.setLabel('left', '<html><head/><body><p>&Omega; [deg] </p></body></html>', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
                
            Om_evol = 0
 
    def plot_energy(self):
        global fit, colors, p20

        p20.plot(clear=True,)  
        
        if len(np.atleast_1d(fit.evol_T_energy)) < 3:
            return
        

        if self.radioButton_energy.isChecked():
 
            p20.plot(fit.evol_T_energy, fit.evol_energy ,pen=dill.copy(fit.colors[0]),symbol=None )    
            p20.setLabel('left', 'Energy', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})    
    
        elif self.radioButton_lx.isChecked():
            p20.plot(fit.evol_T_energy, fit.evol_momentum['lx'] ,pen=dill.copy(fit.colors[0]),symbol=None )    
            p20.setLabel('left', 'Momentum lx', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})    

        elif self.radioButton_ly.isChecked():
            p20.plot(fit.evol_T_energy, fit.evol_momentum['ly'] ,pen=dill.copy(fit.colors[0]),symbol=None )    
            p20.setLabel('left', 'Momentum ly', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
            
        elif self.radioButton_lz.isChecked():
            p20.plot(fit.evol_T_energy, fit.evol_momentum['lz'] ,pen=dill.copy(fit.colors[0]),symbol=None )    
            p20.setLabel('left', 'Momentum lz', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})                    
 
    
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

            if self.a_plot_line.isChecked():
                a_plot_pen={'color':dill.copy(fit.colors[i]), 'width':self.a_line_size.value()}
            else:
                a_plot_pen=None

            if self.a_plot_dot.isChecked():
                a_plot_dot='o'
            else:
                a_plot_dot=None

            p13.plot(fit.evol_T[i], fit.evol_a[i],pen=a_plot_pen,  
        symbol=a_plot_dot,
        symbolPen={'color': dill.copy(fit.colors[i]), 'width': self.a_dot_size.value()},
        symbolSize=1,enableAutoRange=True,viewRect=True,
        symbolBrush=dill.copy(fit.colors[i])
        )    

        if self.orb_evol_auto_range_a.isChecked():
            p13.autoRange()   



         
     
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
            p14.plot(fit.evol_T[i], fit.evol_e[i] ,pen=dill.copy(fit.colors[i]),symbol=None )  

        if self.orb_evol_auto_range_e.isChecked():
            p14.autoRange()   
 
            
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
        symbolPen={'color': dill.copy(fit.colors[i]), 'width': 1.1},
        symbolSize=1,enableAutoRange=True,viewRect=True,
        symbolBrush=dill.copy(fit.colors[i])
        )    
            
        if self.orb_evol_auto_range_p.isChecked():
            p15.autoRange()



    def plot_evol_all(self):
        global fit

        self.plot_evol_a()
        self.plot_evol_e()
        self.plot_evol_p()

        self.per_rat_combo()
        self.delta_omega_combo()
        self.theta_which_pl_combo()
        
        self.plot_per_rat()
        self.plot_delta_omega()
        self.plot_theta()
        self.plot_i_Om()
        self.plot_energy()



    def change_nbody_plot_opt_tab(self):
 
        if self.tabWidget_3.currentIndex() ==  0:
            self.param_tabs.setCurrentWidget(self.plot_options_tabs)
            self.plot_opt_tab.setCurrentWidget(self.tab_149)
            self.tabWidget_nbody_plot_param.setCurrentWidget(self.tab_150)
        elif self.tabWidget_3.currentIndex() ==  1:
            self.param_tabs.setCurrentWidget(self.plot_options_tabs)
            self.plot_opt_tab.setCurrentWidget(self.tab_149)
            self.tabWidget_nbody_plot_param.setCurrentWidget(self.tab_151)
        elif self.tabWidget_3.currentIndex() ==  2:
            self.param_tabs.setCurrentWidget(self.plot_options_tabs)
            self.plot_opt_tab.setCurrentWidget(self.tab_149)
            self.tabWidget_nbody_plot_param.setCurrentWidget(self.tab_152) 
        elif self.tabWidget_3.currentIndex() ==  3:
            self.param_tabs.setCurrentWidget(self.plot_options_tabs)
            self.plot_opt_tab.setCurrentWidget(self.tab_149)
            self.tabWidget_nbody_plot_param.setCurrentWidget(self.tab_153) 
        elif self.tabWidget_3.currentIndex() ==  4:
            self.param_tabs.setCurrentWidget(self.plot_options_tabs)
            self.plot_opt_tab.setCurrentWidget(self.tab_149)
            self.tabWidget_nbody_plot_param.setCurrentWidget(self.tab_154) 


    def change_nbody_plot_opt_tab_res(self):
 
        if self.tabWidget_6.currentIndex() ==  0:
            self.param_tabs.setCurrentWidget(self.plot_options_tabs)
            self.plot_opt_tab.setCurrentWidget(self.tab_149)
            self.tabWidget_nbody_plot_param.setCurrentWidget(self.tab_125)
        elif self.tabWidget_6.currentIndex() ==  1:
            self.param_tabs.setCurrentWidget(self.plot_options_tabs)
            self.plot_opt_tab.setCurrentWidget(self.tab_149)
            self.tabWidget_nbody_plot_param.setCurrentWidget(self.tab_187)
        elif self.tabWidget_6.currentIndex() ==  2:
            self.param_tabs.setCurrentWidget(self.plot_options_tabs)
            self.plot_opt_tab.setCurrentWidget(self.tab_149)
            self.tabWidget_nbody_plot_param.setCurrentWidget(self.tab_plot_opt_res_angles) 


################ Extra Plots (work in progress) ######################

    def update_extra_plots(self):
        global fit

        self.comboBox_extra_plot.clear()

        if fit.npl != 0:
            for i in range(9):
                if not bool(fit.use_planet[i]):
                    continue
                self.comboBox_extra_plot.addItem('pl %s'%(i+1),i+1)
                #self.comboBox_extra_plot.addItem('pl %.3f'%(fit.P[i]),i+1)
                
                #self.comboBox_extra_plot.setCurrentIndex(i) 
            self.phase_plots(1)
        self.comboBox_extra_plot.activated.connect(self.handleActivated)   

 
    def handleActivated(self, index):
        global fit 

        ind = self.comboBox_extra_plot.itemData(index) 

        if ind == -1:
            return

#        elif ind <= fit.npl:
#            self.phase_plots(ind)
 
        else:
            self.phase_plots(ind)
            return



    def update_extra_plots2(self):
        global fit

        self.comboBox_extra_plot_gls.clear()
        if fit.ndset != 0:
            self.comboBox_extra_plot_gls.addItem('RV GLS',fit.npl+1)
            self.comboBox_extra_plot_gls.addItem('RV GLS o-c',fit.npl+2)   

            self.extra_RV_GLS_plots()

        self.comboBox_extra_plot_gls.activated.connect(self.handleActivated2)   

 
    def handleActivated2(self, index):
        global fit

        ind = self.comboBox_extra_plot_gls.itemData(index) 

#        if ind <= fit.npl:
#            self.phase_plots(ind)
        if ind == fit.npl+1: 
            self.extra_RV_GLS_plots()
            #p7.show()
            #pe.setScene(p7.scene())
            #pe.autoRange()


        elif ind == fit.npl+2: 
            self.extra_RV_GLS_o_c_plots()            
            #pe.setScene(p8.scene())
        else:
            return


 
    def phase_plots(self, ind, offset = 0):
        global fit, colors,pe0,pe1,legend_RV_phased

       # pe.plot(clear=True,)
        pe0.plot(clear=True,)
        pe0.setLogMode(False,False)
        pe1.plot(clear=True,)
        pe1.setLogMode(False,False)


        if self.RV_legend.isChecked()==True:
            legend_RV_phased.clear()
            legend_RV_phased.setVisible(True)
        else:
            legend_RV_phased.setVisible(False)

        ######## TBF #############
        if self.radioButton_transit.isChecked():
            return
        ########################

        ph_data = fit.ph_data[ind-1]
        ph_model = fit.ph_model[ind-1]
 
        offset = (self.RV_phase_slider.value()/100.0)*fit.P[ind-1] 

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
 
 

        time_phase = (ph_data[0]-offset)%fit.P[ind-1] 
        rv_data = ph_data[1]
        try:          
            rv_data_o_c = ph_data[4]
        except:
            rv_data_o_c = ph_data[1]



        model_time_phase = np.array((ph_model[0]-offset)%fit.P[ind-1] )

        sort = sorted(range(len(model_time_phase)), key=lambda k: model_time_phase[k])
        model_time_phase  = model_time_phase[sort] 
        ph_model =  ph_model[1][sort] 



        pe0.addLine(x=None, y=0, pen=pg.mkPen('#ff9933', width=0.8))
        pe1.addLine(x=None, y=0, pen=pg.mkPen('#ff9933', width=0.8))



        if self.extra_plot_RVphase_norm.isChecked():
            model_time_phase = (model_time_phase/max(model_time_phase)) - 0.5
            time_phase = ((time_phase)/max(time_phase)) -0.5 # - fit.params.planet_params[7*(ind-1)+1]/2.0) /fit.params.planet_params[7*(ind-1)+1]
 

        model_curve = pe0.plot(model_time_phase,ph_model, pen={'color':  dill.copy(fit.rvs_colors[-1]), 'width': self.rv_model_width.value()+1},
        enableAutoRange=True,viewRect=True)   
 
        model_curve.setZValue(self.RV_model_z.value())

#        for i in range(max(ph_data[3])+1):

#        N_RV_files = len([x for x in range(20) if len(fit.rv_data_sets[x]) != 0])
#        if N_RV_files == 0:
#            rv_filename = ''
#            pe0.plot(clear=True,)
#            pe1.plot(clear=True,)
#            return

        for i in range(20):

#            if len(fit.filelist.files) != 0:
            if len(fit.rv_data_sets[i]) == 0:
                continue
            else:
               # rv_filename = fit.filelist.files[i].name
                rv_filename = fit.rv_data_sets[i][4]


                pe0.plot((time_phase[ph_data[3]==i]),rv_data[ph_data[3]==i],
                pen=None, #{'color': colors[i], 'width': 1.1},
                symbol=dill.copy(fit.pyqt_symbols_rvs[i]),
                symbolPen={'color': dill.copy(fit.rvs_colors[i])+"%02x"%int(dill.copy(fit.pyqt_color_alpha_rvs[i])), 'width': 1.1},
                symbolSize=dill.copy(fit.pyqt_symbols_size_rvs[i]),enableAutoRange=True,viewRect=True,
                symbolBrush=dill.copy(fit.rvs_colors[i])+"%02x"%int(dill.copy(fit.pyqt_color_alpha_rvs[i])), name=rv_filename)

                pe1.plot((time_phase[ph_data[3]==i]),rv_data_o_c[ph_data[3]==i],
                pen=None, #{'color': colors[i], 'width': 1.1},
                symbol=dill.copy(fit.pyqt_symbols_rvs[i]),
                symbolPen={'color': dill.copy(fit.rvs_colors[i])+"%02x"%int(dill.copy(fit.pyqt_color_alpha_rvs[i])), 'width': 1.1},
                symbolSize=dill.copy(fit.pyqt_symbols_size_rvs[i]),enableAutoRange=True,viewRect=True,
                symbolBrush=dill.copy(fit.rvs_colors[i])+"%02x"%int(dill.copy(fit.pyqt_color_alpha_rvs[i])) )


                err_ = pg.ErrorBarItem(x=(time_phase[ph_data[3]==i]), y=rv_data[ph_data[3]==i],
                symbol=dill.copy(fit.pyqt_symbols_rvs[i]), 
                top=error_list[ph_data[3]==i],
                bottom=error_list[ph_data[3]==i],
                beam=0.0, pen=dill.copy(fit.rvs_colors[i])+"%02x"%int(dill.copy(fit.pyqt_color_alpha_rvs[i]))) 

                pe0.addItem(err_)

                err_o_c = pg.ErrorBarItem(x=(time_phase[ph_data[3]==i]), y=rv_data_o_c[ph_data[3]==i],
                symbol=dill.copy(fit.pyqt_symbols_rvs[i]), 
                top=error_list[ph_data[3]==i],
                bottom=error_list[ph_data[3]==i],
                beam=0.0, pen=dill.copy(fit.rvs_colors[i])+"%02x"%int(dill.copy(fit.pyqt_color_alpha_rvs[i]))) 

                pe1.addItem(err_o_c)

                if self.jitter_to_plots.isChecked() and self.split_jitter.isChecked():
     
                    err_2 = pg.ErrorBarItem(x=(time_phase[ph_data[3]==i]), y=rv_data[ph_data[3]==i],
                    symbol=dill.copy(fit.pyqt_symbols_rvs[i]), 
                    top=error_list2[ph_data[3]==i],
                    bottom=error_list2[ph_data[3]==i],
                    beam=0.0, pen=colors_RV_jitter[0])  
                    err_2.setZValue(-10)
                    pe0.addItem(err_2) 

                    err_2_o_c = pg.ErrorBarItem(x=(time_phase[ph_data[3]==i]), y=rv_data_o_c[ph_data[3]==i],
                    symbol=dill.copy(fit.pyqt_symbols_rvs[i]), 
                    top=error_list2[ph_data[3]==i],
                    bottom=error_list2[ph_data[3]==i],
                    beam=0.0, pen=colors_RV_jitter[0])  
                    err_2_o_c.setZValue(-10)
                    pe1.addItem(err_2_o_c) 

        if qso_mode:
            label = "flux [mag.]"
        else:
            label = "RV [m/s]" 

        pe0.setXRange(min(model_time_phase), max(model_time_phase), padding=0.002)
#        pe0.setLabel('bottom', 'phase [days]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})
        #pe0.setLabel('left',   '%s'%label, units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})  

        pe0.setLabel('left', '%s'%label, units='', **{'font-size':'%dpt'%self.plot_font.pointSize()})  

        if self.extra_plot_RVphase_norm.isChecked():
            pe1.setLabel('bottom', 'phase', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})
        else:
            pe1.setLabel('bottom', 'phase [days]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})

        if self.extra_plot_cross_hair.isChecked():
            self.cross_hair(pe0,log=False)   
            self.cross_hair(pe1,log=False)   


    ############### VERY VERY VERY Ugly fix !!!! 

    def extra_RV_GLS_plots(self):
        global fit,  pe2 
 
        pe2.plot(clear=True,)   
        power_levels = np.array([self.gls_fap1.value(),self.gls_fap2.value(),self.gls_fap3.value()])

        if len(fit.fit_results.rv_model.jd) > 5:
            ######################## GLS ##############################
            if self.radioButton_RV_GLS_period.isChecked():
                pe2.setLogMode(True,False)        
                pe2.plot(1/fit.gls.freq, fit.gls.power, pen='r',symbol=None ) 
                pe2.setLabel('bottom', 'period [d]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})    
                pe2.setLabel('left', 'Power', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})    
            else:
                pe2.setLogMode(False,False)
                pe2.plot(fit.gls.freq, fit.gls.power, pen='r',symbol=None ) 
                pe2.setLabel('bottom', 'frequency [1/d]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
                pe2.setLabel('left', 'Power', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})    
    
            if fit.gls.norm == 'ZK':
                [pe2.addLine(x=None, y=fap, pen=pg.mkPen('k', width=0.8, style=QtCore.Qt.PenStyle.DotLine)) for ii,fap in enumerate(fit.gls.powerLevel(np.array(power_levels)))]


            text_peaks, pos_peaks = self.identify_power_peaks(1/fit.gls.freq, fit.gls.power, power_level = power_levels, sig_level = fit.gls.powerLevel(np.array(power_levels)) )           
            self.label_peaks(pe2, pos_peaks, GLS = True, o_c = False)

        if self.extra_plot_cross_hair_gls.isChecked():
            self.cross_hair(pe2,log=self.radioButton_RV_GLS_period.isChecked())   

    def extra_RV_GLS_o_c_plots(self):
        global fit,  pe2 
 
        pe2.plot(clear=True,)
        power_levels = np.array([self.gls_fap1.value(),self.gls_fap2.value(),self.gls_fap3.value()])

        if len(fit.fit_results.rv_model.jd) > 5:
        ######################## GLS o-c ##############################
            if self.radioButton_RV_o_c_GLS_period.isChecked():
                pe2.setLogMode(True,False)        
                pe2.plot(1/fit.gls_o_c.freq, fit.gls_o_c.power, pen='r',symbol=None ) 
                pe2.setLabel('bottom', 'period [d]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})    
                pe2.setLabel('left', 'Power', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})    
            else:
                pe2.setLogMode(False,False)        
                pe2.plot(fit.gls_o_c.freq, fit.gls_o_c.power, pen='r',symbol=None )                    
                pe2.setLabel('bottom', 'frequency [1/d]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
                pe2.setLabel('left', 'Power', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})    
    
    
            if fit.gls.norm == 'ZK':
                [pe2.addLine(x=None, y=fap, pen=pg.mkPen('k', width=0.8, style=QtCore.Qt.PenStyle.DotLine)) for ii,fap in enumerate(fit.gls_o_c.powerLevel(np.array(power_levels)))]

            text_peaks, pos_peaks = self.identify_power_peaks(1/fit.gls_o_c.freq, fit.gls_o_c.power, power_level = power_levels, sig_level = fit.gls_o_c.powerLevel(np.array(power_levels)) )           
            self.label_peaks(pe2, pos_peaks, GLS = True, o_c = True)

        if self.extra_plot_cross_hair_gls.isChecked():
            self.cross_hair(pe2,log=self.radioButton_RV_o_c_GLS_period.isChecked())   




    def update_plots(self):
        global fit

        #self.update_RV_GLS_plots()
        #self.update_RV_o_c_GLS_plots()
        self.run_gls()
        self.run_gls_o_c()
       # self.run_WF(RV = True)

        self.update_RV_plots()
        self.update_extra_plots()
        self.update_extra_plots2()
        self.update_orb_plot()
        #self.change_extra_plot()
        self.update_transit_plots()    
        self.update_ttv_plots()
        self.update_ast_plots()

    def rv_plot_phase_change(self):
        global fit

        ind = self.comboBox_extra_plot.currentIndex()

        if ind == -1:
            return

        #elif ind+1 <= fit.npl:
       #     self.phase_plots(ind+1)
        else:
            self.phase_plots(ind+1)
            return

 

    def set_RV_GP(self):
        global fit

        if self.use_GP_sho_kernel.isChecked():
            fit.gp_kernel = 'SHOKernel'
        elif self.use_GP_rot_kernel.isChecked():
            fit.gp_kernel = 'RotKernel'    
        elif self.use_GP_mat_kernel.isChecked():
            fit.gp_kernel = 'Matern32' 
        elif self.use_GP_drw_kernel.isChecked():
            fit.gp_kernel = 'RealTerm'             
        elif self.use_GP_double_sho_kernel.isChecked():
            fit.gp_kernel = 'dSHOKernel' 

        self.set_link_GP()

    def set_tra_GP(self):
        global fit

        if self.use_tra_GP_sho_kernel.isChecked():
            fit.tra_gp_kernel = 'SHOKernel'
        elif self.use_tra_GP_rot_kernel.isChecked():
            fit.tra_gp_kernel = 'RotKernel'
        elif self.use_tra_GP_mat_kernel.isChecked():
            fit.tra_gp_kernel = 'Matern32'  
        elif self.use_tra_GP_drw_kernel.isChecked():
            fit.tra_gp_kernel = 'RealTerm'            
        elif self.use_tra_GP_double_sho_kernel.isChecked():
            fit.tra_gp_kernel = 'dSHOKernel'

        self.set_link_GP()


    def set_gui_RV_GP(self):
        global fit

        if fit.gp_kernel == 'SHOKernel':
            self.use_GP_sho_kernel.setChecked(True) 
        elif fit.gp_kernel == 'RotKernel':
            self.use_GP_rot_kernel.setChecked(True)
        elif fit.gp_kernel == 'Matern32':
            self.use_GP_mat_kernel.setChecked(True)         
        elif fit.gp_kernel == 'RealTerm':
            self.use_GP_drw_kernel.setChecked(True)              
        elif fit.gp_kernel == 'dSHOKernel':
            self.use_GP_double_sho_kernel.setChecked(True)             
               
    def set_gui_tra_GP(self):
        global fit

        if fit.tra_gp_kernel == 'SHOKernel':
            self.use_tra_GP_sho_kernel.setChecked(True) 
        elif fit.tra_gp_kernel == 'RotKernel':
            self.use_tra_GP_rot_kernel.setChecked(True)
        elif fit.tra_gp_kernel == 'Matern32':
            self.use_tra_GP_mat_kernel.setChecked(True)     
        elif fit.tra_gp_kernel == 'RealTerm':
            self.use_tra_GP_drw_kernel.setChecked(True)   
        elif fit.tra_gp_kernel == 'dSHOKernel':
            self.use_tra_GP_double_sho_kernel.setChecked(True)             


    def set_use_RV_GP(self):
        global fit

        if  self.do_RV_GP.isChecked():
            fit.doGP = True
            self.gls_o_c_GP.setEnabled(True)
        else:
            fit.doGP = False
            self.gls_o_c_GP.setChecked(False)
            self.gls_o_c_GP.setEnabled(False)

    def set_use_tra_GP(self):
        global fit

        if  self.do_tra_GP.isChecked():
            fit.tra_doGP = True
        else:
            fit.tra_doGP = False


#    def set_use_GP(self):
#        global fit
#
#        if  self.do_RV_GP.isChecked():
#            fit.doGP = True
#            self.gls_o_c_GP.setEnabled(True)
#        else:
#            fit.doGP = False
#            self.gls_o_c_GP.setChecked(False)
#            self.gls_o_c_GP.setEnabled(False)

#        if  self.do_tra_GP.isChecked():
#            fit.tra_doGP = True
#        else:
#            fit.tra_doGP = False
            
            
    def set_gui_use_GP(self):
        global fit

        if  fit.doGP:
            self.do_RV_GP.setChecked(True)
        else:
            self.do_RV_GP.setChecked(False)
 
        if  fit.tra_doGP:
            self.do_tra_GP.setChecked(True)
        else:
            self.do_tra_GP.setChecked(False)
 
        self.set_gui_RV_GP()
        self.set_gui_tra_GP()            
       
    def set_link_GP(self):
        global fit       
        
        
        if self.use_tra_GP_rot_kernel.isChecked():
                
            if self.use_tra_GP_rot_kernel_time_sc_link_to_RV.isChecked()  and self.use_GP_rot_kernel.isChecked():
                self.use_GP_rot_kernel_time_sc.setChecked(False)
                self.use_GP_rot_kernel_time_sc.setEnabled(False)
                self.GP_rot_kernel_time_sc.setEnabled(False)
            else:
                self.use_GP_rot_kernel_time_sc.setChecked(True)
                self.use_GP_rot_kernel_time_sc.setEnabled(True)
                self.GP_rot_kernel_time_sc.setEnabled(True)                            
                
            if self.use_tra_GP_rot_kernel_Per_link_to_RV.isChecked()  and self.use_GP_rot_kernel.isChecked():
                self.use_GP_rot_kernel_Per.setChecked(False)
                self.use_GP_rot_kernel_Per.setEnabled(False)
                self.GP_rot_kernel_Per.setEnabled(False)
            else:
                self.use_GP_rot_kernel_Per.setChecked(True)
                self.use_GP_rot_kernel_Per.setEnabled(True)
                self.GP_rot_kernel_Per.setEnabled(True)  
                
            if self.use_tra_GP_rot_kernel_fact_link_to_RV.isChecked()  and self.use_GP_rot_kernel.isChecked():
                self.use_GP_rot_kernel_fact.setChecked(False)
                self.use_GP_rot_kernel_fact.setEnabled(False)
                self.GP_rot_kernel_fact.setEnabled(False)
            else:
                self.use_GP_rot_kernel_fact.setChecked(True)
                self.use_GP_rot_kernel_fact.setEnabled(True)
                self.GP_rot_kernel_fact.setEnabled(True)           

            fit.link_RV_GP = [False,
                              self.use_tra_GP_rot_kernel_time_sc_link_to_RV.isChecked(),
                              self.use_tra_GP_rot_kernel_Per_link_to_RV.isChecked(),
                              self.use_tra_GP_rot_kernel_fact_link_to_RV.isChecked(),False,False]
            
        elif self.use_tra_GP_sho_kernel.isChecked():
                
            if self.use_tra_GP_sho_kernel_Q_link_to_RV.isChecked() and self.use_GP_sho_kernel.isChecked():
                self.use_GP_sho_kernel_Q.setChecked(False)
                self.use_GP_sho_kernel_Q.setEnabled(False)
                self.GP_sho_kernel_Q.setEnabled(False)
            else:
                self.use_GP_sho_kernel_Q.setChecked(True)
                self.use_GP_sho_kernel_Q.setEnabled(True)
                self.GP_sho_kernel_Q.setEnabled(True)                            
                
            if self.use_tra_GP_sho_kernel_omega_link_to_RV.isChecked() and self.use_GP_sho_kernel.isChecked():
                self.use_GP_sho_kernel_omega.setChecked(False)
                self.use_GP_sho_kernel_omega.setEnabled(False)
                self.GP_sho_kernel_omega.setEnabled(False)
            else:
                self.use_GP_sho_kernel_omega.setChecked(True)
                self.use_GP_sho_kernel_omega.setEnabled(True)
                self.GP_sho_kernel_omega.setEnabled(True)         

            fit.link_RV_GP = [False,
                              self.use_tra_GP_sho_kernel_Q_link_to_RV.isChecked(),
                              self.use_tra_GP_sho_kernel_omega_link_to_RV.isChecked(),
                              False,False,False]


        elif self.use_tra_GP_double_sho_kernel.isChecked():
                
            if self.use_tra_GP_double_sho_kernel_P_link_to_RV.isChecked() and self.use_GP_double_sho_kernel.isChecked():
                self.use_GP_double_sho_kernel_P.setChecked(False)
                self.use_GP_double_sho_kernel_P.setEnabled(False)
                self.GP_double_sho_kernel_P.setEnabled(False)
            else:
                self.use_GP_double_sho_kernel_P.setChecked(True)
                self.use_GP_double_sho_kernel_P.setEnabled(True)
                self.GP_double_sho_kernel_P.setEnabled(True)                            
                
            if self.use_tra_GP_double_sho_kernel_Q0_link_to_RV.isChecked() and self.use_GP_double_sho_kernel.isChecked():
                self.use_GP_double_sho_kernel_Q0.setChecked(False)
                self.use_GP_double_sho_kernel_Q0.setEnabled(False)
                self.GP_double_sho_kernel_Q0.setEnabled(False)
            else:
                self.use_GP_double_sho_kernel_Q0.setChecked(True)
                self.use_GP_double_sho_kernel_Q0.setEnabled(True)
                self.GP_double_sho_kernel_Q0.setEnabled(True)          

            if self.use_tra_GP_double_sho_kernel_dQ_link_to_RV.isChecked() and self.use_GP_double_sho_kernel.isChecked():
                self.use_GP_double_sho_kernel_dQ.setChecked(False)
                self.use_GP_double_sho_kernel_dQ.setEnabled(False)
                self.GP_double_sho_kernel_dQ.setEnabled(False)
            else:
                self.use_GP_double_sho_kernel_dQ.setChecked(True)
                self.use_GP_double_sho_kernel_dQ.setEnabled(True)
                self.GP_double_sho_kernel_dQ.setEnabled(True)          


            fit.link_RV_GP = [False,
                              self.use_tra_GP_double_sho_kernel_P_link_to_RV.isChecked(),
                              self.use_tra_GP_double_sho_kernel_Q0_link_to_RV.isChecked(),
                              self.use_tra_GP_double_sho_kernel_dQ_link_to_RV.isChecked(),False,False]


                        
        else:
                
            fit.link_RV_GP = [False, False, False,False,False,False]
            
            
                
################################ RV files #######################################################
    def showDialog_fortran_input_file(self):
        global fit, ses_list
 
        input_files = QtWidgets.QFileDialog.getOpenFileName(self, 'Open session', '', 'Data (*.init)', options=QtWidgets.QFileDialog.Option.DontUseNativeDialog)
        
        if str(input_files[0]) != '':
            fit_new=rv.signal_fit(str(input_files[0]), 'RVmod session',readinputfile=True)

            if len(ses_list) == 1:
                ses_list[0] = fit_new
                fit = fit_new
            else:
                ses_list.append(fit_new)

            self.session_list()
            self.update_use_from_session()
            self.init_fit()
            self.update_RV_file_buttons()


    def showDialog_RVbank_input_file(self):
        global fit, ses_list

        input_files = QtWidgets.QFileDialog.getOpenFileName(self, 'Open RVBank data', '', 'All (*.*);;Data (*.csv)', options=QtWidgets.QFileDialog.Option.DontUseNativeDialog)

        if str(input_files[0]) != '':

            choice = QtWidgets.QMessageBox.information(self, 'Warning!',
                                            "Do you want to split the %s to pre- and post- (if applicable)?"%self.data_str,
                                            QtWidgets.QMessageBox.StandardButton.No | QtWidgets.QMessageBox.StandardButton.Yes)  

            if choice == QtWidgets.QMessageBox.StandardButton.No:
                fit.add_RVbank_dataset(self.file_from_path(input_files[0]), str(input_files[0]), split = False)
            elif choice == QtWidgets.QMessageBox.StandardButton.Yes:            
                fit.add_RVbank_dataset(self.file_from_path(input_files[0]), str(input_files[0]), split = True)
            else:
                return

            fit.type_fit["RV"] = True
            fit.type_fit["Transit"] = False
            self.check_type_fit()
            self.mute_boxes()

            self.init_fit()
            self.update_use_from_session()
            self.update_use()
            self.update_params()
            self.update_RV_file_buttons()
            self.update_act_file_buttons()

    def showDialog_MAROONX_input_file(self):
        global fit, ses_list

        input_files = QtWidgets.QFileDialog.getOpenFileName(self, 'Open MAROON-X data', '', 'All (*.*);;Data (*.csv)', options=QtWidgets.QFileDialog.Option.DontUseNativeDialog)

        if str(input_files[0]) != '':

            choice = QtWidgets.QMessageBox.information(self, 'Warning!',
                                            "Do you want to split the %s to files with different offsets (if applicable)?"%self.data_str,
                                            QtWidgets.QMessageBox.StandardButton.No | QtWidgets.QMessageBox.StandardButton.Yes)  

            if choice == QtWidgets.QMessageBox.StandardButton.No:
                fit.add_MAROONX_dataset(self.file_from_path(input_files[0]), str(input_files[0]), split = False)
            elif choice == QtWidgets.QMessageBox.StandardButton.Yes:            
                fit.add_MAROONX_dataset(self.file_from_path(input_files[0]), str(input_files[0]), split = True)
            else:
                return

            fit.type_fit["RV"] = True
            fit.type_fit["Transit"] = False
            self.check_type_fit()
            self.mute_boxes()
            self.update_bounds()
            self.init_fit()
            self.update_use_from_session()
            self.update_use()
            self.update_params()
            self.update_RV_file_buttons()
            self.update_act_file_buttons()




 
 

    def apply_rvs_data_options(self):
        global fit
        but_ind = self.buttonGroup_apply_rvs_data_options.checkedId()

        print(but_ind)
        if self.rvs_opt[but_ind-1].isChecked() == True and len(fit.rv_data_sets[but_ind-1]) != 0:
            if len(self.RvsWindow.flux_o_c) != 0:
                fit.rv_data_sets[but_ind-1][0] = dill.copy(self.RvsWindow.t_store[but_ind-1])
                fit.rv_data_sets[but_ind-1][1] = dill.copy(self.RvsWindow.flux_o_c_store[but_ind-1])
                fit.rv_data_sets[but_ind-1][2] = dill.copy(self.RvsWindow.flux_err_o_c_store[but_ind-1])
                fit.rv_data_sets[but_ind-1][3] = dill.copy(self.RvsWindow.flux_idset_store[but_ind-1])
                #fit.rv_data_sets[but_ind-1][4] = dill.copy(self.RvsWindow.flux_o_c_store[but_ind-1])

                #fit.act_data_sets[but_ind-1][5] = dill.copy(self.RvsWindow.flux_store[but_ind-1])
                #fit.act_data_sets[but_ind-1][6] = dill.copy(self.RvsWindow.flux_err_store[but_ind-1])
               # fit.act_data_sets[but_ind-1][7] = dill.copy(self.RvsWindow.trend_store[but_ind-1])
                
                #### CHECK FOR NAN ENTRIES AND REMOVE IF ANY, BELOW !!!! ####
               
                bjd    = dill.copy(self.RvsWindow.t_store[but_ind-1])
                flux   = dill.copy(self.RvsWindow.flux_o_c_store[but_ind-1])               
                flux_e = dill.copy(self.RvsWindow.flux_err_o_c_store[but_ind-1])                   
                flux_idset = dill.copy(self.RvsWindow.flux_idset_store[but_ind-1])                   

                #print(but_ind-1)

                for j in range(3): 
                    fit.rv_data_sets[but_ind-1][j] = fit.rv_data_sets[but_ind-1][j][np.isfinite(bjd) & 
                                                                                      np.isfinite(flux) & 
                                                                                      np.isfinite(flux_e)]
                
                if len(bjd) != len(fit.rv_data_sets[but_ind-1][0]):
                    print("WARNING! Something is wrong for these epochs:")          
                    for v in range(len(bjd)):
                        if np.isnan(flux[v]) or np.isnan(flux_e[v]):
                            print("%s  %s   %s"%(bjd[v], flux[v], flux_e[v]))
                    print("These data are not included!")


        #fit.n_rvdata = sum([len(self.rv_data_sets[x][0]) for x in range(20) if len(self.rv_data_sets[x]) != 0])

        fit.model_saved = False 
 
        self.tabWidget_helper.setCurrentWidget(self.tab_info)
        fit.type_fit["RV"] = True
        fit.type_fit["Transit"] = False
        self.check_type_fit()
        self.mute_boxes()
        
        self.init_fit()




    def apply_act_data_options(self):
        global fit
        but_ind = self.buttonGroup_apply_act_data_options.checkedId()
 

        if self.act_opt[but_ind-1].isChecked() == True and len(fit.act_data_sets[but_ind-1]) != 0:
            #rv.transit_data_norm(fit,  file_n = but_ind-1, norm = True)
            if len(self.ActivityWindow.flux_o_c) != 0:
                fit.act_data_sets[but_ind-1][0] = dill.copy(self.ActivityWindow.t_store[but_ind-1])
                fit.act_data_sets[but_ind-1][1] = dill.copy(self.ActivityWindow.flux_o_c_store[but_ind-1])
                fit.act_data_sets[but_ind-1][2] = dill.copy(self.ActivityWindow.flux_err_o_c_store[but_ind-1])
                fit.act_data_sets[but_ind-1][3] = dill.copy(self.ActivityWindow.flux_o_c_store[but_ind-1])
                fit.act_data_sets[but_ind-1][4] = dill.copy(self.ActivityWindow.flux_o_c_store[but_ind-1])

                fit.act_data_sets[but_ind-1][5] = dill.copy(self.ActivityWindow.flux_store[but_ind-1])
                fit.act_data_sets[but_ind-1][6] = dill.copy(self.ActivityWindow.flux_err_store[but_ind-1])
                fit.act_data_sets[but_ind-1][7] = dill.copy(self.ActivityWindow.trend_store[but_ind-1])
                
                #### CHECK FOR NAN ENTRIES AND REMOVE IF ANY, BELOW !!!! ####
               
                bjd    = dill.copy(self.ActivityWindow.t_store[but_ind-1])
                flux   = dill.copy(self.ActivityWindow.flux_o_c_store[but_ind-1])               
                flux_e = dill.copy(self.ActivityWindow.flux_o_c_store[but_ind-1])                   
                #print(but_ind-1)
                for j in range(8): 
                    fit.act_data_sets[but_ind-1][j] = fit.act_data_sets[but_ind-1][j][np.isfinite(bjd) & 
                                                                                      np.isfinite(flux) & 
                                                                                      np.isfinite(flux_e)]
                
                if len(bjd) != len(fit.act_data_sets[but_ind-1][0]):
                    print("WARNING! Something is wrong for these epochs:")          
                    for v in range(len(bjd)):
                        if np.isnan(flux[v]) or np.isnan(flux_e[v]):
                            print("%s  %s   %s"%(bjd[v], flux[v], flux_e[v]))
                    print("These data are not included!")

                #### CHECK IS OVER ####                  
            else:
                rv.sigma_clip(fit, type = 'act', sigma_clip = None, remove_mean =  False, file_n = but_ind-1)
        else:
            rv.sigma_clip(fit, type = 'act', sigma_clip = None, remove_mean =  False, file_n = but_ind-1)

        self.tabWidget_helper.setCurrentWidget(self.tab_info)
        self.update_activity_data_plots(self.comboBox_act_data.currentIndex())
        self.update_activity_gls_plots(but_ind-1)


     
    def apply_tra_dilution(self):
        global fit
        but_ind = self.buttonGroup_apply_tra_dilution.checkedId()


        if self.tra_dilution[but_ind-1][1].isChecked() == True  and len(fit.tra_data_sets[but_ind-1]) != 0:
            fit.tra_dil[but_ind-1] = self.tra_dilution[but_ind-1][0].value()
            fit.tra_data_sets[but_ind-1][8] = fit.tra_dil[but_ind-1]
            
            
        elif self.tra_dilution[but_ind-1][1].isChecked() == False  and len(fit.tra_data_sets[but_ind-1]) != 0:
            fit.tra_dil[but_ind-1] = 1.0
            fit.tra_data_sets[but_ind-1][8] = fit.tra_dil[but_ind-1]
     
        self.update_veiw()

    def apply_tra_data_options(self):
        global fit
        but_ind = self.buttonGroup_apply_tra_data_options.checkedId()
 

        if self.tra_norm[but_ind-1].isChecked() == True and len(fit.tra_data_sets[but_ind-1]) != 0:
            #rv.transit_data_norm(fit,  file_n = but_ind-1, norm = True)
            if len(self.DetrendWindow.flux_o_c) != 0:
                fit.tra_data_sets[but_ind-1][0] = dill.copy(self.DetrendWindow.t_store[but_ind-1])
                fit.tra_data_sets[but_ind-1][1] = dill.copy(self.DetrendWindow.flux_o_c_store[but_ind-1])
                fit.tra_data_sets[but_ind-1][2] = dill.copy(self.DetrendWindow.flux_err_o_c_store[but_ind-1])
                fit.tra_data_sets[but_ind-1][3] = dill.copy(self.DetrendWindow.airmass_store[but_ind-1])
                fit.tra_data_sets[but_ind-1][4] = dill.copy(self.DetrendWindow.flux_o_c_store[but_ind-1])

                fit.tra_data_sets[but_ind-1][5] = dill.copy(self.DetrendWindow.flux_store[but_ind-1])
                fit.tra_data_sets[but_ind-1][6] = dill.copy(self.DetrendWindow.flux_err_store[but_ind-1])
                fit.tra_data_sets[but_ind-1][7] = dill.copy(self.DetrendWindow.trend_store[but_ind-1])
                
                #### CHECK FOR NAN ENTRIES AND REMOVE IF ANY, BELOW !!!! ####
               
                bjd    = dill.copy(self.DetrendWindow.t_store[but_ind-1])
                flux   = dill.copy(self.DetrendWindow.flux_o_c_store[but_ind-1])               
                flux_e = dill.copy(self.DetrendWindow.flux_err_o_c_store[but_ind-1])                   
                
                for j in range(8): 
                    fit.tra_data_sets[but_ind-1][j] = fit.tra_data_sets[but_ind-1][j][np.isfinite(bjd) & np.isfinite(flux) & np.isfinite(flux_e)]
                
                if len(bjd) != len(fit.tra_data_sets[but_ind-1][0]):
                    print("WARNING! Something is wrong for these epochs:")          
                    for v in range(len(bjd)):
                        if np.isnan(flux[v]) or np.isnan(flux_e[v]):
                            print("%s  %s   %s"%(bjd[v], flux[v], flux_e[v]))
                    print("These data are not included!")

                #### CHECK IS OVER ####                  

            else:
                rv.transit_data_norm(fit,  file_n = but_ind-1, norm = False)
        else:
            rv.transit_data_norm(fit,  file_n = but_ind-1, norm = False)

        self.tabWidget_helper.setCurrentWidget(self.tab_info)
        self.update_veiw()



    def showDialog_RV_input_file(self):
        global fit

        but_ind = self.buttonGroup_add_RV_data.checkedId()   
        input_files = QtWidgets.QFileDialog.getOpenFileName(self, 'Open %s'%self.data_str, '', 'All (*.*);;Data (*.vels)', options=QtWidgets.QFileDialog.Option.DontUseNativeDialog)

        if str(input_files[0]) != '':

            ######## This must be in RV_mod once the obsolite code is removed ##########
            try:
                rv_JD_in       = np.genfromtxt("%s"%(str(input_files[0])),skip_header=0, unpack=True,skip_footer=0, usecols = [0])
                rv_data_in     = np.genfromtxt("%s"%(str(input_files[0])),skip_header=0, unpack=True,skip_footer=0, usecols = [1])
                rv_data_sig_in = np.genfromtxt("%s"%(str(input_files[0])),skip_header=0, unpack=True,skip_footer=0, usecols = [2])
                #print(rv_JD_in,rv_data_in, rv_data_sig_in)

                if len(rv_JD_in) != len(rv_data_in) != len(rv_data_sig_in):
                    print("Something is wrong with your %s file! Please provide a valid %s file that contains: BJD RV [m/s] sigma_RV [m/s] "%(self.data_str,self.data_str))
                    return
                 
                rv_JD        = rv_JD_in[      np.isfinite(rv_JD_in) & np.isfinite(rv_data_in) & np.isfinite(rv_data_sig_in)]
                rv_data      = rv_data_in[    np.isfinite(rv_JD_in) & np.isfinite(rv_data_in) & np.isfinite(rv_data_sig_in)]
                rv_data_sig  = rv_data_sig_in[np.isfinite(rv_JD_in) & np.isfinite(rv_data_in) & np.isfinite(rv_data_sig_in)]
                               
                if len(rv_JD_in) ==0 or len(rv_JD_in) != len(rv_JD):
                    print("Something is wrong with your %s file! Perhaps some not all entries are numeric? Please provide a valid %s file that contains: BJD RV [m/s] sigma_RV [m/s] "%(self.data_str,self.data_str))
                    return
     
            except:
                print("Something is wrong with your %s file! Please provide a valid %s file that contains: BJD RV [m/s] sigma_RV [m/s] "%(self.data_str,self.data_str))
                return 
            ################################################################################


            fit.add_dataset(self.file_from_path(input_files[0]), str(input_files[0]),0.0,1.0, index =but_ind-1)
            #### new stuf ####
            #fit.add_rv_dataset('test', str(input_files[0]),rv_idset =but_ind-1)
            ##################
            
            fit.type_fit["RV"] = True
            fit.type_fit["Transit"] = False
            self.check_type_fit()
            self.mute_boxes()
            
            self.init_fit()
            self.update_use_from_session()
            self.update_use()
            self.update_params()
            self.update_RV_file_buttons()

        self.plot_tabs.setCurrentWidget(self.tab_timeseries_RV)


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
        #fit.remove_rv_dataset(but_ind -1)
        #### new stuf ####

        fit.type_fit["RV"] = True
        fit.type_fit["Transit"] = False
        fit.type_fit["AST"] = False
        self.check_type_fit()
        self.mute_boxes()

        self.init_fit()
        #self.update_view()
        self.update_use_from_session()
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

        for i in range(20):
            if len(fit.rv_data_sets[i]) !=0:
                self.buttonGroup_add_RV_data.button(i+1).setStyleSheet("color: %s;"%dill.copy(fit.rvs_colors[i]))
                self.buttonGroup_remove_RV_data.button(i+1).setStyleSheet("color: %s;"%dill.copy(fit.rvs_colors[i]))
                font.setPointSize(9)
                #self.buttonGroup_add_RV_data.button(i+1).setText(fit.filelist.files[i].name) 
                self.buttonGroup_add_RV_data.button(i+1).setText(fit.rv_data_sets[i][4]) 

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
        input_files = QtWidgets.QFileDialog.getOpenFileName(self, 'Open Transit data', '', 'All (*.*);;Data (*.tran)', options=QtWidgets.QFileDialog.Option.DontUseNativeDialog)

        if str(input_files[0]) != '':
            
            
            if input_files[0].endswith(".fits"):
                
                but_n = self.tess_pdc_dialog.get_radio()
                
                if but_n ==1:
                    PDC = False
                elif but_n ==2: 
                    PDC = True
                else:
                    PDC = False
                fit.add_transit_dataset('test', str(input_files[0]),tra_idset =but_ind-1, PDC = PDC)
            else:
                fit.add_transit_dataset('test', str(input_files[0]),tra_idset =but_ind-1)

            self.update_use_from_session()
            self.update_use()
            self.update_gui_params()

            #self.radioButton_transit.setChecked(True)
            #self.worker_transit_fitting(ff=0)

            #self.init_fit()
            
            fit.type_fit["RV"] = False
            fit.type_fit["Transit"] = True
            fit.type_fit["AST"] = False
            self.check_type_fit()
            self.mute_boxes()
            
            self.check_use_tra_GP()
            
            self.update_params()
            self.update_tra_file_buttons()
            #self.buttonGroup_transit_data.button(but_ind).setText(self.file_from_path(input_files[0]))
            self.plot_tabs.setCurrentWidget(self.tab_timeseries_tra)
 

    def remove_tra_file(self):
        global fit

        but_ind = self.buttonGroup_remove_transit_data.checkedId()   
        fit.remove_transit_dataset(but_ind -1)
       # self.init_fit()         
        self.update_use_from_session()   
      #  self.update_use()
      #  self.update_gui_params()
     #   self.update_params()
        fit.type_fit["RV"] = False
        fit.type_fit["Transit"] = True
        fit.type_fit["AST"] = False
        self.check_type_fit()
        self.mute_boxes()
        self.update_tra_file_buttons()


    def update_tra_file_buttons(self):
        global fit, colors          



        for i in range(20):
            if len(fit.tra_data_sets[i]) != 0:
                
                
                self.buttonGroup_transit_data.button(i+1).setStyleSheet("color: %s;"%fit.tra_colors[i])
                self.buttonGroup_remove_transit_data.button(i+1).setStyleSheet("color: %s;"%fit.tra_colors[i])
                self.buttonGroup_transit_data.button(i+1).setText(fit.tra_data_sets[i][-1])

            else:
                self.buttonGroup_transit_data.button(i+1).setStyleSheet("")
                self.buttonGroup_remove_transit_data.button(i+1).setStyleSheet("")
                self.buttonGroup_transit_data.button(i+1).setText("data %s"%(i+1))

                #"background-color: #333399;""background-color: yellow;" "selection-color: yellow;"  "selection-background-color: blue;")               
 
        if len([x for x in range(20) if len(fit.tra_data_sets[x]) != 0]) == 0:
            self.update_transit_plots()
        else:
            self.worker_transit_fitting(ff=0)
#        self.update_transit_plots()
 

################################ activity files #######################################################
        
    def showDialog_act_input_file(self):
        global fit

        but_ind = self.buttonGroup_activity_data.checkedId()   
        input_files = QtWidgets.QFileDialog.getOpenFileName(self, 'Open Activity data', '', 'All (*.*);;Data (*.act)', options=QtWidgets.QFileDialog.Option.DontUseNativeDialog)

        if str(input_files[0]) != '':

            fit.add_act_dataset('test', str(input_files[0]),act_idset =but_ind-1)

            self.update_act_file_buttons()
            self.buttonGroup_activity_data.button(but_ind).setText(self.file_from_path(input_files[0]))

            self.plot_tabs.setCurrentWidget(self.tab_timeseries_act)
            self.comboBox_act_data.setCurrentIndex(but_ind-1)
            self.comboBox_act_data_gls.setCurrentIndex(but_ind-1)

            self.update_activity_data_plots(self.comboBox_act_data.currentIndex())
            self.update_activity_gls_plots(self.comboBox_act_data_gls.currentIndex())


    def remove_act_file(self):
        global fit

        but_ind = self.buttonGroup_remove_activity_data.checkedId()   
        fit.remove_act_dataset(but_ind -1)
       # self.init_fit()         
      #  self.update_use_from_session()   
      #  self.update_use()
      #  self.update_gui_params()
     #   self.update_params()
        self.update_act_file_buttons()

    def update_act_file_buttons(self):
        global fit, colors          

        for i in range(20):
            if len(fit.act_data_sets[i]) != 0:
                self.buttonGroup_activity_data.button(i+1).setStyleSheet("color: %s;"%dill.copy(fit.act_colors[i]))
                self.buttonGroup_remove_activity_data.button(i+1).setStyleSheet("color: %s;"%dill.copy(fit.act_colors[i]))
                self.buttonGroup_activity_data.button(i+1).setText(fit.act_data_sets[i][-1])

            else:
                self.buttonGroup_activity_data.button(i+1).setStyleSheet("")
                self.buttonGroup_remove_activity_data.button(i+1).setStyleSheet("")
                self.buttonGroup_activity_data.button(i+1).setText("data %s"%(i+1))

                #"background-color: #333399;""background-color: yellow;" "selection-color: yellow;"  "selection-background-color: blue;")               
        self.init_correlations_combo()




################################ TTV files #######################################################
        
    def showDialog_ttv_input_file(self):
        global fit

        but_ind = self.buttonGroup_ttv_data.checkedId()   
        input_files = QtWidgets.QFileDialog.getOpenFileName(self, 'Open TTV data', '', 'All (*.*);;Data (*.ttv)', options=QtWidgets.QFileDialog.Option.DontUseNativeDialog)
        
        
        planet_N     = self.ttv_data_to_planet[but_ind-1].value()
        use_planet_N = self.use_ttv_data_to_planet[but_ind-1].isChecked()

        if str(input_files[0]) != '':
 
            fit.add_ttv_dataset('test', str(input_files[0]), ttv_idset =but_ind-1, planet =planet_N, use = use_planet_N )
            self.init_fit()
            #self.update_use_from_session()
            #self.update_use()
            #self.update_params()
            
            
            fit.type_fit["RV"] = False
            fit.type_fit["Transit"] = False
            fit.type_fit["TTV"] = True
            fit.type_fit["AST"] = False
            self.check_type_fit()
            self.mute_boxes()
            
            self.update_params()
            self.update_ttv_file_buttons()
 
            self.plot_tabs.setCurrentWidget(self.tab_timeseries_ttv)


    def remove_ttv_file(self):
        global fit

        but_ind = self.buttonGroup_remove_ttv_data.checkedId()   
        fit.remove_ttv_dataset(but_ind -1)
       # self.init_fit()

        fit.type_fit["RV"] = False
        fit.type_fit["Transit"] = False
        fit.type_fit["TTV"] = True
        self.check_type_fit()
        self.mute_boxes()
        self.update_ttv_file_buttons()
 

    def update_ttv_file_buttons(self):
        global fit, colors          

        for i in range(10):
            if len(fit.ttv_data_sets[i]) != 0:
                self.buttonGroup_ttv_data.button(i+1).setStyleSheet("color: %s;"%dill.copy(fit.colors[i]))
                self.buttonGroup_remove_ttv_data.button(i+1).setStyleSheet("color: %s;"%dill.copy(fit.colors[i]))
                self.buttonGroup_ttv_data.button(i+1).setText(fit.ttv_data_sets[i][5])
                #self.ttv_data_to_planet[i].setValue(fit.ttv_data_sets[i][3])
                #self.use_ttv_data_to_planet[i].setChecked(bool(fit.ttv_data_sets[i][4]))
               # print(i, bool(fit.ttv_data_sets[i][4]))

            else:
                self.buttonGroup_ttv_data.button(i+1).setStyleSheet("")
                self.buttonGroup_remove_ttv_data.button(i+1).setStyleSheet("")
                self.buttonGroup_ttv_data.button(i+1).setText("data %s"%(i+1))
               # self.ttv_data_to_planet[i].setValue(1)
                #self.use_ttv_data_to_planet[i].setChecked(False)
                #self.set_ttv_dataset_to_planet()

                #"background-color: #333399;""background-color: yellow;" "selection-color: yellow;"  "selection-background-color: blue;")               
        #self.init_correlations_combo()



    def ttv_dataset_to_planet(self):
        
        for i in range(10):
            if len(fit.ttv_data_sets[i]) ==0:
                continue
            else:
                fit.ttv_data_sets[i][3] = self.ttv_data_to_planet[i].value()

    def use_ttv_dataset_to_planet(self):

        for i in range(10):
            if len(fit.ttv_data_sets[i]) ==0:
                continue
            else:
                fit.ttv_data_sets[i][4] = self.use_ttv_data_to_planet[i].isChecked()
 

    def set_ttv_dataset_to_planet(self):

        for i in range(10):
            if len(fit.ttv_data_sets[i]) ==0:
                self.ttv_data_to_planet[i].setValue(1)
                self.use_ttv_data_to_planet[i].setChecked(False)
                continue
            else:
                #print(fit.ttv_data_sets[i][3],fit.ttv_data_sets[i][4])                
                self.ttv_data_to_planet[i].setValue(fit.ttv_data_sets[i][3])
                self.use_ttv_data_to_planet[i].setChecked(bool(fit.ttv_data_sets[i][4]))




################################ Astr. files #######################################################
        
    def showDialog_ast_input_file(self):
        global fit

        but_ind = self.buttonGroup_ast_data.checkedId()   
        input_files = QtWidgets.QFileDialog.getOpenFileName(self, 'Open Astrometry data', '', 'All (*.*);;Data (*.ast)', options=QtWidgets.QFileDialog.Option.DontUseNativeDialog)
        
        
        planet_N     = self.ast_data_to_planet[but_ind-1].value()
        use_planet_N = self.use_ast_data_to_planet[but_ind-1].isChecked()

        if str(input_files[0]) != '':
 
            fit.add_ast_dataset('test', str(input_files[0]), ast_idset =but_ind-1, planet =planet_N, use = use_planet_N )
            self.init_fit()
            #self.update_use_from_session()
            #self.update_use()
            #self.update_params()
            
            
            fit.type_fit["RV"] = False
            fit.type_fit["Transit"] = False
            fit.type_fit["TTV"] = False
            fit.type_fit["AST"] = True
            self.check_type_fit()
            self.mute_boxes()
            
            self.update_params()
            self.update_ast_file_buttons()
 
            self.plot_tabs.setCurrentWidget(self.tab_timeseries_ast)


    def remove_ast_file(self):
        global fit

        but_ind = self.buttonGroup_remove_ast_data.checkedId()   
        fit.remove_ast_dataset(but_ind -1)
       # self.init_fit()

        fit.type_fit["RV"] = False
        fit.type_fit["Transit"] = False
        fit.type_fit["TTV"] = False
        fit.type_fit["AST"] = True

        self.check_type_fit()
        self.mute_boxes()
        self.update_ast_file_buttons()
 

    def update_ast_file_buttons(self):
        global fit, colors          

        for i in range(10):
            if len(fit.ast_data_sets[i]) != 0:
                self.buttonGroup_ast_data.button(i+1).setStyleSheet("color: %s;"%dill.copy(fit.colors[i]))
                self.buttonGroup_remove_ast_data.button(i+1).setStyleSheet("color: %s;"%dill.copy(fit.colors[i]))
                self.buttonGroup_ast_data.button(i+1).setText(fit.ast_data_sets[i][7])
                #self.ast_data_to_planet[i].setValue(fit.ast_data_sets[i][3])
                #self.use_ast_data_to_planet[i].setChecked(bool(fit.ast_data_sets[i][4]))
               # print(i, bool(fit.ast_data_sets[i][4]))

            else:
                self.buttonGroup_ast_data.button(i+1).setStyleSheet("")
                self.buttonGroup_remove_ast_data.button(i+1).setStyleSheet("")
                self.buttonGroup_ast_data.button(i+1).setText("data %s"%(i+1))
               # self.ast_data_to_planet[i].setValue(1)
                #self.use_ast_data_to_planet[i].setChecked(False)
                #self.set_ast_dataset_to_planet()

                #"background-color: #333399;""background-color: yellow;" "selection-color: yellow;"  "selection-background-color: blue;")               
        #self.init_correlations_combo()



    def ast_dataset_to_planet(self):
        
        for i in range(10):
            if len(fit.ast_data_sets[i]) ==0:
                continue
            else:
                fit.ast_data_sets[i][5] = self.ast_data_to_planet[i].value()

    def use_ast_dataset_to_planet(self):

        for i in range(10):
            if len(fit.ast_data_sets[i]) ==0:
                continue
            else:
                fit.ast_data_sets[i][6] = self.use_ast_data_to_planet[i].isChecked()
 

    def set_ast_dataset_to_planet(self):

        for i in range(10):
            if len(fit.ast_data_sets[i]) ==0:
                self.ast_data_to_planet[i].setValue(1)
                self.use_ast_data_to_planet[i].setChecked(False)
                continue
            else:
                #print(fit.ast_data_sets[i][3],fit.ast_data_sets[i][4])                
                self.ast_data_to_planet[i].setValue(fit.ast_data_sets[i][5])
                self.use_ast_data_to_planet[i].setChecked(bool(fit.ast_data_sets[i][6]))


################################ Astr. files (Hipp/Gaia) #####################################

    def check_std_ast_solution(self):
        global fit      
        
        for i in range(1): #TBD fixed!!!!
            #print("TEST",fit.ast_data_sets_hipp_gaia[i]) 
            if len(fit.ast_data_sets_hipp_gaia[i]) != 0:

                try:
                    self.ast_std_label_alpha.setText("%.3f +/- %.3f "%(fit.ast_data_sets_hipp_gaia[i][9]["RAdeg"],fit.ast_data_sets_hipp_gaia[i][9]["e_RA"]))                    
                    self.ast_std_label_delta.setText("%.3f +/- %.3f "%(fit.ast_data_sets_hipp_gaia[i][9]["DEdeg"],fit.ast_data_sets_hipp_gaia[i][9]["e_DE"]))                    
                    self.ast_std_label_pi.setText("%.3f +/- %.3f "%(fit.ast_data_sets_hipp_gaia[i][9]["Plx"],fit.ast_data_sets_hipp_gaia[i][9]["e_Plx"]))
                    self.ast_std_label_mu_alpha.setText("%.3f +/- %.3f "%(fit.ast_data_sets_hipp_gaia[i][9]["pm_RA"],fit.ast_data_sets_hipp_gaia[i][9]["e_pmRA"]))                    
                    self.ast_std_label_mu_delta.setText("%.3f +/- %.3f "%(fit.ast_data_sets_hipp_gaia[i][9]["pm_DE"],fit.ast_data_sets_hipp_gaia[i][9]["e_pmDE"]))             
                except:
                    print("No Hipp./Gaia std. solution found in the loaded memory! Make sure the input astrometric file is following the requested format (see READ ME button)")

            else:                       
                self.ast_std_label_alpha.setText("%.3f +/- %.3f "%(0.000,0.000))                    
                self.ast_std_label_delta.setText("%.3f +/- %.3f "%(0.000,0.000))                    
                self.ast_std_label_pi.setText("%.3f +/- %.3f "%(0.000,0.000))
                self.ast_std_label_mu_alpha.setText("%.3f +/- %.3f "%(0.000,0.000))                    
                self.ast_std_label_mu_delta.setText("%.3f +/- %.3f "%(0.000,0.000))                    
                
                               
    def showDialog_ast_input_file_2(self):
        global fit

        but_ind = self.buttonGroup_ast_data_2.checkedId()   
        input_files = QtWidgets.QFileDialog.getOpenFileName(self, 'Open Astrometry (Hipp/Gaia) data', '', 'All (*.*);;Data (*.ast)', options=QtWidgets.QFileDialog.Option.DontUseNativeDialog)
        
        
        planet_N     = self.ast_data_to_planet_2[but_ind-1].value()
        use_planet_N = self.use_ast_data_to_planet_2[but_ind-1].isChecked()

        if str(input_files[0]) != '':
 
            fit.add_ast_dataset_hipp_gaia('test', str(input_files[0]), ast_idset =but_ind-1, planet =planet_N, use = use_planet_N )
 
            
            self.check_std_ast_solution()
            
            self.init_fit()
            #self.update_use_from_session()
            #self.update_use()
            #self.update_params()
            
            
            fit.type_fit["RV"] = False
            fit.type_fit["Transit"] = False
            fit.type_fit["TTV"] = False
            fit.type_fit["AST"] = True
            self.check_type_fit()
            self.mute_boxes()
            
            self.update_params()
            self.update_ast_file_buttons_2()
 
            self.plot_tabs.setCurrentWidget(self.tab_timeseries_ast)


    def remove_ast_file_2(self):
        global fit

        but_ind = self.buttonGroup_remove_ast_data_2.checkedId()   
        fit.remove_ast_dataset_hipp_gaia(but_ind -1)
       # self.init_fit()

        self.check_std_ast_solution()

        fit.type_fit["RV"] = False
        fit.type_fit["Transit"] = False
        fit.type_fit["TTV"] = False
        fit.type_fit["AST"] = True

        self.check_type_fit()
        self.mute_boxes()
        self.update_ast_file_buttons_2()
 

    def update_ast_file_buttons_2(self):
        global fit, colors          

        for i in range(2):
            if len(fit.ast_data_sets_hipp_gaia[i]) != 0:
                self.buttonGroup_ast_data_2.button(i+1).setStyleSheet("color: %s;"%dill.copy(fit.colors[i]))
                self.buttonGroup_remove_ast_data_2.button(i+1).setStyleSheet("color: %s;"%dill.copy(fit.colors[i]))
                self.buttonGroup_ast_data_2.button(i+1).setText(fit.ast_data_sets_hipp_gaia[i][10])
                #self.ast_data_to_planet[i].setValue(fit.ast_data_sets[i][3])
                #self.use_ast_data_to_planet[i].setChecked(bool(fit.ast_data_sets[i][4]))
               # print(i, bool(fit.ast_data_sets[i][4]))

            else:
                self.buttonGroup_ast_data_2.button(i+1).setStyleSheet("")
                self.buttonGroup_remove_ast_data_2.button(i+1).setStyleSheet("")
                self.buttonGroup_ast_data_2.button(i+1).setText("data %s"%(i+1))
               # self.ast_data_to_planet[i].setValue(1)
                #self.use_ast_data_to_planet[i].setChecked(False)
                #self.set_ast_dataset_to_planet()

                #"background-color: #333399;""background-color: yellow;" "selection-color: yellow;"  "selection-background-color: blue;")               
        #self.init_correlations_combo()



    def ast_dataset_to_planet_2(self):
        
        for i in range(2):
            if len(fit.ast_data_sets_hipp_gaia[i]) ==0:
                continue
            else:
                fit.ast_data_sets_hipp_gaia[i][7] = self.ast_data_to_planet_2[i].value()

    def use_ast_dataset_to_planet_2(self):

        for i in range(2):
            if len(fit.ast_data_sets_hipp_gaia[i]) ==0:
                continue
            else:
                fit.ast_data_sets_hipp_gaia[i][8] = self.use_ast_data_to_planet_2[i].isChecked()
 
 

    def set_ast_dataset_to_planet_2(self):

        for i in range(2):
            if len(fit.ast_data_sets_hipp_gaia[i]) ==0:
                self.ast_data_to_planet_2[i].setValue(1)
                self.use_ast_data_to_planet_2[i].setChecked(False)
                continue
            else:
                #print(fit.ast_data_sets[i][3],fit.ast_data_sets[i][4])                
                self.ast_data_to_planet_2[i].setValue(fit.ast_data_sets_hipp_gaia[i][7])
                self.use_ast_data_to_planet_2[i].setChecked(bool(fit.ast_data_sets_hipp_gaia[i][8]))

##################################### Various ################################# 


    def init_fit(self, update_error=True): 
        global fit
        
        
        # A hack when .ses are imported from the Example directory... TBFixed
        fit.cwd = os.getcwd()

        #self.fit_dispatcher(init=True)


        self.check_model_params()
        self.check_mcmc_params()
        self.check_nested_params()
        self.check_settings()

        minimize_fortran=True

        if fit.model_saved == False or len(fit.fit_results.rv_model.jd) != fit.n_rvdata:
            
            fit.model_npoints = self.points_to_draw_model.value()
            fit.model_max = self.model_max_range.value()
            fit.model_min = self.model_min_range.value()
            fit.init_fit = True            
            fit.fitting(outputfiles=[1,1,1], minimize_fortran=minimize_fortran, 
                        doGP=fit.doGP,  fortran_kill=self.dyn_model_to_kill.value(), 
                        timeout_sec=self.master_timeout.value(), minimize_loglik=True,amoeba_starts=0, 
                        print_stat=False, eps=self.dyn_model_accuracy.value(), dt=self.time_step_model.value(), 
                        return_flag=True, npoints=fit.model_npoints, model_max = fit.model_max, model_min = fit.model_min)


            for i in range(9):
                if not bool(fit.use_planet[i]):
                    continue
                try:

                    rv.phase_RV_planet_signal(fit,i)
                except:
                    print("Old session (pre Ver. 0.86)? You are using the old phased plot version!")
                    rv.phase_RV_planet_signal_old(fit,i)  
        
        fit.init_fit = False            
        
        #print(fit.fit_results.rv_model.rv_err)
        self.update_labels()
#        self.update_params()

        self.update_gui_params()
        #print("FFF")
       # if minimize_loglik == True:
        if update_error ==  True:
            self.update_errors() 
            
        self.update_a_mass() 
        
       # self.run_gls()
       # self.run_gls_o_c()
        
        self.update_plots() 
        self.update_transit_plots() 
        self.plot_evol_all()

        self.jupiter_push_vars() 


    def add_jitter(self, errors, ind):
        global fit
 
        errors_with_jitt = np.array([np.sqrt(errors[i]**2 + fit.jitt[ii]**2)  for i,ii in enumerate(ind)])
        return errors_with_jitt






############ MLP ##############################      
       
    def worker_mlp_complete(self, resid = False):
        global fit  

        #start_time = time.time()   

        #if resid == False:     
        #    self.update_RV_MLP_plots() 
        #else:
        #    self.update_mlp_o_c_plots() 

        self.update_RV_MLP_plots()

        self.statusBar().showMessage('')
 
        self.jupiter_push_vars()
        self.calc_MLP.setEnabled(True)
       # self.calc_MLP_o_c.setEnabled(True)

 
    def worker_mlp(self, resid = False):
        global fit  

        self.calc_MLP.setEnabled(False)
        #self.calc_MLP_o_c.setEnabled(False)


        if self.calp_MLP_o_c.isChecked():
            resid = True
        else:
            resid = False
 
        if self.mls_o_c_GP.isChecked():
            resid_GP = True
        else:
            resid_GP = False           
            
     
            
        
        #if z <= 0:
        #    choice = QtWidgets.QMessageBox.information(self, 'Warning!',
        #    "Not possible to look for planets if there are no transit data loaded. Please add your transit data first. Okay?", QtWidgets.QMessageBox.StandardButton.Ok)      
        #    self.calc_TLS.setEnabled(True)         
        #    return   

        self.statusBar().showMessage('Running MLP .... This might take some time!')                 
        worker_mlp_wk = Worker(lambda:  self.mlp_search(resid = resid, resid_GP=resid_GP) )# Any other args, kwargs are passed to the run  
 
        worker_mlp_wk.signals.finished.connect(lambda:  self.worker_mlp_complete(resid = resid))

        self.tabWidget_helper.setCurrentWidget(self.tab_info)

        # worker.signals.result.connect(self.print_output)
        #worker.signals.finished.connect(self.thread_complete)
       # worker.signals.progress.connect(self.progress_fn)
        self.threadpool.start(worker_mlp_wk)       



    def mlp_search(self, resid = False, resid_GP = False):
        global fit

        #omega = 1/ np.logspace(np.log10(self.mlp_min_period.value()), np.log10(self.mlp_max_period.value()), num=int(self.mlp_n_omega.value()))
        ind_norm = self.gls_norm_combo.currentIndex()

        if self.mlp_cross_hair.isChecked():
            gls_like = True
        else:
            gls_like = False            
 
        mlp_N_threads = int(self.mlp_N_threads.value())

        if len(fit.fit_results.rv_model.jd) > 5:  
            
            rv_files_for_mlp = []
            for i in range(fit.ndset):
                
                if resid == True and resid_GP == False:
                    typ = (fit.fit_results.rv_model.jd[fit.fit_results.idset==i],
                       fit.fit_results.rv_model.o_c[fit.fit_results.idset==i], 
                       fit.fit_results.rv_model.rv_err[fit.fit_results.idset==i])
                       
                elif resid_GP == True:
                    
                    if fit.doGP == False:
                        print("No GP model used, thus exiting")
                        return
                
                    typ = (fit.fit_results.rv_model.jd[fit.fit_results.idset==i],
                       fit.fit_results.rv_model.o_c[fit.fit_results.idset==i] - fit.gp_model_data[0][fit.fit_results.idset==i]  , 
                       fit.fit_results.rv_model.rv_err[fit.fit_results.idset==i])                                  
 
                else:
                    typ = (fit.fit_results.rv_model.jd[fit.fit_results.idset==i],
                       fit.fit_results.rv_model.rvs[fit.fit_results.idset==i], 
                       fit.fit_results.rv_model.rv_err[fit.fit_results.idset==i])
                    
                rv_files_for_mlp.append(typ)
 
            RV_per = mlp.Gls(rv_files_for_mlp, fast=True,  verbose=False, nojit=True,ls=gls_like,ncpus= mlp_N_threads,
            ofac=self.mlp_ofac.value(), fbeg=1/self.mlp_max_period.value(), fend=1/self.mlp_min_period.value(), norm='dlnL')

        else:
            return
 
        if resid == True:
            fit.mlp = RV_per  # TB Fixed  
        else:
            fit.mlp = RV_per  # TB Fixed  


    def update_RV_MLP_plots(self):
        global fit, p_mlp 
        
        # a bug fix... must be done smarter.
        if not  hasattr(fit.mlp, 'freq'):
            return
 
        p_mlp.plot(clear=True,)

        self.colors_gls.setStyleSheet("color: %s;"%dill.copy(fit.gls_colors[0]))
        self.colors_alias_mlp.setStyleSheet("color: %s;"%colors_MLP_alias[0])

        power_levels_gui = np.array([self.mlp_fap1.value(),self.mlp_fap2.value(),self.mlp_fap3.value()])

        Delta_f = max(fit.mlp.freq) - min(fit.mlp.freq)
        delta_f = 1.0/(max(fit.fit_results.rv_model.jd) - min(fit.fit_results.rv_model.jd))
        M = Delta_f/delta_f
       # print(M, power_levels)
       # return
        power_levels = np.log(1.0/(power_levels_gui/M))
        #power_levels = np.array([5,10,15])

        gls_model_width = float(self.gls_model_width.value())


        if len(fit.fit_results.rv_model.jd) > 5:

            ######################## GLS ##############################
            if self.radioButton_RV_MLP_period.isChecked():
                p_mlp.setLogMode(True,False)        
                p_mlp.plot(1/fit.mlp.freq, fit.mlp.power,pen={'color': dill.copy(fit.gls_colors[0]), 'width': self.gls_model_width.value()},symbol=None ) 
                p_mlp.setLabel('bottom', 'period [d]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})    
            else:
                p_mlp.setLogMode(False,False)        
                p_mlp.plot(fit.mlp.freq, fit.mlp.power,pen={'color': dill.copy(fit.gls_colors[0]), 'width': self.gls_model_width.value()},symbol=None )                
                p_mlp.setLabel('bottom', 'frequency [1/d]', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 

            if fit.mlp.norm == 'dlnL':
                [p_mlp.addLine(x=None, y=fap, pen=pg.mkPen('k', width=0.8, style=QtCore.Qt.PenStyle.DotLine)) for fap in np.array(power_levels)]
 
#            text_peaks, pos_peaks = self.identify_power_peaks(1/fit.mlp.freq, fit.mlp.power, power_level = power_levels, sig_level = fit.mlp.powerLevel(np.array(power_levels)) )   
            text_peaks, pos_peaks = self.identify_power_peaks(1/fit.mlp.freq, fit.mlp.power, power_level = power_levels_gui, sig_level =np.array(power_levels) )   

            self.label_peaks(p_mlp, pos_peaks, GLS = True, MLP = True)

            self.mlp_print_info.clicked.connect(lambda: self.print_info_for_object(
            fit.mlp.info(stdout=False) + text_peaks))

        if self.mlp_cross_hair.isChecked():
            self.cross_hair(p_mlp,log=self.radioButton_RV_MLP_period.isChecked(), alias=[self.show_alias_MLP.isChecked(), self.alias_days_mlp.value(), colors_MLP_alias[0]])    


############ TLS ##############################      
       
    def worker_tls_complete(self, resid = False):
        global fit  
 
        #start_time = time.time()   
        
        if resid == False:
            self.update_tls_plots() 
        else:
            self.update_tls_o_c_plots() 
                 
        self.statusBar().showMessage('')   
 
        self.jupiter_push_vars()   
        self.calc_TLS.setEnabled(True)         
        self.calc_TLS_o_c.setEnabled(True)  
       # print("--- %s seconds ---" % (time.time() - start_time))     
 
    def worker_tls(self, resid = False):
        global fit  
        
        if tls_not_found==True:
            print("TLS Not found, try to install with 'pip install transitleastsquares'") 
            return

        self.calc_TLS.setEnabled(False)         
        self.calc_TLS_o_c.setEnabled(False)  

        z=0
        for i in range(20):
            if len(fit.tra_data_sets[i]) != 0:
                z=z+1

        if z <= 0:
            choice = QtWidgets.QMessageBox.information(self, 'Warning!',
            "Not possible to look for planets if there are no transit data loaded. Please add your transit data first. Okay?", QtWidgets.QMessageBox.StandardButton.Ok)
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

        N_transit_files = len([x for x in range(20) if len(fit.tra_data_sets[x]) != 0])


        if resid == True:
            lc_data = np.concatenate([fit.transit_results[1][4][x] +1.0 for x in range(20) if len(fit.tra_data_sets[x]) != 0])
        else:
            lc_data = np.concatenate([fit.transit_results[1][1][x] for x in range(20) if len(fit.tra_data_sets[x]) != 0]) 

        lc_time = np.concatenate([fit.transit_results[1][0][x] for x in range(20) if len(fit.tra_data_sets[x]) != 0])

        if self.tls_period_min_use.isChecked():
            tls_period_min = self.tls_min_period.value()
        else:
            tls_period_min = 0

        if self.tls_period_max_use.isChecked():
            tls_period_max = self.tls_max_period.value()
        else:
            tls_period_max = np.inf
 
        
        tls_model = transitleastsquares_ES(lc_time, lc_data)
        tls_results = tls_model.power(period_min = tls_period_min, period_max=tls_period_max, 
                                      oversampling_factor=int(self.tls_ofac.value()),
                                      duration_grid_step=self.tls_grid_step.value())
    
        if resid == True:
            fit.tls_o_c = tls_results  # TB Fixed with an rvmod object (i.e. fit.tls_obj)
        else:
            fit.tls = tls_results  # TB Fixed with an rvmod object (i.e. fit.tls_obj)


    def update_tls_plots(self): 
        global fit, p9, colors

        if len(fit.tls) == 0:
            return

        N_transit_files = len([x for x in range(20) if len(fit.tra_data_sets[x]) != 0])
        SDE_levels = np.array([self.tls_fap_1.value(),self.tls_fap_2.value(),self.tls_fap_3.value()])

        
        p9.plot(clear=True,) 
        
        if self.tls_cross_hair.isChecked():
            self.cross_hair(p9,log=False)
            
        if N_transit_files != 0:

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
            [p9.addLine(x=None, y=fap, pen=pg.mkPen('k', width=0.8, style=QtCore.Qt.PenStyle.DotLine)) for ii,fap in enumerate(SDE_levels)]

            text_peaks, pos_peaks = self.identify_power_peaks(fit.tls.periods,fit.tls.power,  sig_level = SDE_levels   )

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

        N_transit_files = len([x for x in range(20) if len(fit.tra_data_sets[x]) != 0])
        SDE_levels = np.array([self.tls_fap_1.value(),self.tls_fap_2.value(),self.tls_fap_3.value()])


        p10.plot(clear=True,) 

        if self.tls_o_c_cross_hair.isChecked():
            self.cross_hair(p10,log=False) 
            
        if N_transit_files != 0:


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
            [p10.addLine(x=None, y=fap, pen=pg.mkPen('k', width=0.8, style=QtCore.Qt.PenStyle.DotLine)) for ii,fap in enumerate(SDE_levels)]
   
            text_peaks, pos_peaks = self.identify_power_peaks(fit.tls_o_c.periods,fit.tls_o_c.power,  sig_level = SDE_levels )
 
            self.label_peaks(p10, pos_peaks, GLS = False)
            
            self.tls_o_c_print_info.clicked.connect(lambda: self.print_info_for_object(text + text_peaks))
            
            return

        else:
            text_err = pg.TextItem('Nothing to plot',color=(0,0,0))#, anchor=(0,0), border='w',color) #, fill=(0, 0, 255, 100))
            p10.addItem(text_err, ignoreBounds=True)
            self.tls_o_c_print_info.clicked.connect(lambda: self.print_info_for_object(""))
            return


############ transit fitting ##############################      

    def worker_transit_fitting_complete(self):
        global fit  


        self.update_labels()
        self.update_gui_params()
        self.update_errors() 
        self.update_a_mass()
        self.update_transit_combo_phase_pl()


        fit=rv.get_xyz(fit)
                         
        self.statusBar().showMessage('')  
        
       # self.button_fit.setEnabled(True)
        
        if fit.bound_error == True:
            self.get_error_msg(fit.bound_error_msg)
            self.mute_buttons(trigger=True)
            return

        self.update_transit_plots()
        
        if fit.type_fit["RV"] == True:

            for i in range(9):

                if not bool(fit.use_planet[i]):
                    continue
            #for i in range(fit.npl):
                try:
                    rv.phase_RV_planet_signal(fit,i)
                except:
                    print("Old session (pre Ver. 0.86)? You are using the old phased plot version!")
                    rv.phase_RV_planet_signal_old(fit,i)  
           # self.run_gls()
          #  self.run_gls_o_c()
            self.update_plots()  
        self.jupiter_push_vars() 

        self.save_last_session("autosave/auto_save.ses")
        self.mute_buttons(trigger=True)


    def check_model_fact(self):
        global fit  
        fit.tra_model_fact = int(self.tra_model_ndata_fact.value())
        
        

    def worker_transit_fitting(self, ff=1, auto_fit = False ):
        global fit  
        
        #self.button_fit.setEnabled(False)  
        self.mute_buttons(trigger=False)
        
        self.update_params() 
        self.update_use()   
        
        # check if transit data is present
        if self.check_data_and_fit_type(tf="fit") == 0:
            self.mute_buttons(trigger=True)
            return

        if fit.type_fit["RV"] == True :        
            self.statusBar().showMessage('Minimizing Transit + RV parameters.... SciPy in action, please be patient.  ')       
        else:
            self.statusBar().showMessage('Minimizing Transit parameters.... SciPy in action, please be patient. ')       
           
        self.check_ttv_params()
        self.check_ast_params()
        self.set_tra_ld()
        self.check_bounds()
        self.check_priors_nr()
        self.check_priors_jeff()
        self.check_model_fact()

        self.check_scipy_min()
        fit.model_npoints = self.points_to_draw_model.value()
        fit.model_max = self.model_max_range.value()
        fit.model_min = self.model_min_range.value()
          
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


        rv.run_SciPyOp(fit)




    def update_transit_combo_phase_pl(self):
        global fit
        self.comboBox_phase_pl_tran.clear()
        
        #for i in range(fit.npl):

        for i in range(9):

            if not bool(fit.use_planet[i]):
                continue
            self.comboBox_phase_pl_tran.addItem('pl. %s'%str(i+1),i+1) 


        #self.comboBox_phase_pl_tran.setCurrentIndex(0)


 
############ TTV fitting ##############################      

    def worker_ttv_fitting_complete(self):
        global fit  

        self.update_labels()
        self.update_gui_params()
        self.update_errors() 
        self.update_a_mass()
        
        #print(fit.fit_results.mass)

        fit=rv.get_xyz(fit)
                         
        self.statusBar().showMessage('')  
        
        #self.button_fit.setEnabled(True)         
        
        if fit.bound_error == True:
            self.get_error_msg(fit.bound_error_msg)
            self.mute_buttons(trigger=True)
            return

        self.update_transit_plots()
        if fit.type_fit["RV"] == True:


            for i in range(9):

                if not bool(fit.use_planet[i]):
                    continue
           # for i in range(fit.npl):
                try:
                    rv.phase_RV_planet_signal(fit,i)
                except:
                    print("Old session (pre Ver. 0.86)? You are using the old phased plot version!")
                    rv.phase_RV_planet_signal_old(fit,i)  
           # self.run_gls()
           # self.run_gls_o_c()
        self.update_plots()
        self.jupiter_push_vars()
        
        self.save_last_session("autosave/auto_save.ses")
        self.mute_buttons(trigger=True)


    def worker_ttv_fitting(self, ff=1, auto_fit = False ):
        global fit  

        #self.button_fit.setEnabled(False)
        self.update_params() 
        self.update_use()   
        self.mute_buttons(trigger=False)

        if self.check_data_and_fit_type(tf="fit") == 0:
            self.mute_buttons(trigger=True)
            return

        if fit.type_fit["RV"] == True :
            self.statusBar().showMessage('Minimizing TTV + RV parameters.... SciPy in action, please be patient.  ')
        else:
            self.statusBar().showMessage('Minimizing TTV parameters.... SciPy in action, please be patient. ')

        self.check_ttv_params()
        self.check_ast_params()
        self.set_tra_ld()
        self.check_bounds()
        self.check_priors_nr()   
        self.check_priors_jeff()   
        self.check_scipy_min()

        worker_ttv = Worker(lambda:  self.ttv_fit(ff=ff ) )# Any other args, kwargs are passed to the run  
 
        worker_ttv.signals.finished.connect(self.worker_ttv_fitting_complete)
        
        # worker.signals.result.connect(self.print_output)
        #worker.signals.finished.connect(self.thread_complete)
       # worker.signals.progress.connect(self.progress_fn)
        self.threadpool.start(worker_ttv)


    def ttv_fit(self, ff=0 ):
        global fit

        if ff ==0:
            fit.init_fit = True
        else:
            fit.init_fit = False

        rv.run_SciPyOp(fit)


    def check_ttv_params(self):
        global fit

        if fit.type_fit["TTV"] == True:
            fit.epoch_ttv = self.Epoch_ttv.value()
            fit.ttv_dt = self.time_step_model_ttv.value() 
            fit.epoch_ttv_end = self.Epoch_ttv_end.value()
        else:
            fit.epoch_ttv = self.Epoch.value()
            fit.ttv_dt = self.time_step_model.value()
            fit.epoch_ttv_end = self.Epoch.value()+1000          
        
        fit.ttv_times = [fit.epoch_ttv,fit.ttv_dt,fit.epoch_ttv_end]


    def worker_Nbody_complete(self):
        global fit, colors, p13, p14, p15  
          


        self.plot_evol_all()
#        self.update_orb_plot()

        self.plot_tabs.setCurrentWidget(self.tab_Orbital_evol)
            
             
        self.button_orb_evol.setEnabled(True)       
        self.statusBar().showMessage('')      

        self.save_last_session("autosave/auto_save.ses")

 

############ Ast fitting ##############################      

    def worker_ast_fitting_complete(self):
        global fit  

        self.update_labels()
        self.update_gui_params()
        self.update_errors() 
        self.update_a_mass()
        
        #print(fit.fit_results.mass)

        fit=rv.get_xyz(fit)
                         
        self.statusBar().showMessage('')  
        
        #self.button_fit.setEnabled(True)         
        
        if fit.bound_error == True:
            self.get_error_msg(fit.bound_error_msg)
            self.mute_buttons(trigger=True)
            return

        self.update_transit_plots()
        if fit.type_fit["RV"] == True:
           # for i in range(fit.npl):

            for i in range(9):

                if not bool(fit.use_planet[i]):
                    continue
                try:
                    rv.phase_RV_planet_signal(fit,i)
                except:
                    print("Old session (pre Ver. 0.86)? You are using the old phased plot version!")
                    rv.phase_RV_planet_signal_old(fit,i)  
           # self.run_gls()
           # self.run_gls_o_c()
        self.update_plots()
        self.jupiter_push_vars()
        
        self.save_last_session("autosave/auto_save.ses")
        self.mute_buttons(trigger=True)
 
    def worker_ast_fitting(self, ff=1, auto_fit = False ):
        global fit  

        #self.button_fit.setEnabled(False)
        self.update_params() 
        self.update_use()   
        self.mute_buttons(trigger=False)

        if self.check_data_and_fit_type(tf="fit") == 0:
            self.mute_buttons(trigger=True)
            return

        if fit.type_fit["RV"] == True :
            self.statusBar().showMessage('Minimizing Astrometry + RV parameters.... SciPy in action, please be patient.  ')
        else:
            self.statusBar().showMessage('Minimizing Astrometry parameters.... SciPy in action, please be patient. ')

        self.check_ast_params()
        self.set_tra_ld()
        self.check_bounds()
        self.check_priors_nr()   
        self.check_priors_jeff()   
        self.check_scipy_min()

        worker_ast = Worker(lambda:  self.ast_fit(ff=ff ) )# Any other args, kwargs are passed to the run  
 
        worker_ast.signals.finished.connect(self.worker_ast_fitting_complete)
        
        # worker.signals.result.connect(self.print_output)
        #worker.signals.finished.connect(self.thread_complete)
       # worker.signals.progress.connect(self.progress_fn)
        self.threadpool.start(worker_ast)


    def ast_fit(self, ff=0 ):
        global fit

        if ff ==0:
            fit.init_fit = True
        else:
            fit.init_fit = False

        rv.run_SciPyOp(fit)


    def check_ast_params(self):
        global fit

        if fit.type_fit["AST"] == True:
            fit.epoch_ast = self.Epoch_ast.value()
            #fit.ast_dt = self.time_step_model_ast.value() 
            fit.epoch_ast_end = self.Epoch_ast_end.value()
        else:
            fit.epoch_ast = self.Epoch.value()
            fit.ast_dt = self.time_step_model.value()
            fit.epoch_ast_end = self.Epoch.value()+1000          
        
        fit.ast_times = [fit.epoch_ast,fit.ast_dt,fit.epoch_ast_end]
        #print(fit.ast_times)




    def worker_Nbody(self):
        global fit  

        self.button_orb_evol.setEnabled(False)         
        npl_to_fit = np.atleast_1d(fit.fit_results.mass)
 
        # check if any fits where performed, and tus planets present
        if len(np.atleast_1d(npl_to_fit)) == 1 and npl_to_fit[0] <= 0:
            choice = QtWidgets.QMessageBox.information(self, 'Warning!',
            "Not possible to integrate a fit that does not exist. First perform an orbital fitting and then test the orbital stability. Okay?", QtWidgets.QMessageBox.StandardButton.Ok)      
            self.button_orb_evol.setEnabled(True)         
            return        

        if fit.npl < 2:
            choice = QtWidgets.QMessageBox.information(self, 'Warning!'," With less than two planets this makes no sense. Okay?",
                                            QtWidgets.QMessageBox.StandardButton.Ok) 
            self.button_orb_evol.setEnabled(True)                    
            return
        
        ######## this is a fix in case one adds another planet without initialize it 
        #if fit.npl != len(fit.fit_results.mass):  
        #    fit.model_saved = False 
       #     self.init_fit()

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

        #print(fit.npl_arb)

        if fit.npl_arb < 2:
            choice = QtWidgets.QMessageBox.information(self, 'Warning!'," With less than two planets this makes no sense. Okay?",
                                            QtWidgets.QMessageBox.StandardButton.Ok) 
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
             
        start_time = time.time()        
       # fit.run_stability_last_fit_params(timemax=self.max_time_of_evol.value(), timestep=self.time_step_of_evol.value(), integrator=integrator)      
       
        if arbitary == True:
            fit = rv.run_stability_arb(fit, timemax=self.max_time_of_evol.value(), timestep=self.time_step_of_evol.value(), integrator=integrator)      
        else:         
            fit = rv.run_stability(fit, timemax=self.max_time_of_evol.value(), timestep=self.time_step_of_evol.value(), integrator=integrator)      

        self.jupiter_push_vars()         
        
        print("--- %s seconds ---" % (time.time() - start_time))          



############################# Fortran fitting ###############################        
        
    def worker_RV_fitting_complete(self):
        global fit  
        
       # start_time =time.time()   
        fit=rv.get_xyz(fit)
        fit.init_fit= False 
        
        self.update_labels()
        self.update_gui_params()
         #else:
        self.update_errors()      
            
            
            
            
        self.update_a_mass()
        self.update_plots()
        self.jupiter_push_vars() 
        #print("--- %s seconds ---" % (time.time() - start_time))     
        self.statusBar().showMessage('')
        #self.console_widget.print_text(str(fit.print_info(short_errors=False))) 
        #self.button_fit.setEnabled(True)  
        if fit.bound_error == True:
            self.get_error_msg(fit.bound_error_msg)
            self.mute_buttons(trigger=True)
            return
        
        self.run_gls()
        self.run_gls_o_c()

        self.save_last_session("autosave/auto_save.ses")
        self.mute_buttons(trigger=True)
#        print("--- %s seconds ---" % (time.time() - start_time))     

    def worker_RV_fitting(self, ff=20, m_ln=True, auto_fit = False, init = False):
        global fit  
        
        #self.button_fit.setEnabled(False)

        self.mute_buttons(trigger=False)
        # check if RV data is present
        if self.check_data_and_fit_type(tf="fit") == 0:
            self.mute_buttons(trigger=True)
            return

        self.check_model_params()

        self.check_bounds()
        self.check_priors_nr()   
        self.check_priors_jeff()   
        self.mute_buttons(trigger=False)
     
        fit.model_npoints = self.points_to_draw_model.value()
        fit.model_max = self.model_max_range.value()
        fit.model_min = self.model_min_range.value()
        
        self.check_settings()

        
        #self.tabWidget_helper.setCurrentWidget(self.tab_info)

        if init == True:
            fit.init_fit= True
            ff = 0
            doGP=False
        else:
            doGP=self.do_RV_GP.isChecked()
            fit.init_fit= False  
       
        if self.radioButton_fortran77.isChecked() and not self.do_RV_GP.isChecked():
            self.statusBar().showMessage('Minimizing parameters....')    
             # Pass the function to execute
            worker2 = Worker(lambda:  self.optimize_fit(ff=ff,  minimize_fortran=True, m_ln=m_ln, auto_fit = auto_fit)) # Any other args, kwargs are passed to the run  
 
            #worker2 = Worker(lambda:  fit.fitting(outputfiles=[1,1,1], doGP=False,   minimize_fortran=True,  fortran_kill=30, timeout_sec=self.master_timeout.value(),minimize_loglik=True,amoeba_starts=ff, print_stat=False, eps=self.dyn_model_accuracy.value(), dt=self.time_step_model.value(), npoints=self.points_to_draw_model.value(), model_max= self.model_max_range.value(), model_min= self.model_min_range.value()))
 
 
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

        if self.radioButton_scipy.isChecked():
            for i in range(20):     
                self.use_data_jitter_gui[i].setEnabled(True)
        else:
            if self.amoeba_radio_button.isChecked():
                for i in range(20):     
                    self.use_data_jitter_gui[i].setEnabled(True)
            else:
                for i in range(20):     
                    self.use_data_jitter_gui[i].setEnabled(False)
                
  
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
 
        old_err4 = dill.copy(fit.param_errors.GP_params_errors)
 
    
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
 

            fit.fitting(fileinput=self.fortran_debug.isChecked(),outputfiles=[1,1,1], doGP=doGP,  kernel_id=gp_kernel_id,  minimize_fortran=minimize_fortran,  fortran_kill=f_kill, timeout_sec=self.master_timeout.value(),minimize_loglik=True,amoeba_starts=ff, print_stat=False, eps=self.dyn_model_accuracy.value(), dt=self.time_step_model.value(), npoints=self.points_to_draw_model.value(), model_max= self.model_max_range.value(), model_min= self.model_min_range.value())

        elif m_ln == True and doGP == True:       
            fit.fitting(fileinput=self.fortran_debug.isChecked(),outputfiles=[1,1,1], doGP=doGP,  kernel_id=gp_kernel_id,  minimize_fortran=minimize_fortran,  fortran_kill=f_kill, timeout_sec=self.master_timeout.value(),minimize_loglik=True,amoeba_starts=ff,  print_stat=False, eps=self.dyn_model_accuracy.value(), dt=self.time_step_model.value(), npoints=self.points_to_draw_model.value(), model_max= self.model_max_range.value(), model_min= self.model_min_range.value())

        elif m_ln == False and  minimize_fortran==False:
            fit.fitting(fileinput=self.fortran_debug.isChecked(),outputfiles=[1,1,1], doGP=doGP,  kernel_id=gp_kernel_id,  minimize_fortran=minimize_fortran, fortran_kill=f_kill, timeout_sec=self.master_timeout.value(),minimize_loglik=True,amoeba_starts=0, print_stat=False,eps=self.dyn_model_accuracy.value(), dt=self.time_step_model.value(), npoints=self.points_to_draw_model.value(), model_max= self.model_max_range.value(), model_min= self.model_min_range.value())

        else:      
            fit.fitting(fileinput=self.fortran_debug.isChecked(),outputfiles=[1,1,1], doGP=doGP,  kernel_id=gp_kernel_id,  minimize_fortran=minimize_fortran, fortran_kill=f_kill, timeout_sec=self.master_timeout.value(),minimize_loglik=m_ln,amoeba_starts=ff, print_stat=False,eps=self.dyn_model_accuracy.value(), dt=self.time_step_model.value(), npoints=self.points_to_draw_model.value(), model_max= self.model_max_range.value(), model_min= self.model_min_range.value())

        if fit.doGP == True:
            rv.get_RV_gps_model(fit,get_lnl=fit.get_GP_lnl)

        #for i in range(fit.npl):

        self.update_gui_params()


        #sorted_periods = sorted([k for k, v in fit.P.items() if fit.use_planet[k]], key=fit.P.get) 
        for i in range(9):

            if not bool(fit.use_planet[i]):
                continue
            try:
                #print(fit.P[i], i)
                rv.phase_RV_planet_signal(fit,i)
            except:
                print("Old session (pre Ver. 0.86)? You are using the old phased plot version!")
                rv.phase_RV_planet_signal_old(fit,i)              

        if self.reset_errors_at_init.isChecked() == False and fit.init_fit == True:
            print("Warning: 'Reset the errors to 0 when Initialize' is set to 'False', thus errors are not overwritten!")
        #    fit.param_errors.planet_params_errors = dill.copy(old_err1)
        #    fit.param_errors.offset_errors = dill.copy(old_err2)
        #    fit.param_errors.jitter_errors = dill.copy(old_err3)
            fit.param_errors.GP_params_errors = dill.copy(old_err4)
            
          
            
            self.update_params()

        if auto_fit:
            self.update_labels()
            self.update_gui_params()
            self.update_errors() 
            self.update_a_mass() 

            #self.run_gls()
           # self.run_gls_o_c()
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
            np.sqrt(stat.chi2.ppf(0.95,1)),2.0,
            np.sqrt(stat.chi2.ppf(0.99,1)),3.0,
            np.sqrt(stat.chi2.ppf(0.999,1)),4.0,5.0  ]
        
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




    def print_more_stat(self):
        global fit 
        
        #self.dialog.statusBar().showMessage('Ready')
        #self.dialog_more_info.setFixedSize(500,450)
        self.dialog_more_info.setGeometry(1,1, 490, 385)
        self.dialog_more_info.setWindowTitle('Fit stat. info')  
        #self.dialog.setGeometry(300, 300, 800, 800)
        #self.dialog_credits.acceptRichText(True)
        text = ''
        self.dialog_more_info.text.setText(text)         
        ################## text generator #################
        #text_info = """ (Work in progress) """
       # self.dialog_more_info.text.setText(text_info)     
         
        text_info ="""
Fit quality
----------------------------------------------  
max lnL = %.5f
chi^2   = %.5f
red. chi^2   = %.5f
----------------------------------------------"""%(fit.loglik,fit.fit_results.chi2,fit.fit_results.reduced_chi2)   

        self.dialog_more_info.text.append(text_info)

        if fit.type_fit["RV"] == True:

            text_info = """
            
%s rms/wrms 
----------------------------------------------"""%(self.data_str)
            self.dialog_more_info.text.append(text_info)
            
            if fit.ndset != 0:
                #for i in range(max(fit.fit_results.idset)+1):

                for i in range(len(fit.rv_data_sets)):
                    if len(fit.rv_data_sets[i]) == 0:
                        continue         

                    text_info = """ """   
                    self.dialog_more_info.text.append(text_info)    
                    rms = np.sqrt(np.average(fit.fit_results.o_c[fit.fit_results.idset==i]**2))
                    text_wrm = "%s    rms = %.5f m/s"%(fit.rv_data_sets[i][4],rms)       
                    self.dialog_more_info.text.append(text_wrm)
                    wrms = np.sqrt(np.average(fit.fit_results.o_c[fit.fit_results.idset==i]**2, weights=1/fit.fit_results.rv_err[fit.fit_results.idset==i]))
                    text_wrm = "%s wrms = %.5f m/s"%(fit.rv_data_sets[i][4],wrms)       
                    self.dialog_more_info.text.append(text_wrm)        

        if fit.type_fit["Transit"] == True:
 
            
            text_info = """

Transit data rms/wrms 
----------------------------------------------"""   
            self.dialog_more_info.text.append(text_info)

            for i in range(20):
                if len(fit.tra_data_sets[i]) == 0:
                    continue
                else:
                    
                    text_info = """ """   
                    self.dialog_more_info.text.append(text_info)                
                    rms = np.sqrt(np.average(fit.tra_data_sets[i][4]**2))
                    text_wrm = "%s    rms = %.5f rel. flux"%(fit.tra_data_sets[i][-1],rms)       
                    self.dialog_more_info.text.append(text_wrm)
           
                    wrms = np.sqrt(np.average(fit.tra_data_sets[i][4]**2,
                                              weights=1/fit.tra_data_sets[i][2]))
                    text_wrm = "%s wrms = %.5f rel. flux"%(fit.tra_data_sets[i][-1],wrms)       
                    self.dialog_more_info.text.append(text_wrm)        

 
        self.dialog_more_info.text.setReadOnly(True)
   
        self.dialog_more_info.show()



    def print_GP_info(self):
        self.dialog_GP_help.setFixedSize(600, 600)
        self.dialog_GP_help.setWindowTitle('GP modeling help')

        text = ''
        self.dialog_GP_help.text.setText(text) 


        text = """
<br>
<br>
The GP parameters are explained in the <a href='https://ui.adsabs.harvard.edu/abs/2017AJ....154..220F/abstract'> 'Celerite' paper</a>
<br>  
<br>
<br> If you made the use of these GP kernels for your paper, please also cite: <a href='https://ui.adsabs.harvard.edu/abs/2017AJ....154..220F/abstract'> Foreman-Mackey et al. (2017)</a>
"""

        self.dialog_GP_help.text.append(text)
        self.dialog_GP_help.text.setReadOnly(True)
        self.dialog_GP_help.show()



    def print_mcmc_info(self):
        self.dialog_mcmc_help.setFixedSize(600, 600)
        self.dialog_mcmc_help.setWindowTitle('MCMC help')

        text = ''
        self.dialog_mcmc_help.text.setText(text) 


        text = """
<br>
<br>
The MCMC is done via the 'emcee' packade. The parameters and options are explained in the 'emcee: The MCMC Hammer' <a href='https://emcee.readthedocs.io/en/stable/'> user guide</a>.
<br>  
<br>
<br> If you made the use of the Exo-Striker's MCMC for your paper, please cite: <a href='https://ui.adsabs.harvard.edu/abs/2013PASP..125..306F/abstract'> Foreman-Mackey et al. (2013)</a>, and references therein.
"""

        self.dialog_mcmc_help.text.append(text)
        self.dialog_mcmc_help.text.setReadOnly(True)
        self.dialog_mcmc_help.show()

    def print_ns_info(self):
        self.dialog_ns_help.setFixedSize(600, 600)
        self.dialog_ns_help.setWindowTitle('Nested Sampling help')

        text = ''
        self.dialog_ns_help.text.setText(text) 


        text = """
<br>
<br>
The Nested Sampling is done via the 'dynesty' packade. The parameters and options are explained in the 'dynesty' <a href='https://dynesty.readthedocs.io/en/latest/'> user guide</a>.
<br>  
<br>
<br> If you made the use of the Exo-Striker's Nested Sampling for your paper, please also cite: <a href='https://ui.adsabs.harvard.edu/abs/2020MNRAS.493.3132S/abstract'> Speagle (2019)</a>, and references therein.
"""

        self.dialog_ns_help.text.append(text)
        self.dialog_ns_help.text.setReadOnly(True)
        self.dialog_ns_help.show()


    def print_ast_info(self, image=False):
        #self.dialog.statusBar().showMessage('Ready')
        self.dialog_ast_help.setFixedSize(800, 800)
        self.dialog_ast_help.setWindowTitle('Ast. modeling help')  
        #self.dialog.setGeometry(300, 300, 800, 800)
        #self.dialog_credits.acceptRichText(True)
        
        text = ''
        self.dialog_ast_help.text.setText(text) 
        
        text = """
This module is experimental and under construction. Atrometry and Atrometry+RV modeling was 
successfully tested on several exoplanet systems, but you must consider the notes below:
 
* (TBD)

If something is still unclear, or you are experiencing problems you are welcome to open an "issue" 
in https://github.com/3fon3fonov/exostriker

"""
        self.dialog_ast_help.text.append(text)
        self.dialog_ast_help.text.setReadOnly(True)
        #self.dialog.setWindowIcon (QtGui.QIcon('logo.png'))        
        self.dialog_ast_help.show()



    def print_ast_info_2(self, image=False):
        #self.dialog.statusBar().showMessage('Ready')
        self.dialog_ast_help.setFixedSize(800, 800)
        self.dialog_ast_help.setWindowTitle('Ast. modeling help')  
        #self.dialog.setGeometry(300, 300, 800, 800)
        #self.dialog_credits.acceptRichText(True)
        
        text = ''
        self.dialog_ast_help.text.setText(text) 
        
        text = """
This module is experimental and under construction. Atrometry and Atrometry+RV modeling was 
successfully tested on several exoplanet systems, but you must consider the notes below:
 
* (TBD, Hipparcos/Gaia astrometry)

If something is still unclear, or you are experiencing problems you are welcome to open an "issue" 
in https://github.com/3fon3fonov/exostriker

"""
        self.dialog_ast_help.text.append(text)
        self.dialog_ast_help.text.setReadOnly(True)
        #self.dialog.setWindowIcon (QtGui.QIcon('logo.png'))        
        self.dialog_ast_help.show()



    def print_TTV_info(self, image=False):
        #self.dialog.statusBar().showMessage('Ready')
        self.dialog_ttv_help.setFixedSize(800, 800)
        self.dialog_ttv_help.setWindowTitle('TTV modeling help')  
        #self.dialog.setGeometry(300, 300, 800, 800)
        #self.dialog_credits.acceptRichText(True)
        
        text = ''
        self.dialog_ttv_help.text.setText(text) 
        
        text = """
This module is experimental and under construction. TTVs and TTVs+RV modeling was 
successfully tested on several exoplanet systems, but you must consider the notes below:

* The TTV file format should be as follows (e.g.):


#   N transit             t0 [BJD]             sigma t0 [d]
          1                    2458000.5                0.022
          2                    2458020.5                0.023
          4                    2458060.5                0.021
.......................................................................
          8                    2458140.5                0.026

    Anything else than the format above will likely not work and is even possible to crash the GUI!


* The selected epoch MUST be always slightly before the time of the first observed transit,
Using the example above, the epoch should be earlier than 2458000.5, e.g., 2458000.0. 
Otherwise, the TTV model is likely to skip the first transit and start from the next!


* When RV+TTVs are modeled the epoch is ALWAYS chosen to be the epoch of the RV model.

    To change the RV epoch go to:

    Models param. --> Models --> RV Model

    Then, uncheck "first RV" and add whatever epoch you like, as long as it is slightly before the time 
    of the first transit in your TTV input file. 


* Make sure that the time baseline of "End of Model" - "Epoch" >  last t0 - first t0 in your TTV 
input file.


* Chose wisely the "time step in dynamical model" !

If something is still unclear, or you are experiencing problems you are welcome to open an "issue" 
in https://github.com/3fon3fonov/exostriker

"""
        self.dialog_ttv_help.text.append(text)
        self.dialog_ttv_help.text.setReadOnly(True)
        #self.dialog.setWindowIcon (QtGui.QIcon('logo.png'))        
        self.dialog_ttv_help.show()



    def print_info_credits(self, image=False):
 
        
        #self.dialog.statusBar().showMessage('Ready')
        self.dialog_credits.setFixedSize(900, 900)
        self.dialog_credits.setWindowTitle('Credits')  
        #self.dialog.setGeometry(300, 300, 800, 800)
        #self.dialog_credits.acceptRichText(True)
        
        if qso_mode:
            String_es = "QSO"
        else:
            String_es = "Exo"            

        text = ''
        self.dialog_credits.text.setText(text) 
        
        text = "You are using 'The %s-Striker' (ver. %s) \n developed by Trifon Trifonov"%(String_es,self.es_version)
        
        self.dialog_credits.text.append(text)

        text = "\n"*15 +"CREDITS:"+"\n"*2 + "This tool uses the publically \n available packages: \n" 
        self.dialog_credits.text.append(text)
        
        text = "* " + "<a href='https://github.com/pyqtgraph/pyqtgraph'>pyqtgraph</a>"
        self.dialog_credits.text.append(text)
        
        text = "* " + "<a href='https://github.com/mzechmeister/python'>GLS and MLP periogograms</a>"
        self.dialog_credits.text.append(text)
 
        text = "* " + "<a href='https://github.com/dfm/emcee'>emcee</a>" 
        self.dialog_credits.text.append(text) 

        text = "* " + "<a href='https://github.com/dfm/celerite'>celerite</a>" 
        self.dialog_credits.text.append(text)
                
        text = "* " + "<a href='https://dynesty.readthedocs.io/en/latest/'>dynesty</a>" 
        self.dialog_credits.text.append(text)
               
        if qso_mode == False:
            text = "* " + "<a href='https://github.com/lkreidberg/batman'>batman-package</a>" 
            self.dialog_credits.text.append(text)
        
            text = "* " + "<a href='https://github.com/hippke/tls'>transitleastsquares</a>" 
            self.dialog_credits.text.append(text)

            text = "* " + "<a href='https://github.com/hippke/wotan'>wotan</a>" 
            self.dialog_credits.text.append(text)        

            text = "* " + "<a href='https://github.com/mindriot101/ttvfast-python'>ttvfast-python</a>" 
            self.dialog_credits.text.append(text)

            text = "* " + "<a href='https://www.boulder.swri.edu/~hal/swift.html'>swift</a>" 
            self.dialog_credits.text.append(text)

        text = "* " + "<a href='https://github.com/mfitzp/15-minute-apps/tree/master/wordprocessor'>megasolid idiom</a>" 
        self.dialog_credits.text.append(text)  

        text = "(And many 'standard' Python libraries like \n PyQt6, matplotlib, numpy, scipy, pathos, \n jupyter, qtconsole) \n" 
        self.dialog_credits.text.append(text)   


        #self.dialog_credits.text.setText(text)
        #self.dialog_credits.text.insertHtml(text)
        
        
        text = "\n"*5 + """Note:
Please keep in mind that this software is developed 
mostly for my needs and for fun. I hope, however, 
that you may find it capable of solving your scientific 
problems, too.

Feedback and help in further development 
will be highly appreciated!
"""
        self.dialog_credits.text.append(text)
        self.dialog_credits.text.setReadOnly(True)

        if qso_mode:
            self.dialog_credits.setStyleSheet(" QTextEdit{border-image: url(./lib/UI/33_striker_qso.png) 0 0 0 0 stretch stretch;} ")
        else:
            self.dialog_credits.setStyleSheet(" QTextEdit{border-image: url(./lib/UI/33_striker.png) 0 0 0 0 stretch stretch;} ")

        #self.dialog.setWindowIcon (QtGui.QIcon('logo.png'))        
        self.dialog_credits.show()


    def find_planets(self):
        global fit 

        # check if RV data is present
        if fit.ndset <= 0:
             choice = QtWidgets.QMessageBox.information(self, 'Warning!',
             "Not possible to look for planets if there are no %s loaded. Please add your %s first. Okay?"%(self.data_str,self.data_str), QtWidgets.QMessageBox.StandardButton.Ok)      
             #self.mute_buttons(trigger=True)       
             return        

        self.button_auto_fit.setEnabled(False)   

        # the first one on the data GLS
        if fit.gls.power.max() <= fit.gls.powerLevel(self.auto_fit_FAP_level.value()):
             choice = QtWidgets.QMessageBox.information(self, 'Warning!',
             "No significant power on the GLS. Therefore no planets to fit OK?", QtWidgets.QMessageBox.StandardButton.Ok)      
             #self.mute_buttons(trigger=True)                                                          
             return
        
        else:
            if fit.npl !=0:
               # for j in range(fit.npl):
                for i in range(9):

                    if not bool(fit.use_planet[i]):
                        continue

                    fit.remove_planet(i)

            mean_anomaly_from_gls = np.degrees((((fit.epoch - float(fit.gls.hpstat["T0"]) )% (fit.gls.hpstat["P"]) )/ (fit.gls.hpstat["P"]) ) * 2*np.pi)
             
            fit.add_planet(fit.gls.hpstat["amp"],fit.gls.hpstat["P"],0.0,0.0,mean_anomaly_from_gls -90.0,90.0,0.0)

            if qso_mode==True:
                fit.P_use[0] = True
                fit.K_use[0] = True
                fit.e_use[0] = False
                fit.w_use[0] = False
                fit.M0_use[0] = True
                fit.i_use[0]  = False
                fit.Node_use[0] = False     

            else:
                fit.P_use[0] = True
                fit.K_use[0] = True
                fit.e_use[0] = True
                fit.w_use[0] = True
                fit.M0_use[0] = True
                fit.i_use[0]  = False
                fit.Node_use[0] = False   

            self.update_use_from_session()   
            self.update_use()                     
            #self.optimize_fit(20,m_ln=self.amoeba_radio_button.isChecked(),auto_fit = True)
            self.optimize_fit(20,m_ln=False,auto_fit = True) # first time -- LM
            
            
            #now inspect the residuals

            for i in range(1,int(self.auto_fit_N_planets.value())):
                
                if fit.gls_o_c.power.max() <= fit.gls_o_c.powerLevel(self.auto_fit_FAP_level.value()):
                    #for j in range(fit.npl):

                    for j in range(9):

                        if not bool(fit.use_planet[j]):
                            continue

                        fit.P_use[j] = True
                        fit.K_use[j] = True
                        fit.e_use[j] = True
                        fit.w_use[j] = True
                        fit.M0_use[j] = True
                        fit.i_use[j]  = False
                        fit.Node_use[j] = False 
                    self.update_use_from_session()   
                    self.update_use()

                   # fit.sort_by_period(reverse=False)
                    self.optimize_fit(20,m_ln=False,auto_fit = True) 
                    #self.button_auto_fit.setEnabled(True)     
                    return
                else:    
                    mean_anomaly_from_gls = np.degrees((((fit.epoch - float(fit.gls_o_c.hpstat["T0"]) )% (fit.gls_o_c.hpstat["P"]) )/ (fit.gls_o_c.hpstat["P"]) ) * 2*np.pi)
             
                    fit.add_planet(fit.gls_o_c.hpstat["amp"],fit.gls_o_c.hpstat["P"],0.0,0.0,mean_anomaly_from_gls -90.0,90.0,0.0)
 
                    if qso_mode==True:
                        fit.P_use[i] = True
                        fit.K_use[i] = True
                        fit.e_use[i] = False
                        fit.w_use[i] = False
                        fit.M0_use[i] = True
                        fit.i_use[i]  = False
                        fit.Node_use[i] = False      
                    else:
                        fit.P_use[i] = True
                        fit.K_use[i] = True
                        fit.e_use[i] = True
                        fit.w_use[i] = True
                        fit.M0_use[i] = True
                        fit.i_use[i]  = False
                        fit.Node_use[i] = False  
                    self.update_use_from_session()   
                    self.update_use()

                    fit.sort_by_period(reverse=False)
                    self.optimize_fit(20,m_ln=False,auto_fit = True)  

                #else:
                 #   continue

           # for j in range(fit.npl):

            for i in range(9):

                if not bool(fit.use_planet[i]):
                    continue

                if qso_mode==True:
                    fit.P_use[i] = True
                    fit.K_use[i] = True
                    fit.e_use[i] = False
                    fit.w_use[i] = False
                    fit.M0_use[i] = True
                    fit.i_use[i]  = False
                    fit.Node_use[i] = False  
                else:
                    fit.P_use[i] = True
                    fit.K_use[i] = True
                    fit.e_use[i] = True
                    fit.w_use[i] = True
                    fit.M0_use[i] = True
                    fit.i_use[i]  = False
                    fit.Node_use[i] = False  


            self.update_use_from_session()
            self.update_use()
            fit.sort_by_period(reverse=False)

            self.optimize_fit(20,m_ln=self.amoeba_radio_button.isChecked(),auto_fit = True)
 
        self.button_auto_fit.setEnabled(True)   


    def run_auto_fit(self):
        global fit 

 
        self.radioButton_Keplerian.setChecked(True) # this is to be fixed! Only with keplerian fitting th autofit works fine so far.
        #self.button_auto_fit.setEnabled(False)
        self.mute_buttons(trigger=False)    
        
        if fit.npl != 0:        
            choice = QtWidgets.QMessageBox.information(self, 'Warning!',
                                            "Planets already exist. Do you want to overwrite the analysis?",
                                            QtWidgets.QMessageBox.StandardButton.No | QtWidgets.QMessageBox.StandardButton.Yes)  

            if choice == QtWidgets.QMessageBox.StandardButton.No:
                #self.button_auto_fit.setEnabled(True)
                self.mute_buttons(trigger=True)

                return
            elif choice == QtWidgets.QMessageBox.StandardButton.Yes:
                #for j in range(fit.npl):
                #fit.npl = 0
                for j in range(9):

                    if not bool(fit.use_planet[j]):
                        continue

                    fit.remove_planet(j)
                
                self.find_planets()
        else:

           #time.sleep(0.5)
            self.find_planets()

        fit.type_fit["RV"] = True
        fit.type_fit["Transit"] = False
        fit.type_fit["TTV"] = False
        self.check_type_fit()
        self.mute_boxes()
        self.plot_tabs.setCurrentWidget(self.tab_timeseries_RV)
        self.mute_buttons(trigger=True)
        


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

        self.K7.minimize_signal.connect(lambda: fit.minimize_one_param_K(6)) #TBD!
        self.K7.minimize_signal.connect(self.init_fit) #TBD!       
        self.P7.minimize_signal.connect(lambda: fit.minimize_one_param_P(6)) #TBD!
        self.P7.minimize_signal.connect(self.init_fit) #TBD!     
        self.e7.minimize_signal.connect(lambda: fit.minimize_one_param_e(6)) #TBD!
        self.e7.minimize_signal.connect(self.init_fit) #TBD!  
        self.om7.minimize_signal.connect(lambda: fit.minimize_one_param_w(6)) #TBD!
        self.om7.minimize_signal.connect(self.init_fit) #TBD!  
        self.ma7.minimize_signal.connect(lambda: fit.minimize_one_param_M0(6)) #TBD!
        self.ma7.minimize_signal.connect(self.init_fit) #TBD!   

        self.K8.minimize_signal.connect(lambda: fit.minimize_one_param_K(7)) #TBD!
        self.K8.minimize_signal.connect(self.init_fit) #TBD!       
        self.P8.minimize_signal.connect(lambda: fit.minimize_one_param_P(7)) #TBD!
        self.P8.minimize_signal.connect(self.init_fit) #TBD!     
        self.e8.minimize_signal.connect(lambda: fit.minimize_one_param_e(7)) #TBD!
        self.e8.minimize_signal.connect(self.init_fit) #TBD!  
        self.om8.minimize_signal.connect(lambda: fit.minimize_one_param_w(7)) #TBD!
        self.om8.minimize_signal.connect(self.init_fit) #TBD!  
        self.ma8.minimize_signal.connect(lambda: fit.minimize_one_param_M0(7)) #TBD!
        self.ma8.minimize_signal.connect(self.init_fit) #TBD!   

        self.K9.minimize_signal.connect(lambda: fit.minimize_one_param_K(8)) #TBD!
        self.K9.minimize_signal.connect(self.init_fit) #TBD!       
        self.P9.minimize_signal.connect(lambda: fit.minimize_one_param_P(8)) #TBD!
        self.P9.minimize_signal.connect(self.init_fit) #TBD!     
        self.e9.minimize_signal.connect(lambda: fit.minimize_one_param_e(8)) #TBD!
        self.e9.minimize_signal.connect(self.init_fit) #TBD!  
        self.om9.minimize_signal.connect(lambda: fit.minimize_one_param_w(8)) #TBD!
        self.om9.minimize_signal.connect(self.init_fit) #TBD!  
        self.ma9.minimize_signal.connect(lambda: fit.minimize_one_param_M0(8)) #TBD!
        self.ma9.minimize_signal.connect(self.init_fit) #TBD!   


    def jupiter_push_vars(self):
        global fit
        self.console_widget.push_vars({'fit':fit})
        
        #self.console_widget.push_vars({'pg':pg})    

        #self.console_widget.clear()         
        #self.console_widget.print_text(str("Welcome!"+"\n")) 


########################## Sessions ##################################
 
    def getNewses(self):
        global fit, ses_list  
        
        text, okPressed = QtWidgets.QInputDialog.getText(self, "New session","Name session: (No special characters!)", QtWidgets.QLineEdit.Normal, "")
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
#            self.select_session(self.comboBox_select_ses.count()-1)
#            self.select_session(0)        
            #print(self.comboBox_select_ses.count() )
            #self.change_session_label()
           
    def rem_ses(self):
        global fit, ses_list  
        
        ind = self.comboBox_select_ses.currentIndex()
        
        choice = QtWidgets.QMessageBox.information(self, 'Warning!',
        "Do you really want to remove Session %s"%(ses_list[ind].name),
                                            QtWidgets.QMessageBox.StandardButton.No | QtWidgets.QMessageBox.StandardButton.Yes)  
         
        if choice == QtWidgets.QMessageBox.StandardButton.No:
            return
        elif choice == QtWidgets.QMessageBox.StandardButton.Yes: # and ind <=0:
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

        input_file = QtWidgets.QFileDialog.getOpenFileName(self, 'Open session', '', 'Data (*.ses)', options=QtWidgets.QFileDialog.Option.DontUseNativeDialog)

        if str(input_file[0]) != '':

            try:
                file_pi = open(input_file[0], 'rb')
                fit_new = dill.load(file_pi) #, encoding='latin1'
                file_pi.close()     
            except (UnicodeDecodeError, ImportError, KeyError) as e:
                py3_ses = rv.convert_Session_to_Py3(input_file[0])
                
                file_pi = open(py3_ses, 'rb')
                fit_new = dill.load(file_pi) #, encoding='latin1'
                file_pi.close()     


            #print(fit.bounds.GP_params_bounds, fit.bounds.planet_params_bounds, fit.bounds.jitter_bounds, fit.bounds.stellar_mass_bounds, fit.bounds.linear_trend_bounds)

            fit_new = rv.check_for_missing_instances(fit,fit_new)
            #self.check_for_missing_instances(fit_new)

            ses_list.append(fit_new)

            self.check_std_ast_solution()

            #self.check_settings()
            #rv.check_temp_RV_file(fit_new)

            self.session_list()
            
            self.select_session(-1)


    def save_session(self):
        global fit
        
        output_file = QtWidgets.QFileDialog.getSaveFileName(self, 'Save session', '%s.ses'%fit.name, 'Data (*.ses)', options=QtWidgets.QFileDialog.Option.DontUseNativeDialog)


        if str(output_file[0]) != '':
            file_pi = open(output_file[0], 'wb')
            dill.dump(fit, file_pi) #,protocol=2
            file_pi.close()

    def save_last_session(self,session_name):
        global fit


        fit.ses_notes = self.text_editor.editor.toPlainText()

        if str(session_name) != '':
            file_pi = open(session_name, 'wb')
            dill.dump(fit, file_pi) #,protocol=2
            file_pi.close()


    def open_sessions(self):
        global fit, ses_list

        input_file = QtWidgets.QFileDialog.getOpenFileName(self, 'Open session', '', 'Data (*.mses)', options=QtWidgets.QFileDialog.Option.DontUseNativeDialog)

        if str(input_file[0]) != '':

                
            file_pi = open(input_file[0], 'rb')
            fit2 = dill.load(file_pi)
            file_pi.close()

            for jj in fit2:
                jj = rv.check_for_missing_instances(fit,jj)

                #self.check_for_missing_instances(jj)


            choice = QtWidgets.QMessageBox.information(self, 'Warning!',
                                            "Do you want to overwrite the current sessions? If you choose 'No' will add the session, 'Cancel' will exit",
                                            QtWidgets.QMessageBox.StandardButton.Cancel | QtWidgets.QMessageBox.StandardButton.No | QtWidgets.QMessageBox.StandardButton.Yes)  

            if choice == QtWidgets.QMessageBox.StandardButton.No:
                ses_list = ses_list + fit2
            elif choice == QtWidgets.QMessageBox.StandardButton.Yes:
                ses_list = fit2
            elif choice == QtWidgets.QMessageBox.StandardButton.Cancel:
                return


            self.update_settings()  
            self.session_list()
            self.select_session(-1)

    def save_sessions(self):
        global fit, ses_list

        output_file = QtWidgets.QFileDialog.getSaveFileName(self, 'Save multi-session', 'my_sessions.mses', 'Data (*.mses)', options=QtWidgets.QFileDialog.Option.DontUseNativeDialog)

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
            

        self.check_type_fit()

        
        self.update_settings()
        self.update_bounds()

        self.set_type_fit_options()
        #self.mute_boxes()

        #self.init_fit()

        #fit.update_rv_params()

        if len(fit.tra_data_sets) == 10:
            fit = rv.fix_old_to_session_tra(fit)

        
        self.update_use_from_session()
        self.update_use()
        self.update_gui_params()
        self.update_errors() 
        self.update_a_mass()  
        self.update_plots() 
        self.plot_evol_all()
        self.update_labels()
        self.check_cornerplot_samples()          

        self.update_params()
        self.update_RV_file_buttons() 
        
        #if fit.type_fit["Transit"] == True:
        self.update_tra_file_buttons()
        self.update_ttv_file_buttons()
        self.update_act_file_buttons()
        self.update_ast_file_buttons()
        self.update_ast_file_buttons_2()        
        self.update_color_picker()
 
        self.set_gui_use_GP()
        self.check_use_tra_GP()
        


        #self.init_fit()
        #self.update_use_from_session()
        #self.update_use()
        #self.update_gui_params()
       # self.update_params()
       # self.update_RV_file_buttons()
       # self.update_tra_file_buttons()
        #self.update_act_file_buttons()

       # self.fit_dispatcher(init=True)
        self.init_plot_corr()
        self.update_plot_corr()

        self.update_GUI_tra_reg()
        self.update_GUI_mcmc_params()
        self.update_GUI_ns_params()
        self.update_GUI_St_params()
        self.check_cornerplot_samples()
        self.jupiter_push_vars()

        self.text_editor.editor.setText(fit.ses_notes)
        
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

        self.update_labels()
        self.update_gui_params()
        self.update_errors() 
        self.update_a_mass() 

        self.init_plot_corr()

        self.statusBar().showMessage('') 

        if fit.bound_error == True:
            self.get_error_msg(fit.bound_error_msg)
            self.mute_buttons(trigger=True)

            return

        if self.adopt_nest_means_as_par.isChecked() or self.adopt_nest_median_as_par.isChecked() or  self.adopt_nest_best_lnL_as_pars.isChecked() or self.adopt_nest_mode_as_par.isChecked():
            self.init_fit()
        else:
            self.jupiter_push_vars()

        self.save_last_session("autosave/auto_save.ses")
        self.check_cornerplot_samples()
        self.mute_buttons(trigger=True)



    def check_data_and_fit_type(self, tf = ""):
        global fit
        if fit.type_fit["RV"] == True and fit.ndset <= 0:
             choice = QtWidgets.QMessageBox.information(self, 'Warning!',
             "Not possible to run %s if there are no %s loaded. Please add your %s first. Okay?"%(tf, self.data_str,self.data_str), QtWidgets.QMessageBox.StandardButton.Ok)      
            # self.button_nest_samp.setEnabled(True)  
             self.statusBar().showMessage('') 
             return 0

        ntran_data = 0
        for i in range(0,len(fit.tra_data_sets),1):
            ntran_data += len(fit.tra_data_sets[i]) 

        if fit.type_fit["Transit"] == True  and ntran_data == 0:
             choice = QtWidgets.QMessageBox.information(self, 'Warning!',
             "Not possible to run %s if there are no transit data loaded. Please add your transit data first. Okay?"%tf, QtWidgets.QMessageBox.StandardButton.Ok)
           #  self.button_nest_samp.setEnabled(True)  
             self.statusBar().showMessage('') 

             return 0

        nttv_data = 0
        for i in range(0,10,1):
            nttv_data += len(fit.ttv_data_sets[i]) 

        if fit.type_fit["TTV"] == True  and nttv_data == 0:
             choice = QtWidgets.QMessageBox.information(self, 'Warning!',
             "Not possible to run %s if there are no TTV data loaded. Please add your TTV data first. Okay?"%tf, QtWidgets.QMessageBox.StandardButton.Ok)      
             #self.button_nest_samp.setEnabled(True)  
             self.statusBar().showMessage('') 

             return 0
            
        # check if astrometric data is present
        if fit.type_fit["AST"] == True:          
            if self.use_imaging_ast.isChecked():
                z=0
                for i in range(10):
                    if len(fit.ast_data_sets[i]) == 0:
                        continue
                    else:
                        z=z+1
                        if fit.ast_data_sets[i][5] > fit.npl and fit.ast_data_sets[i][6] == True:
                            choice = QtWidgets.QMessageBox.information(self, 'Warning!',"Astrometry dataset %s is set to planet %s but this planet is not included. Okay?"%(i+1,fit.ast_data_sets[i][5]), QtWidgets.QMessageBox.StandardButton.Ok)      
                            #self.button_fit.setEnabled(True)
                            self.mute_buttons(trigger=True)
                            return 0

                if z <= 0:
                    choice = QtWidgets.QMessageBox.information(self, 'Warning!',
                    "Not possible to model planets if there are no imaging astrometry data loaded. Please add your astrometry data first. Okay?", QtWidgets.QMessageBox.StandardButton.Ok)      
                    #self.button_fit.setEnabled(True)
                    self.mute_buttons(trigger=True)
                    return 0
                    
                    
            if self.use_hipp_gaia_ast.isChecked():
                z=0
                for i in range(10):
                    if len(fit.ast_data_sets_hipp_gaia[i]) == 0:
                        continue
                    else:
                        z=z+1
                        if fit.ast_data_sets_hipp_gaia[i][7] > fit.npl and fit.ast_data_sets_hipp_gaia[i][8] == True:
                            choice = QtWidgets.QMessageBox.information(self, 'Warning!',"Astrometry dataset %s is set to planet %s but this planet is not included. Okay?"%(i+1,fit.ast_data_sets_hipp_gaia[i][7]), QtWidgets.QMessageBox.StandardButton.Ok)      
                            #self.button_fit.setEnabled(True)
                            self.mute_buttons(trigger=True)
                            return 0

                if z <= 0:
                    choice = QtWidgets.QMessageBox.information(self, 'Warning!',
                    "Not possible to model planets if there are no Hipparcos/Gaia astrometry data loaded. Please add your astrometry data first. Okay?", QtWidgets.QMessageBox.StandardButton.Ok)      
                    #self.button_fit.setEnabled(True)
                    self.mute_buttons(trigger=True)
                    return 0 
 
            if fit.use_ast_hipp_gaia == False and fit.use_ast_imaging == False: 
            
                choice = QtWidgets.QMessageBox.information(self, 'Warning!',
                "Not possible to start astrometry fit if no Astrometry type is selected. Please select your astrometry data type (and load data) first. Okay?", QtWidgets.QMessageBox.StandardButton.Ok)      
                #self.button_fit.setEnabled(True)
                self.mute_buttons(trigger=True)
                return 0    
                
        return 1 
        
        
                       

    def worker_nest(self):
        global fit  

        self.update_params()
        self.update_use()

        if self.radioButton_RV.isChecked():
            fit.rtg = [True,self.do_RV_GP.isChecked(), False, False]
        elif self.radioButton_transit.isChecked():
            fit.rtg = [False, False, True, self.do_tra_GP.isChecked()]
        elif self.radioButton_transit_RV.isChecked():
            fit.rtg = [True,self.do_RV_GP.isChecked(), True, self.do_tra_GP.isChecked()]

        elif self.radioButton_ttv.isChecked():
            fit.rtg = [False, False, False, False]
        elif self.radioButton_ttv_RV.isChecked():
            fit.rtg = [True,self.do_RV_GP.isChecked(), False, False] 

        fit.nest_mad = self.use_nest_MAD_level.isChecked()

        if self.use_NS_stab_constraints.isChecked():
            fit.NS_AMD_stab   = self.NS_AMD.isChecked()
            fit.NS_Nbody_stab = self.NS_Nbody.isChecked()
        else:
            fit.NS_AMD_stab   = False
            fit.NS_Nbody_stab = False


        #self.button_nest_samp.setEnabled(False)
        self.statusBar().showMessage('Nested Sampling in progress....')
        # check if RV data is present
         
        
        if self.check_data_and_fit_type(tf="Nested Sampling") == 0:
            return


        choice = QtWidgets.QMessageBox.information(self, 'Warning!',
"""
This may take some time. Results are printed in the 'Stdout/Stderr' tab. Okay?


Also, did you setup your priors? By default, the Exo-Striker's priors are WIDELY open! Make sure your priors are 'reasonable' for your science case; otherwise, the Nested Sampling run may take forever!
""",
                                            QtWidgets.QMessageBox.StandardButton.Cancel | QtWidgets.QMessageBox.StandardButton.Ok)       

        if choice == QtWidgets.QMessageBox.StandardButton.Cancel:
            self.statusBar().showMessage('') 
           # self.button_nest_samp.setEnabled(True)
            return

        self.check_ttv_params()
        self.check_ast_params()
        self.set_tra_ld()
        self.check_bounds()
        self.check_priors_nr() 
        fit.model_npoints = self.points_to_draw_model.value()
        fit.model_max = self.model_max_range.value()
        fit.model_min = self.model_min_range.value()
        
        self.tabWidget_helper.setCurrentWidget(self.tab_info)

        if self.use_nest_percentile_level.isChecked():
            fit.nest_percentile_level = self.nest_percentile_level.value()
        else:
            fit.nest_percentile_level = 68.3

        self.mute_buttons(trigger=False)

        # Pass the function to execute
        self.worker_n = Worker(self.run_nest) # Any other args, kwargs are passed to the run  
        # Execute
        self.worker_n.signals.finished.connect(self.worker_nest_complete)

        # worker.signals.result.connect(self.print_output)
        #worker.signals.finished.connect(self.thread_complete)
       # worker.signals.progress.connect(self.progress_fn)
        self.threadpool.start(self.worker_n)

    def stop_nest(self):
        """Stop the worker process."""

        #if hasattr(self, 'worker_n'):
        self.worker_n.stop()
        self.statusBar().showMessage("Stopping Nested Sampling...")

    def run_nest(self):
        global fit

        self.check_model_params() 
        self.check_nested_params()
        self.check_settings()


        if self.run_ns_in_bg.isChecked():
            fit = rv.run_nestsamp_bg(fit)
        else:
            fit = rv.run_nestsamp(fit)

        self.button_nest_samp.setEnabled(True)

    def change_nest_samples_file_name(self):
        global fit
        
        output_file = QtWidgets.QFileDialog.getSaveFileName(self, 'path and name of the nested samples', '', '', options=QtWidgets.QFileDialog.Option.DontUseNativeDialog)

        if output_file[0] != '':
            fit.nest_sample_file = output_file[0] 
            self.nest_samples_change_name.setText(output_file[0])
        else:
            return

    def check_nested_params(self):
        global fit


        fit.live_points_fact = int(self.live_points.value())
        fit.ns_threads=int(self.nest_N_threads.value())
        fit.Dynamic_nest = self.radioButton_dyn_nest_samp.isChecked()
        fit.ns_progress = self.ns_progress.isChecked() 
        fit.stop_crit = self.stop_crit.value()
        fit.ns_fileoutput=self.save_samples_nested.isChecked()
        fit.ns_save_means=self.adopt_nest_means_as_par.isChecked() 
        fit.ns_save_median=self.adopt_nest_median_as_par.isChecked() 
        fit.ns_save_mode=self.adopt_nest_mode_as_par.isChecked() 
        fit.ns_save_maxlnL=self.adopt_nest_best_lnL_as_pars.isChecked() 
        fit.ns_save_sampler=self.save_samples_nested_in_memory.isChecked()

        fit.ns_use_stop = self.ns_use_stop.isChecked()
        fit.ns_maxiter = {0:self.use_ns_maxiter.isChecked(), 1:self.ns_maxiter.value()}
        fit.ns_maxcall = {0:self.use_ns_maxcall.isChecked(), 1:self.ns_maxcall.value()}

        fit.nest_weighted=self.nest_weighted.isChecked()


    def force_nest_check_box(self):
        if self.button_make_cornerplot_nested.isChecked():
            self.save_samples_nested.setChecked(True)
 
    def init_ns_samp_opt_combo(self):
        global fit

        for i in range(len(fit.ns_samp_method_opt)):
            self.comboBox_ns_samp_opt.addItem('%s'%(fit.ns_samp_method_opt[i]),i+1)

        for i in range(len(fit.ns_samp_bound_opt)):
            self.comboBox_ns_bound_opt.addItem('%s'%(fit.ns_samp_bound_opt[i]),i+1)


    def check_ns_samp_opt_combo(self):
        global fit

        ind_ns_opt = self.comboBox_ns_samp_opt.currentIndex()
        fit.ns_samp_method = fit.ns_samp_method_opt[ind_ns_opt]

        ind_ns_bound_opt = self.comboBox_ns_bound_opt.currentIndex()
        fit.ns_samp_bound = fit.ns_samp_bound_opt[ind_ns_bound_opt]

        fit.ns_pfrac = self.nest_pfrac.value()

        self.check_nested_params()



    def update_GUI_ns_params(self):
        global fit
  

        self.live_points.setValue(fit.live_points_fact)
        self.nest_N_threads.setValue(fit.ns_threads)
        self.radioButton_dyn_nest_samp.setChecked(fit.Dynamic_nest)
        self.stop_crit.setValue(fit.stop_crit)
        self.save_samples_nested.setChecked(fit.ns_fileoutput)
        self.adopt_nest_means_as_par.setChecked(fit.ns_save_means) 
        self.adopt_nest_median_as_par.setChecked(fit.ns_save_median) 
        self.adopt_nest_mode_as_par.setChecked(fit.ns_save_mode) 
        self.adopt_nest_best_lnL_as_pars.setChecked(fit.ns_save_maxlnL) 
        self.save_samples_nested_in_memory.setChecked(fit.ns_save_sampler)

        self.ns_use_stop.setChecked(fit.ns_use_stop)
        self.use_ns_maxiter.setChecked(fit.ns_maxiter[0])
        self.ns_maxiter.setValue(fit.ns_maxiter[1])
        self.use_ns_maxcall.setChecked(fit.ns_maxcall[0])
        self.ns_maxcall.setValue(fit.ns_maxcall[1])

        self.comboBox_ns_samp_opt.setCurrentIndex(fit.ns_samp_method_opt.index(fit.ns_samp_method))
        self.comboBox_ns_bound_opt.setCurrentIndex(fit.ns_samp_bound_opt.index(fit.ns_samp_bound))
 
        
        self.nest_pfrac.setValue(fit.ns_pfrac)
        self.nest_weighted.setChecked(fit.nest_weighted)

################################## MCMC #######################################

    def worker_mcmc_complete(self):
        global fit  

        self.update_labels()
        self.update_gui_params()
        self.update_errors() 
        self.update_a_mass() 

        self.init_plot_corr()

        self.statusBar().showMessage('') 

        if fit.bound_error == True:
            self.get_error_msg(fit.bound_error_msg)
            self.mute_buttons(trigger=True)
            return
        
        if self.adopt_mcmc_means_as_par.isChecked() or self.adopt_mcmc_median_as_par.isChecked() or self.adopt_best_lnL_as_pars.isChecked() or self.adopt_mcmc_mode_as_par.isChecked():
            self.init_fit(update_error = False)
        #else:
            
        self.jupiter_push_vars()
          
        #self.fit_dispatcher(init=True)          
        #self.jupiter_push_vars()
            
            

        self.save_last_session("autosave/auto_save.ses")
        self.check_cornerplot_samples()        
        self.mute_buttons(trigger=True)

    def worker_mcmc(self):
        global fit  

        self.update_params()
        self.update_use()

        if self.radioButton_RV.isChecked():
            fit.rtg = [True,self.do_RV_GP.isChecked(), False, False]
        elif self.radioButton_transit.isChecked():
            fit.rtg = [False, False, True, self.do_tra_GP.isChecked()]
        elif self.radioButton_transit_RV.isChecked():
            fit.rtg = [True,self.do_RV_GP.isChecked(), True, self.do_tra_GP.isChecked()]
            
        elif self.radioButton_ttv.isChecked():
            fit.rtg = [False, False, False, False]
        elif self.radioButton_ttv_RV.isChecked():
            fit.rtg = [True,self.do_RV_GP.isChecked(), False, False] 
 
        fit.mcmc_mad = self.use_mcmc_MAD_level.isChecked()
        
        if self.use_mcmc_stab_constraints.isChecked():
            fit.mcmc_AMD_stab   = self.mcmc_AMD.isChecked()
            fit.mcmc_Nbody_stab = self.mcmc_Nbody.isChecked()
        else:
            fit.mcmc_AMD_stab   = False
            fit.mcmc_Nbody_stab = False


        #self.button_MCMC.setEnabled(False)
        self.statusBar().showMessage('MCMC in progress....')

        if self.check_data_and_fit_type(tf="MCMC") == 0:
            return
            
            
        choice = QtWidgets.QMessageBox.question(self, 'Warning!',
                                       """This may take some time. Results are printed in the 'Stdout/Stderr' tab. Okay?
                                       """,
                                       QtWidgets.QMessageBox.StandardButton.Cancel | QtWidgets.QMessageBox.StandardButton.Ok)

        if choice == QtWidgets.QMessageBox.StandardButton.Cancel:
            self.statusBar().showMessage('') 
            self.button_MCMC.setEnabled(True)
            return


        self.check_ttv_params()
        self.check_ast_params()
        self.set_tra_ld()
        self.check_bounds()
        self.check_priors_nr() 
        
        fit.model_npoints = self.points_to_draw_model.value()
        fit.model_max = self.model_max_range.value()
        fit.model_min = self.model_min_range.value()
#        fit.dyn_model_to_kill = 
        
        self.tabWidget_helper.setCurrentWidget(self.tab_info)

        if self.use_percentile_level.isChecked():
            fit.percentile_level = self.percentile_level.value()
        else:
            fit.percentile_level = 68.3

        self.mute_buttons(trigger=False)


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
        self.check_settings()

        if self.run_mcmc_in_bg.isChecked():
            fit = rv.run_mcmc_bg(fit)
        else:
            fit = rv.run_mcmc(fit)

        self.button_MCMC.setEnabled(True)   


    def change_mcmc_samples_file_name(self):
        global fit
        
        output_file = QtWidgets.QFileDialog.getSaveFileName(self, 'path and name of the mcmc samples', '', '', options=QtWidgets.QFileDialog.Option.DontUseNativeDialog)
        
        if output_file[0] != '':
            fit.mcmc_sample_file = output_file[0] 
            self.mcmc_samples_change_name.setText(output_file[0])
        else:
            return

    def check_mcmc_params(self):
        global fit
        #print(fit.mcmc_fileoutput,self.save_samples.isChecked())

        #print(int(self.N_threads.value()))

        fit.gaussian_ball = self.init_gauss_ball.value() 
        fit.nwalkers_fact = int(self.nwalkers_fact.value()) 
        fit.mcmc_burning_ph = self.burning_phase.value() 
        fit.mcmc_ph = self.mcmc_phase.value() 
        fit.mcmc_threads=int(self.N_threads.value())
        fit.mcmc_fileoutput=self.save_samples.isChecked()
        fit.mcmc_save_means=self.adopt_mcmc_means_as_par.isChecked() 
        fit.mcmc_save_median=self.adopt_mcmc_median_as_par.isChecked() 
       # print(fit.mcmc_fileoutput,self.save_samples.isChecked())
        fit.mcmc_progress= self.mcmc_progress.isChecked() 

        fit.mcmc_save_mode=self.adopt_mcmc_mode_as_par.isChecked() 
        fit.mcmc_save_maxlnL=self.adopt_best_lnL_as_pars.isChecked() 
        fit.mcmc_save_sampler=self.save_samples_mcmc_in_memory.isChecked()     
        
    def check_model_params(self):
        global fit
        fit.time_step_model = self.time_step_model.value()
        fit.dyn_model_accuracy = self.dyn_model_accuracy.value()
        fit.master_timeout = self.master_timeout.value()   


    def force_mcmc_check_box(self):
        if self.make_corner_plot.isChecked():
            self.save_samples.setChecked(True)


    def update_GUI_mcmc_params(self):
        global fit
 
        self.init_gauss_ball.setValue(fit.gaussian_ball) 
        self.nwalkers_fact.setValue(fit.nwalkers_fact) 
        self.burning_phase.setValue(fit.mcmc_burning_ph) 
        self.mcmc_phase.setValue(fit.mcmc_ph) 
        self.N_threads.setValue(fit.mcmc_threads)
        #print(fit.mcmc_fileoutput,self.save_samples.isChecked())
       
        self.save_samples.setChecked(bool(fit.mcmc_fileoutput))
        self.adopt_mcmc_means_as_par.setChecked(fit.mcmc_save_means) 
        self.adopt_mcmc_median_as_par.setChecked(fit.mcmc_save_median) 
        #print(fit.mcmc_fileoutput,self.save_samples.isChecked())

        self.adopt_mcmc_mode_as_par.setChecked(fit.mcmc_save_mode) 
        self.adopt_best_lnL_as_pars.setChecked(fit.mcmc_save_maxlnL) 
        self.save_samples_mcmc_in_memory.setChecked(fit.mcmc_save_sampler)   

################################## Cornerplot #######################################

    def switch_to_corenrplot_opt(self):
        global fit  

        self.param_tabs.setCurrentWidget(self.plot_options_tabs)
        self.plot_opt_tab.setCurrentWidget(self.cornerplot_plot_tab)
        
        self.check_cornerplot_samples()

    def remove_mcmc_samples_from_fit(self):
        global fit  
        
        #if isinstance(fit.mcmc_sampler, rv.CustomSampler):
        if len(np.atleast_1d(fit.mcmc_sampler))!=0:          
            choice = QtWidgets.QMessageBox.information(self, 'Warning!',
                                            "Are you sure you want to remove the MCMC samples?",
                                             QtWidgets.QMessageBox.StandardButton.No | QtWidgets.QMessageBox.StandardButton.Yes)

            if choice == QtWidgets.QMessageBox.StandardButton.No:
                return   
            elif choice == QtWidgets.QMessageBox.StandardButton.Yes:
                del fit.mcmc_sampler
                fit.mcmc_sampler = []
                self.check_cornerplot_samples()
                return          
            else: 
                return   
        else:
             choice = QtWidgets.QMessageBox.information(self, 'Warning!',
             "MCMC samples not found.", QtWidgets.QMessageBox.StandardButton.Ok)            

 
    def remove_ns_samples_from_fit(self):
        global fit  

        #if len(fit.ns_sampler)!=0:     
        if len(np.atleast_1d(fit.ns_sampler))!=0:  
        #if isinstance(fit.ns_sampler, dynesty.nestedsamplers.UnitCubeSampler) or isinstance(fit.ns_sampler,dynesty.dynamicsampler.DynamicSampler):
            choice = QtWidgets.QMessageBox.information(self, 'Warning!',
                                            "Are you sure you want to remove the NS samples?",
                                             QtWidgets.QMessageBox.StandardButton.No | QtWidgets.QMessageBox.StandardButton.Yes)

            if choice == QtWidgets.QMessageBox.StandardButton.No:
                return   
            elif choice == QtWidgets.QMessageBox.StandardButton.Yes:
                del fit.ns_sampler
                fit.ns_sampler = []
                self.check_cornerplot_samples()
                return          
            else: 
                return   
        else:
             choice = QtWidgets.QMessageBox.information(self, 'Warning!',
             "NS samples not found.", QtWidgets.QMessageBox.StandardButton.Ok)             

    def check_cornerplot_samples(self):
        global fit  
      
       # import dynesty
        #if isinstance(fit.mcmc_sampler, rv.CustomSampler):
        if len(np.atleast_1d(fit.mcmc_sampler))!=0:        
            MCMC_SAMP_LED = './lib/UI/green_led.png'
            MCMC_SAMP_TXT = "MCMC samples available" 
        else:
            MCMC_SAMP_LED = './lib/UI/red_led.png' 
            MCMC_SAMP_TXT = "No MCMC samples available" 
#dynesty.nestedsamplers.UnitCubeSampler
        #if len(fit.ns_sampler)!=0:
        if len(np.atleast_1d(fit.ns_sampler))!=0:
        #if isinstance(fit.ns_sampler, dynesty.nestedsamplers.UnitCubeSampler) or isinstance(fit.ns_sampler,dynesty.dynamicsampler.DynamicSampler):
            NS_SAMP_LED = './lib/UI/green_led.png'
            NS_SAMP_TXT = "NS samples available" 
        else:
            NS_SAMP_LED = './lib/UI/red_led.png' 
            NS_SAMP_TXT = "No NS samples available" 
            
        self.mcmc_samples_led.setPixmap(QtGui.QPixmap(MCMC_SAMP_LED)) 
        self.mcmc_samples_led_2.setPixmap(QtGui.QPixmap(MCMC_SAMP_LED))         
        self.mcmc_cornerplot_samp_indicator.setText(MCMC_SAMP_TXT)
        
        self.ns_samples_led.setPixmap(QtGui.QPixmap(NS_SAMP_LED)) 
        self.ns_samples_led_2.setPixmap(QtGui.QPixmap(NS_SAMP_LED))         
        self.ns_cornerplot_samp_indicator.setText(NS_SAMP_TXT)        

    def worker_cornerplot_complete(self):
        global fit  
        self.statusBar().showMessage('') 
        self.button_make_mcmc_cornerplot.setEnabled(True)
        self.button_make_nest_cornerplot.setEnabled(True)

        self.tabWidget_helper.setCurrentWidget(self.tab_info)

    def worker_cornerplot(self, type_plot = "mcmc"):
        global fit  


        if type_plot == "mcmc":
           # if isinstance(fit.mcmc_sampler, rv.CustomSampler)==False:
            if len(np.atleast_1d(fit.mcmc_sampler))==0:            
                choice = QtWidgets.QMessageBox.information(self, 'Warning!', "MCMC samples not found.", QtWidgets.QMessageBox.StandardButton.Ok)
                return   
       
        if type_plot == "nest": 
            if len(np.atleast_1d(fit.ns_sampler))==0:
            #if isinstance(fit.ns_sampler, dynesty.nestedsamplers.UnitCubeSampler)==False and isinstance(fit.ns_sampler,dynesty.dynamicsampler.DynamicSampler)==False:
                choice = QtWidgets.QMessageBox.information(self, 'Warning!', "NS samples not found.", QtWidgets.QMessageBox.StandardButton.Ok)
                return   
  

        self.button_make_mcmc_cornerplot.setEnabled(False)
        self.button_make_nest_cornerplot.setEnabled(False)

        self.statusBar().showMessage('Cornerplot in progress....')
        
        
        

        # check if RV data is present
     #   if type_plot == "mcmc":
     #       samp_file = fit.mcmc_sample_file
     #       type_samp = "MCMC"
     #   elif type_plot == "nest":
     #       samp_file = fit.nest_sample_file
     #       type_samp = "Nest. Samp."

     #   if not os.path.exists(samp_file):
     #        choice = QtWidgets.QMessageBox.information(self, 'Warning!',
     #        "%s file not found. Generate one and try again?"%type_samp, QtWidgets.QMessageBox.StandardButton.Ok)
     #        self.button_make_mcmc_cornerplot.setEnabled(True)
     #        self.button_make_nest_cornerplot.setEnabled(True)
     #        self.statusBar().showMessage('')
     #        return


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
        rv.cornerplot(fit, type_plot=type_plot )


    def change_corner_plot_file_name(self, type_plot = "mcmc"):
        global fit
        
        output_file = QtWidgets.QFileDialog.getSaveFileName(self, 'path and name of the corener plot', '', 'All (*.*);;Data (*.pdf);;Data (*.png)', options=QtWidgets.QFileDialog.Option.DontUseNativeDialog)
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

        #print(self.RVBank)   
        
        if self.RVBank == False:
            path = self.inspector_file
            if path == '':
                return 

            filename, file_extension = os.path.splitext(path)

            if file_extension == '.vels' or file_extension == '.dat':
                fit.add_dataset(self.file_from_path(path), str(path),0.0,1.0)
                self.init_fit()
                self.update_use_from_session()
                self.update_use()
                self.update_params()
                self.update_RV_file_buttons()

            elif file_extension == '.act':

                for i in range(20):
                    if len(fit.act_data_sets[i]) == 0:
                        but_ind = i +1
                        fit.add_act_dataset('test', str(path),act_idset =but_ind-1)

                        self.update_act_file_buttons()
                        self.update_activity_gls_plots(but_ind-1)
                        return

            elif  file_extension == '.tran':
                for i in range(20):
                    if len(fit.tra_data_sets[i]) == 0:
                        but_ind = i +1

                        fit.add_transit_dataset('test', str(path),tra_idset =but_ind-1)
                        self.update_use_from_session()
                        self.update_use()
                        self.update_gui_params()
                        self.update_params()
                        self.update_tra_file_buttons()
                        self.buttonGroup_transit_data.button(but_ind).setText(self.file_from_path(path))
                        return
            else: 
                return
        else:


            if self.RVBank_window.data_index < 11 and self.RVBank_window.type_data == "HARPS":
                
                BJD       = self.RVBank_window.x_data 
                rv_data     = self.RVBank_window.y_data
                rv_data_sig = self.RVBank_window.e_y_data
 
                name1 = '%s_pre.dat'%self.RVBank_window.target_name
                name2 = '%s_post.dat'%self.RVBank_window.target_name
                path1 = 'datafiles/%s'%name1
                path2 = 'datafiles/%s'%name2
        
                out1 = open('%s'%path1, 'w')
                out2 = open('%s'%path2, 'w')

                for i in range(len(BJD)):

                    if float(BJD[i]) <= 2457161.5:
                        out1.write('{0:{width}.{precision}f}  {1:{width}.{precision}f}  {2:{width}.{precision}f}  \n'.format(float(BJD[i]), float(rv_data[i]), float(rv_data_sig[i]),  width = 10, precision = 5 )   )
                    elif float(BJD[i]) > 2457161.5:
                        out2.write('{0:{width}.{precision}f}  {1:{width}.{precision}f}  {2:{width}.{precision}f}  \n'.format(float(BJD[i]), float(rv_data[i]), float(rv_data_sig[i]),  width = 10, precision = 5 )   )

                out1.close()
                out2.close()

                if len(BJD[BJD <= 2457161.5]) !=0:
                    fit.add_dataset(name1,path1,0.0,1.0,useoffset=True,usejitter=True)
                if len(BJD[BJD > 2457161.5]) !=0:
                    fit.add_dataset(name2,path2,0.0,1.0,useoffset=True,usejitter=True)
 
    
                fit.type_fit["RV"] = True
                fit.type_fit["Transit"] = False
                self.check_type_fit()
                self.mute_boxes()
    
                self.init_fit()
               # rv.check_temp_RV_file(fit)
                #self.update_veiw()
                self.update_use_from_session()
                self.update_use()
                self.update_params()
                self.update_RV_file_buttons()
                
            elif self.RVBank_window.data_index >= 11 and self.RVBank_window.type_data == "HARPS":
                
                act_JD       = self.RVBank_window.x_data

                if self.RVBank_window.data_index != 53:
#                    act_data     = np.concatenate((
#                self.RVBank_window.y_data[act_JD <= 2457161.5]- np.average(self.RVBank_window.y_data[act_JD <= 2457161.5], weights=1.0/abs(self.RVBank_window.e_y_data[act_JD <= 2457161.5])),  
#                self.RVBank_window.y_data[act_JD > 2457161.5] - np.average(self.RVBank_window.y_data[act_JD  > 2457161.5], weights=1.0/abs(self.RVBank_window.e_y_data[act_JD > 2457161.5]))))


                    # Define a small epsilon to avoid division by zero in weights
                    epsilon = 1e-10

                    # Mask for act_JD <= 2457161.5
                    mask1 = act_JD <= 2457161.5
                    if np.any(mask1):
                        weights1 = 1.0 / (abs(self.RVBank_window.e_y_data[mask1]) + epsilon)
                        avg1 = np.average(self.RVBank_window.y_data[mask1], weights=weights1)
                        data1 = self.RVBank_window.y_data[mask1] - avg1
                    else:
                        data1 = np.array([])  # If no data, create an empty array

                    # Mask for act_JD > 2457161.5
                    mask2 = act_JD > 2457161.5
                    if np.any(mask2):
                        weights2 = 1.0 / (abs(self.RVBank_window.e_y_data[mask2]) + epsilon)
                        avg2 = np.average(self.RVBank_window.y_data[mask2], weights=weights2)
                        data2 = self.RVBank_window.y_data[mask2] - avg2
                    else:
                        data2 = np.array([])  # If no data, create an empty array

                    # Concatenate the two parts
                    act_data = np.concatenate((data1, data2))

 
                
                else:
                    act_data = self.RVBank_window.y_data                    

                act_data_sig = self.RVBank_window.e_y_data

                act_file_name = self.RVBank_window.data_name
                act_data_o_c = dill.copy(act_data)            
                act_data_set = np.array([act_JD,act_data,act_data_sig,act_data_o_c,act_data_o_c,act_data,act_data_sig,act_data_o_c, 1.0, act_file_name],dtype=object) 

                for i in range(20):
                    if len(fit.act_data_sets[i]) == 0:
                        fit.act_data_sets[i]      = dill.copy(act_data_set)
                        fit.act_data_sets_init[i] = dill.copy(fit.act_data_sets[i])
                        break

                self.update_act_file_buttons()



            elif self.RVBank_window.data_index < 5 and self.RVBank_window.type_data == "HIRES":
                
                BJD       = self.RVBank_window.x_data 
                rv_data     = self.RVBank_window.y_data
                rv_data_sig = self.RVBank_window.e_y_data
 
    
                name1 = '%s.dat'%self.RVBank_window.target_name
                path1 = 'datafiles/%s'%name1
                out1 = open('%s'%path1, 'w')

                for i in range(len(BJD)):
     
                    out1.write('{0:{width}.{precision}f}  {1:{width}.{precision}f}  {2:{width}.{precision}f}  \n'.format(float(BJD[i]), float(rv_data[i]), float(rv_data_sig[i]),  width = 10, precision = 5 )   )

                out1.close()
                fit.add_dataset(name1,path1,0.0,1.0,useoffset=True,usejitter=True)


                self.init_fit()

                self.update_use_from_session()
                self.update_use()
                self.update_params()
                self.update_RV_file_buttons()
                

                
            elif self.RVBank_window.data_index >= 5 and self.RVBank_window.type_data == "HIRES":
                
                act_JD       = self.RVBank_window.x_data
                act_data     = self.RVBank_window.y_data
                act_data_sig = self.RVBank_window.e_y_data

                act_file_name = self.RVBank_window.data_name

                #act_data_set = np.array([act_JD,act_data,act_data_sig,act_file_name])
                act_data_o_c = dill.copy(act_data)            
                act_data_set = np.array([act_JD,act_data,act_data_sig,act_data_o_c,act_data_o_c,act_data,act_data_sig,act_data_o_c, 1.0, act_file_name],dtype=object)

                for i in range(10):
                    if len(fit.act_data_sets[i]) == 0:
                        fit.act_data_sets[i]      = dill.copy(act_data_set)
                        fit.act_data_sets_init[i] = dill.copy(fit.act_data_sets[i])
                        break

                self.update_act_file_buttons()

    def save_DataInspector_data(self):
        global fit, pdi 

        if self.RVBank == True:
            print("Exporting data for %s"%self.RVBank_window.target_name)
            pdi_data_name = self.RVBank_window.data_name
            pdi_target_name = self.RVBank_window.target_name
        else:
            path = self.datafiles_window.listview.model().filePath(self.datafiles_window.listview.currentIndex()) 
            pdi_data_name = path
            pdi_target_name = path

        output_file = QtWidgets.QFileDialog.getSaveFileName(self, 'Save to data (ascii) file', '%s_%s.dat'%(pdi_target_name,pdi_data_name), 'All (*.*);;', options=QtWidgets.QFileDialog.Option.DontUseNativeDialog)

        try:         
            x_pdi = dill.copy(pdi.plotItem.items[3].opts["x"])
            y_pdi = dill.copy(pdi.plotItem.items[3].opts["y"])
            e_pdi = dill.copy(pdi.plotItem.items[3].opts["top"])
            
            if len(x_pdi)==0:
                print("Exporting an empty file makes no sense.1")
                return

            if str(output_file[0]) != '':
                f = open(output_file[0], 'w')
                for i in range(len(x_pdi)):
                    f.write('{0:{width}.{precision}f}  {1:{width}.{precision}f}  {2:{width}.{precision}f} \n'.format(
                            float(x_pdi[i]), 
                            float(y_pdi[i]), 
                            float(e_pdi[i]),   
                            width = 14, precision = 7 ))
                f.close()

        except:
            print("Exporting an empty file makes no sense.2")
            return

        return



    def cross_data_inspect(self):
        if self.inpector_plot_cross_hair.isChecked():
            self.cross_hair(pdi,log=False)   
        else:
            self.cross_hair_remove(pdi)
            

    def init_plot_data_inspect(self):
        
        if self.RVBank == True:
            self.plot_data_inspect(0,no_sender=True,RVBank=True)
        else:
            self.plot_data_inspect(0,no_sender=True,RVBank=False)


    def plot_data_inspect(self, index, no_sender=False, RVBank = False):
        global fit, colors, pdi 

        self.RVBank = RVBank
        
        if no_sender==True:
            self.cross_data_inspect()
            return
        else:
            pdi.plot(clear=True,)

        if RVBank == False:

            path = self.datafiles_window.listview.model().filePath(self.datafiles_window.listview.currentIndex()) 
            
            if os.path.exists(path) and os.path.getsize(path) > 0:
                try:
                    x     = np.genfromtxt("%s"%(path),skip_header=0,comments="#", unpack=True,skip_footer=0, usecols = [0])
                    y     = np.genfromtxt("%s"%(path),skip_header=0,comments="#", unpack=True,skip_footer=0, usecols = [1])
                    y_err = np.genfromtxt("%s"%(path),skip_header=0,comments="#", unpack=True,skip_footer=0, usecols = [2])
                    
                    if np.isnan(x).any() or np.isnan(y).any() or np.isnan(y_err).any():
                        print("wrong file format")
                        return
                except (ValueError,Exception) as ex:
                    #print(ex)
                    #pdi.plot(clear=True,)
                    pdi.setLabel('bottom', 'x', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})
                    pdi.setLabel('left',   'y', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})
                    return
            else:
                print("%s is an empty file"%path)
                return
        else:
            
            if self.RVBank_window.url_success == False:
                return
            
            path = self.RVBank_window.path

            try:
                x     = self.RVBank_window.x_data
                y     = self.RVBank_window.y_data  
                y_err = self.RVBank_window.e_y_data 
            except:
                pdi.setLabel('bottom', 'x', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})
                pdi.setLabel('left',   'y', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})
                return
             
        
        #fit.pdi = pdi

        pdi.addLine(x=None, y=np.mean(y), pen=pg.mkPen('#ff9933', width=0.8),name="zero")

        if self.insp_data_size.value() > 2:
            symbolsize = self.insp_data_size.value() -2
        else:
            symbolsize = self.insp_data_size.value() 

        pdi.plot(x,y,
        pen=None, #{'color': colors[i], 'width': 1.1},
        symbol='o',
        symbolPen={'color': dill.copy(fit.colors[0]), 'width': 1.1},
        symbolSize=symbolsize,enableAutoRange=True,viewRect=True,
        symbolBrush=dill.copy(fit.colors[0])
        )

        err_ = pg.ErrorBarItem(x=x, y=y, symbol='o', 
        top = y_err, bottom = y_err,
        #height=y_err, 
        beam=0.0, pen=dill.copy(fit.colors[0]))   
     
        pdi.addItem(err_)
        pdi.autoRange()


        filename, file_extension = os.path.splitext(path)  
            
        if file_extension == '.vels':
            pdi.setLabel('bottom', 'BJD', units='d',  **{'font-size':'%dpt'%self.plot_font.pointSize()})
            pdi.setLabel('left',   'RV', units='m/s',  **{'font-size':'%dpt'%self.plot_font.pointSize()})
 
        elif file_extension == '.act':
            pdi.setLabel('bottom', 'BJD', units='d',  **{'font-size':'%dpt'%self.plot_font.pointSize()})
            pdi.setLabel('left',   'y', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})
            
        elif file_extension == '.tran':      
            pdi.setLabel('bottom', 'BJD', units='d',  **{'font-size':'%dpt'%self.plot_font.pointSize()})
            pdi.setLabel('left',   'flux', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})

        else:
            pdi.setLabel('bottom', 'x', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})
            pdi.setLabel('left',   'y', units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()})

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

        self.check_model_params()
        self.check_mcmc_params()
        self.check_nested_params()
        self.check_settings()

        if self.radioButton_RV.isChecked():
            fit.rtg = [True,self.do_RV_GP.isChecked(),False, False]
            if(init):
                self.worker_RV_fitting(ff=0,m_ln=True, init = init )
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

        elif self.radioButton_ttv.isChecked():

            fit.rtg=[False,False,False,False]
            if(init):
                self.worker_ttv_fitting(ff=0 )
            else:
                self.worker_ttv_fitting()
 
        elif self.radioButton_ttv_RV.isChecked():
            
            fit.rtg = [True,self.do_RV_GP.isChecked(),False, False]
            if(init):
                self.worker_ttv_fitting(ff=0 )  
            else:
                self.worker_ttv_fitting()


        elif self.radioButton_ast.isChecked():

            fit.rtg=[False,False,False,False]
            if(init):
                self.worker_ast_fitting(ff=0 )
            else:
                self.worker_ast_fitting()
   
        elif self.radioButton_ast_RV.isChecked():
            
            fit.rtg = [True,self.do_RV_GP.isChecked(),False, False]
            if(init):
                self.worker_ast_fitting(ff=0 )  
            else:
                self.worker_ast_fitting()

            #return   
###########################  GUI events #############################


    def mute_boxes(self):
        global fit

#        self.get_jupyter_vars()
        
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

        if ttvfast_not_found == True:
            self.radioButton_ttv.setEnabled(False)
            self.radioButton_ttv_RV.setEnabled(False)
            print("""
You dont have TTVfast installed! Therefore, you cannot apply TTV modelling. 
Please install via 'pip install ttvfast'.
""")
            self.tabWidget_helper.setCurrentWidget(self.tab_info)



        ######### TESTS!!!!!!!!!!!###########

        if self.radioButton_RV.isChecked():
            self.label.setEnabled(True)
            self.amoeba_radio_button.setEnabled(True)
            self.lm_radio_button.setEnabled(True)
        else:
            self.label.setEnabled(False)
            self.amoeba_radio_button.setChecked(True)
            self.amoeba_radio_button.setEnabled(False)
            self.lm_radio_button.setEnabled(False)

        if self.radioButton_transit_RV.isChecked():
         
            K_flag = True
            pl_rad_flag = True
            a_sol_flag = True
            fit.type_fit["RV"] = True
            fit.type_fit["Transit"] = True
            fit.type_fit["TTV"] = False
            fit.type_fit["AST"] = False          

            if self.radioButton_Dynamical.isChecked():
                ma_flag = True
                t0_flag = False
                incl_flag = True
                Dom_flag = True
            else:
                ma_flag = False
                t0_flag = True   
                incl_flag = True
                Dom_flag = False
         

        elif self.radioButton_transit.isChecked():

            K_flag = False
            pl_rad_flag = True
            a_sol_flag = True
            fit.type_fit["RV"] = False
            fit.type_fit["Transit"] = True 
            fit.type_fit["TTV"] = False
            fit.type_fit["AST"] = False          

            if self.radioButton_Dynamical.isChecked():
                ma_flag = True
                t0_flag = False
                incl_flag = True
                Dom_flag = True
            else:
                ma_flag = False
                t0_flag = True   
                incl_flag = True
                Dom_flag = False
            
                    
                    
        elif self.radioButton_RV.isChecked():

            K_flag = True
            pl_rad_flag = False
            a_sol_flag = False
            fit.type_fit["RV"] = True
            fit.type_fit["Transit"] = False
            fit.type_fit["TTV"] = False
            fit.type_fit["AST"] = False          
            
            if self.radioButton_Dynamical.isChecked():
                ma_flag = True
                t0_flag = False
                incl_flag = True
                Dom_flag = True
            else:
                ma_flag = True
                t0_flag = False   
                incl_flag = False
                Dom_flag = False
                
        elif self.radioButton_ttv.isChecked() or self.radioButton_ttv_RV.isChecked():

            K_flag = True
            pl_rad_flag = False
            a_sol_flag = False
            if self.radioButton_ttv.isChecked():
                fit.type_fit["RV"] = False
            elif self. radioButton_ttv_RV.isChecked():
                fit.type_fit["RV"] = True
            fit.type_fit["Transit"] = False 
            fit.type_fit["TTV"] = True 
            fit.type_fit["AST"] = False          
            
            if self.radioButton_Dynamical.isChecked():
                ma_flag = True
                t0_flag = False
                incl_flag = True
                Dom_flag = True
            else:
                ma_flag = True
                t0_flag = False
                incl_flag = True
                Dom_flag = True    
      
        elif self.radioButton_ast.isChecked() or self.radioButton_ast_RV.isChecked():

            K_flag = True
            pl_rad_flag = False
            a_sol_flag = False
            if self.radioButton_ast.isChecked():
                fit.type_fit["RV"] = False
            elif self.radioButton_ast_RV.isChecked():
                fit.type_fit["RV"] = True
            fit.type_fit["Transit"] = False 
            fit.type_fit["TTV"] = False 
            fit.type_fit["AST"] = True          
            
            if self.radioButton_Dynamical.isChecked():
                ma_flag = True
                t0_flag = False
                incl_flag = True
                Dom_flag = True
            else:
                ma_flag = True
                t0_flag = False
                incl_flag = True
                Dom_flag = True     

        om_flag = False

        if fit.type_fit["RV"] == True and self.radioButton_Keplerian.isChecked():             
            if self.deattach_omega_dot.isChecked():
                self.set_gr_flag()
                if not self.radioButton_omega_dot_GR.isChecked(): 
                    om_flag = True
 
        for i in range(9):
            self.param_gui_K[i].setEnabled(K_flag)
            self.use_param_gui_K[i].setEnabled(K_flag)
            self.param_gui_ma[i].setEnabled(ma_flag)
            self.use_param_gui_ma[i].setEnabled(ma_flag)

            self.param_gui_incl[i].setEnabled(incl_flag)
            self.use_param_gui_incl[i].setEnabled(incl_flag)
            self.param_gui_Omega[i].setEnabled(Dom_flag)
            self.use_param_gui_Omega[i].setEnabled(Dom_flag)
            self.param_gui_wd[i].setEnabled(om_flag)
            self.use_param_gui_wd[i].setEnabled(om_flag)

            self.param_gui_tr[3*i].setEnabled(t0_flag)
            self.use_param_gui_tr[3*i].setEnabled(t0_flag)
            self.param_gui_tr[3*i + 1].setEnabled(pl_rad_flag)
            self.use_param_gui_tr[3*i + 1].setEnabled(pl_rad_flag)
            self.param_gui_tr[3*i + 2].setEnabled(a_sol_flag)
            self.use_param_gui_tr[3*i + 2].setEnabled(a_sol_flag)

        #self.mute_boxes_dyn()

        if self.force_copl_incl.isChecked()==True and self.radioButton_Dynamical.isChecked():
            incl_flag = False
            Dom_flag = False     
            
            for i in range(9):
                self.param_gui_incl[i].setEnabled(incl_flag)
                self.use_param_gui_incl[i].setEnabled(incl_flag)
                self.param_gui_Omega[i].setEnabled(Dom_flag)
                self.use_param_gui_Omega[i].setEnabled(Dom_flag)
                self.param_gui_wd[i].setEnabled(om_flag)
                self.use_param_gui_wd[i].setEnabled(om_flag)

            self.param_gui_incl[0].setEnabled(True)
            self.use_param_gui_incl[0].setEnabled(True)
            
            
        if self.radioButton_transit_RV.isChecked() or self.radioButton_transit.isChecked():
            
            if self.radioButton_rho.isChecked():
                a_sol_flag = False                
            
                for i in range(9):
                    self.param_gui_tr[3*i + 2].setEnabled(False)
                    self.use_param_gui_tr[3*i + 2].setEnabled(False)            
                self.param_gui_tr[3*0 + 2].setEnabled(True)
                self.use_param_gui_tr[3*0 + 2].setEnabled(True)              
            
        if self.force_sameRV_jit.isChecked()==True:
            print("Use only RV jitter boolean for data 1, and mute the rest!")
        


    ### NOT used anymore 
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
            if self.radioButton_ttv.isChecked()==True or self.radioButton_ttv_RV.isChecked()==True:
                Dom_flag = True
            else:
                Dom_flag = False

        elif self.radioButton_Dynamical.isChecked():

            om_flag = False
            incl_flag = True
            Dom_flag = True 

            if self.radioButton_transit_RV.isChecked() or self.radioButton_transit.isChecked():
                ma_flag = True
                t0_flag = False
            else:
                
                for i in range(9):
                    self.param_gui_ma[i].setEnabled(ma_flag)
                    self.use_param_gui_ma[i].setEnabled(ma_flag)

                    self.param_gui_tr[3*i].setEnabled(t0_flag)
                    self.use_param_gui_tr[3*i].setEnabled(t0_flag)                
                              

        for i in range(9):
            self.param_gui_incl[i].setEnabled(incl_flag)
            self.use_param_gui_incl[i].setEnabled(incl_flag)
            self.param_gui_Omega[i].setEnabled(Dom_flag)
            self.use_param_gui_Omega[i].setEnabled(Dom_flag)

            self.param_gui_wd[i].setEnabled(om_flag)
            self.use_param_gui_wd[i].setEnabled(om_flag)
 


     

    def keyPressEvent(self, event):
        global fit
        
        
        if event.key() in (QtCore.Qt.Key_Enter, QtCore.Qt.Key_Return):
            
            if self.safe_to_init == True:
                self.update_veiw()
            else:
                print("You cannot initialize while a thread is working! This could cause a mess!")
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


    def set_gr_flag(self):
        global fit
        
        if self.radioButton_omega_dot_free.isChecked():
            fit.gr_flag = False
        elif self.deattach_omega_dot.isChecked() and self.radioButton_omega_dot_GR.isChecked():
            fit.gr_flag = True
        else:
            fit.gr_flag = False

############################# Tab selector (not ready) ################################

    def tab_selected(self,ind):

        if ind == 4:
            self.update_activity_data_plots(self.comboBox_act_data.currentIndex())
        if ind == 5:
            self.update_correlations_data_plots()


#############################  Color control ################################  

    ### RV ####

    def update_color_picker(self):
        global fit

        font = QtGui.QFont()
        font.setPointSize(9)
        font.setBold(False)
        #font.setWeight(75)

        for i in range(21):
            self.buttonGroup_color_picker.button(i+1).setStyleSheet("color: %s;"%dill.copy(fit.rvs_colors[i]))
            self.buttonGroup_color_picker.button(i+1).setText("%s"%dill.copy(fit.rvs_colors[i]))
            self.buttonGroup_color_picker.button(i+1).setFont(font)
        for i in range(20):    
            self.buttonGroup_symbol_picker.button(i+1).setStyleSheet("color: %s;"%dill.copy(fit.rvs_colors[i]))  
            self.buttonGroup_symbol_picker.button(i+1).setText(dill.copy(fit.pyqt_symbols_rvs[i])) 
            self.buttonGroup_symbol_picker.button(i+1).setFont(font)
            
          
    def get_color(self):
        global fit

        but_ind = self.buttonGroup_color_picker.checkedId()
        colorz = self.colorDialog.getColor(options=QtWidgets.QColorDialog.DontUseNativeDialog) #|QtWidgets.QColorDialog.ShowAlphaChannel,)

#        print(dir(colorz))
#        print(colorz.name())
#        print(colorz.alpha())
#        print(colorz.rgba())
        #print(colorz.getRgb())

        #QtWidgets.QColorDialog.setOption(QtWidgets.QColorDialog.ShowAlphaChannel,True)
        #colorz = QtWidgets.QColorDialog.getColor()

        if colorz.isValid():
            fit.rvs_colors[but_ind-1]=dill.copy(colorz.name())   #[0]+"B3"+colorz.name()[1:]
            self.update_color_picker()
            self.update_act_file_buttons()
            self.update_RV_file_buttons() 
            self.update_tra_file_buttons() 
            self.update_ttv_file_buttons()

            self.update_RV_plots() 
            self.update_extra_plots()
            self.update_extra_plots2()
            self.update_transit_plots() 
            #self.update_activity_data_plots() 
            #self.update_activity_gls_plots()
        else:
            return

    ### transit ####

    def update_color_picker_tra(self):
        global fit

        font = QtGui.QFont()
        font.setPointSize(9)
        font.setBold(False)
        #font.setWeight(75)

        for i in range(21):
            self.buttonGroup_color_picker_tra.button(i+1).setStyleSheet("color: %s;"%fit.tra_colors[i])
            self.buttonGroup_color_picker_tra.button(i+1).setText("%s"%fit.tra_colors[i])
            self.buttonGroup_color_picker_tra.button(i+1).setFont(font)
        for i in range(20):
            self.buttonGroup_symbol_picker_tra.button(i+1).setStyleSheet("color: %s;"%fit.tra_colors[i])
            self.buttonGroup_symbol_picker_tra.button(i+1).setText(dill.copy(fit.pyqt_symbols_tra[i]))
            self.buttonGroup_symbol_picker_tra.button(i+1).setFont(font)

    def get_color_tra(self):
        global fit

        but_ind = self.buttonGroup_color_picker_tra.checkedId()
        colorz = self.colorDialog.getColor(options=QtWidgets.QColorDialog.DontUseNativeDialog)

        #QtWidgets.QColorDialog.setOption(QtWidgets.QColorDialog.ShowAlphaChannel,True)
        #colorz = QtWidgets.QColorDialog.getColor()
        #print(colorz.name())
        if colorz.isValid():
            fit.tra_colors[but_ind-1]=colorz.name()
            self.update_color_picker_tra()
            self.update_act_file_buttons()
            self.update_RV_file_buttons()
            self.update_tra_file_buttons()
            self.update_ttv_file_buttons()

            self.update_RV_plots()
            self.update_extra_plots()
            self.update_extra_plots2()
            self.update_transit_plots()
            #self.update_activity_data_plots()
            #self.update_activity_gls_plots()
        else:
            return

    ### TTV ####

    def update_color_picker_ttv(self):
        global fit

        font = QtGui.QFont()
        font.setPointSize(11)
        font.setBold(False)
        #font.setWeight(75)

        for i in range(11):
            self.buttonGroup_color_picker_ttv.button(i+1).setStyleSheet("color: %s;"%dill.copy(fit.ttv_colors[i]))
            self.buttonGroup_color_picker_ttv.button(i+1).setFont(font)
        for i in range(10):
            self.buttonGroup_symbol_picker_ttv.button(i+1).setStyleSheet("color: %s;"%dill.copy(fit.ttv_colors[i]))
            self.buttonGroup_symbol_picker_ttv.button(i+1).setText(dill.copy(fit.pyqt_symbols_ttv[i]))
            self.buttonGroup_symbol_picker_ttv.button(i+1).setFont(font)

    def get_color_ttv(self):
        global fit

        but_ind = self.buttonGroup_color_picker_ttv.checkedId()
        colorz = self.colorDialog.getColor(options=QtWidgets.QColorDialog.DontUseNativeDialog)

        #QtWidgets.QColorDialog.setOption(QtWidgets.QColorDialog.ShowAlphaChannel,True)
        #colorz = QtWidgets.QColorDialog.getColor()

        if colorz.isValid():
            fit.ttv_colors[but_ind-1]=colorz.name()
            self.update_color_picker_ttv()
           # self.update_act_file_buttons()
          #  self.update_RV_file_buttons()
         #   self.update_tra_file_buttons()
            self.update_ttv_file_buttons()

         #   self.update_RV_plots()
        #    self.update_extra_plots()
        #    self.update_transit_plots()
            self.update_ttv_plots()

            #self.update_activity_data_plots()
            #self.update_activity_gls_plots()
        else:
            return

    ### AST ####

    def update_color_picker_ast(self):
        global fit

        font = QtGui.QFont()
        font.setPointSize(11)
        font.setBold(False)
        #font.setWeight(75)

        for i in range(11):
            self.buttonGroup_color_picker_ast.button(i+1).setStyleSheet("color: %s;"%dill.copy(fit.ast_colors[i]))
            self.buttonGroup_color_picker_ast.button(i+1).setFont(font)
        for i in range(10):
            self.buttonGroup_symbol_picker_ast.button(i+1).setStyleSheet("color: %s;"%dill.copy(fit.ast_colors[i]))
            self.buttonGroup_symbol_picker_ast.button(i+1).setText(dill.copy(dill.copy(fit.pyqt_symbols_ast[i])))
            self.buttonGroup_symbol_picker_ast.button(i+1).setFont(font)

    def get_color_ast(self):
        global fit

        but_ind = self.buttonGroup_color_picker_ast.checkedId()
        colorz = self.colorDialog.getColor(options=QtWidgets.QColorDialog.DontUseNativeDialog)

        #QtWidgets.QColorDialog.setOption(QtWidgets.QColorDialog.ShowAlphaChannel,True)
        #colorz = QtWidgets.QColorDialog.getColor()

        if colorz.isValid():
            fit.ast_colors[but_ind-1]=colorz.name()
            self.update_color_picker_ast()
           # self.update_act_file_buttons()
          #  self.update_RV_file_buttons()
         #   self.update_tra_file_buttons()
            self.update_ast_file_buttons()
            self.update_ast_file_buttons_2()
         #   self.update_RV_plots()
        #    self.update_extra_plots()
        #    self.update_transit_plots()
            self.update_ast_plots()

            #self.update_activity_data_plots()
            #self.update_activity_gls_plots()
        else:
            return


############################# Symbol controls ################################  

    ### RV ####

    def get_symbol(self):
        global fit
 
        but_ind = self.buttonGroup_symbol_picker.checkedId()   
        but_n = self.dialog_symbols.get_radio()
            
        if but_n != None:
            fit.pyqt_symbols_rvs[but_ind-1] = dill.copy(symbols[but_n-1])
            self.update_color_picker()
            self.update_act_file_buttons()      
            self.update_RV_file_buttons() 
            self.update_tra_file_buttons() 
            self.update_RV_plots() 
            self.update_extra_plots()
            self.update_extra_plots2()
            self.update_transit_plots()     
        else:
            return  



        
      
    def check_RV_symbol_sizes(self):
        global fit
       
       # for i in range(10):
        fit.pyqt_symbols_size_rvs[0] = dill.copy(self.rv_data_size_1.value())
        fit.pyqt_symbols_size_rvs[1] = dill.copy(self.rv_data_size_2.value())
        fit.pyqt_symbols_size_rvs[2] = dill.copy(self.rv_data_size_3.value())
        fit.pyqt_symbols_size_rvs[3] = dill.copy(self.rv_data_size_4.value())
        fit.pyqt_symbols_size_rvs[4] = dill.copy(self.rv_data_size_5.value())
        fit.pyqt_symbols_size_rvs[5] = dill.copy(self.rv_data_size_6.value())
        fit.pyqt_symbols_size_rvs[6] = dill.copy(self.rv_data_size_7.value())
        fit.pyqt_symbols_size_rvs[7] = dill.copy(self.rv_data_size_8.value())
        fit.pyqt_symbols_size_rvs[8] = dill.copy(self.rv_data_size_9.value())
        fit.pyqt_symbols_size_rvs[9] = dill.copy(self.rv_data_size_10.value())
        fit.pyqt_symbols_size_rvs[10] = dill.copy(self.rv_data_size_11.value())
        fit.pyqt_symbols_size_rvs[11] = dill.copy(self.rv_data_size_12.value())
        fit.pyqt_symbols_size_rvs[12] = dill.copy(self.rv_data_size_13.value())
        fit.pyqt_symbols_size_rvs[13] = dill.copy(self.rv_data_size_14.value())
        fit.pyqt_symbols_size_rvs[14] = dill.copy(self.rv_data_size_15.value())
        fit.pyqt_symbols_size_rvs[15] = dill.copy(self.rv_data_size_16.value())
        fit.pyqt_symbols_size_rvs[16] = dill.copy(self.rv_data_size_17.value())
        fit.pyqt_symbols_size_rvs[17] = dill.copy(self.rv_data_size_18.value())
        fit.pyqt_symbols_size_rvs[18] = dill.copy(self.rv_data_size_19.value())
        fit.pyqt_symbols_size_rvs[19] = dill.copy(self.rv_data_size_20.value())

        fit.pyqt_color_alpha_rvs[0] = dill.copy(self.rv_data_alpha_1.value())
        fit.pyqt_color_alpha_rvs[1] = dill.copy(self.rv_data_alpha_2.value())
        fit.pyqt_color_alpha_rvs[2] = dill.copy(self.rv_data_alpha_3.value())
        fit.pyqt_color_alpha_rvs[3] = dill.copy(self.rv_data_alpha_4.value())
        fit.pyqt_color_alpha_rvs[4] = dill.copy(self.rv_data_alpha_5.value())
        fit.pyqt_color_alpha_rvs[5] = dill.copy(self.rv_data_alpha_6.value())
        fit.pyqt_color_alpha_rvs[6] = dill.copy(self.rv_data_alpha_7.value())
        fit.pyqt_color_alpha_rvs[7] = dill.copy(self.rv_data_alpha_8.value())
        fit.pyqt_color_alpha_rvs[8] = dill.copy(self.rv_data_alpha_9.value())
        fit.pyqt_color_alpha_rvs[9] = dill.copy(self.rv_data_alpha_10.value())
        fit.pyqt_color_alpha_rvs[10] = dill.copy(self.rv_data_alpha_11.value())
        fit.pyqt_color_alpha_rvs[11] = dill.copy(self.rv_data_alpha_12.value())
        fit.pyqt_color_alpha_rvs[12] = dill.copy(self.rv_data_alpha_13.value())
        fit.pyqt_color_alpha_rvs[13] = dill.copy(self.rv_data_alpha_14.value())
        fit.pyqt_color_alpha_rvs[14] = dill.copy(self.rv_data_alpha_15.value())
        fit.pyqt_color_alpha_rvs[15] = dill.copy(self.rv_data_alpha_16.value())
        fit.pyqt_color_alpha_rvs[16] = dill.copy(self.rv_data_alpha_17.value())
        fit.pyqt_color_alpha_rvs[17] = dill.copy(self.rv_data_alpha_18.value())
        fit.pyqt_color_alpha_rvs[18] = dill.copy(self.rv_data_alpha_19.value())
        fit.pyqt_color_alpha_rvs[19] = dill.copy(self.rv_data_alpha_20.value())

    ### Transit ###

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
            self.update_extra_plots2()
            self.update_transit_plots()     
        else:
            return    
        
        

        
        
    def check_tra_symbol_sizes(self):
        global fit
       
       # for i in range(10):
        fit.pyqt_symbols_size_tra[0] = dill.copy(self.trans_data_size_1.value())
        fit.pyqt_symbols_size_tra[1] = dill.copy(self.trans_data_size_2.value())
        fit.pyqt_symbols_size_tra[2] = dill.copy(self.trans_data_size_3.value())
        fit.pyqt_symbols_size_tra[3] = dill.copy(self.trans_data_size_4.value())
        fit.pyqt_symbols_size_tra[4] = dill.copy(self.trans_data_size_5.value())
        fit.pyqt_symbols_size_tra[5] = dill.copy(self.trans_data_size_6.value())
        fit.pyqt_symbols_size_tra[6] = dill.copy(self.trans_data_size_7.value())
        fit.pyqt_symbols_size_tra[7] = dill.copy(self.trans_data_size_8.value())
        fit.pyqt_symbols_size_tra[8] = dill.copy(self.trans_data_size_9.value())
        fit.pyqt_symbols_size_tra[9] = dill.copy(self.trans_data_size_10.value())
        fit.pyqt_symbols_size_tra[10] = dill.copy(self.trans_data_size_11.value())
        fit.pyqt_symbols_size_tra[11] = dill.copy(self.trans_data_size_12.value())
        fit.pyqt_symbols_size_tra[12] = dill.copy(self.trans_data_size_13.value())
        fit.pyqt_symbols_size_tra[13] = dill.copy(self.trans_data_size_14.value())
        fit.pyqt_symbols_size_tra[14] = dill.copy(self.trans_data_size_15.value())
        fit.pyqt_symbols_size_tra[15] = dill.copy(self.trans_data_size_16.value())
        fit.pyqt_symbols_size_tra[16] = dill.copy(self.trans_data_size_17.value())
        fit.pyqt_symbols_size_tra[17] = dill.copy(self.trans_data_size_18.value())
        fit.pyqt_symbols_size_tra[18] = dill.copy(self.trans_data_size_19.value())
        fit.pyqt_symbols_size_tra[19] = dill.copy(self.trans_data_size_20.value())

        fit.pyqt_color_alpha_tra[0] = dill.copy(self.trans_data_alpha_1.value())
        fit.pyqt_color_alpha_tra[1] = dill.copy(self.trans_data_alpha_2.value())
        fit.pyqt_color_alpha_tra[2] = dill.copy(self.trans_data_alpha_3.value())
        fit.pyqt_color_alpha_tra[3] = dill.copy(self.trans_data_alpha_4.value())
        fit.pyqt_color_alpha_tra[4] = dill.copy(self.trans_data_alpha_5.value())
        fit.pyqt_color_alpha_tra[5] = dill.copy(self.trans_data_alpha_6.value())
        fit.pyqt_color_alpha_tra[6] = dill.copy(self.trans_data_alpha_7.value())
        fit.pyqt_color_alpha_tra[7] = dill.copy(self.trans_data_alpha_8.value())
        fit.pyqt_color_alpha_tra[8] = dill.copy(self.trans_data_alpha_9.value())
        fit.pyqt_color_alpha_tra[9] = dill.copy(self.trans_data_alpha_10.value())
        fit.pyqt_color_alpha_tra[10] = dill.copy(self.trans_data_alpha_11.value())
        fit.pyqt_color_alpha_tra[11] = dill.copy(self.trans_data_alpha_12.value())
        fit.pyqt_color_alpha_tra[12] = dill.copy(self.trans_data_alpha_13.value())
        fit.pyqt_color_alpha_tra[13] = dill.copy(self.trans_data_alpha_14.value())
        fit.pyqt_color_alpha_tra[14] = dill.copy(self.trans_data_alpha_15.value())
        fit.pyqt_color_alpha_tra[15] = dill.copy(self.trans_data_alpha_16.value())
        fit.pyqt_color_alpha_tra[16] = dill.copy(self.trans_data_alpha_17.value())
        fit.pyqt_color_alpha_tra[17] = dill.copy(self.trans_data_alpha_18.value())
        fit.pyqt_color_alpha_tra[18] = dill.copy(self.trans_data_alpha_19.value())
        fit.pyqt_color_alpha_tra[19] = dill.copy(self.trans_data_alpha_20.value())


    ### TTV ####

    def get_symbol_ttv(self):
        global fit
 
        but_ind = self.buttonGroup_symbol_picker_ttv.checkedId()   
        but_n = self.dialog_symbols.get_radio()
            
        if but_n != None:
            fit.pyqt_symbols_ttv[but_ind-1] = dill.copy(symbols[but_n-1])
            self.update_color_picker_ttv()
        #    self.update_act_file_buttons()      
        #    self.update_RV_file_buttons() 
        #    self.update_tra_file_buttons()
            self.update_ttv_file_buttons() 

        #    self.update_RV_plots() 
        #    self.update_extra_plots()
        #    self.update_transit_plots()
            self.update_ttv_plots()

        else:
            return    
        
    def check_ttv_symbol_sizes(self):
        global fit
       
       # for i in range(10):
        fit.pyqt_symbols_size_ttv[0] = dill.copy(self.ttv_data_size_1.value())
        fit.pyqt_symbols_size_ttv[1] = dill.copy(self.ttv_data_size_2.value())
        fit.pyqt_symbols_size_ttv[2] = dill.copy(self.ttv_data_size_3.value())
        fit.pyqt_symbols_size_ttv[3] = dill.copy(self.ttv_data_size_4.value())
        fit.pyqt_symbols_size_ttv[4] = dill.copy(self.ttv_data_size_5.value())
        fit.pyqt_symbols_size_ttv[5] = dill.copy(self.ttv_data_size_6.value())
        fit.pyqt_symbols_size_ttv[6] = dill.copy(self.ttv_data_size_7.value())
        fit.pyqt_symbols_size_ttv[7] = dill.copy(self.ttv_data_size_8.value())
        fit.pyqt_symbols_size_ttv[8] = dill.copy(self.ttv_data_size_9.value())
        fit.pyqt_symbols_size_ttv[9] = dill.copy(self.ttv_data_size_10.value())

 
    ### AST ####

    def get_symbol_ast(self):
        global fit
 
        but_ind = self.buttonGroup_symbol_picker_ast.checkedId()   
        but_n = self.dialog_symbols.get_radio()
            
        if but_n != None:
            fit.pyqt_symbols_ast[but_ind-1] = dill.copy(symbols[but_n-1])
            self.update_color_picker_ast()
        #    self.update_act_file_buttons()      
        #    self.update_RV_file_buttons() 
        #    self.update_tra_file_buttons()
            self.update_ast_file_buttons() 
            self.update_ast_file_buttons_2()
        #    self.update_RV_plots() 
        #    self.update_extra_plots()
        #    self.update_transit_plots()
            self.update_ast_plots()

        else:
            return    
        
    def check_ast_symbol_sizes(self):
        global fit
       
       # for i in range(10):
        fit.pyqt_symbols_size_ast[0] = dill.copy(self.ast_data_size_1.value())
        fit.pyqt_symbols_size_ast[1] = dill.copy(self.ast_data_size_2.value())
        fit.pyqt_symbols_size_ast[2] = dill.copy(self.ast_data_size_3.value())
        fit.pyqt_symbols_size_ast[3] = dill.copy(self.ast_data_size_4.value())
        fit.pyqt_symbols_size_ast[4] = dill.copy(self.ast_data_size_5.value())
        fit.pyqt_symbols_size_ast[5] = dill.copy(self.ast_data_size_6.value())
        fit.pyqt_symbols_size_ast[6] = dill.copy(self.ast_data_size_7.value())
        fit.pyqt_symbols_size_ast[7] = dill.copy(self.ast_data_size_8.value())
        fit.pyqt_symbols_size_ast[8] = dill.copy(self.ast_data_size_9.value())
        fit.pyqt_symbols_size_ast[9] = dill.copy(self.ast_data_size_10.value())
#        fit.pyqt_symbols_size_ast[10] = dill.copy(self.ast_data_size_11.value())
#        fit.pyqt_symbols_size_ast[11] = dill.copy(self.ast_data_size_12.value())
#        fit.pyqt_symbols_size_ast[12] = dill.copy(self.ast_data_size_13.value())
#        fit.pyqt_symbols_size_ast[13] = dill.copy(self.ast_data_size_14.value())
#        fit.pyqt_symbols_size_ast[14] = dill.copy(self.ast_data_size_15.value())
#        fit.pyqt_symbols_size_ast[15] = dill.copy(self.ast_data_size_16.value())
#        fit.pyqt_symbols_size_ast[16] = dill.copy(self.ast_data_size_17.value())
#        fit.pyqt_symbols_size_ast[17] = dill.copy(self.ast_data_size_18.value())
#        fit.pyqt_symbols_size_ast[18] = dill.copy(self.ast_data_size_19.value())
#        fit.pyqt_symbols_size_ast[19] = dill.copy(self.ast_data_size_20.value())
 
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
       # filename = QtWidgets.QFileDialog.getSaveFileName(self, 'Save image', filter="PNG(*.png);; JPEG(*.jpg)")
        #p.save(filename[0], 'jpg')        
        #label.setPixmap(p)        # just for fun :)
        img, _ = QtWidgets.QFileDialog.getSaveFileName(self,"Save image",
                                            filter="PNG(*.png);; JPEG(*.jpg)", options=QtWidgets.QFileDialog.Option.DontUseNativeDialog)
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

        self.console_widget.print_text("rv.latex_pl_param_table(fit, width = 10, precision = 2, asymmetric = False, file_name='best_fit_param_table.tex',path='./')", before_prompt=False)
#        self.console_widget.execute_command("rv.latex_pl_param_table(fit, width = 10, precision = 2, asymmetric = False, file_name='test.tex',path='./')")  
        self.tabWidget_helper.setCurrentWidget(self.tab_shells) 
        self.terminal_embeded.setCurrentWidget(self.console_widget)

    def get_latex_prior_table(self):
        global fit

        self.console_widget.print_text("rv.latex_prior_table(fit, width = 10, precision = 2,  file_name='prior_table.tex',path='./')", before_prompt=False)  
        self.tabWidget_helper.setCurrentWidget(self.tab_shells) 
        self.terminal_embeded.setCurrentWidget(self.console_widget)

    def get_RV_model(self):
        global fit

        self.console_widget.print_text("rv.export_RV_model(fit, file='RV_model.txt', width = 10, precision = 4)", before_prompt=False)  
        self.tabWidget_helper.setCurrentWidget(self.tab_shells) 
        self.terminal_embeded.setCurrentWidget(self.console_widget)
        
    def get_RV_data(self):
        global fit

        #self.console_widget.print_text("rv.export_RV_data(fit, [0], file='RV_data.txt',  jitter=False, o_c=False, print_data=False, remove_offset = False, width = 10, precision = 3)", before_prompt=False)  
        self.console_widget.print_text("rv.export_RV_data(fit, file='RV_data.txt', delimiter=' ',  print_data=False,  header = True, width = 10, precision = 3)", before_prompt=False)          
        
        self.tabWidget_helper.setCurrentWidget(self.tab_shells) 
        self.terminal_embeded.setCurrentWidget(self.console_widget)   
        
        
    def get_orb_evol(self):
        global fit       

        self.console_widget.print_text("rv.export_orbital_evol(fit, file='planet_N.txt', planet = 1, width = 10, precision = 6)", before_prompt=False)  
        self.tabWidget_helper.setCurrentWidget(self.tab_shells) 
        self.terminal_embeded.setCurrentWidget(self.console_widget)
 
 

################################## Stellar params #######################################

    def update_St_params(self, ind = None):
        global fit

        if ind ==1: fit.stellar_mass     = self.St_mass_input.value()
        if ind ==2: fit.stellar_mass_err = self.err_St_mass_input.value()
        
        if ind ==3: fit.stellar_radius     = self.St_radius_input.value()
        if ind ==4: fit.stellar_radius_err = self.err_St_radius_input.value()
        
        if ind ==5: fit.stellar_luminosity     = self.St_lumin_input.value()
        if ind ==6: fit.stellar_luminosity_err = self.err_St_lumin_input.value()
        
        if ind ==7: fit.stellar_Teff       = self.St_teff_input.value()
        if ind ==8: fit.stellar_Teff_err   = self.err_St_teff_input.value()
 
        if ind ==9: fit.stellar_vsini       = self.St_vsini_input.value()
        if ind ==10: fit.stellar_vsini_err   = self.err_St_vsini_input.value()
        
        st_rot = rv.get_stellar_rotation(fit)
        kb1995 = rv.get_rv_scatter(fit)
        kb2011 = rv.get_rv_scatter(fit,use_kb2011=True)
        
        self.label_St_rot_value.setText("%.4f +/- %.4f [days]"%(st_rot[0],st_rot[1]))
        self.label_kb1995.setText("%.4f +/- %.4f [m/sec]"%(kb1995[0],kb1995[1]))
        self.label_kb2011.setText("%.4f +/- %.4f [m/sec]"%(kb2011[0],kb2011[1]))


    def update_GUI_St_params(self):
        global fit

        self.St_mass_input.setValue(fit.stellar_mass)
        self.err_St_mass_input.setValue(fit.stellar_mass_err)
        
        self.St_radius_input.setValue(fit.stellar_radius)
        self.err_St_radius_input.setValue(fit.stellar_radius_err)
        
        self.St_lumin_input.setValue(fit.stellar_luminosity)
        self.err_St_lumin_input.setValue(fit.stellar_luminosity_err)
        
        self.St_teff_input.setValue(fit.stellar_Teff)
        self.err_St_teff_input.setValue(fit.stellar_Teff_err)
        
        self.St_vsini_input.setValue(fit.stellar_vsini)
        self.err_St_vsini_input.setValue(fit.stellar_vsini_err)
        
 

################################## System #######################################
 
    def set_Win_widget_Style(self, widget):
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create('Windows'))
    def set_Fus_widget_Style(self, widget):
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create('Fusion'))
    def set_Mac_widget_Style(self, widget):
        if sys.platform != "darwin":
            self.tabWidget_helper.setCurrentWidget(self.tab_info)
            print("\n 'Macintosh' window style is only available on MAC OS !!!\n")
            return
        else:
            QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create('Macintosh'))


    def set_widget_font(self, widget):
        #QtWidgets.QFontDialog.setOption(QtWidgets.QFontDialog.DontUseNativeDialog, True)
        font, ok = QtWidgets.QFontDialog.getFont(self.gui_font)

        if ok:
            self.gui_font = font

            QtWidgets.QApplication.setFont(self.gui_font)

            for topLevel in QtWidgets.QApplication.allWidgets():
                topLevel.setFont(self.gui_font)

     ###### Update settings json ########
        GUI_settings = './lib/ES_settings.json'

        with open(GUI_settings, "r") as jsonFile:
            sett = json.load(jsonFile)

        sett["GUI_fonts"]["gui_font"] = str(font.pointSize())

        with open(GUI_settings, "w") as jsonFile:
            json.dump(sett, jsonFile)


 
    def set_plot_font(self):
        #QtWidgets.QFontDialog.setOption(QtWidgets.QFontDialog.DontUseNativeDialog, True)
        font, ok = QtWidgets.QFontDialog.getFont(self.plot_font)

        if ok:
            self.plot_font.setFamily(font.family())
            self.plot_font.setPointSize(font.pointSize())
          
        self.update_font_plots()   


     ###### Update settings json ########
        Plot_settings = './lib/ES_settings.json'

        with open(Plot_settings, "r") as jsonFile:
            sett = json.load(jsonFile)

        sett["fonts"]["plot_font"] = str(font.pointSize())

        with open(Plot_settings, "w") as jsonFile:
            json.dump(sett, jsonFile)

         
    def initialize_font_plot(self): #not working as I want!

        Plot_settings = './lib/ES_settings.json'
        with open(Plot_settings, "r") as jsonFile:
            sett = json.load(jsonFile)

        self.plot_font = QtGui.QFont()
        self.plot_font.setPointSize(int(sett["fonts"]["plot_font"]))
        self.plot_font.setBold(False)
              

    def initialize_font_GUI(self): #not working as I want!

        Plot_settings = './lib/ES_settings.json'
        with open(Plot_settings, "r") as jsonFile:
            sett = json.load(jsonFile)

        self.gui_font = QtGui.QFont()
        self.gui_font.setPointSize(int(sett["GUI_fonts"]["gui_font"]))
        self.gui_font.setBold(False)
              
        for topLevel in QtWidgets.QApplication.allWidgets():
            topLevel.setFont(self.gui_font)


    def closeEvent(self, event):
        choice = QtWidgets.QMessageBox.information(self, 'Warning!',
                                            "Do you want to save the session before you Quit?",
                                            QtWidgets.QMessageBox.StandardButton.Cancel | QtWidgets.QMessageBox.StandardButton.No | QtWidgets.QMessageBox.StandardButton.Yes)

        if choice == QtWidgets.QMessageBox.StandardButton.No:

            if sys.platform[0:5] == "linux":
                self.term_emb.close()
            self.removeEventFilter(self)
            event.accept()

        elif choice == QtWidgets.QMessageBox.StandardButton.Yes:
            self.save_session()
            if sys.platform[0:5] == "linux":
                self.term_emb.close()
            #QtWidgets.QApplication.instance().removeEventFilter(self) 
            self.removeEventFilter(self)
            event.accept()

        elif choice == QtWidgets.QMessageBox.StandardButton.Cancel:
            event.ignore()


    def reset_mid_pannel_buttons(self):
 
        self.button_orb_evol.setEnabled(True)        
        self.mute_buttons(trigger=True)



    def mute_buttons(self, trigger=True):
        global fit
        
        
        self.safe_to_init = trigger
        
        self.button_init_fit.setEnabled(trigger)                
        self.button_MCMC.setEnabled(trigger) 
        self.button_nest_samp.setEnabled(trigger) 
        self.button_auto_fit.setEnabled(trigger) 
        self.button_fit.setEnabled(trigger) 
        self.radioButton_RV.setEnabled(trigger) 
        self.radioButton_transit.setEnabled(trigger) 
        self.radioButton_ttv.setEnabled(trigger) 
        self.radioButton_transit_RV.setEnabled(trigger) 
        self.radioButton_ast.setEnabled(trigger) 
        self.radioButton_ast_RV.setEnabled(trigger) 
        self.radioButton_ttv_RV.setEnabled(trigger) 
        self.comboBox_select_ses.setEnabled(trigger) 
        self.new_ses.setEnabled(trigger) 
        self.copy_ses.setEnabled(trigger) 
        self.remove_ses.setEnabled(trigger) 
        self.radioButton_ewm.setEnabled(trigger)       
        self.radioButton_hkl.setEnabled(trigger) 
        self.radioButton_rho.setEnabled(trigger)       
        self.radioButton_aR.setEnabled(trigger) 

        
        self.radioButton_Dynamical.setEnabled(trigger) 
        self.radioButton_Keplerian.setEnabled(trigger) 
       
        for i in range(9):
            self.planet_checked_gui[i].setEnabled(trigger)
        

        self.Button_RV_data_1.setEnabled(trigger) 
        self.Button_RV_data_2.setEnabled(trigger) 
        self.Button_RV_data_3.setEnabled(trigger) 
        self.Button_RV_data_4.setEnabled(trigger) 
        self.Button_RV_data_5.setEnabled(trigger) 
        self.Button_RV_data_6.setEnabled(trigger) 
        self.Button_RV_data_7.setEnabled(trigger) 
        self.Button_RV_data_8.setEnabled(trigger) 
        self.Button_RV_data_9.setEnabled(trigger) 
        self.Button_RV_data_10.setEnabled(trigger) 
        self.Button_RV_data_11.setEnabled(trigger) 
        self.Button_RV_data_12.setEnabled(trigger) 
        self.Button_RV_data_13.setEnabled(trigger) 
        self.Button_RV_data_14.setEnabled(trigger) 
        self.Button_RV_data_15.setEnabled(trigger) 
        self.Button_RV_data_16.setEnabled(trigger) 
        self.Button_RV_data_17.setEnabled(trigger) 
        self.Button_RV_data_18.setEnabled(trigger) 
        self.Button_RV_data_19.setEnabled(trigger) 
        self.Button_RV_data_20.setEnabled(trigger) 
        
        self.remove_rv_data1.setEnabled(trigger) 
        self.remove_rv_data2.setEnabled(trigger) 
        self.remove_rv_data3.setEnabled(trigger) 
        self.remove_rv_data4.setEnabled(trigger) 
        self.remove_rv_data5.setEnabled(trigger) 
        self.remove_rv_data6.setEnabled(trigger) 
        self.remove_rv_data7.setEnabled(trigger) 
        self.remove_rv_data8.setEnabled(trigger) 
        self.remove_rv_data9.setEnabled(trigger) 
        self.remove_rv_data10.setEnabled(trigger) 
        self.remove_rv_data11.setEnabled(trigger) 
        self.remove_rv_data12.setEnabled(trigger) 
        self.remove_rv_data13.setEnabled(trigger) 
        self.remove_rv_data14.setEnabled(trigger) 
        self.remove_rv_data15.setEnabled(trigger) 
        self.remove_rv_data16.setEnabled(trigger) 
        self.remove_rv_data17.setEnabled(trigger) 
        self.remove_rv_data18.setEnabled(trigger) 
        self.remove_rv_data19.setEnabled(trigger) 
        self.remove_rv_data20.setEnabled(trigger) 

        self.Button_transit_data_1.setEnabled(trigger) 
        self.Button_transit_data_2.setEnabled(trigger) 
        self.Button_transit_data_3.setEnabled(trigger) 
        self.Button_transit_data_4.setEnabled(trigger) 
        self.Button_transit_data_5.setEnabled(trigger) 
        self.Button_transit_data_6.setEnabled(trigger) 
        self.Button_transit_data_7.setEnabled(trigger) 
        self.Button_transit_data_8.setEnabled(trigger) 
        self.Button_transit_data_9.setEnabled(trigger) 
        self.Button_transit_data_10.setEnabled(trigger)
        self.Button_transit_data_11.setEnabled(trigger) 
        self.Button_transit_data_12.setEnabled(trigger) 
        self.Button_transit_data_13.setEnabled(trigger) 
        self.Button_transit_data_14.setEnabled(trigger) 
        self.Button_transit_data_15.setEnabled(trigger) 
        self.Button_transit_data_16.setEnabled(trigger) 
        self.Button_transit_data_17.setEnabled(trigger) 
        self.Button_transit_data_18.setEnabled(trigger) 
        self.Button_transit_data_19.setEnabled(trigger) 
        self.Button_transit_data_20.setEnabled(trigger)

        
        self.remove_transit_data_1.setEnabled(trigger) 
        self.remove_transit_data_2.setEnabled(trigger) 
        self.remove_transit_data_3.setEnabled(trigger) 
        self.remove_transit_data_4.setEnabled(trigger) 
        self.remove_transit_data_5.setEnabled(trigger) 
        self.remove_transit_data_6.setEnabled(trigger) 
        self.remove_transit_data_7.setEnabled(trigger) 
        self.remove_transit_data_8.setEnabled(trigger) 
        self.remove_transit_data_9.setEnabled(trigger) 
        self.remove_transit_data_10.setEnabled(trigger)  
        self.remove_transit_data_11.setEnabled(trigger) 
        self.remove_transit_data_12.setEnabled(trigger) 
        self.remove_transit_data_13.setEnabled(trigger) 
        self.remove_transit_data_14.setEnabled(trigger) 
        self.remove_transit_data_15.setEnabled(trigger) 
        self.remove_transit_data_16.setEnabled(trigger) 
        self.remove_transit_data_17.setEnabled(trigger) 
        self.remove_transit_data_18.setEnabled(trigger) 
        self.remove_transit_data_19.setEnabled(trigger) 
        self.remove_transit_data_20.setEnabled(trigger)  
                
        
        self.Button_ttv_data_1.setEnabled(trigger) 
        self.Button_ttv_data_2.setEnabled(trigger) 
        self.Button_ttv_data_3.setEnabled(trigger) 
        self.Button_ttv_data_4.setEnabled(trigger) 
        self.Button_ttv_data_5.setEnabled(trigger) 
        self.Button_ttv_data_6.setEnabled(trigger) 
        self.Button_ttv_data_7.setEnabled(trigger) 
        self.Button_ttv_data_8.setEnabled(trigger) 
        self.Button_ttv_data_9.setEnabled(trigger) 
        self.Button_ttv_data_10.setEnabled(trigger) 
        
        self.remove_ttv_data_1.setEnabled(trigger) 
        self.remove_ttv_data_2.setEnabled(trigger) 
        self.remove_ttv_data_3.setEnabled(trigger) 
        self.remove_ttv_data_4.setEnabled(trigger) 
        self.remove_ttv_data_5.setEnabled(trigger) 
        self.remove_ttv_data_6.setEnabled(trigger) 
        self.remove_ttv_data_7.setEnabled(trigger) 
        self.remove_ttv_data_8.setEnabled(trigger) 
        self.remove_ttv_data_9.setEnabled(trigger) 
        self.remove_ttv_data_10.setEnabled(trigger)        

        self.Button_ast_data_1.setEnabled(trigger) 
        self.Button_ast_data_2.setEnabled(trigger) 
        self.Button_ast_data_3.setEnabled(trigger) 
        self.Button_ast_data_4.setEnabled(trigger) 
        self.Button_ast_data_5.setEnabled(trigger) 
        self.Button_ast_data_6.setEnabled(trigger) 
        self.Button_ast_data_7.setEnabled(trigger) 
        self.Button_ast_data_8.setEnabled(trigger) 
        self.Button_ast_data_9.setEnabled(trigger) 
        self.Button_ast_data_10.setEnabled(trigger) 

        
        self.remove_ast_data_1.setEnabled(trigger) 
        self.remove_ast_data_2.setEnabled(trigger) 
        self.remove_ast_data_3.setEnabled(trigger) 
        self.remove_ast_data_4.setEnabled(trigger) 
        self.remove_ast_data_5.setEnabled(trigger) 
        self.remove_ast_data_6.setEnabled(trigger) 
        self.remove_ast_data_7.setEnabled(trigger) 
        self.remove_ast_data_8.setEnabled(trigger) 
        self.remove_ast_data_9.setEnabled(trigger) 
        self.remove_ast_data_10.setEnabled(trigger)               
        #self.use_Planet1.setEnabled(trigger) 
        
        self.Button_ast_data_hipp_1.setEnabled(trigger) 
        self.Button_ast_data_hipp_2.setEnabled(trigger) 
        self.Button_ast_data_gaia_1.setEnabled(trigger)       
        self.Button_ast_data_gaia_2.setEnabled(trigger)      
        
        self.remove_ast_data_hipp_1.setEnabled(trigger) 
        self.remove_ast_data_hipp_2.setEnabled(trigger)        
        self.remove_ast_data_gaia_1.setEnabled(trigger) 
        self.remove_ast_data_gaia_2.setEnabled(trigger)              
        
    def file_from_path(self, path):
        head, tail = ntpath.split(path)
        return tail or ntpath.basename(head)
        
       
    def update_settings(self):
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
        self.force_copl_incl.setChecked(fit.copl_incl)
        self.force_sameRV_jit.setChecked(fit.jit_flag)       
        
    def check_settings(self):
        
        fit.time_step_model = self.time_step_model.value()
        fit.dyn_model_accuracy = self.dyn_model_accuracy.value()
        fit.dyn_model_to_kill = self.dyn_model_to_kill.value()
        fit.kep_model_to_kill = self.kep_model_to_kill.value()
        fit.master_timeout = self.master_timeout.value()    
                
        
    def adopt_RV_GLS_param(self):
        global fit   
        
        mean_anomaly_from_gls = np.degrees((((fit.epoch - float(fit.gls_o_c.hpstat["T0"]) )% (fit.gls_o_c.hpstat["P"]) )/ (fit.gls_o_c.hpstat["P"]) ) * 2*np.pi)
 
        fit.add_planet(fit.gls_o_c.hpstat["amp"],fit.gls_o_c.hpstat["P"],0.0,0.0,mean_anomaly_from_gls -90.0,90.0,0.0)
#        fit.use.update_use_planet_params_one_planet(fit.npl-1,True,True,fit.auto_fit_allow_ecc,fit.auto_fit_allow_ecc,True,False,False)   

        fit.P_use[fit.npl] = True
        fit.K_use[fit.npl] = True
        fit.e_use[fit.npl] = fit.auto_fit_allow_ecc
        fit.w_use[fit.npl] = fit.auto_fit_allow_ecc 
        fit.M0_use[fit.npl] = True 
      
        self.update_use_from_session()   
        self.update_use() 
        #self.update_params()  
       # self.update_gui_params()
        self.optimize_fit(0,m_ln=self.amoeba_radio_button.isChecked(),auto_fit = True)
    
 
    def adopt_RV_mlp_param(self):
        global fit   
        
        mean_anomaly_from_gls = np.degrees((((fit.epoch - float(fit.mlp.hpstat["T0"]) )% (fit.mlp.hpstat["P"]) )/ (fit.mlp.hpstat["P"]) ) * 2*np.pi)
 
        fit.add_planet(fit.mlp.hpstat["amp"],fit.mlp.hpstat["P"],0.0,0.0,mean_anomaly_from_gls -90.0,90.0,0.0)
        #fit.use.update_use_planet_params_one_planet(fit.npl-1,True,True,fit.auto_fit_allow_ecc,fit.auto_fit_allow_ecc,True,False,False)   
  

        fit.P_use[fit.npl] = True
        fit.K_use[fit.npl] = True
        fit.e_use[fit.npl] = fit.auto_fit_allow_ecc
        fit.w_use[fit.npl] = fit.auto_fit_allow_ecc 
        fit.M0_use[fit.npl] = True 
     
        self.update_use_from_session()   
        self.update_use() 
        #self.update_params()  
       # self.update_gui_params()
        self.optimize_fit(0,m_ln=self.amoeba_radio_button.isChecked(),auto_fit = True)
    


    def adopt_trans_TLS_param(self):
        global fit   
 
        
        if len(fit.tls_o_c)==0 or np.isnan(fit.tls_o_c.period):
            print("No transit signal that can be adopted as planet")
            return
    
        mean_anomaly_from_tls = np.degrees((((fit.epoch - fit.tls_o_c.transit_times[0] )% (fit.tls_o_c.period) )/ (fit.tls_o_c.period) ) * 2*np.pi)
       
        fit.t0[fit.npl]     = fit.tls_o_c.transit_times[0]
        fit.pl_rad[fit.npl] = fit.stellar_radius*np.sqrt(1.0 - fit.tls_o_c.depth_mean[0])  # alternativly fit.tls_o_c.rp_rs ?
        fit.pl_a[fit.npl]   = 11.44

        fit.add_planet(10.0,fit.tls_o_c.period,0.0,0.0,mean_anomaly_from_tls,90.0,0.0)
        #fit.use.update_use_planet_params_one_planet(fit.npl+1,True,True,True,True,True,False,False)   

        fit.P_use[fit.npl] = True
        fit.K_use[fit.npl] = True
        fit.e_use[fit.npl] = fit.auto_fit_allow_ecc
        fit.w_use[fit.npl] = fit.auto_fit_allow_ecc 
        fit.M0_use[fit.npl] = True 

        self.update_use_from_session()
        self.update_use() 
       # self.update_params()  
        self.update_gui_params()
        self.radioButton_transit.setChecked(True)

        self.worker_transit_fitting(ff=0 )


 
#############################  TEST ZONE ################################  

    def check_use_tra_GP(self):
        global fit
        #but_ind = self.buttonGroup_use_tra_data_GP.checkedId()   

        for j in range(20):
    
            if len(fit.tra_data_sets[j]) == 0:
                continue
            else:
                self.use_tra_data_GP[j].setChecked(bool(fit.tra_data_sets[j][9]))
        #print("Test",but_ind)       
        return

    def set_use_astrometry(self):
        global fit

        fit.use_ast_hipp_gaia = self.use_hipp_gaia_ast.isChecked()
        fit.use_ast_imaging   = self.use_imaging_ast.isChecked() 

    def get_use_astrometry(self):
        global fit

        self.use_hipp_gaia_ast.setChecked(bool(fit.use_ast_hipp_gaia))
        self.use_imaging_ast.setChecked(bool(fit.use_ast_imaging))

    def set_use_tra_GP_data(self):
        global fit
        #but_ind = self.buttonGroup_use_tra_data_GP.checkedId()   

        for j in range(20):
    
            if len(fit.tra_data_sets[j]) == 0:
                continue
            else:
                fit.tra_data_sets[j][9] = self.use_tra_data_GP[j].isChecked()
                
                
        if len([fit.tra_data_sets[j][9] for j in range(10) if len(fit.tra_data_sets[j]) != 0  and fit.tra_data_sets[j][9] ==True]) ==0:
            print("No transit data ready for GP modeling!!! Reverting to 'GP==False'") 
            fit.tra_doGP = False
            self.set_gui_use_GP()
        return    
                

 
    def get_cornerplot_param(self, type_plot = "mcmc"):
        global fit
        
        if type_plot == "mcmc":
            #if isinstance(fit.mcmc_sampler, rv.CustomSampler)==False:
            if len(np.atleast_1d(fit.mcmc_sampler))==0:            
                choice = QtWidgets.QMessageBox.information(self, 'Warning!', "MCMC samples not found.", QtWidgets.QMessageBox.StandardButton.Ok)
                return   
       
        if type_plot == "nest": 
            #if len(fit.ns_sampler)==0:
            if len(np.atleast_1d(fit.ns_sampler))==0:
            #if isinstance(fit.ns_sampler, dynesty.nestedsamplers.UnitCubeSampler)==False and isinstance(fit.ns_sampler,dynesty.dynamicsampler.DynamicSampler)==False:
                choice = QtWidgets.QMessageBox.information(self, 'Warning!', "NS samples not found.", QtWidgets.QMessageBox.StandardButton.Ok)
                return   
 
         
            
        if type_plot == "mcmc":
            self.lables_cornerplot = dill.copy(fit.mcmc_sampler.lbf)
        else:
            self.lables_cornerplot = dill.copy(fit.ns_sampler.lbf)
            
        label_results = self.dialog_select_param_cornerplot.get_labels(self)
        
        if type_plot == "mcmc":
            del fit.mcmc_sampler.lbf
            fit.mcmc_sampler.lbf  = dill.copy(label_results)
        else:
            fit.ns_sampler.lbf  = dill.copy(label_results)
 
        return  
   

    def transit_data_detrend(self):
        global fit
 
        but_ind = self.buttonGroup_detrend_tra.checkedId()
        
        if len(fit.tra_data_sets[but_ind-1]) != 0:
            self.tra_data = dill.copy(fit.tra_data_sets_init[but_ind-1])
            self.tra_data_index = but_ind-1
#           self.DetrendWindow.show()
            #self.DetrendWindow.worker_detrend()
            self.DetrendWindow.reset_data()
#            self.DetrendWindow.plot()
        else:
            
            print("No data ",but_ind)
            
            
    def rvs_data_options(self):
        global fit
 
        but_ind = self.buttonGroup_options_rvs.checkedId()
        
        print(but_ind-1)
        
        if len(fit.rv_data_sets[but_ind-1]) != 0:
            self.rvs_data = dill.copy(fit.rv_data_sets_init[but_ind-1])
            self.rvs_data_index = but_ind-1
            self.RvsWindow.reset_data()
        else:     
            print("No data ",but_ind)
                        

    def activity_data_options(self):
        global fit
 
        but_ind = self.buttonGroup_options_act.checkedId()
        
        print(but_ind-1)
        
        if len(fit.act_data_sets[but_ind-1]) != 0:
            self.act_data = dill.copy(fit.act_data_sets_init[but_ind-1])
            self.act_data_index = but_ind-1
            self.ActivityWindow.reset_data()
        else:     
            print("No data ",but_ind)
                        
            

    # Get variables pushed from the jupyter shell
    def get_jupyter_vars(self):
        global fit  
        fit = dill.copy(self.console_widget.kernel_manager.kernel.shell.user_ns.get('fit'))
 
        self.set_hkl()    
        self.set_aR()     
        self.update_use()
        self.update_gui_params()

        
    #NOT used anymore#
    def check_for_missing_instances(self,fit_new):
        global fit

        for iii in fit.__dict__:
            if iii not in fit_new.__dict__: 
                fit_new.__dict__[iii] = dill.copy(fit.__dict__[iii])

        for iii in fit.fit_results.__dict__:
            if iii not in fit_new.fit_results.__dict__: 
                fit_new.fit_results.__dict__[iii] = dill.copy(fit.fit_results.__dict__[iii])

        fit_new.cwd = dill.copy(fit.cwd)    

        return



    def set_type_fit_options(self):
        global fit 
        
        if self.radioButton_ttv.isChecked() == False:
            self.checkBox_first_TTV_epoch.setEnabled(False)
        else:
            self.checkBox_first_TTV_epoch.setEnabled(True)

        if self.radioButton_ast.isChecked() == False:
            self.checkBox_first_ast_epoch.setEnabled(False)
        else:
            self.checkBox_first_ast_epoch.setEnabled(True)

        
        self.mute_boxes()

    def check_type_fit(self):
        global fit  

        if fit.type_fit["RV"] == True and fit.type_fit["Transit"] == False and fit.type_fit["TTV"]  == False and fit.type_fit["AST"]  == False:
            self.radioButton_RV.setChecked(True)       
        elif fit.type_fit["RV"] == True and fit.type_fit["Transit"] == False and fit.type_fit["TTV"]  == True and fit.type_fit["AST"]  == False:
            self.radioButton_ttv_RV.setChecked(True)         
        elif fit.type_fit["RV"] == True and fit.type_fit["Transit"] == True and fit.type_fit["TTV"]  == False and fit.type_fit["AST"]  == False:
            self.radioButton_transit_RV.setChecked(True)                                 
        elif fit.type_fit["RV"] == False and fit.type_fit["Transit"] == False and fit.type_fit["TTV"]  == True and fit.type_fit["AST"]  == False:
            self.radioButton_ttv.setChecked(True)        
        elif fit.type_fit["RV"] == False and fit.type_fit["Transit"] == True and fit.type_fit["TTV"]  == False and fit.type_fit["AST"]  == False:
            self.radioButton_transit.setChecked(True)
        elif fit.type_fit["RV"] == False and fit.type_fit["Transit"] == False and fit.type_fit["TTV"]  == False and fit.type_fit["AST"]  == True:
            self.radioButton_ast.setChecked(True)            
        elif fit.type_fit["RV"] == True and fit.type_fit["Transit"] == False and fit.type_fit["TTV"]  == False and fit.type_fit["AST"]  == True:
            self.radioButton_ast_RV.setChecked(True)    

    def get_error_msg(self, msg):
        global fit  
        choice = QtWidgets.QMessageBox.information(self, 'Warning!', "%s"%str(msg), QtWidgets.QMessageBox.StandardButton.Ok)
        #return 
        
    def init_plot_corr(self):
        global fit  
  
        self.comboBox_samp_corr_1.clear()
        self.comboBox_samp_corr_2.clear()
            
       # print(self.samp_corr_mcmc.isChecked(), self.samp_corr_ns.isChecked())

        if len(np.atleast_1d(fit.mcmc_sampler))!=0 and self.samp_corr_mcmc.isChecked():       
         
            for i in range(len(fit.mcmc_sampler.lbf['OrigLabels'])):
                self.comboBox_samp_corr_1.addItem(fit.mcmc_sampler.lbf['OrigLabels'][i],i)
                self.comboBox_samp_corr_2.addItem(fit.mcmc_sampler.lbf['OrigLabels'][i],i)
                
        elif len(np.atleast_1d(fit.ns_sampler))!=0 and self.samp_corr_ns.isChecked():       
         
            for i in range(len(fit.ns_sampler.lbf['OrigLabels'])):
                self.comboBox_samp_corr_1.addItem(fit.ns_sampler.lbf['OrigLabels'][i],i)
                self.comboBox_samp_corr_2.addItem(fit.ns_sampler.lbf['OrigLabels'][i],i)


            #self.comboBox_pl_2.setItemText(i, '%s'%(ses_list[i].name))
        self.comboBox_samp_corr_1.setCurrentIndex(0)
        self.comboBox_samp_corr_2.setCurrentIndex(0)
        
        self.update_plot_corr()
        
    def update_plot_corr(self):
        global fit,pcor
 
        corr1_ind = self.comboBox_samp_corr_1.currentIndex()
        corr2_ind = self.comboBox_samp_corr_2.currentIndex()

 
        if corr1_ind ==-1 or corr2_ind ==-1: # or len(fit.e_for_mcmc) ==0 or len(np.atleast_1d(fit.mcmc_sampler))==0: #
            return
            
        if len(np.atleast_1d(fit.mcmc_sampler))!=0 and self.samp_corr_mcmc.isChecked():        
         

            pcor.plot(clear=True,)
            pcor.plot(fit.mcmc_sampler.samples[:,corr1_ind], fit.mcmc_sampler.samples[:,corr2_ind] ,pen=None, #{'color': colors[i], 'width': 1.1},
            symbol='o',
            symbolPen={'color': 'b', 'width': 1},
            symbolSize=1,enableAutoRange=True,viewRect=True,
            symbolBrush='b'
            )  

            pcor.setLabel('bottom', '%s'%fit.mcmc_sampler.lbf['OrigLabels'][corr1_ind], units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
            pcor.setLabel('left', '%s'%fit.mcmc_sampler.lbf['OrigLabels'][corr2_ind], units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
 


        elif len(np.atleast_1d(fit.ns_sampler))!=0 and self.samp_corr_ns.isChecked():            
 
            pcor.plot(clear=True,)
            pcor.plot(fit.ns_sampler.samples[:,corr1_ind], fit.ns_sampler.samples[:,corr2_ind] ,pen=None, #{'color': colors[i], 'width': 1.1},
            symbol='o',
            symbolPen={'color': 'b', 'width': 1},
            symbolSize=1,enableAutoRange=True,viewRect=True,
            symbolBrush='b'
            )  

            #pcor.setLabel('bottom', '%s'%fit.e_for_mcmc[corr1_ind], units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
            #pcor.setLabel('left', '%s'%fit.e_for_mcmc[corr2_ind], units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
            pcor.setLabel('bottom', '%s'%fit.ns_sampler.lbf['OrigLabels'][corr1_ind], units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
            pcor.setLabel('left', '%s'%fit.ns_sampler.lbf['OrigLabels'][corr2_ind], units='',  **{'font-size':'%dpt'%self.plot_font.pointSize()}) 
 
 
        else:
            return



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


    def set_force_copl_incl(self):
        global fit   
        fit.copl_incl = self.force_copl_incl.isChecked()
        self.mute_boxes()
        
    def set_force_sameRV_jit(self):
        global fit   
        fit.jit_flag = self.force_sameRV_jit.isChecked()
        self.mute_boxes()        


    #def update_inspector(self):
   #     self.tree_view_tab.listview.clicked.connect(self.plot_data_inspect)

    def eventFilter(self, obj, event):
        if event.type() == QtCore.QEvent.ToolTip:  # Catch the TouchBegin event.
           # helpEvent = event
           # index = self.itemAt(helpEvent.pos())
           # if index != -1:
           #     QtGui.QToolTip.showText(helpEvent.globalPos(), self.shapeItems[index].toolTip())
           # else:
            #print("TEST")
            QtGui.QToolTip.hideText()
            event.ignore()

            return True
 
        else:
            return super(Exo_striker, self).eventFilter(obj, event)


    def get_credits(self):
        if qso_mode:
            webbrowser.open('https://github.com/3fon3fonov/QSO-Striker')
        else:
            webbrowser.open('https://github.com/3fon3fonov/exostriker')

    def count_cpus(self):

        self.mlp_N_threads.setValue(cpu_count())
        #self.N_threads.setValue(cpu_count()) # somewhat slower with N cpus with Pathos multithreading....
        self.N_threads.setValue(1)        
        self.nest_N_threads.setValue(cpu_count())

 
 
################################################################################################


    

    def __init__(self):
        global fit 
        super().__init__()
    
        self.es_version = __version__

        #self.loading_screen= LoadingScreen()   
 
        QtWidgets.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
       # self.showMaximized()



        self.setupUi(self)



        self.initialize_font_plot()
        self.initialize_font_GUI()

        self.safe_to_init = True
#        self.installEventFilter(self)

        self.param_bounds_gui  = gui_groups.param_bounds_gui(self)
        self.offset_bounds_gui = gui_groups.offset_bounds_gui(self)
        self.jitter_bounds_gui = gui_groups.jitter_bounds_gui(self)
        self.offset_bounds_gui_tra = gui_groups.offset_bounds_gui_tra(self)
        self.jitter_bounds_gui_tra = gui_groups.jitter_bounds_gui_tra(self)

        self.tra_lintr_bounds_gui  = gui_groups.tra_lintr_bounds_gui(self)
        self.tra_quadtr_bounds_gui = gui_groups.tra_quadtr_bounds_gui(self)
        
        self.tra_lin_trend_nr_priors_gui    = gui_groups.tra_lin_trend_nr_priors_gui(self)
        self.tra_lin_trend_jeff_priors_gui  = gui_groups.tra_lin_trend_jeff_priors_gui(self)
        self.tra_quad_trend_nr_priors_gui   = gui_groups.tra_quad_trend_nr_priors_gui(self)
        self.tra_quad_trend_jeff_priors_gui = gui_groups.tra_quad_trend_jeff_priors_gui(self)

#        self.param_gui         = gui_groups.param_gui(self)
#        self.use_param_gui     = gui_groups.use_param_gui(self)
#        self.param_errors_gui  = gui_groups.param_errors_gui(self)

        self.param_gui_K         = gui_groups.param_gui_K(self)
        self.use_param_gui_K     = gui_groups.use_param_gui_K(self)
        self.param_errors_gui_K  = gui_groups.param_errors_gui_K(self)

        self.param_gui_P         = gui_groups.param_gui_P(self)
        self.use_param_gui_P     = gui_groups.use_param_gui_P(self)
        self.param_errors_gui_P  = gui_groups.param_errors_gui_P(self)

        self.param_gui_e         = gui_groups.param_gui_e(self)
        self.use_param_gui_e     = gui_groups.use_param_gui_e(self)
        self.param_errors_gui_e  = gui_groups.param_errors_gui_e(self)

        self.param_gui_om         = gui_groups.param_gui_om(self)
        self.use_param_gui_om     = gui_groups.use_param_gui_om(self)
        self.param_errors_gui_om  = gui_groups.param_errors_gui_om(self)

        self.param_gui_ma         = gui_groups.param_gui_ma(self)
        self.use_param_gui_ma     = gui_groups.use_param_gui_ma(self)
        self.param_errors_gui_ma  = gui_groups.param_errors_gui_ma(self)

        self.param_gui_incl         = gui_groups.param_gui_incl (self)
        self.use_param_gui_incl      = gui_groups.use_param_gui_incl (self)
        self.param_errors_gui_incl   = gui_groups.param_errors_gui_incl (self)

        self.param_gui_Omega         = gui_groups.param_gui_Omega(self)
        self.use_param_gui_Omega     = gui_groups.use_param_gui_Omega(self)
        self.param_errors_gui_Omega  = gui_groups.param_errors_gui_Omega(self)

        self.param_gui_wd      = gui_groups.param_gui_wd(self)
        self.use_param_gui_wd  = gui_groups.use_param_gui_wd(self)
        self.param_errors_gui_wd = gui_groups.param_errors_gui_wd(self)

        ########## Astr.############### 
        self.param_gui_ast      = gui_groups.param_gui_ast(self)
        self.use_param_gui_ast  = gui_groups.use_param_gui_ast(self)
        self.param_errors_gui_ast = gui_groups.param_errors_gui_ast(self)

        self.ast_bounds_gui = gui_groups.ast_bounds_gui(self)
        self.ast_norm_pr_gui = gui_groups.ast_norm_pr_gui(self)
        self.ast_jeff_pr_gui = gui_groups.ast_jeff_pr_gui(self)
        ################################
        
        self.om_dot_bounds_gui = gui_groups.om_dot_bounds_gui(self)
        self.om_dot_norm_pr_gui = gui_groups.om_dot_norm_pr_gui(self)
        self.om_dot_jeff_pr_gui = gui_groups.om_dot_jeff_pr_gui(self)

        self.h_bounds_gui = gui_groups.h_bounds_gui(self)
        self.h_norm_pr_gui = gui_groups.h_norm_pr_gui(self)
        self.h_jeff_pr_gui = gui_groups.h_jeff_pr_gui(self)

        self.k_bounds_gui = gui_groups.k_bounds_gui(self)
        self.k_norm_pr_gui = gui_groups.k_norm_pr_gui(self)
        self.k_jeff_pr_gui = gui_groups.k_jeff_pr_gui(self)

        self.lambda_bounds_gui = gui_groups.lambda_bounds_gui(self)
        self.lambda_norm_pr_gui = gui_groups.lambda_norm_pr_gui(self)
        self.lambda_jeff_pr_gui = gui_groups.lambda_jeff_pr_gui(self)
        
        self.param_gui_tr      = gui_groups.param_gui_tr(self)
        self.use_param_gui_tr  = gui_groups.use_param_gui_tr(self)

        self.err_t0          = gui_groups.err_t0(self)
        self.err_pl_rad      = gui_groups.err_pl_rad(self)
        self.err_a_sol       = gui_groups.err_a_sol(self)
 
        self.rvs_data_gui      = gui_groups.rvs_data_gui(self)
        self.rvs_data_jitter_gui = gui_groups.rvs_data_jitter_gui(self)
        self.use_data_offset_gui = gui_groups.use_data_offset_gui(self)
        self.use_data_jitter_gui = gui_groups.use_data_jitter_gui(self)
        
        self.data_errors_gui            = gui_groups.data_errors_gui(self)
        self.data_errors_jitter_gui     = gui_groups.data_errors_jitter_gui(self)
        
        self.tra_data_gui             = gui_groups.tra_data_gui(self)
        self.tra_data_jitter_gui      = gui_groups.tra_data_jitter_gui(self)
        
        self.use_tra_data_offset_gui  = gui_groups.use_tra_data_offset_gui(self)
        self.use_tra_data_jitter_gui  = gui_groups.use_tra_data_jitter_gui(self)
        
        self.tra_data_errors_gui        = gui_groups.tra_data_errors_gui(self)
        self.tra_data_errors_jitter_gui = gui_groups.tra_data_errors_jitter_gui(self)
        

        self.tra_data_lin_trend_gui       = gui_groups.tra_data_lin_trend_gui(self)
        self.use_tra_data_lin_trend_gui   = gui_groups.use_tra_data_lin_trend_gui(self)
        self.err_tra_data_lin_trend_gui   = gui_groups.err_tra_data_lin_trend_gui(self)
        
        self.tra_data_quad_trend_gui      = gui_groups.tra_data_quad_trend_gui(self)
        self.use_tra_data_quad_trend_gui  = gui_groups.use_tra_data_quad_trend_gui(self)
        self.err_tra_data_quad_trend_gui  = gui_groups.err_tra_data_quad_trend_gui(self)
 

        self.gp_rot_params     = gui_groups.gp_rot_params(self)
        self.use_gp_rot_params = gui_groups.use_gp_rot_params(self)
        self.gp_rot_errors_gui = gui_groups.gp_rot_errors_gui(self)
        
        self.gp_sho_params     = gui_groups.gp_sho_params(self)
        self.use_gp_sho_params = gui_groups.use_gp_sho_params(self)
        self.gp_sho_errors_gui = gui_groups.gp_sho_errors_gui(self)
 
        self.gp_double_sho_params     = gui_groups.gp_double_sho_params(self)
        self.use_gp_double_sho_params = gui_groups.use_gp_double_sho_params(self)
        self.gp_double_sho_errors_gui = gui_groups.gp_double_sho_errors_gui(self)
       
        self.gp_mat_params     = gui_groups.gp_mat_params(self)
        self.use_gp_mat_params = gui_groups.use_gp_mat_params(self)
        self.gp_mat_errors_gui = gui_groups.gp_mat_errors_gui(self)        
        
        self.gp_drw_params     = gui_groups.gp_drw_params(self)
        self.use_gp_drw_params = gui_groups.use_gp_drw_params(self)
        self.gp_drw_errors_gui = gui_groups.gp_drw_errors_gui(self)  
        
        self.tra_gp_rot_params = gui_groups.tra_gp_rot_params(self)
        self.use_tra_gp_rot_params = gui_groups.use_tra_gp_rot_params(self)
        self.tra_gp_rot_errors_gui = gui_groups.tra_gp_rot_errors_gui(self)

        self.tra_gp_sho_params     = gui_groups.tra_gp_sho_params(self)
        self.use_tra_gp_sho_params = gui_groups.use_tra_gp_sho_params(self)
        self.gp_tra_sho_errors_gui = gui_groups.gp_tra_sho_errors_gui(self)
        
        self.tra_gp_mat_params     = gui_groups.tra_gp_mat_params(self)
        self.use_tra_gp_mat_params = gui_groups.use_tra_gp_mat_params(self)        
        self.tra_gp_mat_errors_gui = gui_groups.tra_gp_mat_errors_gui(self)      
        
        self.tra_gp_drw_params     = gui_groups.tra_gp_drw_params(self)
        self.use_tra_gp_drw_params = gui_groups.use_tra_gp_drw_params(self)        
        self.tra_gp_drw_errors_gui = gui_groups.tra_gp_drw_errors_gui(self)          

        self.tra_gp_double_sho_params     = gui_groups.tra_gp_double_sho_params(self)
        self.use_tra_gp_double_sho_params = gui_groups.use_tra_gp_double_sho_params(self)
        self.tra_gp_double_sho_errors_gui = gui_groups.tra_gp_double_sho_errors_gui(self)        

        self.GP_sho_bounds_gui = gui_groups.GP_sho_bounds_gui(self) 
        self.GP_double_sho_bounds_gui = gui_groups.GP_double_sho_bounds_gui(self) 
        self.GP_rot_bounds_gui = gui_groups.GP_rot_bounds_gui(self) 
        self.GP_mat_bounds_gui = gui_groups.GP_mat_bounds_gui(self) 
        self.GP_drw_bounds_gui = gui_groups.GP_drw_bounds_gui(self) 

        self.tra_GP_sho_bounds_gui = gui_groups.tra_GP_sho_bounds_gui(self) 
        self.tra_GP_double_sho_bounds_gui = gui_groups.tra_GP_double_sho_bounds_gui(self) 
        self.tra_GP_rot_bounds_gui = gui_groups.tra_GP_rot_bounds_gui(self) 
        self.tra_GP_mat_bounds_gui = gui_groups.tra_GP_mat_bounds_gui(self) 
        self.tra_GP_drw_bounds_gui = gui_groups.tra_GP_drw_bounds_gui(self) 

        self.GP_sho_nr_priors_gui = gui_groups.GP_sho_nr_priors_gui(self) 
        self.GP_double_sho_nr_priors_gui = gui_groups.GP_double_sho_nr_priors_gui(self) 
        self.GP_rot_nr_priors_gui = gui_groups.GP_rot_nr_priors_gui(self) 
        self.GP_mat_nr_priors_gui = gui_groups.GP_mat_nr_priors_gui(self) 
        self.GP_drw_nr_priors_gui = gui_groups.GP_drw_nr_priors_gui(self) 

        self.tra_GP_sho_nr_priors_gui = gui_groups.tra_GP_sho_nr_priors_gui(self) 
        self.tra_GP_double_sho_nr_priors_gui = gui_groups.tra_GP_double_sho_nr_priors_gui(self) 
        self.tra_GP_rot_nr_priors_gui = gui_groups.tra_GP_rot_nr_priors_gui(self) 
        self.tra_GP_mat_nr_priors_gui = gui_groups.tra_GP_mat_nr_priors_gui(self) 
        self.tra_GP_drw_nr_priors_gui = gui_groups.tra_GP_drw_nr_priors_gui(self) 
        
        self.GP_sho_jeff_priors_gui = gui_groups.GP_sho_jeff_priors_gui(self) 
        self.GP_double_sho_jeff_priors_gui = gui_groups.GP_double_sho_jeff_priors_gui(self) 
        self.GP_rot_jeff_priors_gui = gui_groups.GP_rot_jeff_priors_gui(self) 
        self.GP_mat_jeff_priors_gui = gui_groups.GP_mat_jeff_priors_gui(self) 
        self.GP_drw_jeff_priors_gui = gui_groups.GP_drw_jeff_priors_gui(self) 

        self.tra_GP_sho_jeff_priors_gui = gui_groups.tra_GP_sho_jeff_priors_gui(self) 
        self.tra_GP_double_sho_jeff_priors_gui = gui_groups.tra_GP_double_sho_jeff_priors_gui(self) 
        self.tra_GP_rot_jeff_priors_gui = gui_groups.tra_GP_rot_jeff_priors_gui(self) 
        self.tra_GP_mat_jeff_priors_gui = gui_groups.tra_GP_mat_jeff_priors_gui(self) 
        self.tra_GP_drw_jeff_priors_gui = gui_groups.tra_GP_drw_jeff_priors_gui(self) 
 

        self.param_a_gui      = gui_groups.param_a_gui(self)
        self.param_mass_gui   = gui_groups.param_mass_gui(self)
        self.param_t_peri_gui = gui_groups.param_t_peri_gui(self)
        
        self.planet_checked_gui = gui_groups.planet_checked_gui(self)
        
        self.arb_param_gui     = gui_groups.arb_param_gui(self)
        self.arb_param_gui_use = gui_groups.arb_param_gui_use(self)
        
 #       self.add_rv_error   = gui_groups.add_rv_error(self)
 #       self.rv_sigma_clip  = gui_groups.rv_sigma_clip(self)
 #       self.bin_rv_data    = gui_groups.bin_rv_data(self)
        
        #self.act_sigma_clip  = gui_groups.act_sigma_clip(self)
        #self.act_remove_mean = gui_groups.act_remove_mean(self)

#        self.tra_sigma_clip  = gui_groups.tra_sigma_clip(self)
        self.tra_norm        = gui_groups.tra_norm(self)
        self.tra_dilution    = gui_groups.tra_dilution(self)

        self.rvs_opt         = gui_groups.rvs_opt(self)

        self.act_opt         = gui_groups.act_opt(self)
        self.use_tra_data_GP = gui_groups.use_tra_data_GP(self)

        self.ttv_data_to_planet     = gui_groups.ttv_data_to_planet(self)
        self.use_ttv_data_to_planet = gui_groups.use_ttv_data_to_planet(self)

        self.ast_data_to_planet     = gui_groups.ast_data_to_planet(self)
        self.use_ast_data_to_planet = gui_groups.use_ast_data_to_planet(self)
  
        self.ast_data_to_planet_2     = gui_groups.ast_data_to_planet_2(self)
        self.use_ast_data_to_planet_2 = gui_groups.use_ast_data_to_planet_2(self)

        self.use_uni_ld_models    = gui_groups.use_uni_ld_models(self)
        self.use_lin_ld_models    = gui_groups.use_lin_ld_models(self)
        self.use_quad_ld_models   = gui_groups.use_quad_ld_models(self)
        self.use_nonlin_ld_models = gui_groups.use_nonlin_ld_models(self)

        self.lin_u                = gui_groups.lin_u(self)
        self.quad_u1              = gui_groups.quad_u1(self)
        self.quad_u2              = gui_groups.quad_u2(self)
        self.nonlin_u1            = gui_groups.nonlin_u1(self)
        self.nonlin_u2            = gui_groups.nonlin_u2(self)
        self.nonlin_u3            = gui_groups.nonlin_u3(self)
        self.nonlin_u4            = gui_groups.nonlin_u4(self)

        self.use_lin_u            = gui_groups.use_lin_u(self)
        self.use_quad_u1          = gui_groups.use_quad_u1(self)
        self.use_quad_u2          = gui_groups.use_quad_u2(self)
        self.use_nonlin_u1        = gui_groups.use_nonlin_u1(self)
        self.use_nonlin_u2        = gui_groups.use_nonlin_u2(self)
        self.use_nonlin_u3        = gui_groups.use_nonlin_u3(self)
        self.use_nonlin_u4        = gui_groups.use_nonlin_u4(self)

        self.err_lin_u            = gui_groups.err_lin_u(self)
        self.err_quad_u1          = gui_groups.err_quad_u1(self)
        self.err_quad_u2          = gui_groups.err_quad_u2(self)
        self.err_nonlin_u1        = gui_groups.err_nonlin_u1(self)
        self.err_nonlin_u2        = gui_groups.err_nonlin_u2(self)
        self.err_nonlin_u3        = gui_groups.err_nonlin_u3(self)
        self.err_nonlin_u4        = gui_groups.err_nonlin_u4(self)

        self.ld_u1_bounds_gui     = gui_groups.ld_u1_bounds_gui(self)
        self.ld_u2_bounds_gui     = gui_groups.ld_u2_bounds_gui(self)
        self.ld_u3_bounds_gui     = gui_groups.ld_u3_bounds_gui(self)
        self.ld_u4_bounds_gui     = gui_groups.ld_u4_bounds_gui(self)

        self.ld_u1_norm_pr_gui     = gui_groups.ld_u1_norm_pr_gui(self)
        self.ld_u2_norm_pr_gui     = gui_groups.ld_u2_norm_pr_gui(self)
        self.ld_u3_norm_pr_gui     = gui_groups.ld_u3_norm_pr_gui(self)
        self.ld_u4_norm_pr_gui     = gui_groups.ld_u4_norm_pr_gui(self)

        self.ld_u1_jeff_pr_gui     = gui_groups.ld_u1_jeff_pr_gui(self)
        self.ld_u2_jeff_pr_gui     = gui_groups.ld_u2_jeff_pr_gui(self)
        self.ld_u3_jeff_pr_gui     = gui_groups.ld_u3_jeff_pr_gui(self)
        self.ld_u4_jeff_pr_gui     = gui_groups.ld_u4_jeff_pr_gui(self)
                
        
        self.data_ld_group        = gui_groups.data_ld_group(self)

        self.data_tra_reg_group   = gui_groups.data_tra_reg_group(self)


        self.param_nr_priors_gui = gui_groups.param_nr_priors_gui(self)
        self.param_jeff_priors_gui = gui_groups.param_jeff_priors_gui(self)

        
        self.offset_nr_priors_gui = gui_groups.offset_nr_priors_gui(self)
        self.jitter_nr_priors_gui = gui_groups.jitter_nr_priors_gui(self)
        self.offset_jeff_priors_gui = gui_groups.offset_jeff_priors_gui(self)
        self.jitter_jeff_priors_gui = gui_groups.jitter_jeff_priors_gui(self)


        self.offset_nr_priors_gui_tra = gui_groups.offset_nr_priors_gui_tra(self)
        self.jitter_nr_priors_gui_tra = gui_groups.jitter_nr_priors_gui_tra(self)
        self.offset_jeff_priors_gui_tra = gui_groups.offset_jeff_priors_gui_tra(self)
        self.jitter_jeff_priors_gui_tra = gui_groups.jitter_jeff_priors_gui_tra(self)

        ########### TEMP ##############
        self.TTV_readme_info.clicked.connect(lambda: self.print_TTV_info()) 
        self.ast_readme_info.clicked.connect(lambda: self.print_ast_info()) 
        self.ast_readme_info_2.clicked.connect(lambda: self.print_ast_info_2()) 
        
        
        self.initialize_buttons()
        self.initialize_plots()   
 
        self.initialize_color_dialog()
        #Hill_LED = './lib/UI/grey_led.png'
        #self.Hill_led.setPixmap(QtGui.QPixmap(Hill_LED))
        AMD_LED = './lib/UI/grey_led.png'
        self.AMD_led.setPixmap(QtGui.QPixmap(AMD_LED))
        
        self.check_type_fit()
        
        ###################### Info buttons #############################
    
        self.RV_GP_Rot_readme_info.clicked.connect(self.print_GP_info)
        self.RV_GP_SHO_readme_info.clicked.connect(self.print_GP_info)
        self.RV_GP_Matern_readme_info.clicked.connect(self.print_GP_info)
        self.RV_GP_double_SHO_readme_info .clicked.connect(self.print_GP_info)  
        self.RV_GP_DRW_readme_info.clicked.connect(self.print_GP_info)


        self.Tra_GP_Rot_readme_info.clicked.connect(self.print_GP_info)
        self.Tra_GP_SHO_readme_info.clicked.connect(self.print_GP_info)
        self.Tra_GP_Matern_readme_info.clicked.connect(self.print_GP_info)
        self.Tra_GP_double_SHO_readme_info.clicked.connect(self.print_GP_info)  
        self.Tra_GP_DRW_readme_info.clicked.connect(self.print_GP_info)  
        
        
           
        self.mcmc_readme_info.clicked.connect(self.print_mcmc_info)
        self.ns_readme_info.clicked.connect(self.print_ns_info)
        
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
        #self.console_widget.execute_command("%matplotlib")
        
        self.console_widget = ConsoleWidget_embed(font_size = 10)
        
        self.terminal_embeded.addTab(self.console_widget, "Jupyter")
        
        ###################### Console #############################
        
       
        #self.terminal_embeded.addTab(self.tree_view_tab, "tree")
      
        
        if sys.platform[0:5] == "linux":
            self.term_emb = terminal.mainWindow()
            self.terminal_embeded.addTab(self.term_emb, "Bash shell")

        self.terminal_embeded.addTab(pg_console.ConsoleWidget(), "pqg shell")

        if openai_not_found != True:
            self.chat_widget = ChatWidget()
            self.terminal_embeded.addTab(self.chat_widget, "Exo-Striker ChatBoot (GPT-3)")
        else:
            print("You do not have openAI installed, so You cannot use the GPT-3 chatboot Try $ pip install openai")


        self.text_editor = text_editor_es.MainWindow()
        self.gridLayout_text_editor.addWidget(self.text_editor)
        self.calculator = calc.Calculator()
        self.gridLayout_calculator.addWidget(self.calculator)
        
        #################### data inspector ########################
        self.inspector_file = ""
        self.RVBank = False
        self.datafiles_window = datafiles_window()
        self.RVBank_window = RVBank_window() 
        #self.tree_view_tab = Widget_tree()        
       # self.gridLayout_file_tree.setRowStretch(0, 6)
        #self.gridLayout_file_tree.setRowStretch(1, 4)
        #self.gridLayout_file_tree.addWidget(self.tree_view_tab)


        self.export_DataInspector_data.clicked.connect(self.save_DataInspector_data)
        self.data_insp_load_data.clicked.connect(self.load_data_inspect)
        self.datafiles_window.listview.clicked.connect(self.plot_data_inspect)

        self.dataInspector_thisComputer.clicked.connect(self.datafiles_window.show)
        #self.dataInspector_HARPS_RVBank.clicked.connect(lambda: self.get_error_msg("Still work in progress! <br><br>However, you can inspect the online version of the <a href='http://www.mpia.de/homes/trifonov/HARPS_RVBank.html'>HARPS RVBank</a> and download HARPS data products. Then in the Exo-Striker use:<br><br>File --> Open RVBank file --> <br>(and then load the downoladed .dat file) <br><br>This will load you target's HARPS NZP corrected RVs and all the activity index data. <br><br>If you made use of the HARPS RVBank, please do not forget to cite <a href='https://ui.adsabs.harvard.edu/abs/2020arXiv200105942T/abstract'>Trifonov et al. (2020)</a>"))

        self.RVBank_window.list.clicked.connect(lambda: self.plot_data_inspect(0,RVBank = True))
        self.RVBank_window.list_opt.clicked.connect(lambda: self.plot_data_inspect(0,RVBank = True))
        self.RVBank_window.radio_group.buttonClicked.connect(lambda: self.plot_data_inspect(0,RVBank = True))

        self.dataInspector_HARPS_RVBank.clicked.connect(self.RVBank_window.show)



        #################### stdout pipe ########################

        #if sys.version_info[0] == 3:
        if debug == False:
            self.pipe_text = MyDialog()
        else:
            try:
                import pkg_resources 
                     # list packages to be checked 
                root_packages = [ 
                    'numpy', 
                    'scipy', 
                    'matplotlib',
                    'pyqt6',  
                    'qtconsole',
                    'jupyter-client',
                    'ipykernel',
                    'jupyter', 
                    'pathos', 
                    'emcee',
                    'celerite', 
                    #'corner',
                    'transitleastsquares',
                   # 'dynesty', 
                    'ttvfast',
                    'wotan'] 
                
                
                #root_packages.sort(reverse=True)
                     # print versions, but check if package is imported first 
                for m in pkg_resources.working_set: 
                    if m.project_name.lower() in root_packages: 
                        print("%s==%s"%(m.project_name,m.version))
                    #else:
                    #    print(f"{m.project_name}== Not found!")

            except:
                print(" It seems that you have missing packages!")
                pass

            self.pipe_text = DebugDialog()

        self.gridLayout_stdout.addWidget(self.pipe_text)
        
        
        self.buttonGroup_type_fit.buttonClicked.connect(self.set_type_fit_options)
        
        
        self.buttonGroup_use_tra_data_GP.buttonClicked.connect(self.set_use_tra_GP_data)

        #################### credits  ########################
        
        self.dialog           = print_info(self)
        self.dialog_credits   = print_info(self)
        self.dialog_chi_table = print_info(self)
        self.dialog_ttv_help  = print_info(self)
        self.dialog_ast_help  = print_info(self)
        self.dialog_GP_help   = print_info(self)
        self.dialog_mcmc_help = print_info(self)
        self.dialog_ns_help   = print_info(self)
        
        self.dialog_more_info = print_info(self)



        self.More_info.clicked.connect(self.print_more_stat)

#        self.buttonGroup_apply_rv_data_options.buttonClicked.connect(self.apply_rv_data_options)
        self.buttonGroup_apply_act_data_options.buttonClicked.connect(self.apply_act_data_options)
        self.buttonGroup_apply_rvs_data_options.buttonClicked.connect(self.apply_rvs_data_options)
        self.buttonGroup_apply_tra_data_options.buttonClicked.connect(self.apply_tra_data_options)
        self.buttonGroup_apply_tra_dilution.buttonClicked.connect(self.apply_tra_dilution)


        self.buttonGroup_add_RV_data.buttonClicked.connect(self.showDialog_RV_input_file)
        self.buttonGroup_remove_RV_data.buttonClicked.connect(self.remove_RV_file)
 
        self.buttonGroup_activity_data.buttonClicked.connect(self.showDialog_act_input_file)
        self.buttonGroup_remove_activity_data.buttonClicked.connect(self.remove_act_file)     



        ############### TTV plotting controll ####################     
 
        #self.ttv_data_size.valueChanged.connect(self.update_ttv_plots)

        self.init_scipy_ttv_interpolation()
        self.ttv_model_width.valueChanged.connect(self.update_ttv_plots)
        self.ttv_model_z.valueChanged.connect(self.update_ttv_plots)

        self.ttv_model_interpolate_points.valueChanged.connect(self.update_ttv_plots)
        self.ttv_model_interpolate.stateChanged.connect(self.update_ttv_plots)

        self.comboBox_ttv_model_interpolate.activated.connect(self.update_ttv_plots)

        self.buttonGroup_ttv_data.buttonClicked.connect(self.showDialog_ttv_input_file)
        self.buttonGroup_remove_ttv_data.buttonClicked.connect(self.remove_ttv_file)
        self.buttonGroup_use_ttv_data_to_planet.buttonClicked.connect(self.use_ttv_dataset_to_planet)
        
        self.ttv_data_planet_1.valueChanged.connect(self.ttv_dataset_to_planet)
        self.ttv_data_planet_2.valueChanged.connect(self.ttv_dataset_to_planet)
        self.ttv_data_planet_3.valueChanged.connect(self.ttv_dataset_to_planet)
        self.ttv_data_planet_4.valueChanged.connect(self.ttv_dataset_to_planet)
        self.ttv_data_planet_5.valueChanged.connect(self.ttv_dataset_to_planet)
        self.ttv_data_planet_6.valueChanged.connect(self.ttv_dataset_to_planet)
        self.ttv_data_planet_7.valueChanged.connect(self.ttv_dataset_to_planet)
        self.ttv_data_planet_8.valueChanged.connect(self.ttv_dataset_to_planet)
        self.ttv_data_planet_9.valueChanged.connect(self.ttv_dataset_to_planet)
        self.ttv_data_planet_10.valueChanged.connect(self.ttv_dataset_to_planet)

        self.ttv_pl_combo()
        self.ttv_comboBox_pl.activated.connect(lambda: self.update_ttv_pl_index(o_c =False))
        self.ttv_o_c_comboBox_pl.activated.connect(lambda: self.update_ttv_pl_index(o_c =True))

        self.ttv_plot_cross_hair.stateChanged.connect(self.update_ttv_plots)
        self.ttv_o_c_plot_cross_hair.stateChanged.connect(self.update_ttv_plots)

        self.ttv_o_c.stateChanged.connect(self.update_ttv_plots) 

        self.ttv_apply_mean_period.stateChanged.connect(self.update_ttv_plots)

        self.ttv_subtract_mean.stateChanged.connect(self.update_ttv_plots)
        self.ttv_plot_autorange.stateChanged.connect(self.update_ttv_plots)
        self.ttv_legend.stateChanged.connect(self.update_ttv_plots)


        ############### Astr. plotting controll ####################     

        #self.init_scipy_ast_interpolation()
        self.ast_model_width.valueChanged.connect(self.update_ast_plots)
        self.ast_model_z.valueChanged.connect(self.update_ast_plots)

        self.ast_data_color_gradient.stateChanged.connect(self.update_ast_plots)

        #self.ast_model_interpolate_points.valueChanged.connect(self.update_ast_plots)
        #self.ast_model_interpolate.stateChanged.connect(self.update_ast_plots)

        #self.comboBox_ast_model_interpolate.activated.connect(self.update_ast_plots)
        
        self.ast_model_color_gradient.stateChanged.connect(self.update_ast_plots)

        self.ast_o_c_hipp_gaia.stateChanged.connect(self.update_ast_plots)

        self.buttonGroup_ast_data.buttonClicked.connect(self.showDialog_ast_input_file)
        self.buttonGroup_ast_data_2.buttonClicked.connect(self.showDialog_ast_input_file_2)        
        
        
        self.buttonGroup_remove_ast_data.buttonClicked.connect(self.remove_ast_file)
        self.buttonGroup_remove_ast_data_2.buttonClicked.connect(self.remove_ast_file_2)        
        
        self.buttonGroup_use_ast_data_to_planet.buttonClicked.connect(self.use_ast_dataset_to_planet)
        self.buttonGroup_use_ast_data_to_planet_2.buttonClicked.connect(self.use_ast_dataset_to_planet_2)
                
        self.ast_data_planet_1.valueChanged.connect(self.ast_dataset_to_planet)
        self.ast_data_planet_2.valueChanged.connect(self.ast_dataset_to_planet)
        self.ast_data_planet_3.valueChanged.connect(self.ast_dataset_to_planet)
        self.ast_data_planet_4.valueChanged.connect(self.ast_dataset_to_planet)
        self.ast_data_planet_5.valueChanged.connect(self.ast_dataset_to_planet)
        self.ast_data_planet_6.valueChanged.connect(self.ast_dataset_to_planet)
        self.ast_data_planet_7.valueChanged.connect(self.ast_dataset_to_planet)
        self.ast_data_planet_8.valueChanged.connect(self.ast_dataset_to_planet)
        self.ast_data_planet_9.valueChanged.connect(self.ast_dataset_to_planet)
        self.ast_data_planet_10.valueChanged.connect(self.ast_dataset_to_planet)


        self.ast_data_planet_hipp_1.valueChanged.connect(self.ast_dataset_to_planet_2)
        self.ast_data_planet_hipp_2.valueChanged.connect(self.ast_dataset_to_planet_2)
        self.ast_data_planet_gaia_1.valueChanged.connect(self.ast_dataset_to_planet_2)        
        self.ast_data_planet_gaia_2.valueChanged.connect(self.ast_dataset_to_planet_2)

        self.ast_pl_combo()
        self.ast_comboBox_pl.activated.connect(lambda: self.update_ast_pl_index(o_c =False))
        self.ast_o_c_comboBox_pl.activated.connect(lambda: self.update_ast_pl_index(o_c =True))

        self.ast_plot_cross_hair.stateChanged.connect(self.update_ast_plots)
        self.ast_o_c_plot_cross_hair.stateChanged.connect(self.update_ast_plots)
        self.ast_plot_cross_hair_hipp_gaia.stateChanged.connect(self.update_ast_plots)
        
        self.ast_use_plot_ra_dec.stateChanged.connect(self.update_ast_plots)         
        self.ast_o_c.stateChanged.connect(self.update_ast_plots) 

        self.buttonGroup_use_ast_plot_RA_Dec.buttonClicked.connect(self.update_ast_plots)
        #self.ast_apply_mean_period.stateChanged.connect(self.update_ast_plots)

        #self.ast_subtract_mean.stateChanged.connect(self.update_ast_plots)
        self.ast_plot_autorange.stateChanged.connect(self.update_ast_plots)
        self.ast_legend.stateChanged.connect(self.update_ast_plots)

        self.buttonGroup_use_img_hipp_ast.buttonClicked.connect(self.set_use_astrometry)

        ############################################################


        self.buttonGroup_transit_data.buttonClicked.connect(self.showDialog_tra_input_file)
        self.buttonGroup_remove_transit_data.buttonClicked.connect(self.remove_tra_file)

        self.buttonGroup_use.buttonClicked.connect(self.update_use)
        self.buttonGroup_mixed_fitting.buttonClicked.connect(self.update_mixed_fitting)


        self.button_orb_evol.clicked.connect(self.worker_Nbody) 
        self.button_MCMC.clicked.connect(self.worker_mcmc)
        
        self.buttonGroup_adopt_mcmc_as_par.buttonClicked.connect(self.check_mcmc_params)
        self.buttonGroup_adopt_nest_as_par.buttonClicked.connect(self.check_nested_params)
      

       # self.button_nest_samp.clicked.connect(lambda: self.run_nest_samp())
        self.button_nest_samp.clicked.connect(self.worker_nest)
        
        self.run_orb_evol_arbitary.clicked.connect(self.worker_Nbody_arb) 
 
        self.button_make_mcmc_cornerplot_redir.clicked.connect(self.switch_to_corenrplot_opt)
        self.button_make_mcmc_cornerplot_redir_2.clicked.connect(self.switch_to_corenrplot_opt)

        self.remove_mcmc_samples.clicked.connect(self.remove_mcmc_samples_from_fit)
        self.remove_ns_samples.clicked.connect(self.remove_ns_samples_from_fit)
        
        self.button_make_mcmc_cornerplot.clicked.connect(lambda: self.worker_cornerplot(type_plot = "mcmc"))
        self.button_make_nest_cornerplot.clicked.connect(lambda: self.worker_cornerplot(type_plot = "nest"))
        
        self.mcmc_corner_plot_change_name.clicked.connect(lambda: self.change_corner_plot_file_name(type_plot = "mcmc"))
        self.nest_corner_plot_change_name.clicked.connect(lambda: self.change_corner_plot_file_name(type_plot = "nest"))
        
        self.mcmc_samples_change_name.clicked.connect(self.change_mcmc_samples_file_name)
        self.nest_samples_change_name.clicked.connect(self.change_nest_samples_file_name)
        
        self.comboBox_samp_corr_1.activated.connect(self.update_plot_corr)
        self.comboBox_samp_corr_2.activated.connect(self.update_plot_corr)
        self.buttonGroup_samp_corr.buttonClicked.connect(self.init_plot_corr)
                
        self.comboBox_phase_pl_tran.activated.connect(self.update_transit_plots)

        ########## RV fitting ########################

        self.button_init_fit.clicked.connect(lambda: self.fit_dispatcher(init=True))
        self.button_fit.clicked.connect(lambda: self.fit_dispatcher())
        self.button_auto_fit.clicked.connect(self.run_auto_fit)
        self.minimize_1param()

        self.radioButton_Dynamical.toggled.connect(self.update_dyn_kep_flag)
        self.radioButton_Dynamical.toggled.connect(self.mute_boxes)

        self.radioButton_omega_dot_free.toggled.connect(self.mute_boxes)
#        self.radioButton_omega_dot_GR.toggled.connect(self.set_gr_flag)



        self.radioButton_Keplerian.toggled.connect(self.update_dyn_kep_flag)
        
        self.deattach_omega_dot.stateChanged.connect(self.mute_boxes)
                
        self.amoeba_radio_button.toggled.connect(self.update_RV_jitter_flag)
        self.lm_radio_button.toggled.connect(self.update_RV_jitter_flag)       
        self.radioButton_scipy.toggled.connect(self.update_RV_jitter_flag)       
        self.radioButton_fortran77.toggled.connect(self.update_RV_jitter_flag)       

        self.radioButton_Keplerian.toggled.connect(self.mute_boxes)
        
        
        self.adopt_best_RV__o_c_GLS_per.clicked.connect(self.adopt_RV_GLS_param)
        self.adopt_mpl_param.clicked.connect(self.adopt_RV_mlp_param)        
        self.adopt_tls_o_c_param.clicked.connect(self.adopt_trans_TLS_param)


        ############ Set Widget Style #################

        self.actionWindows.triggered.connect(self.set_Win_widget_Style)
        self.actionMacintosh.triggered.connect(self.set_Mac_widget_Style)
        self.actionLinux_Fusion.triggered.connect(self.set_Fus_widget_Style)

        self.actionSet_GUI_Font.triggered.connect(self.set_widget_font)
        self.actionSet_plots_font.triggered.connect(self.set_plot_font)


        self.actionGet_all_plots.triggered.connect(self.get_pyqtgraph_plots)



        ############ Edit #################

        self.actionRV_Auto_fit.triggered.connect(self.run_auto_fit)
        self.actionReset_Buttons_dev_mode.triggered.connect(self.reset_mid_pannel_buttons)

        ############ View #################

        self.actiongrab_screen.triggered.connect(self.grab_screen) 
        self.actionprint_f_test_FAP.triggered.connect(self.print_f_test_stat)
        self.actionGet_LaTeX_table_with_parameters.triggered.connect(self.get_latex_table)
        self.actionGet_LaTex_table_with_priors.triggered.connect(self.get_latex_prior_table)
        self.actionGet_RV_model.triggered.connect(self.get_RV_model)
        self.actionGet_RV_data.triggered.connect(self.get_RV_data)

        self.actionGet_Orb_evol.triggered.connect(self.get_orb_evol)

        ############ Jupyter  #################

        self.actionpush_Jupyter_to_GUI.triggered.connect(self.get_jupyter_vars) 



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

        if qso_mode:
            self.actionvisit_TRIFON_on_GitHub.setText("&The QSO-Striker on GitHub")

        self.actionvisit_TRIFON_on_GitHub.triggered.connect(self.get_credits)
        self.actionCredits.triggered.connect(self.print_info_credits)
        self.actionConfidence_Intervals_Table.triggered.connect(lambda: self.print_chi_table())

        ############### Orb. Evol. plotting ####################
        self.comboBox_pl_1.activated.connect(self.plot_delta_omega)
        self.comboBox_pl_2.activated.connect(self.plot_delta_omega)
        self.radioButton_dom_180_fold.toggled.connect(self.plot_delta_omega)

        self.domega_esin_ecos.stateChanged.connect(self.plot_delta_omega)
        #self.domega_esin_ecos_e_in.toggled.connect(self.plot_delta_omega)

        self.buttonGroup_ei_eo_sincos.buttonClicked.connect(self.plot_delta_omega)


        self.Dom_plot_line.toggled.connect(self.plot_delta_omega)
        self.Dom_plot_dot.toggled.connect(self.plot_delta_omega)

        self.radioButton_theta_180_fold.toggled.connect(self.plot_theta)
        self.theta_esin_ecos.stateChanged.connect(self.plot_theta)
        self.theta_esin_ecos_e_in.toggled.connect(self.plot_theta)

        self.theta_plot_line.toggled.connect(self.plot_theta)
        self.theta_plot_dot.toggled.connect(self.plot_theta)

        self.per_evol_comboBox_pl_1.activated.connect(self.plot_per_rat)
        self.per_evol_comboBox_pl_2.activated.connect(self.plot_per_rat)

        self.Prat_plot_line.toggled.connect(self.plot_per_rat)
        self.Prat_plot_dot.toggled.connect(self.plot_per_rat)

        self.a_plot_dot.toggled.connect(self.plot_evol_a)
        self.a_plot_line.toggled.connect(self.plot_evol_a)
        self.a_line_size.valueChanged.connect(self.plot_evol_a)
        self.a_dot_size.valueChanged.connect(self.plot_evol_a)   


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
        self.extra_plot_Tra_phase_norm.stateChanged.connect(self.update_transit_plots)
        self.tra_model_width.valueChanged.connect(self.update_transit_plots)
        self.tra_model_z.valueChanged.connect(self.update_transit_plots)
        self.use_rich_tra_model.stateChanged.connect(self.update_transit_plots)
        self.tra_model_ndata_fact.valueChanged.connect(self.update_transit_plots)
        self.tra_label_fontsize.valueChanged.connect(self.update_transit_plots)
       
        self.trans_phase_slider.valueChanged.connect(self.update_transit_plots)       
        self.tra_xaxis_offset.valueChanged.connect(self.update_transit_plots)
        self.tra_plot_add_o_c.stateChanged.connect(self.update_transit_plots)                

        self.plot_transit_GP_model.stateChanged.connect(self.update_transit_plots) 

        self.tra_legend.stateChanged.connect(self.update_transit_plots) 
       

        ############### RV plotting controll ####################      
        self.rv_model_width.valueChanged.connect(self.update_RV_plots)
        self.rv_model_width.valueChanged.connect(self.update_extra_plots)    
        self.RV_model_z.valueChanged.connect(self.update_RV_plots)
        self.RV_model_z.valueChanged.connect(self.update_extra_plots)    
        self.RV_plot_add_o_c.stateChanged.connect(self.update_RV_plots)
        self.extra_plot_RVphase_norm.stateChanged.connect(self.update_extra_plots)    


        self.plot_RV_GP_model.stateChanged.connect(self.update_RV_plots)

        self.RV_legend.stateChanged.connect(self.update_RV_plots)


        ############### RV GLS plotting controll ####################
        self.gls_model_width.valueChanged.connect(self.update_RV_GLS_plots)
        self.gls_o_c_model_width.valueChanged.connect(self.update_RV_o_c_GLS_plots)

        self.N_GLS_peak_to_point.valueChanged.connect(self.update_RV_GLS_plots)
        self.N_GLS_peak_to_point.valueChanged.connect(self.update_RV_o_c_GLS_plots)
        self.avoid_GLS_RV_alias.stateChanged.connect(self.update_RV_GLS_plots)
        self.avoid_GLS_RV_alias.stateChanged.connect(self.update_RV_o_c_GLS_plots)

        self.alias_days_gls.valueChanged.connect(self.update_RV_GLS_plots)
        self.alias_days_gls.valueChanged.connect(self.update_RV_o_c_GLS_plots)
        self.show_alias_GLS.stateChanged.connect(self.update_RV_GLS_plots)
        self.show_alias_GLS.stateChanged.connect(self.update_RV_o_c_GLS_plots)

        self.gls_o_c_GP.stateChanged.connect(self.run_gls_o_c)
        self.gls_o_c_incl_jitter.stateChanged.connect(self.run_gls_o_c)
        self.gls_incl_jitter.stateChanged.connect(self.run_gls)


        self.N_window_peak_to_point.valueChanged.connect(self.update_WF_plots)

        self.N_MLP_peak_to_point.valueChanged.connect(self.update_RV_MLP_plots)
        self.avoid_MLP_RV_alias.stateChanged.connect(self.update_RV_MLP_plots)
        self.alias_days_mlp.valueChanged.connect(self.update_RV_MLP_plots)
        self.show_alias_MLP.stateChanged.connect(self.update_RV_MLP_plots)


        self.N_TLS_peak_to_point.valueChanged.connect(self.update_tls_plots)
        self.N_TLS_peak_to_point.valueChanged.connect(self.update_tls_o_c_plots)

        self.jitter_to_plots.stateChanged.connect(self.update_plots)
        self.split_jitter.stateChanged.connect(self.update_plots)
        self.jitter_color_button.clicked.connect(self.get_RV_jitter_color)



        self.buttonGroup_use_planets.buttonClicked.connect(self.update_veiw)


        self.init_correlations_combo()
        self.init_activity_combo()
        self.init_scipy_combo()
        self.init_gls_norm_combo()

        self.comboBox_scipy_minimizer_1.activated.connect(self.check_scipy_min)
        self.comboBox_scipy_minimizer_2.activated.connect(self.check_scipy_min)
        self.gls_norm_combo.activated.connect(self.update_plots)

        self.init_ns_samp_opt_combo()
        self.comboBox_ns_samp_opt.activated.connect(self.check_ns_samp_opt_combo)
        self.comboBox_ns_bound_opt.activated.connect(self.check_ns_samp_opt_combo)
        
        self.use_ns_maxiter.stateChanged.connect(self.check_ns_samp_opt_combo)
        self.use_ns_maxcall.stateChanged.connect(self.check_ns_samp_opt_combo)
        self.ns_use_stop.stateChanged.connect(self.check_ns_samp_opt_combo)

        self.setWindowIcon(QtGui.QIcon('./lib/UI/33_striker.png'))

        self.radioButton_act_GLS_period.toggled.connect(lambda: self.update_activity_gls_plots(self.comboBox_act_data_gls.currentIndex()))
       
        self.comboBox_act_data_gls.activated.connect(lambda: self.update_activity_gls_plots(self.comboBox_act_data_gls.currentIndex())) 
        self.comboBox_act_data.activated.connect(lambda: self.update_activity_data_plots(self.comboBox_act_data.currentIndex())) 
        self.act_data_size.valueChanged.connect(lambda: self.update_activity_data_plots(self.comboBox_act_data.currentIndex())) 
       
        self.comboBox_corr_1.activated.connect(self.update_correlations_data_plots) 
        self.comboBox_corr_2.activated.connect(self.update_correlations_data_plots) 
        self.plot_corr_err.stateChanged.connect(self.update_correlations_data_plots)
        self.plot_corr_coef.stateChanged.connect(self.update_correlations_data_plots)

        

        self.do_RV_GP.stateChanged.connect(self.set_use_RV_GP)
        self.do_tra_GP.stateChanged.connect(self.set_use_tra_GP)


        ############### Cross hair ####################

        self.gls_cross_hair.stateChanged.connect(self.update_RV_GLS_plots)
        self.gls_o_c_cross_hair.stateChanged.connect(self.update_RV_o_c_GLS_plots)
        self.RV_plot_cross_hair.stateChanged.connect(self.update_RV_plots)

        self.RV_o_c_plot_cross_hair.stateChanged.connect(self.update_RV_plots)
        self.trans_plot_cross_hair.stateChanged.connect(self.update_transit_plots)
        self.trans_o_c_plot_cross_hair.stateChanged.connect(self.update_transit_plots)
        self.tls_cross_hair.stateChanged.connect(self.update_tls_plots)
        self.tls_o_c_cross_hair.stateChanged.connect(self.update_tls_o_c_plots)
        self.gls_act_cross_hair.stateChanged.connect(lambda: self.update_activity_gls_plots(self.comboBox_act_data_gls.currentIndex()))
        self.mlp_cross_hair.stateChanged.connect(self.update_RV_MLP_plots)


        self.extra_plot_cross_hair.stateChanged.connect(self.update_extra_plots)
        self.extra_plot_cross_hair_gls.stateChanged.connect(self.update_extra_plots2)
        self.inpector_plot_cross_hair.stateChanged.connect(self.init_plot_data_inspect)
        #self.inpector_plot_cross_hair.stateChanged.connect(lambda: self.tree_view_tab.listview.connect(self.plot_data_inspect))



        self.color_corr.clicked.connect(self.get_corr_color)
        self.corr_x_label.clicked.connect(self.corr_plot_x_labels)
        self.corr_y_label.clicked.connect(self.corr_plot_y_labels)

        self.colors_gls.clicked.connect(self.get_RV_GLS_plot_color)
        self.colors_gls_o_c.clicked.connect(self.get_RV_o_c_GLS_plot_color)
        self.colors_alias_gls.clicked.connect(self.get_RV_GLS_alias_color)
        self.colors_alias_mlp.clicked.connect(self.get_RV_MLP_alias_color)

        #self.colors_ttv.clicked.connect(self.get_ttv_plot_color)
        #self.ttv_model_color.clicked.connect(self.get_ttv_model_color)


        self.color_delta_om.clicked.connect(self.get_delta_omega_color)
        self.delta_om_x_label.clicked.connect(self.delta_omega_plot_x_labels)
        self.delta_om_y_label.clicked.connect(self.delta_omega_plot_y_labels)

        self.color_per_evol.clicked.connect(self.get_per_rat_color)
        self.per_evol_x_label.clicked.connect(self.per_rat_plot_x_labels)
        self.per_evol_y_label.clicked.connect(self.per_rat_plot_y_labels)

        self.color_theta.clicked.connect(self.get_theta_color)
        self.theta_x_label.clicked.connect(self.theta_plot_x_labels)
        self.theta_y_label.clicked.connect(self.theta_plot_y_labels)


        #self.tab_timeseries_RV.currentChanged.connect(self.tab_selected)

        self.radioButton_RV_o_c_GLS_period.toggled.connect(self.update_RV_o_c_GLS_plots)
        self.radioButton_RV_GLS_period.toggled.connect(self.update_RV_GLS_plots)
        self.radioButton_RV_MLP_period.toggled.connect(self.update_RV_MLP_plots)

        self.mute_boxes()
#        self.radioButton_transit_RV.toggled.connect(self.mute_boxes)
#        self.radioButton_transit.toggled.connect(self.mute_boxes)
#        self.radioButton_RV.toggled.connect(self.mute_boxes)
#        self.radioButton_ttv.toggled.connect(self.mute_boxes)
#        self.radioButton_ttv_RV.toggled.connect(self.mute_boxes)



        self.customize_mcmc_cornerplot.clicked.connect(lambda: self.get_cornerplot_param(type_plot = "mcmc"))
        self.customize_ns_cornerplot.clicked.connect(lambda: self.get_cornerplot_param(type_plot = "nest"))


        self.radioButton_rho.toggled.connect(self.set_aR)
        self.radioButton_ewm.toggled.connect(self.set_hkl)
#        self.radioButton_hkl.toggled.connect(self.set_hkl)
        self.radioButton_KP.toggled.connect(self.set_kp_ma)


        self.radioButton_RV_WF_period.toggled.connect(self.update_WF_plots)

        self.calc_TLS.clicked.connect(self.worker_tls)
        self.calc_TLS_o_c.clicked.connect(lambda: self.worker_tls(resid =True))
        self.calc_MLP.clicked.connect(self.worker_mlp)


        self.actionQuit.triggered.connect(self.close) 

        self.actionopen_RVmod_init_file.triggered.connect(self.showDialog_fortran_input_file)
        self.actionOpen_RVbank_file.triggered.connect(self.showDialog_RVbank_input_file)
        self.actionOpen_MAROONX_file.triggered.connect(self.showDialog_MAROONX_input_file)


        self.jupiter_push_vars()

        self.update_color_picker()

        self.buttonGroup_color_picker.buttonClicked.connect(self.get_color)
        self.update_color_picker_tra()
        self.buttonGroup_color_picker_tra.buttonClicked.connect(self.get_color_tra)
        self.update_color_picker_ttv()
        self.buttonGroup_color_picker_ttv.buttonClicked.connect(self.get_color_ttv)
        self.update_color_picker_ast()
        self.buttonGroup_color_picker_ast.buttonClicked.connect(self.get_color_ast)

        self.tess_pdc_dialog = pdc(self)
        self.dialog_symbols = show_symbols(self)

        self.buttonGroup_symbol_picker.buttonClicked.connect(self.get_symbol) 
        self.buttonGroup_symbol_picker_tra.buttonClicked.connect(self.get_symbol_tra)
        self.buttonGroup_symbol_picker_ttv.buttonClicked.connect(self.get_symbol_ttv)
        self.buttonGroup_symbol_picker_ast.buttonClicked.connect(self.get_symbol_ast)        
        
        self.lables_cornerplot = []
        self.dialog_select_param_cornerplot = show_param_boxes(self)

        ###########  GP control ##########

#        self.set_RV_GP()
 #       self.set_tra_GP()

        self.buttonGroup_use_RV_GP_kernel.buttonClicked.connect(self.set_RV_GP)
        self.buttonGroup_use_tra_GP_kernel.buttonClicked.connect(self.set_tra_GP)
#        self.buttonGroup_use_GP.buttonClicked.connect(self.set_use_GP)

        self.buttonGroup_link_tra_GP_to_RV_GP.buttonClicked.connect(self.set_link_GP)
        #### Transit detrend   ####

        self.buttonGroup_detrend_tra.buttonClicked.connect(self.transit_data_detrend)
        self.DetrendWindow = DetrendWindow(self)

        #### RVS detrend   ####
        self.buttonGroup_options_rvs.buttonClicked.connect(self.rvs_data_options)
        self.RvsWindow = RvsWindow(self)

        #### Activity detrend   ####
        self.buttonGroup_options_act.buttonClicked.connect(self.activity_data_options)
        self.ActivityWindow = ActivityWindow(self)

        ####### LD models #############
 
        self.buttonGroup_use_ld_1.buttonClicked.connect(self.set_tra_ld)
        self.buttonGroup_use_ld_2.buttonClicked.connect(self.set_tra_ld)
        self.buttonGroup_use_ld_3.buttonClicked.connect(self.set_tra_ld)
        self.buttonGroup_use_ld_4.buttonClicked.connect(self.set_tra_ld)
        self.buttonGroup_use_ld_5.buttonClicked.connect(self.set_tra_ld)
        self.buttonGroup_use_ld_6.buttonClicked.connect(self.set_tra_ld)
        self.buttonGroup_use_ld_7.buttonClicked.connect(self.set_tra_ld)
        self.buttonGroup_use_ld_8.buttonClicked.connect(self.set_tra_ld)
        self.buttonGroup_use_ld_9.buttonClicked.connect(self.set_tra_ld)
        self.buttonGroup_use_ld_10.buttonClicked.connect(self.set_tra_ld)
        self.buttonGroup_use_ld_11.buttonClicked.connect(self.set_tra_ld)
        self.buttonGroup_use_ld_12.buttonClicked.connect(self.set_tra_ld)
        self.buttonGroup_use_ld_13.buttonClicked.connect(self.set_tra_ld)
        self.buttonGroup_use_ld_14.buttonClicked.connect(self.set_tra_ld)
        self.buttonGroup_use_ld_15.buttonClicked.connect(self.set_tra_ld)
        self.buttonGroup_use_ld_16.buttonClicked.connect(self.set_tra_ld)
        self.buttonGroup_use_ld_17.buttonClicked.connect(self.set_tra_ld)
        self.buttonGroup_use_ld_18.buttonClicked.connect(self.set_tra_ld)
        self.buttonGroup_use_ld_19.buttonClicked.connect(self.set_tra_ld)
        self.buttonGroup_use_ld_20.buttonClicked.connect(self.set_tra_ld)


        self.Button_apply_set_ld_group.clicked.connect(self.set_tra_gr)

        ####### transit linear models #############

        self.buttonGroup_use_tra_reg_1.buttonClicked.connect(lambda: self.set_tra_reg(ind = 0))
        self.buttonGroup_use_tra_reg_2.buttonClicked.connect(lambda: self.set_tra_reg(ind = 1))
        self.buttonGroup_use_tra_reg_3.buttonClicked.connect(lambda: self.set_tra_reg(ind = 2))
        self.buttonGroup_use_tra_reg_4.buttonClicked.connect(lambda: self.set_tra_reg(ind = 3))
        self.buttonGroup_use_tra_reg_5.buttonClicked.connect(lambda: self.set_tra_reg(ind = 4))
        self.buttonGroup_use_tra_reg_6.buttonClicked.connect(lambda: self.set_tra_reg(ind = 5))
        self.buttonGroup_use_tra_reg_7.buttonClicked.connect(lambda: self.set_tra_reg(ind = 6))
        self.buttonGroup_use_tra_reg_8.buttonClicked.connect(lambda: self.set_tra_reg(ind = 7))
        self.buttonGroup_use_tra_reg_9.buttonClicked.connect(lambda: self.set_tra_reg(ind = 8))
        self.buttonGroup_use_tra_reg_10.buttonClicked.connect(lambda: self.set_tra_reg(ind = 9))
        self.buttonGroup_use_tra_reg_11.buttonClicked.connect(lambda: self.set_tra_reg(ind = 10))
        self.buttonGroup_use_tra_reg_12.buttonClicked.connect(lambda: self.set_tra_reg(ind = 11))
        self.buttonGroup_use_tra_reg_13.buttonClicked.connect(lambda: self.set_tra_reg(ind = 12))
        self.buttonGroup_use_tra_reg_14.buttonClicked.connect(lambda: self.set_tra_reg(ind = 13))
        self.buttonGroup_use_tra_reg_15.buttonClicked.connect(lambda: self.set_tra_reg(ind = 14))
        self.buttonGroup_use_tra_reg_16.buttonClicked.connect(lambda: self.set_tra_reg(ind = 15))
        self.buttonGroup_use_tra_reg_17.buttonClicked.connect(lambda: self.set_tra_reg(ind = 16))
        self.buttonGroup_use_tra_reg_18.buttonClicked.connect(lambda: self.set_tra_reg(ind = 17))
        self.buttonGroup_use_tra_reg_19.buttonClicked.connect(lambda: self.set_tra_reg(ind = 18))
        self.buttonGroup_use_tra_reg_20.buttonClicked.connect(lambda: self.set_tra_reg(ind = 19))


 
       # self.RV_phase_slider.sliderReleased.connect(self.rv_plot_phase_change)
        self.RV_phase_slider.valueChanged.connect(self.rv_plot_phase_change)


        self.threadpool = QtCore.QThreadPool()
        self.threadpool.setMaxThreadCount(int(cpu_count()/2))


        self.select_session(-1)


        self.update_settings()
        self.time_step_model.valueChanged.connect(self.check_settings)
        self.dyn_model_accuracy.valueChanged.connect(self.check_settings)
        self.dyn_model_to_kill.valueChanged.connect(self.check_settings)
        self.kep_model_to_kill.valueChanged.connect(self.check_settings)
        self.points_to_draw_model.valueChanged.connect(self.check_settings)
        


        
        
        self.mute_boxes()
        self.update_RV_jitter_flag()
        self.check_cornerplot_samples()
        
        self.force_copl_incl.stateChanged.connect(self.set_force_copl_incl)
        self.force_sameRV_jit.stateChanged.connect(self.set_force_sameRV_jit)        
 
#        self.set_use_GP()
#        self.set_gui_use_GP()

        self.update_GUI_St_params()
        #self.update_St_params()

        ############### Stellar params ####################      
        self.St_mass_input.valueChanged.connect(lambda: self.update_St_params(ind=1))
        self.St_radius_input.valueChanged.connect(lambda: self.update_St_params(ind=3))
        self.St_lumin_input.valueChanged.connect(lambda: self.update_St_params(ind=5))
        self.St_teff_input.valueChanged.connect(lambda: self.update_St_params(ind=7))
        self.St_vsini_input.valueChanged.connect(lambda: self.update_St_params(ind=9))

        self.err_St_mass_input.valueChanged.connect(lambda: self.update_St_params(ind=2))
        self.err_St_radius_input.valueChanged.connect(lambda: self.update_St_params(ind=4))
        self.err_St_lumin_input.valueChanged.connect(lambda: self.update_St_params(ind=6))
        self.err_St_teff_input.valueChanged.connect(lambda: self.update_St_params(ind=8))
        self.err_St_vsini_input.valueChanged.connect(lambda: self.update_St_params(ind=10))


        #self.plot_opt_tab.tabBarClicked.connect(self.check_cornerplot_samples)
        #self.cornerplot_plot_tab.isVisible.connect(self.check_cornerplot_samples)
        self.tabWidget_3.tabBarDoubleClicked.connect(self.change_nbody_plot_opt_tab)        
        self.tabWidget_6.tabBarDoubleClicked.connect(self.change_nbody_plot_opt_tab_res)        

        self.count_cpus()
        rv.check_swift()


        self.stop_button.clicked.connect(self.stop_nest)

        self.data_str = "RV data"     
        self.QSO_str = "Exo" 
        if qso_mode==True:
            self.set_qso_mode = qso.qso_mode(self)
            self.data_str = "Data"
            self.QSO_str = "QSO"

        print("""You are running a development version of the %s-Striker (ver. %s). 
              
This version is almost full, but there are still some parts of the tool, which are in a 'Work in progress' state. Please, 'git pull' regularly to be up to date with the newest version.
"""%(self.QSO_str,self.es_version))

        if sys.version_info[0] == 2:
            print("""
It seems that you started the '%s-Striker' with Python 2. Please consider Python 3 for your work with the '%s-Striker'.
"""%(self.QSO_str,self.QSO_str))


        print("""Here you can get some more information from the tool's workflow, stdout/strerr, and piped results.""")
       
       
        #self.use_K1.setStyleSheet("color: red")

        #f = open("loooog","w")
        #for j in vars(self):
        #    f.write(j)
        #    f.write('\n')
        #f.close()  
       # print(vars(self), sep='\n')
 
 
  

#class CustomStyle(QtWidgets.QProxyStyle):
#    def drawControl(self, element, option, painter, widget=None):
#        if element in [self.CE_PushButtonLabel, self.CE_RadioButtonLabel] and not option.state & self.State_Enabled:
#            color = QtCore.Qt.red
#            option.palette.setColor(option.palette.ColorGroup.Disabled, option.palette.ColorRole.Text, color)
#        super().drawControl(element, option, painter, widget)


class LoadingScreen(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
 
        self.setFixedSize(1300*0.7, 1250*0.7)

        self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint | QtCore.Qt.CustomizeWindowHint)
        self.label_animation = QtWidgets.QLabel(self)

        self.label_animation.setSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Ignored)
        # self.label_animation.setStyleSheet("""{ background-image: url('./Icons/Pulse-1s-100px.gif');}""")
        self.movie= QtGui.QMovie('./lib/UI/33_striker.png')
        self.label_animation.setMovie(self.movie)
        # self.label_animation.setPixmap(QtGui.QPixmap("./Icons/Pulse-1s-100px.gif"))
 
        timer = QtCore.QTimer
        
        self.startAnimation()
        timer.singleShot(3000, self.stopAniamtion)
        self.show()       
        

    def startAnimation(self):
        self.movie.start()


    def stopAniamtion(self):
        self.movie.stop()
        self.close()




def main():

    app = QtWidgets.QApplication(sys.argv)
    
#    print(app)
    app.setStyle('Fusion') #The available styles depend on your platform but are usually 'Fusion', 'Windows', 'WindowsVista' (Windows only) and 'Macintosh' (Mac only). 
#    app.setStyle(CustomStyle())

    window = Exo_striker()


    screen_resolution = app.primaryScreen().geometry()
    width, height = screen_resolution.width(), screen_resolution.height()
    #print(width, height)
    if height < 920:
        window.setMinimumWidth(int(width*0.6))
        window.setMinimumHeight(int(height*0.6))
        window.resize(int(width*0.8), int(height*0.8))
    else:
        pass
    
#    window.installEventFilter(window)



    window.show()

    sys.exit(app.exec_())
#    sys.exit(app.exec())
if __name__ == '__main__':
    main() 


