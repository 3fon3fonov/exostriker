#!/usr/bin/python

##
######## this is the file where all experiments are happening! 
## the code below is still work in progress. 


from __future__ import print_function
__author__ = 'Trifon Trifonov, Jakub Morawski'

import sys, os
#sys.path.insert(0, '../lib')
sys.path.append('./lib/RV_mod/')
#import gls as gls 
#import prior_functions as pr
import emcee


#import prior_functions as pr
#import rot_kernels
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('SVG') 
 

import time
#import multiprocessing
from pathos import multiprocessing
#from pathos.multiprocessing import ProcessingPool as Pool

#from emcee.utils import MPIPool
import celerite 
from celerite import terms
import dynesty

#import batman

#import copy
import dill
import scipy.optimize as op
from scipy import stats


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




from CustomSampler import CustomSampler
from Warning_log import Warning_log
from parameters import parameters, GP_parameters, parameter_errors, parameter_bounds, use_flags
from fortran_output import fortran_output
from functions import *
from errors import Error, InputError, FittingError
from kernel import kernel, rvmodel, summary
from rv_files import rvfile, rvfile_list

TAU= 2.0*np.pi 
DEFAULT_STELLAR_MASS=1.0
DEFAULT_JITTER=1.0
NPLMAX=20
DEFAULT_PATH='./datafiles/'


import GP_kernels


 

def initiategps(obj,  kernel_id=-1): 
    
    # Prepare objects for Gaussian Processes        
    
    if len(obj.GP_rot_params) != 0:
        obj.params.update_GP_params(obj.GP_rot_params, kernel_id=kernel_id)
        obj.use.update_use_GP_params(obj.GP_rot_use)
    
    #print(obj.gp_kernel)
    if obj.gp_kernel == 'RotKernel':
        kernel = terms.TermSum(GP_kernels.RotationTerm(
                log_amp=np.log(obj.GP_rot_params[0]),
                log_timescale=np.log(obj.GP_rot_params[1]),
                log_period=np.log(obj.GP_rot_params[2]),
                log_factor=np.log(obj.GP_rot_params[3])))
        
    elif obj.gp_kernel == 'SHOKernel':  
        
        kernel = terms.SHOTerm(log_S0=np.log(obj.GP_sho_params[0]), 
                               log_Q=np.log(obj.GP_sho_params[1]), 
                               log_omega0=np.log(obj.GP_sho_params[2]))
    
    gps = celerite.GP(kernel, mean=0.0)
    gps.compute(obj.filelist.time, obj.filelist.rv_err)
    
    obj.gps = gps          
    return
 

def plot_gp(obj, curve=False):
    
    import matplotlib.pyplot as plt

    color="#ff7f0e"
    colors = ['b','g','r']
    
    x     = obj.fit_results.rv_model.jd
    y     = obj.fit_results.rv_model.o_c
    y_err = obj.fit_results.rv_model.rv_err
    idset = obj.filelist.idset
    
    
    if curve==True:
        x_model = np.linspace(min(x), max(x), 5000) #obj.fit_results.model_jd
        mu,var,std = obj.gp_model_curve
        
    else:      
        x_model = x
        mu,var,std = obj.gp_model_data    
    
    #print(mu[0:10])
    #print(y[0:10])
    
    for i in range(obj.filelist.ndset):
        plt.errorbar(x[idset==i],y[idset==i], yerr=y_err[idset==i], fmt=".",color=colors[i],  capsize=0); 
    
    plt.plot(x_model, mu, color = '0.5' );
    plt.fill_between(x_model ,mu+std,  mu-std, color=color, alpha=0.3, edgecolor="none")


def get_gps_model(obj,  kernel_id=-1): 
    
    initiategps(obj,  kernel_id=-1)
    #gp_model_data  = []
    
    ############ DATA ####################
    #for i in range(obj.filelist.ndset):
        #gp.set_parameter_vector(
         
    y = obj.fit_results.rv_model.o_c
    x = obj.fit_results.rv_model.jd
    mu, var = obj.gps.predict(y, x, return_var=True)
    std = np.sqrt(var)

    obj.gp_model_data = [mu,var,std]
        
    ############ MODEL ####################

   # kernel=[]
   # gps=[]
    x = obj.fit_results.model_jd
    #x= np.linspace(min(x2), max(x2), 5000)
    #y = obj.fit_results.model

   # kernel = obj.params.GP_params.rot_kernel
   # gps = celerite.GP(kernel, mean=0.0)
   # gps.compute(x,[0]*len(x))
    #gps.compute(obj.filelist.time, obj.filelist.rv_err)
 
 
    mu, var = obj.gps.predict(y, x, return_var=True)
    std = np.sqrt(var)

    obj.gp_model_curve = [mu,var,std]        

    return



######### Transit GP work in progress ###########    

def initiate_tansit_gps(obj,  kernel_id=-1): 
    
    # Prepare objects for Gaussian Processes        
    
    #print(obj.gp_kernel)
    if obj.tra_gp_kernel == 'RotKernel':
        tra_kernel = terms.TermSum(GP_kernels.RotationTerm(
                log_amp=np.log(obj.tra_GP_rot_params[0]),
                log_timescale=np.log(obj.tra_GP_rot_params[1]),
                log_period=np.log(obj.tra_GP_rot_params[2]),
                log_factor=np.log(obj.tra_GP_rot_params[3])))
        
    elif obj.tra_gp_kernel == 'SHOKernel':   
        tra_kernel = terms.SHOTerm(log_S0=np.log(obj.tra_GP_sho_params[0]), 
                               log_Q=np.log(obj.tra_GP_sho_params[1]), 
                               log_omega0=np.log(obj.tra_GP_sho_params[2]))
    
    tra_gps = celerite.GP(tra_kernel, mean=1.0)
    tra_gps.compute(obj.tra_data_sets[0][0], obj.tra_data_sets[0][2])
    
    obj.tra_gps = tra_gps          
    return    


def plot_transit_gp(obj, curve=False):
    
    import matplotlib.pyplot as plt

    color="#ff7f0e"
    colors = ['b','g','r']
    
    x     = obj.tra_data_sets[0][0]
    y     = obj.tra_data_sets[0][1]
    y_err = obj.tra_data_sets[0][2]
    #idset = obj.filelist.idset
    
    
    if curve==True:
        x_model = np.linspace(min(x), max(x), 5000) #obj.fit_results.model_jd
        mu,var,std = obj.tra_gp_model_curve
        
    else:      
        x_model = x
        mu,var,std = obj.tra_gp_model_data    
    
    #print(mu[0:10])
    #print(y[0:10])
    
    #for i in range(obj.filelist.ndset):
        #plt.errorbar(x[idset==i],y[idset==i], yerr=y_err[idset==i], fmt=".",color=colors[i],  capsize=0); 
    plt.errorbar(x,y, yerr=y_err, fmt=".",color=colors[0],  capsize=0); 
   
    plt.plot(x_model, mu, color = '0.5' );
    plt.fill_between(x_model ,mu+std,  mu-std, color=color, alpha=0.3, edgecolor="none")






def get_transit_ts(obj,  kernel_id=-1): 
                         
    tr_files = []
    
    for i in range(10):
        if len(obj.tra_data_sets[i]) != 0:
            tr_files.append(obj.tra_data_sets[i])
    
    for j in range(len(tr_files)):        
    
    #if len(fit.tra_data_sets[0]) != 0:
        t = np.array(tr_files[j][0])
        flux = np.array(tr_files[j][1] + obj.tra_off[j])
        flux_err = np.sqrt(tr_files[j][2]**2 + obj.tra_jitt[j]**2)
        
        
        
        obj.prepare_for_mcmc(rtg = obj.rtg)    
        par = np.array(obj.parameters)  

        flux_model = np.ones(len(flux))
        m =  {k: [] for k in range(9)}
         
        
        #### a quick fix, TBD! ########
        if obj.rtg[1]:
            if obj.gp_kernel == 'RotKernel':
                rv_gp_npar = 4
            if obj.gp_kernel == 'SHOKernel':
                rv_gp_npar = 3
            #fit.gps = []
        else:
            rv_gp_npar = 0   
            


        obj.tr_params.limb_dark = str(obj.ld_m[j])      #limb darkening model       
        #print(tr_model[0][j], tr_model[1][j] )
        obj.tr_params.u = obj.ld_u[j]
        
        
        for i in range(obj.npl):


            if obj.hkl == True:
                obj.tr_params.ecc = np.sqrt(par[obj.filelist.ndset*2 +7*i+2]**2 + par[obj.filelist.ndset*2 +7*i+3]**2)
                obj.tr_params.w  = np.degrees(np.arctan2(par[obj.filelist.ndset*2 +7*i+2],par[obj.filelist.ndset*2 +7*i+3]))%360
            else:
                obj.tr_params.ecc = par[obj.filelist.ndset*2 +7*i+2] #0.0  
                obj.tr_params.w   = par[obj.filelist.ndset*2 +7*i+3] #90
            
            obj.tr_params.per = par[obj.filelist.ndset*2 +7*i+1] #1.0    #orbital period
            obj.tr_params.inc = par[obj.filelist.ndset*2 +7*i+5]#90. #orbital inclination (in degrees)
                
            obj.tr_params.t0  = par[obj.filelist.ndset*2  +7*obj.npl +2+rv_gp_npar + 3*i]                
            obj.tr_params.a   = par[obj.filelist.ndset*2  +7*obj.npl +2+rv_gp_npar + 3*i+1] #15  #semi-major axis (in units of stellar radii)
            obj.tr_params.rp  = par[obj.filelist.ndset*2  +7*obj.npl +2+rv_gp_npar + 3*i+2] #0.15   #planet radius (in units of stellar radii)
            #print(tr_params.t0)
            #print(tr_params.per, tr_params.ecc,tr_params.w, tr_params.inc, tr_params.t0,tr_params.a,tr_params.rp )
    
            m[i] = batman.TransitModel(obj.tr_params, t)    #initializes model
 
            flux_model = flux_model * m[i].light_curve(obj.tr_params)     
            
 
            
        obj.tra_data_sets[0][4] = flux - flux_model   
        
        
    return  




def get_transit_gps_model(obj,  kernel_id=-1): 
  
    get_transit_ts(obj)      
    initiate_tansit_gps(obj,  kernel_id=-1)
    #gp_model_data  = []

    ############ DATA ####################
    #for i in range(obj.filelist.ndset):
        #gp.set_parameter_vector(
         
    y = obj.tra_data_sets[0][4] #obj.fit_results.rv_model.o_c
    x = obj.tra_data_sets[0][0]
    
    GP_var = False
    if GP_var == True:
        mu, var = obj.tra_gps.predict(y, x, return_var=True) 
        std = np.sqrt(var)
    
        obj.tra_gp_model_data = [mu,var,std]
            
        ############ MODEL ####################
      
        mu, var = obj.tra_gps.predict(y, x, return_var=True)
        std = np.sqrt(var)
    
        obj.tra_gp_model_curve = [mu,var,std]        
        
    else:
        
        mu = obj.tra_gps.predict(y, x, return_cov=False)
        #std = np.sqrt(var)
        obj.tra_gp_model_data = [mu,np.zeros(len(mu)),np.zeros(len(mu))]
            
        ############ MODEL ####################
      
        mu = obj.tra_gps.predict(y, x, return_cov=False)
       # std = np.sqrt(var)
    
        obj.tra_gp_model_curve= [mu,np.zeros(len(mu)),np.zeros(len(mu))]        

    return
 
######### transit GP work in progress ###########   

def transit_loglik(tr_files,vel_files,tr_params,tr_model,par,rv_gp_npar,npl,hkl, rtg,tra_gps  ):
 
    gp_tr_loglik = 0

    for j in range(len(tr_files)):
                 
        t = tr_files[j][0] 
        flux = tr_files[j][1] + par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + len(tr_files)*j]
        #flux_err = np.sqrt(tr_files[j][2]**2 + par[len(vel_files)*2 +7*npl +5 + 3*npl + len(tr_files)*j +1]**2)
        flux_err =  tr_files[j][2] 
        
        flux_model =[1]*len(flux)
        
        m =  {k: [] for k in range(9)}
 
        tr_params.limb_dark = str(tr_model[0][j])      #limb darkening model       
        #print(tr_model[0][j], tr_model[1][j] )
        tr_params.u = tr_model[1][j]   

      #  if str(tr_model[0][j]) == "uniform":    
      #      tr_params.u = tr_model[1][j]   
       # elif  str(tr_model[j]) == "linear":
      #      tr_params.u = [0.1]                  
       # elif  str(tr_model[j]) == "quadratic":
      #      tr_params.u = [0.1,0.3]               
      #  elif  str(tr_model[j]) == "nonlinear":
      #      tr_params.u = [0.5,0.1,0.1,-0.1]   
                 
        for i in range(npl):

            if hkl == True:
                tr_params.ecc = np.sqrt(par[len(vel_files)*2 +7*i+2]**2 + par[len(vel_files)*2 +7*i+3]**2)
                tr_params.w  = np.degrees(np.arctan2(par[len(vel_files)*2 +7*i+2],par[len(vel_files)*2 +7*i+3]))%360
            else:
                tr_params.ecc = par[len(vel_files)*2 +7*i+2] #0.0  
                tr_params.w   = par[len(vel_files)*2 +7*i+3] #90.0   #longitude of periastron (in degrees)               
            
            tr_params.per = par[len(vel_files)*2 +7*i+1] #1.0    #orbital period
            tr_params.inc = par[len(vel_files)*2 +7*i+5]#90. #orbital inclination (in degrees)                

            tr_params.t0  = par[len(vel_files)*2 +7*npl +2 +rv_gp_npar + 3*i]
            tr_params.a   = par[len(vel_files)*2 +7*npl +2 +rv_gp_npar + 3*i+1] #15  #semi-major axis (in units of stellar radii)
            tr_params.rp  = par[len(vel_files)*2 +7*npl +2 +rv_gp_npar + 3*i+2] #0.15   #planet radius (in units of stellar radii)

            m[i] = batman.TransitModel(tr_params, t)    #initializes model
 
            flux_model = flux_model * m[i].light_curve(tr_params)    


        if rtg[3] == False:
            sig2i = 1.0 / (flux_err**2 + par[len(vel_files)*2 +7*npl +2 +rv_gp_npar + 3*npl + len(tr_files)*j+1]**2 )
            #sig2i = sig2i/len(flux_err)
            tr_loglik = -0.5*(np.sum((flux -flux_model)**2 * sig2i - np.log(sig2i / 2./ np.pi))) # - np.log(sig2i / 2./ np.pi)
            #tr_loglik = -0.5*(np.sum((flux -flux_model)**2 * sig2i - np.log(sig2i))) # - np.log(sig2i / 2./ np.pi)
        else:
            
            param_vect = []
            for k in range(len(tra_gps.get_parameter_vector())):
                param_vect.append(np.log(par[len(vel_files)*2  +7*npl  + rv_gp_npar  + 3*npl + len(tr_files)*2 + 2 + k ]))
                #print(par[len(vel_files)*2  +7*npl  + rv_gp_npar  + 3*npl + len(tr_files)*2 + rv_gp_npar + 1 + k  ])
            tra_gps.set_parameter_vector(np.array(param_vect))
        
            tra_gp_pred = tra_gps.predict(flux -flux_model, t, return_cov=False)
            o_c_tra = (flux -flux_model) - tra_gp_pred
            
            sig2i = 1.0 / (flux_err**2 + par[len(vel_files)*2 +7*npl + 2 +rv_gp_npar + 3*npl + len(tr_files)*j+1]**2 )
            gp_tr_loglik = -0.5*(np.sum((o_c_tra)**2 * sig2i - np.log(sig2i / 2./ np.pi))) # - np.log(sig2i / 2./ np.pi)
                                             
    
            tr_loglik =  gp_tr_loglik          
   # print(tr_loglik)        
    return tr_loglik    
 
    


def model_loglik(p, program, par, flags, npl, vel_files, tr_files, tr_model, tr_params, epoch, stmass, gps, tra_gps, rtg, mix_fit, opt, outputfiles = [1,0,0], amoeba_starts=0, prior=0, eps='1.0E-8',dt=864000, when_to_kill=3000, npoints=50, model_max = 100, model_min =0): # generate input string for the fortran code, optionally as a file

    rv_loglik = 0
    gp_rv_loglik = 0
    tr_loglik = 0
    gp_tr_loglik = 0
   
    dt = opt["dt"] 
    eps = opt["eps"]   
    when_to_kill = opt["when_to_kill"] 
    copl_incl = opt["copl_incl"]
    hkl = opt["hkl"]
    cwd = opt["cwd"]
    gr_flag = opt["gr_flag"]
    
    if np.isnan(p).any():
        return -np.inf
    
    for j in range(len(flags)):
        par[flags[j]] = p[j]  
    
    if (rtg[1]):
        outputfiles = [1,1,0]
        rv_gp_npar = len(gps.get_parameter_vector())
    else:
        rv_gp_npar = 0   
 
    if (rtg[3]):
        outputfiles = [1,1,0]
        tra_gp_npar = len(tra_gps.get_parameter_vector())
    else:
        tra_gp_npar = 0          
         


    if rtg[2] == True:
        for i in range(npl): # (per, ecc, om, t_transit, epoch):
            
            if hkl == True:
                ecc_ = np.sqrt(par[len(vel_files)*2 +7*i+2]**2 + par[len(vel_files)*2 +7*i+3]**2)
                om_  = np.degrees(np.arctan2(par[len(vel_files)*2 +7*i+2],par[len(vel_files)*2 +7*i+3]))%360
            else:
                ecc_, om_, = par[len(vel_files)*2 +7*i+2],par[len(vel_files)*2 +7*i+3]
            
            
            par[len(vel_files)*2 +7*i+4] = ma_from_t0(par[len(vel_files)*2 +7*i+1],
                                                      ecc_, om_, par[len(vel_files)*2 +7*npl + 2 +rv_gp_npar + 3*i],epoch)
     
    else:
        for i in range(npl): # (per, ecc, om, ma, epoch):
            if hkl == True:
                ecc_ = np.sqrt(par[len(vel_files)*2 +7*i+2]**2 + par[len(vel_files)*2 +7*i+3]**2)
                om_  = np.degrees(np.arctan2(par[len(vel_files)*2 +7*i+2],par[len(vel_files)*2 +7*i+3]))%360
                Ma_  = (par[len(vel_files)*2 +7*i+4] - om_)%360.0  
            else:
                ecc_, om_, Ma_ = par[len(vel_files)*2 +7*i+2], par[len(vel_files)*2 +7*i+3], par[len(vel_files)*2 +7*i+4]               
            
            par[len(vel_files)*2 +7*npl + 2 +rv_gp_npar + 3*i] = transit_tperi(par[len(vel_files)*2 +7*i+1],
                                                                  ecc_, om_, Ma_ ,epoch)[1]%par[len(vel_files)*2 +7*i+1]
    
    if(rtg[0]):
        
        ppp= '%s << EOF\n%s %f %d %d %d %d %d %d\n%f %d %d %d \n%d\n'%(program, eps,dt,amoeba_starts,when_to_kill,npoints, model_max, model_min, gr_flag, stmass, outputfiles[0], outputfiles[1],outputfiles[2], len(vel_files)) # first three lines of fortran input: precision and timestep for integration, stellar mass and number of datasets
        for i in range(len(vel_files)): 
            # path for each dataset      
            ppp+='%s\n'%(vel_files[i])    
            # offset and jitter information for each dataset
            ppp+='%f\n%d\n'%(par[i],0)
            if (rtg[1]): 
                ppp+='%f\n%d\n'%(0,0)           
            else:
                ppp+='%f\n%d\n'%(par[i + len(vel_files)],0)
         
        # if mixed fitting is requested    
        ppp+='%d\n'%npl
        if mix_fit[0] == True and program == '%s/lib/fr/loglik_dyn+'%cwd:

            for i in range(npl): 
                ppp+='%d\n'%mix_fit[1][i]                
        
        for i in range(npl): # K,P,e,w,M,i,cap0m for each planet, and information which ones we use
            ppp+='%f %f %f %f %f %f %f %f\n'%(par[len(vel_files)*2 + 7*i],
                                               par[len(vel_files)*2 +7*i+1],
                                               par[len(vel_files)*2 +7*i+2],
                                               par[len(vel_files)*2 +7*i+3],
                                               par[len(vel_files)*2 +7*i+4],                                               
                                               par[len(vel_files)*2 +7*i+5],
                                               par[len(vel_files)*2 +7*i+6],
                                               #0
                                               #)
                                               par[len(vel_files)*2  +7*npl  + rv_gp_npar  + 3*npl + len(tr_files)*2 + tra_gp_npar + 2 + i ])
            #print(par[len(vel_files)*2  +7*npl  + rv_gp_npar  + 3*npl + len(tr_files)*2 + rv_gp_npar + 1 + i ]  )
            ppp+='%d %d %d %d %d %d %d %d\n'%(0,0,0,0,0,0,0,0)     
        ppp+='%f\n%d\n'%(par[len(vel_files)*2 +7*npl],0) # information about linear trend
        ppp+='%f\n%d\n'%(par[len(vel_files)*2 +7*npl + 1],0) # information about linear trend
        ppp+='%f\n'%epoch
        ppp+='%d\n'%hkl
        ppp+='EOF' 
 
#        print(ppp)
        
        text,flag=run_command_with_timeout(ppp, when_to_kill, output=True, pipe=True) # running command generated by the fortran_input function 
        fortranoutput=fortran_output(text,npl,len(vel_files),stmass)
        fit_results=fortranoutput.modfit(print_stat=False)
        rv_loglik = float(fit_results.loglik)
        #rv_loglik = float(text[1][0])        
    else:
        rv_loglik = 0
        
    if(rtg[1]): 
        
        gp_rv_loglik = 0
        
        param_vect = []
        for j in range(1,len(gps.get_parameter_vector())+1):
            param_vect.append(np.log(par[len(vel_files)*2  +7*npl +j]))
        
        gps.set_parameter_vector(np.array(param_vect))
    
        gp_pred = gps.predict(fit_results.o_c, fit_results.jd, return_cov=False)
        o_c_kep = fit_results.o_c - gp_pred
        
        for i in range(len(vel_files)):
            sig2i_gp = 1.0 / (fit_results.rv_err[fit_results.idset==i]**2 + par[i + len(vel_files)]**2 )
            
            gp_rv_loglik += -0.5*(np.sum((o_c_kep[fit_results.idset==i])**2 * sig2i_gp - np.log(sig2i_gp / 2./ np.pi)))
                         

        rv_loglik =  gp_rv_loglik 
        #print(rv_loglik)
    if(rtg[2]): 
        
        if len(tr_files[0]) == 0:
            tr_loglik = 0        
        else: 
            tr_loglik = transit_loglik(tr_files,vel_files,tr_params,tr_model,par,rv_gp_npar,npl,hkl,rtg,tra_gps )

    if np.isnan(rv_loglik).any() or np.isnan(tr_loglik).any():
        return -np.inf
    #print(rv_loglik, tr_loglik,rv_loglik + tr_loglik)        
    return rv_loglik + tr_loglik
        

def run_SciPyOp(obj,   threads=1,  kernel_id=-1,  save_means=False, fileoutput=False, save_sampler=False, **kwargs):      
 
    start_time = time.time()    
    rtg = obj.rtg

    check_temp_RV_file(obj)

    vel_files = []
    for i in range(obj.filelist.ndset): 
         vel_files.append(obj.filelist.files[i].path)  
    
    tr_files = []
    tr_mo = []
    tr_ld = []
    
    for i in range(10):
        if len(obj.tra_data_sets[i]) != 0:
            tr_files.append(obj.tra_data_sets[i])
            tr_mo.append(obj.ld_m[i])
            tr_ld.append(obj.ld_u[i])
   
    tr_model = np.array([tr_mo,tr_ld], dtype=object)
    tr_params = obj.tr_params
       
    npl = obj.npl
    epoch = obj.epoch     
    stmass = obj.params.stellar_mass  
 
    mix_fit = obj.mixed_fit    
    
    if (obj.mod_dynamical):
        if mix_fit[0] == True:
            mod='%s/lib/fr/loglik_dyn+'%obj.cwd
            when_to_kill =  obj.dyn_model_to_kill
        else: 
            mod='%s/lib/fr/loglik_dyn'%obj.cwd  
            when_to_kill =  obj.dyn_model_to_kill
    else:
        mod='%s/lib/fr/loglik_kep'%obj.cwd
        when_to_kill =  obj.kep_model_to_kill   


    nll = lambda *args: -lnprob_new(*args)
  
    obj.prepare_for_mcmc(rtg = rtg)    
    pp = obj.par_for_mcmc #.tolist()
    ee = obj.e_for_mcmc #.tolist() 
    bb = np.array(obj.b_for_mcmc)
    pr_nr = np.array(obj.nr_pr_for_mcmc)
    jeff_nr = np.array(obj.jeff_pr_for_mcmc)

    
    flags = obj.f_for_mcmc 
    par = np.array(obj.parameters) 
    
    obj.bound_error = False
    obj.bound_error_msg = ""

    for l in range(len(pp)): 
        if not bb[l,0] <= pp[l] <= bb[l,1]:
            obj.bound_error = True
            obj.bound_error_msg = "Parameter %s is initially out of bounds. Please set the initial parametrs within the parameter limits!"%ee[l]
            return obj
 
    priors = [pr_nr,jeff_nr]
    
    opt = {"eps":obj.dyn_model_accuracy*1e-13,"dt":obj.time_step_model*86400.0,
           "when_to_kill":when_to_kill,"copl_incl":obj.copl_incl,"hkl":obj.hkl,"cwd":obj.cwd, "gr_flag":obj.gr_flag} 
    
    if obj.init_fit == True: 
        flags = []
     
    
#    print(par)
#    print(pp)
    #print(bb)
   # print(flags)  
    #print(rtg)

    gps = []
    if (rtg[1]):
        initiategps(obj, kernel_id=kernel_id) 
        gps = obj.gps
        rv_gp_npar = len(gps.get_parameter_vector())
    else:
        rv_gp_npar = 0 

    tra_gps = []
    if (rtg[3]) and len(tr_files) != 0:
        initiate_tansit_gps(obj)     
        tra_gps = obj.tra_gps   
        tra_gp_npar = len(tra_gps.get_parameter_vector())
    else:
        tra_gp_npar = 0     

    if obj.SciPy_min_use_1 == obj.SciPy_min[0]:
         options1=obj.Simplex_opt  
         fit_bounds = None
    elif obj.SciPy_min_use_1 == obj.SciPy_min[1]:
         options1=obj.Powell_opt   
         fit_bounds = bb              
    elif obj.SciPy_min_use_1 == obj.SciPy_min[2]:
         options1=obj.CG_opt 
         fit_bounds = bb           
    elif obj.SciPy_min_use_1 == obj.SciPy_min[3]:
         options1=obj.BFGS_opt  
         fit_bounds = bb                    
    elif obj.SciPy_min_use_1 == obj.SciPy_min[4]:
         options1=obj.Newton_cg_opt 
         fit_bounds = bb                    
    elif obj.SciPy_min_use_1 == obj.SciPy_min[5]:
         options1=obj.L_BFGS_B_opt  
         fit_bounds = bb                    
    elif obj.SciPy_min_use_1 == obj.SciPy_min[6]:
         options1=obj.TNC_opt     
         fit_bounds = bb                       
   # elif obj.SciPy_min_use_1 == obj.SciPy_min[7]:
   #      options1=obj.COBYLA_opt       
   #      fit_bounds = None
    elif obj.SciPy_min_use_1 == obj.SciPy_min[8]:
         options1=obj.SLSQP_opt 
         fit_bounds = None         
   # elif obj.SciPy_min_use_1 == obj.SciPy_min[9]:
   #      options1=obj.TNC_opt                 
   # elif obj.SciPy_min_use_1 == obj.SciPy_min[10]:
   #      options1=obj.TNC_opt             
    else:
         options1={'disp': True}  


    if obj.SciPy_min_use_2 == obj.SciPy_min[0]:
         options2=obj.Simplex_opt  
         fit_bounds = None
    elif obj.SciPy_min_use_2 == obj.SciPy_min[1]:
         options2=obj.Powell_opt   
         fit_bounds = bb              
    elif obj.SciPy_min_use_2 == obj.SciPy_min[2]:
         options2=obj.CG_opt 
         fit_bounds = bb           
    elif obj.SciPy_min_use_2 == obj.SciPy_min[3]:
         options2=obj.BFGS_opt  
         fit_bounds = bb                    
    elif obj.SciPy_min_use_2 == obj.SciPy_min[4]:
         options2=obj.Newton_cg_opt 
         fit_bounds = bb                    
    elif obj.SciPy_min_use_2 == obj.SciPy_min[5]:
         options2=obj.L_BFGS_B_opt  
         fit_bounds = bb                    
    elif obj.SciPy_min_use_2 == obj.SciPy_min[6]:
         options2=obj.TNC_opt     
         fit_bounds = bb                       
  #  elif obj.SciPy_min_use_2 == obj.SciPy_min[7]:
#         options2=obj.COBYLA_opt       
#         fit_bounds = None
    elif obj.SciPy_min_use_2 == obj.SciPy_min[8]:
         options2=obj.SLSQP_opt 
         fit_bounds = bb    
    else:
         options2={'disp': True}   
    
    
    if len(flags) == 0:
        method1 = 'TNC'
        n1 = 0
        n2 = 0
        
        
    else:
        method1 = obj.SciPy_min_use_1
        method2 = obj.SciPy_min_use_2  
        n1 = obj.SciPy_min_N_use_1
        n2 = obj.SciPy_min_N_use_2
   # print(obj.SciPy_min_use_1,obj.SciPy_min_use_2)
    
    ########################### Primary minimizer #########################
    for k in range(n1): # run at least 3 times the minimizer
        #eps = eps/10.0
       # print('running %s %s %s'%(obj.SciPy_min_use_1, obj.SciPy_min_N_use_1, k))
        #print(par) 
        result = op.minimize(nll,  pp, args=(mod, par,flags, npl,vel_files, tr_files, tr_model, tr_params,  epoch, stmass, bb, priors, gps, tra_gps, rtg, mix_fit, opt ),
                             method=method1,bounds=fit_bounds, options=options1)       
                            #  bounds=bb, tol=None, callback=None, options={'eps': 1e-08, 'scale': None, 'offset': None, 'mesg_num': None, 'maxCGit': -1, 'maxiter': None, 'eta': -1, 'stepmx': 0, 'accuracy': 0, 'minfev': 0, 'ftol': -1, 'xtol': -1, 'gtol': -1, 'rescale': -1, 'disp': True})        
        pp = result["x"]
        #print(par)
        print(method1,' Done!')
       # print("Best fit par.:", result["x"])

    ########################### Secondary minimizer #########################

    for k in range(n2): # run at least 3 times the minimizer
        #print(k,xtol)
      #  print('running %s %s %s'%(obj.SciPy_min_use_2, obj.SciPy_min_N_use_2, k))
        result = op.minimize(nll, pp, args=(mod,par,flags, npl,vel_files, tr_files, tr_model, tr_params, epoch, stmass, bb, priors, gps, tra_gps, rtg, mix_fit, opt ), 
                             method=method2,bounds=fit_bounds, options=options2)
        pp = result["x"]
        print(method2,' Done!')
        
       # print("Best fit par.:", result["x"])
  
    obj.par_for_mcmc = pp  
 #   print(obj.par_for_mcmc)
    newparams = obj.generate_newparams_for_mcmc(obj.par_for_mcmc)   
    obj.overwrite_params(newparams)  

    obj.correct_elements()
    obj.hack_around_rv_params() 
    

    if obj.type_fit["RV"] == True and obj.type_fit["Transit"] == False:
        obj.fitting(minimize_fortran=True, minimize_loglik=True, amoeba_starts=0, npoints= obj.model_npoints, outputfiles=[1,1,1]) # this will help update some things 
    elif obj.type_fit["RV"] == False and obj.type_fit["Transit"] == True:
        obj.loglik = transit_loglik(tr_files,vel_files,tr_params,tr_model,par,rv_gp_npar,obj.npl,obj.hkl, obj.rtg, obj.tra_gps )
    elif obj.type_fit["RV"] == True and obj.type_fit["Transit"] == True:
        obj.fitting(minimize_fortran=True, minimize_loglik=True, amoeba_starts=0, npoints= obj.model_npoints, outputfiles=[1,1,1]) # this will help update some things 
        tr_loglik = transit_loglik(tr_files,vel_files,tr_params,tr_model,par,rv_gp_npar,obj.npl,obj.hkl, obj.rtg , obj.tra_gps ) 
        obj.loglik     =   obj.loglik +  tr_loglik
 
#    obj.loglik = -result["fun"]
       
    errors = [[0.0,0.0] for i in range(len(pp))] 
   
    obj = return_results(obj, pp, ee, par, flags, npl, vel_files, tr_files, tr_model, tr_params, epoch, stmass, bb, priors, gps, tra_gps, rtg, mix_fit, errors)

    obj.init_fit = False
   # print(obj.loglik)
    print("--- %s seconds ---" % (time.time() - start_time))     
    
    return obj





def return_results(obj, pp, ee, par,flags, npl,vel_files, tr_files, tr_model, tr_params, epoch, stmass, bb, pr_nr, gps, tra_gps, rtg, mix_fit, errors):
                
                
    for j in range(len(flags)):
        par[flags[j]] = pp[j] 
        
#    print(par)
#    print(pp)
  #  print(flags)
       
    
    if (rtg[1]):
        if obj.gp_kernel == 'RotKernel':
            for j in range(len(gps.get_parameter_vector())):
                obj.GP_rot_params[j] = par[len(vel_files)*2  +7*npl +2 +j]
                #print(obj.doGP,obj.gp_kernel)
            
        if obj.gp_kernel == 'SHOKernel':
            for j in range(len(gps.get_parameter_vector())):
                obj.GP_sho_params[j] = par[len(vel_files)*2  +7*npl +2 +j]          
        
    if rtg[1]:
        rv_gp_npar = len(gps.get_parameter_vector())
        get_gps_model(obj)  
    else:
        rv_gp_npar = 0 
        
    if rtg[3]:
        tra_gp_npar = len(tra_gps.get_parameter_vector())
        get_transit_gps_model(obj)  
    else:
        tra_gp_npar = 0         
        
    if (rtg[3]) and len(tr_files) != 0:
        if obj.tra_gp_kernel == 'RotKernel':
            for j in range(len(tra_gps.get_parameter_vector())):
                obj.tra_GP_rot_params[j] = par[len(vel_files)*2  +7*npl  + rv_gp_npar  + 3*npl + len(tr_files)*2 + rv_gp_npar + 2 +j]
            
        if obj.tra_gp_kernel == 'SHOKernel':
            for j in range(len(tra_gps.get_parameter_vector())):
                obj.tra_GP_sho_params[j] = par[len(vel_files)*2  +7*npl  + rv_gp_npar  + 3*npl + len(tr_files)*2 + rv_gp_npar + 2 +j]          
            
        
        
    for i in range(npl):   
#        if obj.hkl == True:
           # ecc_ = np.sqrt(par[len(vel_files)*2 +7*i+2]**2 + par[len(vel_files)*2 +7*i+3]**2)
#            om_  = np.degrees(np.arctan2(par[len(vel_files)*2 +7*i+2],par[len(vel_files)*2 +7*i+3]))%360
#            Ma_  = (par[len(vel_files)*2 +7*i+4] - om_)%360.0  
#            obj.params.update_M0(i,Ma_)                 
        if obj.hkl == False:
            obj.params.update_M0(i,par[len(vel_files)*2 +7*i+4]) 
            obj.M0[i] = float(par[len(vel_files)*2 +7*i+4])
        
        obj.t0[i]     = par[len(vel_files)*2 +7*npl +2 +rv_gp_npar + 3*i] #0.0  #time of inferior conjunction
        
        obj.pl_a[i]   = par[len(vel_files)*2 +7*npl +2 +rv_gp_npar + 3*i+1] #15  #semi-major axis (in units of stellar radii)
        obj.pl_rad[i] = par[len(vel_files)*2 +7*npl +2 +rv_gp_npar + 3*i+2] #0.15   #planet radius (in units of stellar radii)   
        
        obj.omega_dot[i] = par[len(vel_files)*2  +7*npl  + rv_gp_npar  + 3*npl + len(tr_files)*2 + tra_gp_npar + 2 + i ]
       # print(obj.t0[i],par[len(vel_files)*2 +7*i+4])
 
    j =0 
    for i in range(10):        
        if len(obj.tra_data_sets[i]) == 0:
            continue
        else:
            obj.tra_off[i] =      par[len(vel_files)*2 +7*npl + 2 +rv_gp_npar + 3*npl + j]
            print(obj.tra_off[i],obj.tra_jitt[i])
            j = j +1
    j =0 
    for i in range(10):        
        if len(obj.tra_data_sets[i]) != 0:
            obj.tra_jitt[i] = abs(par[len(vel_files)*2 +7*npl + 2 +rv_gp_npar + 3*npl + len(tr_files) + j])
            print(obj.tra_off[i],obj.tra_jitt[i])
            j = j +1            
            

    if len(flags) != 0:    
        print("Best lnL: %s"%obj.loglik)
        print("Best fit par.:")  
     
        for j in range(len(pp)):
            #print(ee[j] + "  =  %s"%pp[j])
            #print("{0:{width}s} = {1:{width}.{precision}f}".format(ee[j], pp[j] , width = 10, precision = 4))
            #print("{0:{width}s} = {1:{width}.{precision}f} + {2:{width}.{precision}f} - {3:{width}.{precision}f}".format(ee[j], pp[j],errors[j][0],errors[j][1], width = 10, precision = 4))
            print("{0:{width2}s} = {1:{width}.{precision}f}  -{2:{width}.{precision}f} +{3:{width}.{precision}f}".format(ee[j], pp[j],errors[j][0],errors[j][1], width2 = 16, width = 10, precision = 9))
  
    obj.gps = []
    obj.tra_gps = []    
    
    return obj
 
 

def lnprior(p,b,priors): 
    
    loglik_lnpr = 0
    for j in range(len(p)): 
        if p[j] <= b[j,0] or p[j] >= b[j,1]:
            return -np.inf
        if priors[0][j,2] == True:
            loglik_lnpr = loglik_lnpr + normalprior(p[j],priors[0][j])
        if priors[1][j,2] == True:
            loglik_lnpr = loglik_lnpr + jeffereys_prior(p[j],priors[1][j])  
    return loglik_lnpr    



def jeffereys_prior(p,b):
    
    loglik =   np.log( 1.0/  ( p*(np.log(b[1]) - np.log(b[0]) )) )     
    return loglik

def normalprior(p,b):  
    
    loglik = np.log(1.0/(np.sqrt(2.0*np.pi)*b[1])*np.exp(-(p-b[0])**2.0/(2.0*b[1]**2.0)))
    return loglik

def lnprob_new(p, program, par, flags, npl, vel_files, tr_files, tr_model, tr_params, epoch, stmass, b, priors , gps, tra_gps, rtg, mix_fit, opt):
    
    lp = lnprior(p,b,priors)
    if not np.isfinite(lp):
        return -np.inf
    return lp + model_loglik(p, program, par, flags, npl, vel_files, tr_files, tr_model, tr_params, epoch, stmass, gps, tra_gps, rtg, mix_fit, opt)   




########## Dynesty Work in progress!!! #######
    
def run_nestsamp_bg(obj):
    start_time = time.time()  
     
    target_name = 'ns_run'
    
    print("Nested Sampling is running behind the GUI. For status, see the main terminal.")
    file_ses = open(r"%s.ses"%target_name, 'wb')
    dill.dump(obj, file_ses)
    file_ses.close()    
 
    if sys.version_info[0] == 2:
        os.system("python2 ./lib/run_ns_from_ses.py -ses ./%s.ses %s"%(target_name,target_name))  
    elif sys.version_info[0] == 3:
        os.system("python3 ./lib/run_ns_from_ses.py -ses ./%s.ses %s"%(target_name,target_name))            
          

    file_ses2 = open(r"%s_out.ses"%target_name, 'rb')
    obj = dill.load(file_ses2)
    file_ses2.close()
     
    print("--- %s seconds ---" % (time.time() - start_time))  
    os.system("rm %s.ses"%target_name)
    os.system("rm %s_out.ses"%target_name)
    
    return obj

 

def run_nestsamp(obj, **kwargs):
    
    '''Performs nested sampling and saves results'''
    

    #from contextlib import closing
    #from CustomNestedSampler import CustomNestedSampler
    
   # print("from %s CPS you are using %s CPUs"%(multiprocessing.cpu_count(),threads))
   # if obj.ns_threads == 'max':
   #     threads = multiprocessing.cpu_count()    
      
    
    start_time = time.time()   
    
    rtg = obj.rtg

    check_temp_RV_file(obj)

    #nll = lambda *args: -lnprob_new(*args)
    
    vel_files = []
    for i in range(obj.filelist.ndset): 
        # path for each dataset      
        vel_files.append(obj.filelist.files[i].path)

    tr_files = []
    tr_mo = []
    tr_ld = []
    
    for i in range(10):
        if len(obj.tra_data_sets[i]) != 0:
            tr_files.append(obj.tra_data_sets[i])
            tr_mo.append(obj.ld_m[i])
            tr_ld.append(obj.ld_u[i])
   
    tr_model = np.array([tr_mo,tr_ld], dtype=object)
    tr_params = obj.tr_params
    
    npl = obj.npl
    epoch = obj.epoch
    stmass = obj.params.stellar_mass
    
    mix_fit = obj.mixed_fit
    
    if (obj.mod_dynamical):
        if mix_fit[0] == True:
            mod='%s/lib/fr/loglik_dyn+'%obj.cwd
            when_to_kill =  obj.dyn_model_to_kill
            #print(mix_fit[0],mod)
        else:
            mod='%s/lib/fr/loglik_dyn'%obj.cwd
            when_to_kill =  obj.dyn_model_to_kill

    else:
        mod='%s/lib/fr/loglik_kep'%obj.cwd
        when_to_kill =  obj.kep_model_to_kill

 
    obj.prepare_for_mcmc(rtg = rtg)
    pp = obj.par_for_mcmc #.tolist()
    ee = obj.e_for_mcmc #.tolist()
    bb = np.array(obj.b_for_mcmc)
    pr_nr = np.array(obj.nr_pr_for_mcmc)
    jeff_nr = np.array(obj.jeff_pr_for_mcmc)
    
    flags = obj.f_for_mcmc
    par = np.array(obj.parameters)
    
    obj.bound_error = False
    obj.bound_error_msg = ""

    for l in range(len(pp)):
        if not bb[l,0] <= pp[l] <= bb[l,1]:
            obj.bound_error = True
            obj.bound_error_msg = "Parameter %s is initially out of bounds. Please set the initial parametrs within the parameter limits!"%ee[l]
            return obj
 

 
    
    priors = [pr_nr,jeff_nr]
    
    level = (100.0- obj.nest_percentile_level)/2.0

    opt = {"eps":obj.dyn_model_accuracy*1e-13,"dt":obj.time_step_model*86400.0,
           "when_to_kill":when_to_kill,"copl_incl":obj.copl_incl,"hkl":obj.hkl,
           "cwd":obj.cwd, "gr_flag":obj.gr_flag,"ns_samp_method":obj.ns_samp_method}
    #print(par)
    #print(flags)
   # print(bb)
   # print(pp)
    
    gps = []
    if (rtg[1]):
        initiategps(obj)
        gps = obj.gps

    tra_gps = []
    if (rtg[3]) and len(tr_files) != 0:
        initiate_tansit_gps(obj)
        tra_gps = obj.tra_gps

    ndim, nwalkers = len(pp), len(pp)*obj.live_points_fact
     
    ################## prior TESTS ########################

    def prior_transform(p): 

        u_trans = np.zeros(len(p)) 
        for j in range(len(p)): 

            if priors[0][j,2] == True:
                u_trans[j] = trans_norm(p[j],bb[j][0],bb[j][1])
            elif priors[1][j,2] == True:
                u_trans[j] = trans_loguni(p[j],bb[j][0],bb[j][1]) 
            else:
                u_trans[j] = trans_uni(p[j],bb[j][0],bb[j][1])
        return u_trans   
  

    def trans_norm(p ,mu,sig):
        return stats.norm.ppf(p,loc=mu,scale=sig)

    def trans_uni(p,a,b):
        return a + (b-a)*p
    
    def trans_loguni(p,a,b):
        return np.exp(np.log(a) + p*(np.log(b)-np.log(a)))    
    
    
    def partial_func2(pp):
        loglik = model_loglik(pp, mod, par, flags, npl, vel_files, tr_files, tr_model, tr_params, epoch, stmass, gps, tra_gps, rtg, mix_fit, opt) 
        #loglik = lnprob_new(pp, mod, par, flags, npl, vel_files, tr_files, tr_params, epoch, stmass, bb, priors, gps, rtg, mix_fit)
    #    print(loglik)
        return loglik
    
    ################## TESTS ########################    
    #partial_func = FunctionWrapper(lnprob_new,
   #                 (mod, par, flags, npl, vel_files, tr_files, tr_params, epoch, stmass, bb, priors, gps, rtg, mix_fit) )
  
    partial_func = FunctionWrapper(model_loglik, (mod, par, flags, npl, vel_files, tr_files, tr_model, tr_params, epoch, stmass, gps, tra_gps, rtg, mix_fit, opt) )
    
    #from multiprocessing import Pool, cpu_count    
   # from pathos.multiprocessing import ProcessingPool as Pool
    from pathos.pools import ProcessPool as Pool
#    from multiprocessing import Pool    
#    from contextlib import closing    
 
       
    dynesty_samp = obj.ns_samp_method
    print_progress = obj.std_output #std_output
    threads = int(obj.ns_threads)
    stop_crit = obj.stop_crit
    Dynamic_nest = obj.Dynamic_nest

    thread = Pool(ncpus=threads)

    
    if Dynamic_nest == False:
        print("'Static' Nest. Samp. is running, please wait... (still under tests!)")

        if threads > 1:
            #with closing(Pool(processes=threads)) as thread:
            #    sampler = dynesty.NestedSampler(partial_func, prior_transform, ndim, nlive=nwalkers, pool = thread, 
            #                                    queue_size=threads, sample = dynesty_samp)
     
            #    sampler.run_nested(print_progress=print_progress,dlogz=stop_crit) #dlogz=stop_crit,
            #    thread.close() 
            #    thread.join() 
            #    thread.clear()
            sampler = dynesty.NestedSampler(partial_func, prior_transform, ndim, nlive=nwalkers, pool = thread, 
                                                queue_size=threads, sample = dynesty_samp)
            sampler.run_nested(print_progress=print_progress,dlogz=stop_crit) #dlogz=stop_crit,
            thread.close() 
            thread.join() 
            thread.clear() 
                
        else:
             sampler = dynesty.NestedSampler(partial_func, prior_transform, ndim, nlive=nwalkers, sample = dynesty_samp)
             sampler.run_nested(print_progress=print_progress,dlogz=stop_crit)

 
        obj.dyn_res = sampler.results
        obj.dyn_res.summary()
        
    else:
        print("'Dynamic' Nest. Samp. is running, please wait... (still under tests!)")
        
        if threads > 1:
#            with closing(Pool(processes=threads)) as thread:
#                sampler = dynesty.DynamicNestedSampler(partial_func, prior_transform, ndim, pool = thread,
#                                                       queue_size=threads, sample = dynesty_samp, bound='multi') # nlive=nwalkers, 
     
#                sampler.run_nested(print_progress=print_progress,dlogz_init=stop_crit,nlive_init=nwalkers) #nlive_init=nwalkers, , nlive_batch=1
#                thread.close() 
#                thread.join() 
#                thread.clear()

            sampler = dynesty.DynamicNestedSampler(partial_func, prior_transform, ndim, pool = thread,
                                                   queue_size=threads, sample = dynesty_samp, bound='multi') # nlive=nwalkers, 
 
            sampler.run_nested(print_progress=print_progress,dlogz_init=stop_crit,nlive_init=nwalkers) #nlive_init=nwalkers, , nlive_batch=1
            thread.close() 
            thread.join() 
            thread.clear()            
            
        else:
             sampler = dynesty.DynamicNestedSampler(partial_func, prior_transform, ndim, sample = dynesty_samp, bound='multi')
             sampler.run_nested(print_progress=print_progress,dlogz_init=stop_crit,nlive_init=nwalkers) #nlive_init=nwalkers, 
 
    
        obj.dyn_res = sampler.results

        res = ("niter: {:d}\n"
               "ncall: {:d}\n"
               "eff(%): {:6.3f}\n"
               "logz: {:6.3f} +/- {:6.3f}".format( obj.dyn_res.niter, sum(obj.dyn_res.ncall),
                       obj.dyn_res.eff, obj.dyn_res.logz[-1], obj.dyn_res.logzerr[-1]))

        print('Summary\n=======\n'+res)

    
   # print("--- %s seconds ---" % (time.time() - start_time))  
    ln = np.hstack(sampler.results.logl)
    samples = np.array(sampler.results.samples)
 
        
    if obj.ns_fileoutput == True:
       # start_time = time.time()   
       # print("Please wait... writing the ascii file")  

        outfile = open(str(obj.nest_sample_file), 'w') # file to save samples
        for j in range(len(samples)):
            outfile.write("%s  " %(ln[j]))
            for z in range(len(pp)):
                outfile.write("%s  " %(samples[j,z]))
            outfile.write("\n")
        outfile.close()        

 
    obj.nest_stat["mean"] = get_mean_of_samples(sampler.results.samples,len(pp))
    samp_maxlnl, maxlnl = get_best_lnl_of_samples(sampler.results.samples,ln, len(pp))
    obj.nest_stat["best"] = samp_maxlnl
    obj.nest_stat["mode"] = get_mode_of_samples(sampler.results.samples,len(pp))
 
    
    if (obj.ns_save_means):
        obj.par_for_mcmc = obj.nest_stat["mean"] 
        pp = obj.nest_stat["mean"]  
        
    elif (obj.ns_save_maxlnL):
        obj.par_for_mcmc = obj.nest_stat["best"]  
        pp =  obj.nest_stat["best"]
              
    elif (obj.ns_save_mode):
        obj.par_for_mcmc = obj.nest_stat["mode"]  
        pp =  obj.nest_stat["mode"]  
    # else:
   #     pp = obj.par_for_mcmc
 
    sampler.samples = sampler.results.samples
       
    new_par_errors = [[float(obj.par_for_mcmc[i] - np.percentile(sampler.results.samples[:,i], [level])), float(np.percentile(sampler.results.samples[:,i], [100.0-level])-obj.par_for_mcmc[i])] for i in range(len(obj.par_for_mcmc))] 

    newparams = obj.generate_newparams_for_mcmc(obj.par_for_mcmc)        
   
    obj.fitting(minimize_loglik=True, amoeba_starts=0, npoints=obj.model_npoints, outputfiles=[1,1,1]) # this will help update some things 

    obj.update_with_mcmc_errors(new_par_errors)
    
    obj.overwrite_params(newparams)
    
    obj.hack_around_rv_params() 
 
 
    if (obj.ns_save_means):
        obj.loglik = maxlnl
        
    elif (obj.ns_save_maxlnL):
        obj.loglik = maxlnl

        
    obj = return_results(obj, pp, ee, par, flags, npl,vel_files, tr_files,  tr_model, tr_params, epoch, stmass, bb, priors, gps, tra_gps, rtg, mix_fit, new_par_errors)

 
    if(obj.ns_save_sampler):
        obj.sampler=sampler             
        obj.sampler_saved=True 
  #      sampler.reset()     
  #  else:   
    sampler.reset()
    
    # To avoid memory leak
    if (rtg[1]):    
        obj.gps = []
        gps = []
       # del gps
       # del obj.gps 
    if (rtg[3]):        
       # del tra_gps
       # del obj.tra_gps
       obj.tra_gps = []    
       tra_gps = []

    print("--- %s seconds ---" % (time.time() - start_time))     
    
    return obj
 
def run_mcmc_bg(obj):
    start_time = time.time()  
     
    target_name = 'mcmc_run'
    
    print("MCMC is running behind the GUI. For status, see the main terminal.")
    file_ses = open(r"%s.ses"%target_name, 'wb')
    dill.dump(obj, file_ses)
    file_ses.close()    
 
    if sys.version_info[0] == 2:    
        os.system("python2 ./lib/run_mcmc_from_ses.py -ses ./%s.ses %s"%(target_name,target_name))  
    elif sys.version_info[0] == 3:    
        os.system("python3 ./lib/run_mcmc_from_ses.py -ses ./%s.ses %s"%(target_name,target_name))            
          

    file_ses2 = open(r"%s_out.ses"%target_name, 'rb')
    obj = dill.load(file_ses2)
    file_ses2.close()
     
    print("--- %s seconds ---" % (time.time() - start_time))
    os.system("rm %s.ses"%target_name)
    os.system("rm %s_out.ses"%target_name)
    
    return obj

def run_mcmc(obj, **kwargs):
    
    '''Performs MCMC and saves results'''  
    
    #if threads == 'max':
    #    threads = multiprocessing.cpu_count()    
    
    start_time = time.time()   
    
    rtg = obj.rtg
    check_temp_RV_file(obj)

    #nll = lambda *args: -lnprob_new(*args)
    
    vel_files = []
    for i in range(obj.filelist.ndset): 
        # path for each dataset      
        vel_files.append(obj.filelist.files[i].path)  

    tr_files = []
    tr_mo = []
    tr_ld = []
    
    for i in range(10):
        if len(obj.tra_data_sets[i]) != 0:
            tr_files.append(obj.tra_data_sets[i])
            tr_mo.append(obj.ld_m[i])
            tr_ld.append(obj.ld_u[i])
   
    tr_model = np.array([tr_mo,tr_ld], dtype=object)    
    tr_params = obj.tr_params
    
    npl = obj.npl
    epoch = obj.epoch     
    stmass = obj.params.stellar_mass    

    mix_fit = obj.mixed_fit    
    
    if (obj.mod_dynamical):
        if mix_fit[0] == True:
            mod='%s/lib/fr/loglik_dyn+'%obj.cwd
            when_to_kill =  obj.dyn_model_to_kill
            #print(mix_fit[0],mod) 
        else: 
            mod='%s/lib/fr/loglik_dyn'%obj.cwd  
            when_to_kill =  obj.dyn_model_to_kill

    else:
        mod='%s/lib/fr/loglik_kep'%obj.cwd
        when_to_kill =  obj.kep_model_to_kill   
 
   # print(mod)
    #program='./lib/fr/%s_%s'%(minimized_value,mod) 
 
    obj.prepare_for_mcmc(rtg = rtg)    
    pp = obj.par_for_mcmc #.tolist()
    ee = obj.e_for_mcmc #.tolist() 
    bb = np.array(obj.b_for_mcmc)

    pr_nr = np.array(obj.nr_pr_for_mcmc)
    jeff_nr = np.array(obj.jeff_pr_for_mcmc)
    
    flags = obj.f_for_mcmc 
    par = np.array(obj.parameters)  
    
    obj.bound_error = False
    obj.bound_error_msg = ""

    for l in range(len(pp)): 
        if not bb[l,0] <= pp[l] <= bb[l,1]:
            obj.bound_error = True
            obj.bound_error_msg = "Parameter %s is initially out of bounds. Please set the initial parametrs within the parameter limits!"%ee[l]
            return obj

    #for k in range(len(pp)):
    #    print(ee[k],pp[k],bb[k],pr_nr[k],jeff_nr[k])


    priors = [pr_nr,jeff_nr]
    level = (100.0- obj.percentile_level)/2.0
 
    
    
    opt = {"eps":obj.dyn_model_accuracy*1e-13,"dt":obj.time_step_model*86400.0,
           "when_to_kill":when_to_kill,"copl_incl":obj.copl_incl,"hkl":obj.hkl,"cwd":obj.cwd, "gr_flag":obj.gr_flag}    

    gps = []
    if (rtg[1]):
        initiategps(obj)     
        gps = obj.gps
        
    tra_gps = []
    if (rtg[3]) and len(tr_files) != 0:
        initiate_tansit_gps(obj)     
        tra_gps = obj.tra_gps        


    #from pathos.multiprocessing import ProcessingPool as Pool
    from pathos.pools import ProcessPool as Pool
    #from pathos.threading import ThreadPool as Pool
    pool=Pool(ncpus=obj.mcmc_threads)
    #pool=Pool(processes=threads-1)
    #from pathos.threading import ThreadPool
    #import mkl
   # mkl.set_num_threads(1)    
    #from multiprocessing import Pool 
    #pool=Pool(threads)  
    
    #from schwimmbad import MPIPool

    #with MPIPool() as pool:

   #     if not pool.is_master():
   #         pool.wait()
   #         sys.exit(0) 
            
    ndim, nwalkers = len(pp), len(pp)*obj.nwalkers_fact

    pos = [pp + obj.gaussian_ball*np.random.rand(ndim) for i in range(nwalkers)]
 
    sampler = CustomSampler(nwalkers, ndim, lnprob_new, args=(mod, par, flags, npl, vel_files, tr_files,  tr_model, tr_params, epoch, stmass, bb, priors, gps, tra_gps, rtg, mix_fit,opt), pool=pool)
 

    #sampler = CustomSampler(nwalkers, ndim, lnprob_new, args=(mod, par, flags, npl, vel_files, tr_files,  tr_model, tr_params, epoch, stmass, bb, priors, gps, tra_gps, rtg, mix_fit,opt), threads = threads)


    # burning phase
    pos, prob, state  = sampler.run_mcmc(pos,obj.mcmc_burning_ph)
    sampler.reset()  
    
    # now perform the MCMC
    pos, prob, state  = sampler.run_mcmc(pos,obj.mcmc_ph)
    #ln = np.hstack(sampler.lnprobability)
    sampler.save_samples(obj.f_for_mcmc,obj.filelist.ndset,obj.npl)

    pool.close()
    pool.join() 
    pool.clear()   

 #  print("--- %s seconds ---" % (time.time() - start_time))     

    #print(type(sampler.samples), len(sampler.samples))
    #print(type(sampler.samples[:,0]), len(sampler.samples[:,0]))

            
    fileoutput = obj.mcmc_fileoutput
    if (fileoutput):
     #   start_time = time.time()   
    #    print("Please wait... writing the ascii file")          
        
        outfile = open(str(obj.mcmc_sample_file), 'w') # file to save samples
        for j in range(len(sampler.samples)):
            outfile.write("%s  " %(sampler.lnL[j]))        #BUG here!!!!!
            for z in range(len(pp)):
                outfile.write("%s  " %(sampler.samples[j,z]))
            outfile.write("\n")
        outfile.close()    
   #     print("--- Done for ---")           
   #     print("--- %s seconds ---" % (time.time() - start_time))     
     
   # start_time = time.time()   
 
    obj.mcmc_stat["mean"] = sampler.means
    obj.mcmc_stat["best"] = sampler.minlnL
    obj.mcmc_stat["mode"] = get_mode_of_samples(sampler.samples,len(pp))


    # Now we will save new parameters and their errors (different + and - errors in this case). Flag save_means determines if we want to take means as new best fit parameters or stick to old ones and calculate errors with respect to that           
    if (obj.mcmc_save_means):
        obj.par_for_mcmc = obj.mcmc_stat["mean"] # we will not need to keep the old parameters in this attribbute, so let's store the means now
        pp = obj.mcmc_stat["mean"]  
        
    elif (obj.mcmc_save_maxlnL):
        obj.par_for_mcmc = obj.mcmc_stat["best"] # we will not need to keep the old parameters in this attribbute, so let's store the means now
        pp =  obj.mcmc_stat["best"]
              
    elif (obj.mcmc_save_mode):
        obj.par_for_mcmc = obj.mcmc_stat["mode"] # we will not need to keep the old parameters in this attribbute, so let's store the means now
        pp =  obj.mcmc_stat["mode"]  
    # else:
   #     pp = obj.par_for_mcmc
        
  
    new_par_errors = [[float(obj.par_for_mcmc[i] - np.percentile(sampler.samples[:,i], [level])),float(np.percentile(sampler.samples[:,i], [100.0-level])-obj.par_for_mcmc[i])] for i in range(len(obj.par_for_mcmc))] 
    
    newparams = obj.generate_newparams_for_mcmc(obj.par_for_mcmc)        
    #print(newparams.GP_params)
    #current_GP_params=newparams.GP_params.gp_par # because calling fitting will overwrite them
   # print(current_GP_params)

   
    obj.fitting(minimize_loglik=True, amoeba_starts=0, npoints=obj.model_npoints, outputfiles=[1,1,1]) # this will help update some things 

    obj.update_with_mcmc_errors(new_par_errors)
    
    obj.overwrite_params(newparams)
   # print(new_par_errors)
    
    obj.hack_around_rv_params() 
 
   # obj.params.update_GP_params(current_GP_params)
 
    if (obj.mcmc_save_means):
        obj.loglik = sampler.lnL_min 
        
    elif (obj.mcmc_save_maxlnL):
        obj.loglik = sampler.lnL_min 

    
    obj = return_results(obj, pp, ee, par, flags, npl,vel_files, tr_files, tr_model, tr_params, epoch, stmass, bb, priors, gps,tra_gps, rtg, mix_fit, new_par_errors)

    if(obj.mcmc_save_sampler):
        obj.sampler=sampler             
        obj.sampler_saved=True           
    else:   
        sampler.reset()
 
#    obj.gps = []
    # To avoid memory leak
    if (rtg[1]):    
        obj.gps = []
        gps = []
       # del gps
       # del obj.gps 
    if (rtg[3]):        
       # del tra_gps
       # del obj.tra_gps
       obj.tra_gps = []    
       tra_gps = []


    
    print("--- %s seconds ---" % (time.time() - start_time))     
    
    return obj





class FunctionWrapper(object):
    """
    This is a hack to make the likelihood function pickleable when ``args``
    or ``kwargs`` are also included.
    """
    def __init__(self, f, args, kwargs=None):
        self.f = f
        self.args = [] if args is None else args
        self.kwargs = {} if kwargs is None else kwargs

    def __call__(self, x):
        try:
            result = self.f(x, *self.args, **self.kwargs)
            #print(x, result)
            return result
        except:  # pragma: no cover
            import traceback
            print("emcee: Exception while calling your likelihood function:")
            print("  params:", x)
            print("  args:", self.args)
            print("  kwargs:", self.kwargs)
            print("  exception:")
            traceback.print_exc()
            raise


                  
class signal_fit(object):
 

    def __init__(self, inputfile='init.init', name='session', readinputfile=False): 
        
        #### Old stuff; some these must be replaced! ####
        
        # saving the name for the inputfile and the information that it has not yet been processed
        self.inputfile = inputfile
        self.inputfile_read=False
        self.name = name # for example the name of a host star of the planetary system, any name with which we will identify a given signal_fit object, preferably the same as the name of the signal_fit object
        self.input_processed = False
        self.never_saved=True # important for print_info function
        self.stellar_mass_known=False 
        self.filelist=rvfile_list(0,[],[])
        self.mod_dynamical=False
        self.epoch=0.0
        self.npl=0
        self.use=use_flags([False]*10,[False]*10,[False]*70,False,False) 
        self.params=parameters([0.0]*10,[0.0]*10,[0.0]*70,0.0,DEFAULT_STELLAR_MASS) 
        self.param_errors=parameter_errors([0.0]*10,[0.0]*10,[0.0]*70,0.0,0.0) 
        self.bounds = parameter_bounds([0.0,0.0]*10,[0.0,0.0]*10,[0.0,0.0]*70,[0.0,0.0],[0.0,0.0]*4,[0.0,0.0])  
       
        
        ########## new stuff ##########
        self.init_pl_params()
        self.init_mcmc_par()
        self.init_nest_par()
                
        self.fit_performed = False
        self.model_saved=False
        self.stat_saved=False      
        self.masses=[0.0]*10
        self.semimajor=[0.0]*10
        self.f_for_mcmc=[]
        self.par_for_mcmc=[]
        self.e_for_mcmc=[]  
        self.b_for_mcmc=[] 
        
        #### TBD here ###########
     
        self.init_St_params()   
        
        self.init_GP()  
        self.init_transit_GP()

        self.init_hkl()  
        self.init_omega_dot()
             
        self.init_RV_jitter()       
        self.init_RV_offset()
        self.init_RV_lintr()
        self.init_RV_quadtr()        
        
        self.init_tra_jitter()       
        self.init_tra_offset()        
        
        self.init_st_mass()
        self.cwd = os.getcwd()

        self.init_pl_arb()
        self.init_orb_evol_arb()
        
        
       

        self.type_fit = {"RV": True,"Transit": False,"TTV":False}
 
        self.gr_flag = False
        self.hkl = False
        self.copl_incl = False        
        self.rtg = [True,False,False,False]
 
        
        if (readinputfile):
            self.read_input_file() # after running this a file should be processed and initial values saved in the object properties, unless there was some error
            self.correct_elements() # to make sure angles are between 0 and 360 etc. 
            #self.prepare_for_mcmc() # this will initiate default parameter bounds, in case the user wants to verify if values make sense he will be able to use the verify_params_with_bounds function anytime 
            self.inputfile_read=True           

        # in this case we need to create a new kernel, but we will only give it this information which is needed for plotting
        self.filelist.read_rvfiles(self.params.offsets)
        self.fit_results=kernel(summary(parameters(0.0,0.0,0.0,0.0,DEFAULT_STELLAR_MASS),parameter_errors([0],[0],[0],0,0)), self.filelist.time,self.filelist.rvs,self.filelist.rv_err,np.zeros(len(self.filelist.time)),np.zeros(len(self.filelist.time)),self.filelist.time,0,0,0,self.filelist.idset,0,0,0,0,0,0)
        self.loglik=0.0
        self.rms=0.0
        self.chi2=0.0
        self.reduced_chi2=0.0
        self.sampler=None
        self.sampler_saved=False
        
        self.init_auto_fit()
        
        
        self.colors = ['#0066ff',  '#ff0000','#66ff66','#00ffff','#cc33ff','#ff9900','#cccc00','#3399ff','#990033','#339933','#808080']
 
      
        self.init_orb_evol()
        
        self.tls = []
        self.tls_o_c = []
       
        self.gls = []
        self.gls_o_c =[]

        #self.mlp = []
        
        self.init_dynfit_settings()    
        
 
        self.ph_data = {k: [] for k in range(9)}
        self.ph_model = {k: [] for k in range(9)}

        self.ph_data_tra = {k: [] for k in range(9)}
        self.ph_model_tra = {k: [] for k in range(9)}
        
        self.parameters = []

        self.ttv_data_sets = {k: [] for k in range(10)}
        self.act_data_sets = {k: [] for k in range(10)}
        self.tra_data_sets = {k: [] for k in range(10)}
        self.rv_data_sets  = {k: [] for k in range(10)}
        
        self.pyqt_symbols_rvs = {k: 'o' for k in range(10)} # ['o','t','t1','t2','t3','s','p','h','star','+','d'] 
        self.pyqt_symbols_act = {k: 'o' for k in range(10)} # ['o','t','t1','t2','t3','s','p','h','star','+','d'] 
        self.pyqt_symbols_tra = {k: 'o' for k in range(10)} # ['o','t','t1','t2','t3','s','p','h','star','+','d'] 
        
        self.pyqt_symbols_size_rvs = {k: 6 for k in range(10)} #[6,6,6,6,6,6,6,6,6,6] # 
        self.pyqt_symbols_size_act = {k: 4 for k in range(10)} #[4,4,4,4,4,4,4,4,4,4] #  
        self.pyqt_symbols_size_tra = {k: 2 for k in range(10)} #[2,2,2,2,2,2,2,2,2,2] # 

        self.act_colors = ['#0066ff',  '#ff0000','#66ff66','#00ffff','#cc33ff','#ff9900','#cccc00','#3399ff','#990033','#339933','#666699']
        self.tra_colors = ['#0066ff',  '#ff0000','#66ff66','#00ffff','#cc33ff','#ff9900','#cccc00','#3399ff','#990033','#339933','#000000']
        self.rvs_colors = ['#0066ff',  '#ff0000','#66ff66','#00ffff','#cc33ff','#ff9900','#cccc00','#3399ff','#990033','#339933','#666699']

        self.gls_colors = ['#ff0000',  '#ff0000','#ff0000','#ff0000','#ff0000','#ff0000','#ff0000','#ff0000','#ff0000','#ff0000','#000000']

         
        self.init_sciPy_minimizer()
        self.init_ld_model()
        
        self.init_xyz()
 
        ############################################ 
        try:
           self.init_transit_params()
        except:
            if batman_not_found == True:
                print("You dont have the 'batman' package")            
            else:
                print("Most likely you have the wrong 'batman' package")
            self.tr_params = []




    def init_ld_model(self):

                
        self.ld_models = ["uniform", "linear", "quadratic", "nonlinear"]
        self.ld_m = ["quadratic", "quadratic",  "quadratic", "quadratic", "quadratic", "quadratic", "quadratic", "quadratic", "quadratic", "quadratic"]    #limb darkening model

        self.ld_u = {k: [0.1, 0.3 ] for k in range(10)}    

        self.ld_u_quad = {k: [0.1, 0.3 ] for k in range(10)}    
        self.ld_u_nonlin = {k: [0.5,0.1, 0.3,-0.1] for k in range(10)}    
        self.ld_u_lin    = {k: [0.3] for k in range(10)}         

        
        self.ld_u_quad_use   = {k: [False, False] for k in range(10)}         
        self.ld_u_nonlin_use = {k: [False, False,False, False] for k in range(10)}         
        self.ld_u_lin_use    = {k: [False] for k in range(10)}         

        self.ld_u_quad_err   = {k: [[0,0], [0,0]] for k in range(10)}         
        self.ld_u_nonlin_err = {k: [[0,0], [0,0],[0,0], [0,0]] for k in range(10)}         
        self.ld_u_lin_err    = {k: [0,0] for k in range(10)}         
        
        self.ld_u_quad_bound      = {k: np.array([[-1,1],[-1,1]]) for k in range(10)}
        self.ld_u_nonlin_bound    = {k: np.array([[-1,1],[-1,1],[-1,1],[-1,1]]) for k in range(10)}
        self.ld_u_lin_bound       = {k: np.array([-1,1]) for k in range(10)} 
 
        self.ld_u_quad_norm_pr    = {k: np.array([[0,1, False],[0,1, False]]) for k in range(10)}
        self.ld_u_nonlin_norm_pr  = {k: np.array([[0,1, False],[0,1, False],[0,1, False],[0,1, False]]) for k in range(10)}
        self.ld_u_lin_norm_pr     = {k: np.array([0.1,0.05, False]) for k in range(10)} 

        self.ld_u_quad_jeff_pr    = {k: np.array([[0,1, False],[0,1, False]]) for k in range(10)}
        self.ld_u_nonlin_jeff_pr  = {k: np.array([[0,1, False],[0,1, False],[0,1, False],[0,1, False]]) for k in range(10)}
        self.ld_u_lin_jeff_pr     = {k: np.array([0.1,0.05, False]) for k in range(10)} 
 
        self.ld_u_quad_str        = {k: [[r'ld-quad-1$_%s$'%str(k+1)],[r'ld-quad-2$_%s$'%str(k+1)]] for k in range(10)} 
        self.ld_u_nonlin_str      = {k: [[r'ld-quad-1$_%s$'%str(k+1)],[r'ld-quad-2$_%s$'%str(k+1)],[r'ld-quad-3$_%s$'%str(k+1)],[r'ld-quad-4$_%s$'%str(k+1)]] for k in range(10)} 
        self.ld_u_lin_str         = {k: [r'ld-quad-1$_%s$'%str(k+1)] for k in range(10)} 


        ############################################


#    def constants(self):   
#        self.AU

 
    def init_auto_fit(self):
        self.auto_fit_max_pl = 2
        self.auto_fit_allow_ecc = False
        self.auto_fit_FAP_level = 0.001
      
    def init_xyz(self):         
        self.xyz_mass = {k: [] for k in range(10)}       
        self.xzy    = {k: [] for k in range(10)}
        self.uvw    = {k: [] for k in range(10)}
        self.rpl    = {k: [] for k in range(10)}
        self.rhill  = {k: [] for k in range(10)}

#        self.hack_around_rv_params()
#        self.calc_hkl()
#        self.calc_ewm()  
  
        
        
    def init_pl_params(self): 

        #### RV #####
        self.K    = {k: 50.0 for k in range(9)}
        self.P    = {k: 100.0 + 50*k for k in range(9)}
        self.e    = {k: 0.0  for k in range(9)}
        self.w    = {k: 0.0 for k in range(9)}
        self.M0   = {k: 0.0 for k in range(9)}
        self.i    = {k: 90.0 for k in range(9)}
        self.Node = {k: 0.0 for k in range(9)}     
        self.w_dot= {k: 0.0 for k in range(9)}     
       
        
        

        self.K_err    = {k: np.array([0,0]) for k in range(9)}
        self.P_err    = {k: np.array([0,0]) for k in range(9)}
        self.e_err    = {k: np.array([0,0]) for k in range(9)}
        self.w_err    = {k: np.array([0,0]) for k in range(9)}
        self.M0_err   = {k: np.array([0,0]) for k in range(9)}
        self.i_err    = {k: np.array([0,0]) for k in range(9)}
        self.Node_err = {k: np.array([0,0]) for k in range(9)}  
        self.w_dot_err= {k: np.array([0,0]) for k in range(9)}  
        

        self.K_use    = {k: False for k in range(9)}
        self.P_use    = {k: False for k in range(9)}
        self.e_use    = {k: False for k in range(9)}
        self.w_use    = {k: False for k in range(9)}
        self.M0_use   = {k: False for k in range(9)}
        self.i_use    = {k: False for k in range(9)}
        self.Node_use = {k: False for k in range(9)}  
        self.w_dot_use= {k: False for k in range(9)}  
               
        
        self.K_bound    = {k: np.array([0,10000]) for k in range(9)}
        self.P_bound    = {k: np.array([0,100000]) for k in range(9)}
        self.e_bound    = {k: np.array([0,0.999]) for k in range(9)}
        self.w_bound    = {k: np.array([0.0, 360.0]) for k in range(9)}
        self.M0_bound   = {k: np.array([0.0, 360.0]) for k in range(9)}
        self.i_bound    = {k: np.array([0.0, 180.0]) for k in range(9)}
        self.Node_bound = {k: np.array([0.0, 360.0]) for k in range(9)}
        self.w_dot_bound= {k: np.array([0.0, 360.0]) for k in range(9)}
        
        
        self.K_norm_pr    = {k: np.array([50,100, False]) for k in range(9)}
        self.P_norm_pr    = {k: np.array([150,30, False]) for k in range(9)}
        self.e_norm_pr    = {k: np.array([0.0,0.1, False]) for k in range(9)}
        self.w_norm_pr    = {k: np.array([0, 90, False]) for k in range(9)}
        self.M0_norm_pr   = {k: np.array([0, 90, False]) for k in range(9)}
        self.i_norm_pr    = {k: np.array([90, 90, False]) for k in range(9)}
        self.Node_norm_pr = {k: np.array([0, 360.0, False]) for k in range(9)} 
        self.w_dot_norm_pr= {k: np.array([0, 360.0, False]) for k in range(9)} 
     
        
        self.K_jeff_pr    = {k: np.array([50,100, False]) for k in range(9)}
        self.P_jeff_pr    = {k: np.array([150,30, False]) for k in range(9)}
        self.e_jeff_pr    = {k: np.array([0.0,0.1, False]) for k in range(9)}
        self.w_jeff_pr    = {k: np.array([0, 90, False]) for k in range(9)}
        self.M0_jeff_pr   = {k: np.array([0, 90, False]) for k in range(9)}
        self.i_jeff_pr    = {k: np.array([90, 90, False]) for k in range(9)}
        self.Node_jeff_pr = {k: np.array([0, 360.0, False]) for k in range(9)}          
        self.w_dot_jeff_pr= {k: np.array([0, 360.0, False]) for k in range(9)}          
        
        
        
        self.K_str    = {k: r'K$_%s$'%chr(98+k)  for k in range(9)}
        self.P_str    = {k: r'P$_%s$'%chr(98+k)  for k in range(9)}
        self.e_str    = {k: r'e$_%s$'%chr(98+k)  for k in range(9)}
        self.w_str    = {k: r'$\omega_%s$'%chr(98+k)  for k in range(9)}
        self.M0_str   = {k: r'MA$_%s$'%chr(98+k)  for k in range(9)}
        self.i_str    = {k: r'i$_%s$'%chr(98+k)  for k in range(9)}
        self.Node_str = {k: r'$\Omega_%s$'%chr(98+k)  for k in range(9)}          
        self.w_dot_str= {k: r'$\dot{\omega_%s}$'%chr(98+k)  for k in range(9)}          
  

        
        #### transit #####       
        self.t0      = {k: 0 for k in range(9)}
        self.pl_a    = {k: 15 for k in range(9)}
        self.pl_rad  = {k: 0.10 for k in range(9)}       

        self.t0_use      = {k: False for k in range(9)}
        self.pl_a_use    = {k: False for k in range(9)}
        self.pl_rad_use  = {k: False for k in range(9)}         

        self.t0_err      = {k: np.array([0,0])  for k in range(9)}
        self.pl_a_err    = {k: np.array([0,0])  for k in range(9)}
        self.pl_rad_err  = {k: np.array([0,0])  for k in range(9)}         

        self.t0_bound      = {k: np.array([-10000,10000]) for k in range(9)}
        self.pl_a_bound    = {k: np.array([0,100]) for k in range(9)}
        self.pl_rad_bound  = {k: np.array([0,10000]) for k in range(9)} 
 
        self.t0_norm_pr      = {k: np.array([0,1, False]) for k in range(9)}
        self.pl_a_norm_pr     = {k: np.array([10,10, False]) for k in range(9)}
        self.pl_rad_norm_pr   = {k: np.array([0.1,0.05, False]) for k in range(9)} 
        
        self.t0_jeff_pr      = {k: np.array([0,1, False]) for k in range(9)}
        self.pl_a_jeff_pr     = {k: np.array([10,10, False]) for k in range(9)}
        self.pl_rad_jeff_pr   = {k: np.array([0.1,0.05, False]) for k in range(9)}         
        
 
        self.t0_str      = {k: r't0$_%s$'%chr(98+k) for k in range(9)} 
        self.pl_a_str    = {k: r'pl_a$_%s$'%chr(98+k) for k in range(9)} 
        self.pl_rad_str  = {k: r'pl_rad$_%s$'%chr(98+k) for k in range(9)} 
        


    def init_hkl(self) : 
        
        #### h,k,l #####             
        self.e_sinw        = {k: 0 for k in range(9)}
        self.e_cosw        = {k: 0 for k in range(9)}
        self.lamb          = {k: 0 for k in range(9)}       

        self.e_sinw_use    = {k: False for k in range(9)}
        self.e_cosw_use    = {k: False for k in range(9)}
        self.lamb_use    = {k: False for k in range(9)}         

        self.e_sinw_err    = {k: np.array([0.0,0.0])  for k in range(9)}
        self.e_cosw_err    = {k: np.array([0.0,0.0])  for k in range(9)}
        self.lamb_err    = {k: np.array([0.0,0.0])  for k in range(9)}         

        self.e_sinw_bound  = {k: np.array([-1.0,1.0]) for k in range(9)}
        self.e_cosw_bound  = {k: np.array([-1.0,1.0]) for k in range(9)}
        self.lamb_bound  = {k: np.array([0.0,360.0]) for k in range(9)} 
 
        self.e_sinw_norm_pr = {k: np.array([0.0,0.1, False]) for k in range(9)}
        self.e_cosw_norm_pr = {k: np.array([0.0,0.1, False]) for k in range(9)}
        self.lamb_norm_pr = {k: np.array([0.0,30.0, False]) for k in range(9)} 
        
        self.e_sinw_jeff_pr = {k: np.array([-1,1, False]) for k in range(9)}
        self.e_cosw_jeff_pr = {k: np.array([-1,1, False]) for k in range(9)}
        self.lamb_jeff_pr = {k: np.array([0.0,360.0, False]) for k in range(9)}         
        
 
        self.e_sinw_str     = {k: r'$e sin(\omega_%s)$'%chr(98+k) for k in range(9)} 
        self.e_cosw_str     = {k: r'$e cos(\omega_%s)$'%chr(98+k) for k in range(9)} 
        self.lamb_str     = {k: r'$\lambda_%s$'%chr(98+k) for k in range(9)} 
                

        ######## derived #####################
        self.t_peri = {k: 0.0 for k in range(9)}
        
        




    def init_omega_dot(self) :       
        
        self.omega_dot       = {k: 0 for k in range(9)}
        self.omega_dot_err   = {k: np.array([0,0]) for k in range(9)}
        self.omega_dot_use   = {k: False for k in range(9)}
        self.omega_dot_str   = {k: r'$\omega_%s dot$'%k for k in range(9)}     
        self.omega_dot_bounds  = {k: np.array([0.0,10000.0] )for k in range(9)} 
        self.omega_dot_norm_pr = {k: np.array([0.0,50.0, False] )for k in range(9)} 
        self.omega_dot_jeff_pr = {k: np.array([0.0,360.0, False] )for k in range(9)} 


    def init_RV_jitter(self) :       
        
        self.jitt      = {k: 0 for k in range(10)}
        self.jitt_err  = {k: np.array([0,0]) for k in range(10)}
        self.jitt_use  = {k: True for k in range(10)}
        self.jitt_str  = {k: r'RV jitt$_%s$'%k for k in range(10)}     
        self.jitt_bounds  = {k: np.array([0.0,10000.0] )for k in range(10)} 
        self.jitt_norm_pr = {k: np.array([1.0,5.0, False] )for k in range(10)} 
        self.jitt_jeff_pr = {k: np.array([1.0,5.0, False] )for k in range(10)} 
        
    def init_RV_offset(self) :       
        
        self.rvoff      = {k: 0 for k in range(10)}
        self.rvoff_err  = {k: np.array([0,0])  for k in range(10)}
        self.rvoff_use  = {k: True for k in range(10)}
        self.rvoff_str  = {k: r'RV off$_%s$'%k for k in range(10)}     
        self.rvoff_bounds  = {k: np.array([-1000000.0,1000000.0] )for k in range(10)}        
        self.rvoff_norm_pr = {k: np.array([0,100.0, False] )for k in range(10)}    
        self.rvoff_jeff_pr = {k: np.array([0,100.0, False] )for k in range(10)}    
       
        
    def init_tra_jitter(self) :       
        
        self.tra_jitt      = {k: 0 for k in range(10)}
        self.tra_jitt_err  = {k: np.array([0,0]) for k in range(10)}
        self.tra_jitt_use  = {k: False for k in range(10)}
        self.tra_jitt_str  = {k: r'transit jitt$_%s$'%k for k in range(10)}     
        self.tra_jitt_bounds  = {k: np.array([-0.2,0.2] )for k in range(10)} 
        self.tra_jitt_norm_pr = {k: np.array([0.00,0.1, False] )for k in range(10)} 
        self.tra_jitt_jeff_pr = {k: np.array([0.00,0.1, False] )for k in range(10)} 
        
    def init_tra_offset(self) :       
        
        self.tra_off      = {k: 0 for k in range(10)}
        self.tra_off_err  = {k: np.array([0,0])  for k in range(10)}
        self.tra_off_use  = {k: False for k in range(10)}
        self.tra_off_str  = {k: r'transit off$_%s$'%k for k in range(10)}     
        self.tra_off_bounds  = {k: np.array([-1,2] )for k in range(10)}        
        self.tra_off_norm_pr = {k: np.array([1,0.1, False] )for k in range(10)}    
        self.tra_off_jeff_pr = {k: np.array([1,0.1, False] )for k in range(10)}    
               
 
    def init_RV_lintr(self) :       
         
        self.rv_lintr      = {k: 0 for k in range(1)}
        self.rv_lintr_err  = {k: np.array([0,0]) for k in range(1)}
        self.rv_lintr_use  = {k: False for k in range(1)}
        self.rv_lintr_str  = {k: r'RV lin.tr' for k in range(1)}     
        self.rv_lintr_bounds  = {k: np.array([-1.0,1.0]) for k in range(1)} 
        self.rv_lintr_norm_pr = {k: np.array([0,0.001, False]) for k in range(1)} 
        self.rv_lintr_jeff_pr = {k: np.array([0,0.001, False]) for k in range(1)} 
        
    def init_RV_quadtr(self) :       
         
        self.rv_quadtr      = 0
        self.rv_quadtr_err  = [0,0]
        self.rv_quadtr_use  = False
        self.rv_quadtr_str  = {k: r'RV quad.tr' for k in range(1)}     
        self.rv_quadtr_bounds  = {k: np.array([-1.0,1.0]) for k in range(1)} 
        self.rv_quadtr_norm_pr = {k: np.array([0,0.001, False]) for k in range(1)} 
        self.rv_quadtr_jeff_pr = {k: np.array([0,0.001, False]) for k in range(1)}         
               
        
    def init_st_mass(self) :       
         
        self.st_mass      = {k: 1 for k in range(1)}
        self.st_mass_err  = {k: np.array([0,0]) for k in range(1)}
        self.st_mass_use  = {k: False for k in range(1)}
        self.st_mass_str  = {k: r'St mass' for k in range(1)}     
        self.st_mass_bounds  = {k: np.array([0.01,100]) for k in range(1)}         
        self.st_mass_norm_pr = {k: np.array([1,0.2, False]) for k in range(1)}         
        self.st_mass_jeff_pr = {k: np.array([1,0.2, False]) for k in range(1)}         


    def init_St_params(self): 

        self.stellar_mass = 1.0
        self.stellar_mass_err = 0.0         
        
        self.stellar_radius = 1.0
        self.stellar_radius_err = 0.0    
        
        self.stellar_luminosity = 1.0
        self.stellar_luminosity_err = 0.0            
        
        self.stellar_Teff = 5777.0
        self.stellar_Teff_err = 0.0             
 
        self.stellar_vsini = 2.0
        self.stellar_vsini_err = 0.0        
        
        self.stellar_rotation = 25.0
        self.stellar_rotation_err = 0.0     
      
     
    def init_GP(self):

        self.doGP = False
        self.gps=[]
        
        self.GP_rot_params = [1,10,15,1]# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
        self.GP_rot_err = [0,0,0,0]
        self.GP_rot_use = [False,False,False,False]  
        self.GP_rot_str = [r'Amp', r't', r'per', r'fact']# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway 
        self.GP_rot_bounds  = {k: np.array([0.0,100000.0]) for k in range(len(self.GP_rot_params))}        
        self.GP_rot_norm_pr = {k: np.array([0.0,10.0, False]) for k in range(len(self.GP_rot_params))}        
        self.GP_rot_jeff_pr = {k: np.array([0.0,10.0, False]) for k in range(len(self.GP_rot_params))}        
                
        
        self.GP_sho_params     = [100,1,0.05]# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
        self.GP_sho_err = [0,0,0]
        self.GP_sho_use = [False,False,False]  
        self.GP_sho_str = [r'S', r'Q', r'omega']# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway 
        self.GP_sho_bounds     = {k: np.array([0.0,100000.0]) for k in range(len(self.GP_sho_params))}        
        self.GP_sho_norm_pr    = {k: np.array([0.0,10.0, False]) for k in range(len(self.GP_sho_params))}        
        self.GP_sho_jeff_pr    = {k: np.array([0.0,10.0, False]) for k in range(len(self.GP_sho_params))}        
                

        self.gp_model_curve = {k: 0 for k in range(10)}
        self.gp_model_data  = {k: 0 for k in range(10)}
        
        self.gp_kernels = ['SHOKernel','RotKernel']
        self.gp_kernel = self.gp_kernels[0]
        
        
    def init_transit_GP(self):

        self.tra_doGP = False
        self.tra_gps=[]
        
        self.tra_GP_rot_params = [1,10,15,1]# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
        self.tra_GP_rot_err = [0,0,0,0]
        self.tra_GP_rot_use = [False,False,False,False]  
        self.tra_GP_rot_str = [r'Amp', r't', r'per', r'fact']# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway 
        self.tra_GP_rot_bounds  = {k: np.array([0.0,100000.0]) for k in range(len(self.tra_GP_rot_params))}        
        self.tra_GP_rot_norm_pr = {k: np.array([0.0,10.0, False]) for k in range(len(self.tra_GP_rot_params))}        
        self.tra_GP_rot_jeff_pr = {k: np.array([0.0,10.0, False]) for k in range(len(self.tra_GP_rot_params))}        
                
        
        self.tra_GP_sho_params     = [100,1,0.05]# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
        self.tra_GP_sho_err = [0,0,0]
        self.tra_GP_sho_use = [False,False,False]  
        self.tra_GP_sho_str = [r'S', r'Q', r'omega']# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway 
        self.tra_GP_sho_bounds     = {k: np.array([0.0,100000.0]) for k in range(len(self.tra_GP_sho_params))}        
        self.tra_GP_sho_norm_pr    = {k: np.array([0.0,10.0, False]) for k in range(len(self.tra_GP_sho_params))}        
        self.tra_GP_sho_jeff_pr    = {k: np.array([0.0,10.0, False]) for k in range(len(self.tra_GP_sho_params))}        
                

        self.tra_gp_model_curve = {k: 0 for k in range(10)}
        self.tra_gp_model_data  = {k: 0 for k in range(10)}
        
        self.tra_gp_kernels = ['SHOKernel','RotKernel']
        self.tra_gp_kernel = self.gp_kernels[0]       
        
       
    def init_transit_params(self): 
        # from the example in github

        self.tr_params = batman.TransitParams()       #object to store transit parameters
       
        # WASP 6
        self.tr_params.t0  = 0.0  #time of inferior conjunction
        self.tr_params.per = 3.36    #orbital period
        self.tr_params.ecc = 0.0      
        self.tr_params.w   = 90.0   #longitude of periastron (in degrees)                  
        self.tr_params.rp  = 0.15   #planet radius (in units of stellar radii)
        self.tr_params.inc = 90. #orbital inclination (in degrees)
        self.tr_params.a   = 15  #semi-major axis (in units of stellar radii)
    
 
        #self.tr_params_use = [False, False,False,False,False,False,False]    

        self.tr_params.limb_dark = "quadratic"      #limb darkening model
        self.tr_params.u =  [0.1, 0.3 ]            
        
        # ld_options = ["uniform", "linear", "quadratic", "nonlinear"]
        #ld_coefficients = [[], [0.3], [0.1, 0.3], [0.5, 0.1, 0.1, -0.1]]       
      
        self.tr_params_use = [False, False,False,False,False,False,False]          
 
 
    def init_sciPy_minimizer(self):
        
       # self.SciPy_min = ['Nelder-Mead','Powell','CG','BFGS','Newton-CG','L-BFGS-B', 'TNC','COBYLA','SLSQP','dogleg','trust-ncg']
        self.SciPy_min = ['Nelder-Mead','Powell','CG','BFGS','Newton-CG','L-BFGS-B', 'TNC','SLSQP']
        
        self.SciPy_min_use_1 = self.SciPy_min[6]
        self.SciPy_min_N_use_1 = 1
        
        self.SciPy_min_use_2 = self.SciPy_min[0]
        self.SciPy_min_N_use_2 = 1      
        
        
        self.Simplex_opt    = {'disp': True, 'maxiter': None, 'return_all': False, 'maxfev': None, 'xtol': 0.0001, 'ftol': 0.0001,'adaptive':True}
        self.Powell_opt     = {'disp': True, 'return_all': False, 'maxiter': None, 'direc': None,  'maxfev': None, 'xtol': 0.0001, 'ftol': 0.0001}
        self.CG_opt         = {'disp': True, 'gtol': 1e-05, 'eps': 1.4901161193847656e-08, 'return_all': False, 'maxiter': None, 'norm': np.inf}
        self.BFGS_opt       = {'disp': True, 'gtol': 1e-05, 'eps': 1.4901161193847656e-08, 'return_all': False, 'maxiter': None, 'norm': np.inf}
        self.Newton_cg_opt  = {'disp': True, 'xtol': 1e-05, 'eps': 1.4901161193847656e-08, 'return_all': False, 'maxiter': None}
        self.L_BFGS_B_opt   = {'disp': True,  'maxcor': 10, 'ftol': 2.220446049250313e-09, 'gtol': 1e-05, 'eps': 1e-08, 'maxfun': 15000, 'maxiter': 15000, 'iprint': -1, 'maxls': 20}
        self.TNC_opt        = {'disp': True, 'eps': 1e-08, 'scale': None, 'offset': None, 'mesg_num': None, 'maxCGit': -1, 'maxiter': None, 'eta': -1, 'stepmx': 0, 'accuracy': 0, 'minfev': 0, 'ftol': -1, 'xtol': -1, 'gtol': -1, 'rescale': -1 }
        #self.COBYLA_opt     = {'disp': True, 'rhobeg': 1.0, 'maxiter': 1000, 'catol': 0.0002 }
        self.SLSQP_opt      = {'disp': True, 'maxiter': 100,  'eps': 1.4901161193847656e-08, 'ftol': 1e-06, 'iprint': 1}
       # self.dogleg_opt    = {'disp': True, 'max_trust_radius': 100,  'eta': 1.4901161193847656e-08, 'gtol': 1e-06 }


   # def init_GLS(self):
        
   #     self.RV_gls_power        = {'power': [], 'base': [] }
   #     self.RV_gls_o_c_power    = {'power': [], 'base': [] }
       
    def init_pl_arb(self): 

        self.arb_st_mass = 1.0
        self.npl_arb = 0
        self.pl_arb_use    = {k: False for k in range(9)}
        
        #### RV #####
        self.K_arb    = {k: 50.0 for k in range(9)}
        self.P_arb    = {k: 300.0 + 50*k for k in range(9)}
        self.e_arb    = {k: 0.0  for k in range(9)}
        self.w_arb    = {k: 0.0 for k in range(9)}
        self.M0_arb   = {k: 0.0 for k in range(9)}
        self.i_arb    = {k: 90.0 for k in range(9)}
        self.Node_arb = {k: 0.0 for k in range(9)}     

        self.mass_arb    = {k: 50.0 for k in range(9)}
        self.a_arb    = {k: 300.0 + 50*k for k in range(9)}

        self.pl_arb_use    = {k: False for k in range(9)}
        self.pl_arb_test = False
 

       
    def init_orb_evol_arb(self):
        
        self.evol_T_arb = {k: [] for k in range(9)}
        self.evol_a_arb = {k: [] for k in range(9)}
        self.evol_e_arb = {k: [] for k in range(9)}
        self.evol_p_arb = {k: [] for k in range(9)}
        self.evol_M_arb = {k: [] for k in range(9)}

        self.evol_T_energy_arb   = 0
        self.evol_energy_arb   = 0
        #self.evol_momentum_arb = 0
        self.evol_momentum_arb = {'lx': [], 'ly': [], 'lz': []}
 
    def init_orb_evol(self):
        
        self.evol_T  = {k: [] for k in range(9)}
        self.evol_a  = {k: [] for k in range(9)}
        self.evol_e  = {k: [] for k in range(9)}
        self.evol_p  = {k: [] for k in range(9)}
        self.evol_M  = {k: [] for k in range(9)}
        self.evol_i  = {k: [] for k in range(9)}
        self.evol_Om = {k: [] for k in range(9)}

        self.evol_T_energy   = 0
        self.evol_energy   = 0
        #self.evol_momentum = 0
        self.evol_momentum = {'lx': [], 'ly': [], 'lz': []}

        self.GR_step = 1000


    def init_mcmc_par(self):     
        self.gaussian_ball = 0.0001        
        self.nwalkers_fact = 4
        
        self.percentile_level = 68.3
        
        self.mcmc_burning_ph = 100 
        self.mcmc_ph = 100 
        self.mcmc_threads= 1
        self.mcmc_fileoutput=True
        self.mcmc_save_means=False 
        self.mcmc_save_mode=False 
        self.mcmc_save_maxlnL=False 
        self.mcmc_save_sampler=True        
                                   
        self.mcmc_sample_file = 'mcmc_samples'
        self.mcmc_corner_plot_file = 'cornerplot.pdf'
        self.mcmc_stat = {"mean": [],"mode": [],"best": []}
        
        
        
    def init_nest_par(self):     
        #self.gaussian_ball = 0.0001        
        self.live_points_fact = 4
        self.nest_percentile_level = 68.3    
        
        self.ns_samp_method_opt = ['slice','unif','rwalk','rstagger','rslice','hslice']
        self.ns_samp_method = self.ns_samp_method_opt[0]     
        
        self.ns_threads=1
        self.Dynamic_nest = False
        self.std_output=False
        self.stop_crit = 0.001
        self.ns_fileoutput=True
        self.ns_save_means=False 
        self.ns_save_mode=False 
        self.ns_save_maxlnL=False 
        self.ns_save_sampler=True    
        
        self.nest_sample_file = 'nest_samp_samples'
        self.nest_corner_plot_file = 'nest_samp_cornerplot.pdf'        
        self.nest_stat = {"mean": [],"mode": [],"best": []}
           
    
    
        
    def init_dynfit_settings(self):
        self.mixed_fit = {0: [False], 1:[1,1,1,1,1,1,1,1,1]}     
        self.fitting_method = 'None'
        self.time_step_model = 10.0
        self.dyn_model_accuracy = 1000.0
        self.dyn_model_to_kill = 86400.0
        self.kep_model_to_kill = 60.0
        self.master_timeout = 86400.0
        
        self.model_npoints = 2000
        
    def update_epoch(self,epoch):
        self.epoch=epoch
        return



    def calc_hkl(self):
 
        for i in range(9):
            self.e_sinw[i] = self.e[i]*np.sin(np.radians(self.w[i]))
            self.e_cosw[i] = self.e[i]*np.cos(np.radians(self.w[i]))
            self.lamb[i]   = (self.w[i] + self.M0[i])%360.0
            
    def calc_ewm(self):            
 
        for i in range(9):
            self.e[i]   = np.sqrt(self.e_sinw[i]**2 + self.e_cosw[i]**2)
            self.w[i]   = np.degrees(np.arctan2(np.radians(self.e_sinw[i]),np.radians(self.e_cosw[i])))
            self.M0[i]  = (self.lamb[i] - self.w[i])%360.0  
            
    def hack_around_rv_params(self):

        
         for i in range(9):        
             self.K[i]    = self.params.planet_params[i*7+0]            
             self.P[i]    = self.params.planet_params[i*7+1]                     
             self.i[i]    = self.params.planet_params[i*7+5]           
             self.Node[i] = self.params.planet_params[i*7+6]   
             
            # print( self.params.planet_params[i*7+2] , self.params.planet_params[i*7+3] )
             if self.hkl ==False:
                 self.e[i]    = self.params.planet_params[i*7+2]            
                 self.w[i]    = self.params.planet_params[i*7+3]            
                 self.M0[i]   = self.params.planet_params[i*7+4] 
                 self.e_sinw[i]    = self.e[i]*np.sin(np.radians(self.w[i]))           
                 self.e_cosw[i]    = self.e[i]*np.cos(np.radians(self.w[i]))            
                 self.lamb[i]      = (self.w[i]  + self.M0[i])%360.0 
 
             if self.hkl ==True:
                 
                 self.e_sinw[i]    = self.params.planet_params[i*7+2]            
                 self.e_cosw[i]    = self.params.planet_params[i*7+3]            
                 self.lamb[i]   = self.params.planet_params[i*7+4] 
                 self.e[i]   = np.sqrt(self.e_sinw[i]**2 + self.e_cosw[i]**2)
                 self.w[i]   = np.degrees(np.arctan2(self.e_sinw[i],self.e_cosw[i]))%360
                 self.M0[i]  = (self.lamb[i] - self.w[i])%360.0     
#                

   
    
        
############################ RV datasets ##########################################      
    def add_rv_dataset(self, name, path, rv_idset = 0):
 
        rv_JD       = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [0])
        rv_data     = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [1])
        rv_data_sig = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [2])
	
        rv_data_set = np.array([rv_JD,rv_data,rv_data_sig])

        #ind = 0 
        #for i in range(10):
       #     if len(self.filelist.idset==i) != 0:
        #        ind += 1
        ####### for now ###########
        self.rv_data_sets[max(self.filelist.idset)+1] =  rv_data_set 
 
        return   


    def remove_rv_dataset(self, rv_idset):
 
        self.rv_data_sets[rv_idset] = []
 
        return   

############################ activity datasets ##########################################      
    def add_act_dataset(self, name, path, act_idset = 0):
 
        act_JD       = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [0])
        act_data     = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [1])
        act_data_sig = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [2])
        act_data_set = np.array([act_JD,act_data,act_data_sig]) 
 
        self.act_data_sets[act_idset] = act_data_set
 
        return   


    def remove_act_dataset(self, act_idset):
 
        self.act_data_sets[act_idset] = []
 
        return   
    
############################ TTV datasets ##########################################      
    def add_ttv_dataset(self, name, path, ttv_idset = 0):
 
        ttv_N        = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [0])
        ttv_data     = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [1])
        ttv_data_sig = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [2])
        ttv_data_set = np.array([ttv_N,ttv_data,ttv_data_sig]) 
 
        self.ttv_data_sets[ttv_idset] = ttv_data_set
 
        return   


    def remove_ttv_dataset(self, ttv_idset):
 
        self.ttv_data_sets[ttv_idset] = []
 
        return
    
############################ transit datasets ##########################################      
    def add_transit_dataset(self, name, path, tra_idset = 0):
 
        tra_JD       = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [0])
        tra_data     = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [1])
        tra_data_sig = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [2])
	
        tra_data_o_c = tra_data
   
    
        tra_data_set = np.array([tra_JD,tra_data,tra_data_sig,tra_data_o_c,tra_data_o_c]) 
 
        self.tra_data_sets[tra_idset] = tra_data_set      
 
        return   


    def remove_transit_dataset(self, tra_idset):
 
        self.tra_data_sets[tra_idset] = []
 
        return   


################ Legacy Code!!! TB removed/updated/replaced ######################


    def add_planet(self,K=50,P=100,e=0,w=0,M0=180,i=90,cap=0,useK=True,useP=True,usee=False,usew=False,useM0=True,usei=False, usecap=False):        
        if(self.npl==10):
            self.params.planet_params=np.concatenate((np.atleast_1d(self.params.planet_params),[0.0]*7)) # to allocate more places in planet_params array
            warnings=Warning_log(['By deafult we assume max 10 planets, to satisfy your request we had to overwrite this rule! More then 20 planets will cause an error in fortran codes, modify them if necessary.'],'Adding a new planet')
            warnings.print_warning_log()
        self.params.update_planet_params_one_planet(self.npl,K,P,e,w,M0,i,cap) #because we count planets from 0, so we use old npl and only then increase it
        self.use.update_use_planet_params_one_planet(self.npl,useK,useP,usee,usew,useM0,usei,usecap)
        self.npl=self.npl+1
        return


    def add_dataset(self,name,path,offset,jitter,useoffset=True,usejitter=True):
        
        
        path =  copy_file_to_datafiles(path)
        
        if(self.filelist.ndset==20):
            self.params.offsets=np.concatenate((np.atleast_1d(self.params.offsets),np.atleast_1d(0.0))) # to allocate more places in offsets array
            self.params.jitters=np.concatenate((np.atleast_1d(self.params.jitters),np.atleast_1d(0.0))) # to allocate more places in offsets array
            warnings=Warning_log(['By deafult we assume max 20 datasets, to satisfy your request we had to overwrite this rule! More then 20 datasets will cause an error in fortran codes, modify them if necessary.'],'Adding a new dataset')
            warnings.print_warning_log()
        flag=self.filelist.add_datafile(name,path) 

        if (flag==1): # flag will be zero if invalid path is provided
            self.params.update_offset(self.filelist.ndset-1,offset)
            self.params.update_jitter(self.filelist.ndset-1,jitter)
            self.use.update_use_offset(self.filelist.ndset-1,useoffset)
            self.use.update_use_jitter(self.filelist.ndset-1,usejitter)
            self.filelist.read_rvfiles(self.params.offsets,justthenewone=True)
            if self.epoch < 1:
                self.update_epoch(self.filelist.first_observation())
        return       
        
    def remove_planet(self,planet):
        if not (planet<self.npl):
            warnings=Warning_log(['Planet index outside of range.'],'Removing planet %d'%planet)
            warnings.print_warning_log()
        else:
            for i in range(planet,self.npl-1):
                self.params.update_planet_params_one_planet(i,self.params.planet_params[7*i+7],self.params.planet_params[7*i+8],self.params.planet_params[7*i+9],self.params.planet_params[7*i+10],self.params.planet_params[7*i+11],self.params.planet_params[7*i+12],self.params.planet_params[7*i+13]) #move planets one index down - not the fastest way, but it ensures we don't change the order of the planets, which is important for user's convenience
                self.use.update_use_planet_params_one_planet(i,self.use.use_planet_params[7*i+7],self.use.use_planet_params[7*i+8],self.use.use_planet_params[7*i+9],self.use.use_planet_params[7*i+10],self.use.use_planet_params[7*i+11],self.use.use_planet_params[7*i+12],self.use.use_planet_params[7*i+13]) #move planets one index down - not the fastest way, but it ensures we don't change the order of the planets, which is important for user's convenience
            self.params.update_planet_params_one_planet(self.npl,0.0,0.0,0.0,0.0,0.0,0.0,0.0) # erasing an unnecessary planet slot - not really necessary, but better to be tidy
            self.use.update_use_planet_params_one_planet(self.npl,False,False,False,False,False,False,False)
            self.npl=self.npl-1
        return

    def remove_dataset(self,number):
        flag=self.filelist.remove_datafile(number) 
        if (flag==1): # flag will be zero if too large dataset index is provided
            for i in range(number,self.filelist.ndset+1):
                self.params.update_offset(i,self.params.offsets[i+1])
                self.params.update_jitter(i,self.params.jitters[i+1])
                self.use.update_use_offset(i,self.use.use_offsets[i+1])
                self.use.update_use_jitter(i,self.use.use_jitters[i+1])
            self.params.update_offset(self.filelist.ndset+1,0.0)
            self.params.update_jitter(self.filelist.ndset+1,0.0)
            self.use.update_use_offset(self.filelist.ndset+1,False)
            self.use.update_use_jitter(self.filelist.ndset+1,False)
        return                                                          

    def update_mod_dynamical(self, mod_dynamical):
        self.mod_dynamical=mod_dynamical
        return     
    
    
    
    
    
    def add_RVbank_dataset(self, name, path, offset=0, jitter= 0, split = False):    

       dirname, basename = os.path.split(path)
        
       BJD = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [0])    
      # indices = np.where(BJD > 2457161.5)
 
       fo = open(path, "r")
       lines = fo.readlines() 
       
       
       name1 = '%s_pre.dat'%basename[:-4]
       name2 = '%s_post.dat'%basename[:-4]    
       path1 = 'datafiles/%s'%name1
       path2 = 'datafiles/%s'%name2
 
       out1 = open('%s'%path1, 'w')
       out2 = open('%s'%path2, 'w')      
       
       
       for i in range(len(lines)):
 
           line = lines[i].split()
           if float(line[0]) <= 2457161.5:
               out1.write(lines[i])
           elif float(line[0]) > 2457161.5:
               out2.write(lines[i])               
 
       out1.close()
       out2.close()       
       
       if split == True:
           if len(BJD[BJD <= 2457161.5]) !=0:
               self.add_dataset(name1,path1,offset,jitter,useoffset=True,usejitter=True)
           if len(BJD[BJD > 2457161.5]) !=0:           
               self.add_dataset(name2,path2,offset,jitter,useoffset=True,usejitter=True)
       else:
           self.add_dataset(name,path,offset,jitter,useoffset=True,usejitter=True)
           
 
       for i in range(5):
           
           act_ind = 11 + (2*i)
           act_data     = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [act_ind])
           act_data = act_data - np.mean(act_data)
           
           act_data_sig = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [act_ind+1])
           act_data_set = np.array([BJD,act_data,act_data_sig]) 
 
           self.act_data_sets[i] = act_data_set   
           #print(act_data[0])
    
       #z = 0
       for ii in range(3):
           
           act_ind = 22 + ii
           act_data     = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [act_ind])
           act_data = act_data - np.mean(act_data)
           if ii == 1:
               act_data = act_data * 1000.0
              # print(act_data[0:10])
               
           act_data_sig = act_data*0.05
           act_data_set = np.array([BJD,act_data,act_data_sig]) 
           i = i +1
          # print(act_data[0])
 
           self.act_data_sets[i] = act_data_set       
    
    def BIC(self):        
        if len(self.fit_results.jd) != 0:    
            BIC = 2*self.loglik + self.fit_results.mfit*np.log(len(self.fit_results.jd))    
        else:
            BIC = 0        
        return BIC
      
    def AIC(self):     
        if len(self.fit_results.jd) != 0:    
            AIC = 2*self.loglik - 2*self.fit_results.mfit     
        else:
            AIC = 0        
        return AIC    
    
    def wrms(self):  
        if len(self.fit_results.o_c) != 0:
            wrms = np.sqrt(np.average(self.fit_results.o_c**2, weights=1/self.fit_results.rv_err)) 
        else:
            wrms = 0
        return wrms        
    
    

    def verify_gp_parameters_number(self):
        # Since parameters, use flags and errors are stored in separate objects, and the number of GP parameters can vary, it can lead to problems. To prevent them, verify things by running this function     	
        gpparamsnumwarnings=Warning_log([],'Verifying consistency with the number of GP parameters')
        npar=self.params.GP_params.npar
        if(len(self.use.use_GP_params)>npar):
            gpparamsnumwarnings.update_warning_list('Too many use flags (%d flags while there are %d parameters)! Additional will be discarded.'%(len(self.use.use_GP_params),npar))
            self.use.use_GP_params=self.use.use_GP_params[:npar]
        if(len(self.use.use_GP_params)<npar):
            gpparamsnumwarnings.update_warning_list('Too few use flags (%d flags while there are %d parameters)! Will assume True for remaining parameters.'%(len(self.use.use_GP_params),npar))
            self.use.use_GP_params=np.concatenate((self.use.use_GP_params,[True]*(npar-len(self.use.use_GP_params))))       
        if(len(self.param_errors.GP_params_errors)>npar):
            gpparamsnumwarnings.update_warning_list('Too many errors (%d errors while there are %d parameters)! Additional will be discarded.'%(len(self.use.use_GP_params),npar))
            self.param_errors.GP_params_errors=self.errors.GP_params_errors[:npar]
        if(len(self.param_errors.GP_params_errors)<npar):
            gpparamsnumwarnings.update_warning_list('Too few use flags (%d errors while there are %d parameters)! Will assume [0.0,0.0] for remaining parameters.'%(len(self.use.use_GP_params),npar))
            self.param_errors.GP_params_errors=np.concatenate((self.param_errors.GP_params_errors,np.repeat([[0.0,0.0]],npar-len(self.use.use_GP_params),axis=0)))   
        gpparamsnumwarnings.print_warning_log()
        return
 
    def overwrite_use(self,useflags, save=True): # overwrite use flags with new ones, but optionally save the old ones so we can return to the later
        oldflags=self.use
        self.use=useflags
        if (save): # save old flags, if requested 
            return oldflags
        else:
            return 
 
    def overwrite_params(self,params, save=True): # overwrite use flags with new ones, but optionally save the old ones so we can return to the later
        oldparams=self.params
        self.params=params
        if (save): # save old flags, if requested 
            return oldparams
        else:
            return 
 
    def verify_input(self): # verify that the input file is built correctly, fix such errors that can be solved automatically (but print out warnings), raise exception if it's not possible
        c=read_file_as_array_of_arrays(self.inputfile) # create an array of arrays storing input data   
        input_warnings=Warning_log([],'Reading the input file') # creating a Warning_log object to store warnings related to reading the file

        # Now we go through each line and check for any incorrect input, saving warning messages whenever we encounter and fix a problem, raising errors if it's not possible to fix the problem automatically

        # Checking the first line (stellar mass)
        if (len(c[0])==0):
            c[0]=[0,DEFAULT_STELLAR_MASS,1]
            input_warnings.update_warning_list('Please provide stellar mass information! Unknown stellar mass will be assumed.')
        elif (c[0][0]!=0 and c[0][0]!=1):
            c[0]=[0,DEFAULT_STELLAR_MASS,1]
            input_warnings.update_warning_list('Stellar_mass_is_known should be a boolean value (0 or 1)! Unknown stellar mass will be assumed.')
        elif (c[0][0]==1 and len(c[0])==1):
            c[0]=[0,DEFAULT_STELLAR_MASS,1]
            input_warnings.update_warning_list('If stellar mass is known please provide the value! Unknown stellar mass will be assumed.')     
        elif (c[0][0]==0 and len(c[0])==1):
            c[0]=[0,DEFAULT_STELLAR_MASS,1]
        elif (c[0][0]==1 and not is_float(c[0][1])):
            c[0]=[0,DEFAULT_STELLAR_MASS,1]
            input_warnings.update_warning_list('Stellar mass should be a real number! Unknown stellar mass will be assumed.')  
        elif (c[0][0]==1 and c[0][1]<=0):
            c[0]=[0,DEFAULT_STELLAR_MASS,1]
            input_warnings.update_warning_list('Stellar mass should be a positive number! Unknown stellar mass will be assumed.')   
        elif (len(c[0])==2):
            c[0]=np.concatenate((c[0],np.atleast_1d(1)))
            input_warnings.update_warning_list('Please provide information if you want to use stellar mass for the fit! Will assume True')
        elif (len(c[0])>3):
            c[0]=c[0][:3]
            input_warnings.update_warning_list('First input line too long, additional entries will be ignored.') 

        # Checking the second line (number of datasets and optional path)
        if (not is_int(c[1][0], bounded=[1,0])): # this will test if number of datasets provided by the user is a positive integer
            c[1][0]=len(c[2])
            input_warnings.update_warning_list('Number of datasets should be a positive integer! The number of file names provided in the third line will be assumed as number of datasets.')          
        if (len(c[1])>2):
            c[1]=c[1][:2]          
            input_warnings.update_warning_list('Second input line too long, additional entries will be ignored.')
         
        # Checking the third line (names of RV files)
        if (len(c[2])<c[1][0]):          
            c[1][0]=len(c[2])
            input_warnings.update_warning_list('Fewer datasets than declared! Will assume an appropriate number of datasets')
        if (len(c[2])>c[1][0]):     
            c[2]=c[2][:int(c[1][0])]          
            input_warnings.update_warning_list('Third input line too long, additional entries will be ignored.')
        # Let's check if all the files provided are there. Better to do it now, so we don't have an error later on, also it's a good moment for a sanity check if there isn't something seriously wrong with the structure of input file
        if (len(c[1])==2):
            if not (c[1][1][-1]=='/'):
                c[1][1]=c[1][1]+'/'
            path=c[1][1]
        else:
            path=DEFAULT_PATH
            
        i=0
        while(i<c[1][0]):
            if not (os.path.isfile(path+c[2][i])): # this will check if all files in the list are valid files, remove those which aren't 
                input_warnings.update_warning_list('%s is not a file, will be skipped! Please make sure the file names (and path, if provided) are spelled correctly.'%(path+c[2][i]))
                if (i==c[1][0]-1):
                    c[2]=c[2][:i]
                else:
                    c[2]=np.concatenate((c[2][:i],c[2][i+1:]))
                c[1][0]=c[1][0]-1 
            else:
                i=i+1        
          
        # Checking the fourth line (offsets)
        if (len(c[3])<c[1][0]):
            c[3].extend(['unknown']*(int(c[1][0]-len(c[3]))))
            input_warnings.update_warning_list('Fewer offset information than datasets! Will assume unknown offsets for remaining datasets.')
        if (len(c[3])>c[1][0]):     
            c[3]=c[3][:int(c[1][0])] 
            input_warnings.update_warning_list('Fourth input line too long, additional entries will be ignored.') 
        # Checking the fifth line (use_offsets)
        if (len(c[4])<c[1][0]):
            c[4].extend([1]*(int(c[1][0]-len(c[4]))))
            input_warnings.update_warning_list('Fewer use_offsets information than datasets! Will assume True for remaining datasets.')
        if (len(c[4])>c[1][0]):      
            c[4]=c[4][:int(c[1][0])] 
            input_warnings.update_warning_list('Fifth input line too long, additional entries will be ignored.')  
 
        # Checking the sixth line (jitters)
        if (len(c[5])<c[1][0]):
            c[5].extend(['unknown']*(int(c[1][0]-len(c[5]))))
            input_warnings.update_warning_list('Fewer jitter information than datasets! Will assume unknown jitter for remaining datasets.')
        if (len(c[5])>c[1][0]):      
            c[5]=c[5][:int(c[1][0])] 
            input_warnings.update_warning_list('Sixth input line too long, additional entries will be ignored.') 
        # Checking the seventh line (use_jitters)
        if (len(c[6])<c[1][0]):
            c[6].extend([1]*(int(c[1][0]-len(c[6]))))
            input_warnings.update_warning_list('Fewer use_jitters information than datasets! Will assume True for remaining datasets.')
        if (len(c[6])>c[1][0]):       
            c[6]=c[6][:int(c[1][0])]          
            input_warnings.update_warning_list('Seventh input line too long, additional entries will be ignored.')         
        # Checking the eighth line (linear trend)
        if not (is_float(c[7][0])):
            c[7]=[0,0]
            input_warnings.update_warning_list('Coefficient for linear trend must be a real number! No linear trend will be assumed')  
        if (len(c[7])==1):
            if (c[7][0]==0):
                c[7].append(0)
                input_warnings.update_warning_list('Please provide information if you want to use linear trend as a fitting parameter or not! Since 0 linear coefficient was provided, we will assume there is no linear trend and won\'t use it for the fit.')           
            else:
                c[7].append(1)
                input_warnings.update_warning_list('Please provide information if you want to use linear trend as a fitting parameter or not! Since non-zero linear coefficient was provided, we will use it as a fitting parameter.')
        if (c[7][1]!=0 and c[7][1]!=1):
            if (c[7][0]==0):
                c[7][1]=0
                input_warnings.update_warning_list('Use_linear_trend should be a boolean value (0 or 1)! Since 0 linear coefficient was provided, we will assume there is no linear trend and won\'t use it for the fit.')               
            else:
                c[7][1]=1
                input_warnings.update_warning_list('Use_linear_trend should be a boolean value (0 or 1)! Since non-zero linear coefficient was provided, we will use it as a fitting parameter.')    
        if (len(c[7])>2): 
            c[7]=c[7][:2] 
            input_warnings.update_warning_list('Eighth input line too long, additional entries will be ignored.')
        #Checking the ninth line (Fitting_in_dynamical_mode and do_mcmc information)           
        if (c[8][0]!=0 and c[8][0]!=1):
            c[8][0]=0
            input_warnings.update_warning_list('Fitting_in_dynamical_mode should be a boolean value: 0 if you want to fin in keplerian mode, and 1 if in dynamical! Keplerian fitting mode will be assumed.')
        if (len(c[8])>1):      
            c[8]=c[8][:1] 
            input_warnings.update_warning_list('Ninth input line too long, additional entries will be ignored.')

        #Checking the tenth line (epoch)
        if (len(c[9])==0): 
            c[9]=['unknown']  
            input_warnings.update_warning_list('Epoch should be a positive number! Unknown epoch will be assumed.')        
        else:
            if (is_float(c[9][0])):
                if (c[9][0]<=0):
                    c[9]=['unknown']  
                    input_warnings.update_warning_list('Epoch should be a positive number! Unknown epoch will be assumed.')         
            if (len(c[9])>1): 
                c[9]=c[9][:1] 
                input_warnings.update_warning_list('Tenth input line too long, additional entries will be ignored.')

        #Checking the eleventh line (number of known planets)           
        if not (is_int(c[10][0], bounded=[1,0], equal=[True,False])): # this will test if number of planets provided by the user is a nonnegative integer
            npl_guess=(len(c)-11)/2  # all the lines after the eleventh one should describe planet data, two lines per planet, so (len(c)-11)/2 should be the number of planets
            c[10][0]=npl_guess
            input_warnings.update_warning_list('The number of planets should be a non-negative integer! %d assumed based on the number of subsequent lines in input file.'%npl_guess)    
        if (len(c[10])>1): 
            c[10]=c[10][:1]         
            input_warnings.update_warning_list('Eleventh input line too long, additional entries will be ignored.')
        if (len(c)>11+4*int(c[10][0])):
            c=c[:11+4*int(c[10][0])] 
            input_warnings.update_warning_list('Input file too long, additional entries (beyond declared number of planets) will be ignored.')
        ppp = ["K","P","e","w","M","i","cap0m"]    
        # Checking the rest of the file - information about planet parameters and whether they should be used as fitting parameters
        i=0 # actual iterator for the upcoming loop
        j=0 # this will refer to a number of the planet in user's input, regardless if it was specified correctly and will be used or not 
        while(i<c[10][0]): # looping over all planets
            j=j+1
            p=11+2*i # index in c referring to considered planet's parameters
            warning_found = False
            for q in p,p+1: # we can check line length for both lines regarding the planet in the same way
                if (len(c[q])>7): 
                    c[q]=c[q][:7] 
                    input_warnings.update_warning_list('Input line %d for planet %d too long, additional entries will be ignored.'%(q-p+1,j))
                if (len(c[q])<5): # in dynamical mode we require the inclinations and lines of nodes to be provided, in keplerian it doesn't matter, the program will not use them anyway
                    c[10][0]=c[10][0]-1
                    # remove these two lines from our consideration 
                    c.pop(p) 
                    c.pop(p)
                    input_warnings.update_warning_list('Input line %d for planet %d too short, please provide initial values for all parameters! This planet will be ignored'%(q-p+1,j))
                    warning_found=True 
                    break
            if (len(c[p])==6):
                c[p]=np.concatenate((c[p],np.atleast_1d(0.0)))
                input_warnings.update_warning_list('Missing line of nodes information for planet %d, will assume 0.0'%j)
            elif (len(c[p])==5):
                c[p]=np.concatenate((c[p],[90.0,0.0]))
                input_warnings.update_warning_list('Missing inclination and line of nodes information for planet %d, will assume 90.0 and 0.0'%j) 
            if(warning_found): # If we found a warning and broke the previous loop we go right into the next iteration of the outer loop. i will not be increased, so we will be referring to appropriate place in c even though we removed some elements from the middle
                continue
            if(c[8][0]==1): 
                l=7
            else: # in the keplerian case we don't care if the values for inclination and cap0m are badly stated
               l=5     
            for k in range(l):
                if not is_float(c[p][k]):          
                    c[10][0]=c[10][0]-1
                    # remove these two lines from our consideration 
                    c.pop(p) 
                    c.pop(p)
                    input_warnings.update_warning_list('Invalid %s value for planet %d, initial values should be real numbers! This planet will be ignored'%(ppp[k],j))
                    warning_found=True
                    break
            if(warning_found):  # If we found a warning and broke the previous loop we go right into the next iteration of the outer loop. i will not be increased, so we will be referring to appropriate place in c even though we removed some elements from the middle
                continue               
            p=p+1 # now we check the next line, with information which parameters should be used as fit parameters
            for k in range(l):
                if (c[p][k]!=0 and c[p][k]!=1):  # this can be easily fixed, so we fix it and do not ignore the planet        
                    c[p][k]=1
                    input_warnings.update_warning_list('Invalid %s use flag for planet %d, information which parameters to use for fitting should be provided as boolean values (0s or 1s)! Default option True will be assumed'%(ppp[k],j))                
            i=i+1 # if everything went well, we increase i
    
        input_warnings.print_warning_log()
        return c
    
            
            
    # this should be the default way to provide input - a file structured according to a manual
    def read_input_file (self):
        if (self.input_processed): 
            raise Exception ('Input file already processed.')

        # file has not been processed, so let's start reading it
        alltheinput = self.verify_input() # store all info from input file in an array of arrays      
        
        ### Reading the first line of the input file - stellar mass 
 
        if (alltheinput[0][0]==0): # is stellar mass known? If not we will calculate planet masses in units of the host star mass, so we set stellar_mass = DEFAULT_STELLAR_MASS=1.0
            self.stellar_mass_known = False
        else: # save stellar mass if we know it's value
            self.stellar_mass_known = True
        stellar_mass = alltheinput[0][1]
        use_stellar_mass=bool(alltheinput[0][2])        

        ### Reading lines 2-7, describing all data sets used for the fit

        # first entry in the second line should be the number of RV data files, we will pass it as a first argument (ndset) to the creation of rvfile_list object
        # the second entry in this line is optional and describes the path to directory with datafiles, if different than default
        if (len(alltheinput[1])==1): # if no path was specified by the user, use default path (./datafiles)
            path = DEFAULT_PATH
        else: 
            path = alltheinput[1][1]
            
       # path =  copy_file_to_datafiles(path) 
 
        #filenames = []

        
        for j in range(len(alltheinput[2])):
            #print(alltheinput[2][j], path)
            self.add_dataset(alltheinput[2][j], path+alltheinput[2][j] ,alltheinput[3][j],alltheinput[4][j])
       
        #filenames.append(copy_file_to_datafiles(path))
        #filenames=np.char.array(alltheinput[2])  # the third line contains a list of RV file names
        
        #print(filenames)
        
        
        #self.filelist=rvfile_list(int(alltheinput[1][0]),filenames,path+filenames) # creating a rvfile_list object
        # arrays for offset and jitter, and their use flags, will be used later to create params and use_flags objects
        use_offsets=[]
        use_jitters=[]
        offsets=[]
        jitters=[]
        for i in range(self.filelist.ndset): # looping through all data sets
            # If offset is unknown it should be indicated in the input file by typing 'unknown' 'NaN' or any other string which does not represent a number.
            # We now check it if offset is provided by the user or unknown (for each dataset independently), so if a given offset is not defined we use mean value of a dataset as a first approximation 
            # the next input line informs whether the offset should be used as a parameter for the fit
            # Analogously for jitter, except if unknown we set DEFAULT_JITTER, and we allow to keep it fixed even if jitter was unknown           
            if (is_float(alltheinput[3][i])):
                offsets=np.concatenate((offsets,np.atleast_1d(alltheinput[3][i]))) 
                use_offsets=np.concatenate((use_offsets,np.atleast_1d(bool(alltheinput[4][i])))) 
            # If offset was unknown then we ignore user's choice on fitting - not fitting, the mean value of dataset is just the first approximation and it is necessary to find a better one by using this offset in the fit  
            else:
                off=self.filelist.files[i].mean_value()
                offsets=np.concatenate((offsets,np.atleast_1d(off))) 
                use_offsets=np.concatenate((use_offsets,np.atleast_1d(True)))       
            if (is_float(alltheinput[5][i])): 
                jitters=np.concatenate((jitters,np.atleast_1d(alltheinput[5][i]))) 
            else:
                jitters=np.concatenate((jitters,np.atleast_1d(DEFAULT_JITTER))) 
            # use_jitter is provided but will always be ignored if Chi^2 fitting method is chosen.
            use_jitters=np.concatenate((use_jitters,np.atleast_1d(bool(alltheinput[6][i])))) 
  
        ### Reading the 8th line, describing linear trend
        
        linear_trend=alltheinput[7][0] # saving initial value for linear trend
        use_linear_trend=bool(alltheinput[7][1]) # saving information whether we want to use linear trend or not 
 
        ### Reading the 9th line, describing fitting mode (keplerian vs dynamical)
        
        self.mod_dynamical=bool(alltheinput[8][0])

        ### Reading the 10th line (epoch)
        
        if (is_float(alltheinput[9][0])): # if epoch is defined properly, we will use it 
            self.epoch = alltheinput[9][0]
        else: # if instead we receive 'no_epoch', 'NaN' or similar string then the epoch was not declared, we will set it as time[0] later on                   
            self.epoch=self.filelist.first_observation()
        ### Reading the rest of the file, describing data for all planets 

        # for each planet the user should provide a list of parameters (K,P,e,w,M,i,cap0m, last two optional in keplerian case) and information which ones to use
        self.npl=int(alltheinput[10][0])
        if (self.npl>NPLMAX):
            raise InputError('Too many planets! Maximum number is %d.'%NPLMAX)
        planet_params=[]
        use_planet_params=[]

 
        for i in range(self.npl):
            planet_params=np.concatenate((planet_params,alltheinput[11+2*i][:7])) # add parameters (K,P,e,w,M) of each planet into the array of planet parameters
            use_planet_params=np.concatenate((use_planet_params,np.array(list(map(bool,alltheinput[11+2*i+1][:7]))))) # saving information which data should be used for the fit            
 
        self.use=use_flags(np.array(list(map(bool,use_offsets))),np.array(list(map(bool,use_jitters))),use_planet_params,use_linear_trend,use_stellar_mass) 
        self.params=parameters(offsets,jitters,planet_params,linear_trend,stellar_mass) 
        # information that all the parameters stored at the moment come from user input, no fit has been performed. We need three flags, which will refer to three optional outputs of the fortran code
        self.fit_performed = False
        self.fitting_method = 'None'
        self.model_saved=False
        self.stat_saved=False
        # information that the file is processed and initial data is already stored accordingly
        self.input_processed = True   
        return 
        
   ###################### correct_elements #########################
    #correct initial parameter values for kepfit so that angles are in the [0,180.0] interval etc.

    def correct_elements(self):

        for i in range(self.filelist.ndset): # if jitter is negative, set it to absolute value
            if self.params.jitters[i]<0.0:
                self.params.jitters[i]=-self.params.jitters[i]
 
        if self.hkl == True:
            return
    
        for i in range(self.npl): # looping over all planets
            j = 7*i # 5 parameters per planet

            if self.params.planet_params[j+1] < 0.0:  # if P<0, set P>0 
                self.params.planet_params[j+1] = -self.params.planet_params[j+1]

            if self.params.planet_params[j+2] < 0.0:  # if e<0, set e>0 and w=w+180, M0=M0-180
                self.params.planet_params[j+2] = -self.params.planet_params[j+2]
                self.params.planet_params[j+3] = self.params.planet_params[j+3] +  180.0
                self.params.planet_params[j+4] = self.params.planet_params[j+4] -  180.0
 
            if self.params.planet_params[j+2] > 0.99:
                self.params.planet_params[j+2]=0.99 
 
            if self.params.planet_params[j] < 0.0:  # if K<0, set K>0 and w = w+180.0
                self.params.planet_params[j+3] = self.params.planet_params[j+3] + 180.0
                self.params.planet_params[j] = -self.params.planet_params[j]

            # here we assure arg of periastron, mean anomaly and line of nodes to be between 0 and 360.0, and inclination between 0 and 180.0
            self.params.planet_params[j+3]  = np.fmod(self.params.planet_params[j+3],360.0) 
            self.params.planet_params[j+3]  = np.where(self.params.planet_params[j+3]<0,self.params.planet_params[j+3]+360.0, self.params.planet_params[j+3])
            self.params.planet_params[j+4]  = np.fmod(self.params.planet_params[j+4],360.0)  
            self.params.planet_params[j+4]  = np.where(self.params.planet_params[j+4]<0,self.params.planet_params[j+4]+360.0, self.params.planet_params[j+4])
            self.params.planet_params[j+5]  = np.fmod(self.params.planet_params[j+5],180.0) 
            self.params.planet_params[j+5]  = np.where(self.params.planet_params[j+5]<0,self.params.planet_params[j+5]+180.0, self.params.planet_params[j+5])
            self.params.planet_params[j+6]  = np.fmod(self.params.planet_params[j+6],360.0) 
            self.params.planet_params[j+6]  = np.where(self.params.planet_params[j+6]<0,self.params.planet_params[j+6]+360.0, self.params.planet_params[j+6])
                    
        return                     
        
    def mass_semimajor(self): # calculates planet masses (in Jupiter masses) and orbit semimajor axes (in astronomical units)
        # declaring some constants which will be useful here
        THIRD = 1.0/3.0
        GMSUN = 1.32712497e20
        AU=1.49597892e11
        # arrays for storing mass and semimajor information
        mass  = np.zeros(self.npl+1)
        ap    = np.zeros(self.npl)
        mtotal = self.params.stellar_mass
        mass[0] = self.params.stellar_mass
        f = 5e-6 # iteration step
        for i in range(self.npl):
            T = self.params.planet_params[7*i+1]*86400.0
            dm = 0 
            # we need innitial guess for each planet mass
            mass[i+1] = abs(self.params.planet_params[7*i])*(T*(self.params.stellar_mass)**2.0/(TAU*GMSUN))**THIRD * np.sqrt(1.0-self.params.planet_params[7*i+2]**2.0)/abs(np.sin(self.params.planet_params[7*i+5]))
            # mtotal is the sum of the stellar mass and masses of all planets computed so far, plus the current estimate for the mass of the planet we want to compute, we increase it as we go along
            mtotal = mtotal + mass[i+1]
            mpold = mass[i+1] # we need to keep track of the current estimate for the planet mass, to check if the new one is higher or not
           
            # This is a simple iteration to solve for planet mass
            while (dm <= 0): # if the new estimate for planet mass is lower than the old one, we finish an iteration
                mass[i+1] = abs(self.params.planet_params[7*i])*(T*mtotal**2.0/(TAU*GMSUN))**THIRD * np.sqrt(1.0-self.params.planet_params[7*i+2]**2.0)/abs(np.sin(self.params.planet_params[7*i+5])) # this will become the planet mass at the last run of the loop,until then it is only used for comparison in the while condition, whereas mpold is the estimate for planet mass. 
                dm = (mpold - mass[i+1]) 
                mpold =  mpold + f
                mtotal = mtotal + f 
                
            mtotal = mtotal-dm # the last part of the sum was mpold, now we want mass[m+1]
            ap[i] = (GMSUN*mtotal*(T/TAU)**2)**THIRD
        # 1 047.92612 = solar mass / jupiter mass, to get output in correct units   
        self.masses = mass[1:]*1047.92612
        self.semimajor = ap/AU
        return       
        
    def print_info(self, short_errors=True,show=True): 
        self.sort_by_period(reverse=False) # so the innermost planets will be first
        
        message_str = """ 
"""        
#        if(self.inputfile_read):
#            message_str = message_str +"""This is the information for signal fit %s, based on input file %s. 
            
#"""%(self.name, self.inputfile) 
#        else:
#            message_str = message_str +"""This is the information for signal fit %s.
            
#"""%self.name          
        # checking if we are dealing with fitted parameters or original data
#        if (self.fit_performed and not self.never_saved):
#            message_str = message_str +"""Presenting optimal parameters fitted using %s fitting method.
#            
#"""%self.fitting_method
#            message_str = message_str +"""Fit properties: \n chi^2: %f \n reduced chi^2: %f \n rms: %f \n loglik: %f
#            
#"""%(self.fit_results.chi2,self.fit_results.reduced_chi2,self.fit_results.rms,self.fit_results.loglik)
#        else:
#            message_str = message_str +"""No fit has yet been conducted (or parameters have never been saved), presenting parameters from user\'s original input
#
#"""

        message_str = message_str +"""Fit properties: \n chi^2: %f \n reduced chi^2: %f \n rms: %f \n loglik: %f
            
"""%(self.fit_results.chi2,self.fit_results.reduced_chi2,self.fit_results.rms,self.fit_results.loglik)


        # Now let's print information about RV files
        if (self.filelist.ndset==1): # because word file has a singular and plural form
            message_str = message_str +"""RV signal data was provided in 1 file. We expect it to have following offset and jitter:
            
"""
        else:
            message_str = message_str +"""RV signal data was provided in %d files. We expect these files to have following offsets and jitters:

"""%self.filelist.ndset
        if ((not (self.fit_performed) or self.fitting_method.startswith('loglik')) or self.never_saved): # if there was no fitting done (or values weren't saved) we don't give the errors. Same if the fitting was done with the Symplex method (if loglik was minimized), this method doesn't provide the errors
            for i in range(self.filelist.ndset): # neatly print all names of RV data files with corresponding offsets and jitters
                message_str = message_str +"""\n {0:15} \n offset: {1:>7.4f}  \n jitter: {2:>7.4f}            
""".format(self.filelist.files[i].name,self.params.offsets[i],self.params.jitters[i])
        else: # if the fit was done we provide the errors
            if (short_errors): # assume same + and - error, equal to maximum of the two
                for i in range(self.filelist.ndset): # neatly print all names of RV data files with corresponding offsets and jitters
                    message_str = message_str +"""\n {0:15} \n offset: {1:>7.4f} +/- {3:>7.4f} \n jitter: {2:>7.4f} +/- {4:>7.4f}            
""".format(self.filelist.files[i].name,self.params.offsets[i],self.params.jitters[i],max(self.param_errors.offset_errors[i]),max(self.param_errors.jitter_errors[i]))

            else:
                for i in range(self.filelist.ndset): # neatly print all names of RV data files with corresponding offsets and jitters
                    message_str = message_str +"""\n {0:15} \n offset: {1:>7.4f} + {3:>7.4f} - {5:>7.4f} \n jitter: {2:>7.4f} + {4:>7.4f} - {6:>7.4f}
""".format(self.filelist.files[i].name,self.params.offsets[i],self.params.jitters[i],self.param_errors.offset_errors[i][1],self.param_errors.jitter_errors[i][1],self.param_errors.offset_errors[i][0],self.param_errors.jitter_errors[i][0])                

        # Printing information about stellar activity, if fitting was done with GP
        if(self.fitting_method.startswith('GP')):
            message_str = message_str +"""\nStellar activity was modeled using Gaussian Processes. The resulting parameters are as follows
"""
            if(short_errors):
                for i in range(self.params.GP_params.npar):
                    message_str = message_str +"""\n {2:s} = {0:>7.4f} +/- {1:>7.4f}""".format(self.GP_params.gp_par[i],max(self.param_errors.GP_params_errors[i]),self.params.GP_params.rot_kernel.parameter_labels()[i])         
            else:
                for i in range(self.params.GP_params.npar):
                    message_str = message_str +"""\n {3:s} = {0:>7.4f} + {1:>7.4f} - {2:>7.4f}""".format(self.GP_params.gp_par[i],self.param_errors.GP_params_errors[i][1],self.param_errors.GP_params_errors[i][0],self.params.GP_params.rot_kernel.parameter_labels()[i])         

                          
        # Printing information about linear trend, if any
        if (self.params.linear_trend!=0): # Information about a linear trend only if it exists
            if ((not (self.fit_performed) or self.fitting_method.startswith('loglik')) or self.never_saved): # if there was no fitting done (or values weren't saved) we don't give the errors. Same if the fitting was done with the Symplex method (if loglik was minimized), this method doesn't provide the errors
                message_str = message_str +"""\nThere is a non-zero linear trend in the data, linear coefficient is expected to be equal %f.              
"""%self.params.linear_trend
            else:
                if (short_errors):
                    message_str = message_str +"""\nThere is a non-zero linear trend in the data, linear coefficient is expected to be equal {0:>7.4f} +/- {1:>7.4f}.                 
""".format(self.params.linear_trend,max(self.param_errors.linear_trend_error))
                else:
                    message_str = message_str +"""\nThere is a non-zero linear trend in the data, linear coefficient is expected to be equal {0:>7.4f} + {1:>7.4f} - {2:>7.4f}.
""".format(self.params.linear_trend,self.param_errors.linear_trend_error[1],self.param_errors.linear_trend_error[0])
 
        # Printing information about stellar mass 
        if (not (self.fit_performed) or not ((self.fitting_method=='mcmc_dyn' or self.fitting_method=='GP_dyn') and (self.use.use_stellar_mass))):
            if (self.stellar_mass_known):
                message_str = message_str +"""\nThe mass of the host star is {0:5.4f} solar masses.
""".format(self.params.stellar_mass)
            #else:
            #    message_str = message_str +"""\nThe mass of the host star is not known. Default value of %f is assumed.
#"""%DEFAULT_STELLAR_MASS
        else:
            if(short_errors):  
                message_str = message_str +"""\nThe mass of the host star is expected to be {0:5.4f} +/- {1:>7.4f} solar masses.
""".format(self.params.stellar_mass,max(self.param_errors.stellar_mass_error))       
            else:
                message_str = message_str +"""\nThe mass of the host star is expected to be  {0:>7.4f} + {1:>7.4f} - {2:>7.4f}.
""".format(self.params.stellar_mass,self.param_errors.stellar_mass_error[1],self.param_errors.stellar_mass_error[0])
 
        # Printing planet parameters:
        if (self.npl==0): # if there are zero planets
             message_str = message_str +"""\nThe system has no planets. 
             """
        elif (self.npl==1): # because word planet has a singular and plural form
             message_str = message_str +"""\nThe system has 1 planet. 
             """
        else:
             message_str = message_str +"""\nThe system has %d planets. 
"""%self.npl
        
        if (self.npl>0):
            if not (self.stat_saved):
                self.mass_semimajor() # update mass and semimajor axes values, if they weren't taken from the fortran code 
            message_str = message_str +"""\nKnown planets are expected to have following properties (mean anomalies for epoch {0:7.2f}):
""".format(self.epoch)
            if ((not (self.fit_performed) or self.fitting_method.startswith('loglik')) or self.never_saved): # again, if no fit or loglik fit, no errors
                for i in range(self.npl):
                    message_str = message_str +"""\n Planet {0:2d}: \n signal semiamplitude = {1:5.4f} m/s \n period = {2:5.4f} days \n orbit eccentricity = {3:5.4f} \n argument of periastron = {4:5.4f} \n mean anomally = {5:5.4f} \n inclination = {6:5.4f} \n line of nodes = {7:5.4f} \n mass  = {8:5.4f} M_Jup \n orbit semimajor axis = {9:5.4f} AU
""".format(i+1,self.params.planet_params[7*i],self.params.planet_params[7*i+1],self.params.planet_params[7*i+2],self.params.planet_params[7*i+3],self.params.planet_params[7*i+4],self.params.planet_params[7*i+5],self.params.planet_params[7*i+6],self.masses[i],self.semimajor[i])
            else:
                if (short_errors):
                    for i in range(self.npl):
                        message_str = message_str +"""\n Planet {0:2d}: \n signal semiamplitude = {1:5.4f} +/- {10:5.4f} m/s \n period = {2:5.4f} +/- {11:5.4f}  days \n orbit eccentricity = {3:5.4f} +/- {12:5.4f} \n argument of periastron = {4:5.4f} +/- {13:5.4f} \n mean anomally = {5:5.4f} +/- {14:5.4f} \n inclination = {6:5.4f} +/- {15:5.4f} \n line of nodes = {7:5.4f} +/- {16:5.4f} \n mass  = {8:5.4f} M_Jup \n orbit semimajor axis = {9:5.4f} AU

""".format(i+1,self.params.planet_params[7*i],self.params.planet_params[7*i+1],self.params.planet_params[7*i+2],self.params.planet_params[7*i+3],self.params.planet_params[7*i+4],self.params.planet_params[7*i+5],self.params.planet_params[7*i+6],self.masses[i],self.semimajor[i],max(self.param_errors.planet_params_errors[7*i]),max(self.param_errors.planet_params_errors[7*i+1]),max(self.param_errors.planet_params_errors[7*i+2]),max(self.param_errors.planet_params_errors[7*i+3]),max(self.param_errors.planet_params_errors[7*i+4]),max(self.param_errors.planet_params_errors[7*i+5]),max(self.param_errors.planet_params_errors[7*i+6]))
                else:
                    for i in range(self.npl):
                        message_str = message_str +"""\n Planet {0:2d} \n signal semiamplitude = {1:5.4f} + {10:5.4f} - {17:5.4f} m/s \n period = {2:5.4f} + {11:5.4f} - {18:5.4f} days \n orbit eccentricity = {3:5.4f} + {12:5.4f} - {19:5.4f} \n argument of periastron = {4:5.4f} + {13:5.4f} - {20:5.4f} \n mean anomally = {5:5.4f} + {14:5.4f} -{21:5.4f} \n inclination = {6:5.4f} + {15:5.4f} - {22:5.4f} \n line of nodes = {7:5.4f} + {16:5.4f} - {23:5.4f} \n mass  = {8:5.4f} M_Jup \n orbit semimajor axis = {9:5.4f} AU

""".format(i+1,self.params.planet_params[7*i],self.params.planet_params[7*i+1],self.params.planet_params[7*i+2],self.params.planet_params[7*i+3],self.params.planet_params[7*i+4],self.params.planet_params[7*i+5],self.params.planet_params[7*i+6],self.masses[i],self.semimajor[i],self.param_errors.planet_params_errors[7*i][1],self.param_errors.planet_params_errors[7*i+1][1],self.param_errors.planet_params_errors[7*i+2][1],self.param_errors.planet_params_errors[7*i+3][1],self.param_errors.planet_params_errors[7*i+4][1],self.param_errors.planet_params_errors[7*i+5][1],self.param_errors.planet_params_errors[7*i+6][1],self.param_errors.planet_params_errors[7*i][0],self.param_errors.planet_params_errors[7*i+1][0],self.param_errors.planet_params_errors[7*i+2][0],self.param_errors.planet_params_errors[7*i+3][0],self.param_errors.planet_params_errors[7*i+4][0],self.param_errors.planet_params_errors[7*i+5][0],self.param_errors.planet_params_errors[7*i+6][0])
 
        if show:
            print(message_str)
            return
        return message_str  
        
  
    def fortran_input(self, program='chi2_kep', fileinput=False, filename='debug_input', amoeba_starts=0, outputfiles=[1,1,1],eps='1.0E-8',dt=864000, when_to_kill=300, npoints=50, model_max = 100, model_min =0): # generate input string for the fortran code, optionally as a file


        ### ppp will be the input string. Depending on fileinput parameter we either save it in a file or save it directly 
     
       # print(program, self.mixed_fit[1])
         
        if not (fileinput): # if we want to save input in a file we don't want this line in the input string    
            ppp = '%s << EOF\n'%program
        else:
            ppp = '' 
        ppp+= '%s %f %d %d %d %d %d %d\n%f %d %d %d \n%d\n'%(eps,dt,amoeba_starts,when_to_kill,npoints, model_max, model_min, self.gr_flag, self.params.stellar_mass,outputfiles[0], outputfiles[1],outputfiles[2],self.filelist.ndset) # first three lines of fortran input: precision and timestep for integration, stellar mass and number of datasets
        for i in range(self.filelist.ndset): 
            # path for each dataset      
            ppp+='%s\n'%(self.filelist.files[i].path)    
            # offset and jitter information for each dataset
            ppp+='%f\n%d\n'%(self.params.offsets[i],int(self.use.use_offsets[i]))
            ppp+='%f\n%d\n'%(self.params.jitters[i],int(self.use.use_jitters[i]))  
        ppp+='%d\n'%self.npl
        
        # if mixed fitting is requested        
        if program == '%s/lib/fr/loglik_dyn+'%(self.cwd):       
            for i in range(self.npl):             
                ppp+='%d\n'%self.mixed_fit[1][i] 
                
        for i in range(self.npl): # K,P,e,w,M,i,cap0m for each planet, and information which ones we use
            if self.hkl:
                ppp+='%f %.8f %.8f %f %f %f %f %f\n'%(self.params.planet_params[7*i],self.params.planet_params[7*i+1],self.e_sinw[i],self.e_cosw[i],self.lamb[i],self.params.planet_params[7*i+5],self.params.planet_params[7*i+6],self.omega_dot[i])
                ppp+='%d %d %d %d %d %d %d %d\n'%(int(self.use.use_planet_params[7*i]),int(self.use.use_planet_params[7*i+1]),int(self.use.use_planet_params[7*i+2]),int(self.use.use_planet_params[7*i+3]),int(self.use.use_planet_params[7*i+4]),int(self.use.use_planet_params[7*i+5]),int(self.use.use_planet_params[7*i+6]),int(self.omega_dot_use[i]))     
            else:
                ppp+='%f %.8f %.8f %f %f %f %f %f\n'%(self.params.planet_params[7*i],self.params.planet_params[7*i+1],self.params.planet_params[7*i+2],self.params.planet_params[7*i+3],self.params.planet_params[7*i+4],self.params.planet_params[7*i+5],self.params.planet_params[7*i+6],self.omega_dot[i])
                ppp+='%d %d %d %d %d %d %d %d\n'%(int(self.use.use_planet_params[7*i]),int(self.use.use_planet_params[7*i+1]),int(self.use.use_planet_params[7*i+2]),int(self.use.use_planet_params[7*i+3]),int(self.use.use_planet_params[7*i+4]),int(self.use.use_planet_params[7*i+5]),int(self.use.use_planet_params[7*i+6]),int(self.omega_dot_use[i]))                                
                
        ppp+='%.15f\n%d\n'%(self.params.linear_trend,int(self.use.use_linear_trend)) # information about linear trend
        ppp+='%.15f\n%d\n'%(self.rv_quadtr,int(bool(self.rv_quadtr_use))) # information about linear trend
              
        ppp+='%.5f\n'%self.epoch
        ppp+='%s\n'%int(self.hkl)  
        
       # if program == '%s/lib/fr/loglik_dyn+'%(self.cwd):       
#        print(int(self.hkl),self.npl)    
        # prepare final version of the ppp command to be returned by this function
        if not (fileinput):
            ppp+='EOF' # end of the command to run in the case of saving input directly
        else: # here's what we do if we want to generate a file as well 
            # first we save the ppp string in a file (by default 'Kep_input')
            file_kep = open("%s/lib/fr/%s"%(self.cwd,filename), 'w')
            file_kep.write('%s'%ppp) 
            file_kep.close()
            # then we overwrite ppp with the command to pass this file as input for the fortran code
            ppp='./%s < %s'%(program,filename)
        
        #print(ppp)
        return ppp 

    # sort planets by one of the parameters (K,P,e,w,M0)
    def sort_by_param(self,param,reverse=True):
          if (self.npl>0): # no point in sorting if there's no planets
            sort = convert_array_to_int(np.array(sorted(range(self.npl), key=lambda k: self.params.planet_params[7*k+param], reverse=reverse))) # find a permutation sorting planets by chosen parameter
            # permutate all parameters which refer to planets, according to the permutation found above
            planets_to_sort=np.array([self.params.planet_params[n:n+7] for n in range(0,7*self.npl,7)]) # group parameters for a single planet together for sorting 
            planets_to_sort=planets_to_sort[sort] # sort grouped array according to the found permutation
            self.params=parameters(self.params.offsets,self.params.jitters,np.concatenate(planets_to_sort),self.params.linear_trend,self.params.stellar_mass) # now come back to 1d array      
            # now do the same with use flags
            use_planets_to_sort=np.array([self.use.use_planet_params[n:n+7] for n in range(0,7*self.npl,7)]) # group parameters for a single planet together for sorting
            use_planets_to_sort=use_planets_to_sort[sort] # sort grouped array according to the found permutation
            self.use=use_flags(self.use.use_offsets,self.use.use_jitters,np.concatenate(use_planets_to_sort),self.use.use_linear_trend,self.use.use_stellar_mass) # now come back to 1d array
            # now do the same with errors, only if they are established (fit performed)
            if (self.fit_performed):
                planet_errors_to_sort=np.array([self.param_errors.planet_params_errors[n:n+7] for n in range(0,7*self.npl,7)]) # group parameters for a single planet together for sorting
                planet_errors_to_sort=planet_errors_to_sort[sort] # sort grouped array according to the found permutation
                self.param_errors=parameter_errors([i[0] for i in self.param_errors.offset_errors],[i[0] for i in self.param_errors.jitter_errors],[i[0] for i in np.concatenate(planet_errors_to_sort)],self.param_errors.linear_trend_error[0],self.param_errors.stellar_mass_error[0]) # now come back to 1d array
          return 
                
    def sort_by_semiamplitude(self,reverse=True):
        self.sort_by_param(0,reverse=reverse) 
        return 

    def sort_by_period(self,reverse=True):
        self.sort_by_param(1,reverse=reverse) 
        return 

    def sort_by_eccentricity(self,reverse=True):
        self.sort_by_param(2,reverse=reverse) 
        return 

    def sort_by_arg_of_periastron(self,reverse=True):
        self.sort_by_param(3,reverse=reverse) 
        return 

    def sort_by_mean_anomaly(self,reverse=True):
        self.sort_by_param(4,reverse=reverse) 
        return 
 
    # update planet parameters of signal_fit etc. with fit results.
    def update_with_fit_results(self):
        self.chi2=self.fit_results.chi2
        self.reduced_chi2=self.fit_results.reduced_chi2
        self.rms=self.fit_results.rms
        self.loglik=self.fit_results.loglik     
        if not (self.fit_performed and self.stat_saved): # if we haven't run a fit or didn't save parameters, we can do no more
            pass
        else: # update parameters and their errors based on the recent fit
            self.params=self.fit_results.stat.params
            self.param_errors=self.fit_results.stat.param_errors
            self.masses=self.fit_results.mass
            self.semimajor=self.fit_results.a
           # print(self.fit_results.omega_dot,self.fit_results.omega_dot_err)
            self.omega_dot = self.fit_results.omega_dot
            self.omega_dot_err = self.fit_results.omega_dot_err
            self.rv_quadtr = self.fit_results.rv_quadtr
            self.rv_quadtr_err = self.fit_results.rv_quadtr_err            
            
        return
        
         
    ### this function is a wrapper calling a fortran program to fit parameters in keplerian mode by minimizing chi^2   WORK IN PROGRESS ON THAT ONE! 
        
    def fitting(self, minimize_fortran=True, minimize_loglik=False, fileinput=False, doGP=False, kernel_id=-1, filename='Kep_input', outputfiles=[1,1,1], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 5,model_min =0): # run the fit which will either minimize chi^2 or loglik.
        '''       
         eps, dt - accuracy and step for the integration in dynamical case, in keplerian case ignored, just given for input consistency
         which value to minimize (used for calling an appropriate fortran program)
        '''      
        
        eps = eps * 1e-13 
        dt  = dt  * 86400.0        
        
        check_temp_RV_file(self)

        # bug fix. TBF
        self.bound_error = False
        self.bound_error_msg = ""
        #########################

        if(minimize_loglik):
            minimized_value='loglik'
        else:
            minimized_value='chi2'
        # which mode to use, analogously as above
             
        if(self.mod_dynamical):
            if self.mixed_fit[0] == True and minimized_value=='loglik':            
                mod='dyn+'
            else:
                mod='dyn'                
        else:
            mod='kep'
        #print(mod, self.mixed_fit)
        if minimize_fortran == True and doGP ==False:   
             program='%s/lib/fr/%s_%s'%(self.cwd,minimized_value,mod) 
             #print(program)
             text,flag=run_command_with_timeout(self.fortran_input(program=program, fileinput=fileinput, filename=filename, outputfiles=outputfiles,amoeba_starts=amoeba_starts,eps=eps,dt=dt, when_to_kill=fortran_kill, npoints=npoints, model_max = model_max, model_min =model_min), timeout_sec, output=True,pipe=(not bool(outputfiles[2]))) # running command generated by the fortran_input function 
           
        else:
            #self.fitting_SciPyOp(doGP=doGP, gp_par=gp_par, kernel_id=kernel_id, use_gp_par=use_gp_par)  
            self = run_SciPyOp(self, kernel_id=kernel_id)           
            program='%s/lib/fr/%s_%s'%(self.cwd,minimized_value,mod) 
            text,flag=run_command_with_timeout(self.fortran_input(program=program, fileinput=fileinput, filename=filename, outputfiles=outputfiles,amoeba_starts=0,eps=eps,dt=dt, when_to_kill=fortran_kill, npoints=npoints, model_max = model_max, model_min =model_min), timeout_sec, output=True,pipe=(not bool(outputfiles[2]))) # running command generated by the fortran_input function 
            #print(text)
            
        if (flag==1):
            fortranoutput=fortran_output(text,self.npl,self.filelist.ndset,self.params.stellar_mass) # create an object for fortran output, which we need to split into part
            self.fit_results=fortranoutput.modfit(print_stat=print_stat)
            
            #print(self.fit_results.omega_dot)
            
            
            
            self.stat_saved=self.fit_results.stat_array_saved
            if (self.stat_saved):
                self.never_saved=False
            self.model_saved=bool(outputfiles[2])
            self.fit_performed=True
            if(self.fit_results.stat_array_saved): 
                self.fitting_method=program
            self.update_with_fit_results()
            self.correct_elements() #because amoeba might make things wrong here
            
#        if self.rtg[1]:
#            get_gps_model(self)  
            
        ##################### New stuff to be added here ######################    
        for i in range(self.npl):
            self.t_peri[i] = (float(self.epoch) - (np.radians(self.params.planet_params[7*i + 4])/(2*np.pi))*self.params.planet_params[7*i + 1] ) # epoch  - ((ma/TWOPI)*a[1])
          
        #####################               
        self.hack_around_rv_params()  

            
        if (return_flag):
            return flag
        else:
            return

  
    # if you want to temporarily change use flags and run a new fit, use this wrapper
    def quick_overwrite_use_and_fit(self,useflags,minimize_loglik=False, fileinput=False, filename='Kep_input', outputfiles=[1,1,0], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 500, model_min=0):
        oldflags=self.overwrite_use(useflags,save=True)
        self.fitting(minimize_loglik=minimize_loglik,fileinput=fileinput,filename=filename,
        amoeba_starts=amoeba_starts,outputfiles=outputfiles, eps=eps,dt=dt,timeout_sec=timeout_sec,print_stat=print_stat, fortran_kill=fortran_kill )     
        self.overwrite_use(oldflags,save=False)
        return
        
    # if you want to temporarily change initial parameters and run a new fit, use this wrapper
    def quick_overwrite_params_and_fit(self,params,minimize_loglik=True, fileinput=False, filename='Kep_input',outputfiles=[1,1,0], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 500, model_min=0):
        oldparams=self.overwrite_params(params,save=True)
        if (return_flag):
            flag=self.fitting(minimize_loglik=minimize_loglik,fileinput=fileinput,filename=filename,amoeba_starts=amoeba_starts,outputfiles=outputfiles,eps=eps,dt=dt,timeout_sec=timeout_sec,print_stat=print_stat, return_flag=True, fortran_kill=fortran_kill)     
            self.overwrite_params(oldparams,save=False)
            return flag
        else:
            self.fitting(minimize_loglik=minimize_loglik,fileinput=fileinput,filename=filename,amoeba_starts=amoeba_starts,outputfiles=outputfiles,eps=eps,dt=dt,timeout_sec=timeout_sec,print_stat=print_stat, fortran_kill=fortran_kill)     
            self.overwrite_params(oldparams,save=False)
            return                                   
            
    def minimize_one_param_K(self,planet,minimize_loglik=True, fileinput=False, filename='Kep_input', outputfiles=[1,0,0], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 500, model_min =0):
        useflags=use_flags([False]*20,[False]*20,[False]*70,False,False) 
        useflags.update_use_K(planet,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat, fortran_kill=fortran_kill)
        return    

    def minimize_one_param_P(self,planet,minimize_loglik=True, fileinput=False, filename='Kep_input', outputfiles=[1,0,0], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 500, model_min =0):
        useflags=use_flags([False]*20,[False]*20,[False]*70,False,False) 
        useflags.update_use_P(planet,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat, fortran_kill=fortran_kill)
        return    
        
    def minimize_one_param_e(self,planet,minimize_loglik=True, fileinput=False, filename='Kep_input', outputfiles=[1,0,0], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 500, model_min =0):
        useflags=use_flags([False]*20,[False]*20,[False]*70,False,False) 
        useflags.update_use_e(planet,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat, fortran_kill=fortran_kill)
        return    
                
    def minimize_one_param_w(self,planet,minimize_loglik=True, fileinput=False, filename='Kep_input', outputfiles=[1,0,0], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 500, model_min =0):
        useflags=use_flags([False]*20,[False]*20,[False]*70,False,False) 
        useflags.update_use_w(planet,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat, fortran_kill=fortran_kill)
        return    
        
    def minimize_one_param_M0(self,planet,minimize_loglik=True, fileinput=False, filename='Kep_input', outputfiles=[1,0,0], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 500, model_min =0):
        useflags=use_flags([False]*20,[False]*20,[False]*70,False,False) 
        useflags.update_use_M0(planet,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat, fortran_kill=fortran_kill)
        return    
        
    def minimize_one_param_inclination(self,planet,minimize_loglik=True, fileinput=False, filename='Kep_input', outputfiles=[1,0,0], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 500, model_min =0):
        useflags=use_flags([False]*20,[False]*20,[False]*70,False,False) 
        useflags.update_use_inclination(planet,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat, fortran_kill=fortran_kill)
        return    
        
    def minimize_one_param_lineofnodes(self,planet,minimize_loglik=True, fileinput=False, filename='Kep_input', outputfiles=[1,0,0], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 500, model_min =0):
        useflags=use_flags([False]*20,[False]*20,[False]*70,False,False) 
        useflags.update_use_lineofnodes(planet,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat, fortran_kill=fortran_kill)
        return    
             
    def minimize_one_param_offset(self,dataset,minimize_loglik=True, fileinput=False, filename='Kep_input', outputfiles=[1,0,0], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 500, model_min =0):
        useflags=use_flags([False]*20,[False]*20,[False]*70,False,False) 
        useflags.update_use_offset(dataset,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat, fortran_kill=fortran_kill)
        return                 
  
    def minimize_one_param_jitter(self,dataset,minimize_loglik=True, fileinput=False, filename='Kep_input', outputfiles=[1,0,0], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 500, model_min =0):
        useflags=use_flags([False]*20,[False]*20,[False]*70,False,False) 
        useflags.update_use_jitter(dataset,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat, fortran_kill=fortran_kill)
        return                   
              
    def minimize_one_param_linear_trend(self,minimize_loglik=True, fileinput=False, filename='Kep_input', outputfiles=[1,0,0], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 500, model_min =0):
        useflags=use_flags([False]*20,[False]*20,[False]*70,False,False) 
        useflags.update_use_linear_trend(True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat, fortran_kill=fortran_kill)
        return       
                    
#### prep. #########
    def prepare_for_mcmc(self, rtg=[True,False,False,False], customdatasetlabels=[]):
        '''Prepare bounds and par array needed for the MCMC''' 
 
        preparingwarnings=Warning_log([],'Preparing for MCMC')  
        # put together bounds for K,P,e,w,M0, so they are in an order we are used to
        
        rtg=self.rtg
        
        #print(rtg)
    
        par = []
        flag = []
        par_str = []
        bounds = []
        prior_nr = []
        prior_jeff = []
  
 
        for i in range(self.filelist.ndset):           
            par.append(self.params.offsets[i]) #
            par_str.append(self.rvoff_str[i])
            bounds.append(self.rvoff_bounds[i])        
            prior_nr.append(self.rvoff_norm_pr[i])
            prior_jeff.append(self.rvoff_jeff_pr[i])
           
            if rtg == [False,False,True,True]:
                flag.append(False) #
            elif rtg == [False,False,True,False]:
                flag.append(False) # 
            else:   
                flag.append(self.use.use_offsets[i]) #

            
        for i in range(self.filelist.ndset):              
            par.append(self.params.jitters[i]) #
            par_str.append(self.jitt_str[i]) #
            bounds.append(self.jitt_bounds[i])   
            prior_nr.append(self.jitt_norm_pr[i])
            prior_jeff.append(self.jitt_jeff_pr[i])
            
            if rtg == [False,False,True,True]:
                flag.append(False) #
            elif rtg == [False,False,True,False]:
                flag.append(False) #   
            else:   
                flag.append(self.use.use_jitters[i])
 
        
        for i  in range(self.npl):
            
            par.append(self.params.planet_params[7*i]) #           
            par.append(self.params.planet_params[7*i+1]) #
            par.append(self.params.planet_params[7*i+2]) #
            par.append(self.params.planet_params[7*i+3]) #
            par.append(self.params.planet_params[7*i+4]) #
            par.append(self.params.planet_params[7*i+5]) #
            par.append(self.params.planet_params[7*i+6]) #
            
            if rtg == [False,False,True,True]:
                flag.append(False) #
            elif rtg == [False,False,True,False]:
                flag.append(False) # 
            else:
                flag.append(self.use.use_planet_params[7*i])
            flag.append(self.use.use_planet_params[7*i+1])
            flag.append(self.use.use_planet_params[7*i+2])
            flag.append(self.use.use_planet_params[7*i+3])
            
            if rtg == [False,False,True,True]:
                flag.append(False) #
            elif rtg == [False,False,True,False]:
                flag.append(False) # 
            else:
                flag.append(self.use.use_planet_params[7*i+4])
                
            flag.append(self.use.use_planet_params[7*i+5])
            
            if rtg == [False,False,True,True]:
                flag.append(False) #
            elif rtg == [False,False,True,False]:
                flag.append(False) # 
            else:
                flag.append(self.use.use_planet_params[7*i+6])
            
            
            par_str.append(self.K_str[i])
            par_str.append(self.P_str[i])         
            
            bounds.append(self.K_bound[i])        
            bounds.append(self.P_bound[i]) 

            prior_nr.append(self.K_norm_pr[i])
            prior_nr.append(self.P_norm_pr[i])
           
            prior_jeff.append(self.K_jeff_pr[i])
            prior_jeff.append(self.P_jeff_pr[i])

               
            if self.hkl == False:
                par_str.append(self.e_str[i])
                par_str.append(self.w_str[i])
                par_str.append(self.M0_str[i])            
    
                bounds.append(self.e_bound[i]) 
                bounds.append(self.w_bound[i])  
                bounds.append(self.M0_bound[i])
    
                prior_nr.append(self.e_norm_pr[i])
                prior_nr.append(self.w_norm_pr[i])
                prior_nr.append(self.M0_norm_pr[i]) 
                
                prior_jeff.append(self.e_jeff_pr[i])
                prior_jeff.append(self.w_jeff_pr[i])
                prior_jeff.append(self.M0_jeff_pr[i])    
            else:   
                par_str.append(self.e_sinw_str[i])
                par_str.append(self.e_cosw_str[i])
                par_str.append(self.lamb_str[i])            
    
                bounds.append(self.e_sinw_bound[i]) 
                bounds.append(self.e_cosw_bound[i])  
                bounds.append(self.lamb_bound[i])
    
                prior_nr.append(self.e_sinw_norm_pr[i])
                prior_nr.append(self.e_cosw_norm_pr[i])
                prior_nr.append(self.lamb_norm_pr[i]) 
                
                prior_jeff.append(self.e_sinw_jeff_pr[i])
                prior_jeff.append(self.e_cosw_jeff_pr[i])
                prior_jeff.append(self.lamb_jeff_pr[i])              
        


            par_str.append(self.i_str[i])
            par_str.append(self.Node_str[i] )           

            bounds.append(self.i_bound[i])
            bounds.append(self.Node_bound[i])

            prior_nr.append(self.i_norm_pr[i])
            prior_nr.append(self.Node_norm_pr[i]) 
            
            prior_jeff.append(self.i_jeff_pr[i])
            prior_jeff.append(self.Node_jeff_pr[i])



        par.append(self.params.linear_trend)
        #flag.append(self.use.use_linear_trend)
        par_str.append(self.rv_lintr_str[0])
        bounds.append(self.rv_lintr_bounds[0])
        prior_nr.append(self.rv_lintr_norm_pr[0])
        prior_jeff.append(self.rv_lintr_jeff_pr[0])
        
        
        par.append(self.rv_quadtr)
       # flag.append(self.rv_quadtr_use)
        par_str.append(self.rv_quadtr_str[0])
        bounds.append(self.rv_quadtr_bounds[0])
        prior_nr.append(self.rv_quadtr_norm_pr[0])
        prior_jeff.append(self.rv_quadtr_jeff_pr[0])   
        
        if rtg == [False,False,True,True]:
            flag.append(False) #
            flag.append(False) #            
        elif rtg == [False,False,True,False]:
            flag.append(False) # 
            flag.append(False) #            
        else:   
            flag.append(self.use.use_linear_trend) #        
            flag.append(self.rv_quadtr_use)

       
        if rtg[1] == True:
            if self.gp_kernel == 'RotKernel':
                for i in range(4):  
                    par.append(self.GP_rot_params[i])
                    flag.append(self.GP_rot_use[i])
                    par_str.append(self.GP_rot_str[i])
                    bounds.append(self.GP_rot_bounds[i])
                    prior_nr.append(self.GP_rot_norm_pr[i])
                    prior_jeff.append(self.GP_rot_jeff_pr[i])
                
            elif self.gp_kernel == 'SHOKernel':         
                for i in range(3):  
                    par.append(self.GP_sho_params[i])
                    flag.append(self.GP_sho_use[i])
                    par_str.append(self.GP_sho_str[i])
                    bounds.append(self.GP_sho_bounds[i])
                    prior_nr.append(self.GP_sho_norm_pr[i])
                    prior_jeff.append(self.GP_sho_jeff_pr[i])
           
 
        for i  in range(self.npl):            
            par.append(self.t0[i])
            par.append(self.pl_a[i])
            par.append(self.pl_rad[i])
                        
            par_str.append(self.t0_str[i])
            par_str.append(self.pl_a_str[i])
            par_str.append(self.pl_rad_str[i])
            
            bounds.append(self.t0_bound[i])
            bounds.append(self.pl_a_bound[i])
            bounds.append(self.pl_rad_bound[i])   
            
            prior_nr.append(self.t0_norm_pr[i])
            prior_nr.append(self.pl_a_norm_pr[i])
            prior_nr.append(self.pl_rad_norm_pr[i])   
            
            prior_jeff.append(self.t0_jeff_pr[i])
            prior_jeff.append(self.pl_a_jeff_pr[i])
            prior_jeff.append(self.pl_rad_jeff_pr[i])               
            
            if rtg[2] == [False]:
                flag.append(False) #
                flag.append(False) #
                flag.append(False)           
            else:
                flag.append(self.t0_use[i])
                flag.append(self.pl_a_use[i])
                flag.append(self.pl_rad_use[i])   
                
              
                
        for i in range(10):
            if len(self.tra_data_sets[i]) != 0:
                
                par.append(self.tra_off[i]) #
                par_str.append(self.tra_off_str[i])
                bounds.append(self.tra_off_bounds[i])        
                prior_nr.append(self.tra_off_norm_pr[i])
                prior_jeff.append(self.tra_off_jeff_pr[i])
                
                if rtg == [True, False,False,True]:
                    flag.append(False) #
                elif rtg == [True,False,False,False]:
                    flag.append(False) #
                else:   
                    flag.append(self.tra_off_use[i]) #


        for i in range(10):         
            if len(self.tra_data_sets[i]) != 0:
                par.append(self.tra_jitt[i]) #
                par_str.append(self.tra_jitt_str[i]) #
                bounds.append(self.tra_jitt_bounds[i])   
                prior_nr.append(self.tra_jitt_norm_pr[i])
                prior_jeff.append(self.tra_jitt_jeff_pr[i])
                
                
                if rtg == [True, False,False,True]:
                    flag.append(False) #
                elif rtg == [True,False,False,False]:
                    flag.append(False) #
                else:   
                    flag.append(self.tra_jitt_use[i])                  
                
               

        if rtg[3] == True:
            if self.tra_gp_kernel == 'RotKernel':
                for i in range(4):  
                    par.append(self.tra_GP_rot_params[i])
                    flag.append(self.tra_GP_rot_use[i])
                    par_str.append(self.tra_GP_rot_str[i])
                    bounds.append(self.tra_GP_rot_bounds[i])
                    prior_nr.append(self.tra_GP_rot_norm_pr[i])
                    prior_jeff.append(self.tra_GP_rot_jeff_pr[i])
                
            elif self.tra_gp_kernel == 'SHOKernel':         
                for i in range(3):  
                    par.append(self.tra_GP_sho_params[i])
                    flag.append(self.tra_GP_sho_use[i])
                    par_str.append(self.tra_GP_sho_str[i])
                    bounds.append(self.tra_GP_sho_bounds[i])
                    prior_nr.append(self.tra_GP_sho_norm_pr[i])
                    prior_jeff.append(self.tra_GP_sho_jeff_pr[i])
           

 
        for i  in range(self.npl):  
            par.append(self.omega_dot[i]) #
            flag.append(self.omega_dot_use[i])
            par_str.append(self.omega_dot_str[i]) #
            bounds.append(self.omega_dot_bounds[i])   
            prior_nr.append(self.omega_dot_norm_pr[i])
            prior_jeff.append(self.omega_dot_jeff_pr[i])


                
        par.append(self.params.stellar_mass)
        flag.append(self.use.use_stellar_mass)
        par_str.append(self.st_mass_str[0])
        bounds.append(self.st_mass_bounds[0])
        prior_nr.append(self.st_mass_norm_pr[0])   
        prior_jeff.append(self.st_mass_jeff_pr[0])   
 
        #print(par)    
#        print(flag)    
       # print(par_str)    
       # print(bounds)    
        
 
        self.f_for_mcmc = [idx for idx in range(len(flag)) if flag[idx] ==1 ] # indices for fitted parameters

        self.par_for_mcmc = []  # self par_for_mcmc are the fitted parameters   
        self.e_for_mcmc = [] # labels for fitted parameters only
        self.b_for_mcmc = [] # labels for fitted parameters only
        self.nr_pr_for_mcmc = [] # labels for fitted parameters only
        self.jeff_pr_for_mcmc = [] # labels for fitted parameters only
       
        self.parameters = par

       # print(el_str)

        for j in range(len(par)):
            #print(flag[j])
            if flag[j] > 0:
                self.par_for_mcmc.append(float(par[j]))
                self.e_for_mcmc.append(par_str[j])
                self.b_for_mcmc.append(bounds[j])
                self.nr_pr_for_mcmc.append(prior_nr[j])
                self.jeff_pr_for_mcmc.append(prior_jeff[j])
        
        
        
        #np.asarray(self.par_for_mcmc)
        #np.asarray(self.f_for_mcmc)
        #np.asarray(self.e_for_mcmc)
        #np.asarray(self.b_for_mcmc)
        #np.asarray(self.nr_pr_for_mcmc)
        #np.asarray(self.jeff_pr_for_mcmc)
        # self.par_for_mcmc = np.array(self.par_for_mcmc )
      #  self.f_for_mcmc = np.array(self.f_for_mcmc )
       # print(self.b_for_mcmc)
        self.par_for_mcmc = np.array(self.par_for_mcmc)
        self.f_for_mcmc = np.array(self.f_for_mcmc)
        self.e_for_mcmc = np.array(self.e_for_mcmc)
        self.b_for_mcmc = np.array(self.b_for_mcmc)
        self.nr_pr_for_mcmc = np.array(self.nr_pr_for_mcmc)
        self.jeff_pr_for_mcmc = np.array(self.jeff_pr_for_mcmc)
 
        preparingwarnings.print_warning_log()
        return                  
                      
  

    def verify_params_with_bounds(self):

        '''verify if all planet parameters are within allowed bounds'''     
        verification=True
        if not (verify_array_with_bounds(self.params.offsets,self.bounds.offset_bounds)):
            verification = False 
        elif not (verify_array_with_bounds(self.params.jitters,self.bounds.jitter_bounds)):
            verification = False
        elif not (verify_array_with_bounds(self.params.planet_params,self.bounds.planet_params_bounds)):
            verification = False
        elif not (verify_array_with_bounds(self.params.GP_params.gp_par,self.bounds.GP_params_bounds)):
            verification = False
        elif not (verify_array_with_bounds(np.atleast_1d(self.params.linear_trend),self.bounds.linear_trend_bounds)):
            verification = False 
        elif not (verify_array_with_bounds(np.atleast_1d(self.params.stellar_mass),self.bounds.stellar_mass_bounds)):
            verification = False 

        return verification           
        
    def generate_newparams_for_mcmc(self,p):
        newparams=self.params # we will now modify this object to represent values of new parameters 
        # now update these parameters which should be updated (use flag True, their indices are indicated by f_for_mcmc array)
        i=0
        for idx in self.f_for_mcmc:
            if (idx<self.filelist.ndset):
                newparams.update_offset(idx,p[i])           
                i=i+1             
            elif (idx<2*self.filelist.ndset):
                newparams.update_jitter(idx-self.filelist.ndset,p[i])
                i=i+1    
            elif (idx<2*self.filelist.ndset+7*self.npl):
                nr=idx-2*self.filelist.ndset
                x=int(nr/7)
                if (np.mod(nr,7)==0):
                    newparams.update_K(x,p[i])
                    i=i+1  
                elif (np.mod(nr,7)==1):
                    newparams.update_P(x,p[i])
                    i=i+1                           
                elif (np.mod(nr,7)==2):
                    newparams.update_e(x,p[i])
                    i=i+1       
                elif (np.mod(nr,7)==3):
                    newparams.update_w(x,p[i])
                    i=i+1       
                elif (np.mod(nr,7)==4):
                    newparams.update_M0(x,p[i])
                    i=i+1       
                elif (np.mod(nr,7)==5):
                    newparams.update_inclination(x,p[i])
                    i=i+1  
                elif (np.mod(nr,7)==6):
                    newparams.update_lineofnodes(x,p[i])
                    i=i+1  
            elif (idx<2*self.filelist.ndset+7*self.npl+1):
                newparams.update_linear_trend(p[i]) 
                i=i+1       
            elif (idx<2*self.filelist.ndset+7*self.npl+1+self.params.GP_params.npar):
                newparams.update_GP_param_value(idx-2*self.filelist.ndset-7*self.npl-1,p[i])
                i=i+1 
                
                
            #else:         
            #    newparams.update_stellar_mass(p[i])      
        return newparams                
                
    def update_with_mcmc_errors(self,p):
 
        '''Substitute normal errors with mcmc errors, where + and - errors can be different''' 
 
        #if (not (self.fit_performed) or self.never_saved): # just in case the user calls this function in a wrong moment
        #    print(self.fit_performed, self.never_saved)
       #     return
        #else:
       # print(p)
        i=0
        for idx in self.f_for_mcmc:
            if (idx<self.filelist.ndset):
                self.param_errors.update_offset_error(idx,p[i])
                i=i+1             
            elif (idx<2*self.filelist.ndset):
                self.param_errors.update_jitter_error(idx-self.filelist.ndset,p[i])
                i=i+1    
            elif (idx<2*self.filelist.ndset+7*self.npl):
                nr=idx-2*self.filelist.ndset
                x=int(nr/7)
                if (np.mod(nr,7)==0):
                    self.param_errors.update_Kerror(x,p[i])
                    i=i+1  
                elif (np.mod(nr,7)==1):
                    self.param_errors.update_Perror(x,p[i])
                    i=i+1                           
                elif (np.mod(nr,7)==2):
                    self.param_errors.update_eerror(x,p[i])
                    i=i+1       
                elif (np.mod(nr,7)==3):
                    self.param_errors.update_werror(x,p[i])
                    i=i+1       
                elif (np.mod(nr,7)==4):
                    self.param_errors.update_M0error(x,p[i])
                    i=i+1       
                elif (np.mod(nr,7)==5):
                    self.param_errors.update_inclination_error(x,p[i])
                    i=i+1  
                elif (np.mod(nr,7)==6):
                    self.param_errors.update_lineofnodes_error(x,p[i])
                    i=i+1                          
            elif (idx<2*self.filelist.ndset+7*self.npl+1):
                self.param_errors.update_linear_trend_error(p[i]) 
                i=i+1       
            elif (idx<2*self.filelist.ndset+7*self.npl+2+self.params.GP_params.npar):
                self.param_errors.update_GP_param_errors(idx-2*self.filelist.ndset-7*self.npl-2,p[i])
                #print(idx-2*self.filelist.ndset-7*self.npl-2,p[i])  #TBD!!!
                
                i=i+1
            else:    
                self.param_errors.update_stellar_mass_error(p[i])            
        return                           

  

     
        
    def run_stability_one_sample(self,sample_num,fileinput=False, filename='samples_kep', fileinputgetinit=False, filenamegetinit='geninit_j_input', warnings=None, timeout_sec=1000.0,timemax=3000.0, timestep=10):

        '''Runs stability analysis on one sample. Optionally samples can be read from a file.'''

        print_at_end=False

        if (warnings==None): # this maens we are running this function stand alone and no warnings object was passed, so we need to create a new warnings object and print warnings at the end of this function
            warnings=Warning_log([],'Stability analisys for sample %d'%sample_num)
            print_at_end=True

        if(fileinput):
            samples=read_file_as_array_of_arrays_mcmc(filename)
        elif(self.sampler_saved):
            samples=self.sampler.samples
        else:
            raise Exception ('Please run mcmc and save sampler or provide a valid samples file!')     
            
        if (sample_num>len(samples)):  
            warnings.update_warning_list('There are only %d samples but %d was provided as sample number, cannot analyse, will be skipped!'%(len(samples),sample_num))
        else:
            newparams=self.generate_newparams_for_mcmc(samples[sample_num]) # based on the sample create a parameters object which we will work on
        
            # we need to update masses and semimajor axes according to the parameters
            oldparams=self.params
            self.overwrite_params(newparams)
            self.mass_semimajor()
            self.overwrite_params(oldparams)

            custom_param_file_for_stability(timemax,timestep)

            # runnning fortran codes
            result, flag = run_command_with_timeout(newparams.getinit_input(self.masses,self.semimajor, fileinput=fileinputgetinit, filename=filenamegetinit), timeout_sec)         
            result, flag = run_command_with_timeout('./swift_symba5_j << EOF \nparam.in \npl.in \n1e-40 \nEOF', timeout_sec)                  
 
            for k in range(self.npl):
                result, flag = run_command_with_timeout('./follow_symba2 << EOF \nparam.in \npl.in \n%s \nEOF'%(k+2),timeout_sec)
                result, flag = run_command_with_timeout('mv follow_symba.out pl_%s.out'%(k+1),timeout_sec) 

        return
         
    def run_stability_last_fit_params_new(self, timemax=3000.0, timestep=10, timeout_sec=1000.0):

        custom_param_file_for_stability(timemax,timestep)

        result, flag = run_command_with_timeout(self.params.getinit_input(self.masses,self.semimajor, fileinput=fileinputgetinit, filename=filenamegetinit), timeout_sec)         
        result, flag = run_command_with_timeout('./swift_symba5_j << EOF \nparam.in \npl.in \n1e-40 \nEOF', timeout_sec)                  
 
        
        for k in range(self.npl):
            result, flag = run_command_with_timeout('./follow_symba2 << EOF \nparam.in \npl.in \n%s \nEOF'%(k+2),timeout_sec)

            result, flag = run_command_with_timeout('mv follow_symba.out pl_%s.out'%(k+1),timeout_sec) 

            self.evol_T[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [0]) /  365.2425
            self.evol_a[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [2])
            self.evol_e[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [3])
            self.evol_p[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [6])      
            self.evol_M[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [7])

        return        
        
        
    def run_stability_last_fit_params(self, timemax=3000.0, timestep=10, timeout_sec=1000.0, stab_save_dir = './', integrator='symba'):

#if not os.path.exists(directory):
#    os.makedirs(directory)

        if integrator=='symba':
            os.chdir('./stability/symba/')
        elif integrator=='mvs':
            os.chdir('./stability/mvs/')
        elif integrator=='mvs_gr':
            os.chdir('./stability/mvs_gr/')
    
        print("running stability with: %s"%integrator) 
        ##### crate the param.in file (change only the "t_max" and the "dt" for now) ######
        param_file = open('param.in', 'wb') 
        
        max_time = float(timemax)*365.2425 # make it is days
 
        param_file.write(b"""0.0d0 %s %s
%s %s
        
F T T T T F
0.0001 50.0 50.0 -1. T
bin.dat
unknown
"""%(bytes(str(max_time).encode()), 
     bytes(str(timestep).encode()),
     bytes(str(max_time/1e4).encode()), 
     bytes(str(max_time/1e3).encode())  ))
 
        param_file.close()
        
        #os.system("cp param.in test_param.in__")
        
        
        getin_file = open('geninit_j.in', 'wb') 
        getin_file.write(b"""1 
%s
%s
1.d0
pl.in
    """%(bytes(str(self.params.stellar_mass).encode()), bytes(str(self.npl).encode() ) ))
        
        
  
        for j in range(self.npl):
            getin_file.write(b'%s \n'%bytes(str(self.fit_results.mass[j]/1047.70266835).encode())) 
            getin_file.write(b'%s %s %s %s %s %s \n'%(bytes(str(self.fit_results.a[j]).encode()),
                                                     bytes(str(self.params.planet_params[7*j + 2]).encode()),
                                                     bytes(str(self.params.planet_params[7*j + 5]).encode()),
                                                     bytes(str(self.params.planet_params[7*j + 3]).encode()),
                                                     bytes(str(self.params.planet_params[7*j + 6]).encode()),
                                                     bytes(str(self.params.planet_params[7*j + 4]).encode() )) ) 
              
        getin_file.close()

        # runnning fortran codes
        result, flag = run_command_with_timeout('./geninit_j3_in_days < geninit_j.in', timeout_sec)         

        if integrator=='symba':
            result, flag = run_command_with_timeout('./swift_symba5_j << EOF \nparam.in \npl.in \n1e-40 \nEOF', timeout_sec)                  
        elif integrator=='mvs':
            result, flag = run_command_with_timeout('./swift_mvs_j << EOF \nparam.in \npl.in \nEOF', timeout_sec)                          
        elif integrator=='mvs_gr':
            result, flag = run_command_with_timeout('./swift_mvs_j_GR << EOF \nparam.in \npl.in \nEOF', timeout_sec)          
                 
        
        for k in range(self.npl):
        
            if integrator=='symba':
                result, flag = run_command_with_timeout('./follow_symba2 << EOF \nparam.in \npl.in \n%s \nEOF'%(k+2),timeout_sec)
                result, flag = run_command_with_timeout('mv follow_symba.out pl_%s.out'%(k+1),timeout_sec) 
            elif integrator=='mvs' or integrator=='mvs_gr': 
                result, flag = run_command_with_timeout('./follow2 << EOF \nparam.in \npl.in \n-%s \nEOF'%(k+2),timeout_sec)
                result, flag = run_command_with_timeout('mv follow2.out pl_%s.out'%(k+1),timeout_sec)                 

            self.evol_T[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [0]) /  365.2425
            self.evol_a[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [2])
            self.evol_e[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [3])
            self.evol_p[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [6])      
            self.evol_M[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [7])
        
        try:
            os.system('rm *.out *.dat *.in') 
        except OSError:
            pass
        
        os.chdir('../../')
        
        print("stability with: %s done!"%integrator) 
        
        return  
        
              
