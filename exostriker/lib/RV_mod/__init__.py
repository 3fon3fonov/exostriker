#!/usr/bin/python3

##
######## this is the file where all experiments are happening!
## the code below is still work in progress.
from __future__ import print_function


import sys, os
#sys.path.insert(0, '../lib')
#sys.path.append('./lib/RV_mod/')
#import prior_functions as pr
#import emcee
import emcee_ES as emcee
import numpy as np
#import matplotlib.pyplot as plt
#plt.switch_backend('SVG')

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import time
#import multiprocessing
from pathos import multiprocessing

import celerite
from celerite import terms
import dill
dill.settings['fmode']

#import gc
import dynesty_2_0 as dynesty
#import dynesty
import dynesty_patch
dynesty.results =  dynesty_patch 


import scipy.optimize as op
from scipy import stats

####### fixes https://github.com/3fon3fonov/exostriker/issues/80   ??? ####
import wrapt_ES
import tqdm.std

methods = ['__del__', 'close']
for method_name in methods:
    @wrapt_ES.patch_function_wrapper(tqdm.std.tqdm, method_name)
    def new_del(wrapped, instance, args, kwargs):
        try:
            return wrapped(*args, **kwargs)
        except AttributeError:
            pass
##########################################################################


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
    


try:
    import ttvfast as ttvfast
    ttvfast_not_found = False
except (ImportError, KeyError) as e:
    ttvfast_not_found = True
    pass

try:
    import astropy.io.fits as pyfits
    
    pyfits_not_found = False
except (ImportError, KeyError) as e:
    pyfits_not_found = True
    pass


from .parameters import *
from .CustomSampler import CustomSampler
from .Warning_log import Warning_log
from .functions import *

from .errors import Error, InputError, FittingError
from .rv_files import rvfile, rvfile_list

from .astrometry_ES import * 

TAU= 2.0*np.pi

from . import GP_kernels
from .rvmod import *

import ExTRA as ex
 

def initiate_RV_gps(obj,  kernel_id=-1):
    """Short summary.

    Parameters
    ----------
    obj : type
        Description of parameter `obj`.
    kernel_id : type
        Description of parameter `kernel_id`.

    Returns
    -------
    type
        Description of returned object.

    """

    # Prepare objects for Gaussian Processes

    if len(obj.GP_rot_params) != 0:
        obj.params.update_GP_params(obj.GP_rot_params, kernel_id=kernel_id)
        obj.use.update_use_GP_params(obj.GP_rot_use)

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

    elif obj.gp_kernel == 'Matern32':
        kernel = terms.Matern32Term(log_sigma=np.log(obj.GP_mat_params[0]), 
                                    log_rho=np.log(obj.GP_mat_params[1]), 
                                    eps=obj.GP_mat_params[2])    
        
    elif obj.gp_kernel == 'RealTerm':
        kernel = terms.RealTerm(log_a=np.log(obj.GP_drw_params[0]), 
                                log_c=np.log(obj.GP_drw_params[1]))   

    elif obj.gp_kernel == 'dSHOKernel':
        kernel =GP_kernels.double_SHOTerm(
                sigma_dSHO=obj.GP_double_sho_params[0],
                period_dSHO=obj.GP_double_sho_params[1],
                Q0_dSHO=obj.GP_double_sho_params[2],
                dQ_dSHO=obj.GP_double_sho_params[3], 
                f_dSHO=obj.GP_double_sho_params[4])

 
    gps = celerite.GP(kernel, mean=0.0)
    #gps.compute(obj.filelist.time, obj.filelist.rv_err)

    obj.gps = gps
    return


def init_dSHOKernel(params):
    kernel =GP_kernels.double_SHOTerm(
            sigma_dSHO=params[0],
            period_dSHO=params[1],
            Q0_dSHO=params[2],
            dQ_dSHO=params[3], 
            f_dSHO=params[4])


    gps = celerite.GP(kernel, mean=0.0)
    return gps

def get_RV_gps_model(obj,  kernel_id=-1, get_lnl=False):
 
    
    initiate_RV_gps(obj,  kernel_id=-1)
    #gp_model_data  = []

    ############ DATA ####################
    #for i in range(obj.filelist.ndset):
        #gp.set_parameter_vector(

    y = obj.fit_results.rv_model.o_c
    x = obj.fit_results.rv_model.jd

    errors_with_jitt = np.array([np.sqrt(obj.fit_results.rv_err[i]**2 + obj.jitt[ii]**2)  for i,ii in enumerate(obj.fit_results.idset)])

    obj.gps.compute(x,yerr=errors_with_jitt)    
   #gps.compute(obj.filelist.time, obj.filelist.rv_err)

    
    mu, var = obj.gps.predict(y, x, return_var=True)
    std = np.sqrt(var)

    obj.gp_model_data = [mu,var,std]

    ############ MODEL ####################

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
    
 
    if get_lnl == True:
        
        gp_rv_loglik = obj.gps.log_likelihood(obj.fit_results.o_c)
        
        gp_pred = obj.gps.predict(obj.fit_results.o_c,obj.fit_results.jd, return_cov=False)
        o_c_kep = dill.copy(obj.fit_results.o_c) - gp_pred
        
        gp_rv_chi = 0

        for i in range(max(obj.fit_results.idset)+1):
            sig2i_gp = 1.0 / (obj.fit_results.rv_err[obj.fit_results.idset==i]**2 + obj.jitt[i]**2 )
       #     gp_rv_loglik += -0.5*(np.sum((o_c_kep[obj.fit_results.idset==i])**2 * sig2i_gp - np.log(sig2i_gp / 2./ np.pi)))
            gp_rv_chi    += np.sum((o_c_kep[obj.fit_results.idset==i])**2 * sig2i_gp)


            
        obj.fit_results.chi2 = gp_rv_chi
        obj.fit_results.reduced_chi2 = gp_rv_chi /(len(o_c_kep) - len(obj.par_for_mcmc))

        obj.fit_results.rms = np.sqrt(np.average(o_c_kep**2))
        obj.fit_results.wrms =  np.sqrt(np.average(o_c_kep**2, weights=1/obj.fit_results.rv_err))

        if obj.gp_kernel == 'RotKernel':
            N_gp_pars_used = len([i for i in range(4) if obj.GP_rot_use[i] == True])
        elif obj.gp_kernel == 'SHOKernel':
            N_gp_pars_used = len([i for i in range(3) if obj.GP_sho_use[i] == True])
        elif obj.gp_kernel == 'Matern32':
            N_gp_pars_used = len([i for i in range(2) if obj.GP_mat_use[i] == True])    
        elif obj.gp_kernel == 'RealTerm':
            N_gp_pars_used = len([i for i in range(2) if obj.GP_drw_use[i] == True])            
        elif obj.gp_kernel == 'dSHOKernel':
            N_gp_pars_used = len([i for i in range(5) if obj.GP_double_sho_use[i] == True])

        obj.fit_results.stat.dof = obj.fit_results.stat.dof - N_gp_pars_used
        obj.loglik = gp_rv_loglik

    return



def initiate_tansit_gps(obj,  kernel_id=-1):

    # Prepare objects for Gaussian Processes

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

    elif obj.tra_gp_kernel == 'Matern32':
        tra_kernel = terms.Matern32Term(log_sigma=np.log(obj.tra_GP_mat_params[0]), 
                                    log_rho=np.log(obj.tra_GP_mat_params[1]), 
                                    eps=obj.tra_GP_mat_params[2])     

    elif obj.tra_gp_kernel == 'RealTerm':
        tra_kernel = terms.RealTerm(log_a=np.log(obj.tra_GP_drw_params[0]), 
                                    log_c=np.log(obj.tra_GP_drw_params[1]))   

    elif obj.tra_gp_kernel == 'dSHOKernel':
        tra_kernel =GP_kernels.double_SHOTerm(
                sigma_dSHO=obj.tra_GP_double_sho_params[0],
                period_dSHO=obj.tra_GP_double_sho_params[1],
                Q0_dSHO=obj.GP_double_sho_params[2],
                dQ_dSHO=obj.tra_GP_double_sho_params[3], 
                f_dSHO=obj.tra_GP_double_sho_params[4])

    tra_gps = celerite.GP(tra_kernel, mean=0.0)


    if len([obj.tra_data_sets[j][9] for j in range(60) if len(obj.tra_data_sets[j]) != 0 and obj.tra_data_sets[j][9] ==True]) ==0:
        print("No transit data ready for GP modeling!!! Reverting to 'GP==False'")
        obj.tra_gps = tra_gps
        obj.tra_doGP = False
        #obj.hkl[3] = False
        return        

    else:
        y = np.concatenate([obj.tra_data_sets[j][2] for j in range(60) if len(obj.tra_data_sets[j]) != 0 and obj.tra_data_sets[j][9] ==True])
        x = np.concatenate([obj.tra_data_sets[j][0] for j in range(60) if len(obj.tra_data_sets[j]) != 0 and obj.tra_data_sets[j][9] ==True])
        
        
        tra_gps.compute(x, y)
    
        obj.tra_gps = tra_gps
        return




def get_transit_gps_model(obj, x_model = [], y_model = [],  kernel_id=-1):

    get_transit_ts(obj)
    initiate_tansit_gps(obj,  kernel_id=-1)
    #gp_model_data  = []

    ############ DATA ####################

    y = np.concatenate([obj.tra_data_sets[j][4] for j in range(60) if len(obj.tra_data_sets[j]) != 0 and obj.tra_data_sets[j][9] ==True])
    x = np.concatenate([obj.tra_data_sets[j][0] for j in range(60) if len(obj.tra_data_sets[j]) != 0 and obj.tra_data_sets[j][9] ==True])

   # y_no_gp = np.concatenate([obj.tra_data_sets[j][4] for j in range(10) if len(obj.tra_data_sets[j]) != 0 and obj.tra_data_sets[j][9] !=True])


    if len(x_model) == 0 or len(y_model) == 0:
        x_model = dill.copy(x)
        y_model = dill.copy(y)

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
       # obj.tra_gps.compute(x_model, y_model)

        mu = obj.tra_gps.predict(y, x_model, return_cov=False)
       # std = np.sqrt(var)
       

        obj.tra_gp_model_curve= [mu,np.zeros(len(mu)),np.zeros(len(mu))]


    return 



def get_transit_ts(obj,  kernel_id=-1):

    tr_files = []


    for j in range(60):

        if len(obj.tra_data_sets[j]) == 0:
            continue
        else:
    
            t = np.array(obj.tra_data_sets[j][0])
            flux = np.array(obj.tra_data_sets[j][1] + obj.tra_off[j])
            flux_err = np.sqrt(obj.tra_data_sets[j][2]**2 + obj.tra_jitt[j]**2)



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
            if obj.gp_kernel == 'Matern32':
                rv_gp_npar = 3   
            if obj.gp_kernel == 'RealTerm':
                rv_gp_npar = 2                   
            if obj.gp_kernel == 'dSHOKernel':
                rv_gp_npar = 6                   
            #fit.gps = []
        else:
            rv_gp_npar = 0

        obj.tr_params.limb_dark = str(obj.ld_m[j])      #limb darkening model
        obj.tr_params.u = obj.ld_u[j]


        for i in range(9):
            if not bool(obj.use_planet[i]):
                continue

            if obj.hkl == True:
                obj.tr_params.ecc = np.sqrt(par[obj.ndset*2 +7*i+2]**2 + par[obj.ndset*2 +7*i+3]**2)
                obj.tr_params.w  = np.degrees(np.arctan2(par[obj.ndset*2 +7*i+2],par[obj.ndset*2 +7*i+3]))%360
            else:
                obj.tr_params.ecc = par[obj.ndset*2 +7*i+2] #0.0
                obj.tr_params.w   = par[obj.ndset*2 +7*i+3] #90

            obj.tr_params.per = par[obj.ndset*2 +7*i+1] #1.0    #orbital period
            obj.tr_params.inc = par[obj.ndset*2 +7*i+5]#90. #orbital inclination (in degrees)

            obj.tr_params.t0  = par[obj.ndset*2  +7*obj.npl +2+rv_gp_npar + 3*i]
            obj.tr_params.a   = par[obj.ndset*2  +7*obj.npl +2+rv_gp_npar + 3*i+1] #15  #semi-major axis (in units of stellar radii)
            obj.tr_params.rp  = par[obj.ndset*2  +7*obj.npl +2+rv_gp_npar + 3*i+2] #0.15   #planet radius (in units of stellar radii)

            m[i] = batman.TransitModel(obj.tr_params, t)    #initializes model

            flux_model = flux_model * m[i].light_curve(obj.tr_params)

        obj.tra_data_sets[j][4] = flux - flux_model

    return



def ttvs_mod(par,vel_files,npl, stellar_mass, times, planet_N, hkl, fit_results=False):

    planets = []
    for i in range(npl):
        
        if hkl == True:
            ecc_ = np.sqrt(par[len(vel_files)*2 +7*i+2]**2 + par[len(vel_files)*2 +7*i+3]**2)
            om_  = np.degrees(np.arctan2(par[len(vel_files)*2 +7*i+2],par[len(vel_files)*2 +7*i+3]))%360
            Ma_  = (par[len(vel_files)*2 +7*i+4] - om_)%360.0
            
            #om_  = (om_+180.0)%360.0
            #Ma_  = (Ma_-180.0)%360.0 
                       
        else:
            ecc_, om_, Ma_ = par[len(vel_files)*2 +7*i+2], par[len(vel_files)*2 +7*i+3], par[len(vel_files)*2 +7*i+4]        
 
            #om_  = (om_+180.0)%360.0
            #Ma_  = (Ma_-180.0)%360.0  
       
        if fit_results == False:
            pl_mass,ap = mass_a_from_Kepler_fit([par[len(vel_files)*2 + 7*i],
                                            par[len(vel_files)*2 +7*i+1],
                                            ecc_,
                                            om_,
                                            Ma_],1,stellar_mass) ##################TB FIXED! these are not dynamical masses!
        else:
            pl_mass = float(fit_results.mass[i])
        pl_params = [pl_mass/1047.5654817267318,
                                            par[len(vel_files)*2 +7*i+1],
                                            ecc_,
                                            par[len(vel_files)*2 +7*i+5],
                                            par[len(vel_files)*2 +7*i+6] ,
                                            #par[len(vel_files)*2 +7*i+3]%360.0,
                                            om_,
                                            #par[len(vel_files)*2 +7*i+4]%360.0]
                                            Ma_]

        planet = ttvfast.models.Planet(*pl_params)
        planets.append(planet)

    results = ttvfast.ttvfast(planets, stellar_mass,times[0],times[1],times[2],input_flag=0)
    result_rows = list(zip(*results['positions']))

    n1   = [item[0] for item in result_rows]
 
    n2              = np.array([item[1]+1 for i, item in enumerate(result_rows) if n1[i] == planet_N])
    transits_calc   = np.array([item[2]   for i, item in enumerate(result_rows) if n1[i] == planet_N])    

    rsky   = np.array([item[3]   for i, item in enumerate(result_rows) if n1[i] == planet_N])[np.where(transits_calc > 0)]
    vsky   = np.array([item[4]   for i, item in enumerate(result_rows) if n1[i] == planet_N])[np.where(transits_calc > 0)]

    #print(rsky)
    #print(vsky)

    n2               = n2[np.where(transits_calc > 0)]
    transits_calc    = transits_calc[np.where(transits_calc > 0)]
    rsky             = rsky[np.where(transits_calc > 0)]
    vsky             = vsky[np.where(transits_calc > 0)]

    calc_model = [n2,transits_calc,rsky,vsky]
 
 
    return calc_model


def ast_loglik(par,vel_files, ast_files,npl,stellar_mass, times, hkl, fit_results = False , return_model = False):

 
    loglik_ast = 0
    calc_data   = {kx: [] for kx in range(10)}
    calc_model  = {kx: [] for kx in range(10)}
    for i in range(npl):
        
        if hkl == True:
            ecc_ = np.sqrt(par[len(vel_files)*2 +7*i+2]**2 + par[len(vel_files)*2 +7*i+3]**2)
            om_  = np.degrees(np.arctan2(par[len(vel_files)*2 +7*i+2],par[len(vel_files)*2 +7*i+3]))%360
            Ma_  = (par[len(vel_files)*2 +7*i+4] - om_)%360.0
            # reject e > 0 
            if ecc_ >= 1.0: 
                return -np.inf  
        else:
            ecc_, om_, Ma_ = par[len(vel_files)*2 +7*i+2], par[len(vel_files)*2 +7*i+3], par[len(vel_files)*2 +7*i+4]         

        t0 = transit_tperi(par[len(vel_files)*2 +7*i+1],ecc_, om_, Ma_ ,times[0])[0] #%par[len(vel_files)*2 +7*i+1] 

        om_ = (om_ - 180.0)%360.0
        Ma_ = (Ma_ + 180.0)%360.0
    
        #t0    = times[0]  - ((np.radians(Ma_)/2.0*np.pi)*par[len(vel_files)*2 + 7*i])
        #print(t0,t0_o, times[0])
        if fit_results == False:
            pl_mass,pl_a = mass_a_from_Kepler_fit([par[len(vel_files)*2 + 7*i],
                                            par[len(vel_files)*2 +7*i+1],
                                            ecc_,
                                            om_,
                                            Ma_],1,stellar_mass) ##################TB FIXED! these are not dynamical masses!
        else:
            pl_mass = float(fit_results.mass[i])
            pl_a    = float(fit_results.a[i])


        for x in range(10):
            if len(ast_files[x]) == 0 or ast_files[x][6] == False or ast_files[x][5] != i+1:
                continue
            else:

                results = final_ast(ast_files[x][1],ast_files[x][2],0,
                                              ast_files[x][3],ast_files[x][4],0, 
                                                par[len(vel_files)*2 +7*i+1],
                                                ecc_,
                                                np.radians(om_),
                                                np.radians(par[len(vel_files)*2 +7*i+5]),
                                                np.radians(par[len(vel_files)*2 +7*i+6]),
                                                t0,pl_a,
                                                ast_files[x][0])
                loglik_ast = loglik_ast + results[0]
 

        if return_model == True:
            #times2 = np.linspace(min(ast_files[0][0]),max(ast_files[0][0]),1000)
            times2 = np.linspace(times[0],times[2],1000)

            ast_x_model, ast_y_model = ast_coords(par[len(vel_files)*2 +7*i+1], ecc_, om_, np.radians(par[len(vel_files)*2 +7*i+5]),
                                                    np.radians(par[len(vel_files)*2 +7*i+6]),
                                                    t0,pl_a,times2) 
            ast_x_, ast_y_ = ast_coords(par[len(vel_files)*2 +7*i+1], ecc_, om_, np.radians(par[len(vel_files)*2 +7*i+5]),
                                                    np.radians(par[len(vel_files)*2 +7*i+6]),
                                                    t0,pl_a,ast_files[0][0]) 

            calc_data[i] = np.array([ast_x_, ast_y_])         
            calc_model[i] = np.array([ast_x_model, ast_y_model])

 

    if return_model == True:
        return [loglik_ast, calc_data,calc_model]
    else:
        return loglik_ast


def ast_loglik_hipp(par,vel_files, ast_files_hipp,npl,stellar_mass, times, hkl, fit_results = False , return_model = False):

 
    loglik_ast_hipp = 0
    calc_data   = {kx: [] for kx in range(10)}
    calc_model  = {kx: [] for kx in range(10)}
    for i in range(npl):
        
        if hkl == True:
            ecc_ = np.sqrt(par[len(vel_files)*2 +7*i+2]**2 + par[len(vel_files)*2 +7*i+3]**2)
            om_  = np.degrees(np.arctan2(par[len(vel_files)*2 +7*i+2],par[len(vel_files)*2 +7*i+3]))%360
            Ma_  = (par[len(vel_files)*2 +7*i+4] - om_)%360.0
            # reject e > 0 
            if ecc_ >= 1.0: 
                return -np.inf  
        else:
            ecc_, om_, Ma_ = par[len(vel_files)*2 +7*i+2], par[len(vel_files)*2 +7*i+3], par[len(vel_files)*2 +7*i+4]         

        T0 = transit_tperi(par[len(vel_files)*2 +7*i+1],ecc_, om_, Ma_ ,times[0])[0] #%par[len(vel_files)*2 +7*i+1] 

        om_ = (om_ - 180.0)%360.0
        Ma_ = (Ma_ + 180.0)%360.0
        
        incl = par[len(vel_files)*2 +7*i+5]
        Om_ = (par[len(vel_files)*2 +7*i+6] + 180.0)%360.0
 
        for x in range(10):
            if len(ast_files_hipp[x]) == 0 or ast_files_hipp[x][8] == False or ast_files_hipp[x][7] != i+1:
                continue
            else:

                A_6=ast_files_hipp[x][3]*(ast_files_hipp[x][1])
                A_7=ast_files_hipp[x][4]*(ast_files_hipp[x][1])
 
                hip_ad=np.array([ast_files_hipp[x][3],
                                 ast_files_hipp[x][4],
                                 ast_files_hipp[x][2],
                                 A_6,A_7,
                                 ast_files_hipp[x][5],
                                 ast_files_hipp[x][6]])                
  
                t_HIP=ex.hip_JD(hip_ad)
 
                #J2016=2457389.0
                #J1991=2448349.0625
                
                #gaia=np.array([par[-6],par[-5],par[-4],par[-3],par[-2]])
                #hip_standard=np.zeros(5)
                #correction=np.array([0,0,0,0,0]) 
                
                hip_stand = [
                ast_files_hipp[x][9]["RAdeg"], 
                ast_files_hipp[x][9]["DEdeg"], 
                ast_files_hipp[x][9]["Plx"], 
                ast_files_hipp[x][9]["pm_RA"], 
                ast_files_hipp[x][9]["pm_DE"], 
                ]
                
                parallax=hip_stand[2]+par[-4] #corrected parallax
                
                correction=np.array([par[-6],par[-5],par[-4],par[-3],par[-2]])

                a = a_from_K_in_mas(
                par[len(vel_files)*2 +7*i+1],
                par[len(vel_files)*2 +7*i+0],
                ecc_,
                incl,
                parallax)

                results = -1*ex.L_hip(hip_ad,hip_stand,hip_stand, #Data
                         correction,
                         par[len(vel_files)*2 +7*i+1], 
                         ecc_, 
                         np.radians(om_), 
                         np.radians(incl),
                         np.radians(Om_),                         
                         T0,a,Sepoch=ex.J1991(),s_hip=0)

                loglik_ast_hipp = loglik_ast_hipp + results
                
                
                
        if return_model == True:

            #### TB fixed!
            #corrected_stand = hip_stand
            corrected_stand=ex.stand_correct(hip_stand,correction)

            #times3 = np.linspace(min(t_HIP),min(t_HIP)+par[len(vel_files)*2 +7*i+1]+1,1000)
            times3 = np.linspace(min(t_HIP),max(t_HIP)+1,1000)
            orb_params = [par[len(vel_files)*2 +7*i+1], 
                         ecc_, 
                         np.radians(om_), 
                         np.radians(incl),
                         np.radians(Om_),                         
                         T0,a]
 

            ast_x_model,ast_y_model=ex.orbit(*orb_params,times3)
            
            hip_data = ex.hip_2d(hip_ad) #Rotates the HIP measurements into the RA,Dec frame.
 
 
            hip_data_res_ =ex.hip_residuals(hip_ad,hip_stand,corrected_stand,orb_params, Sepoch=ex.J1991())
 
            hip_data_res  =ex.res_to_orbit(hip_data_res_,hip_ad,orb_params)
 
            calc_data[i] = np.array([hip_data, hip_data_res ,t_HIP],dtype=object)         
            calc_model[i] = np.array([ast_x_model,ast_y_model,times3])
 
    
    if return_model == True:
        return [loglik_ast_hipp, calc_data, calc_model]
    else:
        return loglik_ast_hipp


def ttvs_loglik(par,vel_files,ttv_files,npl,stellar_mass,times, hkl, fit_results = False , return_model = False):


    planets = []
    for i in range(npl):
        
        if hkl == True:
            ecc_ = np.sqrt(par[len(vel_files)*2 +7*i+2]**2 + par[len(vel_files)*2 +7*i+3]**2)
            om_  = np.degrees(np.arctan2(par[len(vel_files)*2 +7*i+2],par[len(vel_files)*2 +7*i+3]))%360
            Ma_  = (par[len(vel_files)*2 +7*i+4] - om_)%360.0
            
           # om_ = (om_+180.0)%360.0
           # Ma_ = (Ma_-180.0)%360.0            
            
            # reject e > 0 
            if ecc_ >= 1.0: 
                return -np.inf  
        else:
            ecc_, om_, Ma_ = par[len(vel_files)*2 +7*i+2], par[len(vel_files)*2 +7*i+3], par[len(vel_files)*2 +7*i+4]         

            #om_ = (om_+180.0)%360.0
            #Ma_ = (Ma_-180.0)%360.0           
        
        if fit_results == False:
            pl_mass,ap = mass_a_from_Kepler_fit([par[len(vel_files)*2 + 7*i],
                                            par[len(vel_files)*2 +7*i+1],
                                            ecc_,
                                            om_,
                                            Ma_],1,stellar_mass) ##################TB FIXED! these are not dynamical masses!
        else:
            pl_mass = float(fit_results.mass[i])
            #print(pl_mass)
        pl_params = [pl_mass/1047.5654817267318,
                                            par[len(vel_files)*2 +7*i+1],
                                            ecc_,
                                            par[len(vel_files)*2 +7*i+5],
                                            par[len(vel_files)*2 +7*i+6],
                                            om_,
                                            Ma_]
        planet = ttvfast.models.Planet(*pl_params)
        planets.append(planet)

    #rv_times  = list(np.linspace(times[0],times[2], 1000)) rv_times=rv_times,
 
    results = ttvfast.ttvfast(planets, stellar_mass, times[0],times[1],times[2], input_flag=0)
    result_rows = list(zip(*results['positions']))
    #result_rows_rv = list(results['rv'])
 


    n1   = [item[0] for item in result_rows]

    loglik_ttv = 0


    calc_data   = {kx: [] for kx in range(10)}
    calc_model  = {kx: [] for kx in range(10)}

    for x in range(10):
        if len(ttv_files[x]) == 0 or ttv_files[x][4] == False:
            continue
        else:
        
            calc_n     = []
            calk_tran  = []
        
            #n2   = np.array([item[1]+1 for i, item in enumerate(result_rows) if n1[i] == 0])
            #transits_calc   = np.array([item[2] for i, item in enumerate(result_rows) if n1[i] == 0])
            n2              = np.array([item[1]+1 for i, item in enumerate(result_rows) if n1[i] == int(ttv_files[x][3])-1])# if ttv_files[0][4] == True])
            transits_calc   = np.array([item[2]   for i, item in enumerate(result_rows) if n1[i] == int(ttv_files[x][3])-1])# if ttv_files[0][4] == True])        
            rsky            = np.array([item[3]   for i, item in enumerate(result_rows) if n1[i] == int(ttv_files[x][3])-1])
            vsky            = np.array([item[4]   for i, item in enumerate(result_rows) if n1[i] == int(ttv_files[x][3])-1])


 
            # A bug fix???
            if len(n2) == 0 or max(n2) < max(ttv_files[x][0]):
                if return_model == True:
                    print("Number of calc. transits:",len(n2))
                    print("Number of obs.  transits:",max(ttv_files[x][0]) )
                    return None
                else:
                    return -np.inf
                
            for i in range(len(ttv_files[x][1])):
        
                calc_n.append(n2[int(ttv_files[x][0][i] -1)])
                calk_tran.append(transits_calc[int(ttv_files[x][0][i])-1])
        
                sig2i_ttv = 1.0 / (ttv_files[x][2][i])**2
                loglik_ttv += -0.5*(np.sum(((ttv_files[x][1][i] - calk_tran[i])**2 * sig2i_ttv - np.log(sig2i_ttv / 2./ np.pi))))
        
            n2  = n2[np.where(transits_calc > 0)]
            transits_calc  = transits_calc[np.where(transits_calc > 0)]
            rsky             = rsky[np.where(transits_calc > 0)]
            vsky             = vsky[np.where(transits_calc > 0)]      
                    
        calc_data[x] = [calc_n,calk_tran,ttv_files[x][2]]
        calc_model[x] = [n2,transits_calc,rsky,vsky]

    if return_model == True:
        return [loglik_ttv, [calc_n,calk_tran],[n2,transits_calc,rsky,vsky],calc_data,calc_model]
    else:
        return loglik_ttv


def transit_loglik(program, tr_files,vel_files,tr_params,tr_model,par,rv_gp_npar,tra_gp_npar,npl,hkl, rtg, tra_gps, stmass, ttv_times, epoch, get_TTVs, opt, return_model = False, tra_model_fact = 10 ,fit_results=False):

    tr_loglik     = 0
    flux_model    = []
    flux          = []
    t             = []
    flux_err      = []
    sig2i         = []
    flux_o_c      = []
    flux_o_c_gp   = []
    tra_gp_model  = []
    use_gp_model  = []

    t_rich  = []
    flux_model_rich  = []

    N_transit_files = len([x for x in range(60) if len(tr_files[x]) != 0])

    if opt['rho'] == True:
        
        rho = par[len(vel_files)*2 +7*npl +2 +rv_gp_npar + 3*0+1]
        
    l = 0

    for j in range(60):

        if len(tr_files[j]) == 0:
            flux_model.append([])
            flux.append([])
            flux_err.append([])
            t.append([])
            sig2i.append([])
            flux_o_c.append([])
            flux_o_c_gp.append([])
            tra_gp_model.append([])
            use_gp_model.append([])
            continue

        t_        = tr_files[j][0]
        flux_     = tr_files[j][1] + par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + l] #* tr_files[j][8]

        sig2i_    = 1./(tr_files[j][2]**2 + par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files + l]**2)
        flux_err_ = np.sqrt(tr_files[j][2]**2 + par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files + l]**2)
        flux_model_ =np.ones(len(flux_))
       


        m  =  {z: [] for z in range(9)}
        tr_params.limb_dark = str(tr_model[0][j])

        
        if tr_model[0][j] == "linear":
            tr_params.u = [par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*4 + tra_gp_npar + tr_model[3][j] +  npl]]
#            k += 1
        elif tr_model[0][j] == "quadratic":
            tr_params.u = [par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*4 + tra_gp_npar + tr_model[3][j] + npl], 
                           par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*4 + tra_gp_npar + tr_model[3][j] + 1+ npl]]
#            k += 2
        elif tr_model[0][j] == "nonlinear":
            tr_params.u = [par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*4 + tra_gp_npar + tr_model[3][j] + npl], 
                           par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*4 + tra_gp_npar + tr_model[3][j] + 1+ npl],
                           par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*4 + tra_gp_npar + tr_model[3][j] + 2+ npl], 
                           par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*4 + tra_gp_npar + tr_model[3][j] + 3+ npl]]
#            k += 4
        else:
            tr_params.u = []


        for i in range(npl):


            if opt['rho'] == True:
                tr_params.a   = get_aR_from_rho(rho, par[len(vel_files)*2 +7*i+1])      
            else:      
                tr_params.a   = par[len(vel_files)*2 +7*npl +2 +rv_gp_npar + 3*i+1]

            if hkl == True:
                tr_params.ecc = np.sqrt(par[len(vel_files)*2 +7*i+2]**2 + par[len(vel_files)*2 +7*i+3]**2)
                tr_params.w  = np.degrees(np.arctan2(par[len(vel_files)*2 +7*i+2],par[len(vel_files)*2 +7*i+3]))%360
            else:
                tr_params.ecc = par[len(vel_files)*2 +7*i+2]
                tr_params.w   = par[len(vel_files)*2 +7*i+3]

            tr_params.per = par[len(vel_files)*2 +7*i+1]
            tr_params.inc = par[len(vel_files)*2 +7*i+5]

            tr_params.t0  = par[len(vel_files)*2 +7*npl +2 +rv_gp_npar + 3*i]
            tr_params.rp  = par[len(vel_files)*2 +7*npl +2 +rv_gp_npar + 3*i+2]


            if program.endswith("dyn+") or program.endswith("dyn"):

                t_os = ttvs_mod(par,vel_files,npl, stmass, [epoch,ttv_times[1],max(t_)], i, hkl, fit_results=fit_results)
                
                for tran_t0 in t_os[1]:
                    tr_params.t0  = float(tran_t0)
                    
                    
                    

                    m[i] = batman.TransitModel(tr_params, t_)
                    tr_ind = np.where(np.logical_and(t_ >= tran_t0-0.3, t_ <= tran_t0+0.3))
                    flux_model_[tr_ind] = m[i].light_curve(tr_params)[tr_ind]


            elif get_TTVs[0] == True:
                
                t_os = par[get_TTVs[1][i][0]:get_TTVs[1][i][1]]
                
                for tran_t0 in t_os:
                    tr_params.t0  = float(tran_t0)

                    m[i] = batman.TransitModel(tr_params, t_)
                    tr_ind = np.where(np.logical_and(t_ >= tran_t0-0.3, t_ <= tran_t0+0.3))
                    flux_model_[tr_ind] = m[i].light_curve(tr_params)[tr_ind]    
            else:
                m[i] = batman.TransitModel(tr_params, t_)
                flux_model_ = flux_model_ * m[i].light_curve(tr_params) 


        if tr_files[j][10] == True:
            flux_model_ = get_airmass_model(tr_files[j][3],flux_model_,0.0,
                                                    par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*2 + tra_gp_npar + l],
                                                    par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*3 + tra_gp_npar + l])           
        else:
            flux_model_ = get_quad_model(t_,flux_model_,0.0,
                                                    par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*2 + tra_gp_npar + l],
                                                    par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*3 + tra_gp_npar + l])



        flux_model_ = flux_model_*tr_files[j][8]  + (1.0 - tr_files[j][8]) 
        

        #flux_model_ = (flux_model_*tr_files[j][8]  + (1.0 - tr_files[j][8])) /  (1+ tr_files[j][8])

        #flux_model_ =  (flux_model_*tr_files[j][8]  + (1.0 - tr_files[j][8]))/  (1+ tr_files[j][8]*par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + l])
        #flux_model_ = flux_model_ + par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + l]

        l +=1

        flux_o_c_ = np.array(flux_) - np.array(flux_model_)

        ###### TBD, GP for each transit dataset #####
          
        #tra_gp_model.append(flux_model_)

        flux_model.append(flux_model_)
        flux.append(flux_)
        flux_err.append(flux_err_)
        t.append(t_)
        sig2i.append(sig2i_)
        flux_o_c.append(flux_o_c_)
        use_gp_model.append(tr_files[j][9])

        if rtg[3] == True and return_model == True:


            param_vect = []

            if opt["tra_GP_kernel"] == "dSHOKernel":
                for k in range(len(tra_gps.get_parameter_vector())-1):
                    param_vect.append(par[len(vel_files)*2  +7*npl  + rv_gp_npar  + 3*npl + N_transit_files*2 + 2 + k])

                tra_gps = init_dSHOKernel(param_vect)

            else:

                for k in range(len(tra_gps.get_parameter_vector())):
                    param_vect.append(np.log(par[len(vel_files)*2  +7*npl  + rv_gp_npar  + 3*npl + N_transit_files*2 + 2 + k]))
                tra_gps.set_parameter_vector(np.array(param_vect))
  

            tra_gps.compute(t_,yerr=flux_err_)
            flux_model_gp   = tra_gps.predict(flux_o_c_, t_ , return_cov=False)
            flux_o_c_gp_  = flux_o_c_  - flux_model_gp  
 
            flux_o_c_gp.append(flux_o_c_gp_)
            tra_gp_model.append(flux_model_gp)
        else:
            flux_o_c_gp.append(flux_o_c_)
            tra_gp_model.append(flux_model_)

    
    flux_model_all  = np.concatenate(flux_model)#[flux_model_ < 1]
    flux_all        = np.concatenate(flux)#[flux_model_ < 1]
    flux_err_all    = np.concatenate(flux_err)#[flux_model_ < 1]
    t_all           = np.concatenate(t)#[flux_model_ < 1]
    sig2i_all       = np.concatenate(sig2i)#[flux_model_ < 1]

    flux_o_c_all   = flux_all -flux_model_all
 
    if rtg[3] == False:

        tr_loglik = tr_loglik -0.5*(np.sum((flux_all  -flux_model_all)**2 * sig2i_all  - np.log(sig2i_all  / 2./ np.pi)))
        
        if return_model == True:
            tra_gp_model_all = flux_model_all
            flux_o_c_gp_all = flux_o_c_all
    else:
        

        param_vect = []

        if opt["tra_GP_kernel"] == "dSHOKernel":
            for k in range(len(tra_gps.get_parameter_vector())-1):
                param_vect.append(par[len(vel_files)*2  +7*npl  + rv_gp_npar  + 3*npl + N_transit_files*2 + 2 + k])

            tra_gps = init_dSHOKernel(param_vect)

        else:

            for k in range(len(tra_gps.get_parameter_vector())):
                param_vect.append(np.log(par[len(vel_files)*2  +7*npl  + rv_gp_npar  + 3*npl + N_transit_files*2 + 2 + k]))
            tra_gps.set_parameter_vector(np.array(param_vect))

        ##### data that will be modeled with a GP ######
        flux_err_all_gp = np.concatenate([flux_err[x] for x in range(len(flux_err)) if use_gp_model[x] == True])
        flux_o_c_all_gp = np.concatenate([flux_o_c[x] for x in range(len(flux_o_c)) if use_gp_model[x] == True])
        t_all_gp        = np.concatenate([t[x] for x in range(len(t)) if use_gp_model[x] == True])
        flux_model_all_ = np.concatenate([flux_model[x] for x in range(len(flux_model)) if use_gp_model[x] == True])     

        tra_gps.compute(t_all_gp,yerr=flux_err_all_gp)
        tr_loglik_gp = tr_loglik + tra_gps.log_likelihood(flux_o_c_all_gp)
 
        
        ##### data that will NOT be modeled with a GP ######

        flux_all_no_gp       = np.concatenate([flux[x] for x in range(len(flux)) if use_gp_model[x] != True])
        flux_model_all_no_gp = np.concatenate([flux_model[x] for x in range(len(flux_model)) if use_gp_model[x] != True])
        sig2i_all_no_gp      = np.concatenate([sig2i[x] for x in range(len(sig2i)) if use_gp_model[x] != True])
       
        flux_o_c_all_no_gp   = flux_all_no_gp - flux_model_all_no_gp
 
        tr_loglik_no_gp = tr_loglik -0.5*(np.sum((flux_o_c_all_no_gp)**2 * sig2i_all_no_gp  - np.log(sig2i_all_no_gp  / 2./ np.pi)))
                    

        tr_loglik = tr_loglik_gp + tr_loglik_no_gp

#        tra_gp_model_all = tra_gps.predict(flux_o_c_all , t_all , return_cov=False)
#        flux_o_c_gp_all = flux_o_c_all  - tra_gp_model_all       
#        tr_loglik2 = -0.5*(np.sum((flux_o_c_gp_all)**2           * sig2i_all  - np.log(sig2i_all  / 2./ np.pi))) # - np.log(sig2i / 2./ np.pi)
        if return_model == True:        
            flux_model_all_gp   = tra_gps.predict(flux_o_c_all_gp , t_all_gp , return_cov=False)
            flux_o_c_gp_all_gp  = flux_o_c_all_gp  - flux_model_all_gp  

            tra_gp_model_all    = np.concatenate([flux_model_all_gp+flux_model_all_, flux_model_all_no_gp])
            flux_o_c_gp_all     = np.concatenate([flux_o_c_gp_all_gp,  flux_o_c_all_no_gp])        
            

    if return_model == True:

        t_all = np.concatenate([tr_files[x][0] for x in range(60) if len(tr_files[x]) != 0])
        t_rich =np.linspace(min(t_all),max(t_all),len(t_all)*tra_model_fact)
        flux_model_rich = np.ones(len(t_rich))

        m  =  {k: [] for k in range(9)}

        tr_params.limb_dark = str(tr_model[0][j])
        tr_params.u = tr_model[1][j]


        for i in range(npl):

            if hkl == True:
                tr_params.ecc = np.sqrt(par[len(vel_files)*2 +7*i+2]**2 + par[len(vel_files)*2 +7*i+3]**2)
                tr_params.w  = np.degrees(np.arctan2(par[len(vel_files)*2 +7*i+2],par[len(vel_files)*2 +7*i+3]))%360
            else:
                tr_params.ecc = par[len(vel_files)*2 +7*i+2]
                tr_params.w   = par[len(vel_files)*2 +7*i+3]

            tr_params.per = par[len(vel_files)*2 +7*i+1]
            tr_params.inc = par[len(vel_files)*2 +7*i+5]

            tr_params.t0  = par[len(vel_files)*2 +7*npl +2 +rv_gp_npar + 3*i]
            tr_params.a   = par[len(vel_files)*2 +7*npl +2 +rv_gp_npar + 3*i+1]
            tr_params.rp  = par[len(vel_files)*2 +7*npl +2 +rv_gp_npar + 3*i+2]
            
            
            if program.endswith("dyn+") or program.endswith("dyn"):

                t_os = ttvs_mod(par,vel_files,npl, stmass, [epoch,ttv_times[1],max(t_rich)], i,hkl, fit_results=fit_results)
    
                for tran_t0 in t_os[1]:
                    tr_params.t0  = float(tran_t0)

                    
                    m[i] = batman.TransitModel(tr_params, t_rich)
                    tr_ind = np.where(np.logical_and(t_rich >= tran_t0-0.3, t_rich <= tran_t0+0.3))
                    flux_model_rich[tr_ind] = m[i].light_curve(tr_params)[tr_ind]
                    
            elif get_TTVs[0] == True:
                
                t_os = par[get_TTVs[1][i][0]:get_TTVs[1][i][1]]
                
                for tran_t0 in t_os:
                    tr_params.t0  = float(tran_t0)

                    m[i] = batman.TransitModel(tr_params, t_rich)
                    tr_ind = np.where(np.logical_and(t_rich >= tran_t0-0.3, t_rich <= tran_t0+0.3))
                    flux_model_rich[tr_ind] = m[i].light_curve(tr_params)[tr_ind]  
                    
            else:
                m[i] = batman.TransitModel(tr_params, t_rich)
                flux_model_rich = flux_model_rich * m[i].light_curve(tr_params) 


#        if tr_files[j][10] == True:
#            np.linspace(min(t_all),max(t_all),len(t_all)*tra_model_fact)
#            flux_model_rich = get_airmass_model(tr_files[j][3],flux_model_rich,0.0,
#                                                    par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*2 + tra_gp_npar + l],
#                                                    par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*3 + tra_gp_npar + l])           
#        else:
#            flux_model_rich = get_quad_model(flux_model_rich,flux_model_rich,0.0,
#                                                    par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*2 + tra_gp_npar + l],
#                                                    par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*3 + tra_gp_npar + l])



#        flux_model_ = flux_model_*tr_files[j][8]  + (1.0 - tr_files[j][8])             
 
        l = 0
        for j in range(60):

            if len(tr_files[j]) == 0:
                continue
            else:
                tr_ind = np.where(np.logical_and(t_rich >= min(tr_files[j][0]), t_rich <= max(tr_files[j][0])))

                flux_model_rich[tr_ind] =  (flux_model_rich[tr_ind]*tr_files[j][8]  + (1.0 - tr_files[j][8]))

                l +=1

        rich_model = np.array([t_rich,flux_model_rich],dtype=object)
        sep_data = np.array([t, flux,flux_err,flux_model, flux_o_c, flux_o_c_gp,tra_gp_model],dtype=object)
        all_data = np.array([t_all, flux_all,flux_err_all,flux_model_all, flux_o_c_all, flux_o_c_gp_all,tra_gp_model_all],dtype=object)

        tr_chi     = np.sum((flux_o_c_gp_all)**2 * sig2i_all ) # - np.log(sig2i / 2./ np.pi)
        tr_chi_red = tr_chi/len(flux_o_c_gp_all)

        tr_rms = np.sqrt(np.average(flux_o_c_gp_all**2))
        tr_wrms =  np.sqrt(np.average(flux_o_c_gp_all**2, weights=1/flux_err_all))
        tr_Ndata = len(flux_o_c_gp_all)

        tr_stat = [tr_chi,tr_chi_red,tr_rms,tr_wrms, tr_Ndata]

        return np.array([tr_loglik, sep_data, all_data,rich_model,tr_stat],dtype=object)
    else:
        return tr_loglik

 
#from threadpoolctl import threadpool_limits
#@threadpool_limits.wrap(limits=1, user_api='blas')
def model_loglik(p, program, par, flags, npl, vel_files, tr_files, tr_model, tr_params, epoch, stmass, gps, tra_gps, rtg, mix_fit, opt, outputfiles = [1,1,0], amoeba_starts=0, prior=0, eps='1.0E-8',dt=864000, when_to_kill=10, npoints=5000, model_max = 100, model_min =0):

    rv_loglik = 0
    gp_rv_loglik = 0
    tr_loglik = 0
    gp_tr_loglik = 0
    ttv_loglik = 0
    astr_loglik = 0

    dt = opt["dt"]
    eps = opt["eps"]
    when_to_kill = opt["when_to_kill"]
    copl_incl = opt["copl_incl"]
    jit_flag = opt["jit_flag"]    
    hkl = opt["hkl"]
    cwd = opt["cwd"]
    gr_flag = opt["gr_flag"]

    rvs_files = opt["RVS_files"]
    ttv_files = opt["TTV_files"]
    ast_files = opt["AST_files"]
    ast_files_hipp_gaia = opt["AST_files_hipp"]    
    ttv_times = opt["TTV_times"]
    ast_times = opt["AST_times"]
    get_TTVs  = opt["get_TTVs"]
#    re = opt["re"]

   
    N_transit_files = len([x for x in range(60) if len(tr_files[x]) != 0])

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

#    if(opt["TTV"]):
#        outputfiles = [1,1,1]        

    if hkl == True:
        for i in range(npl): 
            ecc_ = np.sqrt(par[len(vel_files)*2 +7*i+2]**2 + par[len(vel_files)*2 +7*i+3]**2)
            if ecc_ >= 0.95 or ecc_ < 0.0: 
                return -np.inf  


    if rtg[2] == True:
        

        for i in range(npl): 

            if hkl == True:
                ecc_ = np.sqrt(par[len(vel_files)*2 +7*i+2]**2 + par[len(vel_files)*2 +7*i+3]**2)
                om_  = np.degrees(np.arctan2(par[len(vel_files)*2 +7*i+2],par[len(vel_files)*2 +7*i+3]))%360
            else:
                ecc_, om_, = par[len(vel_files)*2 +7*i+2],par[len(vel_files)*2 +7*i+3]

            #par[len(vel_files)*2 +7*i+4] = ma_from_t0(par[len(vel_files)*2 +7*i+1],
            #                                      ecc_, om_, par[len(vel_files)*2 +7*i+6], par[len(vel_files)*2 +7*npl + 2 +rv_gp_npar + 3*i],epoch)
            #if program.endswith("kep"):
            #    par[len(vel_files)*2 +7*i+4] = get_m0(par[len(vel_files)*2 +7*i+1], ecc_, om_ ,par[len(vel_files)*2 +7*npl + 2 +rv_gp_npar + 3*i], epoch)

            if program.endswith("kep"):
                Ma_= get_m0(    par[len(vel_files)*2 +7*i+1], ecc_, om_ ,par[len(vel_files)*2 +7*npl + 2 +rv_gp_npar + 3*i], epoch)
                #Ma_= ma_from_t0(par[len(vel_files)*2 +7*i+1], ecc_, om_, par[len(vel_files)*2 +7*npl + 2 +rv_gp_npar + 3*i],epoch)
                if hkl == True:              
                    par[len(vel_files)*2 +7*i+4] = (Ma_ + om_)%360.0  
                else:
                    par[len(vel_files)*2 +7*i+4] =  Ma_
    
    else:
        for i in range(npl):
            if hkl == True:
                ecc_ = np.sqrt(par[len(vel_files)*2 +7*i+2]**2 + par[len(vel_files)*2 +7*i+3]**2)
                om_  = np.degrees(np.arctan2(par[len(vel_files)*2 +7*i+2],par[len(vel_files)*2 +7*i+3]))%360
                Ma_  = (par[len(vel_files)*2 +7*i+4] - om_)%360.0
            else:
                ecc_, om_, Ma_ = par[len(vel_files)*2 +7*i+2], par[len(vel_files)*2 +7*i+3], par[len(vel_files)*2 +7*i+4]

            par[len(vel_files)*2 +7*npl + 2 +rv_gp_npar + 3*i] = transit_tperi(par[len(vel_files)*2 +7*i+1],
                                                                  ecc_, om_, Ma_ ,epoch)[1]%par[len(vel_files)*2 +7*i+1]
 
###################################################################
###################################################################
###################################################################
    if(rtg[0]):

        rvmod = Rvfit()

        rv_off = []
        rv_jit = []
        for i in range(len(vel_files)):

            if (rtg[1]):
                rv_off.append([par[i],0])
                rv_jit.append([0,0])
            else:
                rv_off.append([par[i],0])
                rv_jit.append([par[i+ len(vel_files)],0])


######## Ugly fix with many limitations! incilation of planet 1 is given to all. TBD in fortran!!!! #############
#        if copl_incl == True:
#            incl_c = par[len(vel_files)*2 +7*i+5]
#            for i in range(npl):
#                par[len(vel_files)*2 +7*i+5] = incl_c 


        array_npl = []
        
        for i in range(npl): # K,P,e,w,M,i,cap0m for each planet, and information which ones we use
            array_npl.append([[par[len(vel_files)*2 + 7*i], 0],
                               [par[len(vel_files)*2 +7*i+1], 0],
                               [par[len(vel_files)*2 +7*i+2], 0],
                               [par[len(vel_files)*2 +7*i+3], 0],
                               [par[len(vel_files)*2 +7*i+4], 0],
                               [par[len(vel_files)*2 +7*i+5], 0],
                               [par[len(vel_files)*2 +7*i+6], 0],
                               [par[len(vel_files)*2 +7*npl  + rv_gp_npar  + 3*npl + N_transit_files*4 + tra_gp_npar + 2 + i ], 0],
                               [0, 0],
                               [0, 0],
                               [0, 0],
                               [0, 0],
                               [0, 0],
                               [0, 0],
                               [par[len(vel_files)*2 +7*i+2], 0],
                               [par[len(vel_files)*2 +7*i+3], 0],
                               [par[len(vel_files)*2 +7*i+4], 0]])

        
        rvmod.init_args(rvs_files, array_npl, epoch, hkl,
                     get_RV=outputfiles[0], get_best_par=outputfiles[1], get_fit_model=outputfiles[2],
                     rv_jitt=rv_jit, rv_ofset=rv_off,
                     lin_trend_in=[par[len(vel_files) * 2 + 7 * npl], 0],
                     quad_trend_in=[par[len(vel_files) * 2 + 7 * npl + 1], 0],
                     dyn_eps=eps, dyn_dt=dt, amoeba_iter=amoeba_starts, timeout=when_to_kill,
                     rv_model_npoints=npoints, rv_model_max=model_max, rv_model_min=model_min,
                     rv_gr_flag=int(gr_flag), stellar_mass=stmass,
                     ndset=None, ndata=None, nplanet=None,
                     dyn_planets=None,coplar_inc=int(copl_incl),jit_flag = int(jit_flag))

        if program == '%s/lib/fr/loglik_kep'%cwd:
            rvmod.run_amoeba("kep")         
        elif program == '%s/lib/fr/loglik_dyn'%cwd:
            rvmod.run_amoeba("dyn")

        rv_loglik = float(rvmod.loglik)
         
        if(rtg[1]):

            gp_rv_loglik = 0

            param_vect = []

            if opt["RV_GP_kernel"] == "dSHOKernel":
               # print(len(gps.get_parameter_vector()))
                for j in range(len(gps.get_parameter_vector())-1):
                    if opt["link_RV_GP"][j]==True and rtg[3] == True:
                        param_vect.append(par[len(vel_files)*2  +7*npl  + rv_gp_npar  + 3*npl + N_transit_files*2 + 2 + j])
                    else:
                        param_vect.append(par[len(vel_files)*2  +7*npl +2 +j])

                gps = init_dSHOKernel(param_vect)

            else:

                for j in range(len(gps.get_parameter_vector())):
                    if opt["link_RV_GP"][j]==True and rtg[3] == True:
                        param_vect.append(np.log(par[len(vel_files)*2  +7*npl  + rv_gp_npar  + 3*npl + N_transit_files*2 + 2 + j]))
                    else:
                        param_vect.append(np.log(par[len(vel_files)*2  +7*npl +2 +j]))
         

                gps.set_parameter_vector(np.array(param_vect))

            errors_with_jitt = np.array([np.sqrt(rvmod.rv_err[i]**2 + par[int(ii) + len(vel_files)]**2)  for i,ii in enumerate(rvmod.idset)])
            gps.compute(rvmod.jd,yerr=errors_with_jitt)
            gp_rv_loglik = gp_rv_loglik + gps.log_likelihood(rvmod.o_c)

            rv_loglik =  gp_rv_loglik
        
         
        
    else:
        rv_loglik = 0


###################################################################
 
    if(rtg[2]):
        
        if N_transit_files == 0:
            tr_loglik = 0
        elif rtg[0] ==False:      
            tr_loglik = transit_loglik(program, tr_files,vel_files,tr_params,tr_model,par,rv_gp_npar,tra_gp_npar,npl,hkl,rtg,tra_gps,stmass,ttv_times,epoch,get_TTVs,opt,fit_results=False )
        else:
            tr_loglik = transit_loglik(program, tr_files,vel_files,tr_params,tr_model,par,rv_gp_npar,tra_gp_npar,npl,hkl,rtg,tra_gps,stmass,ttv_times,epoch,get_TTVs,opt,fit_results=rvmod )

 
    if(opt["AST"]):
        if(rtg[0])==False:
            astr_loglik = ast_loglik(par,vel_files,ast_files,npl,stmass,ast_times,hkl,fit_results=False, return_model = False)
            astr_loglik = astr_loglik + ast_loglik_hipp(par,vel_files,ast_files_hipp_gaia,npl,stmass,ast_times,hkl,fit_results=False, return_model = False)  
        else:
            astr_loglik = ast_loglik(par,vel_files,ast_files,npl,stmass,ast_times,hkl,fit_results=rvmod, return_model = False)
            astr_loglik = astr_loglik + ast_loglik_hipp(par,vel_files,ast_files_hipp_gaia,npl,stmass,ast_times,hkl,fit_results=rvmod, return_model = False)              


    if(opt["TTV"]):
        if(rtg[0])==False:
            ttv_loglik = ttvs_loglik(par,vel_files,ttv_files,npl,stmass,ttv_times,hkl,fit_results=False, return_model = False)
        else:
            ttv_loglik = ttvs_loglik(par,vel_files,ttv_files,npl,stmass,ttv_times,hkl,fit_results=rvmod, return_model = False)
 
    if opt["AMD_stab"] == True and npl >=2:
               
        for i in range(npl - 1):

            pl_mass_in,ap_in = mass_a_from_Kepler_fit([par[len(vel_files)*2 + 7*i],
                                par[len(vel_files)*2 +7*i+1],
                                par[len(vel_files)*2 +7*i+2],
                                par[len(vel_files)*2 +7*i+3],
                                par[len(vel_files)*2 +7*i+4]],1,stmass) ##################TB FIXED! 

            pl_mass_out,ap_out = mass_a_from_Kepler_fit([par[len(vel_files)*2 + 7*(i+1)],
                                par[len(vel_files)*2 +7*(i+1)+1],
                                par[len(vel_files)*2 +7*(i+1)+2],
                                par[len(vel_files)*2 +7*(i+1)+3],
                                par[len(vel_files)*2 +7*(i+1)+4]],1,stmass) ##################TB FIXED! 

            alpha    = ap_in/ap_out
            gamma    = pl_mass_in/pl_mass_out
            epsilon  = (pl_mass_in + pl_mass_out)/(stmass*1047.5654817267318)

            AMD = gamma*np.sqrt(alpha)*(1.-np.sqrt(1.- par[len(vel_files)*2 +7*i+2]**2)) + 1.- np.sqrt(1.- par[len(vel_files)*2 + 7*(i+1) +2]**2)

            AMD_Hill = gamma*np.sqrt(alpha) + 1. - (1.+gamma)**1.5 * np.sqrt(alpha/(gamma+alpha) * (1.+(3.**(4./3.)*epsilon**(2./3.)*gamma)/((1.+gamma)**2)))
 
            if AMD >= AMD_Hill:
                return (rv_loglik + tr_loglik + ttv_loglik  + astr_loglik)* 2.0*np.exp(1.0 - AMD_Hill/AMD)
 


    if rtg[0]:
        del rvmod
            
    #print(   rv_loglik, tr_loglik,ttv_loglik,astr_loglik     )
    if np.isnan(rv_loglik).any() or np.isnan(tr_loglik).any():
        return -np.inf
    return rv_loglik + tr_loglik + ttv_loglik +  astr_loglik


def run_Fort(obj,mod,minimize_loglik=True):


    import threading
 
    start_time = time.time()

    if(minimize_loglik):
        t1 = threading.Thread(target=lambda:obj.rvmod.run_amoeba(mod))
    else:
        t1 = threading.Thread(target=lambda:obj.rvmod.run_lm(mod))

    t1.start()
    t1.join()

    print("--- %s seconds ---" % (time.time() - start_time))
    obj.flag = 1
    return obj
    
    
def run_SciPyOp(obj,   threads=1,  kernel_id=-1,  save_means=False, fileoutput=False, save_sampler=False, **kwargs):

    start_time = time.time()
    rtg = obj.rtg

#    check_temp_RV_file(obj)

    vel_files = [0]*obj.ndset

    N_transit_files = len([x for x in range(60) if len(obj.tra_data_sets[x]) != 0])

    tr_files = obj.tra_data_sets
    tr_mo    = obj.ld_m
    tr_ld    = obj.ld_u
    tr_gr    = obj.ld_gr
    tr_gr_ind= obj.ld_gr_ind
    
    tr_model = np.array([tr_mo,tr_ld,tr_gr,tr_gr_ind], dtype=object)

    tr_params = obj.tr_params


    final_array = []
    jd, rvs, sig, ids = [], [], [], []


    for z in range(len(obj.rv_data_sets)):
        if len(obj.rv_data_sets[z]) == 0:
            continue
 
        jd = np.append(jd, obj.rv_data_sets[z][0], axis = 0)
        rvs = np.append(rvs, obj.rv_data_sets[z][1], axis = 0)
        sig = np.append(sig, obj.rv_data_sets[z][2], axis = 0)
        ids = np.append(ids, obj.rv_data_sets[z][3]+1, axis = 0)
    final_array = np.array([jd, rvs, sig, ids]).T

    rvs_files  = final_array 
    ttv_files  = obj.ttv_data_sets
    ast_files  = obj.ast_data_sets
    ast_files2 = obj.ast_data_sets_hipp_gaia

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
    pp = obj.par_for_mcmc
    ee = obj.e_for_mcmc 
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
    
#    re = Rvfit()
    opt = {"eps":obj.dyn_model_accuracy*1e-13,
           "dt":obj.time_step_model*86400.0,
           "when_to_kill":when_to_kill,
           "copl_incl":obj.copl_incl,
           "jit_flag":obj.jit_flag,
           "hkl":obj.hkl,
           "rho":obj.rho,                      
           "cwd":obj.cwd, 
           "gr_flag":obj.gr_flag,
           "TTV":obj.type_fit["TTV"],
           "AST":obj.type_fit["AST"],
           "RVS_files":rvs_files,
           "TTV_files":ttv_files, 
           "AST_files":ast_files, 
           "AST_files_hipp":ast_files2,            
           "TTV_times":obj.ttv_times,
           "AST_times":obj.ast_times,
           "AMD_stab":obj.optim_AMD_stab, 
           "Nbody_stab":obj.optim_Nbody_stab,
           "get_TTVs":obj.get_TTVs,
           "link_RV_GP":obj.link_RV_GP,
           "tra_GP_kernel":obj.tra_gp_kernel,
           "RV_GP_kernel":obj.gp_kernel}
#            "re":re}

    if obj.init_fit == True:
        flags = []

    gps = []
    if (rtg[1]):
        initiate_RV_gps(obj, kernel_id=kernel_id)
        gps = obj.gps
        rv_gp_npar = len(gps.get_parameter_vector())
    else:
        rv_gp_npar = 0

    tra_gps = []
    if (rtg[3]) and N_transit_files != 0:
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
    elif obj.SciPy_min_use_1 == obj.SciPy_min[7]:
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
    elif obj.SciPy_min_use_2 == obj.SciPy_min[7]:
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

    ############ TESTS begin ####################
    #bounds = [(-1.00e+05, 1.00e+05),(-1.00e+05, 1.00e+05),(-1.00e+05, 1.00e+05),(-1.00e+05, 1.00e+05),
  #            (-1.00e+05, 1.00e+05),(-1.00e+05, 1.00e+05),(-1.00e+05, 1.00e+05),(-1.00e+05, 1.00e+05),
  #            (-1.00e+05, 1.00e+05),(-1.00e+05, 1.00e+05),(-1.00e+05, 1.00e+05),(-1.00e+05, 1.00e+05)]
    
    #from mystic.samplers import BuckshotSampler
    #solution = BuckshotSampler(bounds, nll,disp=1,retall=1,args=(mod, par,flags, npl,vel_files, 
    #                                                         tr_files, tr_model, tr_params,  epoch, stmass, bb, priors, gps, tra_gps, rtg, mix_fit, opt ))
    #obj.mys = solution                         
   # import mystic
   # from mystic.solvers import fmin
   # solution = fmin(nll,pp,disp=0,retall=1,args=(mod, par,flags, npl,vel_files, tr_files, tr_model, tr_params,  epoch, stmass, bb, priors, gps, tra_gps, rtg, mix_fit, opt ),
    #                         method=method1,bounds=fit_bounds, options=options1)
    #allvecs = solution[-1]
    
   # print(allvecs[])
   ############ TESTS end ####################


#    print(pp)
#    print(mod, par,flags, npl,vel_files, tr_files, tr_model, tr_params,  epoch, stmass, bb, priors, gps, tra_gps, rtg, mix_fit, opt )

    ########################### Primary minimizer #########################
    for k in range(n1): # run at least 3 times the minimizer
        #eps = eps/10.0
       # print('running %s %s %s'%(obj.SciPy_min_use_1, obj.SciPy_min_N_use_1, k))
        result = op.minimize(nll,  pp, args=(mod, par,flags, npl,vel_files, tr_files, tr_model, tr_params,  epoch, stmass, bb, priors, gps, tra_gps, rtg, mix_fit, opt ),
                             method=method1,bounds=fit_bounds, options=options1)
                            #  bounds=bb, tol=None, callback=None, options={'eps': 1e-08, 'scale': None, 'offset': None, 'mesg_num': None, 'maxCGit': -1, 'maxiter': None, 'eta': -1, 'stepmx': 0, 'accuracy': 0, 'minfev': 0, 'ftol': -1, 'xtol': -1, 'gtol': -1, 'rescale': -1, 'disp': True})
        pp = result["x"]
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

    #print(pp)
    errors = [[0.0,0.0] for i in range(len(pp))]

    obj.par_for_mcmc = pp
#    newparams = obj.generate_newparams_for_mcmc(obj.par_for_mcmc)

#    if obj.copl_incl == True:
#        incl_c = par[len(vel_files)*2 +7*0+5]
#        for i in range(npl):
#            par[len(vel_files)*2 +7*i+5] = incl_c 
#            obj.i[i] = incl_c

#    obj.overwrite_params(newparams)
#    obj.correct_elements()
#    obj.hack_around_rv_params()

    obj = return_results(obj, pp, ee, par, flags, npl, vel_files, tr_files, tr_model, tr_params, epoch, stmass, bb, priors, gps, tra_gps, rtg, mix_fit, errors, mod,opt)

  #  obj.init_fit = False
    print("--- %s seconds ---" % (time.time() - start_time))

    return obj





def return_results(obj, pp, ee, par,flags, npl,vel_files, tr_files, tr_model, tr_params, epoch, stmass, bb, pr_nr, gps, tra_gps, rtg, mix_fit, errors, mod,opt):


    N_transit_files = len([x for x in range(60) if len(tr_files[x]) != 0]) 
    
    get_TTVs = opt["get_TTVs"]
    
    e_par = [[0.0,0.0]]*len(par)
    for j in range(len(flags)):
        par[flags[j]]   = pp[j]
        e_par[flags[j]] = [errors[j][0],errors[j][1]]


    for i in range(len(vel_files)):

        obj.rvoff[i] = par[i]
        obj.jitt[i] = par[i+ len(vel_files)]


#        if obj.copl_incl == True:
#            incl_c = par[len(vel_files)*2 +7*0+5]
#            obj.i[i] = incl_c
 

    for i in range(obj.npl):

        obj.K[i]  = par[len(vel_files)*2 + 7*i]
        obj.P[i]  = par[len(vel_files)*2 + 7*i+1]
        obj.i[i]  = par[len(vel_files)*2 + 7*i+5]
        obj.Node[i] = par[len(vel_files)*2 + 7*i+6]

        if obj.hkl == True:

            obj.e_sinw[i]  = par[len(vel_files)*2 + 7*i+2]
            obj.e_cosw[i]  = par[len(vel_files)*2 + 7*i+3]
            obj.lamb[i]    = par[len(vel_files)*2 + 7*i+4]

            obj.e[i]   = np.sqrt(obj.e_sinw[i]**2 + obj.e_cosw[i]**2)
            obj.w[i]   = np.degrees(np.arctan2(np.radians(obj.e_sinw[i]),np.radians(obj.e_cosw[i])))
            obj.M0[i]  = (obj.lamb[i] - obj.w[i])%360.0
            

        else:

            obj.e[i]  = par[len(vel_files)*2 + 7*i+2]
            obj.w[i]  = par[len(vel_files)*2 + 7*i+3]
            obj.M0[i] = par[len(vel_files)*2 + 7*i+4]
            obj.e_sinw[i] = obj.e[i]*np.sin(np.radians(obj.w[i]))
            obj.e_cosw[i] = obj.e[i]*np.cos(np.radians(obj.w[i]))
            obj.lamb[i]   = (obj.w[i] + obj.M0[i])%360.0
 

        obj.t_peri[i] = transit_tperi(obj.P[i],obj.e[i], obj.w[i], obj.M0[i] ,obj.epoch)[0]

    obj.rv_lintr = par[len(vel_files)*2  +7*npl]
    obj.rv_quadtr = par[len(vel_files)*2  +7*npl +1]
 

    if (rtg[1]):
        if obj.gp_kernel == 'RotKernel':
            
            rv_gp_npar =4 
                
            for j in range(len(gps.get_parameter_vector())):
        
                if opt["link_RV_GP"][j]==True and rtg[3] == True and obj.gp_kernel == obj.tra_gp_kernel:
                    obj.GP_rot_params[j] = par[len(vel_files)*2  +7*npl  + rv_gp_npar  + 3*npl + N_transit_files*2 + 2 + j]
                else:
                    obj.GP_rot_params[j] = par[len(vel_files)*2  +7*npl +2 +j]
 
    
        elif obj.gp_kernel == 'SHOKernel':
            
            rv_gp_npar =3 
            
            for j in range(len(gps.get_parameter_vector())):
    
                if opt["link_RV_GP"][j]==True and rtg[3] == True and obj.gp_kernel == obj.tra_gp_kernel:
                    obj.GP_sho_params[j] = par[len(vel_files)*2  +7*npl  + rv_gp_npar  + 3*npl + N_transit_files*2 + 2 + j]
                else:
                    obj.GP_sho_params[j] = par[len(vel_files)*2  +7*npl +2 +j]
    
        elif obj.gp_kernel == 'Matern32':
            for j in range(len(gps.get_parameter_vector())):
                obj.GP_mat_params[j] = par[len(vel_files)*2  +7*npl +2 +j]
            rv_gp_npar =3 
             
        elif obj.gp_kernel == 'RealTerm':
            for j in range(len(gps.get_parameter_vector())):
                obj.GP_drw_params[j] = par[len(vel_files)*2  +7*npl +2 +j]
            rv_gp_npar =2                

        elif obj.gp_kernel == 'dSHOKernel':
            rv_gp_npar =6
            
            for j in range(len(gps.get_parameter_vector())-1):
    
                if opt["link_RV_GP"][j]==True and rtg[3] == True and obj.gp_kernel == obj.tra_gp_kernel:
                    obj.GP_double_sho_params[j] = par[len(vel_files)*2  +7*npl  + rv_gp_npar  + 3*npl + N_transit_files*2 + 2 + j]
                else:
                    obj.GP_double_sho_params[j] = par[len(vel_files)*2  +7*npl +2 +j]
    else:
        rv_gp_npar = 0



    if (rtg[3]) and N_transit_files != 0:
        if obj.tra_gp_kernel == 'RotKernel':
            for j in range(len(tra_gps.get_parameter_vector())):
                obj.tra_GP_rot_params[j] = par[len(vel_files)*2  +7*npl  + rv_gp_npar  + 3*npl + N_transit_files*2 + 2 +j]
            tra_gp_npar = 4
        
        elif obj.tra_gp_kernel == 'SHOKernel':
            for j in range(len(tra_gps.get_parameter_vector())):
                obj.tra_GP_sho_params[j] = par[len(vel_files)*2  +7*npl  + rv_gp_npar  + 3*npl + N_transit_files*2 + 2 +j]
            tra_gp_npar = 3

        elif obj.tra_gp_kernel == 'Matern32':
            for j in range(len(tra_gps.get_parameter_vector())):
                obj.tra_GP_mat_params[j] = par[len(vel_files)*2  +7*npl  + rv_gp_npar  + 3*npl + N_transit_files*2 + 2 +j]
            tra_gp_npar = 3

        elif obj.tra_gp_kernel == 'RealTerm':
            for j in range(len(tra_gps.get_parameter_vector())):
                obj.tra_GP_drw_params[j] = par[len(vel_files)*2  +7*npl  + rv_gp_npar  + 3*npl + N_transit_files*2 + 2 +j]
            tra_gp_npar = 2

        elif obj.tra_gp_kernel == 'dSHOKernel':
            for j in range(len(tra_gps.get_parameter_vector())-1):
                obj.tra_GP_double_sho_params[j] = par[len(vel_files)*2  +7*npl  + rv_gp_npar  + 3*npl + N_transit_files*2 + 2 +j]
            tra_gp_npar = 6
        #tra_gp_npar = len(tra_gps.get_parameter_vector())
    else:
        tra_gp_npar = 0


    if opt['rho'] == True:
        
        rho = par[len(vel_files)*2 +7*npl +2 +rv_gp_npar + 3*0+1]
        
        for i in range(npl):  
            aR = get_aR_from_rho(rho, par[len(vel_files)*2 +7*i+1])
 
            print("Planet %s, aR =%s"%(i+1,aR))
            obj.t0[i]     = par[len(vel_files)*2 +7*npl +2 +rv_gp_npar + 3*i]
            obj.pl_a[i]   = rho
            obj.pl_rad[i] = par[len(vel_files)*2 +7*npl +2 +rv_gp_npar + 3*i+2]
             
            
            
    else:
        for i in range(obj.npl):
     
            obj.t0[i]     = par[len(vel_files)*2 +7*npl +2 +rv_gp_npar + 3*i]
            obj.pl_a[i]   = par[len(vel_files)*2 +7*npl +2 +rv_gp_npar + 3*i+1]
            obj.pl_rad[i] = par[len(vel_files)*2 +7*npl +2 +rv_gp_npar + 3*i+2]
 
 
    j = 0
    for i in range(60):
        if len(obj.tra_data_sets[i]) == 0:
            continue
        else:
            obj.tra_off[i] =      par[len(vel_files)*2 +7*npl + 2 +rv_gp_npar + 3*npl + j]
            obj.tra_jitt[i] = abs(par[len(vel_files)*2 +7*npl + 2 +rv_gp_npar + 3*npl + N_transit_files + j])
             
            obj.tra_lintr[i] =      par[len(vel_files)*2 +7*npl + 2 +rv_gp_npar + 3*npl + N_transit_files*2 + tra_gp_npar + j]
            obj.tra_quadtr[i] =     par[len(vel_files)*2 +7*npl + 2 +rv_gp_npar + 3*npl + N_transit_files*3 + tra_gp_npar + j] 
 
            j = j +1



        if tr_model[0][i] == "linear":
            obj.ld_u_lin[i] = [par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*4 + tra_gp_npar + tr_model[3][i] +  npl]]


        elif tr_model[0][i] == "quadratic":
            obj.ld_u_quad[i] = [par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*4 + tra_gp_npar + tr_model[3][i] + npl], 
                                par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*4 + tra_gp_npar + tr_model[3][i] + 1+ npl]]
 

        elif tr_model[0][j] == "nonlinear":
            obj.ld_u_nonlin[i] = [par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*4 + tra_gp_npar + tr_model[3][i] + npl], 
                                  par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*4 + tra_gp_npar + tr_model[3][i] + 1+ npl],
                                  par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*4 + tra_gp_npar + tr_model[3][i] + 2+ npl], 
                                  par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*4 + tra_gp_npar + tr_model[3][i] + 3+ npl]]
 


    if obj.type_fit["AST"] and obj.use_ast_hipp_gaia == True: 

        obj.ast_alpha[0]    =    par[-6]
        obj.ast_delta[0]    =    par[-5]
        obj.ast_pi[0]       =    par[-4]       
        obj.ast_mu_alpha[0] =    par[-3]       
        obj.ast_mu_delta[0] =    par[-2]

    #print(par)

#################### Get the models for plotting ############################



    if obj.type_fit["RV"] == True and obj.type_fit["Transit"] == False and obj.type_fit["TTV"] == False and obj.type_fit["AST"] == False:

        obj.fitting(outputfiles=[1,1,1], minimize_fortran=True, minimize_loglik=True, amoeba_starts=0, doGP=False, npoints= obj.model_npoints, eps=float(opt["eps"])/1e-13, dt=float(opt["dt"])/86400.0)

        if rtg[1]:
            get_RV_gps_model(obj, get_lnl=True)
 

    elif obj.type_fit["RV"] == False and obj.type_fit["Transit"] == True:

        for i in range(npl): 


            if obj.hkl == True:
                ecc_ = np.sqrt(par[len(vel_files)*2 +7*i+2]**2 + par[len(vel_files)*2 +7*i+3]**2)
                om_  = np.degrees(np.arctan2(par[len(vel_files)*2 +7*i+2],par[len(vel_files)*2 +7*i+3]))%360
            else:
                ecc_, om_, = par[len(vel_files)*2 +7*i+2],par[len(vel_files)*2 +7*i+3]

            if mod.endswith("kep"):
                Ma_= get_m0(par[len(vel_files)*2 +7*i+1], ecc_, om_ ,par[len(vel_files)*2 +7*npl + 2 +rv_gp_npar + 3*i], epoch)
                #Ma_= ma_from_t0(par[len(vel_files)*2 +7*i+1], ecc_, om_, par[len(vel_files)*2 +7*npl + 2 +rv_gp_npar + 3*i],epoch)
                if obj.hkl == True:              
                    par[len(vel_files)*2 +7*i+4] = (Ma_ + om_)%360.0  
                else:
                    par[len(vel_files)*2 +7*i+4] =  Ma_



        obj.transit_results = transit_loglik(mod, tr_files,vel_files,tr_params,tr_model,par,rv_gp_npar,tra_gp_npar,obj.npl,obj.hkl, obj.rtg, obj.tra_gps,stmass, obj.ttv_times, obj.epoch,get_TTVs,opt,return_model = True, tra_model_fact=obj.tra_model_fact)

        obj.loglik = obj.transit_results[0]

        obj.fit_results.Ndata    = obj.transit_results[4][4]
        obj.fit_results.stat.dof = obj.transit_results[4][4] - len(pp)
        
        obj.fit_results.chi2 = obj.transit_results[4][0]
       # obj.fit_results.reduced_chi2 = obj.transit_results[4][1]
        obj.fit_results.reduced_chi2 = obj.fit_results.chi2 / obj.fit_results.stat.dof
        
        obj.fit_results.rms = obj.transit_results[4][2]
        obj.fit_results.wrms = obj.transit_results[4][3]



    elif obj.type_fit["RV"] == True and obj.type_fit["Transit"] == True:

        for i in range(npl): 

            if obj.hkl == True:
                ecc_ = np.sqrt(par[len(vel_files)*2 +7*i+2]**2 + par[len(vel_files)*2 +7*i+3]**2)
                om_  = np.degrees(np.arctan2(par[len(vel_files)*2 +7*i+2],par[len(vel_files)*2 +7*i+3]))%360
            else:
                ecc_, om_, = par[len(vel_files)*2 +7*i+2],par[len(vel_files)*2 +7*i+3]

            if mod.endswith("kep"):
                Ma_= get_m0(par[len(vel_files)*2 +7*i+1], ecc_, om_ ,par[len(vel_files)*2 +7*npl + 2 +rv_gp_npar + 3*i], epoch)
                #Ma_= ma_from_t0(par[len(vel_files)*2 +7*i+1], ecc_, om_, par[len(vel_files)*2 +7*npl + 2 +rv_gp_npar + 3*i],epoch)

                if obj.hkl == True:              
                    par[len(vel_files)*2 +7*i+4] = (Ma_ + om_)%360.0  
                else:
                    par[len(vel_files)*2 +7*i+4] =  Ma_

                #obj.params.update_M0(i,par[len(vel_files)*2 +7*i+4])
                obj.M0[i] = Ma_


        obj.fitting(outputfiles=[1,1,1], minimize_fortran=True, minimize_loglik=True, amoeba_starts=0, doGP=False, npoints= obj.model_npoints, eps=float(opt["eps"])/1e-13, dt=float(opt["dt"])/86400.0)
        
        obj.transit_results = transit_loglik(mod, tr_files,vel_files,tr_params,tr_model,par,rv_gp_npar,tra_gp_npar,obj.npl,obj.hkl, obj.rtg , obj.tra_gps, stmass, obj.ttv_times, obj.epoch, get_TTVs, opt, return_model = True, tra_model_fact=obj.tra_model_fact, fit_results=False)
        
        if rtg[1]:
            get_RV_gps_model(obj, get_lnl=True)

       # if rtg[3]:
       #     get_transit_gps_model(obj)
 
        obj.loglik     =   obj.loglik +  obj.transit_results[0]

 
    elif obj.type_fit["RV"] == True and obj.type_fit["TTV"] == True:
        
        obj.fitting(outputfiles=[1,1,1], minimize_fortran=True, minimize_loglik=True, amoeba_starts=0, doGP=False, npoints= obj.model_npoints, eps=float(opt["eps"])/1e-13, dt=float(opt["dt"])/86400.0)
 
        ttv_loglik = ttvs_loglik(obj.parameters,vel_files, obj.ttv_data_sets,obj.npl,obj.params.stellar_mass,obj.ttv_times,obj.hkl,
                                 fit_results=obj.fit_results, return_model = False)
        #print(obj.fit_results.mass, obj.ttv_times, obj.parameters, obj.loglik, ttv_loglik)


        print(obj.loglik +  ttv_loglik,obj.loglik, ttv_loglik)
        if rtg[1]:
            get_RV_gps_model(obj, get_lnl=True)

        obj.loglik     =   obj.loglik +  ttv_loglik

    elif obj.type_fit["RV"] == True and obj.type_fit["AST"] == True:
        
        obj.fitting(outputfiles=[1,1,1], minimize_fortran=True, minimize_loglik=True, amoeba_starts=0, doGP=False, npoints= obj.model_npoints, eps=float(opt["eps"])/1e-13, dt=float(opt["dt"])/86400.0)

        astr_loglik2 = ast_loglik(obj.parameters,vel_files, obj.ast_data_sets,obj.npl,obj.params.stellar_mass,obj.ast_times,obj.hkl, fit_results=obj.fit_results, return_model = False)
        astr_loglik2 = astr_loglik2 + ast_loglik_hipp(obj.parameters,vel_files, obj.ast_data_sets_hipp_gaia,obj.npl,obj.params.stellar_mass,obj.ast_times,obj.hkl, fit_results=obj.fit_results, return_model = False)

        #astr_loglik = ast_loglik(par,           vel_files, ast_files,            npl,                 stmass,    ttv_times,    hkl, fit_results=False, return_model = False)
        #print(astr_loglik2, obj.fit_results.mass, obj.ttv_times, obj.parameters, obj.loglik)

        if rtg[1]:
            get_RV_gps_model(obj, get_lnl=True)

        obj.loglik     =   obj.loglik +  astr_loglik2


    elif obj.type_fit["RV"] == False and obj.type_fit["Transit"] == False and obj.type_fit["TTV"] == True:
        ttv_loglik = ttvs_loglik(par,vel_files,obj.ttv_data_sets,npl,stmass,obj.ttv_times,obj.hkl,fit_results=False, return_model = False)
        obj.loglik     =  ttv_loglik


    elif obj.type_fit["RV"] == False and obj.type_fit["Transit"] == False and obj.type_fit["AST"] == True:
        astr_loglik = ast_loglik(par,vel_files,obj.ast_data_sets,npl,stmass,obj.ast_times,obj.hkl,fit_results=False, return_model = False)
        astr_loglik = astr_loglik + ast_loglik_hipp(par,vel_files,obj.ast_data_sets_hipp_gaia,npl,stmass,obj.ast_times,obj.hkl,fit_results=False, return_model = False)
        obj.loglik     =  astr_loglik
        #print(obj.fit_results.mass, obj.ttv_times, obj.parameters, obj.loglik, astr_loglik)



############# Errors here ##############

    for i in range(len(vel_files)):

        obj.rvoff_err[i] = e_par[i]
        obj.jitt_err[i] = e_par[i+ len(vel_files)]


#        if obj.copl_incl == True:
#            incl_err = e_par[len(vel_files)*2 +7*0+5]
#            obj.i_err[i] = incl_err
 


 
    for i in range(obj.npl):


        obj.K_err[i] = e_par[len(vel_files)*2 + 7*i]
        obj.P_err[i] = e_par[len(vel_files)*2 + 7*i+1]
        obj.i_err[i] = e_par[len(vel_files)*2 + 7*i+5]
        obj.Node_err[i] = e_par[len(vel_files)*2 + 7*i+6]

        obj.e_err[i] = e_par[len(vel_files)*2 + 7*i+2]
        obj.w_err[i] = e_par[len(vel_files)*2 + 7*i+3]
        obj.M0_err[i] = e_par[len(vel_files)*2 + 7*i+4]
 
        if obj.hkl == True:

            obj.e_sinw_err[i]  = e_par[len(vel_files)*2 + 7*i+2]
            obj.e_cosw_err[i]  = e_par[len(vel_files)*2 + 7*i+3]
            obj.lamb_err[i]    = e_par[len(vel_files)*2 + 7*i+4]

            obj.e_err[i]   = np.array([0.0,0.0]) #.sqrt(obj.e_sinw[i]**2 + obj.e_cosw[i]**2)
            obj.w_err[i]   = np.array([0.0,0.0])   # np.degrees(np.arctan2(np.radians(obj.e_sinw[i]),np.radians(obj.e_cosw[i])))
            obj.M0_err[i]  = np.array([0.0,0.0])   #(obj.lamb[i] - obj.w[i])%360.0

        else:

            obj.e_err[i] = e_par[len(vel_files)*2 + 7*i+2]
            obj.w_err[i] = e_par[len(vel_files)*2 + 7*i+3]
            obj.M0_err[i] = e_par[len(vel_files)*2 + 7*i+4]

            obj.e_sinw_err[i] = np.array([0.0,0.0]) #obj.e[i]*np.sin(np.radians(obj.w[i]))
            obj.e_cosw_err[i] = np.array([0.0,0.0])  #obj.e[i]*np.cos(np.radians(obj.w[i]))
            obj.lamb_err[i]   = np.array([0.0,0.0])  #(obj.w[i] + obj.M0[i])%360.0




    obj.rv_lintr_err = e_par[len(vel_files)*2  +7*npl]
    obj.rv_quadtr_err = e_par[len(vel_files)*2  +7*npl +1]     

 
 
    for i in range(npl):
         
        obj.t0_err[i]     = e_par[len(vel_files)*2 +7*npl +2 +rv_gp_npar + 3*i]
        obj.pl_a_err[i]   = e_par[len(vel_files)*2 +7*npl +2 +rv_gp_npar + 3*i+1]
        obj.pl_rad_err[i] = e_par[len(vel_files)*2 +7*npl +2 +rv_gp_npar + 3*i+2]

        obj.omega_dot[i] = par[len(vel_files)*2  +7*npl  + rv_gp_npar  + 3*npl + N_transit_files*4 + tra_gp_npar + 2 + i ]

    j = 0

    for i in range(60):
        if len(obj.tra_data_sets[i]) == 0:
            continue
        else:
            
            obj.tra_off_err[i] =      e_par[len(vel_files)*2 +7*npl + 2 +rv_gp_npar + 3*npl + j]
            obj.tra_jitt_err[i] =     e_par[len(vel_files)*2 +7*npl + 2 +rv_gp_npar + 3*npl + N_transit_files + j]
             
            obj.tra_lintr_err[i] =      e_par[len(vel_files)*2 +7*npl + 2 +rv_gp_npar + 3*npl + N_transit_files*2 + tra_gp_npar + j]
            obj.tra_quadtr_err[i] =     e_par[len(vel_files)*2 +7*npl + 2 +rv_gp_npar + 3*npl + N_transit_files*3 + tra_gp_npar + j]            

            j = j +1



        if tr_model[0][i] == "linear":
            obj.ld_u_lin_err[i] = [e_par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*4 + tra_gp_npar + tr_model[3][i] +  npl]]

        elif tr_model[0][i] == "quadratic":
            obj.ld_u_quad_err[i] = [e_par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*4 + tra_gp_npar + tr_model[3][i] + npl], 
                                    e_par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*4 + tra_gp_npar + tr_model[3][i] + 1+ npl]]

        elif tr_model[0][j] == "nonlinear":
            obj.ld_u_nonlin_err[i] = [e_par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*4 + tra_gp_npar + tr_model[3][i] + npl], 
                                  e_par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*4 + tra_gp_npar + tr_model[3][i] + 1+ npl],
                                  e_par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*4 + tra_gp_npar + tr_model[3][i] + 2+ npl], 
                                  e_par[len(vel_files)*2 +7*npl + 2 + rv_gp_npar + 3*npl + N_transit_files*4 + tra_gp_npar + tr_model[3][i] + 3+ npl]]



    if obj.type_fit["AST"] and obj.use_ast_hipp_gaia == True: 

        obj.ast_alpha_err[0]    =    e_par[-6]
        obj.ast_delta_err[0]    =    e_par[-5]
        obj.ast_pi_err[0]       =    e_par[-4]       
        obj.ast_mu_alpha_err[0] =    e_par[-3]       
        obj.ast_mu_delta_err[0] =    e_par[-2]

   #####################

    if len(flags) != 0:
        print("Best lnL: %s"%obj.loglik)
        print("Best fit par.:")

        for j in range(len(pp)):
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
    file_ses = open("%s.ses"%target_name, 'wb')
    dill.dump(obj, file_ses)
    file_ses.close()

    #if sys.version_info[0] == 2:
   #     os.system("python2 ./lib/run_ns_from_ses.py -ses ./%s.ses %s"%(target_name,target_name))
   # elif sys.version_info[0] == 3:
    #    os.system("python3 ./lib/run_ns_from_ses.py -ses ./%s.ses %s"%(target_name,target_name))

    os.system("python%s.%s ./lib/run_ns_from_ses.py -ses ./%s.ses %s"%(
        sys.version_info[0],sys.version_info[1],target_name,target_name))

    file_ses2 = open("%s_out.ses"%target_name, 'rb')
    obj = dill.load(file_ses2)
    file_ses2.close()

    print("--- %s seconds ---" % (time.time() - start_time))
    os.system("rm %s.ses"%target_name)
    os.system("rm %s_out.ses"%target_name)

    return obj



def run_nestsamp(obj, **kwargs):

    '''Performs nested sampling and saves results. Work in progress.....'''

    start_time = time.time()
    rtg = obj.rtg
    check_temp_RV_file(obj)


#    vel_files = []
#    for i in range(obj.ndset):
#        # path for each dataset
#        vel_files.append(obj.filelist.files[i].path)

    vel_files = [0]*obj.ndset

    N_transit_files = len([x for x in range(60) if len(obj.tra_data_sets[x]) != 0])

    tr_files = obj.tra_data_sets
    tr_mo    = obj.ld_m
    tr_ld    = obj.ld_u
    tr_gr    = obj.ld_gr
    tr_gr_ind= obj.ld_gr_ind
    
    tr_model = np.array([tr_mo,tr_ld,tr_gr,tr_gr_ind], dtype=object)
    tr_params = obj.tr_params


    final_array = []
    jd, rvs, sig, ids = [], [], [], []


    for z in range(len(obj.rv_data_sets)):
        if len(obj.rv_data_sets[z])==0:
            continue
 
        jd = np.append(jd, obj.rv_data_sets[z][0], axis = 0)
        rvs = np.append(rvs, obj.rv_data_sets[z][1], axis = 0)
        sig = np.append(sig, obj.rv_data_sets[z][2], axis = 0)
        ids = np.append(ids, obj.rv_data_sets[z][3]+1, axis = 0)
    final_array = np.array([jd, rvs, sig, ids]).T



    rvs_files  = final_array
    ttv_files  = obj.ttv_data_sets
    ast_files  = obj.ast_data_sets
    ast_files2 = obj.ast_data_sets_hipp_gaia
    
    
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


    obj.prepare_for_mcmc(rtg = rtg)
    pp = obj.par_for_mcmc
    ee = obj.e_for_mcmc
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

#    re = Rvfit()
    opt = {"eps":obj.dyn_model_accuracy*1e-13,
           "dt":obj.time_step_model*86400.0,
           "when_to_kill":when_to_kill,
           "copl_incl":obj.copl_incl,
           "jit_flag":obj.jit_flag,
           "hkl":obj.hkl,
           "rho":obj.rho,                      
           "cwd":obj.cwd, 
           "gr_flag":obj.gr_flag,
           "ns_samp_method":obj.ns_samp_method,
           "TTV":obj.type_fit["TTV"],
           "AST":obj.type_fit["AST"],
           "RVS_files":rvs_files,
           "TTV_files":ttv_files,
           "AST_files":ast_files,
           "AST_files_hipp":ast_files2,                       
           "TTV_times":obj.ttv_times,
           "AST_times":obj.ast_times,
           "AMD_stab":obj.NS_AMD_stab, 
           "Nbody_stab":obj.NS_Nbody_stab,
           "get_TTVs":obj.get_TTVs,
           "link_RV_GP":obj.link_RV_GP,
           "tra_GP_kernel":obj.tra_gp_kernel,
           "RV_GP_kernel":obj.gp_kernel,
           "nest_weighted":obj.nest_weighted}
#            "re":re}

    gps = []
    if (rtg[1]):
        initiate_RV_gps(obj)
        gps = obj.gps

    tra_gps = []
    if (rtg[3]) and N_transit_files != 0:
        initiate_tansit_gps(obj)
        tra_gps = obj.tra_gps

    ndim, nwalkers = len(pp), len(pp)*obj.live_points_fact

    ################## prior TESTS ########################

    def prior_transform(p):

        u_trans = np.zeros(len(p))
        for j in range(len(p)):

            
            #u_trans[j] = trans_uni(p[j],bb[j][0],bb[j][1])
            
            if priors[0][j,2] == True:
               # u_trans[j] = trans_norm(p[j],priors[0][j,0],priors[0][j,1])
                u_trans[j] = truncated_trans_norm(p[j],priors[0][j,0],priors[0][j,1],bb[j][0],bb[j][1])                
                
            elif priors[1][j,2] == True:
                #u_trans[j] = trans_loguni(p[j],priors[1][j,0],priors[1][j,1])
                u_trans[j] = truncated_trans_loguni(p[j],priors[1][j,0],priors[1][j,1],bb[j][0],bb[j][1])                   
                
            else:
                u_trans[j] = trans_uni(p[j],bb[j][0],bb[j][1])

            #if p[j] <= bb[j][0] or p[j] >= bb[j][1]:
            #    u_trans[j] = 0  
 
        return u_trans


    def truncated_trans_norm(p, loc, scale, lower, upper):
        """
        Computes the inverse CDF (PPF) of a normal distribution truncated within [lower, upper].
        
        Parameters:
            p (float): Probability value, can be any number.
            loc (float): Mean of the normal distribution.
            scale (float): Standard deviation of the normal distribution.
            lower (float): Lower bound of truncation.
            upper (float): Upper bound of truncation.
            
        Returns:
            float: Truncated value of the inverse CDF.
        """
        # Calculate the CDF of the bounds
        cdf_lower = stats.norm.cdf(lower, loc=loc, scale=scale)
        cdf_upper = stats.norm.cdf(upper, loc=loc, scale=scale)

        # Rescale p to the truncated range, clamping if p is out of bounds
        p_clamped = max(0.0, min(1.0, p))  # Ensure p is within [0, 1]
        p_truncated = cdf_lower + p_clamped * (cdf_upper - cdf_lower)

        # Compute the inverse CDF (PPF) for the truncated range
        return stats.norm.ppf(p_truncated, loc=loc, scale=scale)

    def truncated_trans_loguni(p, a, b, flat_lower, flat_upper):
        """
        Transforms a value `p` to a log-uniform distribution within bounds [a, b],
        with additional constraints from a flat prior [flat_lower, flat_upper].
        
        Parameters:
            p (float): Input value, can be any number.
            a (float): Lower bound of the log-uniform distribution (must be positive).
            b (float): Upper bound of the log-uniform distribution (must be positive).
            flat_lower (float): Lower bound of the flat prior.
            flat_upper (float): Upper bound of the flat prior.
            
        Returns:
            float: Transformed value within the combined bounds.
        """
        # Ensure bounds for log-uniform are positive
        if a <= 0 or b <= 0:
            return (flat_lower + flat_upper) / 2.0  # Default to midpoint of flat prior if log bounds are invalid

        # Clamp p to [0, 1]
        p_clamped = max(0.0, min(1.0, p))
        
        # Transform to log-uniform within [a, b]
        log_uni_result = np.exp(np.log(a) + p_clamped * (np.log(b) - np.log(a)))

        # Clamp to flat prior bounds [flat_lower, flat_upper]
        return max(flat_lower, min(flat_upper, log_uni_result))

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
    from pathos.multiprocessing import ProcessingPool as Pool
    #from pathos.pools import ProcessPool as Pool
#    from multiprocessing import Pool
#    from contextlib import closing

    #import multiprocessing as mp
    #from multiprocessing import Pool
   # from pathos.pools import ThreadPool as Pool

    dynesty_samp = obj.ns_samp_method
    print_progress = obj.ns_progress #std_output
    threads = int(obj.ns_threads)
    stop_crit = obj.stop_crit
    Dynamic_nest = obj.Dynamic_nest
    ns_bound = obj.ns_samp_bound
    ns_pfrac = obj.ns_pfrac
    ns_use_stop = obj.ns_use_stop 
    


    if obj.ns_maxiter[0] == False:
         ns_maxiter = None
    else:
         ns_maxiter = obj.ns_maxiter[1]

    if obj.ns_maxcall[0] == False:
         ns_maxcall = None
    else:
         ns_maxcall = obj.ns_maxcall[1]
 
    #from multiprocessing import get_context
    #ctx = get_context("fork")

    pool = Pool(ncpus=threads)#, context=ctx)
    #thread = Pool(threads)#, context=ctx)

    #pool = Pool(nodes=threads)

    if Dynamic_nest == False:
        print("'Static' Nest. Samp. is running, please wait...")

        if threads > 1:
 
            sampler = dynesty.NestedSampler(partial_func, prior_transform, ndim, nlive=nwalkers, pool = pool,
                                                queue_size=threads, sample = dynesty_samp, bound = ns_bound)
            sampler.run_nested(print_progress=print_progress,dlogz=stop_crit, 
            maxiter = ns_maxiter, maxcall = ns_maxcall ) #dlogz=stop_crit,
            #pool.close()
            #pool.join()
            #pool.clear()

        else:
             sampler = dynesty.NestedSampler(partial_func, prior_transform, ndim, nlive=nwalkers, sample = dynesty_samp, bound = ns_bound)
             sampler.run_nested(print_progress=print_progress,dlogz=stop_crit, 
             maxiter = ns_maxiter, maxcall = ns_maxcall )


        #obj.dyn_res = sampler.results
        sampler.results.summary()

    else:
        print("'Dynamic' Nest. Samp. is running, please wait...")

        if threads > 1:
 
 
            if obj.new_mp == True:
    
                import contextlib
                print("This is a test!!!!")
                with contextlib.closing(pool) as pool:
                    
                    sampler = dynesty.DynamicNestedSampler(partial_func, prior_transform, ndim, pool = pool,
                                                               queue_size=threads, sample = dynesty_samp, bound = ns_bound) # nlive=nwalkers,
            
                    sampler.run_nested(print_progress=print_progress,dlogz_init=stop_crit,nlive_init=nwalkers, 
                    maxiter = ns_maxiter, maxcall = ns_maxcall,use_stop = ns_use_stop, wt_kwargs={'pfrac': ns_pfrac})   #nlive_batch=1
                    #pool.close()
                    #pool.join()
                    #pool.clear()
                                
            else:
                
                sampler = dynesty.DynamicNestedSampler(partial_func, prior_transform, ndim, pool = pool,
                                                       queue_size=threads, sample = dynesty_samp, bound = ns_bound) # nlive=nwalkers,
    
                sampler.run_nested(print_progress=print_progress,dlogz_init=stop_crit,nlive_init=nwalkers, 
                maxiter = ns_maxiter, maxcall = ns_maxcall,use_stop = ns_use_stop, wt_kwargs={'pfrac': ns_pfrac})   #nlive_batch=1
                #pool.close()
                #pool.join()
                #pool.clear()
               
            

        else:
            sampler = dynesty.DynamicNestedSampler(partial_func, prior_transform, ndim, sample = dynesty_samp, bound = ns_bound)
            sampler.run_nested(print_progress=print_progress,dlogz_init=stop_crit,nlive_init=nwalkers,  
            maxiter = ns_maxiter, maxcall = ns_maxcall,use_stop = ns_use_stop, wt_kwargs={'pfrac': ns_pfrac} ) 

        # just in case
        pool.close()
        pool.join()
        pool.clear()

       # obj.dyn_res = sampler.results

        res = ("niter: {:d}\n"
               "ncall: {:d}\n"
               "eff(%): {:6.3f}\n"
               "logz: {:6.3f} +/- {:6.3f}".format(sampler.results.niter, sum(sampler.results.ncall),
                       sampler.results.eff, sampler.results.logz[-1], sampler.results.logzerr[-1]))

        print('Summary\n=======\n'+res)

 

   # print("--- %s seconds ---" % (time.time() - start_time))




    ln = np.hstack(sampler.results.logl)
    add_ns_samples(obj,sampler)

 
    if opt["nest_weighted"] == True:
        weighted = np.exp(sampler.results.logwt - sampler.results.logz[-1])
        obj.ns_sampler.samples  =  dill.copy(dynesty.utils.resample_equal(sampler.results.samples, weighted))   
  

    if obj.ns_fileoutput == True:
       # start_time = time.time()
       # print("Please wait... writing the ascii file")
        dirname, basename = os.path.split(obj.nest_sample_file)
        if os.path.exists(dirname):
            ns_file = obj.nest_sample_file
        else:
            print("%s does not exist! Sample file will be saved in %s/ns_sample_file"%(dirname,obj.cwd))
            ns_file = "%s/ns_sample_file"%obj.cwd

        samples = np.array(obj.ns_sampler.samples)
        outfile = open(str(ns_file), 'w') # file to save samples
        for j in range(len(samples)):
            outfile.write("%s  " %(ln[j]))
            for z in range(len(pp)):
                outfile.write("%s  " %(samples[j,z]))
            outfile.write("\n")
        outfile.close()


    obj.nest_stat["mean"] = get_mean_of_samples(obj.ns_sampler.samples,len(pp))
    obj.nest_stat["median"] = get_median_of_samples(obj.ns_sampler.samples,len(pp))
    samp_maxlnl, maxlnl = get_best_lnl_of_samples(obj.ns_sampler.samples,ln, len(pp))
    obj.nest_stat["best"] = samp_maxlnl
    obj.nest_stat["mode"] = get_mode_of_samples(obj.ns_sampler.samples,len(pp))
    obj.nest_stat["MAD"]  = get_MAD_of_samples(obj.ns_sampler.samples,len(pp))


    if (obj.ns_save_means):
        obj.par_for_mcmc = obj.nest_stat["mean"]
        pp = obj.nest_stat["mean"]

    elif (obj.ns_save_median):
        obj.par_for_mcmc = obj.nest_stat["median"]
        pp =  obj.nest_stat["median"]

    elif (obj.ns_save_maxlnL):
        obj.par_for_mcmc = obj.nest_stat["best"]
        pp =  obj.nest_stat["best"]

    elif (obj.ns_save_mode):
        obj.par_for_mcmc = obj.nest_stat["mode"]
        pp =  obj.nest_stat["mode"]
    # else:
   #     pp = obj.par_for_mcmc

    #sampler.samples = obj.ns_sampler.samples


    if obj.nest_mad == False:
        new_par_errors = [[float(obj.par_for_mcmc[i] - np.percentile(obj.ns_sampler.samples[:,i], [level])), float(np.percentile(obj.ns_sampler.samples[:,i], [100.0-level])-obj.par_for_mcmc[i])] for i in range(len(obj.par_for_mcmc))]
    else:
        new_par_errors = [[float(obj.nest_stat["MAD"][i]),float(obj.nest_stat["MAD"][i])] for i in range(len(obj.par_for_mcmc))]


#    newparams = obj.generate_newparams_for_mcmc(obj.par_for_mcmc)


    #obj.fitting(minimize_loglik=True, amoeba_starts=0, npoints=obj.model_npoints, outputfiles=[1,1,1]) # this will help update some things

    obj = return_results(obj, pp, ee, par, flags, npl,vel_files, tr_files,  tr_model, tr_params, epoch, stmass, bb, priors, gps, tra_gps, rtg, mix_fit, new_par_errors,mod,opt)

#    obj.update_with_mcmc_errors(new_par_errors)

#    obj.overwrite_params(newparams)

#    obj.hack_around_rv_params()


    if (obj.ns_save_means):
        obj.loglik = maxlnl

    elif (obj.ns_save_maxlnL):
        obj.loglik = maxlnl


    bestfit_labels      = ["median","mean","mode","best_samp","best_gui","none","mass","semimajor","radius"]
    bestfit_labels_bool = [obj.ns_save_median,obj.ns_save_means,obj.ns_save_mode, obj.ns_save_maxlnL,False,False,False,False,False]

 
 
    if obj.ns_save_sampler != True:
        del obj.ns_sampler        
        obj.ns_sampler = []

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
    file_ses = open("%s.ses"%target_name, 'wb')
    dill.dump(obj, file_ses)
    file_ses.close()

   # if sys.version_info[0] == 2:
#        os.system("python2 ./lib/run_mcmc_from_ses.py -ses ./%s.ses %s"%(target_name,target_name))
#    elif sys.version_info[0] == 3:
  #      os.system("python3 ./lib/run_mcmc_from_ses.py -ses ./%s.ses %s"%(target_name,target_name))

    #if sys.version_info[0] == 2:
    os.system("python%s.%s ./lib/run_mcmc_from_ses.py -ses ./%s.ses %s"%(
        sys.version_info[0],sys.version_info[1],target_name,target_name))
 
    file_ses2 = open("%s_out.ses"%target_name, 'rb')
    obj = dill.load(file_ses2)
    file_ses2.close()

    print("--- %s seconds ---" % (time.time() - start_time))
    os.system("rm %s.ses"%target_name)
    os.system("rm %s_out.ses"%target_name)

    return obj

def run_mcmc(obj, **kwargs):

    '''Performs MCMC and saves results'''

    start_time = time.time()

    rtg = obj.rtg
    check_temp_RV_file(obj)

#    vel_files = []
#    for i in range(obj.ndset):
#        vel_files.append(obj.filelist.files[i].path)

    vel_files = [0]*obj.ndset

    N_transit_files = len([x for x in range(60) if len(obj.tra_data_sets[x]) != 0])

    tr_files = obj.tra_data_sets
    tr_mo    = obj.ld_m
    tr_ld    = obj.ld_u
    tr_gr    = obj.ld_gr
    tr_gr_ind= obj.ld_gr_ind
    
    tr_model = np.array([tr_mo,tr_ld,tr_gr,tr_gr_ind], dtype=object)
    tr_params = obj.tr_params

    final_array = []
    jd, rvs, sig, ids = [], [], [], []


    for z in range(len(obj.rv_data_sets)):
        if len(obj.rv_data_sets[z])==0:
            continue
 
        jd = np.append(jd, obj.rv_data_sets[z][0], axis = 0)
        rvs = np.append(rvs, obj.rv_data_sets[z][1], axis = 0)
        sig = np.append(sig, obj.rv_data_sets[z][2], axis = 0)
        ids = np.append(ids, obj.rv_data_sets[z][3]+1, axis = 0)
    final_array = np.array([jd, rvs, sig, ids]).T


    rvs_files  = final_array #np.array([obj.fit_results.jd,obj.fit_results.rvs,obj.fit_results.rv_err,obj.fit_results.idset +1]).T
    ttv_files  = obj.ttv_data_sets
    ast_files  = obj.ast_data_sets
    ast_files2 = obj.ast_data_sets_hipp_gaia
    
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

    obj.prepare_for_mcmc(rtg = rtg)
    pp = obj.par_for_mcmc 
    ee = obj.e_for_mcmc 
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
    level = (100.0- obj.percentile_level)/2.0


#    re = Rvfit()
    opt = {"eps":obj.dyn_model_accuracy*1e-13,
           "dt":obj.time_step_model*86400.0,
           "when_to_kill":when_to_kill,
           "copl_incl":obj.copl_incl,
           "jit_flag":obj.jit_flag,
           "hkl":obj.hkl,
           "rho":obj.rho,           
           "cwd":obj.cwd,
           "gr_flag":obj.gr_flag,
           "TTV":obj.type_fit["TTV"],
           "AST":obj.type_fit["AST"],
           "RVS_files":rvs_files,
           "TTV_files":ttv_files,
           "AST_files":ast_files,
           "AST_files_hipp":ast_files2,                      
           "TTV_times":obj.ttv_times,
           "AST_times":obj.ast_times,
           "AMD_stab":obj.mcmc_AMD_stab, 
           "Nbody_stab":obj.mcmc_Nbody_stab,
           "get_TTVs":obj.get_TTVs,
           "link_RV_GP":obj.link_RV_GP,
           "tra_GP_kernel":obj.tra_gp_kernel,
           "RV_GP_kernel":obj.gp_kernel}
#            "re":re}

    gps = []
    if (rtg[1]):
        initiate_RV_gps(obj)
        gps = obj.gps
#        rv_gp_npar = len(gps.get_parameter_vector())
#    else:
#        rv_gp_npar = 0
        
    tra_gps = []
    if (rtg[3]) and N_transit_files != 0:
        initiate_tansit_gps(obj)
        tra_gps = obj.tra_gps

    from pathos.multiprocessing import ProcessingPool as Pool
    #from pathos.pools import ProcessPool as Pool
    #from pathos.threading import ThreadPool as Pool
    pool=Pool(ncpus=obj.mcmc_threads)
    #pool=Pool(processes=threads-1)
    #from pathos.threading import ThreadPool
    #import mkl
   # mkl.set_num_threads(1)
    #from multiprocessing import Pool
    #pool=Pool(threads)

    #from schwimmbad import MPIPool

   # with MPIPool() as pool:

  #      if not pool.is_master():
   #         pool.wait()
   #         sys.exit(0)



    ndim, nwalkers = len(pp), len(pp)*obj.nwalkers_fact

    pos = [pp + obj.gaussian_ball*np.random.rand(ndim) for i in range(nwalkers)]
    
    sampler = CustomSampler(nwalkers, ndim, lnprob_new, args=(mod, par, flags, npl, vel_files, tr_files,  tr_model, tr_params, epoch, stmass, bb, priors, gps, tra_gps, rtg, mix_fit,opt), pool=pool)

    #sampler = CustomSampler(nwalkers, ndim, lnprob_new, args=(mod, par, flags, npl, vel_files, tr_files,  tr_model, tr_params, epoch, stmass, bb, priors, gps, tra_gps, rtg, mix_fit,opt), threads = threads)

    emcee_version = emcee.__version__
    
    if int(emcee_version[0]) == 3:
        # burning phase
        if obj.mcmc_burning_ph > 0:
            pos, prob, state  = sampler.run_mcmc(pos,int(obj.mcmc_burning_ph), progress= obj.mcmc_progress)
            sampler.reset()
        # now perform the MCMC
        pos, prob, state  = sampler.run_mcmc(pos,int(obj.mcmc_ph), progress= obj.mcmc_progress)
        
    else:
        print("Please upgrade 'emcee' to Ver. 3! E.g. 'sudo pip install emcee -U'")
        # burning phase
        if obj.mcmc_burning_ph > 0:
            pos, prob, state  = sampler.run_mcmc(pos,int(obj.mcmc_burning_ph))
            sampler.reset()
        # now perform the MCMC
        pos, prob, state  = sampler.run_mcmc(pos,int(obj.mcmc_ph))
        
        
    #ln = np.hstack(sampler.lnprobability)
   # sampler.save_samples(obj.f_for_mcmc,obj.ndset,obj.npl)
    sampler.save_samples(obj.f_for_mcmc,obj.ndset,obj.npl,obj.hkl)

    pool.close()
    pool.join()
    pool.clear()

 #  print("--- %s seconds ---" % (time.time() - start_time))


    fileoutput = obj.mcmc_fileoutput
    if (fileoutput):
     #   start_time = time.time()
    #    print("Please wait... writing the ascii file")
        dirname, basename = os.path.split(obj.mcmc_sample_file)
        if os.path.exists(dirname):
            mcmc_file = obj.mcmc_sample_file
        else:
            print("%s does not exist! Sample file will be saved in %s/mcmc_sample_file"%(dirname,obj.cwd))
            mcmc_file = "%s/mcmc_sample_file"%obj.cwd

        outfile = open(str(mcmc_file), 'w') # file to save samples
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
    obj.mcmc_stat["median"] = sampler.median
    obj.mcmc_stat["best"] = sampler.maxlnL
    obj.mcmc_stat["mode"] = get_mode_of_samples(sampler.samples,len(pp))
    obj.mcmc_stat["MAD"]  = get_MAD_of_samples(sampler.samples,len(pp))


    # Now we will save new parameters and their errors (different + and - errors in this case). Flag save_means determines if we want to take means as new best fit parameters or stick to old ones and calculate errors with respect to that
    if (obj.mcmc_save_means):
        obj.par_for_mcmc = obj.mcmc_stat["mean"] # we will not need to keep the old parameters in this attribbute, so let's store the means now
        pp = obj.mcmc_stat["mean"]
        
    elif (obj.mcmc_save_median):
        obj.par_for_mcmc = obj.mcmc_stat["median"] #
        pp =  obj.mcmc_stat["median"]

    elif (obj.mcmc_save_maxlnL):
        obj.par_for_mcmc = obj.mcmc_stat["best"] # we will not need to keep the old parameters in this attribbute, so let's store the means now
        pp =  obj.mcmc_stat["best"]

    elif (obj.mcmc_save_mode):
        obj.par_for_mcmc = obj.mcmc_stat["mode"] # we will not need to keep the old parameters in this attribbute, so let's store the means now
        pp =  obj.mcmc_stat["mode"]
    # else:
   #     pp = obj.par_for_mcmc

    if obj.mcmc_mad == False:
        new_par_errors = [[float(obj.par_for_mcmc[i] - np.percentile(sampler.samples[:,i], [level])),float(np.percentile(sampler.samples[:,i], [100.0-level])-obj.par_for_mcmc[i])] for i in range(len(obj.par_for_mcmc))]
    else:
        new_par_errors = [[float(obj.mcmc_stat["MAD"][i]),float(obj.mcmc_stat["MAD"][i])] for i in range(len(obj.par_for_mcmc))]

   # newparams = obj.generate_newparams_for_mcmc(obj.par_for_mcmc)

    #obj.fitting(minimize_loglik=True, amoeba_starts=0, npoints=obj.model_npoints, outputfiles=[1,1,1]) # this will help update some things


#    if obj.copl_incl == True:
#        incl_c = par[len(vel_files)*2 +7*0+5]
#        for i in range(npl):
#            par[len(vel_files)*2 +7*i+5] = incl_c 
#            obj.i[i] = incl_c

    obj = return_results(obj, pp, ee, par, flags, npl,vel_files, tr_files, tr_model, tr_params, epoch, stmass, bb, priors, gps,tra_gps, rtg, mix_fit, new_par_errors,mod,opt)

 
    if (obj.mcmc_save_means):
        obj.loglik = sampler.lnL_max

    elif (obj.mcmc_save_maxlnL):
        obj.loglik = sampler.lnL_max

 

    if(obj.mcmc_save_sampler):     
        add_mcmc_samples(obj,sampler)       
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
    (from emcee?)
    """
    def __init__(self, f, args, kwargs=None):
        self.f = f
        self.args = [] if args is None else args
        self.kwargs = {} if kwargs is None else kwargs

    def __call__(self, x):
        try:
            result = self.f(x, *self.args, **self.kwargs)
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


    def __init__(self, inputfile='init.init', name='session'):

        #### Old stuff; some these must be replaced! ####
        self.rvmod = dill.copy(Rvfit())

        # saving the name for the inputfile and the information that it has not yet been processed
        self.name = name # for example the name of a host star of the planetary system
        self.input_processed = False

        self.stellar_mass_known=False
        self.filelist=rvfile_list(0,[],[])
        self.mod_dynamical=False
        self.epoch=0.0
        self.npl=0
        self.use=use_flags([False]*60,[False]*60,[False]*70,False,False)
        self.params=PlanetParamsWrap([0.0]*60,[0.0]*60,[0.0]*70,0.0,1.0)
        self.param_errors=PlanetParamsErrorsWrap(self.rvmod.offset_errors, self.rvmod.jitter_errors, self.rvmod.planet_params_errors,
                                                 [0.0, 0.0], 0.0)
        self.bounds = parameter_bounds([0.0,0.0]*60,[0.0,0.0]*60,[0.0,0.0]*70,[0.0,0.0],[0.0,0.0]*4,[0.0,0.0])
        
        self.use_planet = [0,0,0,0,0,0,0,0,0]

        self.ses_notes = ''

        self.ndset = 0

        #del re

        ########## new stuff ##########
        self.init_pl_params()
        self.init_mcmc_par()
        self.init_nest_par()

        self.fit_performed = False
        self.model_saved=False
        self.stat_saved=False
 
        self.f_for_mcmc=[]
        self.par_for_mcmc=[]
        self.e_for_mcmc=[]
        self.b_for_mcmc=[]

        self.init_St_params()
        self.init_st_mass()
        self.init_mass_a()
        
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
        self.init_tra_dilution()
        self.init_tra_lintr()
        self.init_tra_quadtr()
        
        self.cwd = os.getcwd()

        self.init_pl_arb()
        self.init_orb_evol_arb()

        self.type_fit = {"RV": True,"Transit": False,"TTV":False , "AST":False}

        self.gr_flag = False
        self.hkl = False
        self.rho = False        
        self.copl_incl = False
        self.jit_flag = False
        self.rtg = [True,False,False,False]
        self.link_RV_GP = [False,False,False,False,False,False]                     
                           
        self.ttv_data_sets = {k: [] for k in range(10)}
        self.ast_data_sets = {k: [] for k in range(10)}
        self.act_data_sets = {k: [] for k in range(60)}
        self.tra_data_sets = {k: [] for k in range(60)}
        self.rv_data_sets  = {k: [] for k in range(60)}
        self.ast_data_sets_hipp_gaia = {k: [] for k in range(10)}               
        self.act_data_sets_init = {k: [] for k in range(60)}
        self.tra_data_sets_init = {k: [] for k in range(60)}
        self.rv_data_sets_init  = {k: [] for k in range(60)}
        
        # in this case we need to create a new kernel, but we will only give it this information which is needed for plotting

        self.fit_results = dill.copy(Rvfit())
        self.loglik=0.0
        self.rms=0.0
        self.chi2=0.0
        self.reduced_chi2=0.0
        self.sampler=None
        self.sampler_saved=False

        self.init_auto_fit()


        self.model_npoints = 2000
        self.model_max = 5
        self.model_min =0


        self.init_orb_evol()

        self.tls = []
        self.tls_o_c = []

        self.gls = []
        self.gls_o_c =[]

        self.mlp = []

        self.optim_AMD_stab   = False
        self.optim_Nbody_stab = False
        self.init_dynfit_settings()


        self.ph_data = {k: [] for k in range(9)}
        self.ph_model = {k: [] for k in range(9)}

        self.ph_data_tra = {k: [] for k in range(9)}
        self.ph_model_tra = {k: [] for k in range(9)}

        self.parameters = []

        self.pyqt_symbols_rvs = {k: 'o' for k in range(60)} # ['o','t','t1','t2','t3','s','p','h','star','+','d']
        self.pyqt_symbols_act = {k: 'o' for k in range(60)} # ['o','t','t1','t2','t3','s','p','h','star','+','d']
        self.pyqt_symbols_tra = {k: 'o' for k in range(60)} # ['o','t','t1','t2','t3','s','p','h','star','+','d']
        self.pyqt_symbols_ttv = {k: 'o' for k in range(60)} # ['o','t','t1','t2','t3','s','p','h','star','+','d']
        self.pyqt_symbols_ast = {k: 'o' for k in range(60)} # ['o','t','t1','t2','t3','s','p','h','star','+','d']

        self.pyqt_symbols_size_rvs = {k: 6 for k in range(60)} #[6,6,6,6,6,6,6,6,6,6] #
        self.pyqt_symbols_size_act = {k: 4 for k in range(60)} #[4,4,4,4,4,4,4,4,4,4] #
        self.pyqt_symbols_size_tra = {k: 2 for k in range(60)} #[2,2,2,2,2,2,2,2,2,2] #
        self.pyqt_symbols_size_ttv = {k: 4 for k in range(60)} #[2,2,2,2,2,2,2,2,2,2] #
        self.pyqt_symbols_size_ast = {k: 4 for k in range(60)} #[2,2,2,2,2,2,2,2,2,2] #

        self.pyqt_color_alpha_rvs = {k: 255 for k in range(60)} #[6,6,6,6,6,6,6,6,6,6] #
        self.pyqt_color_alpha_act = {k: 255 for k in range(60)} #[4,4,4,4,4,4,4,4,4,4] #
        self.pyqt_color_alpha_tra = {k: 255 for k in range(60)} #[2,2,2,2,2,2,2,2,2,2] #
        self.pyqt_color_alpha_ttv = {k: 255 for k in range(60)} #[2,2,2,2,2,2,2,2,2,2] #        
        self.pyqt_color_alpha_ast = {k: 255 for k in range(60)} #[2,2,2,2,2,2,2,2,2,2] #        


        self.colors = ['#0066ff',  '#ff0000','#00aa00','#00ffff','#cc33ff','#ff9900','#cccc00','#3399ff','#990033','#339933','#808080']

        self.act_colors = ['#0066ff',  '#ff0000','#00aa00','#00ffff','#cc33ff','#ff9900','#cccc00','#3399ff','#990033','#339933','#0066ff', '#ff0000','#00aa00','#00ffff','#cc33ff','#ff9900','#cccc00','#3399ff','#990033','#339933','#000000']
        self.tra_colors = ['#0066ff',  '#ff0000','#00aa00','#00ffff','#cc33ff','#ff9900','#cccc00','#3399ff','#990033','#339933','#0066ff', '#ff0000','#00aa00','#00ffff','#cc33ff','#ff9900','#cccc00','#3399ff','#990033','#339933','#000000','#0066ff',  '#ff0000','#00aa00','#00ffff','#cc33ff','#ff9900','#cccc00','#3399ff','#990033','#339933','#0066ff', '#ff0000','#00aa00','#00ffff','#cc33ff','#ff9900','#cccc00','#3399ff','#990033','#339933','#000000','#0066ff',  '#ff0000','#00aa00','#00ffff','#cc33ff','#ff9900','#cccc00','#3399ff','#990033','#339933','#0066ff', '#ff0000','#00aa00','#00ffff','#cc33ff','#ff9900','#cccc00','#3399ff','#990033','#339933','#000000']
        self.rvs_colors = ['#0066ff',  '#ff0000','#00aa00','#00ffff','#cc33ff','#ff9900','#cccc00','#3399ff','#990033','#339933','#666699', '#ff0000','#00aa00','#00ffff','#cc33ff','#ff9900','#cccc00','#3399ff','#990033','#339933','#000000']

        self.gls_colors = ['#ff0000',  '#ff0000','#ff0000','#ff0000','#ff0000','#ff0000','#ff0000','#ff0000','#ff0000','#ff0000','#000000']
        self.ttv_colors = ['#0066ff',  '#ff0000','#00aa00','#00ffff','#cc33ff','#ff9900','#cccc00','#3399ff','#990033','#339933','#000000']
        self.ast_colors = ['#0066ff',  '#ff0000','#00aa00','#00ffff','#cc33ff','#ff9900','#cccc00','#3399ff','#990033','#339933','#000000']

        self.init_sciPy_minimizer()
        self.init_ld_model()
        
        self.init_TTVs()
        self.get_TTVs = [False,{k: [] for k in range(10)}]

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



        self.new_mp = False

#######################################
 

    def init_ld_model(self):


        self.ld_models = ["uniform", "linear", "quadratic", "nonlinear"]
        self.ld_m = ["quadratic"]*60    #limb darkening model

        self.ld_u = {k: [0.12, 0.35 ] for k in range(60)}

        self.ld_u_lin    = {k: [0.35] for k in range(60)}
        self.ld_u_quad = {k: [0.12, 0.35 ] for k in range(60)}
        self.ld_u_nonlin = {k: [0.55,0.12, 0.35,-0.11] for k in range(60)}

        self.ld_u_lin_use    = {k: [False] for k in range(60)}
        self.ld_u_quad_use   = {k: [False, False] for k in range(60)}
        self.ld_u_nonlin_use = {k: [False, False,False, False] for k in range(60)}

        self.ld_u_lin_err    = {k: [[0.0,0.0]] for k in range(60)}
        self.ld_u_quad_err   = {k: [[0.0,0.0], [0.0,0.0]] for k in range(60)}
        self.ld_u_nonlin_err = {k: [[0.0,0.0], [0.0,0.0],[0.0,0.0], [0.0,0.0]] for k in range(60)}

        self.ld_u_lin_bound       = {k: np.array([[-1.0,1.0]]) for k in range(60)}
        self.ld_u_quad_bound      = {k: np.array([[-1.0,1.0],[-1.0,1.0]]) for k in range(60)}
        self.ld_u_nonlin_bound    = {k: np.array([[-1.0,1.0],[-1.0,1.0],[-1.0,1.0],[-1.0,1.0]]) for k in range(60)}

        self.ld_u_lin_norm_pr     = {k: np.array([[0.1,0.05, False]]) for k in range(60)}
        self.ld_u_quad_norm_pr    = {k: np.array([[0.0,1.0, False],[0.0,1.0, False]]) for k in range(60)}
        self.ld_u_nonlin_norm_pr  = {k: np.array([[0.0,1.0, False],[0.0,1.0, False],[0.0,1.0, False],[0.0,1.0, False]]) for k in range(60)}

        self.ld_u_lin_jeff_pr     = {k: np.array([[0.1,0.05, False]]) for k in range(60)}
        self.ld_u_quad_jeff_pr    = {k: np.array([[0.0,1.0, False],[0.0,1.0, False]]) for k in range(60)}
        self.ld_u_nonlin_jeff_pr  = {k: np.array([[0.0,1.0, False],[0.0,1.0, False],[0.0,1.0, False],[0.0,1.0, False]]) for k in range(60)}

        self.ld_u_lin_str         = {k: [r'ld-quad-1$_%s$'%str(k+1)] for k in range(60)}
        self.ld_u_quad_str        = {k: [r'ld-quad-1$_%s$'%str(k+1),r'ld-quad-2$_%s$'%str(k+1)] for k in range(60)}
        self.ld_u_nonlin_str      = {k: [r'ld-quad-1$_%s$'%str(k+1),r'ld-quad-2$_%s$'%str(k+1),r'ld-quad-3$_%s$'%str(k+1),r'ld-quad-4$_%s$'%str(k+1)] for k in range(60)}

        self.ld_gr     = list(range(60))
        self.ld_gr_ind = [0]* len(self.ld_gr)

        ############################################


#    def constants(self):
#        self.AU


    def init_auto_fit(self):
        self.auto_fit_max_pl = 2
        self.auto_fit_allow_ecc = True
        self.auto_fit_FAP_level = 0.001

    def init_xyz(self):
        self.xyz_mass = {k: [] for k in range(10)}
        self.xzy    = {k: [] for k in range(10)}
        self.uvw    = {k: [] for k in range(10)}
        self.rpl    = {k: [] for k in range(10)}
        self.rhill  = {k: [] for k in range(10)}


    def init_mass_a(self):
        self.masses    =[]
        self.semimajor = []
        
        for i in range(9):
            masses,semimajor  = mass_a_from_Kepler_fit([self.K[i], self.P[i], self.e[i], self.w[i], self.M0[i]],1,self.st_mass[0])
            self.masses.append(masses[0]),self.semimajor.append(semimajor[0])
            
    def init_pl_params(self):

        #### RV #####
        self.K    = {k: 60.0 for k in range(9)}
        self.P    = {k: 100.0 + 150.0*k for k in range(9)}
        self.e    = {k: 0.0  for k in range(9)}
        self.w    = {k: 0.0 for k in range(9)}
        self.M0   = {k: 0.0 for k in range(9)}
        self.i    = {k: 90.0 for k in range(9)}
        self.Node = {k: 0.0 for k in range(9)}
        self.w_dot= {k: 0.0 for k in range(9)}

        self.K_err    = {k: np.array([0.0,0.0]) for k in range(9)}
        self.P_err    = {k: np.array([0.0,0.0]) for k in range(9)}
        self.e_err    = {k: np.array([0.0,0.0]) for k in range(9)}
        self.w_err    = {k: np.array([0.0,0.0]) for k in range(9)}
        self.M0_err   = {k: np.array([0.0,0.0]) for k in range(9)}
        self.i_err    = {k: np.array([0.0,0.0]) for k in range(9)}
        self.Node_err = {k: np.array([0.0,0.0]) for k in range(9)}
        self.w_dot_err= {k: np.array([0.0,0.0]) for k in range(9)}

        self.K_use    = {k: False for k in range(9)}
        self.P_use    = {k: False for k in range(9)}
        self.e_use    = {k: False for k in range(9)}
        self.w_use    = {k: False for k in range(9)}
        self.M0_use   = {k: False for k in range(9)}
        self.i_use    = {k: False for k in range(9)}
        self.Node_use = {k: False for k in range(9)}
        self.w_dot_use= {k: False for k in range(9)}

        self.K_bound    = {k: np.array([0.0,10000.0]) for k in range(9)}
        self.P_bound    = {k: np.array([0.0,100000.0]) for k in range(9)}
        self.e_bound    = {k: np.array([0.0,0.999]) for k in range(9)}
        self.w_bound    = {k: np.array([0.0, 360.0]) for k in range(9)}
        self.M0_bound   = {k: np.array([0.0, 360.0]) for k in range(9)}
        self.i_bound    = {k: np.array([0.0, 180.0]) for k in range(9)}
        self.Node_bound = {k: np.array([0.0, 360.0]) for k in range(9)}
        self.w_dot_bound= {k: np.array([0.0, 360.0]) for k in range(9)}

        self.K_norm_pr    = {k: np.array([50.0,100.0, False]) for k in range(9)}
        self.P_norm_pr    = {k: np.array([150.0,30.0, False]) for k in range(9)}
        self.e_norm_pr    = {k: np.array([0.0,0.1, False]) for k in range(9)}
        self.w_norm_pr    = {k: np.array([0.0, 90.0, False]) for k in range(9)}
        self.M0_norm_pr   = {k: np.array([0.0, 90.0, False]) for k in range(9)}
        self.i_norm_pr    = {k: np.array([90.0, 90.0, False]) for k in range(9)}
        self.Node_norm_pr = {k: np.array([0.0, 360.0, False]) for k in range(9)}
        self.w_dot_norm_pr= {k: np.array([0.0, 360.0, False]) for k in range(9)}

        self.K_jeff_pr    = {k: np.array([50.0,100.0, False]) for k in range(9)}
        self.P_jeff_pr    = {k: np.array([150.0,30.0, False]) for k in range(9)}
        self.e_jeff_pr    = {k: np.array([0.0,0.1, False]) for k in range(9)}
        self.w_jeff_pr    = {k: np.array([0.0, 90.0, False]) for k in range(9)}
        self.M0_jeff_pr   = {k: np.array([0.0, 90.0, False]) for k in range(9)}
        self.i_jeff_pr    = {k: np.array([90.0, 90.0, False]) for k in range(9)}
        self.Node_jeff_pr = {k: np.array([0.0, 360.0, False]) for k in range(9)}
        self.w_dot_jeff_pr= {k: np.array([0.0, 360.0, False]) for k in range(9)}

        self.K_str    = {k: r'K$_%s$'%chr(98+k)  for k in range(9)}
        self.P_str    = {k: r'P$_%s$'%chr(98+k)  for k in range(9)}
        self.e_str    = {k: r'e$_%s$'%chr(98+k)  for k in range(9)}
        self.w_str    = {k: r'$\omega_%s$'%chr(98+k)  for k in range(9)}
        self.M0_str   = {k: r'MA$_%s$'%chr(98+k)  for k in range(9)}
        self.i_str    = {k: r'i$_%s$'%chr(98+k)  for k in range(9)}
        self.Node_str = {k: r'$\Omega_%s$'%chr(98+k)  for k in range(9)}
        self.w_dot_str= {k: r'$\dot{\omega_%s}$'%chr(98+k)  for k in range(9)}


        #### transit #####
        self.t0      = {k: 2458357.39 for k in range(9)}
        self.pl_a    = {k: 11.0 for k in range(9)}
        self.pl_rad  = {k: 0.145 for k in range(9)}

        self.t0_use      = {k: False for k in range(9)}
        self.pl_a_use    = {k: False for k in range(9)}
        self.pl_rad_use  = {k: False for k in range(9)}

        self.t0_err      = {k: np.array([0.0,0.0])  for k in range(9)}
        self.pl_a_err    = {k: np.array([0.0,0.0])  for k in range(9)}
        self.pl_rad_err  = {k: np.array([0.0,0.0])  for k in range(9)}

        self.t0_bound      = {k: np.array([-10000.0,2600000.0]) for k in range(9)}
        self.pl_a_bound    = {k: np.array([0.0,100.0]) for k in range(9)}
        self.pl_rad_bound  = {k: np.array([0.0,10000.0]) for k in range(9)}

        self.t0_norm_pr      = {k: np.array([0.0,1.0, False]) for k in range(9)}
        self.pl_a_norm_pr     = {k: np.array([10.0,10.0, False]) for k in range(9)}
        self.pl_rad_norm_pr   = {k: np.array([0.1,0.05, False]) for k in range(9)}

        self.t0_jeff_pr      = {k: np.array([0.0,1.0, False]) for k in range(9)}
        self.pl_a_jeff_pr     = {k: np.array([10.0,10.0, False]) for k in range(9)}
        self.pl_rad_jeff_pr   = {k: np.array([0.1,0.05, False]) for k in range(9)}

        self.t0_str      = {k: r't0 $%s$'%chr(98+k) for k in range(9)}
        self.pl_a_str    = {k: r'a/$R_\star$ $%s$'%chr(98+k) for k in range(9)}
        self.pl_rad_str  = {k: r'R/$R_\star$ $%s$'%chr(98+k) for k in range(9)}


        #### Astrometric #####
        
        self.use_ast_hipp_gaia = False
        self.use_ast_imaging   = False 
        
 
        self.ast_alpha     = {k: 0.0 for k in range(1)}
        self.ast_delta      = {k: 0.0 for k in range(1)}
        self.ast_pi        = {k: 0.0 for k in range(1)}
        self.ast_mu_alpha  = {k: 0.0 for k in range(1)}        
        self.ast_mu_delta  = {k: 0.0 for k in range(1)}   
        
        self.ast_alpha_use     = {k: False for k in range(1)}
        self.ast_delta_use      = {k: False for k in range(1)}
        self.ast_pi_use        = {k: False for k in range(1)}
        self.ast_mu_alpha_use  = {k: False for k in range(1)}        
        self.ast_mu_delta_use  = {k: False for k in range(1)}          
        
        self.ast_alpha_err     = {k: np.array([0.0,0.0])  for k in range(1)}
        self.ast_delta_err      = {k: np.array([0.0,0.0])  for k in range(1)}
        self.ast_pi_err        = {k: np.array([0.0,0.0])  for k in range(1)}
        self.ast_mu_alpha_err  = {k: np.array([0.0,0.0])  for k in range(1)}        
        self.ast_mu_delta_err  = {k: np.array([0.0,0.0])  for k in range(1)}          
                  
        self.ast_alpha_bound      = {k: np.array([0.0,360.0])  for k in range(1)}
        self.ast_delta_bound       = {k: np.array([-90.0,90.0])  for k in range(1)}
        self.ast_pi_bound         = {k: np.array([0.0,10000.0])  for k in range(1)}
        self.ast_mu_alpha_bound   = {k: np.array([0.0,10000.0])  for k in range(1)}        
        self.ast_mu_delta_bound   = {k: np.array([-100.0,100.0])  for k in range(1)}          
                             
        self.ast_alpha_norm_pr     = {k: np.array([0.0,3.0, False])  for k in range(1)}
        self.ast_delta_norm_pr      = {k: np.array([0.0,3.0, False])  for k in range(1)}
        self.ast_pi_norm_pr        = {k: np.array([0.0,100.0, False])  for k in range(1)}
        self.ast_mu_alpha_norm_pr  = {k: np.array([0.0,100.0, False])  for k in range(1)}        
        self.ast_mu_delta_norm_pr  = {k: np.array([-100.0,100.0, False])  for k in range(1)}      

        self.ast_alpha_jeff_pr     = {k: np.array([0.0,3.0, False])  for k in range(1)}
        self.ast_delta_jeff_pr      = {k: np.array([0.0,3.0, False])  for k in range(1)}
        self.ast_pi_jeff_pr        = {k: np.array([0.0,100.0, False])  for k in range(1)}
        self.ast_mu_alpha_jeff_pr  = {k: np.array([0.0,100.0, False])  for k in range(1)}        
        self.ast_mu_delta_jeff_pr  = {k: np.array([-100.0,100.0, False])  for k in range(1)}    
          
        self.ast_alpha_str           = {k: r'$\alpha$' for k in range(1)}
        self.ast_delta_str           = {k: r'$\delta$' for k in range(1)}
        self.ast_pi_str              = {k: r'$\pi$' for k in range(1)}
        self.ast_mu_alpha_str        = {k: r'$\mu\alpha$' for k in range(1)}       
        self.ast_mu_delta_str        = {k: r'$\mu\delta$' for k in range(1)}   

        

    def init_hkl(self) :

        #### h,k,l #####
        self.e_sinw        = {k: 0.0 for k in range(9)}
        self.e_cosw        = {k: 0.0 for k in range(9)}
        self.lamb          = {k: 0.0 for k in range(9)}

        self.e_sinw_use    = {k: False for k in range(9)}
        self.e_cosw_use    = {k: False for k in range(9)}
        self.lamb_use    = {k: False for k in range(9)}

        self.e_sinw_err    = {k: np.array([0.0,0.0])  for k in range(9)}
        self.e_cosw_err    = {k: np.array([0.0,0.0])  for k in range(9)}
        self.lamb_err      = {k: np.array([0.0,0.0])  for k in range(9)}

        self.e_sinw_bound  = {k: np.array([-1.0,1.0]) for k in range(9)}
        self.e_cosw_bound  = {k: np.array([-1.0,1.0]) for k in range(9)}
        self.lamb_bound    = {k: np.array([0.0,360.0]) for k in range(9)}

        self.e_sinw_norm_pr = {k: np.array([0.0,0.1, False]) for k in range(9)}
        self.e_cosw_norm_pr = {k: np.array([0.0,0.1, False]) for k in range(9)}
        self.lamb_norm_pr   = {k: np.array([0.0,30.0, False]) for k in range(9)}

        self.e_sinw_jeff_pr = {k: np.array([-1.0,1.0, False]) for k in range(9)}
        self.e_cosw_jeff_pr = {k: np.array([-1.0,1.0, False]) for k in range(9)}
        self.lamb_jeff_pr   = {k: np.array([0.0,360.0, False]) for k in range(9)}

        self.e_sinw_str     = {k: r'$e sin(\omega_%s)$'%chr(98+k) for k in range(9)}
        self.e_cosw_str     = {k: r'$e cos(\omega_%s)$'%chr(98+k) for k in range(9)}
        self.lamb_str       = {k: r'$\lambda_%s$'%chr(98+k) for k in range(9)}

        ######## derived #####################
        self.t_peri = {k: 0.0 for k in range(9)} 

    def init_omega_dot(self) :

        self.omega_dot       = {k: 0.0 for k in range(9)}
        self.omega_dot_err   = {k: np.array([0.0,0.0]) for k in range(9)}
        self.omega_dot_use   = {k: False for k in range(9)}
        self.omega_dot_str   = {k: r'$\omega_%s dot$'%k for k in range(9)}
        self.omega_dot_bounds  = {k: np.array([0.0,10000.0] )for k in range(9)}
        self.omega_dot_norm_pr = {k: np.array([0.0,50.0, False] )for k in range(9)}
        self.omega_dot_jeff_pr = {k: np.array([0.0,360.0, False] )for k in range(9)}


    def init_RV_jitter(self) :

        self.jitt      = {k: 0.0 for k in range(60)}
        self.jitt_err  = {k: np.array([0.0,0.0]) for k in range(60)}
        self.jitt_use  = {k: False for k in range(60)}
        self.jitt_str  = {k: r'RV jitt$_%s$'%str(k+1) for k in range(60)}
        self.jitt_bounds  = {k: np.array([0.0,10000.0] )for k in range(60)}
        self.jitt_norm_pr = {k: np.array([1.0,5.0, False] )for k in range(60)}
        self.jitt_jeff_pr = {k: np.array([1.0,5.0, False] )for k in range(60)}


    def init_RV_offset(self) :

        self.rvoff      = {k: 0.0 for k in range(60)}
        self.rvoff_err  = {k: np.array([0.0,0.0])  for k in range(60)}
        self.rvoff_use  = {k: False for k in range(60)}
        self.rvoff_str  = {k: r'RV off$_%s$'%str(k+1) for k in range(60)}
        self.rvoff_bounds  = {k: np.array([-1000000.0,1000000.0] )for k in range(60)}
        self.rvoff_norm_pr = {k: np.array([0.0,100.0, False] )for k in range(60)}
        self.rvoff_jeff_pr = {k: np.array([0.0,100.0, False] )for k in range(60)}

        self.n_rvdata = 0


    def init_tra_jitter(self) :

        self.tra_jitt      = {k: 0.0 for k in range(60)}
        self.tra_jitt_err  = {k: np.array([0.0,0.0]) for k in range(60)}
        self.tra_jitt_use  = {k: False for k in range(60)}
        self.tra_jitt_str  = {k: r'transit jitt$_%s$'%str(k+1) for k in range(60)}
        self.tra_jitt_bounds  = {k: np.array([-0.2,0.2] )for k in range(60)}
        self.tra_jitt_norm_pr = {k: np.array([0.0,0.1, False] )for k in range(60)}
        self.tra_jitt_jeff_pr = {k: np.array([0.0,0.1, False] )for k in range(60)}


    def init_tra_offset(self) :

        self.tra_off      = {k: 0.0 for k in range(60)}
        self.tra_off_err  = {k: np.array([0.0,0.0])  for k in range(60)}
        self.tra_off_use  = {k: False for k in range(60)}
        self.tra_off_str  = {k: r'transit off$_%s$'%str(k+1) for k in range(60)}
        self.tra_off_bounds  = {k: np.array([-1.0,2.0] )for k in range(60)}
        self.tra_off_norm_pr = {k: np.array([1.0,0.1, False] )for k in range(60)}
        self.tra_off_jeff_pr = {k: np.array([1.0,0.1, False] )for k in range(60)}

    def init_tra_dilution(self) :

        self.tra_dil     = {k: 1.0 for k in range(60)}
        self.tra_dil_err  = {k: np.array([0.0,0.0])  for k in range(60)}
        self.tra_dil_use  = {k: False for k in range(60)}
        self.tra_dil_str  = {k: r'tr. data dilution$_%s$'%str(k+1) for k in range(60)}
        self.tra_dil_bounds  = {k: np.array([0.0,1.0] )for k in range(60)}
        self.tra_dil_norm_pr = {k: np.array([1.0,0.1, False] )for k in range(60)}
        self.tra_dil_jeff_pr = {k: np.array([1.0,0.1, False] )for k in range(60)}

    def init_RV_lintr(self) :

        self.rv_lintr      = 0.0
        self.rv_lintr_err  = [0.0,0.0]
        self.rv_lintr_use  = False
        self.rv_lintr_str  = {k: r'RV lin.tr' for k in range(1)}
        self.rv_lintr_bounds  = {k: np.array([-1.0,1.0]) for k in range(1)}
        self.rv_lintr_norm_pr = {k: np.array([0,0.001, False]) for k in range(1)}
        self.rv_lintr_jeff_pr = {k: np.array([0,0.001, False]) for k in range(1)}


    def init_RV_quadtr(self) :

        self.rv_quadtr      = 0.0
        self.rv_quadtr_err  = [0.0,0.0]
        self.rv_quadtr_use  = False
        self.rv_quadtr_str  = {k: r'RV quad.tr' for k in range(1)}
        self.rv_quadtr_bounds  = {k: np.array([-1.0,1.0]) for k in range(1)}
        self.rv_quadtr_norm_pr = {k: np.array([0,0.001, False]) for k in range(1)}
        self.rv_quadtr_jeff_pr = {k: np.array([0,0.001, False]) for k in range(1)}


    def init_tra_lintr(self) :

        self.tra_lintr      = {k: 0.0 for k in range(60)}
        self.tra_lintr_err  = {k: np.array([0.0,0.0])  for k in range(60)}
        self.tra_lintr_use  = {k: False for k in range(60)}
        self.tra_lintr_str  = {k: r'tra lin.tr$_%s$'%str(k+1) for k in range(60)}
        self.tra_lintr_bounds  = {k: np.array([-1.0,1.0]) for k in range(60)}
        self.tra_lintr_norm_pr = {k: np.array([0,0.001, False]) for k in range(60)}
        self.tra_lintr_jeff_pr = {k: np.array([-0.001,0.001, False]) for k in range(60)}


    def init_tra_quadtr(self) :

        self.tra_quadtr      = {k: 0.0 for k in range(60)}
        self.tra_quadtr_err  = {k: np.array([0.0,0.0])  for k in range(60)}
        self.tra_quadtr_use  = {k: False for k in range(60)}
        self.tra_quadtr_str  = {k: r'tra quad.tr$_%s$'%str(k+1) for k in range(60)}
        self.tra_quadtr_bounds  = {k: np.array([-1.0,1.0]) for k in range(60)}
        self.tra_quadtr_norm_pr = {k: np.array([0,0.001, False]) for k in range(60)}
        self.tra_quadtr_jeff_pr = {k: np.array([-0.001,0.001, False]) for k in range(60)}

 
    def init_TTVs(self):

        # Static

        self.tra_ttv      = {k: [] for k in range(10)}
        self.tra_ttv_err  = {k: [] for k in range(10)}
        self.tra_ttv_use  = {k: [] for k in range(10)}
        self.tra_ttv_str  = {k: [] for k in range(10)}
        self.tra_ttv_bounds  = {k: [] for k in range(10)}
        self.tra_ttv_norm_pr = {k: [] for k in range(10)}
        self.tra_ttv_jeff_pr = {k: [] for k in range(10)}

        # Dynamic
        
        self.epoch_ttv = 2458000.0
        self.epoch_ttv_end = 2459000.0
        self.ttv_dt = 0.02
        self.ttv_times = [self.epoch_ttv,self.ttv_dt,self.epoch_ttv_end] 

        self.epoch_ast = 2458000.0
        self.epoch_ast_end = 2459000.0
        self.ast_dt = 0.02
        self.ast_times = [self.epoch_ttv,self.ttv_dt,self.epoch_ttv_end] 
    
    def init_TTV_set(self):
        
        if len(self.transit_results) <=1:
            print("No transit results stored!")
            return
        
       # for i in range(self.npl):

        for i in range(9):
            if not bool(self.use_planet[i]):
                continue

            tran_times = get_transit_times(self.transit_results, self.P[i], self.t0[i], precise = True)
            
            self.tra_ttv[i] = tran_times[1]
            self.tra_ttv_err[i]  = {k: np.array([0.0,0.0])  for k in range(len(tran_times[1]))}            
            self.tra_ttv_use[i]  = {k: True for k in range(len(tran_times[1]))} 
            self.tra_ttv_str[i]  = {k: r't$_%s$'%str(tran_times[0][k])  for k in range(len(tran_times[1]))}
            self.tra_ttv_bounds[i]  = {k: np.array([tran_times[1][k]-0.2,tran_times[1][k]+0.2])  for k in range(len(tran_times[1]))} 
            self.tra_ttv_norm_pr[i] = {k: np.array([tran_times[1][k],0.1, False])  for k in range(len(tran_times[1]))} 
            self.tra_ttv_jeff_pr[i] = {k: np.array([tran_times[1][k]-0.1,tran_times[1][k]+0.1, False])  for k in range(len(tran_times[1]))}        

        self.get_TTVs[0] = True

    def init_st_mass(self) :

        self.st_mass      = {k: 1.0 for k in range(1)}
        self.st_mass_err  = {k: np.array([0.0,0.0]) for k in range(1)}
        self.st_mass_use  = {k: False for k in range(1)}
        self.st_mass_str  = {k: r'St mass' for k in range(1)}
        self.st_mass_bounds  = {k: np.array([0.01,100.0]) for k in range(1)}
        self.st_mass_norm_pr = {k: np.array([1.0,0.2, False]) for k in range(1)}
        self.st_mass_jeff_pr = {k: np.array([1.0,0.2, False]) for k in range(1)}
        
 
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
        self.get_GP_lnl = True
        self.gps=[]

        self.GP_rot_params = [1.0,10.0,15.0,1.0]# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
        self.GP_rot_err = [0.0,0.0,0.0,0.0]
        self.GP_rot_use = [False,False,False,False]
        self.GP_rot_str = [r'RV GP$_{\rm Rot.}$ Amp.', r'RV GP$_{\rm Rot.}$ timescale', r'RV GP$_{\rm Rot.}$ Period', r'RV GP$_{\rm Rot.}$ fact.']# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
        self.GP_rot_bounds  = {k: np.array([0.0,100000.0]) for k in range(len(self.GP_rot_params))}
        self.GP_rot_norm_pr = {k: np.array([0.0,10.0, False]) for k in range(len(self.GP_rot_params))}
        self.GP_rot_jeff_pr = {k: np.array([0.0,10.0, False]) for k in range(len(self.GP_rot_params))}

        self.GP_sho_params     = [100.0,1.0,0.05]# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
        self.GP_sho_err = [0.0,0.0,0.0]
        self.GP_sho_use = [False,False,False]
        self.GP_sho_str = [r'RV GP$_{\rm SHO}$ $S$', r'RV GP$_{\rm SHO}$ $Q$', r'RV GP$_{\rm SHO}$ $\omega$']# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
        self.GP_sho_bounds     = {k: np.array([0.0,100000.0]) for k in range(len(self.GP_sho_params))}
        self.GP_sho_norm_pr    = {k: np.array([0.0,10.0, False]) for k in range(len(self.GP_sho_params))}
        self.GP_sho_jeff_pr    = {k: np.array([0.0,10.0, False]) for k in range(len(self.GP_sho_params))}

        self.GP_mat_params     = [1.0,1.0,0.001]# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
        self.GP_mat_err = [0.0,0.0,0.0]
        self.GP_mat_use = [False,False,False]
        self.GP_mat_str = [r'RV GP$_{\rm mat 3/2}$ $\sigma$', r'RV GP$_{\rm mat 3/2}$ $\rho$', r'RV GP$_{\rm mat 3/2}$ $eps$']# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
        self.GP_mat_bounds     = {k: np.array([0.0,10.0]) for k in range(len(self.GP_mat_params))}
        self.GP_mat_norm_pr    = {k: np.array([0.0,10.0, False]) for k in range(len(self.GP_mat_params))}
        self.GP_mat_jeff_pr    = {k: np.array([0.0,10.0, False]) for k in range(len(self.GP_mat_params))}

        self.GP_drw_params     = [1.0,1.0]# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
        self.GP_drw_err = [0.0,0.0]
        self.GP_drw_use = [False,False]
        self.GP_drw_str = [r'RV GP$_{\rm DRW}$ $a$', r'RV GP$_{\rm DRW}$ $c$']# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
        self.GP_drw_bounds     = {k: np.array([0.0,10.0]) for k in range(len(self.GP_drw_params))}
        self.GP_drw_norm_pr    = {k: np.array([0.0,10.0, False]) for k in range(len(self.GP_drw_params))}
        self.GP_drw_jeff_pr    = {k: np.array([0.0,10.0, False]) for k in range(len(self.GP_drw_params))}

        self.GP_double_sho_params     = [1.5,3.45,1.3,1.05,0.5, 999]# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
        self.GP_double_sho_err = [0.0,0.0,0.0,0.0,0.0,0.0]
        self.GP_double_sho_use = [False,False,False,False,False,False]
        self.GP_double_sho_str = [r'RV GP$_{\rm dSHO}$ $\sigma$', r'RV GP$_{\rm dSHO}$ $P$', r'RV GP$_{\rm dSHO}$ $Q0$', r'RV GP$_{\rm dSHO}$ $dQ$', r'RV GP$_{\rm dSHO}$ $f$', r'RV GP$_{\rm dSHO}$ $dummy$']# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
        self.GP_double_sho_bounds     = {k: np.array([0.0,100000.0]) for k in range(len(self.GP_double_sho_params))}
        self.GP_double_sho_norm_pr    = {k: np.array([0.0,10.0, False]) for k in range(len(self.GP_double_sho_params))}
        self.GP_double_sho_jeff_pr    = {k: np.array([0.0,10.0, False]) for k in range(len(self.GP_double_sho_params))}

        self.gp_model_curve = {k: 0.0 for k in range(60)}
        self.gp_model_data  = {k: 0.0 for k in range(60)}

        self.gp_kernels = ['SHOKernel','RotKernel','Matern32','dSHOKernel','RealTerm']
        self.gp_kernel = self.gp_kernels[0]


    def init_transit_GP(self):

        self.tra_doGP = False
        self.tra_gps=[]
        self.tra_GP_rot_params = [1.0,10.0,15.0,1.0]# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
        self.tra_GP_rot_err = [0.0,0.0,0.0,0.0]
        self.tra_GP_rot_use = [False,False,False,False]
        self.tra_GP_rot_str = [r'Transit GP$_{\rm Rot.}$ Amp.', r'Transit GP$_{\rm Rot.}$ timescale', r'Transit GP$_{\rm Rot.}$ Period', r'Transit GP$_{\rm Rot.}$ fact.']#we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
        self.tra_GP_rot_bounds  = {k: np.array([0.0,100000.0]) for k in range(len(self.tra_GP_rot_params))}
        self.tra_GP_rot_norm_pr = {k: np.array([0.0,10.0, False]) for k in range(len(self.tra_GP_rot_params))}
        self.tra_GP_rot_jeff_pr = {k: np.array([0.0,10.0, False]) for k in range(len(self.tra_GP_rot_params))}

        self.tra_GP_sho_params     = [100.0,1.0,0.05]# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
        self.tra_GP_sho_err = [0.0,0.0,0.0]
        self.tra_GP_sho_use = [False,False,False]
        self.tra_GP_sho_str = [r'Transit GP$_{\rm SHO}$ $S$', r'Transit GP$_{\rm SHO}$ $Q$', r'Transit GP$_{\rm SHO}$ $\omega$']# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
        self.tra_GP_sho_bounds     = {k: np.array([0.0,100000.0]) for k in range(len(self.tra_GP_sho_params))}
        self.tra_GP_sho_norm_pr    = {k: np.array([0.0,10.0, False]) for k in range(len(self.tra_GP_sho_params))}
        self.tra_GP_sho_jeff_pr    = {k: np.array([0.0,10.0, False]) for k in range(len(self.tra_GP_sho_params))}

        self.tra_GP_mat_params     = [1.0,1.0,0.001]# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
        self.tra_GP_mat_err = [0.0,0.0,0.0]
        self.tra_GP_mat_use = [False,False,False]
        self.tra_GP_mat_str = [r'Transit GP$_{\rm mat 3/2}$ $\sigma$', r'Transit GP$_{\rm mat 3/2}$ $\rho$', r'Transit GP$_{\rm mat 3/2}$ $eps$']# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
        self.tra_GP_mat_bounds     = {k: np.array([0.0,10.0]) for k in range(len(self.tra_GP_mat_params))}
        self.tra_GP_mat_norm_pr    = {k: np.array([0.0,10.0, False]) for k in range(len(self.tra_GP_mat_params))}
        self.tra_GP_mat_jeff_pr    = {k: np.array([0.0,10.0, False]) for k in range(len(self.tra_GP_mat_params))}

        self.tra_GP_drw_params     = [1.0,1.0]# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
        self.GP_drw_err = [0.0,0.0]
        self.tra_GP_drw_use = [False,False]
        self.tra_GP_drw_str = [r'Transit GP$_{\rm DRW}$ $a$', r'Transit GP$_{\rm DRW}$ $c$']# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
        self.tra_GP_drw_bounds     = {k: np.array([0.0,10.0]) for k in range(len(self.tra_GP_drw_params))}
        self.tra_GP_drw_norm_pr    = {k: np.array([0.0,10.0, False]) for k in range(len(self.tra_GP_drw_params))}
        self.tra_GP_drw_jeff_pr    = {k: np.array([0.0,10.0, False]) for k in range(len(self.tra_GP_drw_params))}

        self.tra_GP_double_sho_params     = [1.5,3.45,1.3,1.05,0.5, 999]#  we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
        self.tra_GP_double_sho_err = [0.0,0.0,0.0,0.0,0.0,0.0]
        self.tra_GP_double_sho_use = [False,False,False,False,False,False]
        self.tra_GP_double_sho_str = [r'Transit GP$_{\rm dSHO}$ $\sigma$', r'Transit GP$_{\rm dSHO}$ $P$', r'Transit GP$_{\rm dSHO}$ $Q0$', r'Transit GP$_{\rm dSHO}$ $dQ$', r'Transit GP$_{\rm dSHO}$ $f$', r'Transit GP$_{\rm dSHO}$ $dummy$']# we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
        self.tra_GP_double_sho_bounds     = {k: np.array([0.0,100000.0]) for k in range(len(self.tra_GP_double_sho_params))}
        self.tra_GP_double_sho_norm_pr    = {k: np.array([0.0,10.0, False]) for k in range(len(self.tra_GP_double_sho_params))}
        self.tra_GP_double_sho_jeff_pr    = {k: np.array([0.0,10.0, False]) for k in range(len(self.tra_GP_double_sho_params))}

        self.tra_gp_model_curve = {k: 0.0 for k in range(60)}
        self.tra_gp_model_data  = {k: 0.0 for k in range(60)}

        self.tra_gp_kernels = ['SHOKernel','RotKernel','Matern32','dSHOKernel','RealTerm']
        self.tra_gp_kernel = self.gp_kernels[0]


    def init_transit_params(self):
 

        self.tr_params = batman.TransitParams()       #object to store transit parameters

        # WASP 6
        self.tr_params.t0  = 0.0  #time of inferior conjunction
        self.tr_params.per = 3.36    #orbital period
        self.tr_params.ecc = 0.0
        self.tr_params.w   = 90.0   #longitude of periastron (in degrees)
        self.tr_params.rp  = 0.15   #planet radius (in units of stellar radii)
        self.tr_params.inc = 90. #orbital inclination (in degrees)
        self.tr_params.a   = 15.0  #semi-major axis (in units of stellar radii)
 
        self.tr_params.limb_dark = "quadratic"      #limb darkening model
        self.tr_params.u =  [0.1, 0.3 ]

        # ld_options = ["uniform", "linear", "quadratic", "nonlinear"]
        #ld_coefficients = [[], [0.3], [0.1, 0.3], [0.5, 0.1, 0.1, -0.1]]

        self.tr_params_use = [False, False,False,False,False,False,False]

        self.transit_results = [0.0,0.0,0.0]
        self.tra_model_fact = 10.0

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
        self.evol_Per_arb = {k: [] for k in range(9)}
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
        self.evol_Per = {k: [] for k in range(9)}
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
        self.mcmc_fileoutput=False
        self.mcmc_save_means=False
        self.mcmc_save_median=False
        self.mcmc_save_mode=False
        self.mcmc_save_maxlnL=False
        self.mcmc_save_sampler=True
        self.mcmc_sampler=[]


        self.mcmc_mad = False
        self.mcmc_AMD_stab   = False
        self.mcmc_Nbody_stab = False

        self.mcmc_sample_file = 'mcmc_samples'
        self.mcmc_corner_plot_file = 'cornerplot.pdf'
        self.mcmc_stat = {"mean": [],"median": [],"mode": [],"best": [],"MAD":[]}


    def init_nest_par(self):
        #self.gaussian_ball = 0.0001
        self.live_points_fact = 4
        self.nest_percentile_level = 68.3

        self.ns_samp_method_opt = ['rwalk','slice','unif','rstagger','rslice','hslice']
        self.ns_samp_method = self.ns_samp_method_opt[0]
        
        self.ns_samp_bound_opt = ['multi','single','balls','rslice','cubes','none']
        self.ns_samp_bound = self.ns_samp_bound_opt[0]
        self.ns_pfrac = 1.0

        self.ns_threads=1
        self.Dynamic_nest = False
        self.std_output=False
        self.stop_crit = 0.001
        self.ns_fileoutput=False
        self.ns_save_means=False
        self.ns_save_median=False
        self.ns_save_mode=False
        self.ns_save_maxlnL=False
        self.ns_save_sampler=True
        self.ns_sampler=[]
        
        self.ns_use_stop = True
        self.ns_maxiter = {0:False, 1:10000000}
        self.ns_maxcall = {0:False, 1:10000000}

        self.nest_mad = False
        self.NS_AMD_stab   = False
        self.NS_Nbody_stab = False
        self.nest_weighted = True

        self.nest_sample_file = 'nest_samp_samples'
        self.nest_corner_plot_file = 'nest_samp_cornerplot.pdf'
        self.nest_stat = {"mean": [],"median": [],"mode": [],"best": [],"MAD":[]}


    def init_dynfit_settings(self):
        self.mixed_fit = {0: [False], 1:[1,1,1,1,1,1,1,1,1]}
        self.fitting_method = 'None'
        self.time_step_model = 10.0
        self.dyn_model_accuracy = 1000.0
        self.dyn_model_to_kill = 20.0
        self.kep_model_to_kill = 10.0
        self.master_timeout = 20.0 #86400.0


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

 
 

############################ RV datasets ##########################################
    def add_rv_dataset(self, name, path, rv_idset = 0):
 
 
        rv_JD       = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [0])
        rv_data     = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [1])
        rv_data_sig = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [2])
 

        rv_data_set = np.array([rv_JD,rv_data,rv_data_sig,np.array([rv_idset]*len(rv_JD)),name,path],dtype=object)
 
        ####### for now ###########
        self.rv_data_sets[rv_idset] =  rv_data_set
        self.rv_data_sets_init[rv_idset] = dill.copy(self.rv_data_sets[rv_idset])


        self.ndset = len([x for x in range(60) if len(self.rv_data_sets[x]) != 0])

        self.n_rvdata = sum([len(self.rv_data_sets[x][0]) for x in range(60) if len(self.rv_data_sets[x]) != 0])

        return


    def remove_rv_dataset(self, rv_idset):

        self.rv_data_sets[rv_idset]      = []
        self.rv_data_sets_init[rv_idset] = []

        self.ndset = len([x for x in range(60) if len(self.rv_data_sets[x]) != 0])
        self.n_rvdata = sum([len(self.rv_data_sets[x][0]) for x in range(60) if len(self.rv_data_sets[x]) != 0])

        return

############################ activity datasets ##########################################
    def add_act_dataset(self, name, path, act_idset = 0):


        try:
 
            act_JD_       = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [0])
            act_data_     = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [1])
            act_data_sig_ = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [2])
 
            if len(act_JD_) != len(act_data_) != len(act_data_sig_):
                print("Something is wrong with your activity data file! Please provide a valid activity data file that contains: BJD  y  sigma_y")
                return
            act_JD        = act_JD_[      np.isfinite(act_JD_) & np.isfinite(act_data_) & np.isfinite(act_data_sig_)]
            act_data      = act_data_[    np.isfinite(act_JD_) & np.isfinite(act_data_) & np.isfinite(act_data_sig_)]
            act_data_sig  = act_data_sig_[np.isfinite(act_JD_) & np.isfinite(act_data_) & np.isfinite(act_data_sig_)]

            if len(act_JD) <= 1:
                print("Something is wrong with your activity data file! Perhaps some not all entries are numeric? Please provide a valid activity data file that contains: BJD  y  sigma_y")
                return
        except:
            print("Something is wrong with your activity data file! Please provide a valid activity data file that contains: BJD  y  sigma_y")
            return

        act_data_o_c = act_data            
        act_file_name = file_from_path(path)

        act_data_set = np.array([act_JD,act_data,act_data_sig,act_data_o_c,act_data_o_c,act_data,act_data_sig,act_data_o_c, 1.0, act_file_name],dtype=object)
 
        self.act_data_sets[act_idset]      = act_data_set
        self.act_data_sets_init[act_idset] = dill.copy(self.act_data_sets[act_idset])
 
        return


    def remove_act_dataset(self, act_idset):

        self.act_data_sets[act_idset] = []
        self.act_data_sets_init[act_idset] = []

        return

############################ TTV datasets ##########################################
    def add_ttv_dataset(self, name, path, ttv_idset = 0, planet = 0, use = False):


        try:
            ttv_N_        = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [0])
            ttv_data_     = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [1])
            ttv_data_sig_ = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [2])

            if len(ttv_N_) != len(ttv_data_) != len(ttv_data_sig_):
                print("Something is wrong with your ttv data file! Please provide a ttv data file that contains: N  t_N  sigma_t_N")
                return

            ttv_N         = ttv_N_[       np.isfinite(ttv_N_) & np.isfinite(ttv_data_) & np.isfinite(ttv_data_sig_)]
            ttv_data      = ttv_data_[    np.isfinite(ttv_N_) & np.isfinite(ttv_data_) & np.isfinite(ttv_data_sig_)]
            ttv_data_sig  = ttv_data_sig_[np.isfinite(ttv_N_) & np.isfinite(ttv_data_) & np.isfinite(ttv_data_sig_)]

            if len(ttv_N) == 0:
                print("Something is wrong with your ttv data file! Perhaps some not all entries are numeric? Please provide a ttv data file that contains: N  t_N  sigma_t_N")
                return
        except:
            print("Something is wrong with your ttv data file! Please provide a ttv data file that contains: N  t_N  sigma_t_N")
            return

        ttv_file_name = file_from_path(path)
        ttv_data_set = np.array([ttv_N,ttv_data,ttv_data_sig,planet,use,ttv_file_name],dtype=object)

        self.ttv_data_sets[ttv_idset] = ttv_data_set
        return


    def remove_ttv_dataset(self, ttv_idset):
        self.ttv_data_sets[ttv_idset] = []
        return

############################ Ast datasets ##########################################
    def add_ast_dataset(self, name, path, ast_idset = 0, planet = 0, use = False):
 
        try:
            ast_BJD_        = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [0])
            ast_data_x_     = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [1])
            ast_data_x_sig_ = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [2])
            ast_data_y_     = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [3])
            ast_data_y_sig_ = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [4])


            if len(ast_BJD_) != len(ast_data_x_) != len(ast_data_x_sig_) != len(ast_data_y_) != len(ast_data_y_sig_):
                print("Something is wrong with your Astromtry data file! Please provide an Astromtry data file that contains: BJD  x  x_sigma  y  y_sigma ")
                return

            ast_BJD         = ast_BJD_[       np.isfinite(ast_BJD_) & np.isfinite(ast_data_x_) & np.isfinite(ast_data_y_) & np.isfinite(ast_data_x_sig_) & np.isfinite(ast_data_y_sig_)]
            ast_data_x      = ast_data_x_[    np.isfinite(ast_BJD_) & np.isfinite(ast_data_x_) & np.isfinite(ast_data_y_) & np.isfinite(ast_data_x_sig_) & np.isfinite(ast_data_y_sig_)]
            ast_data_x_sig  = ast_data_x_sig_[np.isfinite(ast_BJD_) & np.isfinite(ast_data_x_) & np.isfinite(ast_data_y_) & np.isfinite(ast_data_x_sig_) & np.isfinite(ast_data_y_sig_)]
            ast_data_y      = ast_data_y_[    np.isfinite(ast_BJD_) & np.isfinite(ast_data_x_) & np.isfinite(ast_data_y_) & np.isfinite(ast_data_x_sig_) & np.isfinite(ast_data_y_sig_)]
            ast_data_y_sig  = ast_data_y_sig_[np.isfinite(ast_BJD_) & np.isfinite(ast_data_x_) & np.isfinite(ast_data_y_) & np.isfinite(ast_data_x_sig_) & np.isfinite(ast_data_y_sig_)]

            if len(ast_BJD) == 0:
                print("Something is wrong with your Astromtry data file! Perhaps some not all entries are numeric? Please provide a Astromtry data file that contains: BJD  x  x_sigma  y  y_sigma")
                return
        except:
            print("Something is wrong with your Astromtry data file! Please provide a Astromtry data file that contains: BJD  x  x_sigma  y  y_sigma")
            return

        ast_file_name = file_from_path(path)
        ast_data_set = np.array([ast_BJD,ast_data_x,ast_data_x_sig,ast_data_y,ast_data_y_sig, planet, use, ast_file_name],dtype=object)

        self.ast_data_sets[ast_idset] = ast_data_set
        return


    def remove_ast_dataset(self, ast_idset):
        self.ast_data_sets[ast_idset] = []
        return

############################ Ast datasets Hipp/Gaia ##########################################
    def add_ast_dataset_hipp_gaia(self, name, path, ast_idset = 0, planet = 0, use = False):
 
# [np.arange(0,len(A8)),frac,A5,A3,A4,A8,A9])

        def read_header_lines(file_path):
            """
            Extract all lines starting with '#' from a file.

            Args:
                file_path (str): Path to the file.

            Returns:
                list: A list of strings containing the header lines.
            """
            header_lines = []
            with open(file_path, 'r') as file:
                for line in file:
                    if line.strip().startswith("#"):
                        header_lines.append(line.strip("#").strip())  # Remove '#' and extra whitespace
            return header_lines



        def parse_hipparcos_data(lines):
            """
            Parses the Hipparcos data lines and returns a structured dictionary
            containing all relevant blocks of headers and values.

            Args:
                lines (list): List of strings containing headers and data rows.

            Returns:
                dict: Dictionary mapping all headers to their corresponding values.
            """
            combined_data = {}  # Final dictionary to store all parsed data
            headers = []
            values = []

            for line in lines:
                if line.strip() == '':
                    # End of a block, process it
                    if headers and values:
                        block_dict = dict(zip(headers, values))
                        combined_data.update(block_dict)
                    headers, values = [], []  # Reset for the next block
                    continue

                # Split the line into parts
                parts = line.split()
                
                # Check if the line is data (numeric or "---") or headers
                if all(x.replace('.', '', 1).replace('-', '', 1).isdigit() or x == "---" for x in parts):
                    values = [float(x) if x.replace('.', '', 1).replace('-', '', 1).isdigit() else x for x in parts]
                else:
                    headers = parts

            # Process the last block if not empty
            if headers and values:
                block_dict = dict(zip(headers, values))
                combined_data.update(block_dict)

            return combined_data

                    
        header = read_header_lines(path)
        
        #print(header)
        hipp_header = parse_hipparcos_data(header)

        #self.header = hipp_header


        try:
        
            ast_BJD_        = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [0])
            ast_data_x_     = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [1])
            ast_data_x_sig_ = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [2])
            ast_data_y_     = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [3])
            ast_data_y_sig_ = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [4])
            ast_data_z_     = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [5])
            ast_data_z_sig_ = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [6])

            if len(ast_BJD_) != len(ast_data_x_) != len(ast_data_x_sig_) != len(ast_data_y_) != len(ast_data_y_sig_):
                print("Something is wrong with your Astromtry data file! Please provide an Astromtry data file that contains the ")
                return

            ast_BJD         = ast_BJD_[       np.isfinite(ast_BJD_) & np.isfinite(ast_data_x_) & np.isfinite(ast_data_y_) & np.isfinite(ast_data_x_sig_) & np.isfinite(ast_data_y_sig_)]
            ast_data_x      = ast_data_x_[    np.isfinite(ast_BJD_) & np.isfinite(ast_data_x_) & np.isfinite(ast_data_y_) & np.isfinite(ast_data_x_sig_) & np.isfinite(ast_data_y_sig_)]
            ast_data_x_sig  = ast_data_x_sig_[np.isfinite(ast_BJD_) & np.isfinite(ast_data_x_) & np.isfinite(ast_data_y_) & np.isfinite(ast_data_x_sig_) & np.isfinite(ast_data_y_sig_)]
            ast_data_y      = ast_data_y_[    np.isfinite(ast_BJD_) & np.isfinite(ast_data_x_) & np.isfinite(ast_data_y_) & np.isfinite(ast_data_x_sig_) & np.isfinite(ast_data_y_sig_)]
            ast_data_y_sig  = ast_data_y_sig_[np.isfinite(ast_BJD_) & np.isfinite(ast_data_x_) & np.isfinite(ast_data_y_) & np.isfinite(ast_data_x_sig_) & np.isfinite(ast_data_y_sig_)]
            ast_data_z      = ast_data_z_[    np.isfinite(ast_BJD_) & np.isfinite(ast_data_x_) & np.isfinite(ast_data_y_) & np.isfinite(ast_data_x_sig_) & np.isfinite(ast_data_y_sig_)]
            ast_data_z_sig  = ast_data_z_sig_[np.isfinite(ast_BJD_) & np.isfinite(ast_data_x_) & np.isfinite(ast_data_y_) & np.isfinite(ast_data_x_sig_) & np.isfinite(ast_data_y_sig_)]

            if len(ast_BJD) == 0:
                print("Something is wrong with your Astromtry data file! Perhaps some not all entries are numeric? Please provide a Astromtry data file that contains the standard Hipparcos data input")
                return
        except:
            print("Something is wrong with your Astromtry data file! Please provide a Astromtry data file that contains the standard Hipparcos data input")
            return


        ast_file_name = file_from_path(path)
        ast_data_set = np.array([ast_BJD,ast_data_x,ast_data_x_sig,ast_data_y,ast_data_y_sig,ast_data_z,ast_data_z_sig, planet, use, hipp_header, ast_file_name],dtype=object)

        self.ast_data_sets_hipp_gaia[ast_idset] = ast_data_set
        return


    def remove_ast_dataset_hipp_gaia(self, ast_idset):
        self.ast_data_sets_hipp_gaia[ast_idset] = []
        return


############################ transit datasets ##########################################
    def add_transit_dataset(self, name, path, tra_idset = 0, PDC = False):
     
  
        if path.endswith(".fits"):

            if pyfits_not_found == True:
                print("""

'astropy.io.fits' not found! Please install 'astropy', and try again with a new start of the 'Exo-Striker'.

""")
                return
            try:
                sc = pyfits.open(path) 
                dd =  sc[1].data
                
                if sc[1].header['TELESCOP'] == 'TESS':
 
                
                    if PDC == True:
                        tra_JD        = dd['TIME'][np.isfinite(dd['TIME']) & np.isfinite(dd['PDCSAP_FLUX']) & np.isfinite(dd['PDCSAP_FLUX_ERR'])]
                        tra_data      = dd['PDCSAP_FLUX'][np.isfinite(dd['TIME']) & np.isfinite(dd['PDCSAP_FLUX']) & np.isfinite(dd['PDCSAP_FLUX_ERR'])]
                        tra_data_sig  = dd['PDCSAP_FLUX_ERR'][np.isfinite(dd['TIME']) & np.isfinite(dd['PDCSAP_FLUX']) & np.isfinite(dd['PDCSAP_FLUX_ERR'])]
                    else:
                        tra_JD        = dd['TIME'][np.isfinite(dd['TIME']) & np.isfinite(dd['SAP_FLUX']) & np.isfinite(dd['SAP_FLUX_ERR'])]
                        tra_data      = dd['SAP_FLUX'][np.isfinite(dd['TIME']) & np.isfinite(dd['SAP_FLUX']) & np.isfinite(dd['SAP_FLUX_ERR'])]
                        tra_data_sig  = dd['SAP_FLUX_ERR'][np.isfinite(dd['TIME']) & np.isfinite(dd['SAP_FLUX']) & np.isfinite(dd['SAP_FLUX_ERR'])]
    
                        print("TESS light curve detected. Automatically adding 2457000.0 to BJD_TBD")
                    
                    tra_JD = tra_JD + 2457000.0
                    tra_airmass_ = np.zeros(len(tra_JD))

                elif sc[1].header['TELESCOP'] == 'CHEOPS':
 
                    tra_JD        = dd['BJD_TIME'][np.isfinite(dd['BJD_TIME']) & np.isfinite(dd['FLUX']) & np.isfinite(dd['FLUXERR'])]
                    tra_data      = dd['FLUX'][np.isfinite(dd['BJD_TIME']) & np.isfinite(dd['FLUX']) & np.isfinite(dd['FLUXERR'])]
                    tra_data_sig  = dd['FLUXERR'][np.isfinite(dd['BJD_TIME']) & np.isfinite(dd['FLUX']) & np.isfinite(dd['FLUXERR'])]

                    print("CHEOPS light curve detected.")
                    tra_airmass_ = np.zeros(len(tra_JD))

                else:
  
                    if PDC == True:
                        tra_JD        = dd['TIME'][np.isfinite(dd['TIME']) & np.isfinite(dd['PDCSAP_FLUX']) & np.isfinite(dd['PDCSAP_FLUX_ERR'])]
                        tra_data      = dd['PDCSAP_FLUX'][np.isfinite(dd['TIME']) & np.isfinite(dd['PDCSAP_FLUX']) & np.isfinite(dd['PDCSAP_FLUX_ERR'])]
                        tra_data_sig  = dd['PDCSAP_FLUX_ERR'][np.isfinite(dd['TIME']) & np.isfinite(dd['PDCSAP_FLUX']) & np.isfinite(dd['PDCSAP_FLUX_ERR'])]
                    else:
                        tra_JD        = dd['TIME'][np.isfinite(dd['TIME']) & np.isfinite(dd['SAP_FLUX']) & np.isfinite(dd['SAP_FLUX_ERR'])]
                        tra_data      = dd['SAP_FLUX'][np.isfinite(dd['TIME']) & np.isfinite(dd['SAP_FLUX']) & np.isfinite(dd['SAP_FLUX_ERR'])]
                        tra_data_sig  = dd['SAP_FLUX_ERR'][np.isfinite(dd['TIME']) & np.isfinite(dd['SAP_FLUX']) & np.isfinite(dd['SAP_FLUX_ERR'])]
    
                    print("Kepler/K2 or other telescope light curve detected.")
                    
                    tra_JD = tra_JD + 2457000.0
                    tra_airmass_ = np.zeros(len(tra_JD))


            except:
                print("Unknown type of .fits file! Please provide a TESS *lc.fits, or a CHEOPS *SCI_COR*.fits file")
                return            

        else:
            try:
                tra_JD_       = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [0])
                tra_data_     = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [1])
                tra_data_sig_ = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [2])
                
                if len(tra_JD_) != len(tra_data_) != len(tra_data_sig_):
                    print("Something is wrong with your transit data file! Please provide a valid transit data that contains: BJD  flux  sigma_flux")
                    return
                    
                tra_JD        = tra_JD_[      np.isfinite(tra_JD_) & np.isfinite(tra_data_) & np.isfinite(tra_data_sig_)]
                tra_data      = tra_data_[    np.isfinite(tra_JD_) & np.isfinite(tra_data_) & np.isfinite(tra_data_sig_)]
                tra_data_sig  = tra_data_sig_[np.isfinite(tra_JD_) & np.isfinite(tra_data_) & np.isfinite(tra_data_sig_)]

                if len(tra_JD) <= 5:
                    print("Something is wrong with your transit data file! Perhaps some not all entries are numeric? Please provide a valid transit data that contains: BJD  flux  sigma_flux")
                    return
            except:
                print("Something is wrong with your transit data file! Please provide a valid transit data that contains: BJD  flux  sigma_flux")
                return

            try:
                tra_airmass_ = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [3])
                
                if len(tra_airmass_) != len(tra_JD_) != len(tra_data_) != len(tra_data_sig_):
                    print("Something is wrong with your transit data file in particular the airmass column! Please provide a valid transit data that contains: BJD  flux  sigma_flux ans airmass")
                else:
                    print("Airmass column detected? ")
            except:
                tra_airmass_ = np.zeros(len(tra_JD))
                pass

        tra_airmass = tra_airmass_ 
        tra_data_o_c = tra_data
        tra_file_name = file_from_path(path)


        tra_data_set = np.array([tra_JD,tra_data,tra_data_sig,tra_airmass,tra_data_o_c,tra_data,tra_data_sig,tra_data_o_c, 1.0, True, False, tra_file_name],dtype=object)
        

        self.tra_off_use[tra_idset] = True
        self.tra_jitt_use[tra_idset] = True     

 
        self.tra_data_sets[tra_idset] = tra_data_set
        self.tra_data_sets_init[tra_idset] = dill.copy(self.tra_data_sets[tra_idset])

        return


    def remove_transit_dataset(self, tra_idset):

        self.tra_data_sets[tra_idset] = []
        self.tra_data_sets_init[tra_idset] = []

        self.tra_off_use[tra_idset] = False
        self.tra_jitt_use[tra_idset] = False     


        return


################ Legacy Code!!! TB removed/updated/replaced ######################


    def add_planet(self,K=50,P=100,e=0,w=0,M0=180,i=90,cap=0,useK=True,useP=True,usee=False,usew=False,useM0=True,usei=False, usecap=False, index = None):

        pl_ind = [0,1,2,3,4,5,6,7,8,9]

       # print(index)
        if index != None and int(index) in pl_ind:
            ind = index
        else:
            ind = self.npl       
    
        self.P[ind] = P
        self.K[ind] = K
        self.e[ind] = e
        self.w[ind] = w
        self.M0[ind] = M0
        self.i[ind] = i
        self.Node[ind] = cap

        self.P_use[ind] = useP
        self.K_use[ind] = useK
        self.e_use[ind] = usee
        self.w_use[ind] = usew
        self.M0_use[ind] = useM0
        self.i_use[ind] = usei
        self.Node_use[ind] = usecap

        self.use_planet[ind] = 1

        self.npl=self.npl+1
        return



    def add_dataset(self,name,path,offset,jitter,useoffset=True,usejitter=True, index =None):



        if index == None:
            empt_rvs_files = min([x for x in range(60) if len(self.rv_data_sets[x]) == 0], default=0) 
            index = int(empt_rvs_files)

        self.add_rv_dataset(name, path, rv_idset = int(index))

        self.rvoff_use[index] = True
        self.jitt_use[index] = True     

        if self.epoch == 0:
            self.update_epoch(min(self.rv_data_sets[index][0]))

        return


    def remove_planet_old(self,planet):
        if not (planet<self.npl):
            warnings=Warning_log(['Planet index outside of range.'],'Removing planet %d'%planet)
            warnings.print_warning_log()
        else:


            self.P_use[planet] = False
            self.K_use[planet] = False
            self.e_use[planet] = False
            self.w_use[planet] = False
            self.M0_use[planet] = False
            self.i_use[planet] = False
            self.Node_use[planet] = False

            self.use_planet[planet] = 0
            self.npl=self.npl-1
        return


    def remove_planet(self,planet):
 

        self.P_use[planet] = False
        self.K_use[planet] = False
        self.e_use[planet] = False
        self.w_use[planet] = False
        self.M0_use[planet] = False
        self.i_use[planet] = False
        self.Node_use[planet] = False

        self.use_planet[planet] = 0
        if not self.npl <= 0:
            self.npl=self.npl-1
        return



    def remove_dataset(self,number):
 
        #### new stuff, TB fixed!             
        self.remove_rv_dataset(number)
        self.rvoff[number] = 0.0
        self.rvoff_use[number] = False
        self.rvoff_err[number] = np.array([0.0,0.0])
        self.rvoff_str[number] = r'RV off$_%s$'%str(number+1)
        self.rvoff_bounds[number] = np.array([-1000000.0,1000000.0] )
        self.rvoff_norm_pr[number] = np.array([0.0,100.0, False])
        self.rvoff_jeff_pr[number] = np.array([0.0,100.0, False] ) 

        self.jitt[number] = 0.0
        self.jitt_use[number] = False
        self.jitt_err[number] = np.array([0.0,0.0])
        self.jitt_str[number] = r'RV jitt$_%s$'%str(number+1)
        self.jitt_bounds[number] = np.array([0.0,10000.0] )
        self.jitt_norm_pr[number] = np.array([1.0,5.0, False])
        self.jitt_jeff_pr[number] = np.array([1.0,5.0, False] ) 
 
#        for i in range(number,9):
#            self.rv_data_sets[i] = self.rv_data_sets[i+1]        
        return


    def update_mod_dynamical(self, mod_dynamical):
        self.mod_dynamical=mod_dynamical
        return


    def add_MAROONX_dataset(self, name, path, offset=0, jitter= 0, split = False):

       dirname, basename = os.path.split(path)

#       MAROONX_data = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0 )
#       num_cols,num_rows = MAROONX_data.shape 
      # print(num_rows, num_cols)
    #   return
#       if num_cols !=50:
#           print("Not a MAROON-X file!") 
#           return                   

       BJD = np.genfromtxt("%s"%(path), skip_header=0, unpack=True, skip_footer=0, usecols = [0])
      # indices = np.where(BJD > 2457161.5)

       fo = open(path, "r")
       lines = fo.readlines()
       fo.close()


       seasons = [[2459358.5,2459368.5],[2459402.5,2459450.5],[2459514.5,2459541.5],[2459541.5001,2460442.5]]

       name_file = []
       path_file = []

       for z in range(len(seasons)):

           name_file.append('%s_S_%s.dat'%(basename[:-5],z+1))
           path_file.append('datafiles/%s'%name_file[z])


           out1 = open('%s'%path_file[z], 'w')

           for i in range(len(lines)):

               line = lines[i].split()
               if len(line) == 0:
                   continue
               if seasons[z][0] <= float(line[0]) <= seasons[z][1]:
                   out1.write(lines[i])
 

           out1.close()
 

           if split == True and os.stat(path_file[z]).st_size != 0:
               self.add_dataset(name_file[z],path_file[z],offset,jitter,useoffset=True,usejitter=True)
               self.rvoff_bounds[len(self.filelist.files)-1] = [-100.0,100.0]
               self.rvoff_norm_pr[len(self.filelist.files)-1] = [0.0,10.0,True]
               self.jitt_bounds[len(self.filelist.files)-1] = [0.0,10.0]
               #self.jitt_norm_pr[len(self.filelist.files)-1] = [0.0,10.0,True]

           elif split == False and z == 0:
               self.add_dataset(name,path,offset,jitter,useoffset=True,usejitter=True)
               self.rvoff_bounds[len(self.filelist.files)-1] = [-100.0,100.0]
               self.rvoff_norm_pr[len(self.filelist.files)-1] = [0.0,10.0,True]
               self.jitt_bounds[len(self.filelist.files)-1] = [0.0,10.0] 

    def add_RVbank_dataset(self, name, path, offset=0, jitter= 0, split = False):

       dirname, basename = os.path.split(path)

       HARPS_RVBank_data = np.genfromtxt("%s"%(path),skip_header=1, unpack=True,skip_footer=0 )
       num_cols,num_rows = HARPS_RVBank_data.shape 
     #  print(num_rows, num_cols)
    #   return
       if num_cols !=57:
           print("Not a HARPS_RVBank ver02 file!") 
           return                   

       BJD = np.genfromtxt("%s"%(path), skip_header=1, unpack=True, skip_footer=0, usecols = [3])
      # indices = np.where(BJD > 2457161.5)

       fo = open(path, "r")
       lines = fo.readlines()[1:]
       fo.close()

       name0 = '%s_all.dat'%basename[:-4]
       name1 = '%s_pre.dat'%basename[:-4]
       name2 = '%s_post.dat'%basename[:-4]
       path0 = 'datafiles/%s'%name0
       path1 = 'datafiles/%s'%name1
       path2 = 'datafiles/%s'%name2

       out0 = open('%s'%path0, 'w')
       out1 = open('%s'%path1, 'w')
       out2 = open('%s'%path2, 'w')


       for i in range(len(lines)):

           line = lines[i].split()
           if len(line) == 0:
               continue

           out0.write(" ".join(line[3:]) + "\n")

           if float(line[3]) <= 2457161.5:
               out1.write(" ".join(line[3:]) + "\n")
           elif float(line[3]) > 2457161.5:
               out2.write(" ".join(line[3:]) + "\n")

       out0.close()
       out1.close()
       out2.close()

       if split == True:
           if len(BJD[BJD <= 2457161.5]) !=0:
               self.add_dataset(name1,path1,offset,jitter,useoffset=True,usejitter=True)
           if len(BJD[BJD > 2457161.5]) !=0:
               self.add_dataset(name2,path2,offset,jitter,useoffset=True,usejitter=True)
       else:
           self.add_dataset(name0,path0,offset,jitter,useoffset=True,usejitter=True)
       
       data_set_name = ['CRX-%s'%basename[:-4],'dLW-%s'%basename[:-4],'Halpha-%s'%basename[:-4],'NaDI-%s'%basename[:-4],'NaDII-%s'%basename[:-4],'BIS-%s'%basename[:-4],'Contrast-%s'%basename[:-4],'FWHM-%s'%basename[:-4],"RHKp-%s"%basename[:-4]]


       for i in range(5):

           act_ind = 3+ 11 + (2*i)
           act_data     = np.genfromtxt("%s"%(path),skip_header=1, unpack=True,skip_footer=0, usecols = [act_ind])
           act_data = act_data - np.mean(act_data)

           act_data_sig = np.genfromtxt("%s"%(path),skip_header=1, unpack=True,skip_footer=0, usecols = [act_ind+1])
           #act_data_set = np.array([BJD,act_data,act_data_sig,data_set_name[i]])
           act_data_o_c = act_data                       
           act_data_set = np.array([BJD,act_data,act_data_sig,act_data_o_c,act_data_o_c,act_data,act_data_sig,act_data_o_c, 1.0, data_set_name[i]],dtype=object)


           self.act_data_sets[i] = act_data_set
           self.act_data_sets_init[i] = act_data_set

       #z = 0
       for ii in range(3):

           act_ind = 3+ 22 + ii
           act_data     = np.genfromtxt("%s"%(path),skip_header=1, unpack=True,skip_footer=0, usecols = [act_ind])
           act_data = act_data - np.mean(act_data)
           if ii == 1:
               act_data = act_data * 1000.0
              # print(act_data[0:10])

           act_data_sig = act_data*0.01
           i = i +1
           #act_data_set = np.array([BJD,act_data,act_data_sig,data_set_name[i]])
           act_data_o_c = act_data                       
           act_data_set = np.array([BJD,act_data,act_data_sig,act_data_o_c,act_data_o_c,act_data,act_data_sig,act_data_o_c, 1.0, data_set_name[i]],dtype=object)

           self.act_data_sets[i] = act_data_set
           self.act_data_sets_init[i] = act_data_set

       for ii in range(1):

           act_ind = 53 + (2*ii)
           act_data     = np.genfromtxt("%s"%(path),skip_header=1, unpack=True,skip_footer=0, usecols = [act_ind])
           act_data = act_data - np.mean(act_data)
           act_data_sig = np.genfromtxt("%s"%(path),skip_header=1, unpack=True,skip_footer=0, usecols = [act_ind+1])

           i = i +1
           #act_data_set = np.array([BJD,act_data,act_data_sig,data_set_name[i]])
           act_data_o_c = act_data                       
           act_data_set = np.array([BJD,act_data,act_data_sig,act_data_o_c,act_data_o_c,act_data,act_data_sig,act_data_o_c, 1.0, data_set_name[i]],dtype=object)

           self.act_data_sets[i] = act_data_set


           self.act_data_sets[i] = act_data_set
           self.act_data_sets_init[i] = act_data_set    

    def BIC(self):
        if len(self.fit_results.jd) != 0:
            BIC = self.fit_results.mfit*np.log(len(self.fit_results.jd)) -2*self.loglik
        else:
            BIC = 0
        return BIC

    def AIC(self):
        if len(self.fit_results.jd) != 0:
            AIC = -2*self.loglik - 2*self.fit_results.mfit
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

 
 
 

    def mass_semimajor2(self, mass_type="J"): # calculates planet masses (in Jupiter masses) and orbit semimajor axes (in astronomical units)
        # declaring some constants which will be useful here
        THIRD = 1.0/3.0
        GMSUN = 1.32712440018e20
        AU=1.49597892e11

        # arrays for storing mass and semimajor information
        mass  = np.zeros(self.npl+1)
        ap    = np.zeros(self.npl)
        mtotal = self.params.stellar_mass
        mass[0] = self.params.stellar_mass

        pl_mass = np.zeros(self.npl)
        mpold = pl_mass
        f = 5e-5 # iteration step



        for i in range(9):
            if not bool(self.use_planet[i]):
                continue
        #for i in range(self.npl):

            T = self.P[i]*86400.0

            # we need innitial guess for each planet mass
            dm = 2

            mass[i+1] = abs(self.K[i])*(T*(self.params.stellar_mass)**2.0/(TAU*GMSUN))**THIRD * np.sqrt(1.0-self.e[i]**2.0)/abs(np.sin(np.radians(self.i[i])))

            # mtotal is the sum of the stellar mass and masses of all planets computed so far, plus the current estimate for the mass of the planet we want to compute, we increase it as we go along
           # mtotal = mtotal + mass[i+1]
           # mpold = mass[i+1] # we need to keep track of the current estimate for the planet mass, to check if the new one is higher or not

            mass[0] = mtotal
            mpold[i] = 0
            while (dm >= f):
 
                if i == 0:
                    mtotal = self.params.stellar_mass
                    mass[i+1] = abs(self.K[i])*(T*(self.params.stellar_mass + mpold[i])**2.0/(TAU*GMSUN))**THIRD * np.sqrt(1.0-self.e[i]**2.0)/abs(np.sin(np.radians(self.i[i])))
                else:
                    mtotal = self.params.stellar_mass
                    for j in range(i):
                        mtotal = mtotal + mass[j+1]
                    mass[i+1] = abs(self.K[i])*(T*(mtotal + mpold[i])**2.0/(TAU*GMSUN))**THIRD * np.sqrt(1.0-self.e[i]**2.0)/abs(np.sin(np.radians(self.i[i])))

 
                dm = abs((mass[i+1] - mpold[i])/mass[i+1] )
                mpold[i] =  mass[i+1]


            ap[i] = (GMSUN*(mtotal + mass[i+1])*(T/2.0*np.pi)**2)**THIRD



            ap[i] = ap[i]/AU # to be in AU
            if mass_type=="J":
                pl_mass[i] = mass[i+1]*1047.5654817267318 # to be in Jup. masses
            elif  mass_type=="E":
                pl_mass[i] = mass[i+1]*1047.5654817267318 * 317.82838 # to be in Earth. masses
            else:
                pl_mass[i] = mass[i+1]

        self.masses = pl_mass
        self.semimajor = ap

        return

    def mass_semimajor(self): # calculates planet masses (in Jupiter masses) and orbit semimajor axes (in astronomical units)
        # declaring some constants which will be useful here
        THIRD = 1.0/3.0
        GMSUN = 1.32712440018e20
        AU=1.49597892e11

        # arrays for storing mass and semimajor information
        mass  = np.zeros(self.npl+1)
        ap    = np.zeros(self.npl)
        mtotal = self.params.stellar_mass
        mass[0] = self.params.stellar_mass
        f = 5e-5 # iteration step
        for i in range(self.npl):
            T = self.P[i]*86400.0
            dm = 0
            # we need innitial guess for each planet mass
            mass[i+1] = abs(self.K[i])*(T*(self.params.stellar_mass)**2.0/(TAU*GMSUN))**THIRD * np.sqrt(1.0-self.e[i]**2.0)/abs(np.sin(np.radians(self.i[i])))
            # mtotal is the sum of the stellar mass and masses of all planets computed so far, plus the current estimate for the mass of the planet we want to compute, we increase it as we go along
            mtotal = mtotal + mass[i+1]
            mpold = mass[i+1] # we need to keep track of the current estimate for the planet mass, to check if the new one is higher or not

            # This is a simple iteration to solve for planet mass
            while (dm <= 0): # if the new estimate for planet mass is lower than the old one, we finish an iteration
                mass[i+1] = abs(self.K[i])*(T*mtotal**2.0/(TAU*GMSUN))**THIRD * np.sqrt(1.0-self.e[i]**2.0)/abs(np.sin(np.radians(self.params.planet_params[7*i+5]))) # this will become the planet mass at the last run of the loop,until then it is only used for comparison in the while condition, whereas mpold is the estimate for planet mass.
                dm = (mpold - mass[i+1])
                mpold =  mpold + f
                mtotal = mtotal + f

            mtotal = mtotal-dm # the last part of the sum was mpold, now we want mass[m+1]
            ap[i] = (GMSUN*mtotal*(T/TAU)**2)**THIRD
        # 1 047.92612 = solar mass / jupiter mass, to get output in correct units
        self.masses = mass[1:]*1047.5654817267318
        self.semimajor = ap/AU
        #print(mtotal,mpold,dm,THIRD,GMSUN,AU,self.npl)
        return
 
 
    # sort planets by one of the parameters (K,P,e,w,M0)
    def sort_by_param(self,param,reverse=True):
          if (self.npl>0): # no point in sorting if there's no planets
            sort = convert_array_to_int(np.array(sorted(range(self.npl), key=lambda k: self.params.planet_params[7*k+param], reverse=reverse))) # find a permutation sorting planets by chosen parameter
            # permutate all parameters which refer to planets, according to the permutation found above
            planets_to_sort=np.array([self.params.planet_params[n:n+7] for n in range(0,7*self.npl,7)]) # group parameters for a single planet together for sorting
            planets_to_sort=planets_to_sort[sort] # sort grouped array according to the found permutation
            for i in range(10-self.npl):
                planets_to_sort = np.concatenate([planets_to_sort, [np.zeros(7)]], axis=0)
            self.params=PlanetParamsWrap(self.params.offsets,self.params.jitters,np.concatenate(planets_to_sort),self.params.linear_trend,self.params.stellar_mass) # now come back to 1d array
            # now do the same with use flags
            use_planets_to_sort=np.array([self.use.use_planet_params[n:n+7] for n in range(0,7*self.npl,7)]) # group parameters for a single planet together for sorting
            use_planets_to_sort=use_planets_to_sort[sort] # sort grouped array according to the found permutation
            self.use=use_flags(self.rvoff_use,self.jitt_use,np.concatenate(use_planets_to_sort),self.use.use_linear_trend,self.use.use_stellar_mass) # now come back to 1d array
            # now do the same with errors, only if they are established (fit performed)
 
            if (self.fit_performed):
                planet_errors_to_sort=np.array([self.param_errors.planet_params_errors[n:n+7] for n in range(0,7*self.npl,7)]) # group parameters for a single planet together for sorting
                planet_errors_to_sort=planet_errors_to_sort[sort] # sort grouped array according to the found permutation
                for i in range(10 - self.npl):
                    planet_errors_to_sort = np.concatenate([planet_errors_to_sort, [np.zeros([7, 2])]], axis=0)
                self.param_errors=PlanetParamsErrorsWrap(np.array([[i[0] for i in self.param_errors.offset_errors]]*2).T,
                                                         np.array([[i[0] for i in self.param_errors.jitter_errors]]*2).T,
                                                         np.array([[i[0] for i in np.concatenate(planet_errors_to_sort)]]*2).T,
                                                         np.array([self.param_errors.linear_trend_error[0]]*2).T,
                                                         self.param_errors.stellar_mass_error) # now come back to 1d array
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

            self.omega_dot = self.fit_results.omega_dot
            self.omega_dot_err = self.fit_results.omega_dot_err


            
            self.rv_lintr = self.fit_results.rv_lintr 
            self.rv_lintr_err = self.fit_results.rv_lintr_err

            self.rv_quadtr = self.fit_results.rv_quadtr
            self.rv_quadtr_err = self.fit_results.rv_quadtr_err


            for i in range(9):
                if not bool(self.use_planet[i]):
                    continue

                self.P[i] = self.fit_results.P[i][0]
                self.K[i] = self.fit_results.K[i][0]
 

                if self.hkl == True:
                    self.e_sinw[i]  = self.fit_results.e[i][0]
                    self.e_cosw[i]  = self.fit_results.w[i][0]
                    self.lamb[i]    =  self.M0[i] = self.fit_results.M0[i][0]

                    self.e[i]   = np.sqrt(self.e_sinw[i]**2 + self.e_cosw[i]**2)
                    self.w[i]   = np.degrees(np.arctan2(np.radians(self.e_sinw[i]),np.radians(self.e_cosw[i])))
                    self.M0[i]  = (self.lamb[i] - self.w[i])%360.0
                else:

                    self.e[i] = self.fit_results.e[i][0]
                    self.w[i] = self.fit_results.w[i][0]
                    self.M0[i] = self.fit_results.M0[i][0]

                    self.e_sinw[i] = self.e[i]*np.sin(np.radians(self.w[i]))
                    self.e_cosw[i] = self.e[i]*np.cos(np.radians(self.w[i]))
                    self.lamb[i]   = (self.w[i] + self.M0[i])%360.0


                self.i[i] = self.fit_results.incl[i][0]
                self.Node[i] = self.fit_results.cap0m[i][0]


                self.P_err[i] = np.array([self.fit_results.P[i][1],self.fit_results.P[i][1]])
                self.K_err[i] = np.array([self.fit_results.K[i][1],self.fit_results.K[i][1]])

                if self.hkl == True:
                    self.e_sinw_err[i] = np.array([self.fit_results.e[i][1],self.fit_results.e[i][1]])
                    self.e_cosw_err[i] = np.array([self.fit_results.w[i][1],self.fit_results.w[i][1]])
                    self.lamb_err[i] = np.array([self.fit_results.M0[i][1],self.fit_results.M0[i][1]])

                    self.e_err[i] = np.array([0,0]) # np.array([self.fit_results.e[i][1],self.fit_results.e[i][1]])
                    self.w_err[i] = np.array([0,0]) # np.array([self.fit_results.w[i][1],self.fit_results.w[i][1]])
                    self.M0_err[i] = np.array([0,0]) # np.array([self.fit_results.M0[i][1],self.fit_results.M0[i][1]])
                else:
                    self.e_err[i] = np.array([self.fit_results.e[i][1],self.fit_results.e[i][1]])
                    self.w_err[i] = np.array([self.fit_results.w[i][1],self.fit_results.w[i][1]])
                    self.M0_err[i] = np.array([self.fit_results.M0[i][1],self.fit_results.M0[i][1]])

                    self.e_sinw_err[i] = np.array([0,0]) # np.array([self.fit_results.e[i][1],self.fit_results.e[i][1]])
                    self.e_cosw_err[i] = np.array([0,0]) # np.array([self.fit_results.w[i][1],self.fit_results.w[i][1]])
                    self.lamb_err[i] = np.array([0,0]) # np.array([self.fit_results.M0[i][1],self.fit_results.M0[i][1]])

                self.i_err[i] = np.array([self.fit_results.incl[i][1],self.fit_results.incl[i][1]])
                self.Node_err[i] = np.array([self.fit_results.cap0m[i][1],self.fit_results.cap0m[i][1]])

                self.masses[i] =self.fit_results.mass[i]
                self.semimajor[i] =self.fit_results.a[i]


            for i in range(20):
 
                self.rvoff[i] = self.fit_results.offsets[i]  
                self.jitt[i] = abs(self.fit_results.jitters[i])  
 
                self.rvoff_err[i] = self.fit_results.offset_errors[i]    
                self.jitt_err[i] = self.fit_results.jitter_errors[i]    


        ##################### New stuff to be added here ######################
            for i in range(9):
                #self.t_peri[i] = (float(self.epoch) - (np.radians(self.M0[i])/(2*np.pi))*self.P[i] ) # epoch  - ((ma/TWOPI)*a[1])
                self.t_peri[i] = transit_tperi(self.P[i],self.e[i], self.w[i], self.M0[i] ,self.epoch)[0]

        return


    ### this function is a wrapper calling a fortran program to fit parameters in keplerian mode by minimizing chi^2   WORK IN PROGRESS ON THAT ONE!

    def fitting(self, minimize_fortran=True, minimize_loglik=False, fileinput=False, doGP=False, kernel_id=-1, filename='Kep_input', outputfiles=[1,1,1], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 5,model_min =0): # run the fit which will either minimize chi^2 or loglik.
        '''
         eps, dt - accuracy and step for the integration in dynamical case, in keplerian case ignored, just given for input consistency
         which value to minimize (used for calling an appropriate fortran program)
        '''

        eps = eps * 1e-13
        dt  = dt  * 86400.0

#        check_temp_RV_file(self)

        # bug fix. TBF
        self.bound_error = False
        self.bound_error_msg = ""
        #########################

        if (minimize_loglik):
            minimized_value = 'loglik'
        else:
            minimized_value = 'chi2'

        # which mode to use, analogously as above

        if(self.mod_dynamical):
            mod='dyn'
        else:
            mod='kep'

        hkl = int(self.hkl)
#        from .rvmod import Rvfit
        #rvmod = Rvfit()
 
        def create_args():
            rv_off = []
            rv_jit = []

            for i in range(len(self.rv_data_sets)):
                if len(self.rv_data_sets[i]) == 0:
                    continue
                rv_off.append([self.rvoff[i], int(self.rvoff_use[i])])
                rv_jit.append([self.jitt[i], int(self.jitt_use[i])])

            array_npl = []

            #for i in range(self.npl): # K,P,e,w,M,i,cap0m for each planet, and information which ones we use
            for i in range(9):
                if not bool(self.use_planet[i]):
                    continue 
                array_npl.append([[self.K[i],    int(self.K_use[i])],
                                  [self.P[i],    int(self.P_use[i])],
                                  [self.e[i],    int(self.e_use[i])],
                                  [self.w[i],    int(self.w_use[i])],
                                  [self.M0[i],   int(self.M0_use[i])],
                                  [self.i[i],    int(self.i_use[i])],
                                  [self.Node[i],    int(self.Node_use[i])],
                           [self.omega_dot[i], int(self.omega_dot_use[i])],
                           [0, 0],
                           [0, 0],
                           [0, 0],
                           [0, 0],
                           [0, 0],
                           [0, 0],
                           [self.e_sinw[i], int(self.e_sinw_use[i])],
                           [self.e_cosw[i], int(self.e_cosw_use[i])],
                           [self.lamb[i], int(self.lamb_use[i])]])


            rvs_files = []
            jd, rvs, sig, ids = [], [], [], []
 
            for z in range(len(self.rv_data_sets)):
    
                if len(self.rv_data_sets[z]) == 0:
                    continue
                else:
                    jd = np.append(jd, self.rv_data_sets[z][0], axis = 0)
                    rvs = np.append(rvs, self.rv_data_sets[z][1], axis = 0)
                    sig = np.append(sig, self.rv_data_sets[z][2], axis = 0)
                    ids = np.append(ids, self.rv_data_sets[z][3]+1, axis = 0)
 
            rvs_files = np.array([jd, rvs, sig, ids]).T
 
            self.rvmod.init_args(rvs_files, array_npl, self.epoch, hkl,
                         get_RV=outputfiles[0], get_best_par=outputfiles[1], get_fit_model=outputfiles[2],
                         rv_jitt=rv_jit, rv_ofset=rv_off, npl_pos = self.use_planet,
                         lin_trend_in=[self.rv_lintr,int(self.rv_lintr_use)],
                         quad_trend_in=[self.rv_quadtr,int(self.rv_quadtr_use)],
                         dyn_eps=eps, dyn_dt=dt, amoeba_iter=amoeba_starts, timeout=fortran_kill,
                         rv_model_npoints=self.model_npoints, rv_model_max=self.model_max, rv_model_min=self.model_min,
                         rv_gr_flag=int(self.gr_flag), stellar_mass=self.params.stellar_mass,
                         ndset=self.ndset, ndata=None, nplanet=None,
                         dyn_planets=None,coplar_inc=int(self.copl_incl),jit_flag = int(self.jit_flag))

        program = '%s/lib/fr/%s_%s' % (self.cwd, minimized_value, mod)
        if minimize_fortran == True and doGP ==False:

            create_args()


            self = run_Fort(self,mod,minimize_loglik=minimize_loglik)

        else:
            self = run_SciPyOp(self, kernel_id=kernel_id)
            create_args()
            flag = self.rvmod.run_amoeba(mod)

 
        if self.flag==1: # or self.rtg[0] == True

            self.fit_results = dill.copy(self.rvmod)
            self.stat_saved=self.fit_results.stat_array_saved

            self.model_saved=bool(outputfiles[2])
            self.fit_performed=True
            if(self.fit_results.stat_array_saved):
                self.fitting_method=program
            self.update_with_fit_results()
#            self.correct_elements() #because amoeba might make things wrong here
 

#        del re, Rvfit
#        gc.collect()

        if (return_flag):
            return self.flag
        else:
            return


    # if you want to temporarily change use flags and run a new fit, use this wrapper
    def quick_overwrite_use_and_fit(self,useflags,minimize_loglik=False, fileinput=False, filename='Kep_input', outputfiles=[1,1,0], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 500, model_min=0):
        oldflags=self.overwrite_use(useflags,save=True)
        
        if self.type_fit["Transit"] == True or self.type_fit["TTV"] == True or self.type_fit["AST"] == True:
            minimize_fortran = False
        else:
            minimize_fortran = True
            
        
        self.fitting(minimize_loglik=minimize_loglik,fileinput=fileinput,filename=filename, minimize_fortran=minimize_fortran,
        amoeba_starts=amoeba_starts,outputfiles=outputfiles, eps=eps,dt=dt,timeout_sec=timeout_sec,print_stat=print_stat, fortran_kill=fortran_kill )
        self.overwrite_use(oldflags,save=False)
        return

    # if you want to temporarily change initial parameters and run a new fit, use this wrapper
    def quick_overwrite_params_and_fit(self,params,minimize_loglik=True, fileinput=False, filename='Kep_input',outputfiles=[1,1,0], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False):
        oldparams=self.overwrite_params(params,save=True)
        
        if self.type_fit["Transit"] == True or self.type_fit["TTV"] == True or self.type_fit["AST"] == True:
            minimize_fortran = False
        else:
            minimize_fortran = True        
        
        if (return_flag):
            flag=self.fitting(minimize_loglik=minimize_loglik,fileinput=fileinput,filename=filename, minimize_fortran=minimize_fortran,
            amoeba_starts=amoeba_starts,outputfiles=outputfiles,eps=eps,dt=dt,timeout_sec=timeout_sec,print_stat=print_stat, return_flag=True, fortran_kill=fortran_kill, npoints=self.model_npoints, model_max = self.model_max, model_min = self.model_min)
            self.overwrite_params(oldparams,save=False)
            return flag
        else:
            self.fitting(minimize_loglik=minimize_loglik,fileinput=fileinput,filename=filename,amoeba_starts=amoeba_starts,outputfiles=outputfiles,eps=eps,dt=dt,timeout_sec=timeout_sec,print_stat=print_stat, fortran_kill=fortran_kill, npoints=self.model_npoints, model_max = self.model_max, model_min = self.model_min)
            self.overwrite_params(oldparams,save=False)
            return

    def minimize_one_param_K(self,planet,minimize_loglik=True, fileinput=False, filename='Kep_input', outputfiles=[1,0,0], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 600, model_min =0):
        useflags=use_flags([False]*60,[False]*60,[False]*70,False,False)
        useflags.update_use_K(planet,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat, fortran_kill=fortran_kill)
        return

    def minimize_one_param_P(self,planet,minimize_loglik=True, fileinput=False, filename='Kep_input', outputfiles=[1,0,0], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 600, model_min =0):
        useflags=use_flags([False]*60,[False]*60,[False]*70,False,False)
        useflags.update_use_P(planet,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat, fortran_kill=fortran_kill)
        return

    def minimize_one_param_e(self,planet,minimize_loglik=True, fileinput=False, filename='Kep_input', outputfiles=[1,0,0], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 600, model_min =0):
        useflags=use_flags([False]*60,[False]*60,[False]*70,False,False)
        useflags.update_use_e(planet,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat, fortran_kill=fortran_kill)
        return

    def minimize_one_param_w(self,planet,minimize_loglik=True, fileinput=False, filename='Kep_input', outputfiles=[1,0,0], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 600, model_min =0):
        useflags=use_flags([False]*60,[False]*60,[False]*70,False,False)
        useflags.update_use_w(planet,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat, fortran_kill=fortran_kill)
        return

    def minimize_one_param_M0(self,planet,minimize_loglik=True, fileinput=False, filename='Kep_input', outputfiles=[1,0,0], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 600, model_min =0):
        useflags=use_flags([False]*60,[False]*60,[False]*70,False,False)
        useflags.update_use_M0(planet,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat, fortran_kill=fortran_kill)
        return

    def minimize_one_param_inclination(self,planet,minimize_loglik=True, fileinput=False, filename='Kep_input', outputfiles=[1,0,0], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 600, model_min =0):
        useflags=use_flags([False]*60,[False]*60,[False]*70,False,False)
        useflags.update_use_inclination(planet,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat, fortran_kill=fortran_kill)
        return

    def minimize_one_param_lineofnodes(self,planet,minimize_loglik=True, fileinput=False, filename='Kep_input', outputfiles=[1,0,0], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 600, model_min =0):
        useflags=use_flags([False]*60,[False]*60,[False]*70,False,False)
        useflags.update_use_lineofnodes(planet,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat, fortran_kill=fortran_kill)
        return

    def minimize_one_param_offset(self,dataset,minimize_loglik=True, fileinput=False, filename='Kep_input', outputfiles=[1,0,0], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 600, model_min =0):
        useflags=use_flags([False]*60,[False]*60,[False]*70,False,False)
        useflags.update_use_offset(dataset,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat, fortran_kill=fortran_kill)
        return

    def minimize_one_param_jitter(self,dataset,minimize_loglik=True, fileinput=False, filename='Kep_input', outputfiles=[1,0,0], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 600, model_min =0):
        useflags=use_flags([False]*60,[False]*60,[False]*70,False,False)
        useflags.update_use_jitter(dataset,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat, fortran_kill=fortran_kill)
        return

    def minimize_one_param_linear_trend(self,minimize_loglik=True, fileinput=False, filename='Kep_input', outputfiles=[1,0,0], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 600, model_min =0):
        useflags=use_flags([False]*60,[False]*60,[False]*70,False,False)
        useflags.update_use_linear_trend(True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat, fortran_kill=fortran_kill)
        return


#### prep. #########
    def prepare_for_mcmc(self, rtg=[True,False,False,False], customdatasetlabels=[]):
        '''Prepare bounds and par array needed for the MCMC'''

        preparingwarnings=Warning_log([],'Preparing for MCMC')
        # put together bounds for K,P,e,w,M0, so they are in an order we are used to

        rtg=self.rtg

        par = []
        flag = []
        par_str = []
        bounds = []
        prior_nr = []
        prior_jeff = []

        for i in range(len(self.rv_data_sets)):
            if len(self.rv_data_sets[i]) == 0:
                continue
            par.append(self.rvoff[i]) #
            par_str.append(self.rvoff_str[i])
            bounds.append(self.rvoff_bounds[i])
            prior_nr.append(self.rvoff_norm_pr[i])
            prior_jeff.append(self.rvoff_jeff_pr[i])

            if rtg == [False,False,True,True]:
                flag.append(False) #
            elif rtg == [False,False,True,False]:
                flag.append(False) #
            elif rtg == [False,False,False,False]:
                flag.append(False) #
            else:
                flag.append(self.rvoff_use[i]) #


        for i in range(len(self.rv_data_sets)):
            if len(self.rv_data_sets[i]) == 0:
                continue
            par.append(self.jitt[i]) #
            par_str.append(self.jitt_str[i]) #
            bounds.append(self.jitt_bounds[i])
            prior_nr.append(self.jitt_norm_pr[i])
            prior_jeff.append(self.jitt_jeff_pr[i])

            if rtg == [False,False,True,True]:
                flag.append(False) #
            elif rtg == [False,False,True,False]:
                flag.append(False) #
            elif rtg == [False,False,False,False]:
                flag.append(False) #
            else:
                flag.append(self.jitt_use[i])


       # for i  in range(self.npl):
        for i in range(9):
            if not bool(self.use_planet[i]):
                continue

 
            par.append(self.K[i]) #
            par.append(self.P[i]) #

            if self.hkl == True:
                par.append(self.e_sinw[i]) #
                par.append(self.e_cosw[i]) #
                par.append(self.lamb[i]) #
            else:
                par.append(self.e[i]) #
                par.append(self.w[i]) #
                par.append(self.M0[i]) #

            par.append(self.i[i]) #
            par.append(self.Node[i]) #

            if rtg == [False,False,True,True]:
                flag.append(False) #
            elif rtg == [False,False,True,False]:
                flag.append(False) #
            else:
                flag.append(self.K_use[i])                

            flag.append(self.P_use[i]) #

            if self.hkl == True:
                flag.append(self.e_sinw_use[i]) #
                flag.append(self.e_cosw_use[i]) #
            else:
                flag.append(self.e_use[i]) #
                flag.append(self.w_use[i]) #

            if rtg == [False,False,True,True]:
                flag.append(False) #
            elif rtg == [False,False,True,False]:
                flag.append(False) #
            else:
                if self.hkl == True:
                    flag.append(self.lamb_use[i])
                else:
                    flag.append(self.M0_use[i])

            flag.append(self.i_use[i])

            if rtg == [False,False,True,True]:
                flag.append(False) #
            elif rtg == [False,False,True,False]:
                flag.append(False) #
            else:
                flag.append(self.Node_use[i])

            par_str.append(self.K_str[i])
            par_str.append(self.P_str[i])

            bounds.append(self.K_bound[i])
            bounds.append(self.P_bound[i])

            prior_nr.append(self.K_norm_pr[i])
            prior_nr.append(self.P_norm_pr[i])

            prior_jeff.append(self.K_jeff_pr[i])
            prior_jeff.append(self.P_jeff_pr[i])


            if self.hkl == True:

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

            else:
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

            par_str.append(self.i_str[i])
            par_str.append(self.Node_str[i] )

            bounds.append(self.i_bound[i])
            bounds.append(self.Node_bound[i])

            prior_nr.append(self.i_norm_pr[i])
            prior_nr.append(self.Node_norm_pr[i])

            prior_jeff.append(self.i_jeff_pr[i])
            prior_jeff.append(self.Node_jeff_pr[i])



        par.append(self.rv_lintr)
        #flag.append(self.rv_lintr_use)
        par_str.append(self.rv_lintr_str[0])
        bounds.append(self.rv_lintr_bounds[0])
        prior_nr.append(self.rv_lintr_norm_pr[0])
        prior_jeff.append(self.rv_lintr_jeff_pr[0])


        par.append(self.rv_quadtr)
        #flag.append(self.rv_quadtr_use)
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
            flag.append(self.rv_lintr_use) #
            flag.append(self.rv_quadtr_use)


        if rtg[1] == True:# and self.type_fit['RV'] == True:
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

            elif self.gp_kernel == 'Matern32':
                for i in range(3):
                    par.append(self.GP_mat_params[i])
                    flag.append(self.GP_mat_use[i])
                    par_str.append(self.GP_mat_str[i])
                    bounds.append(self.GP_mat_bounds[i])
                    prior_nr.append(self.GP_mat_norm_pr[i])
                    prior_jeff.append(self.GP_mat_jeff_pr[i])
                    
            elif self.gp_kernel == 'RealTerm':
                for i in range(2):
                    par.append(self.GP_drw_params[i])
                    flag.append(self.GP_drw_use[i])
                    par_str.append(self.GP_drw_str[i])
                    bounds.append(self.GP_drw_bounds[i])
                    prior_nr.append(self.GP_drw_norm_pr[i])
                    prior_jeff.append(self.GP_drw_jeff_pr[i])
 

            elif self.gp_kernel == 'dSHOKernel':
                for i in range(6):
                    par.append(self.GP_double_sho_params[i])
                    flag.append(self.GP_double_sho_use[i])
                    par_str.append(self.GP_double_sho_str[i])
                    bounds.append(self.GP_double_sho_bounds[i])
                    prior_nr.append(self.GP_double_sho_norm_pr[i])
                    prior_jeff.append(self.GP_double_sho_jeff_pr[i])



        #for i  in range(self.npl):

        for i in range(9):
            if not bool(self.use_planet[i]):
                continue

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

            if rtg[2] == False:
                flag.append(False) #
                flag.append(False) #
                flag.append(False)
            else:
                flag.append(self.t0_use[i])
                flag.append(self.pl_a_use[i])
                flag.append(self.pl_rad_use[i])


       # if rtg[3] == True: # and self.type_fit['Transit'] == True:


        for i in range(len(self.tra_data_sets)):
            if len(self.tra_data_sets[i]) != 0:

                par.append(self.tra_off[i]) #
                par_str.append(self.tra_off_str[i])
                bounds.append(self.tra_off_bounds[i])
                prior_nr.append(self.tra_off_norm_pr[i])
                prior_jeff.append(self.tra_off_jeff_pr[i])

                if rtg[2] == True:
                    flag.append(self.tra_off_use[i]) #
                else:
                    flag.append(False) #                    



        for i in range(len(self.tra_data_sets)):
            if len(self.tra_data_sets[i]) != 0:
                par.append(self.tra_jitt[i]) #
                par_str.append(self.tra_jitt_str[i]) #
                bounds.append(self.tra_jitt_bounds[i])
                prior_nr.append(self.tra_jitt_norm_pr[i])
                prior_jeff.append(self.tra_jitt_jeff_pr[i])


                if rtg[2] == True:
                    flag.append(self.tra_jitt_use[i])
                else:
                    flag.append(False) #



        if rtg[3] == True:
            if self.tra_gp_kernel == 'RotKernel':
                for i in range(4):
                    par.append(self.tra_GP_rot_params[i])
                    par_str.append(self.tra_GP_rot_str[i])
                    bounds.append(self.tra_GP_rot_bounds[i])
                    prior_nr.append(self.tra_GP_rot_norm_pr[i])
                    prior_jeff.append(self.tra_GP_rot_jeff_pr[i])
                    
                    if self.type_fit['Transit'] == True:
                        flag.append(self.tra_GP_rot_use[i])  
                    else:
                        flag.append(False)                  

            elif self.tra_gp_kernel == 'SHOKernel':
                for i in range(3):
                    par.append(self.tra_GP_sho_params[i])
                    par_str.append(self.tra_GP_sho_str[i])
                    bounds.append(self.tra_GP_sho_bounds[i])
                    prior_nr.append(self.tra_GP_sho_norm_pr[i])
                    prior_jeff.append(self.tra_GP_sho_jeff_pr[i])
                    
                    if self.type_fit['Transit'] == True:
                        flag.append(self.tra_GP_sho_use[i])
                    else:
                        flag.append(False)     

            elif self.tra_gp_kernel == 'Matern32':
                for i in range(3):
                    par.append(self.tra_GP_mat_params[i])
                    par_str.append(self.tra_GP_mat_str[i])
                    bounds.append(self.tra_GP_mat_bounds[i])
                    prior_nr.append(self.tra_GP_mat_norm_pr[i])
                    prior_jeff.append(self.tra_GP_mat_jeff_pr[i])

                    if self.type_fit['Transit'] == True:
                        flag.append(self.tra_GP_mat_use[i])
                    else:
                        flag.append(False)   


            elif self.tra_gp_kernel == 'RealTerm':
                for i in range(2):
                    par.append(self.tra_GP_drw_params[i])
                    par_str.append(self.tra_GP_drw_str[i])
                    bounds.append(self.tra_GP_drw_bounds[i])
                    prior_nr.append(self.tra_GP_drw_norm_pr[i])
                    prior_jeff.append(self.tra_GP_drw_jeff_pr[i])

                    if self.type_fit['Transit'] == True:
                        flag.append(self.tra_GP_drw_use[i])
                    else:
                        flag.append(False) 

            elif self.tra_gp_kernel == 'dSHOKernel':
                for i in range(6):
                    par.append(self.tra_GP_double_sho_params[i])
                    par_str.append(self.tra_GP_double_sho_str[i])
                    bounds.append(self.tra_GP_double_sho_bounds[i])
                    prior_nr.append(self.tra_GP_double_sho_norm_pr[i])
                    prior_jeff.append(self.tra_GP_double_sho_jeff_pr[i])

                    if self.type_fit['Transit'] == True:
                        flag.append(self.tra_GP_double_sho_use[i])
                    else:
                        flag.append(False) 


        for i in range(len(self.tra_data_sets)):
            if len(self.tra_data_sets[i]) != 0:
                par.append(self.tra_lintr[i]) #
                par_str.append(self.tra_lintr_str[i]) #
                bounds.append(self.tra_lintr_bounds[i])
                prior_nr.append(self.tra_lintr_norm_pr[i])
                prior_jeff.append(self.tra_lintr_jeff_pr[i])


                if rtg[2] == True:
                    flag.append(self.tra_lintr_use[i])
                else:
                    flag.append(False) #


        for i in range(len(self.tra_data_sets)):
            if len(self.tra_data_sets[i]) != 0:
                par.append(self.tra_quadtr[i]) #
                par_str.append(self.tra_quadtr_str[i]) #
                bounds.append(self.tra_quadtr_bounds[i])
                prior_nr.append(self.tra_quadtr_norm_pr[i])
                prior_jeff.append(self.tra_quadtr_jeff_pr[i])

 
                if rtg[2] == True:
                    flag.append(self.tra_quadtr_use[i])
                else:
                    flag.append(False) #

        for i in range(9):
            if not bool(self.use_planet[i]):
                continue

            par.append(self.omega_dot[i]) #
            flag.append(self.omega_dot_use[i])
            par_str.append(self.omega_dot_str[i]) #
            bounds.append(self.omega_dot_bounds[i])
            prior_nr.append(self.omega_dot_norm_pr[i])
            prior_jeff.append(self.omega_dot_jeff_pr[i])


        #print(len(self.ld_gr))
        for i in range(len(self.tra_data_sets)):
            if len(self.tra_data_sets[i]) == 0 or self.ld_gr[i] != i:
                continue
          #  elif self.npl==0:
          #      continue
            else:
                if self.ld_m[i] == "linear":
                    par.append(self.ld_u_lin[i][0])
                    par_str.append(self.ld_u_lin_str[i][0])
                    bounds.append(self.ld_u_lin_bound[i][0])
                    prior_nr.append(self.ld_u_lin_norm_pr[i][0])
                    prior_jeff.append(self.ld_u_lin_jeff_pr[i][0])
                    
                    if rtg[2] == False:
                        flag.append(False) #
                    else:
                        flag.append(self.ld_u_lin_use[i][0])

                    
                elif self.ld_m[i] ==  "quadratic":
                    for x in range(2):
                        par.append(self.ld_u_quad[i][x])
                        par_str.append(self.ld_u_quad_str[i][x])
                        bounds.append(self.ld_u_quad_bound[i][x])
                        prior_nr.append(self.ld_u_quad_norm_pr[i][x])
                        prior_jeff.append(self.ld_u_quad_jeff_pr[i][x])

                        if rtg[2] == False:
                            flag.append(False) #
                        else:
                            flag.append(self.ld_u_quad_use[i][x])

                elif self.ld_m[i] ==  "nonlinear":
                    for x in range(4):
                        par.append(self.ld_u_nonlin[i][x])
                        par_str.append(self.ld_u_nonlin_str[i][x])
                        bounds.append(self.ld_u_nonlin_bound[i][x])
                        prior_nr.append(self.ld_u_nonlin_norm_pr[i][x])
                        prior_jeff.append(self.ld_u_nonlin_jeff_pr[i][x])

                        if rtg[2] == False:
                            flag.append(False) #
                        else:
                            flag.append(self.ld_u_nonlin_use[i][x])
                            
                            


        if self.get_TTVs[0] == True:
           # for i in range(self.npl):

            for i in range(9):
                if not bool(self.use_planet[i]):
                    continue
                self.get_TTVs[1][i] = [len(par),len(par)+len(self.tra_ttv[i])]                         

                for z in range(len(self.tra_ttv[i])):
                
                    par.append(self.tra_ttv[i][z]) #
                    flag.append(self.tra_ttv_use[i][z])
                    par_str.append(self.tra_ttv_str[i][z]) #
                    bounds.append(self.tra_ttv_bounds[i][z])
                    prior_nr.append(self.tra_ttv_norm_pr[i][z])
                    prior_jeff.append(self.tra_ttv_jeff_pr[i][z])   
               
       #### Astr. & Stellar parameters #####
       

        par.append(self.ast_alpha[0])
        par_str.append(self.ast_alpha_str[0])
        bounds.append(self.ast_alpha_bound[0])
        prior_nr.append(self.ast_alpha_norm_pr[0])
        prior_jeff.append(self.ast_alpha_jeff_pr[0])
        
        par.append(self.ast_delta[0])
        par_str.append(self.ast_delta_str[0])
        bounds.append(self.ast_delta_bound[0])
        prior_nr.append(self.ast_delta_norm_pr[0])
        prior_jeff.append(self.ast_delta_jeff_pr[0])
        
        
        par.append(self.ast_pi[0])
        par_str.append(self.ast_pi_str[0])
        bounds.append(self.ast_pi_bound[0])
        prior_nr.append(self.ast_pi_norm_pr[0])
        prior_jeff.append(self.ast_pi_jeff_pr[0])        

        par.append(self.ast_mu_alpha[0])
        par_str.append(self.ast_mu_alpha_str[0])
        bounds.append(self.ast_mu_alpha_bound[0])
        prior_nr.append(self.ast_mu_alpha_norm_pr[0])
        prior_jeff.append(self.ast_mu_alpha_jeff_pr[0])
        
        par.append(self.ast_mu_delta[0])
        par_str.append(self.ast_mu_delta_str[0])
        bounds.append(self.ast_mu_delta_bound[0])
        prior_nr.append(self.ast_mu_delta_norm_pr[0])
        prior_jeff.append(self.ast_mu_delta_jeff_pr[0])
        
                      

        if self.type_fit["AST"] and self.use_ast_hipp_gaia == True: 
     
            flag.append(self.ast_alpha_use[0])   
            flag.append(self.ast_delta_use[0])
            flag.append(self.ast_pi_use[0])        
            flag.append(self.ast_mu_alpha_use[0])        
            flag.append(self.ast_mu_delta_use[0])
        else:
            #flag.append(False)   
            [flag.append(False) for _ in range(5)]        
         
         
         
                
        par.append(self.params.stellar_mass)
        flag.append(self.use.use_stellar_mass)
        par_str.append(self.st_mass_str[0])
        bounds.append(self.st_mass_bounds[0])
        prior_nr.append(self.st_mass_norm_pr[0])
        prior_jeff.append(self.st_mass_jeff_pr[0])
        
#        print(par)
#        print(flag)
#        print(par_str)
#        print(bounds)


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
       # #print(self.b_for_mcmc)
        self.par_for_mcmc = np.array(self.par_for_mcmc)
        self.f_for_mcmc = np.array(self.f_for_mcmc)
        self.e_for_mcmc = np.array(self.e_for_mcmc)
        self.b_for_mcmc = np.array(self.b_for_mcmc)
        self.nr_pr_for_mcmc = np.array(self.nr_pr_for_mcmc)
        self.jeff_pr_for_mcmc = np.array(self.jeff_pr_for_mcmc)

        preparingwarnings.print_warning_log()
        return




    def generate_newparams_for_mcmc(self,p):
        newparams=self.params # we will now modify this object to represent values of new parameters
        # now update these parameters which should be updated (use flag True, their indices are indicated by f_for_mcmc array)
        i=0
        for idx in self.f_for_mcmc:
            if (idx<self.ndset):
                newparams.update_offset(idx,p[i])
                i=i+1
            elif (idx<2*self.ndset):
                newparams.update_jitter(idx-self.ndset,p[i])
                i=i+1
            elif (idx<2*self.ndset+7*self.npl):
                nr=idx-2*self.ndset
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
            elif (idx<2*self.ndset+7*self.npl+1):
                newparams.update_linear_trend(p[i])
                i=i+1
            elif (idx<2*self.ndset+7*self.npl+1+self.params.GP_params.npar):
                newparams.update_GP_param_value(idx-2*self.ndset-7*self.npl-1,p[i])
                i=i+1


            #else:
            #    newparams.update_stellar_mass(p[i])
        return newparams



    def update_with_mcmc_errors(self,p):

        '''Substitute normal errors with mcmc errors, where + and - errors can be different'''

 
       # print(p)
        i=0
        for idx in self.f_for_mcmc:
            if (idx<self.ndset):
                self.param_errors.update_offset_error(idx,p[i])
                i=i+1
            elif (idx<2*self.ndset):
                self.param_errors.update_jitter_error(idx-self.ndset,p[i])
                i=i+1
            elif (idx<2*self.ndset+7*self.npl):
                nr=idx-2*self.ndset
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
            elif (idx<2*self.ndset+7*self.npl+1):
                self.param_errors.update_linear_trend_error(p[i])
                i=i+1
            elif (idx<2*self.ndset+7*self.npl+2+self.params.GP_params.npar):
                self.param_errors.update_GP_param_errors(idx-2*self.ndset-7*self.npl-2,p[i])
                #print(idx-2*self.filelist.ndset-7*self.npl-2,p[i])  #TBD!!!

                i=i+1
            else:
                self.param_errors.update_stellar_mass_error(p[i])
        return


 

 



def mcmc_stab(obj,samp_i):

    fit2 = dill.copy(obj)

    ncpus = 40
    random_dir_str = 6
    
    timemax=10000 
    timestep=10
    timeout_sec=10.0
    integrator='mvs'
    
    
    a_threshold = 0.2 # i.e. 10 % in semi-major axis
    e_threshold = 0.95 # i.e. above e>0.95 the system is defined as unstable
    
    
    rand_samp = 1000 # if None then all the samples are integrated
    
    cwd = os.getcwd()
    
    pl1_ind = 0
    pl2_ind = 1
    Per_2 = 0
    Per_1 = 5
    
    newparams = fit2.generate_newparams_for_mcmc(samp_i)
    #fit2.overwrite_params(newparams)
 
    random_dir = randomString(random_dir_str)

  
    run_stability(fit2, timemax=timemax, timestep=timestep, timeout_sec=timeout_sec,  integrator=integrator,stab_save_dir=random_dir)
 

    last_stab = []
    if np.argwhere(np.isnan(fit2.evol_M[pl1_ind])).size > 0:
        last_stab.append(min(min(np.argwhere(np.isnan(fit2.evol_M[pl1_ind])))))
    if np.argwhere(np.isnan(fit2.evol_M[pl2_ind])).size > 0:
        last_stab.append(min(min(np.argwhere(np.isnan(fit2.evol_M[pl2_ind])))))
 
    
    last_stab.append(min(np.argwhere( fit2.evol_a[pl1_ind] > fit2.evol_a[pl1_ind] +fit2.evol_a[pl1_ind]*a_threshold ),default=len(fit2.evol_a[pl1_ind])))            
    last_stab.append(min(np.argwhere( fit2.evol_a[pl2_ind] > fit2.evol_a[pl2_ind] +fit2.evol_a[pl2_ind]*a_threshold ),default=len(fit2.evol_a[pl2_ind])) )    
    last_stab.append(min(np.argwhere( fit2.evol_a[pl1_ind] < fit2.evol_a[pl1_ind] -fit2.evol_a[pl1_ind]*a_threshold ),default=len(fit2.evol_a[pl1_ind])) )           
    last_stab.append(min(np.argwhere( fit2.evol_a[pl2_ind] < fit2.evol_a[pl2_ind] -fit2.evol_a[pl2_ind]*a_threshold ),default=len(fit2.evol_a[pl2_ind])) )  
    
    last_stab.append(min(np.argwhere( fit2.evol_e[pl1_ind] > e_threshold ),default=len(fit2.evol_e[pl1_ind])) )
    last_stab.append(min(np.argwhere( fit2.evol_e[pl2_ind] > e_threshold ),default=len(fit2.evol_e[pl2_ind])) )
 

    last_stable = int(min(last_stab)-1)


    stable = fit2.evol_T[pl1_ind][last_stable]
 

    
    Prat = fit2.evol_Per[pl2_ind][0:last_stable] / fit2.evol_Per[pl1_ind][0:last_stable]        

    lambda1  = (fit2.evol_M[pl1_ind][0:last_stable]    + fit2.evol_p[pl1_ind][0:last_stable]   + 0)%360
    lambda2  = (fit2.evol_M[pl2_ind][0:last_stable]    + fit2.evol_p[pl2_ind][0:last_stable]   + 0)%360

 
 
    mean_e1 = np.mean(fit2.evol_e[pl1_ind][0:last_stable])
    min_e1 = min(fit2.evol_e[pl1_ind][0:last_stable], default=0)
    max_e1 = max(fit2.evol_e[pl1_ind][0:last_stable], default=0.99)
    ampl_e1 = (max_e1 - min_e1) / 2.0

    mean_e2 = np.mean(fit2.evol_e[pl2_ind][0:last_stable])
    min_e2 = min(fit2.evol_e[pl2_ind][0:last_stable], default=0)
    max_e2 = max(fit2.evol_e[pl2_ind][0:last_stable], default=0.9)
    ampl_e2 = (max_e2 - min_e2) / 2.0
    
    theta      = {k: [ ] for k in range(10)}  
    theta_ampl = {k: [ ] for k in range(10)}  
    theta_mean = {k: [ ] for k in range(10)}  
    
    
    coef1 = Per_2 +1
    coef2 = Per_1 +1
    order = abs(coef2 - coef1)
 
    for i in range(order+1):
 
        theta[i] = (coef1*lambda1%360 - coef2*lambda2%360 )%360 + (((coef2 -coef1) -i)*fit2.evol_p[pl1_ind][0:last_stable] + i*fit2.evol_p[pl2_ind][0:last_stable])%360 
        theta[i] = theta[i]%360

        if 90 > circ_mean_np(theta[i] ,azimuth=True) or circ_mean_np(theta[i] ,azimuth=True) > 270:
            theta[i][theta[i] >=180.0] -= 360.0
            theta_mean[i] = circ_mean_np(theta[i],azimuth=False)
        else:
            theta_mean[i] = circ_mean_np(theta[i],azimuth=True)        

        print(max(theta[i]),min(theta[i]))
        theta_ampl[i] = (max(theta[i], default=360) - min(theta[i], default=0))/2.0
#        theta_circ_mean = circ_mean_np(theta[i])


    dom  = (fit2.evol_p[pl2_ind][0:last_stable] - fit2.evol_p[pl1_ind][0:last_stable])%360
    
    
    print(circ_mean_np(dom,azimuth=True))
    if 90 > circ_mean_np(dom,azimuth=True) or circ_mean_np(dom,azimuth=True) > 270:
        dom[dom>=180.0] -= 360.0
        mean_dom = circ_mean_np(dom,azimuth=False)
    else:
        mean_dom = circ_mean_np(dom,azimuth=True)
        
 

    dom_amp = (max(dom, default=360) - min(dom, default=0))/2.0


 
    
    del fit2 
 



    return [
    np.mean(np.array(Prat)), 
    mean_e1, ampl_e1,
    mean_e2, ampl_e2,
    mean_dom,
    dom_amp, 
    theta_mean[0],
    theta_ampl[0],
    theta_mean[1],
    theta_ampl[1],
    stable]


#################### Junk here ##################################
  

