#!/usr/bin/python

from __future__ import print_function
__author__ = 'Trifon Trifonov, Jakub Morawski'

import sys, os
#sys.path.insert(0, '../lib')
sys.path.append('./lib/RV_mod/')
import jac2astrocen
import gls as gls 
import prior_functions as pr


#import prior_functions as pr
#import rot_kernels
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('SVG') 
 

import time
import multiprocessing
from scipy.signal import argrelextrema

#from emcee.utils import MPIPool
import corner
import celerite 
from celerite import terms

import batman

#import copy
import dill
import scipy.optimize as op

from CustomSampler import CustomSampler
from Warning_log import Warning_log
from parameters import parameters, GP_parameters, parameter_errors, parameter_bounds, use_flags
from fortran_output import fortran_output
from functions import *
from errors import Error, InputError, FittingError
from kernel import kernel, rvmodel, summary


TAU=6.2831853071
DEFAULT_STELLAR_MASS=1.0
DEFAULT_JITTER=1.0
NPLMAX=20
NDSETMAX=20
DEFAULT_PATH='./datafiles/'


###################### Fit transits (work in progress) #######################



def run_SciPyOp_transit(obj):      
 
    
    nll = lambda *args: -compute_loglik_transit(*args)

    flag_ind = [idx for idx in range(len(obj.tr_params_use)) if obj.tr_params_use[idx] == True ]
    
    p = []  #'p" are the fitted parameters
    b = []
    e = []
    
    for j in range(len(obj.tr_par)):
        if obj.tr_params_use[j]:
            p.append(obj.tr_par[j])
            b.append(obj.tr_bounds[j])
            e.append(obj.tr_el_str[j])
    
    b = np.array(b)
    
    
    minimzers = ['Nelder-Mead','Powell','CG','BFGS','Newton-CG','L-BFGS-B', 'TNC','COBYLA','SLSQP','dogleg','trust-ncg']
 
    #----------------- one more time using the Simplex method ---------------------#
    xtol = 1e-3
    for k in range(3): # run at least 3 times the minimizer
        xtol = xtol/10.0 
        if len(p) ==0:
            #print("Transit fitting not possible: All parameters are fixed! ")
            loglik_tra = compute_loglik_transit(p,obj,flag_ind,b,e)
            return
        else:
     
            result = op.minimize(nll, p, args=(obj,flag_ind,b,e), method=minimzers[0], options={'xtol': xtol, 'disp': True, 'maxiter':30000, 'maxfev':30000 })
            p = result["x"]

    print("Best fit par.:")  
 
    for j in range(len(p)):
        print(e[j] + "  =  %s"%p[j])
   
 
    
def compute_loglik_transit(p,copied_obj,flag_ind,b,e):
    #newparams=copied_obj.generate_newparams_for_mcmc(p)
  #  oldparams=signalfit.params
              
 
    for j in range(len(p)):
 
        if p[j] <= b[j,0] or p[j] >= b[j,1]:
            return -np.inf
        else:
            copied_obj.tr_par[flag_ind[j]] = p[j]  
            
    #print(copied_obj.tr_par)
    if len(copied_obj.tra_data_sets[0]) == 0:
        return -np.inf        
    else: 
        t = copied_obj.tra_data_sets[0][0] 
        flux = copied_obj.tra_data_sets[0][1] 
        flux_err = copied_obj.tra_data_sets[0][2] 
        
        
        #print(copied_obj.tr_par[0],copied_obj.tr_par[1])
        copied_obj.tr_params.t0  = copied_obj.tr_par[0] #0.0  #time of inferior conjunction
        copied_obj.tr_params.per = copied_obj.tr_par[1] #1.0    #orbital period
        copied_obj.tr_params.ecc = copied_obj.tr_par[2] #0.0  
        copied_obj.tr_params.w   = copied_obj.tr_par[3] #90.0   #longitude of periastron (in degrees)               
        copied_obj.tr_params.rp  = copied_obj.tr_par[4] #0.15   #planet radius (in units of stellar radii)
        copied_obj.tr_params.inc = copied_obj.tr_par[5] #90. #orbital inclination (in degrees)
        copied_obj.tr_params.a   = copied_obj.tr_par[6] #15  #semi-major axis (in units of stellar radii)

       
        m = batman.TransitModel(copied_obj.tr_params, t)    #initializes model
 
        flux_model = m.light_curve(copied_obj.tr_params)          #calculates light curve  
        tr_o_c = flux -flux_model
        S = 0
        for i in range(len(tr_o_c)):
            S= S + (((tr_o_c[i]**2)/(2*flux_err[i]**2) ) - 0.5*(np.log(TAU*(flux_err[i]**2))))
 
        return -S
  
 
###################### Fit transits END (work in progress) #######################
    
    
    
def run_SciPyOp(obj,  Kbounds=[[0.0,100000.0]],Pbounds=[[0.0,100000.0]],ebounds=[[-0.99,0.99]],wbounds=[[-2.0*360.0, 2.0*360.0]],Mbounds=[[-2.0*360.0, 2.0*360.0]],ibounds=[[-2.0*180.0, 2.0*180.0]],capbounds=[[-2.0*360.0, 2.0*360.0]],offbounds=[[-100000.0,100000.0]],jitbounds=[[0.0,10000.0]],lintrbounds=[[-10.0,10.0]], GPbounds=[[0.0,100000.0]], stmassbounds=[[0.01,1000.0]], prior=0, samplesfile='', level=(100.0-68.3)/2.0, threads=1, doGP=False, gp_par=None, kernel_id=-1, use_gp_par=[False,False,False,False], save_means=False, fileoutput=False, save_sampler=False,burning_ph=20, mcmc_ph=20, **kwargs):      
 
 
    if (doGP):
        obj.initiategps(gp_par=gp_par, use_gp_par=use_gp_par, kernel_id=kernel_id)
        nll = lambda *args: -compute_loglik_SciPyOp2(*args)
    else:
        nll = lambda *args: -compute_loglik_SciPyOp(*args)
        GPbounds=[[x-10.0,x+10.0] for x in obj.params.GP_params.gp_par] # just to make sure GPbounds don't cause lnprob return -infinity when we don't do GP (now all GP_params will be within bounds for sure)
     
 
    obj.prepare_for_mcmc(Kbounds=Kbounds,Pbounds=Pbounds,ebounds=ebounds,wbounds=wbounds,Mbounds=Mbounds,ibounds=ibounds,capbounds=capbounds,offbounds=offbounds,jitbounds=jitbounds,lintrbounds=lintrbounds, GPbounds=GPbounds, stmassbounds=stmassbounds)    
    pp = obj.par_for_mcmc
 
    
   # b = np.array(obj.bounds.offset_bounds,obj.bounds.jitter_bounds,obj.bounds.planet_params_bounds,
   #                    obj.bounds.linear_trend_bounds,obj.bounds.GP_params_bounds,obj.bounds.stellar_mass_bounds)
   # b.flatten()
   # print(b)
   # bounds=b,
    minimzers = ['Nelder-Mead','Powell','CG','BFGS','Newton-CG','L-BFGS-B', 'TNC','COBYLA','SLSQP','dogleg','trust-ncg']

   # for k in range(2): # run at least 3 times the minimizer
   #     result = op.minimize(nll, pp, args=(obj), method=minimzers[6], bounds=None, options={'xtol': 1e-6, 'disp': True })
   #     pp = result["x"]
        
  #  print("Best fit par.:", result["x"])
#----------------- one more time using the Simplex method ---------------------#
    xtol = 1e-3
    for k in range(3): # run at least 3 times the minimizer
        xtol = xtol/10.0 
        print(k,xtol)
     
        result = op.minimize(nll, pp, args=(obj), method=minimzers[6], options={'xtol': xtol, 'disp': True, 'maxiter':30000, 'maxfev':30000 })
        pp = result["x"]

    print("Best fit par.:", result["x"])
    
    obj.par_for_mcmc = result["x"]
       

def compute_loglik_SciPyOp(p,copied_obj):
    newparams=copied_obj.generate_newparams_for_mcmc(p)
  #  oldparams=signalfit.params
    copied_obj.overwrite_params(newparams)
    if not (copied_obj.verify_params_with_bounds()):
        return -np.inf # if parameters are outside of bounds we return -infinity
    else:               
        flag=copied_obj.fitting(fileinput=False, filename='Kep_input', minimize_loglik=True, amoeba_starts=0, outputfiles=[1,0,0],return_flag=True)

    return copied_obj.loglik 
 


def compute_loglik_SciPyOp2(p,copied_obj):
    newparams=copied_obj.generate_newparams_for_mcmc(p)
  #  oldparams=signalfit.params
    copied_obj.overwrite_params(newparams)
    if not (copied_obj.verify_params_with_bounds()):
        return -np.inf # if parameters are outside of bounds we return -infinity
   # elif not (copied_obj.fitting_method.startswith('GP')):               
   #     flag=copied_obj.fitting(fileinput=False, filename='Kep_input', minimize_loglik=True, amoeba_starts=0, outputfiles=[1,0,0],return_flag=True)

   #     return copied_obj.loglik 
    
   # elif(copied_obj.fitting_method.startswith('GP')):
    else:
        flag=copied_obj.fitting(fileinput=False, filename='Kep_input', minimize_loglik=True, amoeba_starts=0, outputfiles=[1,1,0],return_flag=True)
        #newparams=copied_obj.generate_newparams_for_mcmc(p)
        #copied_obj.overwrite_params(newparams)   
        S=0
        for i in range(copied_obj.filelist.ndset):
            copied_obj.gps[i].set_parameter_vector(np.array(list(map(np.log,np.concatenate((copied_obj.params.GP_params.gp_par,np.atleast_1d(copied_obj.params.jitters[i]))))))) 
            S = S + copied_obj.gps[i].log_likelihood(copied_obj.fit_results.rv_model.o_c[copied_obj.fit_results.idset==i])
                      
       #print(copied_obj.fit_results.loglik, compute_loglik_TT(p,copied_obj),  S)#,copied_obj.params.GP_params[:4],copied_obj.params.planet_params[:14])   
        return  S  



  
def lnprob(p,copied_obj,prior):
    '''Compute logarithmic probability for given parameters'''     
    # now we need to compute loglikelihood using the fortran code on new parameters, and then add it to lp
    #copied_obj=signalfit
    newparams=copied_obj.generate_newparams_for_mcmc(p)
   # oldparams=copied_obj.params
    copied_obj.overwrite_params(newparams)
    if not (copied_obj.verify_params_with_bounds()):
        return -np.inf # if parameters are outside of bounds we return -infinity
    else:      

       # print(copied_obj.params.jitters[:3],copied_obj.params.GP_params[:4],copied_obj.params.planet_params[:14])        
       # print(copied_obj.use.use_jitters[:3])    
        flag=copied_obj.fitting(fileinput=False, filename='Kep_input', minimize_loglik=True, amoeba_starts=0, outputfiles=[0,0,0],return_flag=True)
        #print(copied_obj.fit_results.loglik, flag)         
        #print(copied_obj.fit_results.loglik, compute_loglik_TT(p,copied_obj) )#,copied_obj.params.GP_params[:4],copied_obj.params.planet_params[:14])   
        #copied_obj.overwrite_params(oldparams)
        #now self.fit_results.loglik is the loglikelihood corresponding to new parameters
        if (flag==1):
            #print(pr.choose_prior(p,prior)+copied_obj.fit_results.loglik)
            
            return pr.choose_prior(p,prior)+copied_obj.fit_results.loglik
        else:
            return -np.inf
 

def compute_loglik_TT(p,signalfit):
    #newparams=signalfit.generate_newparams_for_mcmc(p)
    #oldparams=signalfit.params
   # signalfit.overwrite_params(newparams)
    S=0
    for i in range(len(signalfit.fit_results.rv_model.o_c)):
        S=S-(signalfit.fit_results.rv_model.o_c[i]**2)/(2*(signalfit.fit_results.rv_model.rv_err[i]**2+signalfit.params.jitters[signalfit.fit_results.idset[i]]**2))-0.5*np.log(TAU*(signalfit.fit_results.rv_model.rv_err[i]**2+signalfit.params.jitters[signalfit.fit_results.idset[i]]**2))
    return S
     

      
def compute_loglik(p,signalfit):
    newparams=signalfit.generate_newparams_for_mcmc(p)
    oldparams=signalfit.params
    signalfit.overwrite_params(newparams)
    if not (signalfit.verify_params_with_bounds()):
        return -np.inf # if parameters are outside of bounds we return -infinity
    else:               
        flag=signalfit.fitting(fileinput=True, filename='Kep_input', minimize_loglik=True, amoeba_starts=0, outputfiles=[0,1,0],return_flag=True)

        #now self.fit_results.loglik is the loglikelihood corresponding to new parameters
        if not (flag==1):
            signalfit.overwrite_params(oldparams)         
            return -np.inf        
        else: 
            S=0
            for i in range(len(signalfit.fit_results.rv_model.o_c)):
                S=S-(signalfit.fit_results.rv_model.o_c[i]**2)/(2*(signalfit.fit_results.rv_model.rv_err[i]**2+signalfit.params.jitters[signalfit.fit_results.idset[i]]**2))-0.5*np.log(TAU*(signalfit.fit_results.rv_model.rv_err[i]**2+signalfit.params.jitters[signalfit.fit_results.idset[i]]**2))
            return S

def lnprobGP(p,signalfit,prior):
    '''Compute logarithmic probability for given parameters'''     
    # now we need to compute loglikelihood using the fortran code on new parameters, and then add it to lp
    newparams=signalfit.generate_newparams_for_mcmc(p)
    oldparams=signalfit.params
    signalfit.overwrite_params(newparams)
    if not (signalfit.verify_params_with_bounds()):
        return -np.inf # if parameters are outside of bounds we return -infinity
    else:               
        flag=signalfit.fitting(fileinput=True, filename='Kep_input', minimize_loglik=True, amoeba_starts=0, outputfiles=[0,1,0],return_flag=True)

        #now self.fit_results.loglik is the loglikelihood corresponding to new parameters
        if not (flag==1):
            signalfit.overwrite_params(oldparams)
            return -np.inf        
        else:
            S=0

            for i in range(signalfit.filelist.ndset):
                
                #print(len(signalfit.fit_results.rv_model.o_c[signalfit.fit_results.idset==i]),signalfit.filelist.ndset,max(signalfit.fit_results.idset))
                
                signalfit.gps[i].set_parameter_vector(np.array(list(map(np.log,np.concatenate((signalfit.params.GP_params.gp_par,np.atleast_1d(signalfit.params.jitters[i])))))))
                S+=signalfit.gps[i].log_likelihood(signalfit.fit_results.rv_model.o_c[signalfit.fit_results.idset==i])
                signalfit.overwrite_params(oldparams)
                #print(S)
 
                #print("something")
                #wait = raw_input("PRESS ENTER TO CONTINUE.")
                #print("something")
                
            return pr.choose_prior(p,prior)+S
            
def lnprobGP2(p,copied_obj,prior):
    '''Compute logarithmic probability for given parameters'''     
    # now we need to compute loglikelihood using the fortran code on new parameters, and then add it to lp
    
    newparams=copied_obj.generate_newparams_for_mcmc(p)
    copied_obj.overwrite_params(newparams)
    
    if not (copied_obj.verify_params_with_bounds()):
        return -np.inf # if parameters are outside of bounds we return -infinity
    else:  
        flag=copied_obj.fitting(fileinput=False, filename='Kep_input', minimize_loglik=True, amoeba_starts=0, outputfiles=[0,1,0],return_flag=True)
        #now self.fit_results.loglik is the loglikelihood corresponding to new parameters
        if not (flag==1):
            return -np.inf  
            
        elif copied_obj.fit_results.loglik == 0:
            return -np.inf             
                  
        else:
            S=0
            for i in range(copied_obj.filelist.ndset):

                copied_obj.gps[i].set_parameter_vectorset_parameter_vector(np.array(list(map(np.log,np.concatenate((copied_obj.params.GP_params.gp_par,np.atleast_1d(signalfit.params.jitters[i])))))))    
                S = S + copied_obj.gps[i].log_likelihood(copied_obj.fit_results.rv_model.o_c[copied_obj.fit_results.idset==i])

            return pr.choose_prior(p,prior)+S   


def lnprobGP3(p,copied_obj,prior): 
    '''Compute logarithmic probability for given parameters'''     
 
    fit_log = lnprob(p,copied_obj,prior)  
 
    if fit_log == -np.inf or fit_log == 0:
        return -np.inf # if parameters are outside of bounds we return -infinity
    else:  
       newparams=copied_obj.generate_newparams_for_mcmc(p)
       copied_obj.overwrite_params(newparams)
       
       S=0
       for i in range(copied_obj.filelist.ndset):
           copied_obj.gps[i].set_parameter_vector(np.array(list(map(np.log,np.concatenate((copied_obj.params.GP_params.gp_par,np.atleast_1d(copied_obj.params.jitters[i]))))))) 
           S = S + copied_obj.gps[i].log_likelihood(copied_obj.fit_results.rv_model.o_c[copied_obj.fit_results.idset==i])
                      
       #print(copied_obj.fit_results.loglik, compute_loglik_TT(p,copied_obj),  S)#,copied_obj.params.GP_params[:4],copied_obj.params.planet_params[:14])   
       return pr.choose_prior(p,prior)+S    
 

def run_mcmc(obj,  Kbounds=[[0.0,100000.0]],Pbounds=[[0.0,100000.0]],ebounds=[[-0.99,0.99]],wbounds=[[-2.0*360.0, 2.0*360.0]],Mbounds=[[-2.0*360.0, 2.0*360.0]],ibounds=[[-2.0*180.0, 2.0*180.0]],capbounds=[[-2.0*360.0, 2.0*360.0]],offbounds=[[-100000.0,100000.0]],jitbounds=[[0.0,10000.0]],lintrbounds=[[-10.0,10.0]], GPbounds=[[0.0,100000.0]], stmassbounds=[[0.01,1000.0]], prior=0, samplesfile='', level=(100.0-68.3)/2.0, threads=1, doGP=False, gp_par=None, gp_kernel_id=-1, use_gp_par=[True,True,True,True], save_means=False, fileoutput=False, save_sampler=False,burning_ph=10, mcmc_ph=10, **kwargs):      

    '''Performs MCMC and saves results'''  
    
    if threads == 'max':
        threads = multiprocessing.cpu_count()
 
    
    # Let's prepare things depending if we want to do GP or not

    if (doGP):
        obj.initiategps(gp_par=gp_par, use_gp_par=use_gp_par, kernel_id=gp_kernel_id)
        fun=lnprobGP3                            
    else:
        fun=lnprob                           
        GPbounds=[[x-10.0,x+10.0] for x in obj.params.GP_params.gp_par] # just to make sure GPbounds don't cause lnprob return -infinity when we don't do GP (now all GP_params will be within bounds for sure)
 
    # we will need this later
    if (obj.mod_dynamical):
        mod='dyn'
    else:
        mod='kep'

    obj.prepare_for_mcmc(Kbounds=Kbounds,Pbounds=Pbounds,ebounds=ebounds,wbounds=wbounds,Mbounds=Mbounds,ibounds=ibounds,capbounds=capbounds,offbounds=offbounds,jitbounds=jitbounds,lintrbounds=lintrbounds, GPbounds=GPbounds, stmassbounds=stmassbounds)    
 

    start_time = time.time()

    ndim, nwalkers = len(obj.par_for_mcmc), len(obj.par_for_mcmc)*4

    #par_for_mcmc_ = obj.par_for_mcmc
    fff = dill.copy(obj)
    pos = [obj.par_for_mcmc + 1e-3*np.random.rand(ndim) for i in range(nwalkers)]

    sampler = CustomSampler(nwalkers, ndim, fun, args=[fff,prior], threads = threads)
 
 
    
    # burning phase
    pos, prob, state  = sampler.run_mcmc(pos,burning_ph)

    sampler.reset()
 
    # now perform the MCMC

    pos, prob, state  = sampler.run_mcmc(pos,mcmc_ph)
    
    
              
    print("--- %s seconds ---" % (time.time() - start_time))  
 
    ln = np.hstack(sampler.lnprobability)
    sampler.save_samples(obj.f_for_mcmc,obj.filelist.ndset,obj.npl)
            
    
    if (fileoutput):
        #if (samplesfile==''): # that means no file name for samples file has been provided, so we generate a default one
       #     samplesfile='samples_%s'%mod
        #obj.mcmc_sample_file = 'mcmc_samples'+'_%s'%mod
        #obj.corner_plot_file = 'cornerplot.png'
        outfile = open(str(obj.mcmc_sample_file), 'w') # file to save samples

        for j in range(len(sampler.samples)):
            outfile.write("%s  " %(ln[j]))        
            for z in range(len(obj.par_for_mcmc)):
                outfile.write("%s  " %(sampler.samples[j,z]))
            outfile.write("\n")

        outfile.close()        
            
    # Now we will save new parameters and their errors (different + and - errors in this case). Flag save_means determines if we want to take means as new best fit parameters or stick to old ones and calculate errors with respect to that           
    if (save_means):
        obj.par_for_mcmc = sampler.means # we will not need to keep the old parameters in this attribbute, so let's store the means now
        
    new_par_errors = [[float(obj.par_for_mcmc[i] - np.percentile(sampler.samples[:,i], [level])),float(np.percentile(sampler.samples[:,i], [100.0-level])-obj.par_for_mcmc[i])] for i in range(len(obj.par_for_mcmc))] 
    newparams = obj.generate_newparams_for_mcmc(obj.par_for_mcmc)        
    #print(newparams.GP_params)
    current_GP_params=newparams.GP_params.gp_par # because calling fitting will overwrite them
   # print(current_GP_params)

    obj.overwrite_params(newparams)
    obj.fitting(minimize_loglik=True, amoeba_starts=0, outputfiles=[1,1,1]) # this will help update some things 
    obj.update_with_mcmc_errors(new_par_errors)
    obj.params.update_GP_params(current_GP_params)

    #print(current_GP_params)

    if (doGP):
        obj.fitting_method='GP_%s'%mod
    else:
        obj.fitting_method='mcmc_%s'%mod  

    if (doGP):
        loglik_to_save = lnprobGP(obj.par_for_mcmc,obj,prior)
        obj.loglik=loglik_to_save         
        obj.fit_results.loglik=loglik_to_save  
        
        
        ########### Ugly things are happening here! to be fixed! ##########
        
        fitted_GP = sampler.means[-len(obj.params.GP_params.gp_par[obj.use.use_GP_params==True]):]
 
        z = 0
        for k in range(obj.params.GP_params.npar):
            obj.param_errors.GP_params_errors[k] = [0,0]
            if not obj.use.use_GP_params[k]:
                continue
            else:
                obj.params.GP_params.gp_par[k] = fitted_GP[z]
               # obj.param_errors.GP_params_errors[k] = [float(sampler.means[-obj.params.GP_params.npar+k] - np.percentile(sampler.means[-obj.params.GP_params.npar+k], [level])),float(np.percentile(sampler.means[-obj.params.GP_params.npar+k], [100.0-level])-sampler.means[-obj.params.GP_params.npar+k])]
                obj.param_errors.GP_params_errors[k] = [float(fitted_GP[z] - np.percentile(sampler.samples[:,-len(fitted_GP)+z],[level])), float(np.percentile(sampler.samples[:,-len(fitted_GP)+z], [100.0-level])-fitted_GP[z])]
                #print(obj.param_errors.GP_params_errors[k], fitted_GP[z], np.percentile(fitted_GP[z],level),level)
            z = z+1  

        obj.GP_params=GP_parameters(len(obj.params.GP_params.gp_par),obj.params.GP_params.gp_par,kernel_id=0) # we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
         
        message_str =""" """
        if(obj.fitting_method.startswith('GP')):
            message_str = message_str +"""\nStellar activity was modeled using Gaussian Processes. The resulting GP parameters (means) are as follows:
"""
            message_str = message_str +"""\n A = {0:>7.4f} +/- {4:>7.4f}\n t = {1:>7.4f} +/- {5:>7.4f}\n P = {2:>7.4f} +/- {6:>7.4f}\n f = {3:>7.4f} +/- {7:>7.4f}
""".format(obj.GP_params.gp_par[0],obj.GP_params.gp_par[1],obj.GP_params.gp_par[2],obj.GP_params.gp_par[3],max(obj.param_errors.GP_params_errors[0]),max(obj.param_errors.GP_params_errors[1]),max(obj.param_errors.GP_params_errors[2]),max(obj.param_errors.GP_params_errors[3]))         
            message_str = message_str +"""
(The GP_params are printed here as these are still not part of the "params" structure. TBD!)                
"""
       # print(message_str)
        ###################################################################
        

    ###############  This is not working! you cannot save the sampler as an atribute and call it back later!
    ###############  See https://github.com/dfm/emcee/issues/148
    if(save_sampler):
        obj.sampler=sampler             
        obj.sampler_saved=True           
        
    sampler.reset()

    return obj
 
    
 

def cornerplot(obj, fileinput=False, level=(100.0-68.3)/2.0, **kwargs): 

    #obj = dill.copy(copied_obj)
    '''Generates a corner plot visualizing the mcmc samples. Optionally samples can be read from a file.'''
    #self.mcmc_sample_file = 'mcmc_samples'+'_%s'%mod
    #self.corner_plot_file = 'cornerplot.png'
    if(fileinput):
        samples=read_file_as_array_of_arrays_mcmc(obj.mcmc_sample_file)
   # elif(obj.sampler_saved):
   #     samples=obj.sampler.samples
    else:
        raise Exception ('Please run mcmc and save sampler or provide a valid samples file!')
    
    fig = corner.corner(samples,bins=25, color="k", reverse=True, upper= True, labels=obj.e_for_mcmc, quantiles=[level/100.0, 1.0-level/100.0],levels=(0.6827, 0.9545,0.9973), smooth=1.0, smooth1d=1.0, plot_contours= True, show_titles=True, truths=obj.par_for_mcmc, dpi = 300, pad=15, labelpad = 50 ,truth_color ='r', title_kwargs={"fontsize": 12}, scale_hist=True,  no_fill_contours=True, plot_datapoints=True, kwargs=kwargs)
    fig.savefig(obj.corner_plot_file)  

   # obj = None
    #fig.clf()
    return      
    
        
def fitting(obj, minimize_fortran=True, minimize_loglik=False, fileinput=False, doGP=False, gp_par=None, kernel_id=-1, use_gp_par=[False,False,False,False], filename='Kep_input', outputfiles=[1,1,1], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 500,model_min =0): # run the fit which will either minimize chi^2 or loglik.
    '''       
     eps, dt - accuracy and step for the integration in dynamical case, in keplerian case ignored, just given for input consistency
     which value to minimize (used for calling an appropriate fortran program)
    '''      
    
    eps = eps * 1e-13 
    dt  = dt  * 86400.0         
    
    if(minimize_loglik):
        minimized_value='loglik'
    else:
        minimized_value='chi2'
    # which mode to use, analogously as above
    if(obj.mod_dynamical):
        mod='dyn'
    else:
        mod='kep'
        
    if(minimize_fortran):   
        program='./lib/fr/%s_%s'%(minimized_value,mod) 
        text,flag=run_command_with_timeout(obj.fortran_input(program=program, fileinput=fileinput, filename=filename, outputfiles=outputfiles,amoeba_starts=amoeba_starts,eps=eps,dt=dt, when_to_kill=fortran_kill, npoints=npoints, model_max = model_max, model_min =model_min), timeout_sec, output=True,pipe=(not bool(outputfiles[2]))) # running command generated by the fortran_input function 
       
    else:
        obj.fitting_SciPyOp(doGP=doGP, gp_par=gp_par, kernel_id=kernel_id, use_gp_par=use_gp_par)         
        program='./lib/fr/%s_%s'%(minimized_value,mod) 
        text,flag=run_command_with_timeout(obj.fortran_input(program=program, fileinput=fileinput, filename=filename, outputfiles=outputfiles,amoeba_starts=0,eps=eps,dt=dt, when_to_kill=fortran_kill, npoints=npoints, model_max = model_max, model_min =model_min), timeout_sec, output=True,pipe=(not bool(outputfiles[2]))) # running command generated by the fortran_input function 
       
              
    if (flag==1):
        fortranoutput=fortran_output(text,obj.npl,obj.filelist.ndset,obj.params.stellar_mass) # create an object for fortran output, which we need to split into part
        obj.fit_results=fortranoutput.modfit(print_stat=print_stat)
        obj.stat_saved=obj.fit_results.stat_array_saved
        if (obj.stat_saved):
            obj.never_saved=False
        obj.model_saved=bool(outputfiles[1])
        obj.fit_performed=True
        if(obj.fit_results.stat_array_saved): 
            obj.fitting_method=program
        obj.update_with_fit_results()
        obj.correct_elements() #because amoeba might make things wrong here
    if (return_flag):
        return flag
    else:
        return obj



def custom_param_file_for_stability(max_time,time_step):

        ##### create the param.in file (change only the "t_max" and the "dt" for now) ######
    param_file = open('param.in', 'wb') 
        
    max_time = float(max_time)*365.2425 # make it is days
 
    param_file.write("""0.0d0 %s %s
%s %s
        
F T T T T F
0.001 50.0 50.0 -1. T
 bin.dat
unknown
"""%(max_time, time_step, max_time/1e4, max_time/1e3 ))
 
    param_file.close()
    return


def phase_planet_signal(obj,planet):

    if obj.npl ==0 or len(obj.fit_results.rv_model.jd) ==0:
        return [-1], [-1] #[[0],[0]], [[0],[0],[0],[0]]
    else:
        #copied_obj = copy.deepcopy(obj) 
         
        copied_obj = dill.copy(obj) 
   
        index = planet - 1
        ############################################      
        ######### and here is the trick!  ##########
        ############################################        
        pp0 =  copied_obj.params.planet_params[7*index+0]  # we define a variable to be the planet amplitude Kj   
        #print(pp0)
        copied_obj.params.planet_params[7*index+0] = 0 # then we set Kj to be 0, i.e. remove the j-th planet signal
        copied_obj.fitting(fileinput=False, filename='Kep_input', minimize_loglik=True, amoeba_starts=0, 
                           outputfiles=[0,1,1],return_flag=False, npoints=int(len(obj.fit_results.model)), 
                           model_max=int(max(obj.fit_results.model_jd)-max(copied_obj.fit_results.rv_model.jd)),
                           model_min=int(min(copied_obj.fit_results.rv_model.jd)-min(obj.fit_results.model_jd)))
        # and we create the static Nplanet model for the data and the model curve 
        # now this model residuals will contain ONLY the j-th planet signal + the best fit residuals
       
        copied_obj.params.planet_params[7*index+0] = pp0 # we restore Kj to its best fit value.
        ############################################      
        #########      trick is over      ##########
        ############################################  
        
    
        ############ phase fold fix for sparse model ######use_flags
        model_time_phase = np.array( (copied_obj.fit_results.model_jd -copied_obj.fit_results.rv_model.jd[0] )% copied_obj.params.planet_params[7*index+1] )
             
        sort = sorted(range(len(model_time_phase)), key=lambda k: model_time_phase[k])                        
        model_time_phase  = model_time_phase[sort]
        phased_model      = obj.fit_results.model[sort] - copied_obj.fit_results.model[sort]
    
        ############ phase data ######
        data_time_phase = np.array( (copied_obj.fit_results.rv_model.jd - copied_obj.fit_results.rv_model.jd[0])% copied_obj.params.planet_params[7*index+1] )
             
        sort = sorted(range(len(data_time_phase)), key=lambda k: data_time_phase[k])                        
        data_time_phase      = data_time_phase[sort]
        phased_data          = copied_obj.fit_results.rv_model.o_c[sort]#  - copied_obj.fit_results.rv_model.rvs[sort] 
        phased_data_err      = copied_obj.fit_results.rv_model.rv_err[sort]  
        phased_data_idset    = copied_obj.fit_results.idset[sort]  
     
        
        
        ##################### 
              
        #phased_data =  rvs[i]-signal1[i] 
     
        model = [model_time_phase,  phased_model]
        data  = [data_time_phase,  phased_data, phased_data_err, phased_data_idset]
       
        
        return data, model



def planet_orbit_xyz(obj, planet):

    u1 = obj.params.stellar_mass * (4*np.pi*np.pi)/(365.25*365.25)
    mean_orb = np.linspace(0,2.0*np.pi, 360.0)
    
    x = np.zeros(len(mean_orb))
    y = np.zeros(len(mean_orb))
    z = np.zeros(len(mean_orb))
    u = np.zeros(len(mean_orb))
    v = np.zeros(len(mean_orb))
    w = np.zeros(len(mean_orb))
    
    dist =  np.zeros(len(mean_orb))
        
    q = (1.0 - obj.params.planet_params[2 + int(planet)*7])*float(obj.fit_results.a[int(planet)])
    
    
    #this need to be fixed to work with arrays
    for f in range(len(mean_orb)):
        x[f],y[f],z[f],u[f],v[f],w[f] = jac2astrocen.mco_el2x(u1,q,
                                                       obj.params.planet_params[2 + int(planet)*7],
                                                       np.radians(obj.params.planet_params[5 + int(planet)*7]-90.0),
                                                       np.radians(obj.params.planet_params[3 + int(planet)*7]) - np.radians(obj.params.planet_params[6 + int(planet)*7]),
                                                       np.radians(obj.params.planet_params[6 + int(planet)*7] ), mean_orb[f])    
                                                       
        dist[f] =  np.sqrt(x[f]**2.0 + y[f]**2.0 + z[f]**2.0)                                        
    
    x_p,y_p,z_p,u_p,v_p,w_p = jac2astrocen.mco_el2x(u1,q,
                                                       obj.params.planet_params[2 + int(planet)*7],
                                                       np.radians(obj.params.planet_params[5 + int(planet)*7] -90.0),
                                                       np.radians(obj.params.planet_params[3 + int(planet)*7]) - np.radians(obj.params.planet_params[6 + int(planet)*7]),
                                                       np.radians(obj.params.planet_params[6 + int(planet)*7]), np.radians(obj.params.planet_params[4 + int(planet)*7]))    
 

    min_index = np.unravel_index(np.argmin(dist, axis=None), dist.shape)                                                    
    max_index = np.unravel_index(np.argmax(dist, axis=None), dist.shape)                                                                                                           
                                                       
    return np.array([x,y,z,u,v,w]), np.array([x_p,y_p,z_p,u_p,v_p,w_p]), np.array([x[min_index],y[min_index],z[min_index],u[min_index],v[min_index],w[min_index]]), np.array([x[max_index],y[max_index],z[max_index],u[max_index],v[max_index],w[max_index]])



            

# Defining some custom exceptions and warnings


# class for signal to run gls and plotting functions on (it can be both the entire signal and O-C, depends what we pass as rvs in the __ini__ function)
# This class must completely refurbished!
class signal_data(object): 

    def __init__(self,jd,rvs,rv_error, sig_for_gls=np.array([0.1,0.01,0.001])):
        self.jd=np.array(list(map(float,jd))) # for some reason sometimes we get strings instead of floats, so...
        self.rvs=np.array(list(map(float,rvs)))
        self.rv_error=np.array(list(map(float,rv_error)))
        self.sig=sig_for_gls # list of significances for gls, so the algorithm will do the computation for all three of them and include relevant levels in the z array outputed by lomb_scargle
        self.sig=np.array(self.sig)
        # It is convenient to store significance levels in decreasing order, so let's sort them
       
        sort = convert_array_to_int(np.array(sorted(range(len(self.sig)), key=lambda k: self.sig[k], reverse=True))) # find the permutation in which they are sorted in decreasing order 
        # sorting based on the permutation sort calculated above
        self.sig = self.sig[sort]             
        
  
    def gls(self, ind_sig=2,sig_for_gls=np.array([0.1,0.01,0.001])): # ind_sig will be the index of the significance level which we assume to be sufficient to declare a planet candidate

        gls_warnings=Warning_log([],'Running gls')
        if not (is_int(str(ind_sig),bounded=[True,True],bounds=[0,len(self.sig)],equal=[True,False])):
            ind_sig=len(self.sig)-1
            gls_warnings.update_warning_list('ind_sig for gls must be a non-negative integer smaller than the size of sig_for_gls array! Default value of len(sig)-1 will be assumed.')

        ### Compute the Lomb-Scargle Periodogram
        try:

            #self.gls_range = (float(max(self.jd))-float(min(self.jd)))*2 # range for gls, twice the time range     
            #self.periods = np.linspace(1, self.gls_range, 2000) # abscissas for the period range
            #omega = TAU / self.periods # converting the periods array into the frequency range for the Lomb-Scargle Periodogram evaluation
           # omega = 1 / self.periods # converting the periods array into the frequency range for the Lomb-Scargle Periodogram evaluation
             #omega = 1/ np.logspace(-0.05, 4, num=1000)
            omega = 1/ np.logspace(-0.05, 4, num=1000)        
        
            RV_gls = gls.Gls((self.jd, self.rvs, self.rv_error), fast=True,  verbose=False, norm= "ZK",ofac=5, fbeg=omega[999], fend=omega[ 0],)
        
            #self.P_G, self.z = lomb_scargle(self.jd, self.rvs, self.rv_error, omega, generalized=True, significance=self.sig) # Lomb-Scargle for the RV signal    
            self.P_G = RV_gls.power 
            self.z = RV_gls.powerLevel(sig_for_gls)
            self.periods = 1/RV_gls.freq
            self.gls_range = (float(max(self.jd))-float(min(self.jd)))*2 # range for gls, twice the time range     
            
            per_ind = argrelextrema(self.P_G, np.greater) # generates an array with the indices of all the relative extrema in the periodogram
 
 
    
            self.best_per = self.periods[per_ind] # periods corresponding to the indices calculated above
            self.best_peaks = self.P_G[per_ind] # peak heights of these extrema
            if (len(self.best_peaks)>0): # don't sort if there's no peaks    
                sort = convert_array_to_int(np.array(sorted(range(len(self.best_peaks)), key=lambda k: self.best_peaks[k], reverse=True))) # find the permutation of the peaks in which they are sorted by height in decreasing order 
                # sorting both arrays (periods and peak heights) by peak height, based on the permutation sort calculated above
                self.best_peaks = self.best_peaks[sort]
                self.best_per = self.best_per[sort]  
        
            # now we save the number of peaks which exceed the significance requirement          
                for i in range(len(self.best_peaks)):
                    if not (self.best_peaks[i]>self.z[ind_sig]): # z[ind_sig] will be the minimum height of a peak according to the chosen significance level 
                        break # since best_peaks and best_per are already sorted, if we break at this point i will be the index of the first peak which is too low, and the rest is then too low also
                self.number_of_significant_peaks=i 
            else:
                self.number_of_significant_peaks=0        
        except (RuntimeError, ValueError, TypeError):
            gls_warnings.update_warning_list('Failed to conduct gls, assuming no peaks')
            self.number_of_significant_peaks=0
            self.P_G=[]
            self.z=[]
            self.best_per=[]
            self.best_peaks=[]
       ############### Sorting the best peaks ##########################
 
           
        gls_warnings.print_warning_log()
        return

 
        
# classes for processing RV files

class rvfile(object):
    
    def __init__(self,name,path):
        # let's make sure all the files are sorted, it won't take long and it will be needed for mean_value and first_observation functions
        #command = 'echo "$(sort %s)" > %s'%(path,path) # command to sort the file in place
        #text,flag=run_command_with_timeout(command, 3) # run the command  
        # now save the properties
        #path =  copy_file_to_datafiles(path)

        self.name=name
        self.path=path
        self.reading_in_progress=open(str(self.path),'r') # to be able to read file line by line using read_one_line function
        self.reading_in_progress.close() # for now let's keep it closed, though  

    # function to read just one line
    def read_one_line(self,offset=0.0): 
        if (self.reading_in_progress.closed): # if the file is currently closed we treat it the same way as end of file
            isline=False
            comment=False
            x=0
            y=0
            y_error=0        
        else:
            line=self.reading_in_progress.readline()
            if not line: # end of file, only isline value matters in this case but we need to declare all for output
                isline=False
                comment=False
                x=0
                y=0
                y_error=0
            else:
                isline=True
                if line[0] == '#': #skip comments, in this case x, y, y_error values are not important, just needed for function output
                    comment=True
                    x=0
                    y=0
                    y_error=0
                else:
                    line = line.strip()
                    line = line.split()
                    if is_float(line[1]): #omit lines with invalid RV data, for example 'NaN', '-Inf' etc.
                        comment=False
                        x=float(line[0])
                        y=float(line[1])-offset # subtracting the offset at this point already
                        y_error=float(line[2])    
                    else: # lines with invalid RV data will also be treated as comments
                        comment=True
                        x=0
                        y=0
                        y_error=0
        return [x,y,y_error,isline,comment]

    # numerical integration of RV data to find the mean value as a first proxy for offset
    def mean_value(self):
        self.reading_in_progress=open(self.path,'r') # to be able to read file line by line using read_one_line function
        comment = True
        # In case there is only one data point we want to return it's abscissa, so we will add a tiny rectangle around this point to our integral, but we want it to be tiny, so it normally won't affect the results 
        dt=0.00001
        while(comment):
            x0,y0,y_error,isline,comment = self.read_one_line()
            if not (isline):
                raise InputError('RV file %s contains no data!'%self.name)
            else: # S will be our integral, t0 is going to be the beginning of the interval over which we calculate our integral
                S=y0*dt
                t0=x0-dt              
        while 1: # now we have the starting point, let's go through the rest of the file
            x,y,y_error,isline,comment = self.read_one_line()
            if not (isline):
                break # the file is over, we finish
            elif (comment):
                continue
            else:
                S=S+(x-x0)*(y+y0)/2 # add a trapezoid to the integral 
                x0=x
                y0=y              
        self.reading_in_progress.close()
        return S/(x0-t0) # divide integral by interval length

    # time of first observation in a given dataset
    def first_datapoint(self):
        self.reading_in_progress=open(self.path,'r') # to be able to read file line by line using read_one_line function
        comment = True
        while(comment):
            x0,y0,y_error,isline,comment = self.read_one_line()
            if not (isline):
                raise InputError('RV file %s contains no data!'%self.name)
        self.reading_in_progress.close()
        return x0            

class rvfile_list(object): # this will store a list of rvfile objects
    def __init__(self,ndset,names,paths): # in initialization we use default offsets and jitters

        if (ndset>NDSETMAX):
            raise InputError('Too many data sets! Maximum number is %d'%NDSETMAX)
        self.ndset=ndset
        self.files=[]
        self.time=[]
        self.rvs=[]
        self.rv_err=[]
        self.idset=[]
        for i in range(self.ndset):
            self.files.append(rvfile(names[i],paths[i])) # creating a rvfile object for each data set
        
         
    def add_datafile(self,name,path): # add another data file to the list
        #path =  copy_file_to_datafiles(path)

        if (os.path.isfile(path)):
            self.files.append(rvfile(name,path))
            self.ndset=self.ndset+1    
            flag=1     
        else:
            warnings=Warning_log(['Path %s does not correspond to a valid file! Request ignored'%path], 'Adding a dataset %s'%name)
            warnings.print_warning_log()
            flag=0
        return flag
         
    def remove_datafile(self,number): # add another data file to the list
        if not (number<self.ndset):
            warnings=Warning_log(['Dataset index outside of range'], 'Removing a dataset %d'%number)
            warnings.print_warning_log()
            flag=0
        else:
            self.ndset-=1
            self.files.pop(number)
            self.idset  = self.idset[np.where(self.idset != number)]
            self.time   = self.time[np.where(self.idset != number)]
            self.rvs    = self.rvs[np.where(self.idset != number)]
            self.rv_err = self.rv_err[np.where(self.idset != number)] 
            self.idset[self.idset > number] -= 1  
            flag=1     
 
        return flag     
         
    def first_observation (self): # search for the earliest observation
        if (len(self.files)==0): # shouldn't happen, but just in case
            raise Exception ('No RV files provided!')
        else: 
            # each file is by now sorted, so the first line of the file is it's earliest Julian Day
            firstone = self.files[0].first_datapoint()
            for i in range(1,len(self.files)): # search all datafiles for earliest JD of all
                anotherone = self.files[i].first_datapoint() # to compare with current minimum
                if (anotherone<firstone): 
                    firstone=anotherone
        return firstone
        
   # def read_rvfiles(self,offsets,justthenewone=False):
        
        
    #    x =
   #     y = 
    #    y_error = 
    #    data_set =       
   


    def read_rvfiles(self,offsets,justthenewone=False):
        i = 0
        x = []
        y = []
        y_error = []
        data_set = []
        if (justthenewone): # this means we had already read rvfiles, but we're calling a function to read a new freshly added one
            i=0#self.ndset-1 #not working!!!!
        else:
            i=0
        while i<self.ndset:
            self.files[i].reading_in_progress=open(self.files[i].path,'r')
            while 1:
                xx, yy, yy_error, isline, comment = self.files[i].read_one_line(offsets[i])
                if not (isline):
                    break
                elif (comment):
                    continue
                else:
                    x.append(xx)
                    y.append(yy)
                    y_error.append(yy_error)
                    data_set.append(i)       
            self.files[i].reading_in_progress.close()   
                   
            i += 1
        #casting all lists on arrays which work faster 
        x = np.array(x)
        y = np.array(y)
        y_error = np.array(y_error)
        data_set = np.array(data_set)
       

        # saving sorted arrays to appropriate object attributes        
        self.time = x #np.concatenate((self.time,x))
        self.rvs = y#np.concatenate((self.rvs,y))
        self.rv_err = y_error#np.concatenate((self.rv_err,y_error))
        self.idset = data_set
        
        dataset_ind = np.concatenate((self.idset,data_set)) #there is a problem with np.concatenate() it always returns floats not integers
        self.idset = dataset_ind.astype(int) # this is a quick fix
        #print(self.idset)

        #sorting data by time
        ind_sort = np.argsort(self.time)
        self.time = self.time[ind_sort]
        self.rvs = self.rvs[ind_sort]
        self.rv_err = self.rv_err[ind_sort]
        self.idset = self.idset[ind_sort]
        return 
         
                    
class signal_fit(object):
 
    def __init__(self, inputfile='init.init', name='', readinputfile=False): 
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
        self.GP_params=GP_parameters(4,[1,10,15,1],kernel_id=0) # we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
       
        ########## new stuff ##########
        self.init_bounds()
        
        self.fit_performed = False
        self.fitting_method = 'None'
        self.model_saved=False
        self.stat_saved=False      
        self.masses=[0.0]*10
        self.semimajor=[0.0]*10
        self.f_for_mcmc=[]
        self.par_for_mcmc=[]
        self.e_for_mcmc=[]  
        self.b_for_mcmc=[] 
        
        self.gps=[]
        if (readinputfile):
            self.read_input_file() # after running this a file should be processed and initial values saved in the object properties, unless there was some error
            self.correct_elements() # to make sure angles are between 0 and 360 etc. 
            self.prepare_for_mcmc() # this will initiate default parameter bounds, in case the user wants to verify if values make sense he will be able to use the verify_params_with_bounds function anytime 
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
        
        self.init_transit_params()
        #self.tr_params = batman.TransitParams() 
        self.colors = ['#0066ff',  '#ff0000','#66ff66','#00ffff','#cc33ff','#ff9900','#cccc00','#3399ff','#990033','#339933','#666699']
       
        self.mcmc_sample_file = 'mcmc_samples'
        self.corner_plot_file = 'cornerplot.pdf'
      
        
        #self.print_info=()        
        ##### this is how I wanted the kernel parameters to be 
        ##### initially defined: as a python dictionary and only to filled in by functions! 
        ##### It will be a lot of work to change those above now, but it has to be done as some point.
        
        self.evol_T = {k: [] for k in range(9)}
        self.evol_a = {k: [] for k in range(9)}
        self.evol_e = {k: [] for k in range(9)}
        self.evol_p = {k: [] for k in range(9)}
        self.evol_M = {k: [] for k in range(9)}
        
        
        self.act_data_sets = {k: [] for k in range(10)}
        self.tra_data_sets = {k: [] for k in range(10)}
        self.rad_data_sets = {k: [] for k in range(10)}
        


    def init_bounds(self): 

        self.K_bound   = {k: np.array([0,10000]) for k in range(9)}
        self.P_bound    = {k: np.array([0,100000]) for k in range(9)}
        self.e_bound    = {k: np.array([-0.999,0.999]) for k in range(9)}
        self.w_bound    = {k: np.array([-2.0*360.0, 2.0*360.0]) for k in range(9)}
        self.M0_bound   = {k: np.array([-2.0*360.0, 2.0*360.0]) for k in range(9)}
        self.i_bound    = {k: np.array([-180.0, 180.0]) for k in range(9)}
        self.Node_bound = {k: np.array([-2.0*360.0, 2.0*360.0]) for k in range(9)}
        
        self.rvoff_bounds = {k: np.array([-1000000.0,1000000.0]) for k in range(10)} 
        self.jitt_bounds  = {k: np.array([0.0,10000.0] )for k in range(10)} 
        
        self.lintr_bounds  = {k: np.array([-1.0,1.0]) for k in range(1)} 
        self.st_mass_bounds  = {k: np.array([0.01,100]) for k in range(1)} 


        self.GP_bounds  = {k: np.array([0.0,100000.0]) for k in range(4)} 
 

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
    

        self.tr_params.limb_dark = "quadratic"        #limb darkening model
        self.tr_params.u = [0.1, 0.3 ]           
      
        self.tr_params_use = [True, True,False,False,True,False,True]    
        #self.tr_params_use = [False, False,False,False,False,False,False]    
       
        
        
        self.tr_par = [self.tr_params.t0,  
        self.tr_params.per, 
        self.tr_params.ecc, 
        self.tr_params.w,             
        self.tr_params.rp,
        self.tr_params.inc, 
        self.tr_params.a,        
        ]

        self.tr_el_str  = [r't0', r'P', r'e',r'omega [deg]',r'rp[Rsol]', r'i [deg]' ,r'a [Rsol]']     
 
    
        self.tr_bounds = [[-2460000.0, 2460000.0],[0.5, 4.0],[0.0, 0.999],[0.0, 359.9],[0, 0.20],[84.0, 96.0001],[1, 100]] # planet 1
      
 
        
         
    def update_epoch(self,epoch):
        self.epoch=epoch
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
############################ activity datasets END ##########################################      

############################ transit datasets ##########################################      
    def add_transit_dataset(self, name, path, tra_idset = 0):
 
        tra_JD       = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [0])
        tra_data     = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [1])
        tra_data_sig = np.genfromtxt("%s"%(path),skip_header=0, unpack=True,skip_footer=0, usecols = [2])
	
        tra_data_set = np.array([tra_JD,tra_data,tra_data_sig]) 
 
        self.tra_data_sets[tra_idset] = tra_data_set
 
        return   


    def remove_transit_dataset(self, tra_idset):
 
        self.tra_data_sets[tra_idset] = []
 
        return   
############################ activity datasets END ##########################################      





    def add_planet(self,K=50,P=100,e=0,w=0,M0=180,i=90,cap=0,useK=True,useP=True,usee=False,usew=False,useM0=True,usei=True, usecap=True):        
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
        filenames=np.char.array(alltheinput[2])  # the third line contains a list of RV file names
        self.filelist=rvfile_list(int(alltheinput[1][0]),filenames,path+filenames) # creating a rvfile_list object
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
        if(self.inputfile_read):
            message_str = message_str +"""This is the information for signal fit %s, based on input file %s. 
            
"""%(self.name, self.inputfile) 
        else:
            message_str = message_str +"""This is the information for signal fit %s.
            
"""%self.name          
        # checking if we are dealing with fitted parameters or original data
        if (self.fit_performed and not self.never_saved):
            message_str = message_str +"""Presenting optimal parameters fitted using %s fitting method.
            
"""%self.fitting_method
            message_str = message_str +"""Fit properties: \n chi^2: %f \n reduced chi^2: %f \n rms: %f \n loglik: %f
            
"""%(self.fit_results.chi2,self.fit_results.reduced_chi2,self.fit_results.rms,self.fit_results.loglik)
        else:
            message_str = message_str +"""No fit has yet been conducted (or parameters have never been saved), presenting parameters from user\'s original input

"""

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
            else:
                message_str = message_str +"""\nThe mass of the host star is not known. Default value of %f is assumed.
"""%DEFAULT_STELLAR_MASS
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
 
        return message_str  
        
  
    def fortran_input(self, program='chi2_kep', fileinput=False, filename='Kep_input', amoeba_starts=1, outputfiles=[1,1,1],eps='1.0E-8',dt=864000, when_to_kill=300, npoints=50, model_max = 100, model_min =0): # generate input string for the fortran code, optionally as a file


        ### ppp will be the input string. Depending on fileinput parameter we either save it in a file or save it directly 
     
        if not (fileinput): # if we want to save input in a file we don't want this line in the input string    
            ppp = './%s << EOF\n'%program
        else:
            ppp = '' 
        ppp+= '%s %f %d %d %d %d %d\n%f %d %d %d \n%d\n'%(eps,dt,amoeba_starts,when_to_kill,npoints, model_max, model_min, self.params.stellar_mass,outputfiles[0], outputfiles[1],outputfiles[2],self.filelist.ndset) # first three lines of fortran input: precision and timestep for integration, stellar mass and number of datasets
        for i in range(self.filelist.ndset): 
            # path for each dataset      
            ppp+='%s\n'%(self.filelist.files[i].path)    
            # offset and jitter information for each dataset
            ppp+='%f\n%d\n'%(self.params.offsets[i],int(self.use.use_offsets[i]))
            ppp+='%f\n%d\n'%(self.params.jitters[i],int(self.use.use_jitters[i]))  
        ppp+='%d\n'%self.npl
        for i in range(self.npl): # K,P,e,w,M,i,cap0m for each planet, and information which ones we use
            ppp+='%f %f %f %f %f %f %f\n'%(self.params.planet_params[7*i],self.params.planet_params[7*i+1],self.params.planet_params[7*i+2],self.params.planet_params[7*i+3],self.params.planet_params[7*i+4],self.params.planet_params[7*i+5],self.params.planet_params[7*i+6])
            ppp+='%d %d %d %d %d %d %d\n'%(int(self.use.use_planet_params[7*i]),int(self.use.use_planet_params[7*i+1]),int(self.use.use_planet_params[7*i+2]),int(self.use.use_planet_params[7*i+3]),int(self.use.use_planet_params[7*i+4]),int(self.use.use_planet_params[7*i+5]),int(self.use.use_planet_params[7*i+6]))     
        ppp+='%f\n%d\n'%(self.params.linear_trend,int(self.use.use_linear_trend)) # information about linear trend
        ppp+='%f\n'%self.epoch
        # prepare final version of the ppp command to be returned by this function
        if not (fileinput):
            ppp+='EOF' # end of the command to run in the case of saving input directly
        else: # here's what we do if we want to generate a file as well 
            # first we save the ppp string in a file (by default 'Kep_input')
            file_kep = open(filename, 'w')
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
        return
        
        
        
    def fitting_SciPyOp(self,  Kbounds=[[0.0,100000.0]],Pbounds=[[0.0,100000.0]],ebounds=[[-0.99,0.99]],wbounds=[[-2.0*360.0, 2.0*360.0]],Mbounds=[[-2.0*360.0, 2.0*360.0]],ibounds=[[-2.0*180.0, 2.0*180.0]],capbounds=[[-2.0*360.0, 2.0*360.0]],offbounds=[[-100000.0,100000.0]],jitbounds=[[0.0,10000.0]],lintrbounds=[[-10.0,10.0]], GPbounds=[[0.0,100000.0]], stmassbounds=[[0.01,1000.0]], prior=0, samplesfile='', level=(100.0-68.3)/2.0, threads=1, doGP=False, gp_par=None, kernel_id=-1, use_gp_par=[False,False,False,False], save_means=False, fileoutput=False, save_sampler=False,burning_ph=20, mcmc_ph=20, **kwargs):      

         #def fittingSciPyOp(, doGP=False, gp_par=None, kernel_id=-1, use_gp_par=[False,False,False,False]):    
         ####### find the -LogLik "minimum" using the "truncated Newton" method ######### 
    
        if (doGP):
            self.initiategps(gp_par=gp_par, use_gp_par=use_gp_par, kernel_id=kernel_id)
            nll = lambda *args: -compute_loglik_SciPyOp2(*args)
        else:
            nll = lambda *args: -compute_loglik_SciPyOp(*args)
            GPbounds=[[x-10.0,x+10.0] for x in self.params.GP_params.gp_par] # just to make sure GPbounds don't cause lnprob return -infinity when we don't do GP (now all GP_params will be within bounds for sure)
     
 
        self.prepare_for_mcmc(Kbounds=Kbounds,Pbounds=Pbounds,ebounds=ebounds,wbounds=wbounds,Mbounds=Mbounds,ibounds=ibounds,capbounds=capbounds,offbounds=offbounds,jitbounds=jitbounds,lintrbounds=lintrbounds, GPbounds=GPbounds, stmassbounds=stmassbounds)    
        pp = self.par_for_mcmc
 
    
        # b = np.array(obj.bounds.offset_bounds,obj.bounds.jitter_bounds,obj.bounds.planet_params_bounds,
        #                    obj.bounds.linear_trend_bounds,obj.bounds.GP_params_bounds,obj.bounds.stellar_mass_bounds)
        # b.flatten()
        # print(b)
        # bounds=b,
        minimzers = ['Nelder-Mead','Powell','CG','BFGS','Newton-CG','L-BFGS-B', 'TNC','COBYLA','SLSQP','dogleg','trust-ncg']

        # for k in range(2): # run at least 3 times the minimizer
        #     result = op.minimize(nll, pp, args=(obj), method=minimzers[6], bounds=None, options={'xtol': 1e-6, 'disp': True })
        #     pp = result["x"]
         
        #  print("Best fit par.:", result["x"])
        #----------------- one more time using the Simplex method ---------------------#
        xtol = 1e-3
        for k in range(3): # run at least 3 times the minimizer
            xtol = xtol/10.0 

            result = op.minimize(nll, pp, args=(self), method=minimzers[0], options={'xtol': xtol, 'disp': True, 'maxiter':10000, 'maxfev':10000 })
            pp = result["x"]

        print("Best fit par.:", result["x"])
    
        self.par_for_mcmc = result["x"]        
        
        #self.initiategps = None
        self.gps= None
                
    ### this function is a wrapper calling a fortran program to fit parameters in keplerian mode by minimizing chi^2   WORK IN PROGRESS ON THAT ONE! 
        
    def fitting(self, minimize_fortran=True, minimize_loglik=False, fileinput=False, doGP=False, gp_par=None, kernel_id=-1, use_gp_par=[False,False,False,False], filename='Kep_input', outputfiles=[1,1,1], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 500,model_min =0): # run the fit which will either minimize chi^2 or loglik.
        '''       
         eps, dt - accuracy and step for the integration in dynamical case, in keplerian case ignored, just given for input consistency
         which value to minimize (used for calling an appropriate fortran program)
        '''      
        
        eps = eps * 1e-13 
        dt  = dt  * 86400.0         
        
        if(minimize_loglik):
            minimized_value='loglik'
        else:
            minimized_value='chi2'
        # which mode to use, analogously as above
        if(self.mod_dynamical):
            mod='dyn'
        else:
            mod='kep'
            
        if(minimize_fortran):   
            program='./lib/fr/%s_%s'%(minimized_value,mod) 
            text,flag=run_command_with_timeout(self.fortran_input(program=program, fileinput=fileinput, filename=filename, outputfiles=outputfiles,amoeba_starts=amoeba_starts,eps=eps,dt=dt, when_to_kill=fortran_kill, npoints=npoints, model_max = model_max, model_min =model_min), timeout_sec, output=True,pipe=(not bool(outputfiles[2]))) # running command generated by the fortran_input function 
           
        else:
            self.fitting_SciPyOp(doGP=doGP, gp_par=gp_par, kernel_id=kernel_id, use_gp_par=use_gp_par)         
            program='./lib/fr/%s_%s'%(minimized_value,mod) 
            text,flag=run_command_with_timeout(self.fortran_input(program=program, fileinput=fileinput, filename=filename, outputfiles=outputfiles,amoeba_starts=0,eps=eps,dt=dt, when_to_kill=fortran_kill, npoints=npoints, model_max = model_max, model_min =model_min), timeout_sec, output=True,pipe=(not bool(outputfiles[2]))) # running command generated by the fortran_input function 
           
                  
        if (flag==1):
            fortranoutput=fortran_output(text,self.npl,self.filelist.ndset,self.params.stellar_mass) # create an object for fortran output, which we need to split into part
            self.fit_results=fortranoutput.modfit(print_stat=print_stat)
            self.stat_saved=self.fit_results.stat_array_saved
            if (self.stat_saved):
                self.never_saved=False
            self.model_saved=bool(outputfiles[1])
            self.fit_performed=True
            if(self.fit_results.stat_array_saved): 
                self.fitting_method=program
            self.update_with_fit_results()
            self.correct_elements() #because amoeba might make things wrong here
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
                    

    #### WORK IN PROGRESS HERE!
   # def prepare_for_mcmc(self, customdatasetlabels=[]):
    def prepare_for_mcmc(self, Kbounds=[[0.0,100000.0]],Pbounds=[[0.0,100000.0]],ebounds=[[-0.99,0.99]],wbounds=[[-2.0*360.0, 2.0*360.0]],Mbounds=[[-2.0*360.0, 2.0*360.0]],ibounds=[[-2.0*180.0, 2.0*180.0]],capbounds=[[-2.0*360.0, 2.0*360.0]],offbounds=[[-100000.0,100000.0]],jitbounds=[[0.0,10000.0]],lintrbounds=[[-100000.0,100000.0]], GPbounds=[[0.0,100000.0]], stmassbounds=[[0.01,1000.0]], customdatasetlabels=[]):
        '''Prepare bounds and par array needed for the MCMC''' 
 
        preparingwarnings=Warning_log([],'Preparing for MCMC')  
        # put together bounds for K,P,e,w,M0, so they are in an order we are used to
        
        parambounds=[[0.0,100000.0]]*(7*self.npl) # set some initial values just in case
        offbounds  =[[0.0,100000.0]]*(self.filelist.ndset) # set some initial values just in case
        jitbounds  =[[0.0,100000.0]]*(self.filelist.ndset) # set some initial values just in case
        GPbounds  =[[0.0,100000.0]]*(4) # set some initial values just in case
        ltrbound  = [[-1,1.0]] 
        stmass_bound = [[0.01,100.0]]     
 
 
        for i in range(self.npl):
            parambounds[7*i]=self.K_bound[i]        
            parambounds[7*i+1]=self.P_bound[i] 
            parambounds[7*i+2]=self.e_bound[i]  
            parambounds[7*i+3]=self.w_bound[i]   
            parambounds[7*i+4]=self.M0_bound[i] 
            parambounds[7*i+5]=self.i_bound[i]
            parambounds[7*i+6]=self.Node_bound[i]
 
 
        for i in range(self.filelist.ndset):           
            offbounds[i] = self.rvoff_bounds[i]
            jitbounds[i] = self.jitt_bounds[i]
        # Concatenate bounds into one array
        
        for i in range(4):           
            GPbounds[i] = self.GP_bounds[i]
        
        for i in range(1):         
            ltrbound[i] = self.lintr_bounds[i]
            stmass_bound[i] = self.st_mass_bounds[i]
        
        self.bounds = parameter_bounds(offbounds,jitbounds,parambounds,self.lintr_bounds,self.GP_bounds,self.st_mass_bounds)     
        
        # Now prepare parameters, only those which are used         
        par = np.concatenate((self.params.offsets[:self.filelist.ndset],self.params.jitters[:self.filelist.ndset],self.params.planet_params[:7*self.npl],np.atleast_1d(self.params.linear_trend),np.atleast_1d(self.params.GP_params.gp_par),np.atleast_1d(self.params.stellar_mass)))         
        flag = np.concatenate((self.use.use_offsets[:self.filelist.ndset],self.use.use_jitters[:self.filelist.ndset],self.use.use_planet_params[:7*self.npl],np.atleast_1d(self.use.use_linear_trend),np.atleast_1d(self.use.use_GP_params),np.atleast_1d(self.use.use_stellar_mass)))
      
        bounds = np.concatenate((offbounds[:self.filelist.ndset],jitbounds[:self.filelist.ndset],parambounds[:7*self.npl], GPbounds[:4],ltrbound, stmass_bound))   
         
        
       # print(par,flag)
        if not (self.mod_dynamical): # we need to make sure we don't try to fit for inclination in keplerian case
            for i in range(self.npl):
                flag[2*self.filelist.ndset+7*i+5]=False
                flag[2*self.filelist.ndset+7*i+6]=False


        # prepare element names for corner plot labels
        if (len(customdatasetlabels)<self.filelist.ndset):
            if not (customdatasetlabels==[]): # this means the user wanted to provide custom dataset labels, but didn't provide it for all datasets, we need to give a warning
                preparingwarnings.update_warning_list('Too few customdatasetlabels! Will use default labels for remaining.')
            customdatasetlabels=np.concatenate((customdatasetlabels,np.array(list(map(str,np.arange(len(customdatasetlabels),self.filelist.ndset))))))
        if (len(customdatasetlabels)>self.filelist.ndset):
            customdatasetlabels=customdatasetlabels[:self.filelist.ndset]
            preparingwarnings.update_warning_list('Too many customdatasetlabels! Additional will be ignored.')
        el_str=np.concatenate(([r'$\gamma_%s$ [m/s]'%customdatasetlabels[i] for i in range(self.filelist.ndset)],
                               [r'jitt$_%s$ [m/s]'%customdatasetlabels[i] for i in range(self.filelist.ndset)],
                               np.concatenate([[r'K$_%s$ [m/s]'%chr(98+i),
                                                r'P$_%s$ [day]'%chr(98+i),
                                                r'e$_%s$'%chr(98+i),
                                                r'$\omega_%s$ [deg]'%chr(98+i),
                                                r'M$_%s$ [deg]'%chr(98+i), 
                                                r'i$_%s$ [deg]'%chr(98+i), 
                                                r'$\Omega_%s$ [deg]'%chr(98+i)] for i in range(self.npl)]),np.atleast_1d(r'lin.trend [m/s/day]'),[r'Amp', r't', r'per', r'fact'],np.atleast_1d(r'st. mass [$M_\odot$]')))
   
        self.f_for_mcmc = [idx for idx in range(len(flag)) if flag[idx] ==1 ] # indices for fitted parameters

        self.par_for_mcmc = []  # self par_for_mcmc are the fitted parameters   
        self.e_for_mcmc = [] # labels for fitted parameters only
        self.b_for_mcmc = [] # labels for fitted parameters only
        

        for j in range(len(par)):
            if flag[j] > 0:
                self.par_for_mcmc=np.concatenate((self.par_for_mcmc,np.atleast_1d(par[j])))
                self.e_for_mcmc=np.concatenate((self.e_for_mcmc,np.atleast_1d(el_str[j])))
                #self.b_for_mcmc=np.concatenate((self.b_for_mcmc,np.atleast_1d(bounds[j])))
                self.b_for_mcmc.append(bounds[j])

            
 
        preparingwarnings.print_warning_log()
        return                  
        
    #### WORK IN PROGRESS HERE!
    def prepare_for_mcmcold(self, Kbounds=[[0.0,100000.0]],Pbounds=[[0.0,100000.0]],ebounds=[[-0.99,0.99]],wbounds=[[-2.0*360.0, 2.0*360.0]],Mbounds=[[-2.0*360.0, 2.0*360.0]],ibounds=[[-2.0*180.0, 2.0*180.0]],capbounds=[[-2.0*360.0, 2.0*360.0]],offbounds=[[-100000.0,100000.0]],jitbounds=[[0.0,10000.0]],lintrbounds=[[-100000.0,100000.0]], GPbounds=[[0.0,100000.0]], stmassbounds=[[0.01,1000.0]], customdatasetlabels=[]):
 
        '''Prepare bounds and par array needed for the MCMC''' 
 
        preparingwarnings=Warning_log([],'Preparing for MCMC') 
 
        # set default bounds for these parameters which were not provided. Make warnings if there's too many bounds
        if(len(Kbounds)<self.npl):
            if not (Kbounds==[[0.0,100000.0]]): # This means user intentionally provided some bounds, but fewer than there should be
                preparingwarnings.update_warning_list('Too few Kbounds! Assuming default for remaining planets.') 
            Kbounds=np.concatenate((Kbounds,np.repeat([[0.0,100000.0]],self.npl-len(Kbounds),axis=0)))
        elif(len(Kbounds)>self.npl):
            Kbounds=Kbounds[:self.npl]
            preparingwarnings.update_warning_list('Too many Kbounds! Additional will be ignored.')
        if(len(Pbounds)<self.npl):
            if not (Pbounds==[[0.0,100000.0]]):
                preparingwarnings.update_warning_list('Too few Pbounds! Assuming default for remaining planets.')         
            Pbounds=np.concatenate((Pbounds,np.repeat([[0.0,100000.0]],self.npl-len(Pbounds),axis=0)))
        elif (len(Pbounds)>self.npl):
            Pbounds=Pbounds[:self.npl]
            preparingwarnings.update_warning_list('Too many Pbounds! Additional will be ignored.')   
        if(len(ebounds)<self.npl):
            if not (ebounds==[[-0.99,0.99]]):
                preparingwarnings.update_warning_list('Too few ebounds! Assuming default for remaining planets.')         
            ebounds=np.concatenate((ebounds,np.repeat([[-0.99,0.99]],self.npl-len(ebounds),axis=0)))
        elif (len(ebounds)>self.npl):
            ebounds=ebounds[:self.npl]
            preparingwarnings.update_warning_list('Too many ebounds! Additional will be ignored.')  
        if(len(wbounds)<self.npl):
            if not (wbounds==[[-2.0*360.0, 2.0*360.0]]):
                preparingwarnings.update_warning_list('Too few wbounds! Assuming default [-360.0*2.0, 360.0*2.0] for remaining planets.')         
            wbounds=np.concatenate((wbounds,np.repeat([[-2.0*360.0, 2.0*360.0]],self.npl-len(wbounds),axis=0)))
        elif (len(wbounds)>self.npl):
            wbounds=wbounds[:self.npl]
            preparingwarnings.update_warning_list('Too many wbounds! Additional will be ignored.')       
        if(len(Mbounds)<self.npl):
            if not (Mbounds==[[-2.0*360.0, 2.0*360.0]]):
                preparingwarnings.update_warning_list('Too few Mbounds! Assuming default [-360.0*2.0, 360.0*2.0] for remaining planets.')         
            Mbounds=np.concatenate((Mbounds,np.repeat([[-2.0*360.0, 2.0*360.0]],self.npl-len(Mbounds),axis=0)))           
        elif (len(Mbounds)>self.npl):
            Mbounds=Mbounds[:self.npl]
            preparingwarnings.update_warning_list('Too many Mbounds! Additional will be ignored.') 
        if(len(ibounds)<self.npl):
            if not (ibounds==[[-2.0*180.0, 2.0*180.0]]):
                preparingwarnings.update_warning_list('Too few ibounds! Assuming default [-180.0*2.0, 180.0*2.0] for remaining planets.')         
            ibounds=np.concatenate((ibounds,np.repeat([[-2.0*180.0, 2.0*180.0]],self.npl-len(ibounds),axis=0)))
        elif (len(ibounds)>self.npl):
            ibounds=ibounds[:self.npl]
            preparingwarnings.update_warning_list('Too many ibounds! Additional will be ignored.')  
        if(len(capbounds)<self.npl):
            if not (capbounds==[[-2.0*360.0, 2.0*360.0]]):
                preparingwarnings.update_warning_list('Too few capbounds! Assuming default [-360.0*2.0, 360.0*2.0] for remaining planets.')         
            capbounds=np.concatenate((capbounds,np.repeat([[-2.0*360.0, 2.0*360.0]],self.npl-len(capbounds),axis=0))) 
        elif (len(capbounds)>self.npl):
            capbounds=capbounds[:self.npl]
            preparingwarnings.update_warning_list('Too many capbounds! Additional will be ignored.')   
        if(len(offbounds)<self.filelist.ndset):
            if not (offbounds==[[-100000.0,100000.0]]):
                preparingwarnings.update_warning_list('Too few offbounds! Assuming default [-100000.0, 100000.0] for remaining datasets.')           
            offbounds=np.concatenate((offbounds,np.repeat([[-100000.0,100000.0]],self.filelist.ndset-len(offbounds),axis=0)))
        elif (len(offbounds)>self.filelist.ndset):
            offbounds=offbounds[:self.filelist.ndset]
            preparingwarnings.update_warning_list('Too many offbounds! Additional will be ignored.')  
        if(len(jitbounds)<self.filelist.ndset):
            if not (jitbounds==[[0.0,10000.0]]):
                preparingwarnings.update_warning_list('Too few jitbounds! Assuming default [0.0, 10000.0] for remaining datasets.')           
            jitbounds=np.concatenate((jitbounds,np.repeat([[0.0,10000.0]],self.filelist.ndset-len(jitbounds),axis=0)))           
        elif (len(jitbounds)>self.filelist.ndset):
            jitbounds=jitbounds[:self.filelist.ndset]
            preparingwarnings.update_warning_list('Too many jitbounds! Additional will be ignored.')
        if(len(lintrbounds)<1):
            lintrbounds=[[-100000.0,100000.0]] 
        elif (len(lintrbounds)>1):
            lintrbounds=lintrbounds[:1]
            preparingwarnings.update_warning_list('Too many lintrbounds! Additional will be ignored.') 
        if(len(stmassbounds)<1):
            stmassbounds=[[0.01,1000.0]] 
        elif (len(stmassbounds)>1):
            stmassbounds=stmassbounds[:1]
            preparingwarnings.update_warning_list('Too many stmassbounds! Additional will be ignored.') 
        if(len(GPbounds)<self.params.GP_params.npar):
            if not (GPbounds==[[0.0,100000.0]]):
                preparingwarnings.update_warning_list('Too few GPbounds! Assuming default [0.0, 10000.0] for remaining parameters.')           
            GPbounds=np.concatenate((GPbounds,np.repeat([[0.0,10000.0]],self.params.GP_params.npar-len(GPbounds),axis=0)))           
        elif (len(GPbounds)>self.params.GP_params.npar):
            GPbounds=GPbounds[:self.params.GP_params.npar]
            preparingwarnings.update_warning_list('Too many GPbounds! Additional will be ignored.')
        
        # put together bounds for K,P,e,w,M0, so they are in an order we are used to
        
        parambounds=[[0.0,100000.0]]*(7*self.npl) # set some initial values just in case
        for i in range(self.npl):
            parambounds[7*i]=Kbounds[i]        
            parambounds[7*i+1]=Pbounds[i] 
            parambounds[7*i+2]=ebounds[i]  
            parambounds[7*i+3]=wbounds[i]   
            parambounds[7*i+4]=Mbounds[i] 
            parambounds[7*i+5]=ibounds[i]
            parambounds[7*i+6]=capbounds[i]
                              
        # Concatenate bounds into one array
        
        self.bounds = parameter_bounds(offbounds,jitbounds,parambounds,lintrbounds,GPbounds,stmassbounds)     
        
        # Now prepare parameters, only those which are used         

        par = np.concatenate((self.params.offsets[:self.filelist.ndset],self.params.jitters[:self.filelist.ndset],self.params.planet_params[:7*self.npl],np.atleast_1d(self.params.linear_trend),np.atleast_1d(self.params.GP_params.gp_par),np.atleast_1d(self.params.stellar_mass)))
           
        flag = np.concatenate((self.use.use_offsets[:self.filelist.ndset],self.use.use_jitters[:self.filelist.ndset],self.use.use_planet_params[:7*self.npl],np.atleast_1d(self.use.use_linear_trend),np.atleast_1d(self.use.use_GP_params),np.atleast_1d(self.use.use_stellar_mass)))
        
       # print(par,flag)


        if not (self.mod_dynamical): # we need to make sure we don't try to fit for inclination in keplerian case
            for i in range(self.npl):
                flag[2*self.filelist.ndset+7*i+5]=False
                flag[2*self.filelist.ndset+7*i+6]=False


        # prepare element names for corner plot labels
        if (len(customdatasetlabels)<self.filelist.ndset):
            if not (customdatasetlabels==[]): # this means the user wanted to provide custom dataset labels, but didn't provide it for all datasets, we need to give a warning
                preparingwarnings.update_warning_list('Too few customdatasetlabels! Will use default labels for remaining.')
            customdatasetlabels=np.concatenate((customdatasetlabels,np.array(list(map(str,np.arange(len(customdatasetlabels),self.filelist.ndset))))))
        if (len(customdatasetlabels)>self.filelist.ndset):
            customdatasetlabels=customdatasetlabels[:self.filelist.ndset]
            preparingwarnings.update_warning_list('Too many customdatasetlabels! Additional will be ignored.')
        el_str=np.concatenate(([r'$\gamma_%s$ [m/s]'%customdatasetlabels[i] for i in range(self.filelist.ndset)],
                               [r'jitt$_%s$ [m/s]'%customdatasetlabels[i] for i in range(self.filelist.ndset)],
                               np.concatenate([[r'K$_%s$ [m/s]'%chr(98+i),
                                                r'P$_%s$ [day]'%chr(98+i),
                                                r'e$_%s$'%chr(98+i),
                                                r'$\omega_%s$ [deg]'%chr(98+i),
                                                r'M$_%s$ [deg]'%chr(98+i), 
                                                r'i$_%s$ [deg]'%chr(98+i), 
                                                r'$\Omega_%s$ [deg]'%chr(98+i)] for i in range(self.npl)]),np.atleast_1d(r'lin.trend [m/s/day]'),[r'Amp', r't', r'per', r'fact'],np.atleast_1d(r'st. mass [$M_\odot$]')))
   
        self.f_for_mcmc = [idx for idx in range(len(flag)) if flag[idx] ==1 ] # indices for fitted parameters

        self.par_for_mcmc = []  # self par_for_mcmc are the fitted parameters   
        self.e_for_mcmc = [] # labels for fitted parameters only
        

        for j in range(len(par)):
            if flag[j] > 0:
                self.par_for_mcmc=np.concatenate((self.par_for_mcmc,np.atleast_1d(par[j])))
                self.e_for_mcmc=np.concatenate((self.e_for_mcmc,np.atleast_1d(el_str[j])))
 
            
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
                elif (np.mod(nr,7)==6):
                    newparams.update_lineofnodes(x,p[i])
            elif (idx<2*self.filelist.ndset+7*self.npl+1):
                newparams.update_linear_trend(p[i]) 
                i=i+1       
            elif (idx<2*self.filelist.ndset+7*self.npl+1+self.params.GP_params.npar):
                newparams.update_GP_param_value(idx-2*self.filelist.ndset-7*self.npl-1,p[i])
                i=i+1
            else:         
                newparams.update_stellar_mass(p[i])      
        return newparams                
                
    def update_with_mcmc_errors(self,p):
 
        '''Substitute normal errors with mcmc errors, where + and - errors can be different''' 
 
        if (not (self.fit_performed) or self.never_saved): # just in case the user calls this function in a wrong moment
            return
        else:
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
                elif (idx<2*self.filelist.ndset+7*self.npl+1+self.params.GP_params.npar):
                    self.param_errors.update_GP_param_errors(idx-2*self.filelist.ndset-7*self.npl-1,p[i])
                    i=i+1
                else:    
                    self.param_errors.update_stellar_mass_error(p[i])            
        return                           

    def initiategps(self, gp_par=None, use_gp_par=[], kernel_id=-1): 
        
        # Prepare objects for Gaussian Processes        
        
        # Redefine GP parameters if new ones were provided
        
        #if not (gp_par==None):
        if len(gp_par) != 0:
            self.params.update_GP_params(gp_par,kernel_id=kernel_id)
            self.use.update_use_GP_params(use_gp_par)
            self.verify_gp_parameters_number()
        
#        kernel_jitters=[]
        kernels=[]
        gps=[]
        #print(gp_par[0],gp_par[1],gp_par[2],gp_par[3] )
       #copied_obj.gps[i].set_parameter_vector(np.array(list(map(np.log,np.concatenate((copied_obj.params.GP_params.gp_par,np.atleast_1d(copied_obj.params.jitters[i]))))))) 

        for i in range(self.filelist.ndset):
            kernels.append(self.params.GP_params.rot_kernel +terms.JitterTerm(np.log(self.params.jitters[i])))
            gps.append(celerite.GP(kernels[i], mean=0.0))
            gps[i].compute(self.filelist.time[self.filelist.idset==i],self.filelist.rv_err[self.filelist.idset==i])
            #gps[i].compute(self.filelist.time[self.filelist.idset==i])
        self.gps=gps
        
       # print(self.params.GP_params.gp_par)    
       # print(self.use.use_GP_params)    
        
 
                   
        return

    
    def mcmc(self,  Kbounds=[[0.0,100000.0]],Pbounds=[[0.0,100000.0]],ebounds=[[-0.99,0.99]],wbounds=[[-2.0*360.0, 2.0*360.0]],Mbounds=[[-2.0*360.0, 2.0*360.0]],ibounds=[[-2.0*180.0, 2.0*180.0]],capbounds=[[-2.0*360.0, 2.0*360.0]],offbounds=[[-100000.0,100000.0]],jitbounds=[[0.0,10000.0]],lintrbounds=[[-10.0,10.0]], GPbounds=[[0.0,100000.0]], stmassbounds=[[0.01,1000.0]], prior=0, samplesfile='', level=(100.0-68.3)/2.0, threads=1, doGP=False, gp_par=None, gp_kernel_id=-1, use_gp_par=[True,True,True,True], save_means=False, fileoutput=False, save_sampler=False,burning_ph=10, mcmc_ph=10, **kwargs):      

        '''Performs MCMC and saves results'''  
        
        if threads == 'max':
            threads = multiprocessing.cpu_count()
 
        
        # Let's prepare things depending if we want to do GP or not

        if (doGP):
            self.initiategps(gp_par=gp_par, use_gp_par=use_gp_par, kernel_id=gp_kernel_id)
            fun=lnprobGP3                            
        else:
            fun=lnprob                           
            GPbounds=[[x-10.0,x+10.0] for x in self.params.GP_params.gp_par] # just to make sure GPbounds don't cause lnprob return -infinity when we don't do GP (now all GP_params will be within bounds for sure)
 
        # we will need this later
        if (self.mod_dynamical):
            mod='dyn'
        else:
            mod='kep'
    
        self.prepare_for_mcmc(Kbounds=Kbounds,Pbounds=Pbounds,ebounds=ebounds,wbounds=wbounds,Mbounds=Mbounds,ibounds=ibounds,capbounds=capbounds,offbounds=offbounds,jitbounds=jitbounds,lintrbounds=lintrbounds, GPbounds=GPbounds, stmassbounds=stmassbounds)    
 
    
        start_time = time.time()
    
        ndim, nwalkers = len(self.par_for_mcmc), len(self.par_for_mcmc)*4

        #par_for_mcmc_ = self.par_for_mcmc

        pos = [self.par_for_mcmc + 1e-3*np.random.rand(ndim) for i in range(nwalkers)]

        sampler = CustomSampler(nwalkers, ndim, fun, args=[self,prior], threads = threads)
 
        
        # burning phase
        pos, prob, state  = sampler.run_mcmc(pos,burning_ph)

        sampler.reset()
 
        # now perform the MCMC

        pos, prob, state  = sampler.run_mcmc(pos,mcmc_ph)
        
        
                  
        print("--- %s seconds ---" % (time.time() - start_time))  
 
        ln = np.hstack(sampler.lnprobability)
        sampler.save_samples(self.f_for_mcmc,self.filelist.ndset,self.npl)
                
        
        if (fileoutput):
            #if (samplesfile==''): # that means no file name for samples file has been provided, so we generate a default one
           #     samplesfile='samples_%s'%mod
            #self.mcmc_sample_file = 'mcmc_samples'+'_%s'%mod
            #self.corner_plot_file = 'cornerplot.png'
            outfile = open(str(self.mcmc_sample_file), 'w') # file to save samples

            for j in range(len(sampler.samples)):
                outfile.write("%s  " %(ln[j]))        
                for z in range(len(self.par_for_mcmc)):
                    outfile.write("%s  " %(sampler.samples[j,z]))
                outfile.write("\n")
    
            outfile.close()        
                
        # Now we will save new parameters and their errors (different + and - errors in this case). Flag save_means determines if we want to take means as new best fit parameters or stick to old ones and calculate errors with respect to that           
        if (save_means):
            self.par_for_mcmc = sampler.means # we will not need to keep the old parameters in this attribbute, so let's store the means now
            
        new_par_errors = [[float(self.par_for_mcmc[i] - np.percentile(sampler.samples[:,i], [level])),float(np.percentile(sampler.samples[:,i], [100.0-level])-self.par_for_mcmc[i])] for i in range(len(self.par_for_mcmc))] 
        newparams = self.generate_newparams_for_mcmc(self.par_for_mcmc)        
        #print(newparams.GP_params)
        current_GP_params=newparams.GP_params.gp_par # because calling fitting will overwrite them
       # print(current_GP_params)

        self.overwrite_params(newparams)
        self.fitting(minimize_loglik=True, amoeba_starts=0, outputfiles=[1,1,1]) # this will help update some things 
        self.update_with_mcmc_errors(new_par_errors)
        self.params.update_GP_params(current_GP_params)

        #print(current_GP_params)

        if (doGP):
            self.fitting_method='GP_%s'%mod
        else:
            self.fitting_method='mcmc_%s'%mod  

        if (doGP):
            loglik_to_save = lnprobGP(self.par_for_mcmc,self,prior)
            self.loglik=loglik_to_save         
            self.fit_results.loglik=loglik_to_save  
            
            
            ########### Ugly things are happening here! to be fixed! ##########
            
            fitted_GP = sampler.means[-len(self.params.GP_params.gp_par[self.use.use_GP_params==True]):]
 
            z = 0
            for k in range(self.params.GP_params.npar):
                self.param_errors.GP_params_errors[k] = [0,0]
                if not self.use.use_GP_params[k]:
                    continue
                else:
                    self.params.GP_params.gp_par[k] = fitted_GP[z]
                   # self.param_errors.GP_params_errors[k] = [float(sampler.means[-self.params.GP_params.npar+k] - np.percentile(sampler.means[-self.params.GP_params.npar+k], [level])),float(np.percentile(sampler.means[-self.params.GP_params.npar+k], [100.0-level])-sampler.means[-self.params.GP_params.npar+k])]
                    self.param_errors.GP_params_errors[k] = [float(fitted_GP[z] - np.percentile(sampler.samples[:,-len(fitted_GP)+z],[level])), float(np.percentile(sampler.samples[:,-len(fitted_GP)+z], [100.0-level])-fitted_GP[z])]
                    #print(self.param_errors.GP_params_errors[k], fitted_GP[z], np.percentile(fitted_GP[z],level),level)
                z = z+1  

            self.GP_params=GP_parameters(len(self.params.GP_params.gp_par),self.params.GP_params.gp_par,kernel_id=0) # we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
             
            message_str =""" """
            if(self.fitting_method.startswith('GP')):
                message_str = message_str +"""\nStellar activity was modeled using Gaussian Processes. The resulting GP parameters (means) are as follows:
"""
                message_str = message_str +"""\n A = {0:>7.4f} +/- {4:>7.4f}\n t = {1:>7.4f} +/- {5:>7.4f}\n P = {2:>7.4f} +/- {6:>7.4f}\n f = {3:>7.4f} +/- {7:>7.4f}
""".format(self.GP_params.gp_par[0],self.GP_params.gp_par[1],self.GP_params.gp_par[2],self.GP_params.gp_par[3],max(self.param_errors.GP_params_errors[0]),max(self.param_errors.GP_params_errors[1]),max(self.param_errors.GP_params_errors[2]),max(self.param_errors.GP_params_errors[3]))         
                message_str = message_str +"""
(The GP_params are printed here as these are still not part of the "params" structure. TBD!)                
"""
           # print(message_str)
            ###################################################################
            

        ###############  This is not working! you cannot save the sampler as an atribute and call it back later!
        ###############  See https://github.com/dfm/emcee/issues/148
        if(save_sampler):
            self.sampler=sampler             
            self.sampler_saved=True           
            
        sampler.reset()
       
        return
             
    def cornerplot(self, fileinput=False, level=(100.0-68.3)/2.0, **kwargs): 

        '''Generates a corner plot visualizing the mcmc samples. Optionally samples can be read from a file.'''
        #self.mcmc_sample_file = 'mcmc_samples'+'_%s'%mod
        #self.corner_plot_file = 'cornerplot.png'
        if(fileinput):
            samples=read_file_as_array_of_arrays_mcmc(self.mcmc_sample_file)
        elif(self.sampler_saved):
            samples=self.sampler.samples
        else:
            raise Exception ('Please run mcmc and save sampler or provide a valid samples file!')
        
        fig = corner.corner(samples,bins=25, color="k", reverse=True, upper= True, labels=self.e_for_mcmc, quantiles=[level/100.0, 1.0-level/100.0],levels=(0.6827, 0.9545,0.9973), smooth=1.0, smooth1d=1.0, plot_contours= True, show_titles=True, truths=self.par_for_mcmc, dpi = 300, pad=15, labelpad = 50 ,truth_color ='r', title_kwargs={"fontsize": 12}, scale_hist=True,  no_fill_contours=True, plot_datapoints=True, kwargs=kwargs)
        fig.savefig(self.corner_plot_file)  

 
        #fig.clf()
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
 
        param_file.write("""0.0d0 %s %s
%s %s
        
F T T T T F
0.001 50.0 50.0 -1. T
bin.dat
unknown
"""%(max_time, timestep, max_time/1e4, max_time/1e3 ))
 
        param_file.close()
        
        
        getin_file = open('geninit_j.in', 'wb') 
        getin_file.write("""1 
%s
%s
1.d0
pl.in
    """%(str(self.params.stellar_mass),str(self.npl)))
        
        
  
        for j in range(self.npl):
            getin_file.write('%s \n'%str(self.fit_results.mass[j]/1047.70266835)) 
            getin_file.write('%s %s %s %s %s %s \n'%(str(self.fit_results.a[j]),str(self.params.planet_params[7*j + 2]),str(self.params.planet_params[7*j + 5]),
                                                  str(self.params.planet_params[7*j + 3]),str(self.params.planet_params[7*j + 6]),str(self.params.planet_params[7*j + 4])) ) 
              
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
        
              
