#!/usr/bin/python

from __future__ import print_function
__author__ = 'Trifon Trifonov, Jakub Morawski'

import sys #,os
sys.path.insert(0, './addons')
import jac2astrocen
import gls as gls 


import prior_functions as pr
import numpy as np
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
#from pylab import *
import re
#import nmpfit as mp
from matplotlib import gridspec
from subprocess import PIPE, Popen 
import signal, os
import tempfile, shutil

#from emcee.utils import MPIPool
import time
import copy
import multiprocessing
from threading import Thread
#from astroML.time_series import lomb_scargle
#from astroML.time_series import lomb_scargle_bootstrap
from scipy.signal import argrelextrema

import emcee
import corner
import celerite 
from celerite import terms


TAU=6.2831853071
DEFAULT_STELLAR_MASS=1.0
DEFAULT_JITTER=1.0
NPLMAX=20
NDSETMAX=20
DEFAULT_PATH='./datafiles/'



def mut_incl(i1,i2,capOm):
    '''
    Calculates the mutual inclination of two planets   
    
    input parameters:
    
    i1,i2, Delta Omega: inclinations and diffence of the line of nodes in deg.
  
    
    output parameters:    
     
    Delta i: mutual orbital inclination in deg.
    '''
    fb = np.degrees(np.arccos(((np.cos(np.radians(i1))*np.cos(np.radians(i2)))+
    (np.sin(np.radians(i1))*np.sin(np.radians(i2))*np.cos(np.radians(capOm))))))
    return fb


def create_temporary_copy(path): # not really a good idea......
    '''
    creates a temp_ velocity file in the root directory of the GUI.   
    
    input: full path to the file  
    output: temp_name of the file to be loaded
    '''    

    dirname, basename = os.path.split(path)
    temp_dir = './'#tempfile.gettempdir()
    temp_path = os.path.join(temp_dir, '.temp_'+basename)
    shutil.copy2(path, temp_path)
    #temp = tempfile.NamedTemporaryFile(prefix=basename, dir='./',delete=False)
    return temp_path

def copy_file_to_datafiles(path):
    '''
    creates a temp_ velocity file in the root directory of the GUI.   
    
    input: full path to the file  
    output: temp_name of the file to be loaded
    '''    

    dirname, basename = os.path.split(path)
    temp_dir = './datafiles'#tempfile.gettempdir()   
    temp_path = os.path.join(temp_dir, basename)
    
    os.system("cp %s %s"%(path, temp_path))
    #print(temp_path, path)
   # if os.path.exists(temp_path):
    #    try:
   #         shutil.copy2(path, temp_path)
   #     except shutil.Error:
   #         temp_path == path    
 
        #os.system("cp %s %s"%(path, temp_path))
    #temp = tempfile.NamedTemporaryFile(prefix=basename, dir='./',delete=False)
    return temp_path



def run_command_with_timeout(args, secs, output=False, pipe=False): # set output=True if you need to save the output
    '''
    Run a command and kill if it takes too long.
    '''

    if not (pipe):
        text=tempfile.TemporaryFile() # because PIPE usually has too low capacity
        proc = Popen(args, shell=True, preexec_fn=os.setsid, stdout=text, stderr=text)
    else:
        proc = Popen(args, shell=True, preexec_fn=os.setsid, stdout=PIPE, stderr=PIPE)        
    proc_thread = Thread(target=proc.wait)
    proc_thread.start()
    proc_thread.join(secs)
    if proc_thread.is_alive():
        #print (proc.pid)
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        except OSError:
            pass
        print('Process #{} killed after {} seconds'.format(proc.pid, secs))
        flag = -1
        return '',flag
    if not (pipe):
        text.seek(0)
        string_to_output=text.readlines()
    else:
        text=proc.communicate()[0]
        string_to_output=text.splitlines()
    for i in range(len(string_to_output)):
        string_to_output[i]=string_to_output[i].decode('utf-8').split()
    if not (pipe):
        text.close()    
    flag = 1
    if (output):
        return string_to_output,flag # besides the flag which informs about successful termination we also return all the console output in case we want to save it in a variable
    else:
        return '',flag
         
def is_float(n):
    '''
    Given a string n, verify if it expresses a valid float.
    Casting n to string in case an object of type float or similar is given as an argument
    '''
    return re.match(r'^-?\d*(\.\d+)?(E-?\d+)?$', str(n))   

# Given a float or string, verify if it expresses an integer. Possible to introduce upper and lower bounds and if the inequalities on either side should be strong or weak . 
def is_int(s,bounded=[False,False],bounds=[0,0],equal=[False,False]):
    if is_float(s): # if it is an int, it is certainly float as well 
        n=float(s) # we need n as a number, not as a string, for comparisons with bounds later
        is_an_int=float(s).is_integer()
    else:
        is_an_int=False
    # is_an_int now contains an information if s is an int, but without bounds. Let's introduce bounds:
    if(is_an_int): # if it's not an int at all we don't need to check any further
        if(bounded[0]): # if there is a lower bound let's apply it
            if (n<bounds[0] or (not equal[0] and n==bounds[0])):
                is_an_int=False
    if(is_an_int): # if the lower bound returned False we don't need to check any further
        if(bounded[1]): # if there is a lower bound let's apply it
            if (n>bounds[1] or (not equal[1] and n==bounds[1])):
                is_an_int=False
    return is_an_int

# If save_wrong_lines is enabled we will save a string 'wrong_line' instead of this line and save indices at which this occurred, otherwise we will skip this line
def convert_array_to_float(a,save_wrong_lines=False): 
    converting_warnings=Warning_log([],'Converting array to float')
    b=[]
    if (save_wrong_lines):
        wrong_indices=[]
    for i in range(len(a)):
        if not is_float(a[i]):
            if not (save_wrong_lines):
                converting_warnings.update_warning_list('Array passed to convert_array_to_float function should only contain floats! Line %d skipped'%(i+1))
            else:
                b.append('wrong_line')
                wrong_indices=np.concatenate((wrong_indices,np.atleast_1d(i)))
        else:
            b.append(float(a[i]))
            
    converting_warnings.print_warning_log()
    if (save_wrong_lines):
        return np.array(b),wrong_indices  
    else:
        return np.array(b)

def convert_array_to_int(a, save_wrong_lines=False):
    converting_warnings=Warning_log([],'Converting array to int')
    b=[]
    if (save_wrong_lines):
        wrong_indices=[]
    for i in range(len(a)):
        if not is_int(a[i]):
            if not (save_wrong_lines):
                converting_warnings.update_warning_list('Array passed to convert_array_to_int function should only contain ints! Line %d skipped'%(i+1))
            else:
                b.append('wrong_line')
                wrong_indices=np.concatenate((wrong_indices,np.atleast_1d(i)))
        else:
            b.append(int(a[i]))
            
    converting_warnings.print_warning_log()
    if (save_wrong_lines):
        return np.array(b),wrong_indices  
    else:
        return np.array(b)

#for convenient reading of the input file    
def read_file_as_array_of_arrays(inputfile): 
    a=open(inputfile, 'r')
    b=a.readlines() # b as array of strings
    c=[]
    ic=0 # iterator for values in c
    for i in range(len(b)): 
        b[i]=np.atleast_1d(b[i].split()) # turn a row of b into an array of arrays
        c.append([]) # need to make a separate array so every element is of correct type
        # convert each string that represents a float into float
        for j in range(len(b[i])): 
            if (is_float(b[i][j])):
                c[ic].append(float(b[i][j]))
            elif not (b[i][j][-1]==':'): # ignore comments, which can be place by the user as strings which end with a collon, in the comments use underline instead of space or an error will arise
                c[ic].append(b[i][j])
        ic=ic+1
    return c
    
def verify_array_with_bounds(ar,bounds):
    '''Verify if values of array ar fit withind declared bounds, if too many/too few bounds do as much as can be done'''
    if (len(ar)<=len(bounds)):
        num=len(ar) # number of values to check
    elif (len(ar)>len(bounds)):
        num=len(bounds) # number of values to check  
    verification=True # initial value
    for i in range(num):
        # check if some of the values doesn't fit in the bounds, if so return False
        if (ar[i]<bounds[i][0] or ar[i]>bounds[i][1]):
            verification=False
            break       
        
    return verification
  
def lnprob(p,copied_obj,prior):
    '''Compute logarithmic probability for given parameters'''     
    # now we need to compute loglikelihood using the fortran code on new parameters, and then add it to lp
    #copied_obj=signalfit
    newparams=copied_obj.generate_newparams_for_mcmc(p)
    copied_obj.overwrite_params(newparams)
    if not (copied_obj.verify_params_with_bounds()):
        return -np.inf # if parameters are outside of bounds we return -infinity
    else:      

       # print(copied_obj.params.jitters[:3],copied_obj.params.GP_params[:4],copied_obj.params.planet_params[:14])        
       # print(copied_obj.use.use_jitters[:3])    
        flag=copied_obj.fitting(fileinput=False, filename='Kep_input', minimize_loglik=True, amoeba_starts=0, outputfiles=[0,1,0],return_flag=True)
        #print(copied_obj.fit_results.loglik, flag)         
        #print(copied_obj.fit_results.loglik, compute_loglik_TT(p,copied_obj) )#,copied_obj.params.GP_params[:4],copied_obj.params.planet_params[:14])   

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
                
                signalfit.gps[i].set_parameter_vector(np.array(list(map(np.log,np.concatenate((signalfit.params.GP_params,np.atleast_1d(signalfit.params.jitters[i])))))))
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

                copied_obj.gps[i].set_parameter_vector(np.array([np.log(copied_obj.params.GP_params[0]),np.log(copied_obj.params.GP_params[1]) ,np.log(copied_obj.params.GP_params[2]) ,np.log(copied_obj.params.GP_params[3]),copied_obj.params.jitters[i]]))       
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
           copied_obj.gps[i].set_parameter_vector(np.array([np.log(copied_obj.params.GP_params[0]),np.log(copied_obj.params.GP_params[1]) ,np.log(copied_obj.params.GP_params[2]) ,np.log(copied_obj.params.GP_params[3]),copied_obj.params.jitters[i]]))       
           S = S + copied_obj.gps[i].log_likelihood(copied_obj.fit_results.rv_model.o_c[copied_obj.fit_results.idset==i])
                      
       #print(copied_obj.fit_results.loglik, compute_loglik_TT(p,copied_obj),  S)#,copied_obj.params.GP_params[:4],copied_obj.params.planet_params[:14])   
       return pr.choose_prior(p,prior)+S    
     
        
def run_mcmc(signalfit,  Kbounds=[[0.0,100000.0]],Pbounds=[[0.0,100000.0]],ebounds=[[-0.99,0.99]],wbounds=[[-2.0*360.0, 2.0*360.0]],Mbounds=[[-2.0*360.0, 2.0*360.0]],ibounds=[[-2.0*180.0, 2.0*180.0]],capbounds=[[-2.0*360.0, 2.0*360.0]],offbounds=[[-100000.0,100000.0]],jitbounds=[[0.0,10000.0]],lintrbounds=[[-10.0,10.0]], GPbounds=[[0.0,100000.0]], stmassbounds=[[0.01,1000.0]], prior=0, samplesfile='', level=(100.0-68.3)/2.0, threads=1, doGP=False, gp_par=[10.0,10.0,10.0,1.0], use_gp_par=[False,False,False,False], save_means=False, fileoutput=False, save_sampler=False,burning_ph=20, mcmc_ph=20, **kwargs):      

    '''Performs MCMC and saves results'''  
    
    if threads == 'max':
        threads = multiprocessing.cpu_count()
    
    
    # Let's prepare things depending if we want to do GP or not
    if (doGP):
        signalfit.initiategps(gp_par=gp_par)
        fun=lnprobGP3
        signalfit.use.update_use_GP_params(use_gp_par)                            
        signalfit.params.update_GP_params(gp_par)
    else:
        fun=lnprob
        signalfit.use.update_use_GP_params(use_gp_par)                            
        signalfit.params.update_GP_params(gp_par)
        GPbounds=[[x-10.0,x+10.0] for x in gp_par] # just to make sure GPbounds don't cause lnprob return -infinity when we don't do GP (now all GP_params will be within bounds for sure)
    
    # we will need this later
    if (signalfit.mod_dynamical):
        mod='dyn'
    else:
        mod='kep'

    signalfit.prepare_for_mcmc(Kbounds=Kbounds,Pbounds=Pbounds,ebounds=ebounds,wbounds=wbounds,Mbounds=Mbounds,ibounds=ibounds,capbounds=capbounds,offbounds=offbounds,jitbounds=jitbounds,lintrbounds=lintrbounds, GPbounds=GPbounds, stmassbounds=stmassbounds)    
 
    start_time = time.time()


    ndim, nwalkers = len(signalfit.par_for_mcmc), len(signalfit.par_for_mcmc)*4

    #par_for_mcmc_ = self.par_for_mcmc
   #import pathos as pa 
   # pmult = pa.pools._ProcessPool(nodes=4)
    #from pathos.pools import ProcessPool

    #pmult = ProcessPool(nodes=4)
   # import pathos
   # pmult = pathos.multiprocessing.Pool(processes=2)    
    
    
    pos = [signalfit.par_for_mcmc + 1e-3*np.random.rand(ndim) for i in range(nwalkers)]

    sampler = CustomSampler(nwalkers, ndim, fun, args=[signalfit,prior], threads = threads)# pool = pmult)
 
    # burning phase
    pos, prob, state  = sampler.run_mcmc(pos,burning_ph)

    sampler.reset()
 
    # now perform the MCMC

    pos, prob, state  = sampler.run_mcmc(pos,mcmc_ph)
 
    ln = np.hstack(sampler.lnprobability)
    sampler.save_samples(signalfit.f_for_mcmc,signalfit.filelist.ndset,signalfit.npl)
          
    print("--- %s seconds ---" % (time.time() - start_time))    

#    pool.close()
    
    if (fileoutput):
        if (samplesfile==''): # that means no file name for samples file has been provided, so we generate a default one
            samplesfile='samples_%s'%mod
        outfile = open(samplesfile, 'w') # file to save samples

        for j in range(len(sampler.samples)):
            outfile.write("%s  " %(ln[j]))        
            for z in range(len(signalfit.par_for_mcmc)):
                outfile.write("%s  " %(sampler.samples[j,z]))
            outfile.write("\n")

        outfile.close()        
            
    # Now we will save new parameters and their errors (different + and - errors in this case). Flag save_means determines if we want to take means as new best fit parameters or stick to old ones and calculate errors with respect to that           
    if (save_means):
        signalfit.par_for_mcmc = sampler.means # we will not need to keep the old parameters in this attribbute, so let's store the means now
        
    new_par_errors = [[float(signalfit.par_for_mcmc[i] - np.percentile(sampler.samples[:,i], [level])),float(np.percentile(sampler.samples[:,i], [100.0-level])-signalfit.par_for_mcmc[i])] for i in range(len(signalfit.par_for_mcmc))] 
    newparams = signalfit.generate_newparams_for_mcmc(signalfit.par_for_mcmc)        
    print(newparams.GP_params)
    current_GP_params=newparams.GP_params # because calling fitting will overwrite them
    print(current_GP_params)

    signalfit.overwrite_params(newparams)
    signalfit.fitting(minimize_loglik=True, amoeba_starts=0, outputfiles=[1,1,1]) # this will help update some things 
    signalfit.update_with_mcmc_errors(new_par_errors)
    signalfit.params.update_GP_params(current_GP_params)

    if (doGP):
        signalfit.fitting_method='GP_%s'%mod
    else:
        signalfit.fitting_method='mcmc_%s'%mod  

    if (doGP):
        loglik_to_save = lnprobGP(signalfit.par_for_mcmc,signalfit,prior)
        signalfit.loglik=loglik_to_save         
        signalfit.fit_results.loglik=loglik_to_save   
         

    ###############  This is not working! you cannot save the sampler as an atribute and call it back later!
    ###############  See https://github.com/dfm/emcee/issues/148
    if(save_sampler):
        signalfit.sampler=sampler             
        signalfit.sampler_saved=True           
        
    #sampler.reset()

    return signalfit
    
    
    
    
    

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
        copied_obj = copy.deepcopy(obj)     
        index = planet - 1
        ############################################      
        ######### and here is the trick!  ##########
        ############################################        
        pp0 =  copied_obj.params.planet_params[7*index+0]  # we define a variable to be the planet amplitude Kj   
        #print(pp0)
        copied_obj.params.planet_params[7*index+0] = 0 # then we set Kj to be 0, i.e. remove the j-th planet signal
        copied_obj.fitting(fileinput=False, filename='Kep_input', minimize_loglik=True, amoeba_starts=0, outputfiles=[0,1,1],return_flag=False, npoints=int(len(obj.fit_results.model)), model_max=int(max(obj.fit_results.model_jd)-max(copied_obj.fit_results.rv_model.jd)))
        # and we create the static Nplanet model for the data and the model curve 
        # now this model residuals will contain ONLY the j-th planet signal + the best fit residuals
       
        copied_obj.params.planet_params[7*index+0] = pp0 # we restore Kj to its best fit value.
        ############################################      
        #########      trick is over      ##########
        ############################################  
        
    
        ############ phase fold fix for sparse model ######
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







class CustomSampler(emcee.EnsembleSampler):
 
    def unique_rows(self):
    
        '''
        Given an array, remove identical rows and also sort it
        '''
    
        # Perform lex sort and get sorted data
        sorted_idx = np.lexsort(self.flatchain.T)
        sorted_data =  self.flatchain[sorted_idx,:]

        # Get unique row mask
        row_mask = np.append([True],np.any(np.diff(sorted_data,axis=0),1))

        # Get unique rows 
        out = sorted_data[row_mask] 
        self.samples = out
        return

    def correct_rows(self,f,ndset,npl):
    
        '''Corrects angles and eccentricities for all samples'''     
     
        i=0
        self.means=np.zeros(len(f))
        for k in range(len(f)):
            idx=f[k]
            if (idx<2*ndset):          
                self.means[i]=np.mean(self.samples[:,i])         
            elif (idx<2*ndset+7*npl):
                nr=idx-2*ndset
                #x=int(nr/7)
                if (np.mod(nr,7)<2):
                    self.means[i]=np.mean(self.samples[:,i]) 
                elif (np.mod(nr,7)==2): # correct eccentricities
                    for j in range(len(self.samples)):
                        if (self.samples[j,i]<0):
                            self.samples[j,i]=abs(self.samples[j,i])
                            #if(f[k+1]==i+1):     
                            #    self.samples[j,i+1]=self.samples[j,i+1]+180.0
                            #if(f[k+2]==i+2):
                            #    self.samples[j,i+2]=self.samples[j,i+2]+-180.0
                    self.means[i]=np.mean(self.samples[:,i])      
                elif (np.mod(nr,7)==3): # correct w to be in a 360 interval around mean value 
                    self.means[i]=np.mean(self.samples[:,i]) 
                    meanw=self.means[i]                
                    for j in range(len(self.samples)):
                        self.samples[j,i]=np.where(self.samples[j,i]<meanw-180.0,self.samples[j,i]+360.0,self.samples[j,i])    
                        self.samples[j,i]=np.where(self.samples[j,i]>meanw+180.0,self.samples[j,i]-360.0,self.samples[j,i])
                    # now let's make sure meanw is between 0 and 360:
                    newmeanw=np.fmod(meanw,360.0)
                    delta=newmeanw-meanw
                    if not (delta==0):
                        for j in range(len(self.samples)):    
                            self.samples[j,i]=self.samples[j,i]+delta                
                elif (np.mod(nr,7)==4):# correct M to be in a 360 interval around mean value
                    self.means[i]=np.mean(self.samples[:,i]) 
                    meanM=self.means[i]               
                    for j in range(len(self.samples)):
                        self.samples[j,i]=np.where(self.samples[j,i]<meanM-180.0,self.samples[j,i]+360.0,self.samples[j,i])    
                        self.samples[j,i]=np.where(self.samples[j,i]>meanM+180.0,self.samples[j,i]-360.0,self.samples[j,i])
                    # now let's make sure meanw is between 0 and 360:
                    newmeanM=np.fmod(meanM,360.0)
                    delta=newmeanM-meanM
                    if not (delta==0):
                        for j in range(len(self.samples)):    
                            self.samples[j,i]=self.samples[j,i]+delta                
            elif (idx<2*ndset+6*npl):# correct i to be in a 180 interval around mean value
                self.means[i]=np.mean(self.samples[:,i])          
                meani=self.means[i]               
                for j in range(len(self.samples)):
                    self.samples[j,i]=np.where(self.samples[j,i]<meani-90.0,self.samples[j,i]+180.0,self.samples[j,i])    
                    self.samples[j,i]=np.where(self.samples[j,i]>meani+90.0,self.samples[j,i]-180.0,self.samples[j,i])
                # now let's make sure meani is between 0 and 180:
                newmeani=np.fmod(meani,180.0)
                delta=newmeani-meani
                if not (delta==0):
                    for j in range(len(self.samples)):    
                        self.samples[j,i]=self.samples[j,i]+delta                
            elif (idx<2*ndset+7*npl):# correct lineofnodes to be in a 360 interval around mean value 
                self.means[i]=np.mean(self.samples[:,i])
                meancap=self.means[i]              
                for j in range(len(self.samples)):
                    self.samples[j,i]=np.where(self.samples[j,i]<meancap-180.0,self.samples[j,i]+360.0,self.samples[j,i])    
                    self.samples[j,i]=np.where(self.samples[j,i]>meancap+180.0,self.samples[j,i]-360.0,self.samples[j,i])
                # now let's make sure meancap is between 0 and 360:
                newmeancap=np.fmod(meancap,360.0)
                delta=newmeancap-meancap
                if not (delta==0):
                    for j in range(len(self.samples)):    
                        self.samples[j,i]=self.samples[j,i]+delta      
            else:
                self.means[i]=np.mean(self.samples[:,i])           
            i=i+1                 
        return     
                
    def save_samples(self,f,ndset,npl):
        self.unique_rows()
        self.correct_rows(f,ndset,npl)
        return

class parameters(object): # class for all parameters which can be fitted

    def __init__(self,offsets,jitters,planet_params,linear_trend,stellar_mass):
        self.offsets=np.concatenate((np.atleast_1d(offsets),[0.0]*(max(0,20-len(np.atleast_1d(offsets)))))) 
        self.jitters=np.concatenate((np.atleast_1d(jitters),[0.0]*(max(0,20-len(np.atleast_1d(jitters))))))
        self.planet_params=np.concatenate((np.atleast_1d(planet_params),[0.0]*(max(0,70-len(np.atleast_1d(planet_params))))))
        self.linear_trend=linear_trend
        self.GP_params=[10.0]*4 # we always want to hav ethis attribute, but we only use it if we call GP, and then we update it anyway
        self.stellar_mass=stellar_mass
        
    def update_offset(self,dataset,offset): # change offset of one dataset
        self.offsets[dataset]=offset
        return
        
    def update_offsets(self,offsets): # change all offsets
        self.offsets=offsets
        return
        
    def update_jitter(self,dataset,jitter): # change jitter of one dataset
        self.jitters[dataset]=jitter
        return
        
    def update_jitters(self,jitters): # change all jitters 
        self.jitters=jitters
        return
        
    def update_K(self,planet,newK): # update K for a given planet
        self.planet_params[7*planet]=newK
        return
        
    def update_P(self,planet,newP): # update P for a given planet
        self.planet_params[7*planet+1]=newP
        return
        
    def update_e(self,planet,newe): # update e for a given planet
        self.planet_params[7*planet+2]=newe
        return
        
    def update_w(self,planet,neww): # update w for a given planet
        self.planet_params[7*planet+3]=neww
        return
        
    def update_M0(self,planet,newM0): # update M0 for a given planet
        self.planet_params[7*planet+4]=newM0
        return
        
    def update_inclination(self,planet,newi): # update inclination info for one planet
        self.planet_params[7*planet+5]=newi
        return              
        
    def update_lineofnodes(self,planet,newcap): # update lineofnodes info for one planet
        self.planet_params[7*planet+6]=newcap
        return              

    def update_planet_params_one_planet(self,planet,K,P,e,w,M0,i,cap): # update all planet_params for one planet
        self.update_K(planet,K)        
        self.update_P(planet,P)
        self.update_e(planet,e)
        self.update_w(planet,w)
        self.update_M0(planet,M0)
        self.update_inclination(planet,i)
        self.update_lineofnodes(planet,cap)
        return                
                
    def update_planet_params(self,planet_params): # update all planet_params in one go
        self.planet_params=planet_params
        return
                
    def update_linear_trend(self,linear_trend): # update linear trend 
        self.linear_trend=linear_trend       
        return      
        
    def update_GP_param(self,i,newpar):
        self.GP_params[i]=newpar
        return
        
    def update_GP_params(self,newparams):
        self.GP_params=newparams
        return
        
    def update_stellar_mass(self,stellar_mass):
        self.stellar_mass=stellar_mass
        return
        
    def getinit_input(self, masses, semimajor, fileinput=False, filename='geninit_j_input'):
          
        '''Prepares input for a fortran code which calculates Jacobi coordinates based on orbital parameters'''
        
        if not (fileinput): # if we want to save input in a file we don't want this line in the input string    
            ppp = './geninit_j3_in_days << EOF\n'
        else:
            ppp = ''  
        
        npl=len(self.planet_params)/7         
        
        ppp+='1 \n%s \n%s \n1.d0 \npl.in\n'%(str(self.stellar_mass),str(npl))    
        
        for i in range(npl):
            ppp+='%s \n%s %s %s %s %s %s\n'%(str(masses[i]),str(semimajor[i]),str(self.planet_params[7*i+2]),str(self.planet_params[7*i+5]),str(self.planet_params[7*i+3]),str(self.planet_params[7*i+6]),str(self.planet_params[7*i+4]))    
        
        if not (fileinput):
            ppp+='EOF' # end of the command to run in the case of saving input directly
        else: # here's what we do if we want to generate a file as well 
            # first we save the ppp string in a file (by default 'Kep_input')
            file_getin = open(filename, 'wb')
            file_getin.write('%s'%ppp) 
            file_getin.close()
            # then we overwrite ppp with the command to pass this file as input for the fortran code
            ppp='./geninit_j3_in_days < %s'%(program,filename)
                      
        return ppp        
                      
class use_flags(object): # class for all use flags

    def __init__(self,use_offsets,use_jitters,use_planet_params,use_linear_trend,use_stellar_mass):
        self.use_offsets=use_offsets 
        self.use_jitters=use_jitters
        self.use_planet_params=use_planet_params
        self.use_linear_trend=use_linear_trend
        self.use_GP_params=[False]*4 
        self.use_stellar_mass=False
        self.use_offsets=np.concatenate((np.atleast_1d(self.use_offsets),[0.0]*(20-len(np.atleast_1d(self.use_offsets))))) 
        self.use_jitters=np.concatenate((np.atleast_1d(self.use_jitters),[0.0]*(20-len(np.atleast_1d(self.use_jitters)))))
        self.use_planet_params=np.concatenate((np.atleast_1d(self.use_planet_params),[0.0]*(70-len(np.atleast_1d(self.use_planet_params)))))
        
    def update_use_offset(self,dataset,use_offset): # change flag for offset of one dataset
        self.use_offsets[dataset]=use_offset
        return
        
    def update_use_offsets(self,use_offsets): # change flags for all offsets
        self.use_offsets=use_offsets
        return
        
    def update_use_jitter(self,dataset,use_jitter): # change flag for jitter of one dataset
        self.use_jitters[dataset]=use_jitter
        return
        
    def update_use_jitters(self,use_jitters): # change flags for all jitters 
        self.use_jitters=use_jitters
        return
        
    def update_use_K(self,planet,newflag): # update K flag for a given planet
        self.use_planet_params[7*planet]=newflag
        return
        
    def update_use_P(self,planet,newflag): # update P flag for a given planet
        self.use_planet_params[7*planet+1]=newflag
        return
        
    def update_use_e(self,planet,newflag): # update e flag for a given planet
        self.use_planet_params[7*planet+2]=newflag
        return
        
    def update_use_w(self,planet,newflag): # update w flag for a given planet
        self.use_planet_params[7*planet+3]=newflag
        return
        
    def update_use_M0(self,planet,newflag): # update M0 flag for a given planet
        self.use_planet_params[7*planet+4]=newflag
        return
                
    def update_use_planet_params(self,use_planet_params): # update all use_planet_params flags in one go
        self.use_planet_params=use_planet_params
        return
              
    def update_use_inclination(self,planet,newflag): # update use_inclination info for one planet
        self.use_planet_params[7*planet+5]=newflag
        return             

    def update_use_lineofnodes(self,planet,newflag): # update use_lineofnodes info for one planet
        self.use_planet_params[7*planet+6]=newflag
        return                  
        
    def update_use_planet_params_one_planet(self,planet,Kflag,Pflag,eflag,wflag,M0flag,iflag,lineofnodesflag): # update all use_planet_params flags for one planet
        self.update_use_K(planet,Kflag)        
        self.update_use_P(planet,Pflag)
        self.update_use_e(planet,eflag)
        self.update_use_w(planet,wflag)
        self.update_use_M0(planet,M0flag)
        self.update_use_inclination(planet,iflag)
        self.update_use_lineofnodes(planet,lineofnodesflag)        
        return                
                
    def update_use_linear_trend(self,use_linear_trend): # update use_linear_trend flag
        self.use_linear_trend=use_linear_trend       
        return   
        
    def update_use_GP_param(self,i,newflag):
        self.use_GP_params[i]=newflag
        return
        
    def update_use_GP_params(self,use_GP_params):
        self.use_GP_params=use_GP_params
        return             

    def update_use_stellar_mass(self,use_stellar_mass):
        self.use_stellar_mass=use_stellar_mass
        return
        
class parameter_errors(object): # Class for parameter errors. 

    def __init__(self,offset_errors,jitter_errors,planet_params_errors,linear_trend_error,stellar_mass_error):
        '''At initiation we provide single values for each error, and we extrapolate them into 2 element arrays corresponding to + and - errors, which are at this point equal. They can be updated with updating functions found below'''
        # allocating room for up to 20 datasets and 10 planets
        offset_errors=np.concatenate((np.atleast_1d(offset_errors),[0.0]*(max(0,20-len(np.atleast_1d(offset_errors)))))) 
        jitter_errors=np.concatenate((np.atleast_1d(jitter_errors),[0.0]*(max(0,20-len(np.atleast_1d(jitter_errors))))))
        planet_params_errors=np.concatenate((np.atleast_1d(planet_params_errors),[0.0]*max((0,70-len(np.atleast_1d(planet_params_errors))))))

        self.offset_errors=np.array([[offset_errors[i],offset_errors[i]] for i in range(len(offset_errors))])
        self.jitter_errors=np.array([[jitter_errors[i],jitter_errors[i]] for i in range(len(jitter_errors))])
        self.planet_params_errors=np.array([[planet_params_errors[i],planet_params_errors[i]] for i in range(len(planet_params_errors))])
        self.linear_trend_error=np.array([linear_trend_error,linear_trend_error])
        self.GP_params_errors=[[0.0,0.0]]*4
        self.stellar_mass_error=np.array([stellar_mass_error,stellar_mass_error])

    '''In all functions below 'error' should in fact be a [-error,+error] 2 element array'''

    def update_offset_error(self,dataset,offset_error): # change offset error of one dataset
        self.offset_errors[dataset]=offset_error
        return
        
    def update_offset_errors(self,offset_errors): # change all offset errors
        self.offset_errors=offset_errors
        return
        
    def update_jitter_error(self,dataset,jitter_error): # change jitter error for one dataset
        self.jitter_errors[dataset]=jitter_error
        return
        
    def update_jitter_errors(self,jitter_errors): # change all jitters 
        self.jitter_errors=jitter_errors
        return
        
    def update_Kerror(self,planet,newKerror): # update K error for a given planet
        self.planet_params_errors[7*planet]=newKerror
        return
        
    def update_Perror(self,planet,newPerror): # update P error for a given planet
        self.planet_params_errors[7*planet+1]=newPerror
        return
        
    def update_eerror(self,planet,neweerror): # update e error for a given planet
        self.planet_params_errors[7*planet+2]=neweerror
        return
        
    def update_werror(self,planet,newwerror): # update w error for a given planet
        self.planet_params_errors[7*planet+3]=newwerror
        return
        
    def update_M0error(self,planet,newM0error): # update M0 error for a given planet
        self.planet_params_errors[7*planet+4]=newM0error
        return
        
    def update_inclination_errors(self,planet,newierror): # update inclination error for one planet
        self.planet_params_errors[7*planet+5]=newierror
        return              

    def update_lineofnodes_error(self,planet,newcaperror): # update lineofnodes error for one planet
        self.planet_params_errors[7*planet+6]=newcaperror
        return              
        
    def update_planet_param_errors_one_planet(self,planet,Kerror,Perror,eerror,werror,M0error,ierror,lineofnodeserror): # update all planet_params_errors for one planet
        self.update_Kerror(planet,Kerror)        
        self.update_Perror(planet,Perror)
        self.update_eerror(planet,eerror)
        self.update_werror(planet,werror)
        self.update_M0error(planet,M0error)
        self.update_inclination_error(planet,ierror)
        self.update_lineofnodes_error(planet,lineofnodeserror)
        return                
                
    def update_planet_param_errors(self,planet_param_errors): # update all planet_param_errors in one go
        self.planet_param_errors=planet_param_errors
        return
                
    def update_linear_trend_error(self,linear_trend_error): # update linear trend error
        self.linear_trend_error=linear_trend_error       
        return    

    def update_GP_param_errors(self,i,newerror): 
        self.GP_params_errors[i]=newerror
        return

    def update_GP_params_errors(self,GP_params_errors):  
        self.GP_params_errors=GP_params_errors
        return
  
    def update_stellar_mass_error(self,stellar_mass_error):
        self.stellar_mass_error=stellar_mass_error
        return  

class parameter_bounds(object): # Class for parameter bounds. We do not need updating functions here, as updating is only done together with fitting, and this is resolved there

    def __init__(self,offset_bounds,jitter_bounds,planet_params_bounds,linear_trend_bounds,GP_params_bounds,stellar_mass_bounds):
        self.offset_bounds=offset_bounds 
        self.jitter_bounds=jitter_bounds
        self.planet_params_bounds=planet_params_bounds
        self.linear_trend_bounds=linear_trend_bounds
        self.GP_params_bounds=GP_params_bounds
        self.stellar_mass_bounds=stellar_mass_bounds

#class for an object storing all the information about a given planet. Warning: contains less information than original kernel from old versions of RV_mod! Full information is now only stored in signal_fit class object
class kernel(object):
    
    def __init__(self,stat=0, jd=0,rvs=0,rv_err=0,o_c=0, model=0,model_jd=0,npl=0,a=0,mass=0,idset=0,stat_array_saved=0,reduced_chi2=0,chi2=0,rms=0,loglik=0):
        self.stat = stat
        self.rv_model=rvmodel(jd,rvs,rv_err,o_c)
        self.jd = jd
        self.rvs = rvs
        self.rv_err = rv_err
        self.o_c = o_c
        self.model = model
        self.model_jd = model_jd
        self.npl = npl
        self.a = a
        self.mass = mass
        self.idset = idset
        self.stat_array_saved=stat_array_saved
        self.reduced_chi2=reduced_chi2
        self.chi2=chi2
        self.rms=rms
        self.loglik=loglik
        
        
        
    def plot_ker(self, filetosave='time_series.png', plot_title='RV time series', save = True, dpi=200, colors = ['b','r','g','m','c','y','k','w','0.5','pink','lightsalmon','indigo','olive','plum','mediumvioletred','chocolate','blurywood','rosybrown','mistyrose','beige','ivory'], plot_oc=True):

        '''Plotting RV time series'''
        # preparing plot properties        
        if (plot_oc):
            plt.figure(figsize=(8,6.5))
            gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
            gs.update(wspace=0.15)
            plt.subplots_adjust(hspace=0.005)
            ax1 = plt.subplot(gs[0,0])
            ax2 = plt.subplot(gs[1,0])
            ax1.set_title(plot_title)  
        else:
            plt.figure(figsize=(8,10.5))
            gs = gridspec.GridSpec(1, 1) 
            gs.update(wspace=0.15, hspace = 0.2)
            ax1 = plt.subplot(gs[0,0])            

        # setting limits on the ordinate axis to fit all the data with some free space above and below. If statement for the cases of wrong offset value so that all RVs have the same sign 
        model_RV_min=min(self.rv_model.rvs)
        model_RV_max=max(self.rv_model.rvs)
        if model_RV_min>0:
            ax1.set_ylim(model_RV_min*0.45,model_RV_max*1.55)
        elif model_RV_max<0:
            ax1.set_ylim(model_RV_min*1.55,model_RV_max*0.45)
        else: 
            ax1.set_ylim(model_RV_min*1.55,model_RV_max*1.55) 

        # horizontal line at zero level
        zerx = np.array(range(int(min(self.rv_model.jd))-100, int(max(self.rv_model.jd))+100,1))
        zer = [0]*zerx
        ax1.plot(zerx,zer,'--', linewidth=1.0, color='k',markersize = 3.9,  mew=0)
        if (plot_oc):
            ax2.plot(zerx,zer,'--', linewidth=1.0, color='k',markersize = 3.9,  mew=0  )
        
        # plotting every data point separately, so they can have diffent colors, depending on which dataset they come from (idset)
        for z in range(len(self.rv_model.jd)):
            ax1.errorbar(self.rv_model.jd[z],self.rv_model.rvs[z],yerr=self.rv_model.rv_err[z], fmt='o',color=colors[self.idset[z]])
            if (plot_oc):
                ax2.errorbar(self.rv_model.jd[z],self.rv_model.o_c[z],yerr=self.rv_model.rv_err[z], fmt='o',color=colors[self.idset[z]])

        if(len(self.model_jd)==0):  
            self.model_jd=self.rv_model.jd
            self.model=np.zeros(len(self.model_jd))  
    
        ax1.plot(self.model_jd,self.model, color = 'r', label='')
    
        ax1.set_ylabel('RV [m/s]',fontsize=20)
        if (plot_oc):
            ax2.set_ylabel(r'o--c [m/s]',fontsize=20)
            ax2.set_xlabel('JD',fontsize=20)

            ranges2 = max(abs(min(self.rv_model.o_c)), max(self.rv_model.o_c)) # Parameter on which we will base the ordinate range for the O-C plot

        # setting time range for the plot, with some empty space on the right, but none on the left 
        ax1.set_xlim(min(self.rv_model.jd),max(self.rv_model.jd)*1.1-min(self.rv_model.jd)*0.1)    
        if(plot_oc):
            ax2.set_xlim(min(self.rv_model.jd),max(self.rv_model.jd)*1.1-min(self.rv_model.jd)*0.1) 
            # ordinate range for O-C plot      
            ax2.set_ylim(ranges2*-1.55,ranges2*1.55)

            plt.setp(ax2.get_yticklabels(), fontsize=15,weight='bold') #tuning the fontsize of yticklabels

        if save: # save the plot if requested by user (save=True as kep_plot input, deafult option)
            plt.savefig(filetosave, dpi=dpi, bbox_inches='tight') 
        
            
        return

class rvmodel(object):
 
    def __init__(self,jd,rvs,rv_err,o_c):
        self.jd=jd
        self.rvs=rvs
        self.rv_err=rv_err
        self.o_c=o_c

    def plot_gls(self, sig=np.array([0.1,0.01,0.001]),ind_sig=2, filetosave='periodogram.png',linestyle = [':','--','-.','-']): # plot gls for RV and O-C, only if they exist
     
        plot_gls_warnings=Warning_log([],'Plotting gls')
        # first let's check if we have RV and O-C data
        if (len(self.rvs)==0):
            plot_gls_warnings.update_warning_list('No RV data to run gls, gls will not be run!')
        elif (len(self.o_c)==0):
            plot_gls_warnings.update_warning_list('No O-C data to run gls, gls will not be run!')
        else:
            # building objects to run gls on
            rv_data=signal_data(self.jd,self.rvs,self.rv_err,sig_for_gls=sig)
            oc_data=signal_data(self.jd,self.o_c,self.rv_err,sig_for_gls=sig)
            rv_data.gls(ind_sig=ind_sig)
            oc_data.gls(ind_sig=ind_sig)
            
            # preparing plot properties
            plt.rc('text',usetex= True)
            font = {'family' : 'serif','weight' : 'bold','size'   : 18,'serif':['Helvetica']}
            plt.rc('font', **font)    
            plt.figure(figsize=(8,6.5))
            gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1.5]) 
            gs.update(wspace=0.15)
            plt.subplots_adjust(hspace=0.005)
            ax1 = plt.subplot(gs[0,0]) # in ax1 we will plot the RV signal periodogram
            ax2 = plt.subplot(gs[1,0]) # in ax2 we will plot the O-C periodogram     

            ############### Ploting the best peaks ##########################         
            # loop to generate annotations next to highest peaks
            for l in range(min(3,rv_data.number_of_significant_peaks)): # highlight 3 highest peaks, unless fewer pass the significance criteria 
                ax1.annotate(r'%.2f d'%rv_data.best_per[l],(rv_data.best_per[l],rv_data.best_peaks[l]*0.95), fontsize=16, color = 'k')
  
            ax1.plot(rv_data.periods, rv_data.P_G, '-',color ='r', lw=1.3, label='generalized')
            ax2.plot(oc_data.periods, oc_data.P_G, '-',color ='r', lw=1.3, label='generalized')
 
            # limits for abssica axes, lower bound is always 1 day     
            ax1.set_xlim(1, rv_data.gls_range)
            ax2.set_xlim(1, rv_data.gls_range) 

            ############ plot the significance lines. ############

            j = 0 # iterator to go through the linestyle list

            xlim = (rv_data.periods[0], rv_data.periods[-1])

            for zi in rv_data.z:
                ax1.plot(xlim, (zi, zi), linestyle[j], color = '0.25', lw=1)
                j = j +1

            # logscale for abscissa axis only
            ax1.semilogx()
            ax2.semilogx() 
               
            ax2.set_xlabel('days', fontsize=22, color = 'k')
            ax1.set_ylabel(r'Power', fontsize=22, color = 'k')
          
            plt.savefig(filetosave, format='png', dpi=600 ) 
            plot_gls_warnings.print_warning_log()           
            return
            
# class for the stat array in kernel
class summary(object):
 
    def __init__(self,params,param_errors,dof=0):
        self.params=params
        self.param_errors=param_errors
        self.dof=dof

# Defining some custom exceptions and warnings

class Error(Exception):
    pass

class InputError(Error):
 
    def __init__(self, message):
        self.message = message

class FittingError(Error):
 
    def __init__(self, message):
        self.message = message
        
class Warning_log(object): # here we can store all warnings raised when calling a certain function, so we can print them out nicely 
 
    warning_count=0
    warning_messages=[]
    procedure_name='Your query' # what function's warnings are we storing, will be needed for printing warning log
 
    def __init__(self,warning_messages, procedure_name):
        self.warning_count=len(warning_messages)
        self.warning_messages=warning_messages
        self.procedure_name=procedure_name
 
    def update_warning_list(self,warning_message): # add a new warning message to the warning list 
        self.warning_messages.append(warning_message)
        self.warning_count = self.warning_count+1
        return
 
    def print_warning_log(self):
        if(self.warning_count==0): # do not print a warning log if there were no warnings
            pass
        else: # printing a custom warning log with a list of all warnings encountered when running the procedure
            print('')
            print('%s resulted in %d warnings:'%(self.procedure_name,self.warning_count))
            print('')
            for i in range(self.warning_count):
                print('%d %s'%(i+1,self.warning_messages[i]))
            print('') 
        return

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
        
        
    def gls(self, ind_sig=2): # ind_sig will be the index of the significance level which we assume to be sufficient to declare a planet candidate

        gls_warnings=Warning_log([],'Running gls')
        if not (is_int(str(ind_sig),bounded=[True,True],bounds=[0,len(self.sig)],equal=[True,False])):
            ind_sig=len(self.sig)-1
            gls_warnings.update_warning_list('ind_sig for gls must be a non-negative integer smaller than the size of sig_for_gls array! Default value of len(sig)-1 will be assumed.')

        ### Compute the Lomb-Scargle Periodogram
        try:

            self.gls_range = (float(max(self.jd))-float(min(self.jd)))*2 # range for gls, twice the time range     
            self.periods = np.linspace(1, self.gls_range, 2000) # abscissas for the period range
            #omega = TAU / self.periods # converting the periods array into the frequency range for the Lomb-Scargle Periodogram evaluation
            omega = 1 / self.periods # converting the periods array into the frequency range for the Lomb-Scargle Periodogram evaluation
             #omega = 1/ np.logspace(-0.05, 4, num=1000)
        
        
            RV_gls = gls.Gls((self.jd, self.rvs, self.rv_error), fast=True,  verbose=False, norm= "ZK",ofac=5, fbeg=omega[999], fend=omega[ 0],)
        
            #self.P_G, self.z = lomb_scargle(self.jd, self.rvs, self.rv_error, omega, generalized=True, significance=self.sig) # Lomb-Scargle for the RV signal    
            self.P_G = RV_gls.power 
            self.z = RV_gls.powerLevel(sig_for_gls)
            
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
        
# We also need a class for the fortran output. Originally fortran code outputted everything in several files, now we want it all in console output, this causes some problems which will be solved by applying functions of this class

class fortran_output(object):
    # it's better to give initial values for these attributes, just in case the user calls functions in the wrong order
    best_par=[]
    RV_kep=[]
    keplerian_fit=[]
    params=parameters([],[],[],0.0,0.0)
    param_errors=parameter_errors([],[],[],0.0,0.0)
    ndata_str=''
    mfit_str=''
    rms_str=''
    chi_str=''
    epoch_str=''
    masses=[]
    semiM=[]
    JD_model=[]
    model=[]    
    jd=[]
    o_c=[]
    rv_obs=[]
    rv_error=[]
    data_set=[]
    loglik=0 
    reduced_chi2=1
    chi2=1
    rms=1
    
    def __init__(self,text, npl, ndset, stellar_mass):
        self.text=text
        self.npl=npl
        self.ndset=ndset
        self.stellar_mass=stellar_mass
        
    def dismantle_RV_kep(self): # save each column of RV_kep in a separate attribute 
        a=np.transpose(np.array(self.RV_kep)) # so columns are rows now
        if (len(a)==6): # otherwise something is wrong or keplerian_fit is empty, so we do nothing
            # we need to remove lines in which one of these entries was wrong

            indices=np.arange(min(list(map(len,a)))) # we will be deleting some of them as we go along
            self.jd,wrong_indices=convert_array_to_float(a[0],save_wrong_lines=True) 
            indices=[item for item in indices if item not in wrong_indices] # set substraction to exclude indices where thins went wrong
            self.rv_obs,wrong_indices=convert_array_to_float(a[2],save_wrong_lines=True) 
            indices=[item for item in indices if item not in wrong_indices] # set substraction to exclude indices where thins went wrong
            self.o_c,wrong_indices=convert_array_to_float(a[3],save_wrong_lines=True) 
            indices=[item for item in indices if item not in wrong_indices] # set substraction to exclude indices where thins went wrong
            self.rv_error,wrong_indices=convert_array_to_float(a[4],save_wrong_lines=True) 
            indices=[item for item in indices if item not in wrong_indices] # set substraction to exclude indices where thins went wrong
            self.data_set,wrong_indices=convert_array_to_int(a[5],save_wrong_lines=True) 
            indices=[item for item in indices if item not in wrong_indices] # set substraction to exclude indices where thins went wrong
            # and now we filter out wrong lines
            self.jd=self.jd[indices]
            self.rv_obs=self.rv_obs[indices]
            self.o_c=self.o_c[indices]
            self.rv_error=self.rv_error[indices]
            self.data_set=self.data_set[indices]
            self.data_set=np.array([i-1 for i in self.data_set]) # cause fortran names datasets starting from 1, and python from 0
        
    def dismantle_keplerian_fit(self): # save each column of keplerian_fit in a separate attribute 
        a=np.transpose(np.array(self.keplerian_fit, dtype=float)) # so columns are rows now
        if (len(a)==2): # otherwise something is wrong or keplerian_fit is empty, so we do nothing
            self.JD_model=a[0]
            self.model=a[1]
        return    
        
    def sort_out_text (self): # initially self.text will contain all that would normally be in fortran output files, we need to separate it into parts
        # similar action to read_file_as_array_of_arrays, except on a string
        T=self.text
        l=len(T)
        fortran_output_warnings=Warning_log([],'Reading fortran output')
        # now a loop which will do the sorting out
        i=0
        while i<l:
            if(len(T[i])==0): # skip empty lines
                i=i+1
            elif(T[i][0]=='ITMAX'): # error in amoeba
                raise RuntimeError ('Runtime error in amoeba (fortran subroutine). Change initial parameters.')
            elif(T[i][0]=='loglik,'): # these two lines always appear in the output
                self.loglik=float(T[i+1][0]) 
                self.reduced_chi2=float(T[i+1][1])
                self.chi2=float(T[i+1][2])
                self.rms=float(T[i+1][3])
                i=i+2
            elif(T[i][0]=='Best-fit'): # characteristic opening for the best_par.out file
                i0=i # save the index where file begins
                while 1:
                    i=i+1
                    if (i>l-1):
                        raise Exception ('Wrong fortran output')
                    if(len(T[i])==0): # empty lines
                        raise Exception ('Empty line in best_par part of fortran output in a place where it is not allowed!')
                    if (T[i][0]=='Jupiter'): # fourth last line of best_par.out begins this way
                        i=i+4 # so this will be the first line that doesn't belong to this file
                        self.best_par=T[i0:i] # now we can save the part of T which corresponds to best_par
                        break # and leave the while 1 loop

            elif(len(T[i])==6): # other two options are RV.out and keplerian_fit.out, the first one has 6 columns
                i0=i # save the index where file begins
                while 1:
                    i=i+1
                    if (i==l): # that means we reached the end of the file, so let's save what we have and leave
                        self.RV_kep=T[i0:i] # now we can save the part of T which corresponds to RV
                        break # and leave the while 1 loop                      
                    if (len(T[i])!=6 or not is_float(T[i][0])): # that means we went outside of the part of the output corresponding to RV_kep.out, second condition just in case
                        self.RV_kep=T[i0:i] # now we can save the part of T which corresponds to RV
                        break # and leave the while 1 loop     
                    
            elif(len(T[i])==2): # last option is keplerian_fit.out, which has 2 columns
                i0=i # save the index where file begins
                while 1:
                    i=i+1
                    if (i==l): # that means we reached the end of the file, so let's save what we have and leave
                        self.keplerian_fit=T[i0:i] # now we can save the part of T which corresponds to fit
                        break # and leave the while 1 loop                      
                    if (len(T[i])!=2 or not is_float(T[i][0])): # that means we went outside of the part of the output corresponding to RV_kep.out, second condition just in case
                        self.keplerian_fit=T[i0:i] # now we can save the part of T which corresponds to fit
                        break # and leave the while 1 loop                         
            else: # just in case there is something wrong with fortran output, skip this line
                i=i+1
                fortran_output_warnings.update_warning_list('Wrong data format in line %d of fortran output, line skipped.'%i) #i, not i+1, becase we increase it above    
        fortran_output_warnings.print_warning_log()
        return
    
    # this function is a bit of a work around, it puts together the part of the output which corresponds to best_par.out
    def print_stat_array(self):
        string_to_be_printed=''
        for line in self.best_par:
            for element in line:
                string_to_be_printed+='%s '%element
            string_to_be_printed+='\n'
        print(string_to_be_printed)
        return    

    # this one on the other hand reads bes_par.out and saves fitted parameters in object attributes
    def save_stat_array(self):
        planet_params=[]
        planet_params_errors=[]
        offsets=[]
        offset_errors=[]
        jitters=[]
        jitter_errors=[]
        linear_trend=0.0
        linear_trend_error=0.0
        i=0
        l=len(self.best_par)
        if(l==0): # to know that if there was no best_par in the output we do not really save anything here
            self.stat_array_saved=False
        else:
            self.stat_array_saved=True
        fortran_stat_warnings=Warning_log([],'Saving best fit parameters')
        while i<l:
            if (self.best_par[i][1]=='K'):# characteristic of beginning of planet parameters and errors information 
                for k in range(self.npl):
                    planet_params  = np.concatenate((planet_params,np.array(list(map(float,self.best_par[i+1+(k*2)][:7])))))
                    planet_params_errors  = np.concatenate((planet_params_errors,np.array(list(map(float,self.best_par[i+2+(k*2)][:7])))))
                i=i+2*self.npl+1
            elif (self.best_par[i][1]=='V0'):# characteristic of beginning of RV offsets and their errors information 
                for k in range(self.ndset):
                    offsets  = np.concatenate((offsets,np.array(list(map(float,self.best_par[i+1+(k*2)])))))
                    offset_errors  = np.concatenate((offset_errors,np.array(list(map(float,self.best_par[i+2+(k*2)])))))
                i=i+2*self.ndset+1
            elif (self.best_par[i][0]=='Jitters'):# characteristic of beginning of RV offsets and their errors information 
                for k in range(self.ndset):
                    jitters  = np.concatenate((jitters,np.array(list(map(float,self.best_par[i+1+(k*2)])))))
                    jitter_errors  = np.concatenate((jitter_errors,np.array(list(map(float,self.best_par[i+2+(k*2)])))))
                i=i+2*self.ndset+1
            elif (self.best_par[i][0]=='linear'): # characteristic of beginning of linear trend and it's error information
                linear_trend  = float(self.best_par[i+1][0])
                linear_trend_error  = float(self.best_par[i+2][0])
                i=i+3
            elif (self.best_par[i][0]=='ndata'): # rest of the information is below that
                self.ndata_str = self.best_par[i]
                self.mfit_str = self.best_par[i+1]
                self.rms_str   = self.best_par[i+2]
                self.chi_str   = self.best_par[i+3]
                self.epoch_str = self.best_par[i+4]
                self.masses  = map(float,self.best_par[i+6])
                self.semiM = map(float,self.best_par[i+8])
                i=i+9
            else:
                i=i+1
                fortran_stat_warnings.update_warning_list('Wrong data format in line %d of best fit parameter information in fortran output, line skipped. Please check if number of planets and number of datasets is specified correctly!'%i) #i, not i+1, becase we increase it above    
        self.params=parameters(offsets,jitters,planet_params,linear_trend,self.stellar_mass)
        self.param_errors=parameter_errors(offset_errors,jitter_errors,planet_params_errors,linear_trend_error,0.0)
        fortran_stat_warnings.print_warning_log()
        return
     
    def generate_summary(self):
        if(len(self.ndata_str)>2 and len(self.mfit_str)>2):
            dof = float(self.ndata_str[2]) - float(self.mfit_str[2]) # degrees of freedom for summary 
        else:
            dof= 1 # we have to give something y'know
        stat = summary(self.params,self.param_errors,dof=dof)
        return stat
       
    def modfit(self, print_stat=True):
        self.sort_out_text()
        if(print_stat): # if requested by user (passing print_stat=True to modfit) prints out the output file of the fortran code describing best fit parameters (by default True)
            self.print_stat_array()
        self.save_stat_array()
        self.dismantle_keplerian_fit()
        self.dismantle_RV_kep()
        results = kernel(self.generate_summary(), self.jd, self.rv_obs, self.rv_error,self.o_c, self.model, self.JD_model, self.npl,self.semiM,self.masses,self.data_set,self.stat_array_saved,self.reduced_chi2,self.chi2,self.rms,self.loglik) 
        return results
      

          
class RotationTerm(terms.Term):
    parameter_names = ("log_amp", "log_timescale", "log_period", "log_factor")

    def get_real_coefficients(self, params):
        log_amp, log_timescale, log_period, log_factor = params
        f = np.exp(log_factor)
        return (np.exp(log_amp) * (1.0 + f) / (2.0 + f), np.exp(-log_timescale))   


    def get_complex_coefficients(self, params):
        log_amp, log_timescale, log_period, log_factor = params
        f = np.exp(log_factor)
        return (np.exp(log_amp) / (2.0 + f), 0.0, np.exp(-log_timescale), 2*np.pi*np.exp(-log_period))   
                
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
        self.use=use_flags([False]*20,[False]*20,[False]*70,False,False) 
        self.params=parameters([0.0]*20,[0.0]*20,[0.0]*70,0.0,DEFAULT_STELLAR_MASS) 
        self.param_errors=parameter_errors([0.0]*20,[0.0]*20,[0.0]*70,0.0,0.0) 
        self.bounds = parameter_bounds([0.0,0.0]*20,[0.0,0.0]*20,[0.0,0.0]*70,[0.0,0.0],[0.0,0.0]*4,[0.0,0.0])  
        self.fit_performed = False
        self.fitting_method = 'None'
        self.model_saved=False
        self.stat_saved=False      
        self.masses=[0.0]*10
        self.semimajor=[0.0]*10
        self.f_for_mcmc=[]
        self.par_for_mcmc=[]
        self.e_for_mcmc=[]        
        self.gps=[]
        if (readinputfile):
            self.read_input_file() # after running this a file should be processed and initial values saved in the object properties, unless there was some error
            self.correct_elements() # to make sure angles are between 0 and 360 etc. 
            self.prepare_for_mcmc() # this will initiate default parameter bounds, in case the user wants to verify if values make sense he will be able to use the verify_params_with_bounds function anytime 
            self.inputfile_read=True
        # in this case we need to create a new kernel, but we will only give it this information which is needed for plotting
        self.filelist.read_rvfiles(self.params.offsets)
        self.fit_results=kernel(summary(parameters(0.0,0.0,0.0,0.0,DEFAULT_STELLAR_MASS),parameter_errors([0],[0],[0],0,0)), self.filelist.time,self.filelist.rvs,self.filelist.rv_err,np.zeros(len(self.filelist.time)),np.zeros(len(self.filelist.time)),self.filelist.time,0,0,0,self.filelist.idset,0,0,0,0,0)
        self.loglik=0.0
        self.rms=0.0
        self.chi2=0.0
        self.reduced_chi2=0.0
        self.sampler=None
        self.sampler_saved=False
        
        ##### this is how I wanted the kernel parameters to be 
        ##### initially defined: as a python dictionary and only to filled in by functions! 
        ##### It will be a lot of work to change those above now, but it has to be done as some point.
        
        self.evol_T =   {k: [] for k in range(9)}
        self.evol_a =   {k: [] for k in range(9)}
        self.evol_e =   {k: [] for k in range(9)}
        self.evol_p =   {k: [] for k in range(9)}
        self.evol_M =   {k: [] for k in range(9)}
         
    def update_epoch(self,epoch):
        self.epoch=epoch
        return

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
        
    def print_info (self, short_errors=True): 
        self.sort_by_period(reverse=False) # so the innermost planets will be first
        if(self.inputfile_read):
            print ('\nThis is the information for signal fit %s, based on input file %s. '%(self.name, self.inputfile), end="") 
        else:
            print ('\nThis is the information for signal fit %s.'%self.name, end="")            
        # checking if we are dealing with fitted parameters or original data
        if (self.fit_performed and not self.never_saved):
            print ('Presenting optimal parameters fitted using %s fitting method.'%self.fitting_method)
            print ('\n Fit properties: \n chi^2: %f \n reduced chi^2: %f \n rms: %f \n loglik: %f'%(self.fit_results.chi2,self.fit_results.reduced_chi2,self.fit_results.rms,self.fit_results.loglik))
        else:
            print ('No fit has yet been conducted (or parameters have never been saved), presenting parameters from user\'s original input')
        print(' ')  
        
        # Now let's print information about RV files
        if (self.filelist.ndset==1): # because word file has a singular and plural form
            print('RV signal data was provided in 1 file. We expect it to have following offset and jitter:')
        else:
            print('RV signal data was provided in %d files. We expect these files to have following offsets and jitters:'%self.filelist.ndset)
        if ((not (self.fit_performed) or self.fitting_method.startswith('loglik')) or self.never_saved): # if there was no fitting done (or values weren't saved) we don't give the errors. Same if the fitting was done with the Symplex method (if loglik was minimized), this method doesn't provide the errors
            for i in range(self.filelist.ndset): # neatly print all names of RV data files with corresponding offsets and jitters
                print('\n {0:15} \n offset: {1:>7.4f}  \n jitter: {2:>7.4f}'.format(self.filelist.files[i].name,self.params.offsets[i],self.params.jitters[i]))
            print('')  
        else: # if the fit was done we provide the errors
            if (short_errors): # assume same + and - error, equal to maximum of the two
                for i in range(self.filelist.ndset): # neatly print all names of RV data files with corresponding offsets and jitters
                    print('\n {0:15} \n offset: {1:>7.4f} +/- {3:>7.4f} \n jitter: {2:>7.4f} +/- {4:>7.4f}'.format(self.filelist.files[i].name,self.params.offsets[i],self.params.jitters[i],max(self.param_errors.offset_errors[i]),max(self.param_errors.jitter_errors[i])))
                print('')  
            else:
                for i in range(self.filelist.ndset): # neatly print all names of RV data files with corresponding offsets and jitters
                    print('\n {0:15} \n offset: {1:>7.4f} + {3:>7.4f} - {5:>7.4f} \n jitter: {2:>7.4f} + {4:>7.4f} - {6:>7.4f}'.format(self.filelist.files[i].name,self.params.offsets[i],self.params.jitters[i],self.param_errors.offset_errors[i][1],self.param_errors.jitter_errors[i][1],self.param_errors.offset_errors[i][0],self.param_errors.jitter_errors[i][0]))                 

        # Printing information about stellar activity, if fitting was done with GP
        if(self.fitting_method.startswith('GP')):
            print('\nStellar activity was modelled using Gaussian Processes. The resulting parameters are as follows:')
            if(short_errors):
                print('\n A = {0:>7.4f} +/- {4:>7.4f}\n t = {1:>7.4f} +/- {5:>7.4f}\n P = {2:>7.4f} +/- {6:>7.4f}\n f = {3:>7.4f} +/- {7:>7.4f}'.format(self.params.GP_params[0],self.params.GP_params[1],self.params.GP_params[2],self.params.GP_params[3],max(self.param_errors.GP_params_errors[0]),max(self.param_errors.GP_params_errors[1]),max(self.param_errors.GP_params_errors[2]),max(self.param_errors.GP_params_errors[3])))           
            else:
                print('\n A = {0:>7.4f} + {4:>7.4f} - {8:>7.4f}\n t = {1:>7.4f} + {5:>7.4f} - {9:>7.4f}\n P = {2:>7.4f} + {6:>7.4f} - {10:>7.4f}\n f = {3:>7.4f} + {7:>7.4f} - {11:>7.4f}'.format(self.params.GP_params[0],self.params.GP_params[1],self.params.GP_params[2],self.params.GP_params[3],self.param_errors.GP_params_errors[0][1],self.param_errors.GP_params_errors[1][1],self.param_errors.GP_params_errors[2][1],self.param_errors.GP_params_errors[3][1],self.param_errors.GP_params_errors[0][0],self.param_errors.GP_params_errors[1][0],self.param_errors.GP_params_errors[2][0],self.param_errors.GP_params_errors[3][0]))                                                                                  
                    
        # Printing information about linear trend, if any
        if (self.params.linear_trend!=0): # Information about a linear trend only if it exists
            if ((not (self.fit_performed) or self.fitting_method.startswith('loglik')) or self.never_saved): # if there was no fitting done (or values weren't saved) we don't give the errors. Same if the fitting was done with the Symplex method (if loglik was minimized), this method doesn't provide the errors
                print('There is a non-zero linear trend in the data, linear coefficient is expected to be equal %f.'%self.params.linear_trend)
            else:
                if (short_errors):
                    print('There is a non-zero linear trend in the data, linear coefficient is expected to be equal {0:>7.4f} +/- {1:>7.4f}.'.format(self.params.linear_trend,max(self.param_errors.linear_trend_error)))
                else:
                    print('There is a non-zero linear trend in the data, linear coefficient is expected to be equal {0:>7.4f} + {1:>7.4f} - {2:>7.4f}.'.format(self.params.linear_trend,self.param_errors.linear_trend_error[1],self.param_errors.linear_trend_error[0]))
            print('')
        
        # Printing information about stellar mass 
        if (not (self.fit_performed) or not ((self.fitting_method=='mcmc_dyn' or self.fitting_method=='GP_dyn') and (self.use.use_stellar_mass))):
            if (self.stellar_mass_known):
                print('The mass of the host star is {0:5.4f} solar masses.'.format(self.params.stellar_mass))
            else:
                print('The mass of the host star is not known. Default value of %f is assumed.'%DEFAULT_STELLAR_MASS)
        else:
            if(short_errors):  
                print('The mass of the host star is expected to be {0:5.4f} +/- {1:>7.4f} solar masses.'.format(self.params.stellar_mass,max(self.param_errors.stellar_mass_error)))      
            else:
                print('The mass of the host star is expected to be  {0:>7.4f} + {1:>7.4f} - {2:>7.4f}.'.format(self.params.stellar_mass,self.param_errors.stellar_mass_error[1],self.param_errors.stellar_mass_error[0]))

        print('')        
        
        # Printing planet parameters:
        if (self.npl==0): # if there are zero planets
            print('The system has no planets. ', end="")
        elif (self.npl==1): # because word planet has a singular and plural form
            print('The system has 1 planet. ', end="") 
        else:
            print('The system has %d planets. '%self.npl, end="")
        
        if (self.npl>0):
            if not (self.stat_saved):
                self.mass_semimajor() # update mass and semimajor axes values, if they weren't taken from the fortran code 
            print('Known planets are expected to have following properties (mean anomalies for epoch {0:7.2f}):'.format(self.epoch))
            if ((not (self.fit_performed) or self.fitting_method.startswith('loglik')) or self.never_saved): # again, if no fit or loglik fit, no errors
                for i in range(self.npl):
                    print('\n Planet {0:2d}: \n signal semiamplitude = {1:5.4f} m/s \n period = {2:5.4f} days \n orbit eccentricity = {3:5.4f} \n argument of periastron = {4:5.4f} \n mean anomally = {5:5.4f} \n inclination = {6:5.4f} \n line of nodes = {7:5.4f} \n mass  = {8:5.4f} M_Jup \n orbit semimajor axis = {9:5.4f} AU'.format(i+1,self.params.planet_params[7*i],self.params.planet_params[7*i+1],self.params.planet_params[7*i+2],self.params.planet_params[7*i+3],self.params.planet_params[7*i+4],self.params.planet_params[7*i+5],self.params.planet_params[7*i+6],self.masses[i],self.semimajor[i]))
            else:
                if (short_errors):
                    for i in range(self.npl):
                        print('\n Planet {0:2d}: \n signal semiamplitude = {1:5.4f} +/- {10:5.4f} m/s \n period = {2:5.4f} +/- {11:5.4f}  days \n orbit eccentricity = {3:5.4f} +/- {12:5.4f} \n argument of periastron = {4:5.4f} +/- {13:5.4f} \n mean anomally = {5:5.4f} +/- {14:5.4f} \n inclination = {6:5.4f} +/- {15:5.4f} \n line of nodes = {7:5.4f} +/- {16:5.4f} \n mass  = {8:5.4f} M_Jup \n orbit semimajor axis = {9:5.4f} AU'.format(i+1,self.params.planet_params[7*i],self.params.planet_params[7*i+1],self.params.planet_params[7*i+2],self.params.planet_params[7*i+3],self.params.planet_params[7*i+4],self.params.planet_params[7*i+5],self.params.planet_params[7*i+6],self.masses[i],self.semimajor[i],max(self.param_errors.planet_params_errors[7*i]),max(self.param_errors.planet_params_errors[7*i+1]),max(self.param_errors.planet_params_errors[7*i+2]),max(self.param_errors.planet_params_errors[7*i+3]),max(self.param_errors.planet_params_errors[7*i+4]),max(self.param_errors.planet_params_errors[7*i+5]),max(self.param_errors.planet_params_errors[7*i+6])))
                else:
                    for i in range(self.npl):
                        print('\n Planet {0:2d} \n signal semiamplitude = {1:5.4f} + {10:5.4f} - {17:5.4f} m/s \n period = {2:5.4f} + {11:5.4f} - {18:5.4f} days \n orbit eccentricity = {3:5.4f} + {12:5.4f} - {19:5.4f} \n argument of periastron = {4:5.4f} + {13:5.4f} - {20:5.4f} \n mean anomally = {5:5.4f} + {14:5.4f} -{21:5.4f} \n inclination = {6:5.4f} + {15:5.4f} - {22:5.4f} \n line of nodes = {7:5.4f} + {16:5.4f} - {23:5.4f} \n mass  = {8:5.4f} M_Jup \n orbit semimajor axis = {9:5.4f} AU'.format(i+1,self.params.planet_params[7*i],self.params.planet_params[7*i+1],self.params.planet_params[7*i+2],self.params.planet_params[7*i+3],self.params.planet_params[7*i+4],self.params.planet_params[7*i+5],self.params.planet_params[7*i+6],self.masses[i],self.semimajor[i],self.param_errors.planet_params_errors[7*i][1],self.param_errors.planet_params_errors[7*i+1][1],self.param_errors.planet_params_errors[7*i+2][1],self.param_errors.planet_params_errors[7*i+3][1],self.param_errors.planet_params_errors[7*i+4][1],self.param_errors.planet_params_errors[7*i+5][1],self.param_errors.planet_params_errors[7*i+6][1],self.param_errors.planet_params_errors[7*i][0],self.param_errors.planet_params_errors[7*i+1][0],self.param_errors.planet_params_errors[7*i+2][0],self.param_errors.planet_params_errors[7*i+3][0],self.param_errors.planet_params_errors[7*i+4][0],self.param_errors.planet_params_errors[7*i+5][0],self.param_errors.planet_params_errors[7*i+6][0]))
        print('') 
        return        
    
  
    def fortran_input(self, program='chi2_kep', fileinput=False, filename='Kep_input', amoeba_starts=1, outputfiles=[1,1,1],eps='1.0E-8',dt=864000, when_to_kill=300, npoints=50, model_max = 100): # generate input string for the fortran code, optionally as a file


        ### ppp will be the input string. Depending on fileinput parameter we either save it in a file or save it directly 
     
        if not (fileinput): # if we want to save input in a file we don't want this line in the input string    
            ppp = './%s << EOF\n'%program
        else:
            ppp = '' 
        ppp+= '%s %f %d %d %d %d\n%f %d %d %d \n%d\n'%(eps,dt,amoeba_starts,when_to_kill,npoints, model_max, self.params.stellar_mass,outputfiles[0], outputfiles[1],outputfiles[2],self.filelist.ndset) # first three lines of fortran input: precision and timestep for integration, stellar mass and number of datasets
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
                
    ### this function is a wrapper calling a fortran program to fit parameters in keplerian mode by minimizing chi^2   WORK IN PROGRESS ON THAT ONE! 
        
    def fitting(self, minimize_loglik=False, fileinput=False, filename='Kep_input', outputfiles=[1,1,1], amoeba_starts=1, eps=1, dt=1, fortran_kill=300, timeout_sec=600, print_stat=False, return_flag=False, npoints=1000, model_max = 500): # run the fit which will either minimize chi^2 or loglik.
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
        program='./fitting_routines/%s_%s'%(minimized_value,mod) 
        text,flag=run_command_with_timeout(self.fortran_input(program=program, fileinput=fileinput, filename=filename, outputfiles=outputfiles,amoeba_starts=amoeba_starts,eps=eps,dt=dt, when_to_kill=fortran_kill, npoints=npoints, model_max = model_max), timeout_sec, output=True,pipe=(not bool(outputfiles[2]))) # running command generated by the fortran_input function 
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
    def quick_overwrite_use_and_fit(self,useflags,minimize_loglik=False, fileinput=False, filename='Kep_input', amoeba_starts=1, outputfiles=[1,1,0], eps=1, dt=1, timeout_sec=600, print_stat=False):
        oldflags=self.overwrite_use(useflags,save=True)
        self.fitting(minimize_loglik=minimize_loglik,fileinput=fileinput,filename=filename,
        amoeba_starts=amoeba_starts,outputfiles=outputfiles, eps=eps,dt=dt,timeout_sec=timeout_sec,print_stat=print_stat )     
        self.overwrite_use(oldflags,save=False)
        return
        
    # if you want to temporarily change initial parameters and run a new fit, use this wrapper
    def quick_overwrite_params_and_fit(self,params,minimize_loglik=False, fileinput=False, filename='Kep_input', amoeba_starts=1, outputfiles=[1,1,0], eps=1, dt=1, timeout_sec=600, print_stat=False, return_flag=False):
        oldparams=self.overwrite_params(params,save=True)
        if (return_flag):
            flag=self.fitting(minimize_loglik=minimize_loglik,fileinput=fileinput,filename=filename,amoeba_starts=amoeba_starts,outputfiles=outputfiles,eps=eps,dt=dt,timeout_sec=timeout_sec,print_stat=print_stat, return_flag=True)     
            self.overwrite_params(oldparams,save=False)
            return flag
        else:
            self.fitting(minimize_loglik=minimize_loglik,fileinput=fileinput,filename=filename,amoeba_starts=amoeba_starts,outputfiles=outputfiles,eps=eps,dt=dt,timeout_sec=timeout_sec,print_stat=print_stat)     
            self.overwrite_params(oldparams,save=False)
            return                                   
            
    def minimize_one_param_K(self,planet,minimize_loglik=False, fileinput=False, filename='Kep_input', amoeba_starts=5, outputfiles=[1,0,0], eps=1, dt=1, timeout_sec=600, print_stat=False):
        useflags=use_flags([False]*20,[False]*20,[False]*70,False,False) 
        useflags.update_use_K(planet,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat)
        return    

    def minimize_one_param_P(self,planet,minimize_loglik=False, fileinput=False, filename='Kep_input', amoeba_starts=5, outputfiles=[1,0,0], eps=1, dt=1, timeout_sec=600, print_stat=False):
        useflags=use_flags([False]*20,[False]*20,[False]*70,False,False) 
        useflags.update_use_P(planet,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat)
        return    
        
    def minimize_one_param_e(self,planet,minimize_loglik=False, fileinput=False, filename='Kep_input', amoeba_starts=5, outputfiles=[1,0,0], eps='1.0E-8', dt=864000, timeout_sec=60, print_stat=False):
        useflags=use_flags([False]*20,[False]*20,[False]*70,False,False) 
        useflags.update_use_e(planet,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat)
        return    
                
    def minimize_one_param_w(self,planet,minimize_loglik=False, fileinput=False, filename='Kep_input', amoeba_starts=5, outputfiles=[1,0,0], eps=1, dt=1, timeout_sec=60, print_stat=False):
        useflags=use_flags([False]*20,[False]*20,[False]*70,False,False) 
        useflags.update_use_w(planet,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat)
        return    
        
    def minimize_one_param_M0(self,planet,minimize_loglik=False, fileinput=False, filename='Kep_input', amoeba_starts=5, outputfiles=[1,0,0], eps=1, dt=1, timeout_sec=60, print_stat=False):
        useflags=use_flags([False]*20,[False]*20,[False]*70,False,False) 
        useflags.update_use_M0(planet,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat)
        return    
        
    def minimize_one_param_inclination(self,planet,minimize_loglik=False, fileinput=False, filename='Kep_input', amoeba_starts=5, outputfiles=[1,0,0], eps=1, dt=1, timeout_sec=60, print_stat=False):
        useflags=use_flags([False]*20,[False]*20,[False]*70,False,False) 
        useflags.update_use_inclination(planet,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat)
        return    
        
    def minimize_one_param_lineofnodes(self,planet,minimize_loglik=False, fileinput=False, filename='Kep_input', amoeba_starts=5, outputfiles=[1,0,0], eps=1, dt=1, timeout_sec=60, print_stat=False):
        useflags=use_flags([False]*20,[False]*20,[False]*70,False,False) 
        useflags.update_use_lineofnodes(planet,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat)
        return    
             
    def minimize_one_param_offset(self,dataset,minimize_loglik=False, fileinput=False, filename='Kep_input', amoeba_starts=5, outputfiles=[1,0,0], eps=1, dt=1, timeout_sec=60, print_stat=False):
        useflags=use_flags([False]*20,[False]*20,[False]*70,False,False) 
        useflags.update_use_offset(dataset,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat)
        return                 
  
    def minimize_one_param_jitter(self,dataset,minimize_loglik=False, fileinput=False, filename='Kep_input', amoeba_starts=5, outputfiles=[1,0,0], eps=1, dt=1, timeout_sec=60, print_stat=False):
        useflags=use_flags([False]*20,[False]*20,[False]*70,False,False) 
        useflags.update_use_jitter(dataset,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat)
        return                   
              
    def minimize_one_param_linear_trend(self,minimize_loglik=False, fileinput=False, filename='Kep_input', amoeba_starts=5, outputfiles=[1,0,0], eps=1, dt=1, timeout_sec=60, print_stat=False):
        useflags=use_flags([False]*20,[False]*20,[False]*70,False,False) 
        useflags.update_use_linear_trend(True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat)
        return       
        
    def minimize_one_param_stellar_mass(self,dataset,minimize_loglik=False, fileinput=False, filename='Kep_input', amoeba_starts=5, outputfiles=[1,0,0], eps=1, dt=1, timeout_sec=60, print_stat=False):
        useflags=use_flags([False]*20,[False]*20,[False]*70,False,False) 
        useflags.update_use_stellar_mass(dataset,True)
        self.quick_overwrite_use_and_fit(useflags,minimize_loglik=minimize_loglik, fileinput=fileinput, filename=filename, amoeba_starts=amoeba_starts, outputfiles=outputfiles, eps=eps, dt=dt, timeout_sec=timeout_sec, print_stat=print_stat)
        return                               
              
    def plot_periodogram(self, filetosave='periodogram.png'):
        if(self.model_saved): # if we have the model we can run plot_gls on an appropriate kernel
            self.fit_results.rv_model.plot_gls(filetosave=filetosave)
        else: # otherwise we can only plot it for initial data
            object_for_gls=rvmodel(self.filelist.time,self.filelist.rvs,self.filelist.rv_err,[0.0]*len(self.filelist.time))
            object_for_gls.plot_gls(filetosave=filetosave)
        return
        
    def plot_time_series(self, filetosave='time_series.png', plot_title='RV time series', save = True, dpi=200):

        '''Generate a kernel object and plot it's time series'''
        
        if not (self.fit_performed): # to update in case we hadn't performed a fit, but added new datasets         
            self.fit_results=kernel(summary(parameters(0,0,0,0,0),parameter_errors([0],[0],[0],0,0)), self.filelist.time,self.filelist.rvs,self.filelist.rv_err,np.zeros(len(self.filelist.time)),np.zeros(len(self.filelist.time)),self.filelist.time,0,0,0,self.filelist.idset,0,0,0,0,0)
 
        self.fit_results.plot_ker(filetosave=filetosave, plot_title=plot_title, save=save, dpi=dpi)
        return
                   
    def plot_signal_one_planet(self,planet,filetosave='one_planet_signal.png', plot_title='One planet signal', save = True, dpi=200, warnings=None):
    
        '''Plots a signal for a chosen planet, understood as the overall signal minus contribution from all other planets'''     
        print_at_end=False

        if (warnings==None): # this maens we are running this function stand alone, rather than inside plot_signal_each_planet, so we need to create a new warnings object and print warnings at the end of this function
            warnings=Warning_log([],'Plotting one planet signal')
            print_at_end=True

        if (planet>self.npl):
            warnings.update_warning_list('Planet number %d higher than number of planets %d, cannot generate plot!'%(planet,self.npl))
        else:                 
            special_task_copy=copy.deepcopy(self) # create new signal_fit object, which will be the same except semiamplitude for chosen planet will be zero
            normal_copy=copy.deepcopy(self) # create new signal_fit object, which will be the same          
            flag=normal_copy.fitting(fileinput=True, filename='temp_test', minimize_loglik=True, amoeba_starts=0, outputfiles=[0,1,1],return_flag=True) # this should find signal from other planets
            if not (flag==1):
                warnings.update_warning_list('Failed to find signal from all planets!')
            else:
                special_task_copy.params.update_K(planet,0.0)    
                flag=special_task_copy.fitting(fileinput=True, filename='temp_test', minimize_loglik=True, amoeba_starts=0, outputfiles=[0,1,1],return_flag=True) # this should find signal from other planets
                if not (flag==1):
                    warnings.update_warning_list('Failed to find signal from all planets except planet %d.'%planet)
                else:
                    if not (len(self.fit_results.rv_model.rvs)==len(special_task_copy.fit_results.rv_model.rvs)):
                        warnings.update_warning_list('Number of data points is %d, but after removing %d planet the number is different, i. e. %d. Something must have gone wrong.'%(len(normal_copy.fit_results.rv_model.rvs),planet,len(special_task_copy.fit_results.rv_model.rvs)))
                    rvs_to_pass=[]        
                    for i in range(min(len(normal_copy.fit_results.rv_model.rvs),len(special_task_copy.fit_results.rv_model.rvs))):
                        rvs_to_pass.append(normal_copy.fit_results.rv_model.rvs[i]-special_task_copy.fit_results.rv_model.rvs[i])
                    rvs_to_pass=np.array(rvs_to_pass)
                    if not (len(normal_copy.fit_results.model)==len(special_task_copy.fit_results.model)):
                        warnings.update_warning_list('Number of points in model is %d, but after removing %d planet the number is different, i. e. %d. Something must have gone wrong.'%(len(normal_copy.fit_results.model),planet,len(special_task_copy.fit_results.model)))
                    model_to_pass=[]  
                    for i in range(min(len(normal_copy.fit_results.model),len(special_task_copy.fit_results.model))):          
                        model_to_pass.append(float(normal_copy.fit_results.model[i])-float(special_task_copy.fit_results.model[i]))
                    model_to_pass=np.array(model_to_pass)
                    special_task_kernel=kernel(normal_copy.fit_results.stat, normal_copy.fit_results.rv_model.jd,rvs_to_pass,normal_copy.fit_results.rv_model.rv_err,normal_copy.fit_results.rv_model.o_c,model_to_pass,normal_copy.fit_results.model_jd,normal_copy.npl,0.0,0.0,normal_copy.fit_results.idset,normal_copy.fit_results.stat_array_saved,0.0,0.0,0.0,0.0)
                    special_task_kernel.plot_ker(filetosave=filetosave, plot_title=plot_title, save=save, dpi=dpi, plot_oc=False)
        
        if (print_at_end):
            warnings.print_warning_log()
        
        return
       
    def plot_signal_each_planet(self,filestosave=[],plot_titles=[],dpi=200):
        
        planetsignalplotwarnings=Warning_log([],'Plotting individual planet signals')       
        
        if (len(filestosave)<self.npl): # this means the user didn't provide custom filenames for all the plots
            if (len(filestosave)>0): # so he wanted to provide custom names, but didn't provide for all
                planetsignalplotwarnings.update_warning_list('Too few file names for the plots, will assume default for remaining!')
            for i in np.arange(len(filestosave),self.npl,1):
                filestosave.append('signal_plot_pl_%d.png'%(i+1))
            filestosave=np.array(filestosave)        
        if (len(filestosave)>self.npl):
            planetsignalplotwarnings.update_warning_list('Too many file names for the plots, additional will be ignored!')
            filestosave=filestosave[:self.npl]
        if (len(plot_titles)<self.npl): # this means the user didn't provide custom plottitles for all the plots
            if (len(plot_titles)>0): # so he wanted to provide custom names, but didn't provide for all
                planetsignalplotwarnings.update_warning_list('Too few plot titles, will assume default for remaining!')
            for i in np.arange(len(plot_titles),self.npl,1):
                plot_titles.append('signal_pl_%d'%i)
            plot_titles=np.array(plot_titles)        
        if (len(plot_titles)>self.npl):
            planetsignalplotwarnings.update_warning_list('Too many plot tiltes, additional will be ignored!')
            plot_titles=plot_titles[:self.npl]
        
        # Done with the warnings, now let's plot
        for i in range(self.npl):
            self.plot_signal_one_planet(i,filetosave=filestosave[i], plot_title=plot_titles[i], save = True, dpi=dpi, warnings=planetsignalplotwarnings)
                    
        planetsignalplotwarnings.print_warning_log()
        return                   

    #### WORK IN PROGRESS HERE!
    def prepare_for_mcmc(self, Kbounds=[[0.0,100000.0]],Pbounds=[[0.0,100000.0]],ebounds=[[-0.99,0.99]],wbounds=[[-2.0*360.0, 2.0*360.0]],Mbounds=[[-2.0*360.0, 2.0*360.0]],ibounds=[[-2.0*180.0, 2.0*180.0]],capbounds=[[-2.0*360.0, 2.0*360.0]],offbounds=[[-100000.0,100000.0]],jitbounds=[[0.0,10000.0]],lintrbounds=[[-100000.0,100000.0]], GPbounds=[[0.0,100000.0]], stmassbounds=[[0.01,1000.0]], customdatasetlabels=[]):
 
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
        if(len(GPbounds)<4):
            if not (GPbounds==[[0.0,100000.0]]):
                preparingwarnings.update_warning_list('Too few GPbounds! Assuming default [0.0, 10000.0] for remaining parameters.')           
            GPbounds=np.concatenate((GPbounds,np.repeat([[0.0,10000.0]],4-len(GPbounds),axis=0)))           
        elif (len(GPbounds)>4):
            GPbounds=GPbounds[:4]
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

        par = np.concatenate((self.params.offsets[:self.filelist.ndset],self.params.jitters[:self.filelist.ndset],self.params.planet_params[:7*self.npl],np.atleast_1d(self.params.linear_trend),self.params.GP_params,np.atleast_1d(self.params.stellar_mass)))
           
        flag = np.concatenate((self.use.use_offsets[:self.filelist.ndset],self.use.use_jitters[:self.filelist.ndset],self.use.use_planet_params[:7*self.npl],np.atleast_1d(self.use.use_linear_trend),self.use.use_GP_params,np.atleast_1d(self.use.use_stellar_mass)))
        
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
        el_str=np.concatenate(([r'$\gamma_%s$ [m/s]'%customdatasetlabels[i] for i in range(self.filelist.ndset)],[r'jitt$_%s$ [m/s]'%customdatasetlabels[i] for i in range(self.filelist.ndset)],np.concatenate([[r'K$_%s$ [m/s]'%chr(98+i),r'P$_%s$ [day]'%chr(98+i) ,r'e$_%s$'%chr(98+i),r'$\omega_%s$ [deg]'%chr(98+i),r'M$_%s$ [deg]'%chr(98+i)] for i in range(self.npl)]),[r'i$_%s$ [deg]'%chr(98+i) for i in range(self.npl)],[r'$\Omega_%s$ [deg]'%chr(98+i) for i in range(self.npl)],np.atleast_1d(r'lin.trend [m/s/day]'),[r'Amp', r't', r'per', r'fact'],np.atleast_1d(r'st. mass [$M_\odot$]')))
   
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
        elif not (verify_array_with_bounds(self.params.GP_params,self.bounds.GP_params_bounds)):
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
            elif (idx<2*self.filelist.ndset+7*self.npl+5):
                newparams.update_GP_param(idx-2*self.filelist.ndset-7*self.npl-1,p[i])
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
                elif (idx<2*self.filelist.ndset+7*self.npl+5):
                    self.param_errors.update_GP_param_errors(idx-2*self.filelist.ndset-7*self.npl-1,p[i])
                    i=i+1
                else:    
                    self.param_errors.update_stellar_mass_error(p[i])            
        return                           

    def initiategps(self, gp_par=[10.0,10.0,10.0,1.0]): 
        
        # Prepare objects for Gaussian Processes        
        
        rot_kernel = terms.TermSum(RotationTerm(log_amp=np.log(gp_par[0]),log_timescale=np.log(gp_par[1]),log_period=np.log(gp_par[2]), log_factor=np.log(gp_par[3])))
        #rot_kernel = terms.SHOTerm(log_S0=np.log(gp_par[0]), log_Q=np.log(gp_par[1]), log_omega0=np.log(gp_par[2]))        
        
#        kernel_jitters=[]
        kernels=[]
        gps=[]
        #print(gp_par[0],gp_par[1],gp_par[2],gp_par[3] )

        for i in range(self.filelist.ndset):
            kernels.append(rot_kernel+terms.JitterTerm(np.log(self.params.jitters[i])))
            gps.append(celerite.GP(kernels[i], mean=0.0))
            gps[i].compute(self.filelist.time[self.filelist.idset==i],self.filelist.rv_err[self.filelist.idset==i])
            #gps[i].compute(self.filelist.time[self.filelist.idset==i])
        self.gps=gps
        return

    
    def mcmc(self,  Kbounds=[[0.0,100000.0]],Pbounds=[[0.0,100000.0]],ebounds=[[-0.99,0.99]],wbounds=[[-2.0*360.0, 2.0*360.0]],Mbounds=[[-2.0*360.0, 2.0*360.0]],ibounds=[[-2.0*180.0, 2.0*180.0]],capbounds=[[-2.0*360.0, 2.0*360.0]],offbounds=[[-100000.0,100000.0]],jitbounds=[[0.0,10000.0]],lintrbounds=[[-10.0,10.0]], GPbounds=[[0.0,100000.0]], stmassbounds=[[0.01,1000.0]], prior=0, samplesfile='', level=(100.0-68.3)/2.0, threads=1, doGP=False, gp_par=[10.0,10.0,10.0,1.0], use_gp_par=[False,False,False,False], save_means=False, fileoutput=False, save_sampler=False,burning_ph=100, mcmc_ph=1000, **kwargs):      

        '''Performs MCMC and saves results'''  
        
        if threads == 'max':
            threads = multiprocessing.cpu_count()
 
        
        # Let's prepare things depending if we want to do GP or not
        if (doGP):
            self.initiategps(gp_par=gp_par)
            fun=lnprobGP3
            self.use.update_use_GP_params(use_gp_par)                            
            self.params.update_GP_params(gp_par)
        else:
            fun=lnprob
            self.use.update_use_GP_params(use_gp_par)                            
            self.params.update_GP_params(gp_par)
            GPbounds=[[x-10.0,x+10.0] for x in gp_par] # just to make sure GPbounds don't cause lnprob return -infinity when we don't do GP (now all GP_params will be within bounds for sure)
        
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
            if (samplesfile==''): # that means no file name for samples file has been provided, so we generate a default one
                samplesfile='samples_%s'%mod
            outfile = open(samplesfile, 'w') # file to save samples

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
        current_GP_params=newparams.GP_params # because calling fitting will overwrite them
       # print(current_GP_params)

        self.overwrite_params(newparams)
        self.fitting(minimize_loglik=True, amoeba_starts=0, outputfiles=[1,1,1]) # this will help update some things 
        self.update_with_mcmc_errors(new_par_errors)
        self.params.update_GP_params(current_GP_params)

        if (doGP):
            self.fitting_method='GP_%s'%mod
        else:
            self.fitting_method='mcmc_%s'%mod  

        if (doGP):
            loglik_to_save = lnprobGP(self.par_for_mcmc,self,prior)
            self.loglik=loglik_to_save         
            self.fit_results.loglik=loglik_to_save   
             

        ###############  This is not working! you cannot save the sampler as an atribute and call it back later!
        ###############  See https://github.com/dfm/emcee/issues/148
        if(save_sampler):
            self.sampler=sampler             
            self.sampler_saved=True           
            
        #sampler.reset()

        return
             
    def cornerplot(self, cornerplotname='cornerplot.png', fileinput=False, filename='samples_kep'): 

        '''Generates a corner plot visualizing the mcmc samples. Optionally samples can be read from a file.'''

        if(fileinput):
            samples=read_file_as_array_of_arrays(filename)
        elif(self.sampler_saved):
            samples=self.sampler.samples
        else:
            raise Exception ('Please run mcmc and save sampler or provide a valid samples file!')
        fig = corner.corner(samples,bins=25, color="k", reverse=True, upper= True, labels=self.e_for_mcmc, quantiles=[level/100.0, 1.0-level/100.0],levels=(0.6827, 0.9545,0.9973), smooth=1.0, smooth1d=1.0, plot_contours= True, show_titles=True, truths=self.par_for_mcmc, dpi = 300, pad=15, labelpad = 50 ,truth_color ='r', title_kwargs={"fontsize": 12}, scale_hist=True,  no_fill_contours=True, plot_datapoints=True, kwargs=kwargs)
        fig.savefig(cornerplotname)  
 
        return   
        
    def run_stability_one_sample(self,sample_num,fileinput=False, filename='samples_kep', fileinputgetinit=False, filenamegetinit='geninit_j_input', warnings=None, timeout_sec=1000.0,timemax=3000.0, timestep=10):

        '''Runs stability analysis on one sample. Optionally samples can be read from a file.'''

        print_at_end=False

        if (warnings==None): # this maens we are running this function stand alone and no warnings object was passed, so we need to create a new warnings object and print warnings at the end of this function
            warnings=Warning_log([],'Stability analisys for sample %d'%sample_num)
            print_at_end=True

        if(fileinput):
            samples=read_file_as_array_of_arrays(filename)
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
        elif integrator=='mvs_GR':
            os.chdir('./stability/mvs_gr/')
    

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

        result, flag = run_command_with_timeout('./swift_symba5_j << EOF \nparam.in \npl.in \n1e-40 \nEOF', timeout_sec)                  
         
        
        for k in range(self.npl):
            result, flag = run_command_with_timeout('./follow_symba2 << EOF \nparam.in \npl.in \n%s \nEOF'%(k+2),timeout_sec)

            result, flag = run_command_with_timeout('mv follow_symba.out pl_%s.out'%(k+1),timeout_sec) 

            self.evol_T[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [0]) /  365.2425
            self.evol_a[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [2])
            self.evol_e[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [3])
            self.evol_p[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [6])      
            self.evol_M[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [7])


        
        os.chdir('../../')

        return   
        
              
