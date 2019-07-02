#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is one VERY simple example how to use RVmod as a fitting library. 
The RVmod options however, are a lot more and it can do much more than the
Exo-Striker GUI interface (which is warped around the RVmod). 


This example script delas with the Eta Ceti system (the usual demo in 
the Exo-Striker tool). 

1. We add the RV data
2. We find the offsets
3. We run GLS to identify the significant peaks
4. We apply an autofitting routine
5. We make a simple plot showing the the Keplerian model find on the data

 

More detailed examples of how to use the RVmod will be provided 
as Jupyter notebooks in future.


Created on Sun Jul  2 08:23:07 2019

@author: Trifon Trifonov
"""


import sys 
sys.path.append('./lib/') #RV_mod directory must be in your path
import RV_mod as rv
import gls as gls
import numpy as np

# Lets create the RVmod object
fit=rv.signal_fit('Eta Ceti demo',readinputfile=False);

fit.add_dataset("./datafiles/", "./datafiles/hip5364.vels",0.0,0.0)  # the last two entries are initial offset and jitter
fit.add_dataset("./datafiles/", "./datafiles/hip5364_crires.vels",0.0,0.0)

# add the stellar mass
fit.params.stellar_mass = 1.7 # In M sol.

# Lets not fit for jitters now, i.e. keep at the initial value of 0 m/s
fit.use.use_jitters[0] = False
fit.use.use_jitters[1] = False


#  Run it once to find the RV offsets, no planets yet.
fit.fitting(outputfiles=[1,1,1], doGP=False,  minimize_fortran=True, minimize_loglik=False, amoeba_starts=20, print_stat=False)

#lets print the best fit params  
print("Loglik = %s"%fit.loglik)
fit.print_info() #this is an obsolete function call, will be replaced!



def find_planets(obj):
 
    # check if RV data is present
    #if obj.filelist.ndset <= 0:  
   #      return        

    # the first one on the data GLS
    if obj.gls.power.max() <= obj.gls.powerLevel(0.001):                                                       
         return obj
    
    else:
        if obj.npl !=0:
            for j in range(obj.npl):
                obj.remove_planet(obj.npl-(j+1))

        mean_anomaly_from_gls = np.degrees((((obj.epoch - float(obj.gls.hpstat["T0"]) )% (obj.gls.hpstat["P"]) )/ (obj.gls.hpstat["P"]) ) * 2*np.pi)
         
        obj.add_planet(obj.gls.hpstat["amp"],obj.gls.hpstat["P"],0.0,0.0,mean_anomaly_from_gls -90.0,90.0,0.0)
        obj.use.update_use_planet_params_one_planet(0,True,True,False,False,True,False,False)     
               
        obj.fitting(fileinput=False,outputfiles=[1,1,1], doGP=False,   minimize_fortran=True, fortran_kill=3, timeout_sec= 3)
        obj = run_gls_o_c(obj)
        for i in range(obj.npl):
             rv.phase_RV_planet_signal(obj,i+1)          
        #now inspect the residuals
        
        for i in range(1,int(2)):
            
            if obj.gls_o_c.power.max() <= obj.gls_o_c.powerLevel(0.001):
                for j in range(obj.npl):
                    obj.use.update_use_planet_params_one_planet(j,True,True,False,False,True,False,False)     
           
                    obj.fitting(fileinput=False,outputfiles=[1,1,1], doGP=False,   minimize_fortran=True, fortran_kill=3, timeout_sec= 3)
                    obj = run_gls_o_c(obj)   
                return obj
            #elif (1/RV_per_res.hpstat["fbest"]) > 1.5:
            else:    
                mean_anomaly_from_gls = np.degrees((((obj.epoch - float(obj.gls_o_c.hpstat["T0"]) )% (obj.gls_o_c.hpstat["P"]) )/ (obj.gls_o_c.hpstat["P"]) ) * 2*np.pi)
         
                obj.add_planet(obj.gls_o_c.hpstat["amp"],obj.gls_o_c.hpstat["P"],0.0,0.0,mean_anomaly_from_gls -90.0,90.0,0.0)
                obj.use.update_use_planet_params_one_planet(i,True,True,False,False,True,False,False)  
                
                 
                obj.fitting(fileinput=False,outputfiles=[1,1,1], doGP=False,   minimize_fortran=True, fortran_kill=3, timeout_sec= 3)
                obj = run_gls_o_c(obj)

            #else:
             #   continue
                                   
        for j in range(obj.npl):
            obj.use.update_use_planet_params_one_planet(j,True,True,False,False,True,False,False)     

                  
        obj.fitting(fileinput=False,outputfiles=[1,1,1], doGP=False,   minimize_fortran=True, fortran_kill=3, timeout_sec= 3)
        obj = run_gls_o_c(obj) 
    return obj



def run_gls(obj):
             
    omega = 1/ np.logspace(np.log10(0.85), np.log10(5000), num=int(1000))
 
 

    if len(fit.fit_results.rv_model.jd) > 5:      
        RV_per = gls.Gls((obj.fit_results.rv_model.jd, obj.fit_results.rv_model.rvs, obj.fit_results.rv_model.rv_err), 
        fast=True,  verbose=False, norm='ZK',ofac=10, fbeg=omega[-1], fend=omega[0],)
        
        obj.gls = RV_per
    else:
        return obj
    
    return obj

def run_gls_o_c(obj):
                     
    omega = 1/ np.logspace(np.log10(0.85), np.log10(5000), num=int(1000))

 
    if len(obj.fit_results.rv_model.jd) > 5:
        RV_per_res = gls.Gls((obj.fit_results.rv_model.jd, obj.fit_results.rv_model.o_c, obj.fit_results.rv_model.rv_err), 
        fast=True,  verbose=False, norm='ZK', ofac=10, fbeg=omega[-1], fend=omega[ 0],)            

        obj.gls_o_c = RV_per_res        
    else:
        return obj

    return obj




 
 
fit = run_gls(fit)
fit = run_gls_o_c(fit)

 

fit = find_planets(fit)

 
 
 
 
#################  Dynamical fitting #############################
# now lets fit dynamics model starting from the best Keplerian derived above

#first lets copy the Keplerian object, we will need it later for plotting
import dill

kep_fit = dill.copy(fit)

 
 


 
 
################# Plotting #############################
 
# lets make some basic plots with the "fit" object results:

import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl

import numpy as np

###### For nice plotting ##############

mpl.rcParams['axes.linewidth'] = 2.0 #set the value globally
mpl.rcParams['xtick.major.pad']='8'
mpl.rcParams['ytick.major.pad']='2'

# set tick width
mpl.rcParams['xtick.major.size'] = 8
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['xtick.minor.size'] = 5
mpl.rcParams['xtick.minor.width'] = 2

mpl.rcParams['ytick.major.size'] = 8
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['ytick.minor.size'] = 5
mpl.rcParams['ytick.minor.width'] = 2


mpl.rc('text',usetex=True)
font = {'family' : 'normal','weight' : 'bold','size'   : 18,'serif':['Helvetica']}
mpl.rc('font', **font)
########################################



f = plt.figure(0, figsize=(8,6.5))
plt.subplots_adjust(hspace=0.005)
format_im = 'pdf'
dpi = 300

gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
#gs.update(  wspace=0.05)

ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])

 
color = ['b', 'r', 'g', 'r']
symbol = ['o', 'o', 'o', 'o']   
markersize = [6, 6, 6, 6] 
alpha = [1, 1, 1, 1]

model_color = 'k'
model_lw = '1.0'


#### Get the time series (these below are self explanatory) ########     
jd        = kep_fit.fit_results.rv_model.jd
rvs       = kep_fit.fit_results.rv_model.rvs
rv_err    = kep_fit.fit_results.rv_model.rv_err
o_c       = kep_fit.fit_results.rv_model.o_c

data_set  = kep_fit.filelist.idset

# we can add the jitter
add_jitter = True
if add_jitter == True:
    rv_err = np.array([np.sqrt(rv_err[i]**2 + kep_fit.params.jitters[ii]**2)  for i,ii in enumerate(data_set)])



# Kep model time series #
kep_model_x = kep_fit.fit_results.model_jd
kep_model_y = kep_fit.fit_results.model

 

###################################################################



zero_point_T = range((int(min(jd))-10),(int(max(jd))-10),10)
zero_point   = np.zeros(len(zero_point_T))


ax1.plot(kep_model_x, kep_model_y,       '-', linewidth=model_lw, color=model_color)
ax2.plot(zero_point_T,zero_point,'-', linewidth=model_lw, color=model_color)      

 

for i in range(len(data_set)):

        ax1.errorbar(jd[i],rvs[i], yerr=rv_err[i], alpha=alpha[int(data_set[i])], fmt=symbol[int(data_set[i])], linestyle='None', markersize = markersize[int(data_set[i])], color=color[int(data_set[i])], capsize = 0, elinewidth=1,mew=0.1)
        ax2.errorbar(jd[i],o_c[i], yerr=rv_err[i], alpha=alpha[int(data_set[i])], fmt=symbol[int(data_set[i])], linestyle='None', markersize = markersize[int(data_set[i])],color=color[int(data_set[i])], capsize = 0, elinewidth=1,mew=0.1)


 
ax1.set_ylabel(r'RV [m/s]',fontsize=16, rotation = 'vertical') 
ax1.set_xlim(min(jd),max(jd))
 

ax2.set_xlabel(r'JD - 2450000 [day]',fontsize=16)
ax2.set_ylabel(r'o$-$c  [m/s]',fontsize=16, rotation = 'vertical') 
ax2.set_xlim(min(jd),max(jd))

ax2.locator_params(axis="x", nbins=9)
plt.setp( ax2.get_yticklabels(), fontsize=15,weight='bold')
plt.setp( ax2.get_xticklabels(), fontsize=15,weight='bold')
 
# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
 
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False) 


plt.savefig('RV_plot_example_2.%s'%(format_im), format=format_im,dpi=dpi, bbox_inches='tight' )
ax1.cla() 











#####################


# Run MCMC 
#fit = rv.run_mcmc(fit, burning_ph=1000, mcmc_ph=5000, threads=30, output=False, fileoutput=True,save_means=False, save_mode=True, save_maxlnL=False)

# Run Nested sampling 

#fit = rv.run_nestsamp(fit, threads=30, std_output=False, stop_crit = 0.0001, Dynamic_nest = False, live_points = 500, fileoutput=True, save_means=False, save_mode=False, save_maxlnL=True)

# WARNING! setup the bounds/prioirs first. Usually these are wide open and if you dont set them up 
# it my take forever for the Nest. Samp. to finish. Unfortunatly I have to provide another example how to work with the RVmod priors/
# Work in progress....
      



#if you already have a session saved you may try:

#import dill

#file = open("session.ses", 'rb')
#fit = dill.load(file)
#file.close()     

# and then for example
#fit = rv.run_mcmc(fit, burning_ph=1000, mcmc_ph=5000, threads=30, output=False, fileoutput=True, save_means=False, save_mode=True, save_maxlnL=False)



























