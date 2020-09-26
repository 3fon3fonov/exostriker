#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is one VERY simple example how to use RVmod as a fitting library. 
The RVmod options however, are a lot more. Reminder: the Exo-Striker 
GUI interface is warped around the RVmod. 


This example script delas with the Eta Ceti system (the usual demo in 
the Exo-Striker tool) and demonstrates how find planets automatically.
I.e. We only insert the data and we let RVmod to find the planets for us!

1. We add the RV data
2. We find the RV offsets
3. We run GLS to identify the significant peaks
4. We apply the auto-fitting routine
5. We make a simple time series plot showing the the Keplerian model find on the data
6. We plot the phased signals for each planet

More detailed examples of how to use the RVmod will be provided 
as Jupyter notebooks in future.


Created on Sun Jul  2 08:23:07 2019

@author: Trifon Trifonov
"""


import sys 
sys.path.append('../exostriker/lib/') #RV_mod directory must be in your path
import RV_mod as rv
import gls as gls
import numpy as np





# Lets create the RVmod object
fit=rv.signal_fit('Eta Ceti demo',readinputfile=False);

fit.cwd = '../exostriker/' # it is also important that the ES current working directory (cwd) point to the "lib" directory. This will be fixed in future releases 

# add the stellar mass
fit.params.stellar_mass = 1.7 # In M sol.

 
fit.add_dataset("hip5364_lick", "../exostriker/datafiles/hip5364.vels",0.0,0.0)  # the last two entries are initial offset and jitter
fit.add_dataset("hip5364_VLT", "../exostriker/datafiles/hip5364_crires.vels",0.0,0.0)

 
# Lets not fit for jitters now, i.e. keep at the initial value of 0 m/s
fit.use.use_jitters[0] = False
fit.use.use_jitters[1] = False


#  Run it once to find the RV offsets, no planets yet.
fit.fitting(outputfiles=[1,1,1], doGP=False,  minimize_fortran=True, minimize_loglik=False, amoeba_starts=20, print_stat=False)

#lets print the best fit params  
print("Loglik = %s"%fit.loglik)
fit.print_info() #this is an obsolete function call, will be replaced!




# Run GLS once 
rv.run_gls(fit)
rv.run_gls_o_c(fit)

# now lets find the planets in our data!
fit.auto_fit_max_pl = 2
fit.auto_fit_FAP_level = 0.001 # this corresponds to FAP = 0.1%. GLS power with FAP below that level we take as planet candidate.
fit.auto_fit_allow_ecc = True # otherwise will look only for circular planets

fit = rv.find_planets(fit)



 
#Lets print the best fit params  
print("Loglik = %s"%fit.loglik)
fit.print_info() #this is an obsolete function call, will be replaced!

 
# Lets fit one more time with RV jitters modeled
fit.use.use_jitters[0] = True
fit.use.use_jitters[1] = True

fit.fitting(minimize_loglik=True) # This will run the Fortran Simplex, which optimizes the lnL and thus can fit for the jitters (but no errors so far)
fit.fitting(minimize_loglik=False) # Now that the jiters are known we can get to the L-M method to gett Covar.Matrix par. uncertanties.


#Lets print the best fit params  
print("Loglik = %s"%fit.loglik)
fit.print_info() #this is an obsolete function call, will be replaced!

# To get the phased signals:
for i in range(fit.npl):
    rv.phase_RV_planet_signal(fit,i+1)   
 
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



################# Time series #####################



f = plt.figure(0, figsize=(8,6.5))
plt.subplots_adjust(hspace=0.005)
format_im = 'png' #'pdf'
dpi = 300

gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
#gs.update(  wspace=0.05)

ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])

 
color = ['b', 'r', 'g', 'r']
symbol = ['o', 'D', 'o', 'o']   
markersize = [5, 5, 6, 6] 
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

offset_pre  = 250
offset_post = 250

zero_point_T = range((int(min(jd))-offset_pre),(int(max(jd))+offset_post),10)
zero_point   = np.zeros(len(zero_point_T))


ax1.plot(kep_model_x, kep_model_y,       '-', linewidth=model_lw, color=model_color)
ax2.plot(zero_point_T,zero_point,'-', linewidth=model_lw, color=model_color)      

 

for i in range(len(data_set)):

        ax1.errorbar(jd[i],rvs[i], yerr=rv_err[i], alpha=alpha[int(data_set[i])], fmt=symbol[int(data_set[i])], linestyle='None', markersize = markersize[int(data_set[i])], color=color[int(data_set[i])], capsize = 0, elinewidth=1,mew=0.1)
        ax2.errorbar(jd[i],o_c[i], yerr=rv_err[i], alpha=alpha[int(data_set[i])], fmt=symbol[int(data_set[i])], linestyle='None', markersize = markersize[int(data_set[i])],color=color[int(data_set[i])], capsize = 0, elinewidth=1,mew=0.1)



 
ax1.set_ylabel(r'RV [m/s]',fontsize=16, rotation = 'vertical') 
ax1.set_xlim(min(jd)-offset_pre,max(jd)+offset_post)
 

ax2.set_xlabel(r'JD [day]',fontsize=16)
ax2.set_ylabel(r'o$-$c  [m/s]',fontsize=16, rotation = 'vertical') 
ax2.set_xlim(min(jd)-offset_pre,max(jd)+offset_post)

ax2.locator_params(axis="x", nbins=9)
plt.setp( ax2.get_yticklabels(), fontsize=15,weight='bold')
plt.setp( ax2.get_xticklabels(), fontsize=15,weight='bold')
 
# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
 
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False) 


plt.savefig('RV_plot_example_time_series.%s'%(format_im), format=format_im,dpi=dpi, bbox_inches='tight' )
ax1.cla() 
ax2.cla()










################## Phase plotting ###############

format_im = 'png' #'pdf'
dpi = 300

color = ['b', 'r', 'g', 'r']
symbol = ['o', 'D', 'o', 'o']   
markersize = [5, 5, 6, 6] 
alpha = [1, 1, 1, 1]

model_color = 'k'
model_lw = '1.0'

for j in range(fit.npl):

    f = plt.figure(1, figsize=(6,6))
 
    ax1 = plt.subplot(111)
 
    planet = j+1
    data, model = rv.phase_RV_planet_signal(fit,planet)

    jd = data[0]
    rvs = data[1]
    rv_err = data[2]
    data_set = data[3]

    ax1.plot(model[0],model[1], linestyle =  '-', linewidth=model_lw, color=model_color)

    # we can add the jitter
    add_jitter = True
    if add_jitter == True:
        rv_err = np.array([np.sqrt(rv_err[i]**2 + kep_fit.params.jitters[ii]**2)  for i,ii in enumerate(data_set)])

    for i in range(len(data_set)):
        ax1.errorbar(jd[i],rvs[i], yerr=rv_err[i], alpha=alpha[int(data_set[i])], fmt=symbol[int(data_set[i])], linestyle='None', markersize = markersize[int(data_set[i])], color=color[int(data_set[i])], capsize = 0, elinewidth=1,mew=0.1)
       

    ax1.set_xlabel(r'phase [day]',fontsize=16)
    ax1.set_ylabel(r'RV [m/s]',fontsize=16, rotation = 'vertical') 

    plt.savefig('RV_plot_example_phased_planet_%s.%s'%(j+1,format_im), format=format_im,dpi=dpi, bbox_inches='tight' )
    ax1.cla() 
 





 






















