#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is one VERY simple example how to use RVmod as a fitting library. 
The RVmod options however, are a lot more. Reminder: the Exo-Striker 
GUI interface is warped around the RVmod. 


This example script delas with the Eta Ceti system (the usual demo in 
the Exo-Striker tool) and demonstrates how to make .pdf plots using
the saves .ses files from the Exo-Striker and Matplotlib.


 
Created on 19 Feb 2020

@author: Trifon Trifonov
"""

import sys 
sys.path.append('../exostriker/lib/') #RV_mod directory must be in your path
import RV_mod as rv
import dill


 
################# Plotting #############################

import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl
import sys
import numpy as np

###### For nice plotting ##############

mpl.rcParams['axes.linewidth'] = 2.0 #set the value globally
mpl.rcParams['xtick.major.pad']='8'
mpl.rcParams['ytick.major.pad']='6'

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


#if sys.version_info[0] == 2:
#    file_ses = open(r"./Eta_Ceti_py2.ses", 'rb')
#else: 
 


file_ses = open(r"./Eta_Ceti_py3.ses", 'rb')

fit = dill.load(file_ses)
file_ses.close()      
 
fit.cwd = '../exostriker/' # it is also important that the ES current working directory (cwd) point to the "lib" directory. This will be fixed in future releases 
# add the stellar mass


#fit.params.stellar_mass = 1.7 # In M sol. 
 

add_jitter = True
pluss   =  130
minuss  =  130


################## Time series plotting ###############

f = plt.figure(0, figsize=(16,6.5))
plt.subplots_adjust(hspace=0.005)
format_im = 'pdf'
dpi = 300

gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
#gs.update(  wspace=0.05)

ax1 = plt.subplot(gs[:-1, -1])
ax2 = plt.subplot(gs[-1, -1])

 
color = ['#0280E8', 'r', 'g', 'r']
symbol = ['o', 'D', 'o', 'o']   
markersize = [6, 5, 6, 6] 
alpha = [1, 1, 1, 1]

model_color = '0.6'
model_lw = 2.0


#### Get the time series (these below are self explanatory) ########     
jd        = fit.fit_results.rv_model.jd 
rvs       = fit.fit_results.rv_model.rvs
rv_err    = fit.fit_results.rv_model.rv_err
o_c       = fit.fit_results.rv_model.o_c

data_set  = fit.filelist.idset

# we can add the jitter
 
if add_jitter == True:
    rv_err = np.array([np.sqrt(rv_err[i]**2 + fit.params.jitters[ii]**2)  for i,ii in enumerate(data_set)])



# Kep model time series #
kep_model_x = fit.fit_results.model_jd #-2458000.0
kep_model_y = fit.fit_results.model

print(len(kep_model_x))

# Dyn model time series #
#dyn_model_x = dyn_fit.fit_results.model_jd
#dyn_model_y = dyn_fit.fit_results.model

###################################################################



zero_point_T = range((int(min(jd))-minuss),(int(max(jd))+pluss),10)
zero_point   = np.zeros(len(zero_point_T))


ax1.plot(kep_model_x, kep_model_y,       '-', linewidth=model_lw, color=model_color)
ax2.plot(zero_point_T,zero_point,'-', linewidth=model_lw, color=model_color)      

 
for i in range(len(data_set)):
        ax1.errorbar(jd[i],rvs[i], yerr=rv_err[i], alpha=alpha[int(data_set[i])], fmt=symbol[int(data_set[i])], linestyle='None', markersize = markersize[int(data_set[i])], color=color[int(data_set[i])], capsize = 0, elinewidth=1,mew=0.1,zorder=10)
        ax2.errorbar(jd[i],o_c[i], yerr=rv_err[i], alpha=alpha[int(data_set[i])], fmt=symbol[int(data_set[i])], linestyle='None', markersize = markersize[int(data_set[i])],color=color[int(data_set[i])], capsize = 0, elinewidth=1,mew=0.1,zorder=10)
 

#legend =ax1.legend(frameon=True, loc='upper right', ncol=1,handlelength=1)


# Now add the legend with some customizations.
#legend = ax1.legend(loc='best', borderaxespad=0., shadow=False,numpoints=1)

# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
#frame = legend.get_frame()
#frame.set_facecolor('0.90')

# Set the fontsize
#for label in legend.get_texts():
#	label.set_fontsize('xx-small')

#for label in legend.get_lines():
#	label.set_linewidth(2.0)  # the legend line width

 
ax1.set_ylabel(r'RV [m/s]',fontsize=16, rotation = 'vertical') 
ax1.set_xlim(min(jd)-minuss,max(jd)+pluss)
 

ax2.set_xlabel(r'BJD [days]',fontsize=16)
ax2.set_ylabel(r'o$-$c  [m/s]',fontsize=16, rotation = 'vertical') 
ax2.set_xlim(min(jd)-minuss,max(jd)+pluss)

#ax1.set_ylim(-70,50)
#ax2.set_ylim(-15,15)

ax2.locator_params(axis="x", nbins=9)
plt.setp( ax2.get_yticklabels(), fontsize=15,weight='bold')
plt.setp( ax2.get_xticklabels(), fontsize=15,weight='bold')
 
# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
 
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False) 
ax1.tick_params(direction="in")
ax2.tick_params(direction="in")



#x1,x2,y1,y2 = ax1.axis() 
#ax2.annotate(r'wrms = 1.67\,ms$^{-1}$',(2456500.0,4), fontsize=18, color = 'k')


plt.savefig('RV_plot.%s'%(format_im), format=format_im,dpi=dpi, bbox_inches='tight' )
ax1.cla() 
ax2.cla()


 
################## Time series plotting 2 ###############

f = plt.figure(0, figsize=(12, 5))
plt.subplots_adjust(hspace=0.005)
format_im = 'pdf'
dpi = 300

 
#gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
#gs.update(  wspace=0.05)

ax1 = plt.subplot(111)
#ax2 = plt.subplot(gs[-1, -1])

 
color = ['#0280E8', 'r', 'g', 'r']
symbol = ['o', 'D', 'o', 'o']   
markersize = [6, 6, 6, 6] 
alpha = [1, 1, 1, 1]

model_color = '0.7'
model_lw = 1.3


#### Get the time series (these below are self explanatory) ########     
jd        = fit.fit_results.rv_model.jd
rvs       = fit.fit_results.rv_model.rvs
rv_err    = fit.fit_results.rv_model.rv_err
o_c       = fit.fit_results.rv_model.o_c

data_set  = fit.filelist.idset

# we can add the jitter
 
if add_jitter == True:
    rv_err = np.array([np.sqrt(rv_err[i]**2 + fit.params.jitters[ii]**2)  for i,ii in enumerate(data_set)])



# Kep model time series #
kep_model_x = fit.fit_results.model_jd
kep_model_y = fit.fit_results.model

 
###################################################################
 
zero_point_T = range((int(min(jd))-minuss),(int(max(jd))+pluss),10)
zero_point   = np.zeros(len(zero_point_T))


ax1.plot(kep_model_x, kep_model_y,       '-', linewidth=model_lw, color=model_color,zorder=-5)
 

for i in range(len(data_set)):
     ax1.errorbar(jd[i],rvs[i], yerr=rv_err[i], alpha=alpha[int(data_set[i])], fmt=symbol[int(data_set[i])], linestyle='None', markersize = markersize[int(data_set[i])], color=color[int(data_set[i])], capsize = 0, elinewidth=2,mew=0.1,zorder=10,label="")
 
 
 
#legend =ax1.legend(frameon=True, loc='upper right', ncol=1,handlelength=1)


# Now add the legend with some customizations.
#legend = ax1.legend(loc='best', borderaxespad=0., shadow=False,numpoints=1)

# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
#frame = legend.get_frame()
#frame.set_facecolor('0.90')

# Set the fontsize
#for label in legend.get_texts():
#	label.set_fontsize('xx-small')

#for label in legend.get_lines():
#	label.set_linewidth(2.0)  # the legend line width

 
ax1.set_ylabel(r'RV [m/s]',fontsize=20, rotation = 'vertical') 
ax1.set_xlim(min(jd)-minuss,max(jd)+pluss)
 

ax1.set_xlabel(r'BJD [days]',fontsize=20)
 
 
#ax1.set_ylim(-70,55)
 
ax1.locator_params(axis="x", nbins=9)
plt.setp( ax1.get_yticklabels(), fontsize=18,weight='bold')
plt.setp( ax1.get_xticklabels(), fontsize=18,weight='bold')
 
# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
 
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False) 
ax1.tick_params(direction="in")
 



#x1,x2,y1,y2 = ax1.axis() 
#ax2.annotate(r'wrms = 1.67\,ms$^{-1}$',(2456500.0,4), fontsize=18, color = 'k')


plt.savefig('RV_plot_ts.%s'%(format_im), format=format_im,dpi=dpi, bbox_inches='tight' )
ax1.cla()
 

################## Phase plotting ###############

format_im = 'pdf' #'pdf'
dpi = 300
 
model_color = '0.3'
model_lw = 1.5

 

f, ax = plt.subplots(1,2, figsize=(12, 5), facecolor='w', edgecolor='k',gridspec_kw = {'width_ratios':[1,1]})
f.subplots_adjust(hspace = .3, wspace=0.15)
ax = ax.ravel()

 
for j in range(fit.npl):

 
    planet = j+1
    data, model = rv.phase_RV_planet_signal(fit,planet)

    jd = data[0]
    rvs = data[1]
    rv_err = data[2]
    data_set = data[3]
    
    if j ==0:
        offset = -1
    elif j ==1:
        offset = -1

    model_time_phase = np.array((model[0]-offset)%fit.params.planet_params[7*(j)+1] )
                
    sort = sorted(range(len(model_time_phase)), key=lambda k: model_time_phase[k])                        
    model_time_phase  = model_time_phase[sort] 
    ph_model =  model[1][sort] 


    ax[j].plot(model_time_phase,ph_model, linestyle =  '-', linewidth=model_lw, color=model_color,zorder=-2)

    # we can add the jitter
 
    if add_jitter == True:
        rv_err = np.array([np.sqrt(rv_err[i]**2 + fit.params.jitters[ii]**2)  for i,ii in enumerate(data_set)])

    for i in range(len(data_set)):

        ax[j].errorbar((jd[i]-offset)%fit.params.planet_params[7*(j)+1],rvs[i], yerr=rv_err[i], alpha=alpha[int(data_set[i])], fmt=symbol[int(data_set[i])], linestyle='None', markersize = markersize[int(data_set[i])], color=color[int(data_set[i])], capsize = 0, elinewidth=1,mew=0.1)
 

    ax[j].set_xlim(min(model[0])-max(model[0])*0.02,max(model[0])+max(model[0])*0.02)

    ax[j].set_xlabel(r'phase [days]',fontsize=17)
    ax[j].tick_params(direction="in")
    ax[j].tick_params(axis='both', which='major', labelsize=14)    
    if j == 0:
        ax[j].set_xticks(np.arange(min(model[0]),max(model[0]),50))
        #ax[j].set_yticks(np.arange(-60,50,20))
        ax[j].set_ylabel(r'RV [m/s]',fontsize=17, rotation = 'vertical') 
       # ax[j].set_ylim(-58.5,50.5)
        #ax[j].annotate(r'GJ 1148 b',(0.8,8), fontsize=22)    
    elif j == 1:
        ax[j].set_xticks(np.arange(int(min(model[0])),max(model[0]),100))
       # ax[j].set_yticks(np.arange(-100,140,50))
        #ax[j].annotate(r'GJ 1148 c',(3.8,8), fontsize=22)  
       # ax[j].set_ylim(-28.5,28.5)  
plt.savefig('RV_plot_ph.%s'%(format_im), format=format_im,dpi=dpi, bbox_inches='tight' )
#ax.cla() 
 
  









