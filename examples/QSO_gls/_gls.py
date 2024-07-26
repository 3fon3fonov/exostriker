#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is one VERY simple example how to use RVmod as a fitting library. 
The RVmod options however, are a lot more. Reminder: the Exo-Striker 
GUI interface is warped around the RVmod. 


This example script delas with the Eta Ceti system (the usual demo in 
the Exo-Striker tool) and demonstrates how to fit Keplerian 
and Dynamical models to RV data.


1. We add the RV data
2. We find the RV offsets
3. We apply approx. parameters (to be taken from a GLS, for example.)
4. We fit to get the best two-planet Keplerian model
5. We adopt the best Keplerian fit and we include the dynamics into the modeling.
6. We make a simple plot showing the deviation between Keplerian and N-body models.

There are some (commented) examples how one can run mcmc and/or nested sampling
to get errors and/or posterior distributions. 

More detailed examples of how to use the RVmod will be provided 
as Jupyter notebooks in future.


Created on Sun Jun  2 09:30:02 2019

@author: Trifon Trifonov
"""

import sys 
sys.path.append('../lib/') #RV_mod directory must be in your path
import RV_mod as rv
import dill
import gls as gls
 
################# Plotting #############################

import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl
from matplotlib import ticker
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


format_im = 'pdf'
dpi = 300

 
 
color = ['#0280E8', 'r', 'g', 'r']
symbol = ['o', 'D', 'o', 'o']   
markersize = [6, 5, 6, 6] 
alpha = [1, 1, 1, 1]

model_color = '0.6'
model_lw = 2.0

 
file_ses = open(r"./TIC358107516_py2_GLS.ses", 'rb')
fit = dill.load(file_ses)
file_ses.close()     
 


f, ax = plt.subplots(10,1, figsize=(11, 21), facecolor='w', edgecolor='k')
f.subplots_adjust(hspace = .03, wspace=0.4)

ax = ax.ravel()

 
color_gls = ['#0280E8', 'green', '#0280E8', 'green','#800080', '#808000', '#FF0000','#FFA500', '#FF4500', '#00BFFF', '#9ACD32','#FF00FF']
 
for i in range(10):


    if i == 0:
	    x = fit.fit_results.rv_model.jd
	    y = fit.fit_results.rv_model.rvs 
	    e_y = fit.fit_results.rv_model.rv_err  

    elif i == 1:
	    x = fit.fit_results.rv_model.jd
	    y = fit.fit_results.rv_model.o_c  
	    e_y = fit.fit_results.rv_model.rv_err 


    elif i > 1:
		act_data_set = fit.act_data_sets[i-2]
		y_act  = act_data_set[1] - np.mean(act_data_set[1])
 
	 
		x = act_data_set[0] 
		y = y_act   
		e_y = act_data_set[2]

 

    omega = 1/ np.logspace(np.log(0.85), np.log(5000), num=1000)
    power_levels = np.array([0.1000,0.01000,0.001000])

    act_per = gls.Gls((x,y,e_y), fast=True,  verbose=False, norm= "ZK",ofac=10, fbeg=omega[-1], fend=omega[ 0])
    ax2 = ax[i].twiny() 
    ax2.yaxis.set_major_locator(plt.MaxNLocator(3)) 

    factor = (len(x) -1.0) /2.0


 
    ax[i].plot([1/24.75, 1/24.75],[0.1,40], '--',color ='b',  alpha=0.7, lw=1.5, label='')
    ax[i].plot([1/11.91, 1/11.91],[0.1,40], '--',color ='b',  alpha=0.7, lw=1.5, label='') 
    #ax[i].plot([1/36.6, 1/36.6],[0.1,40], '--',color ='m',  alpha=0.7, lw=1.2, label='')     
    #if i >3:
   #     ax[i].plot([1/36.6, 1/36.6],[0.1,90], '--',color ='m',  alpha=0.7, lw=1.5, label='') 
   #     ax[i].plot([1/36.6, 1/36.6],[0.1,90], '--',color ='m',  alpha=0.7, lw=1.5, label='') 

    [ax[i].plot([max(act_per.freq),min(act_per.freq)],[fap*factor,fap*factor], linestyle='--',color='k') for ii,fap in enumerate(act_per.powerLevel(np.array(power_levels)))] 
    ax[i].plot(act_per.freq,act_per.power*factor, linestyle='-',color=color_gls[i], lw=2.2,)                    




    ax[i].set_xlim( 1/600.0,1/10.0)		
   # plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False) 
    ax[i].tick_params(direction="in")
    ax2.tick_params(direction="in") 

    if i == 5:
        ax[i].set_ylabel(r'Power',fontsize=24, rotation = 'vertical', labelpad = 15) 


    if i == 9:
        ax[i].set_xlabel(r'frequency [1/d]',fontsize=24) 
        ax2.set_xticklabels([])
    else:
        ax[i].set_xticklabels([])
        #ax2.set_xticklabels([])


	if i != 0:
		ax2.set_xticklabels([])
		ax2.set_xticks([])
	else:
		ax2.set_xlabel(r'period [d]',fontsize=24) 



		def tick_function(X):
		    V = 1/X
		    return ["%.1f" % z for z in V]
		#ax1Ticks = ax[pp].get_xticks()   
		ax2Ticks = np.array([1/11.91, 1/24.67 ])

                   # print(ax1Ticks)
		ax2.set_xticks(ax2Ticks)
		ax2.set_xbound(ax[i].get_xbound())
		ax2.set_xticklabels(tick_function(ax2Ticks))


	ax[0].set_ylim(0.05, 20.1)
	ax[1].set_ylim(0.05, 20.1)

	ax[2].set_ylim(0.05, 10.1)
	ax[3].set_ylim(0.05, 10.1)
	ax[4].set_ylim(0.05, 10.1)
	ax[5].set_ylim(0.05, 10.1)
	ax[6].set_ylim(0.05, 10.1)
	ax[7].set_ylim(0.05, 10.1)
	ax[8].set_ylim(0.05, 10.1) 
	ax[9].set_ylim(0.05, 15.1) 

    x1,x2,y1,y2 = ax[i].axis() 
    data_set_name = [r'FEROS + HARPS RVs',r'FEROS + HARPS o-c RVs',r'HARPS BIS', r'HARPS Contrast',r'HARPS FWHM', r'HARPS CRX',r'HARPS dLW',r'HARPS H$_\alpha$',r'HARPS NaD1',r'FEROS BIS',]
    ax[i].annotate(r'%s'%data_set_name[i],(1/313.90,y2 - (y2*0.2)), fontsize=18) 
    #ax2.locator_params(axis="x", nbins=9)
    plt.setp( ax[i].get_yticklabels(), fontsize=22,weight='bold')
    plt.setp( ax[i].get_xticklabels(), fontsize=22,weight='bold')
    plt.setp( ax2.get_xticklabels(), fontsize=22,weight='bold')
 


plt.savefig('gls_act.pdf', format=format_im,dpi=dpi, bbox_inches='tight' )















