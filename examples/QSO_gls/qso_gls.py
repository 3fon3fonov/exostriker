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

#import sys 
#sys.path.append('../../lib/') #RV_mod directory must be in your path
 
import dill
import sys,  os
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

 
 

files = [     
"binned_BIS_J00183+440.ccfpar.dat",
"binned_contrast_J00183+440.ccfpar.dat"]

 
 
color_gls = ['0.3']*len(files)



def make_GLS(x,y,e_y, name=''):


    f = plt.figure(1, figsize=(9,7.5))
    plt.subplots_adjust(hspace=0.005)
    format_im = 'png'
    dpi = 150
    ax = plt.subplot(1,1,1)


    omega = 1/ np.logspace(np.log10(9.99), np.log10(10000), num=1000)
    power_levels = np.array([0.1000,0.01000,0.001000])

    act_per = gls.Gls((x,y,e_y), fast=True, verbose=False, norm="ZK",ofac=5, fbeg=1/max(x)*2, fend=1/10.0)
    
    ax2 = ax.twiny() 
    ax2.yaxis.set_major_locator(plt.MaxNLocator(3)) 

    factor = (len(x) -1.0) /2.0
 
    #ax.plot([1/11.44, 1/11.44],[0.1,40], '--',color ='b',  alpha=0.7, lw=1.0, label='')
 
 
    [ax.plot([max(act_per.freq),min(act_per.freq)],[fap*factor,fap*factor], linestyle='--',color='k',zorder=10) for ii,fap in enumerate(act_per.powerLevel(np.array(power_levels)))] 
    ax.plot(act_per.freq,act_per.power*factor, linestyle='-',color=color_gls[i], lw=1.2,)                    


    ax.set_xlim( 1/10000.0,1/10.0)       
    
      
 
    ax.tick_params(direction="in")
    ax2.tick_params(direction="in") 

    ax.set_ylabel(r'GLS Power',fontsize=24, rotation = 'vertical', labelpad = 15) 
    ax.set_xlabel(r'frequency [1/d]',fontsize=24)
    ax2.set_xlabel(r'period [d]',fontsize=24) 



    def tick_function(X):
        V = 1/X
        return ["%.1f" % z for z in V]
    #ax1Ticks = ax[pp].get_xticks()   
    ax2Ticks = np.array([1/10.0, 1/25.0, 1/10000.0, 1/50.0 ])

               # print(ax1Ticks)
    ax2.set_xticks(ax2Ticks)
    ax2.set_xbound(ax.get_xbound())
    ax2.set_xticklabels(tick_function(ax2Ticks))


#    ax[0].set_ylim(0.05, 38.9)
 

    x1,x2,y1,y2 = ax.axis() 
#    data_set_name = [r'GJ~15A HARPS-N BIS']
#    ax.annotate(r'%s'%data_set_name[i],(1/18.0,y2 - (y2*0.2)), fontsize=18) 
    #ax2.locator_params(axis="x", nbins=9)
    plt.setp( ax.get_yticklabels(), fontsize=22,weight='bold')
    plt.setp( ax.get_xticklabels(), fontsize=22,weight='bold')
    plt.setp( ax2.get_xticklabels(), fontsize=22,weight='bold')
 
 
 
    plt.savefig('QSO_%s.pdf'%name, format=format_im,dpi=dpi, bbox_inches='tight' )






for i in range(len(files)):


    x   = np.genfromtxt("%s"%(files[i]),skip_header=0, unpack=True,skip_footer=0, usecols = [0])  #t_day
    y   = np.genfromtxt("%s"%(files[i]),skip_header=0, unpack=True,skip_footer=0, usecols = [1])  # y_mag
    e_y = np.genfromtxt("%s"%(files[i]),skip_header=0, unpack=True,skip_footer=0, usecols = [2])  # y_mag_err

    
    make_GLS(x,y,e_y,name=files[i])
    
    













