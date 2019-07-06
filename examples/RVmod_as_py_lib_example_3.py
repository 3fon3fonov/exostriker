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
sys.path.append('../lib/') #RV_mod directory must be in your path
import RV_mod as rv
import gls as gls
import numpy as np
import glob
import dill
  

targs = sorted(glob.glob("../data_Olga/FEROS/FEROS-RVSPY/*.vels"))



f  = open("test.tex", 'w') # open the file


f.write("""
\documentclass[letterpaper,10pt,oneside]{article}
\usepackage[utf8]{inputenc}
 
\usepackage{graphicx}
\usepackage[margin=0.7in]{geometry}
 
\\begin{document}
  
""")
   

        



 
for j in range(0,5):

   # print(targs[j][28:-5])
    
    rv_JD = np.genfromtxt("./%s"%(targs[j]),skip_header=0, unpack=True,skip_footer=0, usecols = [0])
 
    if len(np.atleast_1d(rv_JD))<6:
        continue

    
    # Lets create the RVmod object
    fit=rv.signal_fit('test',readinputfile=False);

    fit.cwd = '../'  
    
    fit.add_dataset(targs[j][31:],targs[j], 0.0,0.0)  # the last two entries are initial offset and jitter
   # fit.add_dataset("hip5364_crires.vels","~/git/trifon/datafiles2/hip5364_crires.vels", 0.0,0.0)
    
    # add the stellar mass
    fit.params.stellar_mass = 1.0 # In M sol.
    
    # Lets not fit for jitters now, i.e. keep at the initial value of 0 m/s
    fit.use.use_jitters[0] = False
    fit.use.use_jitters[1] = False
    
    
    #  Run it once to find the RV offsets, no planets yet.
    fit.fitting(outputfiles=[1,1,1], doGP=False,  minimize_fortran=True, minimize_loglik=False, amoeba_starts=20, print_stat=False,npoints=3000, model_max=10 )
    
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
    
   # fit.fitting(outputfiles=[1,1,1], doGP=False,  minimize_fortran=True, minimize_loglik=False, amoeba_starts=20, print_stat=False,npoints=3000, model_max=10 )

     
    rv.latex_pl_param_table(fit, width = 10, precision = 2, asymmetric = False, file_name='%s_params.tex'%targs[j][31:-5],path='./')
     
 
    #first lets copy the Keplerian object, we will need it later for plotting
   
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
    
    
    
    fig = plt.figure(0, figsize=(8,6.5))
    plt.subplots_adjust(hspace=0.005)
    format_im = 'pdf'
    dpi = 300
    
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
    #gs.update(  wspace=0.05)
    
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    
     
    color = ['b', 'r', 'g', 'r']
    symbol = ['o', 'o', 'o', 'o']   
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
    
    # we can add the jitter
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
     
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False) 
    
    
    plt.savefig('RV_plot_%s.%s'%(targs[j][31:-5],format_im), format=format_im,dpi=dpi, bbox_inches='tight' )
    ax1.cla() 
    ax2.cla() 
 

#################### GLS #####################
    
    fig = plt.figure(1, figsize=(8,6.5))
    plt.subplots_adjust(hspace=0.005)
    format_im = 'pdf'
    dpi = 300
    
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
    #gs.update(  wspace=0.05)
    
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    
 
 
    ###################################################################
    
  
    
    ax1.plot(1/fit.gls.freq, fit.gls.power, color='r',linestyle= '-')
    ax2.plot(1/fit.gls_o_c.freq, fit.gls_o_c.power, color='r',linestyle= '-')    
    
     
    power_levels = np.array([0.1,0.01,0.001]) 
    
    
              
    if fit.gls.norm == 'ZK':
        [ax1.plot(1/fit.gls.freq, np.ones(len(fit.gls.freq))*fap, '--', linewidth=2, color='k')  for ii,fap in enumerate(fit.gls.powerLevel(np.array(power_levels)))]
        [ax2.plot(1/fit.gls_o_c.freq, np.ones(len(fit.gls_o_c.freq))*fap, '--', linewidth=2, color='k')  for ii,fap in enumerate(fit.gls.powerLevel(np.array(power_levels)))]
  
    
     
    ax1.set_ylabel(r'GLS',fontsize=16, rotation = 'vertical') 
    #ax1.set_xlim(min(jd),max(jd))
    ax1.semilogx()
    ax2.semilogx()
    
    ax2.set_xlabel(r'days',fontsize=16)
    ax2.set_ylabel(r'o$-$c  GLS',fontsize=16, rotation = 'vertical') 
    #ax2.set_xlim(min(jd),max(jd))
    
    #ax2.locator_params(axis="x", nbins=9)
    plt.setp( ax2.get_yticklabels(), fontsize=15,weight='bold')
    plt.setp( ax2.get_xticklabels(), fontsize=15,weight='bold')
     
    # Fine-tune figure; make subplots close to each other and hide x ticks for
    # all but bottom plot.
     
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False) 
    
    
    plt.savefig('GLS_plot_%s.%s'%(targs[j][31:-5],format_im), format=format_im,dpi=dpi, bbox_inches='tight' )
    ax1.cla() 
    ax2.cla() 


    f.write("""
            
\\begin{figure}[ht]
   \\begin{center}$
	\\begin{array}{ccc}    
	\includegraphics[width=0.5\\textwidth]{GLS_plot_%s.%s}        
	\includegraphics[width=0.5\\textwidth]{RV_plot_%s.%s} \\\\
 
"""%(targs[j][31:-5],format_im,targs[j][31:-5],format_im))   

    ################## Phase plotting ###############

    format_im = 'png' #'pdf'
    dpi = 300

    color = ['b', 'r', 'g', 'r']
    symbol = ['o', 'D', 'o', 'o']   
    markersize = [5, 5, 6, 6] 
    alpha = [1, 1, 1, 1]

    model_color = 'k'
    model_lw = '1.0'

    for jj in range(fit.npl):

        fig = plt.figure(1, figsize=(6,6))
     
        ax1 = plt.subplot(111)
     
        planet = jj+1
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

        plt.savefig('%s_phased_planet_%s.%s'%(targs[j][31:-5],jj+1,format_im), format=format_im,dpi=dpi, bbox_inches='tight' )
        ax1.cla() 
        
        if jj%2 == 0:
            sepp = "\\\\"
        else:
            sepp = ""
        f.write("""  
	\includegraphics[width=0.5\\textwidth]{%s_phased_planet_%s.%s} %s       
    """%(targs[j][31:-5],jj+1,format_im,sepp))   


    f.write("""
 
	\end{array} $
  \end{center}
\end{figure}
    
\input{%s_params.tex} 

\clearpage
"""%(targs[j][31:-5]))   


f.write("""
            
\end{document}
""")

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



























