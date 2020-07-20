#!/usr/bin/env python3
# coding: utf-8
 
# coding: utf-8
# Created on Sun Jun  14 2020
# @author: Trifon Trifonov
 
import os, sys 
import numpy as np
#import corner as corner
import corner as corner
import dill
from scipy.ndimage import gaussian_filter
 
import time
sys.path.append('../lib/') #RV_mod directory must be in your path
import RV_mod as rv
 
################# Plotting #############################
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl

from pathos.pools import ProcessPool as Pool
    
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



################# Check script arguments ####################
arguments = len(sys.argv) - 1


if arguments != 0 and sys.argv[1] == '-ses' and os.path.exists(sys.argv[2]):
    try:
        file_pi = open(sys.argv[2], 'rb')
        fit = dill.load(file_pi)
        file_pi.close()   
 
    except (ImportError, KeyError, AttributeError) as e:
        print("You have entered wrong session. %s cannot be recognaized or it does not exist"%sys.argv[2])

else:
    print("You have entered wrong session. %s cannot be recognaized or it does not exist"%sys.argv[2])
    sys.exit(0)

if '-mass' in sys.argv:
    make_mass = True
else:
    make_mass = False

if '-semimajor' in sys.argv:
    make_a = True
else:
    make_a = False
    
if '-stab' in sys.argv:
    make_stab = True
else:
    make_stab = False

if '-best' in sys.argv:
    best = True
else:
    best = False
    
if '-median' in sys.argv:
    median = True
else:
    median = False

if '-mean' in sys.argv:
    mean = True
else:
    mean = False

if '-print_output' in sys.argv:
    print_output = True
else:
    print_output = False

if '-help' in sys.argv:
    help = True
else:
    help = False
 
 
 
################# Work starts here ########################

#### load the samples, labels and lnL values
ln = np.hstack(fit.mcmc_sampler.lnL)
samples = np.array(fit.mcmc_sampler.samples)
labels = fit.e_for_mcmc


if help:
    help_text ="""
Permitted keywords are:

-mass         (calculates and includes the pl. mass)
-semimajor    (calculates and includes the pl. semimajor axis)
-median       (posterior medians showed)
-mean         (posterior medians showed)
-print_output (prints |best|median|mean| parameters and 1sigma errors)
-help         (shows this)
"""
    print(help_text)
    print(" ")
    print("parameter indecies are:")
    for i in range(len(labels)):
        print("%s  %s"%(i,labels[i]))
    sys.exit(0)



if mean:
    best_fit_par = fit.mcmc_stat["mean"] 
    median = False
elif median:
    best_fit_par = fit.mcmc_stat["median"] 
else:
    best_fit_par = fit.mcmc_stat["best"]

 



############### make "Gaussan" samples of the stellar parameters ##############
m_s   = np.random.normal(loc=fit.stellar_mass,      scale=fit.stellar_mass_err,      size=len(samples[:,0]))
r_s   = np.random.normal(loc=fit.stellar_radius,    scale=fit.stellar_radius_err,    size=len(samples[:,0]))
L_s   = np.random.normal(loc=fit.stellar_luminosity,scale=fit.stellar_luminosity_err,size=len(samples[:,0]))
vsini = np.random.normal(loc=fit.stellar_vsini,     scale=fit.stellar_vsini_err,     size=len(samples[:,0]))



######### define new samples, labels and best-fit params to be refilled again
######### with masses, semi-major axes, etc. (to be re-worked).

samp    = []
samp_labels =  []
samp_best_fit_par = []

for i in range(len(labels)):

    ss = np.hstack(samples[:,i])
    samp.append(ss)
    samp_labels.append(labels[i])
    samp_best_fit_par.append(best_fit_par[i])
    
    
letters = ['b','c','d','e'] #... For the planets
    
if make_mass:
    for i in range(fit.npl):
        let = letters[i]
        K   = np.hstack(samples[:,[ii for ii, j in enumerate(labels) if j == 'K$_%s$'%let]])
        P   = np.hstack(samples[:,[ii for ii, j in enumerate(labels) if j == 'P$_%s$'%let]])
        ecc = np.hstack(samples[:,[ii for ii, j in enumerate(labels) if j == 'e$_%s$'%let]])
#        i   = samples[:,[ii for ii, j in enumerate(labels) if j == 'i$_%s$'%let]]

        samp.append(np.array(rv.get_mass(K,P, ecc, 90.0, m_s)))
        samp_labels.append(r'm $\sin i_%s$'%let)
        
        if mean:
            samp_best_fit_par.append(rv.get_mass(np.mean(K),np.mean(P),np.mean(ecc), 90.0, np.mean(m_s)))
        elif median:
            samp_best_fit_par.append(rv.get_mass(np.median(K),np.median(P),np.median(ecc), 90.0, np.median(m_s)))
        else:
            samp_best_fit_par.append(rv.get_mass(K[np.argmax(ln)],P[np.argmax(ln)], ecc[np.argmax(ln)], 90.0, fit.stellar_mass))


if make_a:
    for i in range(fit.npl):
        let = letters[i]
        P   = np.hstack(samples[:,[ii for ii, j in enumerate(labels) if j == 'P$_%s$'%let]])

        samp.append(rv.P_to_a(P,m_s))
        samp_labels.append(r'a$_%s$'%let)


        if mean:
            samp_best_fit_par.append(rv.P_to_a(np.mean(P),np.mean(m_s)))
        elif median:
            samp_best_fit_par.append(rv.P_to_a(np.median(P),np.median(m_s)))
        else:
            samp_best_fit_par.append(rv.P_to_a(P[np.argmax(ln)],fit.stellar_mass))
            
 
        

if make_stab:

    def circ_mean_np(angles,azimuth=True):  
        """find circular mean"""  
        rads = np.radians(angles)  
        av_sin = np.mean(np.sin(rads)) 
        av_cos = np.mean(np.cos(rads))  
        ang_rad = np.arctan2(av_sin,av_cos)  
        ang_deg = np.degrees(ang_rad)  
        if azimuth:  
            ang_deg = np.mod(ang_deg,360.)  
        return  ang_deg  

    def mcmc_satb(obj,samp_i):

        fit2 = dill.copy(obj)
        newparams = fit2.generate_newparams_for_mcmc(samp_i)
        fit2.overwrite_params(newparams) 
        
        random_dir = rv.randomString(5)
        
        rv.run_stability(fit2, timemax=5.0, timestep=0.1, timeout_sec=100.0,  integrator='mvs',stab_save_dir=random_dir)
        last_stable = min(len(fit2.evol_Per[0]),len(fit2.evol_Per[1]))


        pl1_ind = 0
        pl2_ind = 1
        Per_2 = 0
        Per_1 = 1
        
        Prat = fit2.evol_Per[pl2_ind][0:last_stable] / fit2.evol_Per[pl1_ind][0:last_stable]
        dom  = (fit2.evol_p[pl2_ind][0:last_stable] - fit2.evol_p[pl1_ind][0:last_stable])%360
        dom[dom>=180.0] -= 360.0

        lambda1  = (fit2.evol_M[pl1_ind][0:last_stable]    + fit2.evol_p[pl1_ind][0:last_stable]   + 0)%360
        lambda2  = (fit2.evol_M[pl2_ind][0:last_stable]    + fit2.evol_p[pl2_ind][0:last_stable]   + 0)%360
 
 

        
        theta = {k: [ ] for k in range(10)}    
        coef1 = Per_2 +1
        coef2 = Per_1 +1
        order = abs(coef2 - coef1)
 
        for i in range(order+1):
 
            theta[i] = (coef1*lambda1%360 - coef2*lambda2%360 )%360 + (((coef2 -coef1) -i)*fit2.evol_p[pl1_ind][0:last_stable] + i*fit2.evol_p[pl2_ind][0:last_stable])%360 
            theta[i] = theta[i]%360
            theta[i][theta[i]>=180.0] -= 360.0
            
        theta[0] = theta[0] - circ_mean_np(theta[0])
        theta[1] = theta[1] - circ_mean_np(theta[1])
        dom      = dom - circ_mean_np(dom)
        
        theta_1_amp = (max(theta[0]) - min(theta[0]))/2.0
        theta_2_amp = (max(theta[1]) - min(theta[1]))/2.0
        dom_amp = (max(dom) - min(dom))/2.0
        
        
        del fit2 
        
        return [
        np.mean(np.array(Prat)),
        np.mean(np.array(dom_amp)),
        np.mean(np.array(theta_1_amp)),
        np.mean(np.array(theta_2_amp))]


    os.chdir('../')
    
    ncpus = 4
    fit2 = dill.copy(fit)
    del fit2.mcmc_sampler

    def partial_func(samp_ii):
        loglik = mcmc_satb(fit2, samp_ii)
        return loglik

    with Pool(ncpus=ncpus) as thread: 
        results = thread.map(partial_func,   [samp_i for samp_i in samples])
 

    thread.close()
    thread.join()
    thread.clear()

        
    os.chdir('./MCMC/')

    #print(len(results[0]))
    results = np.transpose(results)
   # print(len(results[0]))

    mean_Prat = np.array(results[0])
    samp.append(mean_Prat)
    samp_labels.append(r'P$_{\rm rat}$')
    
    if mean:
        samp_best_fit_par.append(np.mean(mean_Prat))
    elif median:
        samp_best_fit_par.append(np.median(mean_Prat))
    else:
        samp_best_fit_par.append(mean_Prat[np.argmax(ln)])
        
    dom_amp = np.array(results[1])
    samp.append(dom_amp)
    samp_labels.append(r'$\Delta\omega$')
    
    if mean:
        samp_best_fit_par.append(np.mean(dom_amp))
    elif median:
        samp_best_fit_par.append(np.median(dom_amp))
    else:
        samp_best_fit_par.append(dom_amp[np.argmax(ln)])
    
    t1_amp = np.array(results[2])
    samp.append(t1_amp)
    samp_labels.append(r'$\Delta\theta 1$')
    
    if mean:
        samp_best_fit_par.append(np.mean(t1_amp))
    elif median:
        samp_best_fit_par.append(np.median(t1_amp))
    else:
        samp_best_fit_par.append(t1_amp[np.argmax(ln)])
    
    t2_amp = np.array(results[3])
    samp.append(t2_amp)
    samp_labels.append(r'$\Delta\theta 2$')

    if mean:
        samp_best_fit_par.append(np.mean(t2_amp))
    elif median:
        samp_best_fit_par.append(np.median(t2_amp))
    else:
        samp_best_fit_par.append(t2_amp[np.argmax(ln)])



 
############# Here you can make changes to the input structures ################
 
# e.g. to change the labels:
 
#samp_labels[2] = "RV jitt Lick m/s"
#samp_labels[3] = "RV jitt CRIRES m/s"
 


########### if you want to remove parameters from the cornerplot use this: ####

#del samp[2] 
#del samp_labels[2] 
#del samp_best_fit_par[2] 

#del samp[3]               # repeat for another element. (use -help to get the indices)
#del samp_labels[3] 
#del samp_best_fit_par[3] 


######### here index 2 will remove the Lick jitter from the plot. #############
######### (there must be a better way to do this but it works for now ) #######



################### Transpose is needed fro the cornerplot. ###################
samples_ = np.transpose(samp)
labels = samp_labels
best_fit_par =samp_best_fit_par
 


if print_output:


    print("Best fit par. and their 1 sigma errors")
    lines = ""
    for i in range(len(best_fit_par)):
        ci = np.percentile(samples_[:,i], [level, 100.0-level])
        print(labels[i],'=', best_fit_par[i], "- %s"%(best_fit_par[i]-ci[0]), "+ %s"%(ci[1]  - best_fit_par[i] ))
        lines += ("%s   &  %.3f$_{-%.3f}^{+%.3f}$ \\\\ \\noalign{\\vskip 0.9mm} \n"%(labels[i],best_fit_par[i], best_fit_par[i]-ci[0], ci[1]  - best_fit_par[i] ))

    print("")
    print("---------------")
    print(lines)


############# Making the cornerplot #################


fig = corner.corner(samples_,
bins=25, 
color="k", 
reverse=False, 
labels=labels, 
quantiles=[0.1585, 0.8415],
levels=(0.6827, 0.9545,0.9973),
smooth=1.0, 
smooth1d=1.0,
plot_contours= True, 
show_titles=True, 
truths=best_fit_par, 
dpi = 200, 
pad=15, 
labelpad = 50 ,
truth_color ='r', 
title_kwargs={"fontsize": 12}, 
scale_hist=True,  
no_fill_contours=True, 
plot_datapoints=True)


fig.savefig("MCMC_cornerplot.pdf")





 
############### Junk, not used but that might be useful ############


junk = """

for j in range(len(best_fit_par)):
    par[f[j]] = mean_params[j] 


file_ses = open(r"best_fit_mean.pkl", 'wb')
dill.dump(np.array(par), file_ses)
file_ses.close()

for j in range(len(best_fit_par)):
    par[f[j]] = median_params[j] 

file_ses = open(r"best_fit_median.pkl", 'wb')
dill.dump(np.array(par), file_ses)
file_ses.close()


for j in range(len(best_fit_par)):
    par[f[j]] = mode_params[j] 

file_ses = open(r"best_fit_mod.pkl", 'wb')
dill.dump(np.array(par), file_ses)
file_ses.close()

for j in range(len(best_fit_par)):
    par[f[j]] = best_fit_par[j] 

file_ses = open(r"best_fit_NS.pkl", 'wb')
dill.dump(np.array(par), file_ses)
file_ses.close()   

 
#print(len(samples))
#samples = samples[ln  > max(lnZ) - 15.0]
# print(len(samples))

# samples = np.array(samp) #np.transpose(samp)
# samples = np.transpose(samples)

print("Best fit par. and their 1 sigma errors")
for i in range(len(best_fit_par)):
    ci = np.percentile(samples[:,i], [level, 100.0-level])
    print(e[i],'=', best_fit_par[i], "- %s"%(best_fit_par[i]-ci[0]), "+ %s"%(ci[1]  - best_fit_par[i] ))
    #lines += ("%s   &  %.3f$_{-%.3f}^{+%.3f}$ \\\\ \\noalign{\\vskip 0.9mm} \n"%(e[i],best_fit_par[i], best_fit_par[i]-ci[0], ci[1]  - best_fit_par[i] ))

for i in range(len(der_param)):
    ci = np.percentile(der_param[i], [level, 100.0-level])
    print(e_der_param[i],'=', der_param[i][np.argmax(ln)], "- %s"%(der_param[i][np.argmax(ln)]-ci[0]), "+ %s"%(ci[1]  - der_param[i][np.argmax(ln)] ))
    #mean_params.append(der_param[i][np.argmax(ln)])

print("   ")
print("   ")
print("Means and their 1 sigma errors")
for i in range(len(best_fit_par)):
    ci = np.percentile(samples[:,i], [level, 100.0-level])
    print(e[i],'=', np.mean(samples[:,i]), "- %s"%(np.mean(samples[:,i])-ci[0]), "+ %s"%(ci[1]  - np.mean(samples[:,i]) ))
    mean_params.append(np.mean(samples[:,i]))
for i in range(len(der_param)):
    ci = np.percentile(der_param[i], [level, 100.0-level])
    print(e_der_param[i],'=', np.mean(der_param[i]), "- %s"%(np.mean(der_param[i])-ci[0]), "+ %s"%(ci[1]  - np.mean(der_param[i]) ))
    #mean_params.append(np.mean(der_param[i]))


print("   ")
print("   ")
print("Median and their 1 sigma errors")


############# Print the LateX lines (not a full Latex Table) #################


for i in range(len(best_fit_par)):


    if pr_nr[i][2]==True:
        sign,f_arg,s_arg,pow_arg = "$\mathcal{N}$",pr_nr[i][0],pr_nr[i][1],"$^2$"
    elif jeff_nr[i][2]==True:
        sign,f_arg,s_arg,pow_arg = "$\mathcal{J}$",jeff_nr[i][0],jeff_nr[i][1],""
    else:
        sign,f_arg,s_arg,pow_arg = "$\mathcal{U}$",bb[i][0],bb[i][1],""
        
    #    text = text + '''& {0:s}({1:{width}.{precision}f},{2:{width}.{precision}f}{3:s})'''.format(sign, f_arg,s_arg,pow_arg, width = width, precision = precision)
    #text = text + '''\\\\

    ci = np.percentile(samples[:,i], [level, 100.0-level])
    print(e[i],'=', np.median(samples[:,i]), "- %s"%(np.median(samples[:,i])-ci[0]), "+ %s"%(ci[1]  - np.median(samples[:,i]) ))
    median_params.append(np.median(samples[:,i]))
    if i == 3:
        lines += ("%s   &  %.5f$_{-%.5f}^{+%.5f}$  &  %.5f   &  %s(%.5f,%.5f%s) &   \\\\ \\noalign{\\vskip 0.9mm} \n"%(e[i],np.median(samples[:,i]),np.median(samples[:,i])-ci[0], ci[1]  -np.median(samples[:,i]),  best_fit_par[i],sign, f_arg,s_arg,pow_arg ))
    elif i == 5 or i == 7:
        lines += ("%s   &  %.4f$_{-%.4f}^{+%.4f}$  &  %.4f   &  %s(%.4f,%.4f%s) &   \\\\ \\noalign{\\vskip 0.9mm} \n"%(e[i],np.median(samples[:,i]),np.median(samples[:,i])-ci[0], ci[1]  -np.median(samples[:,i]),  best_fit_par[i],sign, f_arg,s_arg,pow_arg ))
        
    else:
        lines += ("%s   &  %.2f$_{-%.2f}^{+%.2f}$  &  %.2f   &  %s(%.2f,%.2f%s) &   \\\\ \\noalign{\\vskip 0.9mm} \n"%(e[i],np.median(samples[:,i]),np.median(samples[:,i])-ci[0], ci[1]  -np.median(samples[:,i]),  best_fit_par[i],sign, f_arg,s_arg,pow_arg ))






for i in range(len(der_param)):
    ci = np.percentile(der_param[i], [level, 100.0-level])
    print(e_der_param[i],'=', np.median(der_param[i]), "- %s"%(np.median(der_param[i])-ci[0]), "+ %s"%(ci[1]  - np.median(der_param[i]) ))
    #median_params.append(np.median(der_param[i]))
    if i == 1:
        lines += ("%s   &  %.5f$_{-%.5f}^{+%.5f}$  &  %.5f   &   $   \\\\ \\noalign{\\vskip 0.9mm} \n"%(e_der_param[i],np.median(der_param[i]),np.median(der_param[i])-ci[0], ci[1]  -np.median(der_param[i]),  der_param[i][np.argmax(ln)]))
    else:
        lines += ("%s   &  %.2f$_{-%.2f}^{+%.2f}$  &  %.2f   &   $   \\\\ \\noalign{\\vskip 0.9mm} \n"%(e_der_param[i],np.median(der_param[i]),np.median(der_param[i])-ci[0], ci[1]  -np.median(der_param[i]),  der_param[i][np.argmax(ln)]))
print("   ")
print("   ")    
print("Mode and their 1 sigma errors")    




for i in range(len(best_fit_par)):
    ci = np.percentile(samples[:,i], [level, 100.0-level])
    #mmm = stats.binned_statistic(np.array([samples[:,i]]), axis=None)
    n, b = np.histogram(samples[:,i], bins=100)
    n = gaussian_filter(n, 1.0)
    x0 = np.array(list(zip(b[:-1], b[1:]))).flatten()
    y0 = np.array(list(zip(n, n))).flatten()    
    k  = np.unravel_index(y0.argmax(),y0.shape)
    mode_params.append(x0[k])
    print(e[i],'=', x0[k], "- %s"%(x0[k]-ci[0]), "+ %s"%(ci[1]  - x0[k] ))
    #lines += ("%s   &  %.3f$_{-%.3f}^{+%.3f}$ \\\\ \\noalign{\\vskip 0.9mm} \n"%(e[i],x0[k], x0[k]-ci[0], ci[1]  - x0[k] ))




 
print("")
print("---------------")
print(lines)




"""



