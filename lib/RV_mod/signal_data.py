 #!/usr/bin/python


__author__ = 'Trifon Trifonov, Jakub Morawski'


#import sys 
#sys.path.insert(0, '../lib')
import numpy as np
import gls as gls 
from functions import *  
from errors import Error, InputError, FittingError
from Warning_log import Warning_log
     


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
        