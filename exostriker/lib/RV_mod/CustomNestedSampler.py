#!/usr/bin/python


__author__ = 'Trifon Trifonov'

import sys 
#sys.path.insert(0, '../lib')
import numpy as np
from .posterior import circular_mean, posterior_summary
#from dynesty import NestedSampler, DynamicNestedSampler
import dynesty


class CustomNestedSampler(dynesty.NestedSampler,dynesty.DynamicNestedSampler):

    
    def __init__(self, *args, **kwargs):
        super(CustomNestedSampler, self).__init__(*args, **kwargs)   
    
    def convert_to_samples(self):
        self.samples = self.results.samples
    
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
        self.results.samples = out
        
        
        # same for LnL
        lnL = np.hstack(self.lnprobability)       
        sorted_lnL =  lnL[sorted_idx]   
        self.lnL = sorted_lnL[row_mask]        
        
        
        lnL_min_idx = np.argmax(self.lnL)
        #print(abs(lnL_min))
        
        # get samples at minimum Lnl
        self.minlnL = out[lnL_min_idx]
        self.lnL_min = self.lnL[lnL_min_idx]

        return


    def correct_rows(self, f, ndset, npl, hkl=False):
        """Summarize angles on a common branch without modifying samples."""
        self.means, self.median = posterior_summary(
            self.results.samples, f, ndset, npl, hkl)
        return

    def circ_mean_np(self, angles, azimuth=True):
        angle = circular_mean(np.asarray(angles))
        return angle if azimuth else (angle + 180.) % 360. - 180.

    def get_meadians(self, f):
        """Retain the historical API and summaries computed by correct_rows."""
        if not hasattr(self, 'median'):
            self.median = np.median(self.results.samples, axis=0)
        return

    def save_samples(self, f, ndset, npl, hkl=False):
        self.correct_rows(f, ndset, npl, hkl)
        return
