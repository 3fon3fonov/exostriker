#!/usr/bin/python3


__author__ = 'Trifon Trifonov, Jakub Morawski'

import sys 
import numpy as np
from .posterior import circular_mean, posterior_summary
#import emcee
import emcee_ES as emcee


class CustomSampler(emcee.EnsembleSampler):
 
    
    def unique_rows(self):

        '''
        Given an array, sort it
        '''

        # Perform lex sort and get sorted data
        sorted_idx = np.lexsort(self.flatchain.T)
        sorted_data =  self.flatchain[sorted_idx,:]
        sorted_lnL  =  self.flatlnprobability[sorted_idx]
        # Get unique row mask
        #row_mask = np.append([True],np.any(np.diff(sorted_data,axis=0),1))

        # Get unique rows 
        self.samples = sorted_data #[row_mask] 
        self.lnL     = sorted_lnL #[row_mask]

        lnL_max_idx = np.argmax(self.lnL)
        #print(abs(lnL_min))

        # get samples at minimum Lnl
        self.maxlnL = self.samples[lnL_max_idx]
        self.lnL_max = self.lnL[lnL_max_idx]

        return


    def correct_rows(self, f, ndset, npl, hkl=False):
        """Summarize angles on a common branch without modifying samples."""
        self.means, self.median = posterior_summary(
            self.samples, f, ndset, npl, hkl)
        return

    def circ_mean_np(self, angles, azimuth=True):
        angle = circular_mean(np.asarray(angles))
        return angle if azimuth else (angle + 180.) % 360. - 180.

    def get_meadians(self, f):
        """Retain the historical API and summaries computed by correct_rows."""
        if not hasattr(self, 'median'):
            self.median = np.median(self.samples, axis=0)
        return

    def save_samples(self, f, ndset, npl, hkl=False):
        self.unique_rows()
        self.correct_rows(f, ndset, npl, hkl)
        return
