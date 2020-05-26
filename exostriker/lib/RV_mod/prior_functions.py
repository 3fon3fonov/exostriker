#!/usr/bin/python

from __future__ import print_function
__author__ = 'Trifon Trifonov, Jakub Morawski'

import sys, os
#sys.path.insert(0, '../lib')


#from pylab import *
import numpy as np

'''Here you can define custom priors. Give each an id, then you will be able to choose it by passing the id to mcmc (prior=id), since mcmc calls the choose_prior function.
   If you don't provide prior id for mcmc id=0, i. e. flatprior, will be chosen.
'''

def choose_prior(p,i):
	if (i==0):
	    return flatprior(p)
	else:
	    return 0.0

def flatprior(p):
    '''id=0'''
    return 0.0


def normalprior(p):    
    loglik = np.log( 1./(np.sqrt(2*np.pi) * b[2]) * np.exp( - (p - b[1])**2 / (2.*b[2]**2) ) )
    return loglik

#prob_ecc1 = np.log(np.exp(-0.5 * (p[j]/SIG_E)**2.0))