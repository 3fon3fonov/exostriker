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
sys.path.append('../../../exostriker/lib/') #RV_mod directory must be in your path
import RV_mod as rv
import dill

from rvmod import *

file_pi = open('session.ses', 'rb')
fit = dill.load(file_pi)
file_pi.close()

fit.mod_dynamical = True

for i in range(200): 

   # re = Rvfit()
    fit.fitting(outputfiles=[1,1,1], doGP=False,  minimize_fortran=True, minimize_loglik=False, amoeba_starts=20, print_stat=False,timeout_sec=6,fortran_kill=6)
    print(i,fit.loglik)
    #flag = re.run_lm("kep")
    #print(flag)
   # del re
 
 
 















