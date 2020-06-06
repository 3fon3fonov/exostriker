#!/usr/bin/python

from __future__ import print_function
__author__ = 'Trifon Trifonov, Jakub Morawski'



        
#class for an object storing all the information about a given planet. Warning: contains less information than original kernel from old versions of RV_mod! Full information is now only stored in signal_fit class object
class kernel(object):
    
    def __init__(self,stat=0, jd=0,rvs=0,rv_err=0,o_c=0, model=0,model_jd=0,npl=0,a=[0,0,0,0,0,0,0,0,0],pl_mass=[0,0,0,0,0,0,0,0,0], idset=0,stat_array_saved=0,
                 reduced_chi2=0,chi2=0,rms=0,wrms=0,loglik=0,mfit=0, omega_dot=0,omega_dot_err=0,quad_trend=0,quad_trend_error=0,Ndata=0):
        self.stat = stat
        self.rv_model=rvmodel(jd,rvs,rv_err,o_c)
        self.jd = jd
        self.rvs = rvs
        self.rv_err = rv_err
        self.o_c = o_c
        self.model = model
        self.model_jd = model_jd
        self.npl = npl
        self.a = a
        self.mass = pl_mass
        self.idset = idset
        self.stat_array_saved=stat_array_saved
        self.reduced_chi2=reduced_chi2
        self.chi2=chi2
        self.rms=rms
        self.wrms=wrms
        self.loglik=loglik
        self.mfit = mfit
        self.omega_dot = omega_dot
        self.omega_dot_err =omega_dot_err
        self.rv_quadtr = quad_trend
        self.rv_quadtr_err= quad_trend_error 
        self.Ndata = len(self.jd)
        
class rvmodel(object):
 
    def __init__(self,jd,rvs,rv_err,o_c):
        self.jd=jd
        self.rvs=rvs
        self.rv_err=rv_err
        self.o_c=o_c
        
        # class for the stat array in kernel
class summary(object):
 
    def __init__(self,params,param_errors,dof=0):
        self.params=params
        self.param_errors=param_errors
        self.dof=dof
