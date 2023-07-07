 #!/usr/bin/python


__author__ = 'Trifon Trifonov, Jakub Morawski'

import sys 
#sys.path.insert(0, '../lib')
import numpy as np
from .Warning_log import Warning_log
#import GP_kernels


class GP_parameters(object): # class for GP process parameters which allows for different kernels with different number of parameters

    def __init__(self,npar,parameters,kernel_id=0):
        gpparameterswarnings=Warning_log([],'generating GP_parameters object')
        self.gp_par=parameters
        if not (npar==len(parameters)):
            npar=len(parameters)
            gpparameterswarnings.update_warning_list('Different number of parameters declared than the number of parameters actually provided! Redefined.')	          
        self.npar=npar
        self.kernel_id=kernel_id
        #self.rot_kernel=rot_kernels.choose_kernel(kernel_id,parameters)
        gpparameterswarnings.print_warning_log()

class parameters(object): # class for all parameters which can be fitted

    def __init__(self,offsets,jitters,planet_params,linear_trend,stellar_mass, GP_params=[10.0]*4, GP_kernel_id=0):
        self.offsets=np.concatenate((np.atleast_1d(offsets),[0.0]*(max(0,20-len(np.atleast_1d(offsets)))))) 
        self.jitters=np.concatenate((np.atleast_1d(jitters),[0.0]*(max(0,20-len(np.atleast_1d(jitters))))))
        self.planet_params=np.concatenate((np.atleast_1d(planet_params),[0.0]*(max(0,70-len(np.atleast_1d(planet_params))))))
        self.linear_trend=linear_trend
        self.GP_params=GP_parameters(len(GP_params),GP_params,kernel_id=GP_kernel_id) # we always want to have this attribute, but we only use it if we call GP, and then we update it anyway
        self.stellar_mass=stellar_mass
        
    def update_offset(self,dataset,offset): # change offset of one dataset
        self.offsets[dataset]=offset
        return
        
    def update_offsets(self,offsets): # change all offsets
        self.offsets=offsets
        return
        
    def update_jitter(self,dataset,jitter): # change jitter of one dataset
        self.jitters[dataset]=jitter
        return
        
    def update_jitters(self,jitters): # change all jitters 
        self.jitters=jitters
        return
        
    def update_K(self,planet,newK): # update K for a given planet
        self.planet_params[7*planet]=newK
        return
        
    def update_P(self,planet,newP): # update P for a given planet
        self.planet_params[7*planet+1]=newP
        return
        
    def update_e(self,planet,newe): # update e for a given planet
        self.planet_params[7*planet+2]=newe
        return
        
    def update_w(self,planet,neww): # update w for a given planet
        self.planet_params[7*planet+3]=neww
        return
        
    def update_M0(self,planet,newM0): # update M0 for a given planet
        self.planet_params[7*planet+4]=newM0
        return
        
    def update_inclination(self,planet,newi): # update inclination info for one planet
        self.planet_params[7*planet+5]=newi
        return              
        
    def update_lineofnodes(self,planet,newcap): # update lineofnodes info for one planet
        self.planet_params[7*planet+6]=newcap
        return              

    def update_planet_params_one_planet(self,planet,K,P,e,w,M0,i,cap): # update all planet_params for one planet
        self.update_K(planet,K)        
        self.update_P(planet,P)
        self.update_e(planet,e)
        self.update_w(planet,w)
        self.update_M0(planet,M0)
        self.update_inclination(planet,i)
        self.update_lineofnodes(planet,cap)
        return                
                
    def update_planet_params(self,planet_params): # update all planet_params in one go
        self.planet_params=planet_params
        return
                
    def update_linear_trend(self,linear_trend): # update linear trend 
        self.linear_trend=linear_trend       
        return      
        
    def update_GP_param_value(self,i,newpar):
        self.GP_params.gp_par[i]=newpar
        return
        
    def update_GP_params(self,newparams, kernel_id=-1):
        #redefine entire GP_params object, if kernel_id=-1 then we do not wish to change kernel type
        if (kernel_id==-1):
            kernel_id=self.GP_params.kernel_id
        #self.GP_params=GP_parameters(newparams,newparams,kernel_id=kernel_id)
        self.GP_params=GP_parameters(len(newparams),newparams,kernel_id=kernel_id)
        return
        
    def update_stellar_mass(self,stellar_mass):
        self.stellar_mass=stellar_mass
        return
        
    def getinit_input(self, masses, semimajor, fileinput=False, filename='geninit_j_input'):
          
        '''Prepares input for a fortran code which calculates Jacobi coordinates based on orbital parameters'''
        
        if not (fileinput): # if we want to save input in a file we don't want this line in the input string    
            ppp = './geninit_j3_in_days << EOF\n'
        else:
            ppp = ''  
        
        npl=len(self.planet_params)/7         
        
        ppp+='1 \n%s \n%s \n1.d0 \npl.in\n'%(str(self.stellar_mass),str(npl))    
        
        for i in range(npl):
            ppp+='%s \n%s %s %s %s %s %s\n'%(str(masses[i]),str(semimajor[i]),str(self.planet_params[7*i+2]),str(self.planet_params[7*i+5]),str(self.planet_params[7*i+3]),str(self.planet_params[7*i+6]),str(self.planet_params[7*i+4]))    
        
        if not (fileinput):
            ppp+='EOF' # end of the command to run in the case of saving input directly
        else: # here's what we do if we want to generate a file as well 
            # first we save the ppp string in a file (by default 'Kep_input')
            file_getin = open(filename, 'wb')
            file_getin.write('%s'%ppp) 
            file_getin.close()
            # then we overwrite ppp with the command to pass this file as input for the fortran code
            ppp='./geninit_j3_in_days < %s'%(program,filename)
                      
        return ppp 
    

class parameter_errors(object): # Class for parameter errors. 

    def __init__(self,offset_errors,jitter_errors,planet_params_errors,linear_trend_error,stellar_mass_error,GP_params_errors=[[0.0,0.0]]*4):
        '''At initiation we provide single values for each error, and we extrapolate them into 2 element arrays corresponding to + and - errors, which are at this point equal. They can be updated with updating functions found below'''
        # allocating room for up to 20 datasets and 10 planets
 
        
        offset_errors=np.concatenate((np.atleast_1d(offset_errors),[0.0]*(max(0,20-len(np.atleast_1d(offset_errors)))))) 
        jitter_errors=np.concatenate((np.atleast_1d(jitter_errors),[0.0]*(max(0,20-len(np.atleast_1d(jitter_errors))))))
        planet_params_errors=np.concatenate((np.atleast_1d(planet_params_errors),[0.0]*max((0,70-len(np.atleast_1d(planet_params_errors))))))

 

        self.offset_errors=np.array([[offset_errors[i],offset_errors[i]] for i in range(len(offset_errors))])
        self.jitter_errors=np.array([[jitter_errors[i],jitter_errors[i]] for i in range(len(jitter_errors))])
        self.planet_params_errors=np.array([[planet_params_errors[i],planet_params_errors[i]] for i in range(len(planet_params_errors))])
        self.linear_trend_error=np.array([linear_trend_error,linear_trend_error])
        self.GP_params_errors=GP_params_errors
        self.stellar_mass_error=np.array([stellar_mass_error,stellar_mass_error])

    '''In all functions below 'error' should in fact be a [-error,+error] 2 element array'''

    def update_offset_error(self,dataset,offset_error): # change offset error of one dataset
        self.offset_errors[dataset]=offset_error
        return
        
    def update_offset_errors(self,offset_errors): # change all offset errors
        self.offset_errors=offset_errors
        return
        
    def update_jitter_error(self,dataset,jitter_error): # change jitter error for one dataset
        self.jitter_errors[dataset]=jitter_error
        return
        
    def update_jitter_errors(self,jitter_errors): # change all jitters 
        self.jitter_errors=jitter_errors
        return
        
    def update_Kerror(self,planet,newKerror): # update K error for a given planet
        self.planet_params_errors[7*planet]=newKerror
        return
        
    def update_Perror(self,planet,newPerror): # update P error for a given planet
        self.planet_params_errors[7*planet+1]=newPerror
        return
        
    def update_eerror(self,planet,neweerror): # update e error for a given planet
        self.planet_params_errors[7*planet+2]=neweerror
        return
        
    def update_werror(self,planet,newwerror): # update w error for a given planet
        self.planet_params_errors[7*planet+3]=newwerror
        return
        
    def update_M0error(self,planet,newM0error): # update M0 error for a given planet
        self.planet_params_errors[7*planet+4]=newM0error
        return
        
    def update_inclination_error(self,planet,newierror): # update inclination error for one planet
        self.planet_params_errors[7*planet+5]=newierror
        return              

    def update_lineofnodes_error(self,planet,newcaperror): # update lineofnodes error for one planet
        self.planet_params_errors[7*planet+6]=newcaperror
        return              
        
    def update_planet_param_errors_one_planet(self,planet,Kerror,Perror,eerror,werror,M0error,ierror,lineofnodeserror): # update all planet_params_errors for one planet
        self.update_Kerror(planet,Kerror)        
        self.update_Perror(planet,Perror)
        self.update_eerror(planet,eerror)
        self.update_werror(planet,werror)
        self.update_M0error(planet,M0error)
        self.update_inclination_error(planet,ierror)
        self.update_lineofnodes_error(planet,lineofnodeserror)
        return                
                
    def update_planet_param_errors(self,planet_params_errors): # update all planet_param_errors in one go
        self.planet_params_errors=planet_params_errors
        return
                
    def update_linear_trend_error(self,linear_trend_error): # update linear trend error
        self.linear_trend_error=linear_trend_error       
        return    

    def update_GP_param_errors(self,i,newerror): 
        self.GP_params_errors[i]=newerror
        return

    def update_GP_params_errors(self,GP_params_errors):  
        self.GP_params_errors=GP_params_errors
        return
  
    def update_stellar_mass_error(self,stellar_mass_error):
        self.stellar_mass_error=stellar_mass_error
        return  

class parameter_bounds(object): # Class for parameter bounds. We do not need updating functions here, as updating is only done together with fitting, and this is resolved there

    def __init__(self,offset_bounds,jitter_bounds,planet_params_bounds,linear_trend_bounds,GP_params_bounds,stellar_mass_bounds):
        self.offset_bounds=offset_bounds 
        self.jitter_bounds=jitter_bounds
        self.planet_params_bounds=planet_params_bounds
        self.linear_trend_bounds=linear_trend_bounds
        self.GP_params_bounds=GP_params_bounds
        self.stellar_mass_bounds=stellar_mass_bounds
    
    
	       
                      
class use_flags(object): # class for all use flags

    def __init__(self,use_offsets,use_jitters,use_planet_params,use_linear_trend,use_stellar_mass,use_GP_params=[False]*4):
        self.use_offsets=use_offsets 
        self.use_jitters=use_jitters
        self.use_planet_params=use_planet_params
        self.use_linear_trend=use_linear_trend
        self.use_GP_params=use_GP_params
        self.use_stellar_mass=False
        self.use_offsets=np.concatenate((np.atleast_1d(self.use_offsets),[0.0]*(20-len(np.atleast_1d(self.use_offsets))))) 
        self.use_jitters=np.concatenate((np.atleast_1d(self.use_jitters),[0.0]*(20-len(np.atleast_1d(self.use_jitters)))))
        self.use_planet_params=np.concatenate((np.atleast_1d(self.use_planet_params),[0.0]*(70-len(np.atleast_1d(self.use_planet_params)))))
        
    def update_use_offset(self,dataset,use_offset): # change flag for offset of one dataset
        self.use_offsets[dataset]=use_offset
        return
        
    def update_use_offsets(self,use_offsets): # change flags for all offsets
        self.use_offsets=use_offsets
        return
        
    def update_use_jitter(self,dataset,use_jitter): # change flag for jitter of one dataset
        self.use_jitters[dataset]=use_jitter
        return
        
    def update_use_jitters(self,use_jitters): # change flags for all jitters 
        self.use_jitters=use_jitters
        return
        
    def update_use_K(self,planet,newflag): # update K flag for a given planet
        self.use_planet_params[7*planet]=newflag
        return
        
    def update_use_P(self,planet,newflag): # update P flag for a given planet
        self.use_planet_params[7*planet+1]=newflag
        return
        
    def update_use_e(self,planet,newflag): # update e flag for a given planet
        self.use_planet_params[7*planet+2]=newflag
        return
        
    def update_use_w(self,planet,newflag): # update w flag for a given planet
        self.use_planet_params[7*planet+3]=newflag
        return
        
    def update_use_M0(self,planet,newflag): # update M0 flag for a given planet
        self.use_planet_params[7*planet+4]=newflag
        return
                
    def update_use_planet_params(self,use_planet_params): # update all use_planet_params flags in one go
        self.use_planet_params=use_planet_params
        return
              
    def update_use_inclination(self,planet,newflag): # update use_inclination info for one planet
        self.use_planet_params[7*planet+5]=newflag
        return             

    def update_use_lineofnodes(self,planet,newflag): # update use_lineofnodes info for one planet
        self.use_planet_params[7*planet+6]=newflag
        return                  
        
    def update_use_planet_params_one_planet(self,planet,Kflag,Pflag,eflag,wflag,M0flag,iflag,lineofnodesflag): # update all use_planet_params flags for one planet
        self.update_use_K(planet,Kflag)        
        self.update_use_P(planet,Pflag)
        self.update_use_e(planet,eflag)
        self.update_use_w(planet,wflag)
        self.update_use_M0(planet,M0flag)
        self.update_use_inclination(planet,iflag)
        self.update_use_lineofnodes(planet,lineofnodesflag)   
        return                
                
    def update_use_linear_trend(self,use_linear_trend): # update use_linear_trend flag
        self.use_linear_trend=use_linear_trend       
        return   
        
    def update_use_GP_param(self,i,newflag):
        self.use_GP_params[i]=newflag
        return
        
    def update_use_GP_params(self,use_GP_params):
        self.use_GP_params=use_GP_params
        return             

    def update_use_stellar_mass(self,use_stellar_mass):
        self.use_stellar_mass=use_stellar_mass
        return   
    

