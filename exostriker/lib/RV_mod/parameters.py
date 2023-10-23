 #!/usr/bin/python


__author__ = 'Trifon Trifonov, Jakub Morawski'

import sys 
#sys.path.insert(0, '../lib')
import numpy as np
#import GP_kernels

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
    

