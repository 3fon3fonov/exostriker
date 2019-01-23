from pylab import *
import numpy as np
import celerite 
from celerite import terms

'''Here you can define custom rotational kernels for Gaussian Processes. Give each an id, then you will be able to choose it by passing the id to the GP_parameters object
   If you don't provide an id then id=0, i. e. DefaultRotationTerm will be chosen.
   When defining a new custom kernel make sure to define a parameter_labels function (see examples here), otherwise there will be trouble with printing
'''

def choose_kernel(i, params):
	if (i==0):
	    #Kernel=DefaultRotationTerm(log_amp=np.log(params[0]),log_timescale=np.log(params[1]),log_period=np.log(params[2]),log_factor=np.log(params[3]))
 	    Kernel=DefaultRotationTerm(log_amp=params[0],log_timescale=params[1],log_period=params[2],log_factor=params[3])
       
	if not (Kernel.verify_number_of_parameters(len(params))):
	    raise Exception('Number of parameters not suitable for this kernel, this would lead to errors later on!!!')
	return Kernel
	    

class RotationTerm(terms.Term):

    ''' to make sure a desired kernel has as many parameters as provided, otherwise there will be problems '''

    def verify_number_of_parameters(self,npar):
        if (self.number_of_parameters==npar):
            return True
        else:
            return False

class DefaultRotationTerm(RotationTerm):
    parameter_names = ("log_amp", "log_timescale", "log_period", "log_factor")
    number_of_parameters=4 
    
    def parameter_labels(self):
        return ['A','t','P','f']

    def get_real_coefficients(self, params):
        log_amp, log_timescale, log_period, log_factor = params
        f = np.exp(log_factor)
        return (np.exp(log_amp) * (1.0 + f) / (2.0 + f), np.exp(-log_timescale))   


    def get_complex_coefficients(self, params):
        log_amp, log_timescale, log_period, log_factor = params
        f = np.exp(log_factor)
        return (np.exp(log_amp) / (2.0 + f), 0.0, np.exp(-log_timescale), 2*np.pi*np.exp(-log_period))
