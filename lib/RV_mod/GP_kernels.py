from pylab import *
import numpy as np
import celerite 
from celerite import terms

'''Here you can define custom rotational kernels for Gaussian Processes. Give each an id, then you will be able to choose it by passing the id to the GP_parameters object
   If you don't provide an id then id=0, i. e. DefaultRotationTerm will be chosen.
   When defining a new custom kernel make sure to define a parameter_labels function (see examples here), otherwise there will be trouble with printing
'''

    
class RotationTerm(terms.Term):
    
    parameter_names = ("log_amp", "log_timescale", "log_period", "log_factor")
   # number_of_parameters=4 
    
   # def parameter_labels(self):
   #     return ['A','t','P','f']

    def get_real_coefficients(self, params):
        log_amp, log_timescale, log_period, log_factor = params
        f = np.exp(log_factor)
        return (np.exp(log_amp) * (1.0 + f) / (2.0 + f), np.exp(-log_timescale))   

    def get_complex_coefficients(self, params):
        log_amp, log_timescale, log_period, log_factor = params
        f = np.exp(log_factor)
        return (np.exp(log_amp) / (2.0 + f), 0.0, np.exp(-log_timescale), 2*np.pi*np.exp(-log_period))



class SHOTerm2(terms.SHOTerm):

    parameter_names = ("Ampl", "Plife", "Prot")
    
    def __init__(self, params):
        Ampl, Plife, Prot = params
        
        S0 = (Ampl*(Prot**2))/(2. * (np.pi**2) * Plife)
        w0 = (2.0*np.pi)/Prot
        Q = Plife(np.pi)/Prot 
        
        kernel = terms.SHOTerm(log_S0=np.log(S0), log_Q=np.log(Q), log_omega0=np.log(w0) )