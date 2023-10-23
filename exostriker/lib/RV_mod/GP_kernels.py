#from pylab import *
import numpy as np
#import celerite 
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


 


class double_SHOTerm(terms.TermSum):
    r"""A mixture of two SHO terms that can be used to model stellar rotation
    This term has two modes in Fourier space: one at ``period`` and one at
    ``0.5 * period``. This can be a good descriptive model for a wide range of
    stochastic variability in stellar time series from rotation to pulsations.
    More precisely, the parameters of the two :class:`SHOTerm` terms are
    defined as
    .. math::
        Q_1 = 1/2 + Q_0 + \delta Q \\
        \omega_1 = \frac{4\,\pi\,Q_1}{P\,\sqrt{4\,Q_1^2 - 1}} \\
        S_1 = \frac{\sigma^2}{(1 + f)\,\omega_1\,Q_1}
    for the primary term, and
    .. math::
        Q_2 = 1/2 + Q_0 \\
        \omega_2 = \frac{8\,\pi\,Q_1}{P\,\sqrt{4\,Q_1^2 - 1}} \\
        S_2 = \frac{f\,\sigma^2}{(1 + f)\,\omega_2\,Q_2}
    for the secondary term.
    Args:
        sigma: The standard deviation of the process.
        period: The primary period of variability.
        Q0: The quality factor (or really the quality factor minus one half;
            this keeps the system underdamped) for the secondary oscillation.
        dQ: The difference between the quality factors of the first and the
            second modes. This parameterization (if ``dQ > 0``) ensures that
            the primary mode alway has higher quality.
        f: The fractional amplitude of the secondary mode compared to the
            primary. This should probably always be ``0 < f < 1``, but that
            is not enforced.
    """
    @staticmethod
    def get_test_parameters():
        return dict(sigma_dSHO=1.5, period_dSHO=3.45, Q0_dSHO=1.3, dQ_dSHO=1.05, f_dSHO=0.5)    
    
    def __init__(self, *, sigma_dSHO, period_dSHO, Q0_dSHO, dQ_dSHO, f_dSHO):
        self.sigma_dSHO = float(sigma_dSHO)
        self.period_dSHO = float(period_dSHO)
        self.Q0_dSHO = float(Q0_dSHO)
        self.dQ_dSHO = float(dQ_dSHO)
        self.f_dSHO = float(f_dSHO) 

        self.amp = self.sigma_dSHO ** 2 / (1 + self.f_dSHO)

        # One term with a period of period
        Q1 = 0.5 + self.Q0_dSHO + self.dQ_dSHO
        w1 = 4 * np.pi * Q1 / (self.period_dSHO * np.sqrt(4 * Q1 ** 2 - 1))
        S1 = self.amp / (w1 * Q1)

        # Another term at half the period
        Q2 = 0.5 + self.Q0_dSHO
        w2 = 8 * np.pi * Q2 / (self.period_dSHO * np.sqrt(4 * Q2 ** 2 - 1))
        S2 = self.f_dSHO * self.amp / (w2 * Q2)

        
        kernel1 = terms.SHOTerm(log_S0=np.log(S1), log_Q=np.log(Q1), log_omega0=np.log(w1))
        kernel2 = terms.SHOTerm(log_S0=np.log(S2), log_Q=np.log(Q2), log_omega0=np.log(w2))

        #kernel = terms.TermSum(kernel1+ kernel2)
 
        super().__init__(
            kernel1, kernel2
        )



class SHOTerm2(terms.SHOTerm):

    parameter_names = ("Ampl", "Plife", "Prot")
    
    def __init__(self, params):
        Ampl, Plife, Prot = params
        
        S0 = (Ampl*(Prot**2))/(2. * (np.pi**2) * Plife)
        w0 = (2.0*np.pi)/Prot
        Q = Plife(np.pi)/Prot 
        
        kernel = terms.SHOTerm(log_S0=np.log(S0), log_Q=np.log(Q), log_omega0=np.log(w0) )
        return kernel   


class Matern32(terms.Matern32Term):

    parameter_names = ("sigma", "rho", "eps")
    
    def __init__(self, params):
        sigma, rho, eps = params
  
        kernel = terms.Matern32Term(log_sigma=np.log(sigma), log_rho=np.log(rho), eps=eps)     
        return kernel    
 
        
        
