from .RVsolo import *
import numpy as np
__all__ = ["RV_comb", "K", "TfromM", "alpha_max","prepare_rv"]

#the RV model, just with i and a as new parameters instead of K, so its consistent with the orbital model
def RV_comb(v0,planet,t):

    P,e,om,i,Om,T0,a=planet

    #print(e)
    
    M=2*np.pi*((t-T0)%P)/P
    E=calc_E(e,M)
    
    f=calc_f(E,e)
    
    #the last factor is a conversion from [mas]/[d] to [m]/[s], since a will be in [mas] and P in [d]
    
    K1=2*np.pi*a*np.sin(i)*1.495978707e11
    K2=(P*86400*(1-e**2)**0.5)
    K=K1/K2

    #i introduce a -1* in the next line so the coordinate systems align. the Z-Axis is pointing at the observer
    #usually this is not the case and its pointing away from the observer, but for astrometry its pointing towards
    #print("K",K)
    v_r=1*(K*(np.cos(f+om)+e*np.cos(om)))+v0
    #v_r=1*(K*(np.cos(f+om)+e*np.cos(om))+v0)

    return v_r


def K(P,e,i,a,parallax):
    K1=2*np.pi*a*np.sin(i)*1.495978707e11/(parallax)
    K2=(P*86400*(1-e**2)**0.5)
    K=K1/K2
    return K


def TfromM(P,M): #M(t)=n⋅(t−Tp​) , n=  2pi /P
    T=M*P/(2*np.pi)
    return T+2500000 #in JD

def alpha_max(K1, P, e, i): #computing a in mas from the RV solution 
    """
    Computes the maximum astrometric signal (alpha_max) using radial velocity amplitude.
    
    Parameters:
        K1 (float): Radial velocity semi-amplitude (in m/s).
        P (float): Orbital period (in days).
        e (float): Orbital eccentricity.
        i (float): Orbital inclination (in degrees).
    
    Returns:
        float: Maximum astrometric signal (alpha_max) in m
    """
    

    
    # Compute alpha_max
    alpha_max_value = (K1 * P  * 86400 * np.sqrt(1 - e**2)) / (2 * np.pi * np.sin(i)*1.496e11)
    
    return alpha_max_value

def prepare_rv(RVs):
    """Run once before fitting: concatenate all instruments, sort by time."""
    t   = np.concatenate([d[0] for d in RVs])
    rv  = np.concatenate([d[1] for d in RVs])
    err = np.concatenate([d[2] for d in RVs])
    idx = np.concatenate([np.full(len(d[0]), i) for i, d in enumerate(RVs)])
    order = np.argsort(t, kind="stable")
    return dict(t=t[order], rv=rv[order], err=err[order], idx=idx[order],
                n_inst=len(RVs))

