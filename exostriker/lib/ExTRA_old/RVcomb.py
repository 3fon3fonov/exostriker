from .RVsolo import *

#the RV model, just with i and a as new parameters instead of K, so its consistent with the orbital model
def RV_comb(v0,P,e,om,i,T0,a,parallax,t):
    
    M=2*np.pi*((t-T0)%P)/P
    E=calc_E(e,M)
    
    f=calc_f(E,e)
    
    #the last factor is a conversion from [mas]/[d] to [m]/[s], since a will be in [mas] and P in [d]
    
    K1=2*np.pi*a*np.sin(i)*1.495978707e11/(parallax)
    K2=(P*86400*(1-e**2)**0.5)
    K=K1/K2

    #i introduce a -1* in the next line so the coordinate systems align. the Z-Axis is pointing at the observer
    #usually this is not the case and its pointing away from the observer, but for astrometry its pointing towards
    
    #v_r=-1*(K*(np.cos(f+om)+e*np.cos(om))+v0)
    v_r=1*(K*(np.cos(f+om)+e*np.cos(om))+v0)

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
