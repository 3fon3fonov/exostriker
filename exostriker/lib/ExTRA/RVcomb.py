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

    #i introduce a -1* in the next line so the coordinate systems aline. the Z-Axis is pointing at the observer
    #usually this is not the case and its pointing away from the observer, but for astrometry its pointing towards
    
    v_r=-1*(K*(np.cos(f+om)+e*np.cos(om))+v0)
    #v_r=1*(K*(np.cos(f+om)+e*np.cos(om))+v0)

    return v_r