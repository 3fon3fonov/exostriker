

import numpy as np
from numba import jit


@jit(nopython=True)
def calc_E(e,M,n=30):
    
   
    final=np.ones(len(M))
    for index, j in enumerate(M):
        if e<=0.8:
            E=j
        if e>0.8:
            E=np.pi
        values=np.ones(n)
        values[0]=E
        i=1
        while i<n:
            values[i]=values[i-1]-(values[i-1]-e*np.sin(values[i-1])-j)/(1-e*np.cos(values[i-1]))
            if abs(values[i]-values[i-1])<1e-8:
                values[-1]=values[i]
                break
            i=i+1
            if i==n:
                return print("error in calculating E")
        final[index] = values[-1]
    
    # final=np.array(final)
            
    return final


#the RV model needs the true anomaly, we calculate it with help of the newton raphson for E

def calc_f(E,e):
    f0=(np.cos(E)-e)
    f1=np.sin(E)*(1-e**2)**0.5
    f=np.arctan2(f1,f0)
    #f=(f%(2*np.pi))
    return f

#the classic RV model
def RV_solo(v0,K,P,e,om,T0,t):
    
    M=2*np.pi*((t-T0)%P)/P
    E=calc_E(e,M)
    
    f=calc_f(E,e)
    
    
    v_r=K*(np.cos(f+om)+e*np.cos(om))+v0
    return v_r