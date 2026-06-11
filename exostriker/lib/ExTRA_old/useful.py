import numpy as np

#for cases like Hipparchos its J1991.25, or 2448349.0625JD
#for gaia it is J2016, or 2456389.0
def J1991():
    return 2448349.0625

def J2016():
    return 2457389.0

def J2017(): #this is J2017.5 as in gaia dr4
    return 2457936.875


def jitter_estimate(residuals,err):
    s=np.mean(abs((residuals**2 -err**2 ))**0.5)
    return s
