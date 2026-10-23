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

def mas_to_deg(angle_mas):
    return angle_mas/3600000.0


# Astronomical Unit in meter, IAU constant and defining length
au_in_meter = 149597870700.0

# AU expressed in mas*pc or muas*kpc
au_mas_parsec = 1000.0

# Number of seconds in Julian year
julian_year_seconds = 365.25 * 86400.0

#astronomical unit in km yr /s
astro_unit=4.740470464

#astronomical unit in mas km yr /s
astro_unit_z=9.777922181e8

#ly
lightyear=	9460730472580800/1000 #in km

# AU expressed in km*yr/s
au_km_year_per_sec = au_in_meter / (julian_year_seconds * 1000.0)



