import numpy as np
from astropy.table import Table

__all__ = [
    "J1991", "J2016", "J2017", "jitter_estimate", "mas_to_deg", "calc_a",
    "signature", "astropy_to_numpy", "mas_to_rad", "pc_in_km", "au_in_meter",
    "au_mas_parsec", "julian_year_seconds", "astro_unit", "astro_unit_z",
    "lightyear", "au_km_year_per_sec", "c_0","angular_distance"
]
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


def calc_a(M_s,P):
    a=(M_s*P**2)**(1/3) #P in years M in solar masses

    return a


def signature(M_s,M_p,a_p,par):
    """M_s"""
    
    s=M_p/M_s *a_p *(par/1000)

    return s



def astropy_to_numpy(table):

    arr = np.vstack([table[col] for col in table.colnames])
    return arr


def angular_distance(alpha1, delta1, alpha2, delta2):
	dalpha = alpha2 - alpha1
	numerator = np.sqrt(
		(np.cos(delta2) * np.sin(dalpha))**2
		+ (np.cos(delta1) * np.sin(delta2) - np.sin(delta1) * np.cos(delta2) * np.cos(dalpha))**2
	)
	denominator = (
		np.sin(delta1) * np.sin(delta2)
		+ np.cos(delta1) * np.cos(delta2) * np.cos(dalpha)
	)
	return np.arctan2(numerator, denominator)





mas_to_rad=1/3600 * 2*np.pi/360 *1/1000
pc_in_km=30856775814913.673
# Astronomical Unit in meter, IAU constant and defining length
au_in_meter = 149597870700.0

# AU expressed in mas*pc or muas*kpc
au_mas_parsec = 1000.0

# Number of seconds in Julian year
julian_year_seconds = 365.25 * 86400.0

#astronomical unit in km yr /s
astro_unit=4.740470464

#astronomical unit in mas km yr /s
astro_unit_z=9.7779222168e8


#ly
lightyear=	9460730472580800/1000 #in km

# AU expressed in km*yr/s
au_km_year_per_sec = au_in_meter / (julian_year_seconds * 1000.0)

c_0=299792.458 #kms



