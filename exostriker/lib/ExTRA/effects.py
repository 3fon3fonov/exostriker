import numpy as np
import scipy.constants
from astropy.time import Time
from .hipparcos import hip_JD, abs_res
from .gaia import gaia_JD
from .useful import *
from .vectorastrometry import spherical_to_cartesian, normal_triad, EpochPropagation, cartesian_to_spherical

__all__ = [
    "mu_to_v", "secular_zech", "secular_lindegren", "secular_acceleration", "secular_shift",
    "ltd_approx", "ltd_accurate", "consistent_models", "secular_correction", "propagate_v"
]

def mu_to_v(parallax,mu):
    return mu*astro_unit/ parallax






## secular formula by zechmeister, it does not need an RV
def secular_zech(sss):
    conv=pc_in_km*mas_to_rad**2 *1/julian_year_seconds *1000
    v_r_dot=conv*(sss[3]**2 + sss[4]**2) *1/sss[2]

    return v_r_dot

def secular_lindegren(sss,v_r):

    delta_mu_asc=-2*v_r*sss[2]*sss[3]/astro_unit_z 
    delta_mu_dec=-2*v_r*sss[2]*sss[4]/astro_unit_z 
    delta_par=-1*sss[2]**2 *v_r/astro_unit_z



    return np.array([delta_par,delta_mu_asc,delta_mu_dec])# in mas/yr , mas/yr^2 , mas/yr^2,km/s/yr



def secular_acceleration(sss,v_r):

    v_r_dot=secular_zech(sss)

    par_dot,mu_asc_dot,mu_dec_dot=secular_lindegren(sss,v_r)


    
    return np.array([par_dot,mu_asc_dot,mu_dec_dot,v_r_dot]) # in mas/yr , mas/yr^2 , mas/yr^2,km/s/yr

def secular_shift(t,sss,v_rad,format="years"):

    


    pc_to_km = 3.085677581e13
    d0_pc=1000/sss[2] #distance in pc using parallax
    d0=d0_pc*pc_to_km #distance in km

    

    r0_cartesian=spherical_to_cartesian(d0,np.radians(sss[0]),np.radians(sss[1])) #cartesian vector to source sss in rad rad 
    v_asc,v_dec=mu_to_v(sss[2],np.array([sss[3],sss[4]])) #mu to v using parallax
    

    if format=="years":
        v0_vec=np.array([v_asc,v_dec,v_rad])*365.25*24*60**2 #km/yr
    if format=="days":
        v0_vec=np.array([v_asc,v_dec,v_rad])*24*60**2 
    

    triad=np.array(normal_triad(np.radians(sss[0]),np.radians(sss[1])))
    v0_cartesian=triad @ v0_vec #cartesian velocity

    asc_shift=[]
    dec_shift=[]
    for time in t:
        r1_cartesian=r0_cartesian+v0_cartesian * time

    
    
        r1_spherical=cartesian_to_spherical(*r1_cartesian)
        asc_shift.append(np.degrees(r1_spherical[1])-sss[0])
        dec_shift.append(np.degrees(r1_spherical[2])-sss[1])
    


    


    return np.array([asc_shift,dec_shift])










def ltd_approx(v,format="years"): #lighttimedifference

    
    if format=="years":
        delta_t=(v*365.25*24*60*60)/lightyear  # in 1/YEARS 
        return delta_t
    if format=="days":

        delta_t_day=(v*24*60*60)/lightyear # in 1/DAYS 
        return delta_t_day



def ltd_accurate(t,sss,v_rad,format="years"):

    pc_to_km = 3.085677581e13
    d0_pc=1000/sss[2] #distance in pc using parallax
    d0=d0_pc*pc_to_km #distance in km

    

    r0_vec=spherical_to_cartesian(d0,sss[0],sss[1]) #cartesian vector to source sss in rad rad 
    r0=np.linalg.norm(r0_vec)

    

    v_asc,v_dec=mu_to_v(sss[2],np.array([sss[3],sss[4]])) #mu to v using parallax
    triad=np.array(normal_triad(np.radians(sss[0]),np.radians(sss[1])))

    if format=="years":
        v0_vec=np.array([v_asc,v_dec,v_rad])*365.25*24*60**2 #spherical velocities in km/year
    if format=="days":
        v0_vec=np.array([v_asc,v_dec,v_rad])*24*60**2
    
    v0_cartesian=triad @ v0_vec #cartesian velocity



   

    tau_ltd=[]
    for value in t:
        tau=1000*(np.linalg.norm(r0_vec+v0_cartesian*value)-r0)/(scipy.constants.c) #c in kms, v in kms r0 in pc
        tau_ltd.append(tau)

    tau_ltd=np.array(tau_ltd)/(julian_year_seconds)
    return tau_ltd

def consistent_models(hip_iad,hip_sss,gaia_sss,v_rad,hip_Sepoch=J1991(),gaia_Sepoch=J2017()):
    propagator=EpochPropagation() 

    hip_Sepoch_jyear=Time(hip_Sepoch,format="jd").jyear
    gaia_Sepoch_jyear=Time(gaia_Sepoch,format="jd").jyear


    #Notice i give gaia model absolute position in RADIANS
    #input of the propagator is [rad],[rad],[mas],[mas/yr],[mas/yr],[km/s],[yr],[yr]
    gaia1991=np.array(propagator.propagate_astrometry(np.radians(gaia_sss[0]),np.radians(gaia_sss[1]),
                                                      gaia_sss[2],gaia_sss[3],gaia_sss[4],
                                                      v_rad,hip_Sepoch_jyear,gaia_Sepoch_jyear))
    gaia1991_sss=gaia1991[:5] #abspos in radians
    #print(gaia1991_sss)
    gaia1991_sss_deg=gaia1991_sss.copy()
    #abs pos in degrees
    gaia1991_sss_deg[:2]=np.degrees(gaia1991_sss[:2])

    #input of my abs res func is [deg],[deg],[mas],[mas/yr],[mas/yr]
    abs_consistent=abs_res(hip_iad[-2],gaia1991_sss_deg,hip_sss,hip_iad)


    mu_1991=np.array([gaia1991[3],gaia1991[4],gaia1991[5]])
    v_1991=mu_to_v(gaia1991[2],mu_1991)
    v_rad_1991=v_1991[-1]

    return abs_consistent,gaia1991_sss_deg,v_rad_1991



def secular_correction(iad,sss,v_rad,Sepoch=J1991(),ltd=True):
    """ applies the secular and light time travel correction onto a dataset WITHIN itself. 
    the radial velocity has to be at the standard epoch!"""


    Sepoch_jyear = Time(Sepoch,format="jd").jyear

    if Sepoch==J1991():
        t_0=hip_JD(iad)
    else:
        t_0=gaia_JD(iad)

    t_jyear_relative = Time(t_0,format="jd").jyear -Sepoch_jyear #timestamps relative to sepoch


    if ltd==True:
        

        
        #ltd correction:

        tau_ltd=[]
        
        for t in t_jyear_relative:

            tau_ltd.append(ltd_accurate(t,sss,v_rad))
        tau_ltd=np.array(tau_ltd)/(julian_year_seconds)
        t_corrected=t_jyear_relative+tau_ltd
    
    else:
        t_corrected=t_jyear_relative
    
    
    shift_asc,shift_dec=secular_shift(sss,v_rad,t_corrected) ###this can be done more accuratly
    #please see that we used the already ltd corrected timestamps
    #removing the shift accordingly:
    abs_secular_corrected=iad[-2]-(shift_asc*iad[0]+shift_dec*iad[1])
    return abs_secular_corrected

def propagate_v(sss,v,epoch_sss,epoch_v):
    """Propagates a radial velocity NOT WORKING YET!"""
    return


        



    




