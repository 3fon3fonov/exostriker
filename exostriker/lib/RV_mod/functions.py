 #!/usr/bin/python3

__author__ = 'Trifon Trifonov'

import sys, os
#sys.path.insert(0, '../lib')
import numpy as np
import jac2astrocen
import corner_ES as corner
from sample_angles import fold_plot_angles
from posterior_display import print_summary
from hkl_plot import derived_hkl_columns
from sample_plot_labels import default_rv_units

import matplotlib
matplotlib.pyplot.switch_backend('Agg')

import re
import subprocess
from subprocess import PIPE, Popen
import signal
import platform
import tempfile, shutil
from threading import Thread
from .Warning_log import Warning_log
from .jacobi import mass_a_from_jacobi_el
import scipy.stats as pdf
import dill
from scipy.signal import argrelextrema
from scipy.ndimage import gaussian_filter
import random
import string
import ntpath
from pathlib import Path, PureWindowsPath


import gls as gls

TAU= 2.0*np.pi
RSUN_AU = 0.00465047
 

G = 6.67430e-11
DAY = 86400.0
AU = 1.49597892e11
MSUN = 1.9884098709677423e30 #1.98847e30 
GMSUN = 1.32712440018e20
REARTH = 6.371e6
RJUP = 7.1492e7

THIRD = 1.0/3.0
PI    = 3.14159265358979e0
TWOPI = 2.0*PI




matplotlib.rcParams['axes.formatter.useoffset'] = False

def migrate_ttv_datasets(fit):
    """Upgrade saved TTV records without reloading their observation files."""
    datasets = getattr(fit, 'ttv_data_sets', None)
    if datasets is None:
        return
    for slot in range(len(datasets)):
        record = datasets[slot]
        if not len(record):
            continue
        values = list(record)
        if len(values) == 6:
            size = len(np.atleast_1d(values[0]))
            values.extend([np.ones(size), np.full(size, .1),
                           np.zeros(size), np.full(size, .1),
                           np.nan, np.nan, False, False])
        elif len(values) == 12:
            # Optional depth/duration columns predate density/radius metadata.
            values = values[:10] + [np.nan, np.nan] + values[10:12]
        else:
            continue
        datasets[slot] = np.array(values, dtype=object)


def check_for_missing_instances(fit,fit_new):

    for iii in fit.__dict__:
        if iii not in fit_new.__dict__:
            fit_new.__dict__[iii] = dill.copy(fit.__dict__[iii])
#        elif iii in fit_new.__dict__ and len(np.atleast_1d(fit_new.__dict__[iii])) != len(np.atleast_1d(fit.__dict__[iii])):
#            fit_new.__dict__[iii] = dill.copy(fit.__dict__[iii])

    for iii in fit.fit_results.__dict__:
        if iii not in fit_new.fit_results.__dict__:
            fit_new.fit_results.__dict__[iii] = dill.copy(fit.fit_results.__dict__[iii])

    #fit.mass_semimajor = dill.copy(fit.mass_semimajor)


    try:
        if len(np.atleast_1d(fit_new.ns_sampler))!=0:
            try:
                fit_new.ns_sampler.lbf = dill.copy(fit_new.ns_sampler.lbf)
            except:
                pass
     
        if len(np.atleast_1d(fit_new.mcmc_sampler))!=0:
            try:
                fit_new.mcmc_sampler.lbf = dill.copy(fit_new.mcmc_sampler.lbf)
            except:
                pass
    except Exception as e:
        print("Warning: Could not load nested/mcmc sampler:", e)
        fit_new.ns_sampler = []
        fit_new.mcmc_sampler = []

    if len(fit_new.type_fit) != 4:
        try:
            fit_new.type_fit = dill.copy(fit.type_fit)
        except:
            pass

    fit_new.bounds  = dill.copy(fit.bounds)

    if len(fit_new.rvoff_bounds) <= 11:
        fit_new.rvoff_bounds = dill.copy(fit.rvoff_bounds)
        fit_new.jitt_bounds = dill.copy(fit.jitt_bounds)
        fit_new.rvoff_norm_pr = dill.copy(fit.rvoff_norm_pr)
        fit_new.jitt_norm_pr = dill.copy(fit.jitt_norm_pr)
        fit_new.rvoff_jeff_pr = dill.copy(fit.rvoff_jeff_pr)
        fit_new.jitt_jeff_pr = dill.copy(fit.jitt_jeff_pr)
        fit_new.rvs_colors = dill.copy(fit.rvs_colors)
        fit_new.pyqt_symbols_rvs  = dill.copy(fit.pyqt_symbols_rvs)




    if len(fit_new.tra_colors) <= 11:
        fit_new.tra_colors = dill.copy(fit.tra_colors)
        fit_new.tra_quadtr_jeff_pr = dill.copy(fit.tra_quadtr_jeff_pr)
        fit_new.tra_jitt_use  = {k: False for k in range(60)}
        fit_new.tra_off_use  = {k: False for k in range(60)}
        fit_new.tra_dil_use  = {k: False for k in range(60)}

        fit_new.init_tra_jitter()
        fit_new.init_tra_offset()
        fit_new.init_tra_dilution()
        fit_new.init_tra_lintr()
        fit_new.init_tra_quadtr()


    migrate_ttv_datasets(fit_new)
    fit_new.cwd = dill.copy(fit.cwd)

    return fit_new





def inc_from_rsky_vsky(rsky_AU, a_over_R, Rstar_Rsun):
    Rstar_AU = Rstar_Rsun * RSUN_AU
    b = rsky_AU / Rstar_AU
    x = b / a_over_R
    x = np.clip(x, -1.0, 1.0)
    inc_deg = np.degrees(np.arccos(x))
    return inc_deg, b


def transit_tperi_old(per, ecc, om, ma, epoch):
    """It derives Time of periatron [tp]
    and time of mid transit [t0]

    Parameters
    ----------
    per : float
        Period of the planet [days].
    ecc : float
        Eccentricity of the orbit.
    om : float
        Argument of periastron [deg]
    ma : float
        Mean anomaly [deg].
    epoch : float
        Epoch for wich the orbital elements are valid [BJD].

    Returns
    -------
    [tp,t0]
        if the epoch in is BJD then tp and t0 are also in BJD.
    """
    om = np.radians(om)
    ma = np.radians(ma)

    E = 2.0*np.arctan( np.sqrt( ( (1.0-ecc)/(1.0+ecc) ) ) * np.tan( (np.pi/4.0)-(om/2.0) ) )
    t_peri    = epoch  - ((ma/TAU)*per)
    t_transit = t_peri + (E + ecc*np.sin(E)) * (per/TAU)

    return t_peri, t_transit


def transit_tperi_old(per, ecc, om, ma, epoch, primary = True):
    """It derives Time of periatron [tp]
    and time of mid transit [t0]
    Parameters
    ----------
    per : float
        Period of the planet [days].
    ecc : float
        Eccentricity of the orbit.
    om : float
        Argument of periastron [deg]
    ma : float
        Mean anomaly [deg].
    epoch : float
        Epoch for wich the orbital elements are valid [BJD].
    Returns
    -------
    [tp,t0]
        if the epoch in is BJD then tp and t0 are also in BJD.
    """

    trueA = np.pi/2.0
    om = np.radians((om)%360)
    ma = np.radians((ma)%360)
    f = trueA - om
    E = 2.0*np.arctan( np.sqrt( (1.0-ecc)/(1.0+ecc) ) * np.tan(f/2.0)  )
    #t_peri    =  epoch - (per/TAU)*(E - ecc*np.sin(E))
    t_peri    = epoch  - ((ma/TAU)*per)
    t_transit = t_peri + (E - ecc*np.sin(E)) * (per/TAU)

    if primary != True:
        trueA = 3.0 * np.pi/2.0
        f = trueA - om
        E = 2.0*np.arctan( np.sqrt( (1.0-ecc)/(1.0+ecc) ) * np.tan(f/2.0)  )
        t_transit = t_peri + (E - ecc*np.sin(E)) * (per/TAU)
        return t_peri, t_transit


    return t_peri, t_transit
    
 

def transit_tperi(per, ecc, om, ma, epoch, primary=True):
    """
    Estimate time of periastron passage and time of conjunction/transit.

    Parameters
    ----------
    per : float
        Orbital period [days].
    ecc : float
        Eccentricity.
    om : float
        Argument of periastron [deg].
    ma : float
        Mean anomaly at epoch [deg].
    epoch : float
        Epoch at which the orbital elements are defined [BJD].
    primary : bool
        If True, compute inferior conjunction, approximate primary transit.
        If False, compute superior conjunction, approximate secondary eclipse.

    Returns
    -------
    t_peri, t_transit : floats
        Time of previous periastron passage and approximate transit/conjunction time.
    """

    if per <= 0.0:
        raise ValueError("Orbital period must be positive.")

    if ecc < 0.0 or ecc >= 1.0:
        raise ValueError("This function is valid only for elliptical orbits, 0 <= ecc < 1.")

    om = np.radians(om % 360.0)
    ma = np.radians(ma % 360.0)

    # Time of previous periastron passage relative to epoch.
    t_peri = epoch - (ma / TAU) * per

    # Approximate conjunction geometry for an edge-on orbit.
    if primary:
        u_conj = 0.5 * np.pi
    else:
        u_conj = 1.5 * np.pi

    # True anomaly at conjunction.
    f = u_conj - om
    f = f % TAU

    # Convert true anomaly to eccentric anomaly using quadrant-safe formula.
    E = 2.0 * np.arctan2(
        np.sqrt(1.0 - ecc) * np.sin(0.5 * f),
        np.sqrt(1.0 + ecc) * np.cos(0.5 * f)
    )
    E = E % TAU

    # Mean anomaly at conjunction.
    M_conj = E - ecc * np.sin(E)
    M_conj = M_conj % TAU

    # Time of conjunction after periastron.
    t_transit = t_peri + (M_conj / TAU) * per

    # Put the transit time close to the input epoch.
    while t_transit < epoch - 0.5 * per:
        t_transit += per

    while t_transit > epoch + 0.5 * per:
        t_transit -= per

    return t_peri, t_transit    


def get_m0(per, ecc, om, t0, epoch):
    '''
    '''
    om = np.radians(om)
    f = np.pi/2.0 - om
    E = 2.0*np.arctan( np.sqrt( (1.0-ecc)/(1.0+ecc) ) * np.tan(f/2.0)  )

#    t_peri    =  epoch  - (per/TAU)*(E - ecc*np.sin(E))
#    print(t_peri)
 #   t0 = 2458334.3166

    #m0 = E - ecc*np.sin(E) + 2.0*np.pi*( (epoch-t0)/per % 1.)
    #ma = np.degrees(ma)%360.0

    ma = E - ecc*np.sin(E)

    m0 = ma + 2.0*np.pi*( (epoch-t0)/per % 1.)
    #m0 = ma -  (t0-epoch)* 2.0*np.pi/per

    m0 = np.degrees(m0)%360.0


    return m0

def ma_from_t0(per, ecc, om, t_transit, epoch):
    '''
    '''
    om = np.radians(om)
    E = 2.0*np.arctan( np.sqrt( ( (1.0-ecc)/(1.0+ecc) ) ) * np.tan( (np.pi/4.0)-(om/2.0) ) )
   # t_transit = epoch  - ((ma/TAU)*per) + (E + ecc*np.sin(E)) * (per/TAU)

    ma =  ((epoch  - t_transit + (E + ecc*np.sin(E)) * (per/TAU))*TAU)/per
    ma = np.degrees(ma)%360.0

    return ma


def ma_from_epoch_old(per, t_peri, epoch):
    '''
    '''
    ma =  np.degrees(2.0*np.pi*( (epoch-t_peri)/per % 1.))

    return ma
    
    
def ma_from_epoch(per, t_peri, epoch):
    """
    Mean anomaly at epoch, in degrees, wrapped to [0, 360).

    Parameters
    ----------
    per : float or array
        Orbital period, same time units as epoch and t_peri.
    t_peri : float or array
        Time of periastron passage.
    epoch : float or array
        Epoch at which the mean anomaly is evaluated.

    Returns
    -------
    ma : float or array
        Mean anomaly in degrees, in the range [0, 360).
    """
    phase = np.mod((epoch - t_peri)/per, 1.0)
    ma = 360.0 * phase
    
    return ma    

def mass_to_K(P,ecc,incl, pl_mass,Stellar_mass):

    '''Returns the RV semi-amplitude K in m/s

    Parameters
    ----------
    P : float
        Period of the planet in [d]
    ecc : float
        eccentricity
    incl: float
        inclination in [deg]
    pl_mass: float
        planet mass in [Msol]
    Stellar_mass: float
        Primary mass in [Msol]

    Returns
        float
    K in [m/s]
    -------
    '''
    THIRD = 1.0/3.0
    GMSUN = 1.32712440018e20
    AU=1.49597892e11

    T = P*86400.0

    #K = ((2.0*np.pi*GMSUN)/T)**THIRD * (pl_mass*np.sin(np.radians(incl)) /
    #    (Stellar_mass+pl_mass)**(2.0/3.0)) * 1.0/np.sqrt(1.0-ecc**2.0)

  #  K = ((2.0*np.pi*GMSUN)/T)**THIRD * (pl_mass*np.sin(np.radians(incl)) /
   #     (Stellar_mass+pl_mass)**(2.0/3.0))
    K = ((2.0*np.pi*GMSUN)/T)**THIRD * (pl_mass*np.sin(np.radians(incl)) /
        (Stellar_mass+pl_mass)**(2.0/3.0)) * 1.0/np.sqrt(1.0-ecc**2.0)

    return K


def uniform_disk_depth_from_b_k_vectorized(b, k):
    """
    Uniform-disk transit depth for projected separation b and radius ratio k.

    b can be scalar or array.
    k must be scalar.

    Returns relative flux loss.
    """

    b = np.asarray(b, dtype=float)
    z = np.abs(b)
    k = float(k)

    depth = np.zeros_like(z, dtype=float)

    if k <= 0.0:
        return depth

    # No overlap
    no_overlap = z >= 1.0 + k

    # Planet fully inside stellar disk
    full_planet = z <= 1.0 - k
    depth[full_planet] = k*k

    # Star fully inside planet disk, mostly irrelevant for exoplanets
    full_star = z <= k - 1.0
    depth[full_star] = 1.0

    # Partial overlap
    partial = ~(no_overlap | full_planet | full_star)

    if np.any(partial):

        zp = z[partial]

        arg1 = (zp*zp + 1.0 - k*k) / (2.0*zp)
        arg2 = (zp*zp + k*k - 1.0) / (2.0*zp*k)

        arg1 = np.clip(arg1, -1.0, 1.0)
        arg2 = np.clip(arg2, -1.0, 1.0)

        area = (
            np.arccos(arg1)
            + k*k*np.arccos(arg2)
            - 0.5*np.sqrt(
                np.maximum(
                    (-zp + 1.0 + k) *
                    ( zp + 1.0 - k) *
                    ( zp - 1.0 + k) *
                    ( zp + 1.0 + k),
                    0.0
                )
            )
        )

        depth[partial] = area / np.pi

    return depth 


def transit_duration_depth_from_rsky_vsky(
        rsky_AU,
        vsky_AU_day,
        rho_kg_m3,
        Mstar_Msun,
        Rp_Rstar):
    """
    Estimate transit durations and uniform-disk transit depths from TTVFast
    rsky and vsky at mid-transit.

    This version accepts either scalars or arrays.

    Input
    -----
    rsky_AU : float or ndarray
        Sky-projected separation at mid-transit [AU].

    vsky_AU_day : float or ndarray
        Sky-projected relative velocity at mid-transit [AU/day].

    rho_kg_m3 : float
        Mean stellar density [kg/m^3].

    Mstar_Msun : float
        Stellar mass [Msun].

    Rp_Rstar : float
        Planet-to-star radius ratio.

    Output
    ------
    T14_day : float or ndarray
        First-to-fourth contact duration [day].

    T23_day : float or ndarray
        Second-to-third contact duration [day].

    depth_uniform : float or ndarray
        Uniform-disk geometric depth.

    b : float or ndarray
        Impact parameter.

    Rstar_AU : float
        Stellar radius [AU].
    """

    scalar_input = np.ndim(rsky_AU) == 0 and np.ndim(vsky_AU_day) == 0

    rsky_AU = np.asarray(rsky_AU, dtype=float)
    vsky_AU_day = np.asarray(vsky_AU_day, dtype=float)

    rsky_AU, vsky_AU_day = np.broadcast_arrays(rsky_AU, vsky_AU_day)

    Mstar_kg = float(Mstar_Msun) * MSUN
    rho_kg_m3 = float(rho_kg_m3)

    Rstar_m = (3.0 * Mstar_kg / (4.0 * np.pi * rho_kg_m3))**(1.0/3.0)
    Rstar_AU = Rstar_m / AU

    k = float(Rp_Rstar)

    vsky_AU_day = np.abs(vsky_AU_day)

    b = rsky_AU / Rstar_AU

    T14_day = np.zeros_like(b, dtype=float)
    T23_day = np.zeros_like(b, dtype=float)
    depth_uniform = np.zeros_like(b, dtype=float)

    # Valid velocity
    valid_v = vsky_AU_day > 0.0

    # Any transit, including grazing
    has_transit = valid_v & (b < 1.0 + k)

    # Full, non-grazing transit
    full_transit = has_transit & (b < abs(1.0 - k))

    # T14
    chord14 = np.zeros_like(b, dtype=float)
    chord14[has_transit] = (
        2.0 * Rstar_AU *
        np.sqrt(np.maximum((1.0 + k)**2 - b[has_transit]**2, 0.0))
    )

    T14_day[has_transit] = chord14[has_transit] / vsky_AU_day[has_transit]

    # T23
    chord23 = np.zeros_like(b, dtype=float)
    chord23[full_transit] = (
        2.0 * Rstar_AU *
        np.sqrt(np.maximum((1.0 - k)**2 - b[full_transit]**2, 0.0))
    )

    T23_day[full_transit] = chord23[full_transit] / vsky_AU_day[full_transit]

    # Uniform-disk depth, including grazing geometry
    depth_uniform = uniform_disk_depth_from_b_k_vectorized(b, k)

    # No transit should have zero depth
    depth_uniform[~has_transit] = 0.0

    if scalar_input:
        return (float(T14_day),
                float(T23_day),
                float(depth_uniform),
                float(b),
                float(Rstar_AU))

    return T14_day, T23_day, depth_uniform, b, Rstar_AU

 

def transit_geom_from_rsky_rho_eomega(
    rsky_AU,
    rho_kg_m3,
    P_day,
    Mstar_Msun,
    ecc,
    omega_deg,
    Rp_over_Rstar=None,
    Rp_Rearth=None,
    Rp_Rjup=None,
):
    """
    Compute transit geometry and simple observables from:
      - sky-projected separation at mid-transit, rsky_AU
      - stellar density, rho_kg_m3
      - period, P_day
      - stellar mass, Mstar_Msun
      - eccentricity, ecc
      - argument of periastron, omega_deg

    Returns:
      inc_deg
      b
      a_over_R
      Rstar_AU
      depth_relflux
      T14_day
      T23_day

    Notes:
      - depth is geometric, no limb darkening, delta = (Rp/Rstar)^2
      - durations are in days
      - for grazing/non-transiting cases, durations are set consistently
    """

    P_sec = P_day * DAY
    Mstar_kg = Mstar_Msun * MSUN
    omega = np.deg2rad(omega_deg)

    # Stellar radius from mean density and stellar mass
    Rstar_m = (3.0 * Mstar_kg / (4.0 * np.pi * rho_kg_m3)) ** (1.0 / 3.0)
    Rstar_AU = Rstar_m / AU

    # Semi-major axis in units of stellar radii
    a_over_R = ((G * rho_kg_m3 * P_sec**2) / (3.0 * np.pi)) ** (1.0 / 3.0)

    # Impact parameter from sky-projected separation at mid-transit
    b = rsky_AU / Rstar_AU

    # Eccentric correction
    fac = (1.0 + ecc * np.sin(omega)) / (1.0 - ecc**2)

    cosi = (b * fac) / a_over_R
    cosi = np.clip(cosi, -1.0, 1.0)
    inc_deg = np.degrees(np.arccos(cosi))
    sini = np.sqrt(np.maximum(0.0, 1.0 - cosi**2))

    # Planet radius ratio
    n_set = sum(x is not None for x in [Rp_over_Rstar, Rp_Rearth, Rp_Rjup])
    if n_set != 1:
        raise ValueError("Provide exactly one of Rp_over_Rstar, Rp_Rearth, or Rp_Rjup.")

    if Rp_over_Rstar is not None:
        k = Rp_over_Rstar
    elif Rp_Rearth is not None:
        Rp_m = Rp_Rearth * REARTH
        k = Rp_m / Rstar_m
    else:
        Rp_m = Rp_Rjup * RJUP
        k = Rp_m / Rstar_m

    # Simple geometric depth
    depth_relflux = k**2

    # Transit existence checks
    # Transit if b < 1 + k
    if b >= (1.0 + k):
        T14_day = 0.0
        T23_day = 0.0
        return inc_deg, b, a_over_R, Rstar_AU, depth_relflux, T14_day, T23_day

    # Common prefactor for eccentric orbit
    prefac = (P_day / np.pi) * (np.sqrt(1.0 - ecc**2) / (1.0 + ecc * np.sin(omega)))

    denom = a_over_R * sini

    # Numerical guard
    if denom <= 0.0:
        T14_day = 0.0
        T23_day = 0.0
        return inc_deg, b, a_over_R, Rstar_AU, depth_relflux, T14_day, T23_day

    # T14
    arg14_num = np.sqrt(np.maximum(0.0, (1.0 + k)**2 - b**2))
    x14 = np.clip(arg14_num / denom, -1.0, 1.0)
    T14_day = prefac * np.arcsin(x14)

    # T23, only if full transit exists
    if b < (1.0 - k):
        arg23_num = np.sqrt(np.maximum(0.0, (1.0 - k)**2 - b**2))
        x23 = np.clip(arg23_num / denom, -1.0, 1.0)
        T23_day = prefac * np.arcsin(x23)
    else:
        T23_day = 0.0

    return inc_deg, b, a_over_R, Rstar_AU, depth_relflux, T14_day, T23_day

 
def inc_from_rsky_rho_eomega(rsky_AU, rho_kg_m3, P_day, Mstar_Msun, ecc, omega_deg):
    #if ecc < 0.0 or ecc >= 1.0:
    #    raise ValueError("ecc must satisfy 0 <= ecc < 1.")

    P_sec = P_day * DAY
    Mstar_kg = Mstar_Msun * MSUN
    omega = np.deg2rad(omega_deg)

    # Stellar radius from mean density and stellar mass
    Rstar_m = (3.0 * Mstar_kg / (4.0 * np.pi * rho_kg_m3))**(1.0/3.0)
    Rstar_AU = Rstar_m / AU

    # Semi-major axis in units of stellar radii (Kepler + rho)
    a_over_R = ((G * rho_kg_m3 * P_sec**2) / (3.0 * np.pi))**(1.0/3.0)

    # Impact parameter from TTVFast sky-projected separation at mid-transit
    b = rsky_AU / Rstar_AU

    # Eccentric correction for transit geometry
    fac = (1.0 + ecc * np.sin(omega)) / (1.0 - ecc**2)

    cosi = (b * fac) / a_over_R
    cosi = np.clip(cosi, -1.0, 1.0)

    inc_deg = np.degrees(np.arccos(cosi))

    return inc_deg, b, a_over_R, Rstar_AU



def inc_from_rsky_rho(rsky_AU, rho_kg_m3, P_day, Mstar_Msun):

    P_sec = P_day * DAY
    Mstar = Mstar_Msun * MSUN

    # stellar radius from density
    Rstar_m = (3.0 * Mstar / (4.0 * np.pi * rho_kg_m3))**(1.0/3.0)
    Rstar_AU = Rstar_m / AU

    # a/R*
    a_over_R = ((G * rho_kg_m3 * P_sec**2) / (3.0*np.pi))**(1.0/3.0)

    # impact parameter
    b = rsky_AU / Rstar_AU

    x = b / a_over_R
    x = np.clip(x, -1.0, 1.0)

    inc_deg = np.degrees(np.arccos(x))

    return inc_deg, b, a_over_R


def get_rho_from_aR(aR, P, aR_e=None, P_e=None):
    """
    Calculate stellar density (rho) from scaled semi-major axis (a/R*) and period.

    Parameters
    ----------
    aR : float
        Scaled semi-major axis (a/R*).

    P : float
        Orbital period in days.

    aR_e : float, optional
        Uncertainty in a/R*.

    P_e : float, optional
        Uncertainty in orbital period in days.

    Returns
    -------
    float or tuple
        If uncertainties are not provided:
            rho : float
                Stellar density in kg/m^3.

        If uncertainties are provided:
            rho : float
            rho_e : float
                Uncertainty in stellar density.
    """

    import numpy as np

    G = 6.67408e-11

    P = P * 86400.0

    rho = (3.0 * np.pi / (G * P**2.0)) * aR**3.0

    if aR_e is not None and P_e is not None:

        P_e = P_e * 86400.0

        # derivatives
        drho_daR = 3.0 * rho / aR
        drho_dP = -2.0 * rho / P

        rho_e = np.sqrt((drho_daR * aR_e)**2 + (drho_dP * P_e)**2)

        return rho, rho_e

    return rho



def get_aR_from_rho(rho, P, rho_e=None, P_e=None):
    """
    Calculate the scaled semi-major axis (a/R*) from stellar density and orbital period.
    
    Parameters
    ----------
    rho : float
        Stellar density in kg/m^3.
    P : float
        Orbital period in days.
        
    rho_e : float, optional
        Uncertainty in stellar density in kg/m^3.
    P_e : float, optional
        Uncertainty in orbital period in seconds.
    
    Returns
    -------
    float or tuple
        If uncertainties are not provided, returns:
            aR : float
                Scaled semi-major axis (a/R*).
        If uncertainties are provided, returns:
            aR : float
                Scaled semi-major axis (a/R*).
            aR_e : float
                Uncertainty in scaled semi-major axis.
    
    Notes
    -----
    The formula used is:
        a/R* = [ (G * rho * P^2) / (3 * pi) ]^(1/3)
    where G is the gravitational constant.
    
    The uncertainty is estimated using standard error propagation.
    """
    G = 6.67408e-11  # gravitational constant in m^3 kg^-1 s^-2
    P = P*86400.0
    aR = ((rho * G * P**2.0) / (3.0 * np.pi))**(1.0/3.0)

    if rho_e is not None and P_e is not None:
        # Partial derivatives
        
        P_e = P_e*86400.0
        
        dadrho = (1.0/3.0) * ((G * P**2.0) / (3.0 * np.pi))**(1.0/3.0) * rho**(-2.0/3.0)
        dadP = (2.0/3.0) * ((G * rho) / (3.0 * np.pi))**(1.0/3.0) * P**(-1.0/3.0)

        # Error propagation
        aR_e = np.sqrt((dadrho * rho_e)**2.0 + (dadP * P_e)**2.0)

        return aR, aR_e

    return aR


def get_mass(K, P, ecc, i, Stellar_mass):
    '''Returns the mass in

    Parameters
    ----------
    K : float
        Semi-amplitude of the RV signal in [m/s]
    P : float
        Period of the planet in [d]
    ecc : float
        eccentricity
    incl: float
        inclination in [deg]
    Stellar_mass: float
        Primary mass in [Msol]

    Returns
        float
        planet mass in [Mj]
    -------
    '''
    T = P*86400.0
    THIRD = 1.0/3.0

    GMSUN = 1.32712440018e20
    msini = (T/(2.0*np.pi*GMSUN))**THIRD * K * Stellar_mass**(2./3) * np.sqrt(1.0-ecc**2.0)

    msini = msini/np.sin(np.radians(i))*1047.5654817267318

    return msini

def calculate_ttv_period(P_in, P_out):
    """
    Calculate the TTV period based on the periods of two near-resonant planets.
    
    Parameters:
    P_in (float): Orbital period of the inner planet.
    P_out (float): Orbital period of the outer planet.
    
    Returns:
    float: TTV period (P_TTV)
    """
    # Estimate the resonance order (j) based on the ratio of periods
    j = round(P_out / P_in)

 
    # Calculate the TTV period
    P_TTV = abs((j - 1) / P_in - j / P_out)**-1
    
    return P_TTV
    
    
def a_from_K_in_mas(P,K,e,i,parallax):
    """
    Calculate the semi-major axis 'a' given K, parallax, period P, eccentricity e, and inclination i.

    Parameters:
    K : float
        Radial velocity amplitude (m/s)
    parallax : float
        Parallax in miliarcseconds
    P : float
        Orbital period in days
    e : float
        Eccentricity
    i : float
        Inclination in degerees

    Returns:
    a : float
        Semi-major axis in mas
    """    
    i = np.radians(i)

    a = (K * P *86400.0 * np.sqrt(1.0-e**2)* parallax) / (2*np.pi*np.sin(i)*1.496e11)

    return a    

def get_gravity(m_p, r_p):
    '''Returns the gravity in

    Parameters
    ----------
    m_p in kg
    r_p in m
    Returns
        float
        g mass in in m/s^2
    -------
    '''
    G = 6.67e-11
    return G*m_p/r_p**2


def a_to_P(a,m0):
    GMSUN = 1.32712440018e20
    AU=1.49597892e11
    T = np.sqrt( (a*AU)**3.0 * (2.0*np.pi)**2.0 /(GMSUN*(m0)))
    T = T /86400.0
    return T

def P_to_a(P,m0):
    GMSUN = 1.32712440018e20
    AU=1.49597892e11
    P = P * 86400.0
    a = ((P**2.0 * (GMSUN*(m0)))/(4.0*(np.pi)**2.0))**(1.0/3.0)

    return a/AU



def get_transit_times(tr_res, p, t0, precise = False, verbose=False):
    '''Returns linear transit times (approximate!)

    Parameters
    ----------
    -------
    '''
    #t_ind = np.where(np.logical_and(t >= t0-0.17, t <= t0+0.17))

    #transits = t0 + p

    t = tr_res[2][0]
    f = tr_res[2][3]
    range_ = int(max(t-t[0])/p)+1

    tr_index = []
    tr_t0    = []
     #tran_t0 = t0
    for i in range(range_):
        tran_t0 = t0 + p*i

        tr_ind = np.where(np.logical_and(t >= tran_t0-0.17, t <= tran_t0+0.17))
        #print(i,tran_t0)
        if len(t[tr_ind]) !=0:
            minn = np.argmin(f[tr_ind])
            tr_index.append(i)
            tr_t0.append(t[tr_ind][minn])

    if precise:
        t = tr_res[3][0]
        f = tr_res[3][1]
        tr_t0    = []

        for i in tr_index:
            tran_t0 = t0 + p*i

            tr_ind = np.where(np.logical_and(t >= tran_t0-0.15, t <= tran_t0+0.15))
            minn = np.argmin(f[tr_ind])
            tr_t0.append(t[tr_ind][minn])

    if verbose == True:
        for i in range(len(tr_index)):
            print(tr_index[i],tr_t0[i])


    return [tr_index, tr_t0]


############################ density ###########################################
def get_density(m_p, r_p, d_m_p=None, d_r_p=None):
    """
    # m_p in Sol
    # r_p in Sol

    Warning! You must do:  rho = get_density(pl_mass_kg*1000.0, pl_radii_m*100.0) to get
    rho in [g cm^-3]
    """

    solarrad2m = 6.957e8
    solarmass2kg = 1.9891e30

    r_p = r_p*solarrad2m
    m_p = m_p*solarmass2kg
    const =  (4.0/3.0) * np.pi

    volume =  const * (r_p ** 3.0)
    rho = m_p/volume

    if d_m_p != None:
        d_m_p = d_m_p*solarmass2kg
        d_r_p = d_r_p*solarrad2m
        delta_rho = np.sqrt( (d_m_p**2*r_p**2) + (9*d_r_p**2 * m_p**2)   ) / (const*r_p**4)
    else:
        delta_rho = 0.0

    return rho,delta_rho



####################### mass_semimajor #########################################
def get_mass_a_samples(K, P, ecc, incl, m_s, mass_type="J"):
    """Physical Jacobi masses/axes for a posterior sample, in caller order.

    Keep the legacy nine-slot output; units are Mjup (J), Mearth (E), or Msun.
    """
    npl = len(K)
    rows = np.column_stack((K, P, ecc, np.zeros(npl), np.zeros(npl), incl, np.zeros(npl)))
    masses, axes = mass_a_from_jacobi_el(rows, npl, m_s)
    if mass_type == "E":
        masses *= 317.82838
    elif mass_type != "J":
        masses /= 1047.5654817267318
    return np.pad(masses, (0, 9 - npl)), np.pad(axes, (0, 9 - npl))

####################### mass_semimajor #########################################
def get_mass_a(obj, mass_type="J"):
    """Physical Jacobi masses/axes in the active planet slots."""
    slots = [i for i, active in enumerate(obj.use_planet) if active]
    masses = np.zeros(9)
    axes = np.zeros(9)
    rows = [[obj.K[i], obj.P[i], obj.e[i], obj.w[i], obj.M0[i], obj.i[i], obj.Node[i]]
            for i in slots]
    physical, jacobi_axes = mass_a_from_jacobi_el(rows, len(slots), obj.st_mass[0])
    if mass_type == "E":
        physical *= 317.82838
    elif mass_type != "J":
        physical /= 1047.5654817267318
    masses[slots], axes[slots] = physical, jacobi_axes
    return masses, axes




def LoadSession(input_file, template_session = None):

    try:
        file_pi = open(input_file, 'rb')
        fit_new = dill.load(file_pi) #, encoding='latin1'
        file_pi.close()
    except (UnicodeDecodeError, ImportError, KeyError) as e:
        py3_ses = convert_Session_to_Py3(input_file)

        file_pi = open(py3_ses, 'rb')
        fit_new = dill.load(file_pi) #, encoding='latin1'
        file_pi.close()

    if template_session != None:
        fit_new = check_for_missing_instances(template_session,fit_new)
    #self.check_for_missing_instances(fit_new)

    #self.check_settings()
    check_temp_RV_file(fit_new)
    return fit_new

def convert_Session_to_Py3(old_ses):
    """
    Convert a Python 2 sesion to Python 3
    """
    # Make a name for the new pickle
    new_ses = os.path.splitext(os.path.basename(old_ses))[0]+"_p3.ses"

    # Convert Python 2 "ObjectType" to Python 3 object
    dill._dill._reverse_typemap["ObjectType"] = object

    # Open the pickle using latin1 encoding
    with open(old_ses, "rb") as f:
        loaded = dill.load(f, encoding="latin1")
    f.close()

    # Re-save as Python 3 pickle
    with open(new_ses, "wb") as outfile:
        dill.dump(loaded, outfile)
    outfile.close()

    return new_ses


def fix_old_to_session_tra(old_ses):

        ################  Some problems with old sessions are fixed here. TB removed later on.# #############

        if len(old_ses.tra_data_sets) == 20:
            return old_ses

        else:
            for i in range(10):
                old_ses.tra_data_sets[10 + i] = []
                old_ses.tra_data_sets_init[10 + i] = []
                old_ses.ld_m = ["quadratic"]*20    #limb darkening model

                old_ses.ld_u[10 + i] =  [0.12, 0.35 ]

                old_ses.ld_u_lin[10 + i] = [0.35]
                old_ses.ld_u_quad[10 + i] =  [0.12, 0.35 ]
                old_ses.ld_u_nonlin[10 + i] =  [0.55,0.12, 0.35,-0.11]

                old_ses.ld_u_lin_use[10 + i] =  [False]
                old_ses.ld_u_quad_use[10 + i] =  [False, False]
                old_ses.ld_u_nonlin_use[10 + i] =  [False, False,False, False]

                old_ses.ld_u_lin_err[10 + i] =  [[0.0,0.0]]
                old_ses.ld_u_quad_err[10 + i] =  [[0.0,0.0], [0.0,0.0]]
                old_ses.ld_u_nonlin_err[10 + i] = [[0.0,0.0], [0.0,0.0],[0.0,0.0], [0.0,0.0]]

                old_ses.ld_u_lin_bound[10 + i] =  np.array([[-1.0,1.0]])
                old_ses.ld_u_quad_bound[10 + i] =  np.array([[-1.0,1.0],[-1.0,1.0]])
                old_ses.ld_u_nonlin_bound[10 + i] =  np.array([[-1.0,1.0],[-1.0,1.0],[-1.0,1.0],[-1.0,1.0]])

                old_ses.ld_u_lin_norm_pr[10 + i] =  np.array([[0.1,0.05, False]])
                old_ses.ld_u_quad_norm_pr[10 + i] = np.array([[0.0,1.0, False],[0.0,1.0, False]])
                old_ses.ld_u_nonlin_norm_pr[10 + i] = np.array([[0.0,1.0, False],[0.0,1.0, False],[0.0,1.0, False],[0.0,1.0, False]])

                old_ses.ld_u_lin_jeff_pr[10 + i] =  np.array([[0.1,0.05, False]])
                old_ses.ld_u_quad_jeff_pr[10 + i] =  np.array([[0.0,1.0, False],[0.0,1.0, False]])
                old_ses.ld_u_nonlin_jeff_pr[10 + i] =  np.array([[0.0,1.0, False],[0.0,1.0, False],[0.0,1.0, False],[0.0,1.0, False]])

                old_ses.ld_u_lin_str[10 + i] =  [r'ld-quad-1$_%s$'%str(10 + i+1)]
                old_ses.ld_u_quad_str[10 + i] =  [r'ld-quad-1$_%s$'%str(10 + i+1),r'ld-quad-2$_%s$'%str(10 + i+1)]
                old_ses.ld_u_nonlin_str[10 + i] =  [r'ld-quad-1$_%s$'%str(10 + i+1),r'ld-quad-2$_%s$'%str(10 + i+1),r'ld-quad-3$_%s$'%str(10 + i+1),r'ld-quad-4$_%s$'%str(10 + i+1)]
                old_ses.ld_gr.append(10 + i)
                old_ses.ld_gr_ind.append(10 + i)


        for i in range(20):
            if len(old_ses.tra_data_sets[i]) != 0:

                if len(old_ses.tra_data_sets[i]) ==11:
                     old_ses.tra_data_sets[i] = np.insert(old_ses.tra_data_sets[i], 9, True)
                     old_ses.tra_data_sets_init[i] = np.insert(old_ses.tra_data_sets_init[i], 9, True)
                elif len(old_ses.tra_data_sets[i]) ==10:
                     old_ses.tra_data_sets[i] = np.insert(old_ses.tra_data_sets[i], 9, True)
                     old_ses.tra_data_sets_init[i] = np.insert(old_ses.tra_data_sets_init[i], 9, False)
                     old_ses.tra_data_sets[i] = np.insert(old_ses.tra_data_sets[i], 10, True)
                     old_ses.tra_data_sets_init[i] = np.insert(old_ses.tra_data_sets_init[i], 10, False)

        return old_ses

         ########################################################################################################


def find_close_elements(a, b, precision = 0.01):
    """Finds close elements in two arrays with diffrent sizes.

    Parameters
    ----------
    a : array of floats with dimention of N
        Description of parameter `a`.
    b : array of floats with dimention of M
        Description of parameter `b`.
    precision : threshold withing which two elements are considered the same.
        Description of parameter `precision`.

    Returns
    -------
    [array, array]
        returns two arrays with the elements that mathched withing the
        precision.
    """
    return [[x for x in a for i in b if abs(x - i) < precision], [x for x in b for i in a if abs(x - i) < precision]]

def custom_param_file_for_stability(max_time,time_step):

        ##### create the param.in file (change only the "t_max" and the "dt" for now) ######
    param_file = open('param.in', 'wb')

    max_time = float(max_time)*365.25 # make it is days

    param_file.write("""0.0d0 %s %s
%s %s

F T T T T F
0.001 50.0 50.0 -1. T
 bin.dat
unknown
"""%(max_time, time_step, max_time/1e4, max_time/1e3 ))

    param_file.close()
    return

def _default_cornerplot_labels(labels, save_flags):

    bestfit_labels      = ["median","mean","mode","best_samp","best_gui","none",
                           "mass","use_Me","use_Mj","use_Ms",
                           "semimajor","radius","use_Re","use_Rj","use_Rs","use_lambda", "use_ppm"]
    bestfit_labels_bool = [*save_flags,False,False,False,
                           True,False,False,False,False,True,False,False,False,False]


    # Reserve room for edited labels while preserving longer initial labels.
    settings             = {k: np.array([label, True], dtype=f"U{max(30, len(label))}")
                               for k, label in enumerate(labels)}
    for k in range(17):
        settings[bestfit_labels[k]] = bestfit_labels_bool[k]

    cornerplot_opt = {"bins":25,
                      "color":"k",
                      "reverse":False,
                      "upper":True,
                      "quantiles":68.3,
                      "levels":(0.6827, 0.9545,0.9973),
                      "smooth":1.0,
                      "smooth1d":1.0,
                      "plot_contours":True,
                      "show_titles":True,
                      "dpi":300,
                      "pad":15,
                      "labelpad":0,
                      "truth_color":'r',
                      "title_kwargs":{"fontsize": 12},
                      "scale_hist":True,
                      "fill_contours":True,
                      "no_fill_contours":False,
                      "plot_datapoints":False,
                      "stab_color":"r",
                      "stab_threshold":100000
                      }

    settings["cornerplot"] = cornerplot_opt
    settings["OrigLabels"] = dill.copy(labels)

    return settings


def ensure_cornerplot_labels(obj, type_plot="mcmc"):
    """Restore missing plot metadata only when usable posterior samples exist."""
    sampler = obj.mcmc_sampler if type_plot == "mcmc" else obj.ns_sampler
    samples = np.asarray(getattr(sampler, 'samples', []))
    if samples.ndim != 2 or not samples.size:
        return False
    if getattr(sampler, 'lbf', None) is None:
        labels = list(getattr(obj, 'e_for_mcmc', []))
        if len(labels) != samples.shape[1]:
            print("Cannot customize corner plot: sample columns do not match parameter labels.")
            return False
        prefix = 'mcmc' if type_plot == 'mcmc' else 'ns'
        flags = [getattr(obj, prefix + '_save_' + name, default)
                 for name, default in (('median', True), ('means', False),
                                       ('mode', False), ('maxlnL', False))]
        sampler.lbf = _default_cornerplot_labels(labels, flags)
    labels = list(sampler.lbf.get('OrigLabels', getattr(obj, 'e_for_mcmc', [])))
    if len(labels) != samples.shape[1]:
        print('Cannot customize corner plot: sample columns do not match original labels.')
        return False
    defaults = _default_cornerplot_labels(labels, [True, False, False, False])
    for key, value in defaults.items():
        sampler.lbf.setdefault(key, value)
    for column, label in enumerate(labels):
        entry = sampler.lbf[column]
        display = default_rv_units(entry[0])
        sampler.lbf[column] = np.array([display, entry[1]], dtype=f'U{max(30, len(display))}')
    return True


def add_mcmc_samples(obj,sampler):
    sampler.lbf = _default_cornerplot_labels(
        obj.e_for_mcmc, [obj.mcmc_save_median, obj.mcmc_save_means,
                        obj.mcmc_save_mode, obj.mcmc_save_maxlnL])
    obj.mcmc_sampler = sampler
    obj.sampler_saved = True


def add_ns_samples(obj,sampler):

    bestfit_labels      = ["median","mean","mode","best_samp","best_gui","none",
                           "mass","use_Me","use_Mj","use_Ms",
                           "semimajor","radius","use_Re","use_Rj","use_Rs","use_lambda","use_ppm"]
    bestfit_labels_bool = [obj.ns_save_median,obj.ns_save_means,obj.ns_save_mode,
                           obj.ns_save_maxlnL,False,False,False,
                           True,False,False,False,False,True,False,False,False,False]



    obj.ns_sampler= dill.copy(sampler.results)
    obj.ns_sampler.lbf     = {k: np.array([label, True], dtype=f"U{max(30, len(label))}")
                             for k, label in enumerate(obj.e_for_mcmc)}

    for k in range(17):
        obj.ns_sampler.lbf[bestfit_labels[k]] = bestfit_labels_bool[k]

    cornerplot_opt = {"bins":25,
                      "color":"k",
                      "reverse":False,
                      "upper":True,
                      "quantiles":68.3,
                      "levels":(0.6827, 0.9545,0.9973),
                      "smooth":1.0,
                      "smooth1d":1.0,
                      "plot_contours":True,
                      "show_titles":True,
                      "dpi":300,
                      "pad":15,
                      "labelpad":0,
                      "truth_color":'r',
                      "title_kwargs":{"fontsize": 12},
                      "scale_hist":True,
                      "fill_contours":True,
                      "no_fill_contours":False,
                      "plot_datapoints":False,
                      "stab_color":"r",
                      "stab_threshold":100000
                      }

    obj.ns_sampler.lbf["cornerplot"] = cornerplot_opt
    obj.ns_sampler.lbf["OrigLabels"] = dill.copy(obj.e_for_mcmc)

    #delattr(obj.ns_sampler, 'rstate')
    obj.sampler_saved=True
    return obj

def get_quad_model(x,y,a1,a2,a3):

    x = x - x[0]
    y = y + a1 + a2*x + a3*x**2
    return y

def get_airmass_model(x,y,a1,a2,a3):

    #x = x - x[0]
    #print(a1,a2,x[0:2],a2*x[0:2])
    y = y + a1 + a2*x + a3*x**2
    return y


def get_mode_of_samples(samples, nsamp):

    mode_samp = []
   # err1_samp = []
  #  err2_samp = []


    for i in range(nsamp):
        #ci = np.percentile(samples[:,i], [level, 100.0-level])
        #mmm = stats.binned_statistic(np.array([samples[:,i]]), axis=None)
        n, b = np.histogram(samples[:,i], bins=100)
        n = gaussian_filter(n, 1.0)
        x0 = np.array(list(zip(b[:-1], b[1:]))).flatten()
        y0 = np.array(list(zip(n, n))).flatten()
        k  = np.unravel_index(y0.argmax(),y0.shape)
        mode_samp.append(x0[k])
        #err1_samp.append(x0[k]- ci[0])
        #err2_samp.append(ci[1]- x0[k])
   # print el_str[i],'=', x0[k], "- %s"%(x0[k]-ci[0]), "+ %s"%(ci[1]  - x0[k] )
    return mode_samp #,err1_samp,err2_samp

def get_mean_of_samples(samples, nsamp):

    mean_samp = []
    for i in range(nsamp):
        mean_samp.append(np.mean(samples[:,i]))
    return mean_samp

def get_median_of_samples(samples, nsamp):

    median_samp = []
    for i in range(nsamp):
        median_samp.append(np.median(samples[:,i]))
    return median_samp

def get_MAD_of_samples(samples, nsamp):

    mad_samp = []
    for i in range(nsamp):
        mad_samp.append(np.mean(np.absolute(samples[:,i] - np.mean(samples[:,i]))))
    return mad_samp



def get_best_lnl_of_samples(samples,lnl, nsamp):

    best_ln_samp = []
    lnL_best_idx = np.argmax(lnl)
    lnL_best = lnl[lnL_best_idx]


    for i in range(nsamp):

        minlnL = samples[lnL_best_idx,i]
        best_ln_samp.append(minlnL)


    return best_ln_samp,lnL_best #,err1_samp,err2_samp




def _corner_stellar_draws(value, uncertainty, count):
    """Gaussian draws for legacy scalar or current [lower, upper] errors.

    Corner-derived quantities use the larger error magnitude as a conservative
    symmetric sigma. Zero errors produce a constant stellar parameter.
    """
    center = np.asarray(value, dtype=float)
    errors = np.asarray(uncertainty, dtype=float).reshape(-1)
    if center.size != 1 or not np.isfinite(center).all():
        raise ValueError('Stellar parameter must be a finite scalar')
    if errors.size not in (1, 2) or not np.isfinite(errors).all():
        raise ValueError('Stellar uncertainty must be a finite scalar or lower/upper pair')
    sigma = float(np.max(np.abs(errors)))
    return np.random.normal(loc=float(center.reshape(-1)[0]), scale=sigma, size=count)


def _corner_lambda_samples(obj, samples, labels, slot):
    """Use sampled angles when present, otherwise the fixed orbital element."""
    letter = 'bcdefghij'[slot]
    labels = list(labels)
    def column(candidates, fixed):
        for label in candidates:
            if label in labels:
                return np.asarray(samples[:, labels.index(label)], dtype=float)
        return np.full(len(samples), float(fixed))
    if obj.hkl:
        return column([getattr(obj, 'lamb_str', {}).get(slot, ''),
                       r'$\lambda_%s$ [deg]' % letter], obj.lamb[slot]) % 360.
    mean_anomaly = column([getattr(obj, 'M0_str', {}).get(slot, ''),
                          'MA$_%s$ [deg]' % letter, 'MA$_%s$' % letter,
                          'MA$_%s [deg]' % letter], obj.M0[slot])
    omega = column([getattr(obj, 'w_str', {}).get(slot, ''),
                    r'$\omega_%s$ [deg]' % letter], obj.w[slot])
    return (mean_anomaly + omega) % 360.


def _corner_mass_axes(obj, samples, labels, stellar_masses):
    slots = [i for i in range(9) if obj.use_planet[i]]
    elements = np.empty((len(samples), len(slots), 7))
    labels = list(labels)
    def values(root, slot):
        label = getattr(obj, root+'_str')[slot]
        return samples[:, labels.index(label)] if label in labels else np.full(len(samples), getattr(obj, root)[slot])
    for packed, slot in enumerate(slots):
        roots = ('K','P','e_sinw','e_cosw','lamb','i','Node') if obj.hkl else ('K','P','e','w','M0','i','Node')
        elements[:,packed,:] = np.column_stack([values(root,slot) for root in roots])
    if obj.hkl:
        h,k,longitude = (elements[:,:,j].copy() for j in (2,3,4))
        elements[:,:,2] = np.hypot(h,k)
        elements[:,:,3] = np.where(elements[:,:,2]>0,np.degrees(np.arctan2(h,k)) % 360.,0.)
        elements[:,:,4] = (longitude-elements[:,:,3]) % 360.
    if obj.copl_incl:
        for row in elements:
            row[:,5] = row[np.argmin(row[:,1]),5]
    masses = np.full((len(samples),len(slots)),np.nan)
    axes = masses.copy()
    for row, (orbit, stellar_mass) in enumerate(zip(elements, stellar_masses)):
        try:
            masses[row], axes[row] = mass_a_from_jacobi_el(orbit,len(slots),stellar_mass)
        except ValueError:
            continue
    return slots,masses,axes


def cornerplot(obj, level=(100.0-68.3)/2.0, type_plot = 'mcmc', **kwargs):

    statistics_only = kwargs.pop('return_statistics', False)

    #obj = dill.copy(copied_obj)
    '''Generates a corner plot visualizing the mcmc samples. Optionally samples can be read from a file.'''
    #self.mcmc_sample_file = 'mcmc_samples'+'_%s'%mod
    #self.corner_plot_file = 'cornerplot.png'


    if not ensure_cornerplot_labels(obj, type_plot):
        return

    if type_plot == 'mcmc':

        #### load the samples, labels and lnL values
        ln      = dill.copy(np.hstack(obj.mcmc_sampler.lnL))
        samples = dill.copy(np.array(obj.mcmc_sampler.samples))
        #labels  = dill.copy(obj.e_for_mcmc)
        labels  = dill.copy(obj.mcmc_sampler.lbf["OrigLabels"])
        mod_labels  = dill.copy(obj.mcmc_sampler.lbf)

        if mod_labels['mean'] ==True:
            best_fit_par = obj.mcmc_stat["mean"]
        elif mod_labels['median']==True:
            best_fit_par = obj.mcmc_stat["median"]
        elif mod_labels['mode']==True:
            best_fit_par = obj.mcmc_stat["mode"]
        elif mod_labels['best_samp']==True:
            best_fit_par = obj.mcmc_stat["best"]
        elif mod_labels['none']==True:
            best_fit_par = np.array([0]*len(labels), dtype=np.float64) #obj.par_for_mcmc

        elif mod_labels['best_gui']==True:
            best_fit_par = obj.par_for_mcmc
        else:
            best_fit_par = obj.par_for_mcmc

        cornerplot_opt = dill.copy(obj.mcmc_sampler.lbf["cornerplot"])


    elif type_plot == 'nest':

        #### load the samples, labels and lnL values
        #ln      = dill.copy(obj.ns_sampler.results.logl)
       # samples = dill.copy(np.array(obj.ns_sampler.results.samples))
        ln      = dill.copy(obj.ns_sampler.logl)
        samples = dill.copy(np.array(obj.ns_sampler.samples))
        #labels  = dill.copy(obj.e_for_mcmc)
        labels  = dill.copy(obj.ns_sampler.lbf["OrigLabels"])
        mod_labels  = dill.copy(obj.ns_sampler.lbf)

        if mod_labels['mean'] ==True:
            best_fit_par = obj.nest_stat["mean"]
        elif mod_labels['median']==True:
            best_fit_par = obj.nest_stat["median"]
        elif mod_labels['mode']==True:
            best_fit_par = obj.nest_stat["mode"]
        elif mod_labels['best_samp']==True:
            best_fit_par = obj.nest_stat["best"]
        elif mod_labels['none']==True:
            best_fit_par = np.array([0]*len(labels), dtype=np.float64) #obj.par_for_mcmc

        elif mod_labels['best_gui']==True:
            best_fit_par = obj.par_for_mcmc
        else:
            best_fit_par = obj.par_for_mcmc

        cornerplot_opt = dill.copy(obj.ns_sampler.lbf["cornerplot"])

    else:
        return




    ############### make "Gaussan" samples of the stellar parameters ##############
    #m_s   = np.random.normal(loc=obj.stellar_mass,      scale=obj.stellar_mass_err,      size=len(samples[:,0]))
    #r_s = _corner_stellar_draws(obj.stellar_radius, obj.stellar_radius_err, len(samples))
   # L_s   = np.random.normal(loc=obj.stellar_luminosity,scale=obj.stellar_luminosity_err,size=len(samples[:,0]))
   # vsini = np.random.normal(loc=obj.stellar_vsini,     scale=obj.stellar_vsini_err,     size=len(samples[:,0]))
   #

    ######### define new samples, labels and best-fit params to be refilled again
    ######### with masses, semi-major axes, etc. (to be re-worked).

    if statistics_only:
        best_fit_par = np.zeros(len(labels))

    m_s     = []
    r_s     = []
    samp    = []
    samp_labels =  []
    samp_best_fit_par = []
    #print(best_fit_par)
    for i in range(len(labels)):

        ss = np.hstack(samples[:,i])

        if mod_labels['use_ppm']:
            if 'transit' in labels[i]:
                ss = ss * 1000000.0

        samp.append(ss)
        samp_labels.append(labels[i])
        samp_best_fit_par.append(best_fit_par[i])


    index_to_keep = []
    for i in range(len(labels)):
        samp_labels[i] = mod_labels[i][0]
        if mod_labels[i][1] == 'True':
            index_to_keep.append(i)


    if statistics_only:
        index_to_keep = list(range(len(labels)))

    samp              = [samp[i] for i in index_to_keep]
    samp_labels       = [samp_labels[i] for i in index_to_keep]
    samp_best_fit_par = [samp_best_fit_par[i] for i in index_to_keep]

    letters = ['b','c','d','e','f','g','h','i','j',] #... For the planets

    for option, slot, label, draws in derived_hkl_columns(obj, samples, labels, mod_labels):
        samp.append(draws)
        samp_labels.append(label)
        if mod_labels['best_gui']:
            h, k = obj.e_sinw[slot], obj.e_cosw[slot]
            omega = np.degrees(np.arctan2(h, k)) % 360.
            value = np.hypot(h, k) if option == 'derive_e' else omega if option == 'derive_omega' else (obj.lamb[slot]-omega) % 360.
        elif mod_labels['mean'] or mod_labels['median']:
            value = (np.nanmean if mod_labels['mean'] else np.nanmedian)(draws)
        else:
            value = draws[np.nanargmax(ln)] if np.isfinite(ln).any() else np.nanmedian(draws)
        samp_best_fit_par.append(value)

    if mod_labels['use_lambda']:
        for i in range(9):
            if not bool(obj.use_planet[i]):
                continue
            let = letters[i]
            lambda_label = r'$\lambda_%s$ [deg]' % let
            # h/k/lambda posteriors already contain the longitude column.
            existing = [lambda_label, getattr(obj, 'lamb_str', {}).get(i, '')]
            if any(label and label in labels for label in existing):
                continue
            lambda_ = _corner_lambda_samples(obj, samples, labels, i)
            samp.append(lambda_)
            samp_labels.append(lambda_label)
            if mod_labels['mean']:
                samp_best_fit_par.append(np.mean(lambda_))
            elif mod_labels['median']:
                samp_best_fit_par.append(np.median(lambda_))
            elif mod_labels['best_gui']:
                samp_best_fit_par.append(float(obj.lamb[i]) if obj.hkl else (obj.M0[i]+obj.w[i]) % 360.)
            else:
                samp_best_fit_par.append(lambda_[np.argmax(ln)])


    if mod_labels['mass'] or mod_labels['semimajor']:
        m_s = _corner_stellar_draws(obj.stellar_mass, obj.stellar_mass_err, len(samples))
        slots, masses, semimajor = _corner_mass_axes(obj, samples, labels, m_s)
        M_fact, mass_lab = (317.82838, r'[M$_\oplus$]') if mod_labels['use_Me'] else ((1./1047.5654817267318, r'[M$_\odot$]') if mod_labels['use_Ms'] else (1., r'[M$_{\rm Jup.}$]'))
        for option, values, factor, symbol, unit in [
                ('mass', masses, M_fact, 'm', mass_lab),
                ('semimajor', semimajor, 1., 'a', '[au]')]:
            if not mod_labels[option]:
                continue
            for packed, slot in enumerate(slots):
                column = values[:, packed]*factor
                samp.append(column)
                samp_labels.append(f'{symbol}$_{letters[slot]}$ {unit}')
                if mod_labels['mean'] or mod_labels['median']:
                    value = (np.nanmean if mod_labels['mean'] else np.nanmedian)(column)
                elif mod_labels['best_gui']:
                    value = (obj.masses if option == 'mass' else obj.semimajor)[slot]*factor
                else:
                    value = column[np.nanargmax(ln)] if np.isfinite(ln).any() else np.nanmedian(column)
                samp_best_fit_par.append(value)


    if mod_labels['radius']:
        r_s = _corner_stellar_draws(obj.stellar_radius, obj.stellar_radius_err, len(samples))

        for i in range(9):
            if not bool(obj.use_planet[i]):
                continue
            let = letters[i]

            if not 'R/$R_\star$ $%s$'%let in labels:
                continue

            if mod_labels['use_Re']:
                R_fact = 109.076
                rad_lab = r'[R$_\oplus$]'
            elif mod_labels['use_Rj']:
                R_fact = 9.955201593
                rad_lab = r'[R$_{\rm Jup.}$]'
            elif mod_labels['use_Rs']:
                R_fact = 1.0
                rad_lab = r'[R$_\odot$]'

            rad = np.hstack(samples[:,[ii for ii, j in enumerate(labels) if j == 'R/$R_\star$ $%s$'%let]]) * R_fact
            samp.append(np.array(rad*r_s))
            samp_labels.append(r'R$_%s$ %s'%(let,rad_lab))


            if mod_labels['mean']:
                samp_best_fit_par.append(np.mean(rad)*np.mean(r_s))
            elif mod_labels['median']:
                samp_best_fit_par.append(np.median(rad)*np.median(r_s))
            else:
                samp_best_fit_par.append(rad[np.argmax(ln)]*obj.stellar_radius)


#    if mod_labels['gravity']:

#        if len(r_s) == 0:
#            r_s = _corner_stellar_draws(obj.stellar_radius, obj.stellar_radius_err, len(samples))


    ################### Transpose is needed for the cornerplot. ###################

    samples_ = np.transpose(samp)
    if statistics_only:
        table_samples, _ = fold_plot_angles(samples_, samp_labels)
        stable = (samples[:, list(labels).index('max. time')] >= cornerplot_opt['stab_threshold']) if 'max. time' in labels else None
        return dict(samples=table_samples, labels=samp_labels, stable=stable, options=cornerplot_opt)
    usable = np.isfinite(samples_).any(axis=0)
    for j in np.flatnonzero(~usable):
        print('Skipping undefined corner parameter: '+samp_labels[j])
    samples_ = samples_[:, usable]
    samp_labels = [label for j,label in enumerate(samp_labels) if usable[j]]
    samp_best_fit_par = np.asarray(samp_best_fit_par)[usable]
    finite_rows = np.isfinite(samples_).all(axis=1)
    if not finite_rows.all():
        print(f'Corner plot: excluding {np.count_nonzero(~finite_rows)} rows with undefined selected parameters.')
    samples_ = samples_[finite_rows]
    samples = samples[finite_rows]
    ln = np.asarray(ln)[finite_rows]
    if not samples_.shape[1] or len(samples_) < 2:
        print('Not enough finite samples/parameters for a corner plot.')
        return
    #labels = samp_labels
    best_fit_par =samp_best_fit_par
    angle_summary = 'mean' if mod_labels['mean'] else 'median' if mod_labels['median'] else None
    samples_, samp_best_fit_par = fold_plot_angles(samples_, samp_labels,
                                                   samp_best_fit_par, angle_summary)
    samp = list(samples_.T)
    best_fit_par = samp_best_fit_par

    ################### Verbose output (TBD). ###################

    verbose = True
    level = (100.0-68.3)/2.0
    quantiles = None
    title_quantiles = None    
    level_q = (100.0-cornerplot_opt["quantiles"])/2.0    
    
    if cornerplot_opt["quantiles"] != None:
        quantiles = [level_q/100.0, 0.5, 1.0-level_q/100.0]

    if cornerplot_opt["show_titles"] != None:
        title_quantiles = [level_q/100.0, 0.5, 1.0-level_q/100.0]    
 


    summary_mode = 'mean' if mod_labels['mean'] else 'median' if mod_labels['median'] else 'best'
    if verbose:
        print_summary(samples_, samp_labels, summary_mode, samp_best_fit_par, None if mod_labels['best_gui'] else ln,
                      rounded=cornerplot_opt.get('round_statistics', True),
                      show_mad=cornerplot_opt.get('show_mad', False))

    range_stab =[]
#    for i in range(len(samp)):
#        range_stab.append([0.0,1.0])

    if mod_labels['none']==True:
        best_fit_par = None
        
        



    if "max. time" in labels:
        samples_stab = []

        samples_ = np.transpose(samples_)
        for i in range(len(samples_)):
            samples_stab.append(samples_[i][np.where(samples[:,-1]>= cornerplot_opt["stab_threshold"])])
#        print(max(samples[:,-1])-1000.0,max(samples[:,-1]))
            range_stab.append([min(samples_[i]),max(samples_[i])])
        N_samp = len(np.atleast_1d(samples_stab[0]))



        if verbose:
            stable_mask = samples[:, -1] >= cornerplot_opt['stab_threshold']
            print_summary(np.asarray(samples_stab).T, samp_labels, summary_mode,
                          samp_best_fit_par, None if mod_labels['best_gui'] else np.asarray(ln)[stable_mask],
                          rounded=cornerplot_opt.get('round_statistics', True),
                          show_mad=cornerplot_opt.get('show_mad', False), stable=True)

        samples_stab = np.transpose(samples_stab)
        samples_ = np.transpose(samples_)



        if N_samp > len(samp_best_fit_par):

            fig = corner.corner(samples_stab,
                            range = range_stab,
                            bins=cornerplot_opt["bins"],
                            color=cornerplot_opt["stab_color"],
                            reverse=cornerplot_opt["reverse"],
                            upper=cornerplot_opt["upper"],
                            labels=samp_labels,
                            quantiles=quantiles,
                            title_quantiles=title_quantiles,
                            levels=(0.6827, 0.9545,0.9973),
                            smooth=cornerplot_opt["smooth"],
                            smooth1d=cornerplot_opt["smooth1d"],
                            plot_contours=cornerplot_opt["plot_contours"],
                            show_titles=cornerplot_opt["show_titles"],
                            truths=best_fit_par,
                            dpi=cornerplot_opt["dpi"],
                            pad=cornerplot_opt["pad"],
                            labelpad=cornerplot_opt["labelpad"],
                            truth_color=cornerplot_opt["truth_color"],
                            title_kwargs={"fontsize": 12},
                            scale_hist=cornerplot_opt["scale_hist"],
                            no_fill_contours=cornerplot_opt["no_fill_contours"],
                            fill_contours=cornerplot_opt["fill_contours"],
                            plot_datapoints=cornerplot_opt["plot_datapoints"],
                       #     data_kwargs={'zorder':10,'color':cornerplot_opt["stab_color"]},
                            contour_kwargs={'color':cornerplot_opt["stab_color"]},
                            hist_kwargs={'color':cornerplot_opt["stab_color"]},
                            kwargs={'zorder':10,'color':cornerplot_opt["stab_color"]}
                            )
        else:
            print(f"Stable overlay omitted: {N_samp} stable rows for {len(samp_best_fit_par)} plotted parameters; plotting the full finite posterior.")
            fig = corner.corner(samples_,
                                #range = range_stab,
                                bins=cornerplot_opt["bins"],
                                color=cornerplot_opt["color"],
                                reverse=cornerplot_opt["reverse"],
                                upper=cornerplot_opt["upper"],
                                labels=samp_labels,
                                quantiles=quantiles,
                                title_quantiles=title_quantiles,
                                levels=(0.6827, 0.9545,0.9973),
                                smooth=cornerplot_opt["smooth"],
                                smooth1d=cornerplot_opt["smooth1d"],
                                plot_contours=cornerplot_opt["plot_contours"],
                                show_titles=cornerplot_opt["show_titles"],
                                truths=best_fit_par,
                                dpi=cornerplot_opt["dpi"],
                                pad=cornerplot_opt["pad"],
                                labelpad=cornerplot_opt["labelpad"],
                                truth_color=cornerplot_opt["truth_color"],
                                title_kwargs={"fontsize": 12},
                                scale_hist=cornerplot_opt["scale_hist"],
                                no_fill_contours=cornerplot_opt["no_fill_contours"],
                                fill_contours=cornerplot_opt["fill_contours"],
                                plot_datapoints=cornerplot_opt["plot_datapoints"],
                                kwargs=kwargs)


        corner.corner(samples_,
                        range = range_stab,
                        bins=cornerplot_opt["bins"],
                        color=cornerplot_opt["color"],
                        reverse=cornerplot_opt["reverse"],
                        upper=cornerplot_opt["upper"],
                        labels=samp_labels,
                        quantiles=quantiles,
                        title_quantiles=title_quantiles,
                        levels=(0.6827, 0.9545,0.9973),
                        smooth=cornerplot_opt["smooth"],
                        smooth1d=cornerplot_opt["smooth1d"],
                        plot_contours=cornerplot_opt["plot_contours"],
                        show_titles=cornerplot_opt["show_titles"],
                        truths=best_fit_par,
                        dpi=cornerplot_opt["dpi"],
                        pad=cornerplot_opt["pad"],
                        labelpad=cornerplot_opt["labelpad"],
                        truth_color=cornerplot_opt["truth_color"],
                        title_kwargs={"fontsize": 12},
                        scale_hist=cornerplot_opt["scale_hist"],
                        no_fill_contours=cornerplot_opt["no_fill_contours"],
                        fill_contours=cornerplot_opt["fill_contours"],
                        plot_datapoints=cornerplot_opt["plot_datapoints"],
                       # data_kwargs={'zorder':-1,'color':cornerplot_opt["color"]},
                        contour_kwargs={'color':cornerplot_opt["color"]},
                        hist_kwargs={'color':cornerplot_opt["color"]},
                        kwargs={'zorder':-1,'color':cornerplot_opt["color"]},
                        fig =fig)

    else:

        fig = corner.corner(samples_,
                            #range = range_stab,
                            bins=cornerplot_opt["bins"],
                            color=cornerplot_opt["color"],
                            reverse=cornerplot_opt["reverse"],
                            upper=cornerplot_opt["upper"],
                            labels=samp_labels,
                            quantiles=quantiles,
                            title_quantiles=title_quantiles,
                            levels=(0.6827, 0.9545,0.9973),
                            smooth=cornerplot_opt["smooth"],
                            smooth1d=cornerplot_opt["smooth1d"],
                            plot_contours=cornerplot_opt["plot_contours"],
                            show_titles=cornerplot_opt["show_titles"],
                            truths=best_fit_par,
                            dpi=cornerplot_opt["dpi"],
                            pad=cornerplot_opt["pad"],
                            labelpad=cornerplot_opt["labelpad"],
                            truth_color=cornerplot_opt["truth_color"],
                            title_kwargs={"fontsize": 12},
                            scale_hist=cornerplot_opt["scale_hist"],
                            no_fill_contours=cornerplot_opt["no_fill_contours"],
                            fill_contours=cornerplot_opt["fill_contours"],
                            plot_datapoints=cornerplot_opt["plot_datapoints"],
                            kwargs=kwargs)

    for ax in fig.get_axes():
          #ax.tick_params(axis='both', which='major', labelsize=14)
          #ax.tick_params(axis='both', which='minor', labelsize=12)
          ax.tick_params(axis='both', labelsize=16)
          ax.xaxis.label.set_size(16)
          ax.yaxis.label.set_size(16)\

    if type_plot == 'mcmc':
        fig.savefig(obj.mcmc_corner_plot_file, bbox_inches='tight')
    if type_plot == 'nest':
        fig.savefig(obj.nest_corner_plot_file, bbox_inches='tight')

    ### to avoid memory leak in loops!
    fig.clf()
    del fig
    del samples_
    del samp
    del samples
    del ln
    print("Cornerplot done!")

    return





def planet_orbit_xyz(obj, planet):

    u1 = obj.params.stellar_mass * (4*np.pi*np.pi)/(365.25*365.25)
    mean_orb = np.linspace(0,2.0*np.pi, 360)

    x = np.zeros(len(mean_orb))
    y = np.zeros(len(mean_orb))
    z = np.zeros(len(mean_orb))
    u = np.zeros(len(mean_orb))
    v = np.zeros(len(mean_orb))
    w = np.zeros(len(mean_orb))

    dist =  np.zeros(len(mean_orb))

    q = (1.0 - obj.e[int(planet)])*float(obj.fit_results.a[int(planet)])


    #this need to be fixed to work with arrays
    for f in range(len(mean_orb)):
        x[f],y[f],z[f],u[f],v[f],w[f] = jac2astrocen.mco_el2x(u1,q, obj.e[int(planet)],
                                                       np.radians(obj.i[int(planet)]-90.0),
                                                       np.radians(obj.w[int(planet)]) - np.radians(obj.Node[int(planet)]),
                                                       np.radians(obj.Node[int(planet)]), mean_orb[f])

        dist[f] =  np.sqrt(x[f]**2.0 + y[f]**2.0 + z[f]**2.0)

    x_p,y_p,z_p,u_p,v_p,w_p = jac2astrocen.mco_el2x(u1,q, obj.e[int(planet)],
                                                       np.radians(obj.i[int(planet)]-90.0),
                                                       np.radians(obj.w[int(planet)]) - np.radians(obj.Node[int(planet)]),
                                                       np.radians(obj.Node[int(planet)]), np.radians(obj.M0[int(planet)]))


    min_index = np.unravel_index(np.argmin(dist, axis=None), dist.shape)
    max_index = np.unravel_index(np.argmax(dist, axis=None), dist.shape)

    return np.array([x,y,z,u,v,w]), np.array([x_p,y_p,z_p,u_p,v_p,w_p]), np.array([x[min_index],y[min_index],z[min_index],u[min_index],v[min_index],w[min_index]]), np.array([x[max_index],y[max_index],z[max_index],u[max_index],v[max_index],w[max_index]])



def get_xyz(obj):

    st_mass = obj.params.stellar_mass * (4*np.pi*np.pi)/(365.25*365.25)
    frho3 = 1.0
    u1 = st_mass
    obj.xyz_mass[0] = u1

    Msun = 1.989e33
    Au = 1.49597892e13

    ##### this is a hack avoiding transit init crash!!!! TB fixed/removed
    if len(np.atleast_1d(obj.fit_results.mass)) == 0:
        return obj
    #####################################################################


    for i in range(9):
        if not bool(obj.use_planet[i]):
            continue

        pl_mass_in_st = float(obj.masses[i])/ 1047.5654817267318

        pl_mass = pl_mass_in_st * (4*np.pi*np.pi)/(365.25*365.25)
        q = (1.0 - obj.e[i])*float(obj.semimajor[i])

        obj.rpl[i+1] = frho3*(1.5*pl_mass_in_st*Msun/(2*np.pi))**0.3333333333/Au
#             rpl(i) =  frho3*(1.5d0*mpl0*MSUN/TWOPI)**0.3333333333d0/AU

        obj.rhill[i+1] = float(obj.semimajor[i])*(pl_mass/(3.0*st_mass))**0.3333333333

        u1 = u1 +pl_mass
        obj.xyz_mass[i+1] = pl_mass

        x_p,y_p,z_p,u_p,v_p,w_p = jac2astrocen.mco_el2x(u1,q,
                                                       obj.e[i],
                                                       np.radians(obj.i[i]),
                                                       np.radians(obj.w[i]) - np.radians(obj.Node[i]),
                                                       np.radians(obj.Node[i]),  np.radians(obj.M0[i]))

        obj.xzy[i+1] =     np.array([x_p,y_p,z_p])
        obj.uvw[i+1] =     np.array([u_p,v_p,w_p])


    return obj

def get_Hill_satb(obj):

    st_mass = float(obj.params.stellar_mass)* 1047.5654817267318

    if len(np.atleast_1d(obj.fit_results.mass)) <=1:
        return False
    #####################################################################

    else:

        Delta_a = (float(obj.fit_results.a[1]) - float(obj.fit_results.a[0]))/float(obj.fit_results.a[0])
        Mu = 2.4*( (float(obj.fit_results.mass[0])/ st_mass) + (float(obj.fit_results.mass[1])/ st_mass) )**(1.0/3.0)

        if Mu >= Delta_a:
            return False
        else:
            return True


def get_AMD_stab(obj):

    st_mass = float(obj.params.stellar_mass)* 1047.5654817267318

    AMD_stable = True

    if len(np.atleast_1d(obj.fit_results.mass)) <=1:
        return False
    else:
        pl_ecc    = np.array([float(obj.e[i]) for i in range(obj.npl)])
        pl_a      = np.array([float(obj.fit_results.a[i]) for i in range(obj.npl)])
        pl_mass   = np.array([float(obj.fit_results.mass[i]) for i in range(obj.npl)])
        inds = pl_a.argsort()
        sorted_pl_ecc = pl_ecc[inds]
        sorted_pl_a = pl_a[inds]
        sorted_pl_mass = pl_mass[inds]

        for i in range(obj.npl - 1):

            alpha    = sorted_pl_a[i]/sorted_pl_a[i+1]
            gamma    = sorted_pl_mass[i]/sorted_pl_mass[i+1]
            epsilon  = (sorted_pl_mass[i]+sorted_pl_mass[i+1])/st_mass

            AMD = gamma*np.sqrt(alpha)*(1.-np.sqrt(1.-sorted_pl_ecc[i]**2)) + 1.-np.sqrt(1.-sorted_pl_ecc[i+1]**2)

            AMD_Hill = gamma*np.sqrt(alpha) + 1. - (1.+gamma)**1.5 * np.sqrt(alpha/(gamma+alpha) * (1.+(3.**(4./3.)*epsilon**(2./3.)*gamma)/((1.+gamma)**2)))
            #print( AMD,AMD_Hill)

            if AMD >= AMD_Hill:
                return False

    return AMD_stable

def loglik_AMD_penalty(pl_a,pl_ecc,pl_mass,st_mass):

    for i in range(len(pl_a) - 1):

        alpha    = pl_a[i]/pl_a[i+1]
        gamma    = pl_mass[i]/pl_mass[i+1]
        epsilon  = (pl_mass[i]+pl_mass[i+1])/st_mass

        AMD = gamma*np.sqrt(alpha)*(1.-np.sqrt(1.- pl_ecc[i]**2)) + 1.-np.sqrt(1.- pl_ecc[i+1]**2)

        AMD_Hill = gamma*np.sqrt(alpha) + 1. - (1.+gamma)**1.5 * np.sqrt(alpha/(gamma+alpha) * (1.+(3.**(4./3.)*epsilon**(2./3.)*gamma)/((1.+gamma)**2)))
        if AMD >= AMD_Hill:
            return -np.inf

    return 0



def randomString(stringLength=5):
    """
    Generate a random string of fixed length
    """
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))

def copy_file_to_datafiles(path):
    '''
    creates a temp_ velocity file in the root directory of the GUI.

    input: full path to the file
    output: temp_name of the file to be loaded
    '''

    dirname, basename = os.path.split(path)
    #temp_dir = './datafiles'#tempfile.gettempdir()

    tmp = '/tmp/es2'


    if platform.system() == 'Darwin':
        if not os.path.exists(tmp):
            os.system("mkdir %s"%tmp)
        else:
            if os.access(tmp, os.W_OK):
                tmp = '/tmp/es2'
            else:
                tmp = '/tmp/%s'%(randomString(3))

        tmp = '%s/%s'%(tmp,randomString(5))
        os.system("mkdir %s"%tmp)
    else:
        tmp = tempfile.mkdtemp()

    temp_path = os.path.join(tmp, basename)

   # os.system("cp %s %s"%(path, temp_path))

    f_in = open(path, "r")
    lines = f_in.readlines()
    f  = open(temp_path, 'wb') # open the file

    for j in range(len(lines)):

        if lines[j].isspace() or lines[j].split()[0].startswith("#"):
            continue

        line = lines[j].split()
        text = b"%s  %s  %s \n"%(bytes(str(lines[j].split()[0]).encode()),bytes(str(lines[j].split()[1]).encode()),bytes(str(lines[j].split()[2]).encode()) )
        f.write(text)
    f.close()
    f_in.close()




    return temp_path

def mut_incl(i1,i2,capOm):
    '''
    Calculates the mutual inclination of two planets

    input parameters:

    i1,i2, Delta Omega: inclinations and diffence of the line of nodes in deg.

    output parameters:

    Delta i: mutual orbital inclination in deg.
    '''
    fb = np.degrees(np.arccos(((np.cos(np.radians(i1))*np.cos(np.radians(i2)))+
    (np.sin(np.radians(i1))*np.sin(np.radians(i2))*np.cos(np.radians(capOm))))))
    return fb


def add_jitter(obj, errors, ind):
    errors_with_jitt = np.array([np.sqrt(errors[i]**2 + obj.jitt[ii]**2)  for i,ii in enumerate(ind)])
    return errors_with_jitt





def get_stellar_rotation(obj, print_output=False):
    '''
    '''
    vsini   = float(obj.stellar_vsini)
    vsini_d = float(obj.stellar_vsini_err)

    R   = float(obj.stellar_radius)
    R_d = float(obj.stellar_radius_err)

    Rot = (2*np.pi*R*695700.0)/ (vsini * 86400.0)

    Delta_Rot = np.sqrt( ( ( (2*np.pi*R*695700.0)**2 * (vsini_d*86400.0)**2) + ((2*np.pi*R_d*695700.0)**2 * (vsini*86400.0)**2) ) /
(vsini*86400.0)**4
)

    if print_output == True:
        print("Stellar Rot. period = %s  +/- %s [days]"%(Rot, Delta_Rot))

    return [Rot,  Delta_Rot]


def get_rv_scatter(obj, print_output=False,use_kb2011=False):
    '''
    '''
    Solar_fact = 0.234
    Solar_fact_d = 0.014

    M   = float(obj.stellar_mass)
    M_d = float(obj.stellar_mass_err)

    L   = float(obj.stellar_luminosity)
    L_d = float(obj.stellar_luminosity_err)

    Teff   = float(obj.stellar_Teff)/5771.0
    Teff_d = float(obj.stellar_Teff_err)/5771.0

    if use_kb2011==True:

        A = (L / ((M**1.5)*(Teff**4.25))) * Solar_fact

        delta_A = 0.25*np.sqrt(
        (L**2.0 * (
        4.0*Teff**2.0 *(
        (4.0 * (M**2.0) * (Solar_fact_d**2.0)) +
        (9.0 * M_d**2.0 * Solar_fact**2.0)) +
        (289.0 * (Teff_d**2.0) * (M**2.0) * (Solar_fact**2.0) )) +
        (16.0 * (L_d**2.0) * (Teff**2.0) * (M**2.0) * (Solar_fact**2.0)) ) / ((Teff**(21.0/2.0)) * (M**5.0)) )

        if print_output == True:
            print("KB2011 jitter = %s +/- %s [m/s]"%(A, delta_A))

    else:
        A = (L/M) * Solar_fact

        delta_A = np.sqrt( (L**2.0 * ((M**2.0)*(Solar_fact_d**2.0) + (M_d**2.0)*(Solar_fact**2.0)  ) +
        ((L_d**2.0) *(M**2.0) * (Solar_fact**2.0) ) )/ M**4.0 )


        if print_output == True:
            print("KB1995 jitter = %s +/- %s [m/s]"%(A, delta_A))


    return [A, delta_A]



def export_transit_data(obj, file='transit_data.dat', delimiter=' ',  print_data=False,  header=True, width=10, precision=6):
    import numpy as np

    if len(obj.transit_results[1]) == 0:
        return

    transit_results_sep = obj.transit_results[1]

    all_BJD = []
    all_flux = []
    all_flux_err = []
    all_flux_model = []
    all_o_c = []
    all_idset = []
    all_flux_gp = []

    for z in range(len(transit_results_sep[0])):
        all_BJD.extend(np.array(transit_results_sep[0][z], dtype=float))
        all_flux.extend(np.array(transit_results_sep[1][z], dtype=float))
        all_flux_err.extend(np.array(transit_results_sep[2][z], dtype=float))
        all_flux_model.extend(np.array(transit_results_sep[3][z], dtype=float))
        all_o_c.extend(np.array(transit_results_sep[4][z], dtype=float))
        all_idset.extend(np.array([z] * len(transit_results_sep[0][z]), dtype=int))
        all_flux_gp.extend(np.array(transit_results_sep[6][z], dtype=float))

    # Stack and sort by BJD
    data = np.vstack((all_BJD, all_flux, all_flux_err, all_flux_model, all_o_c, all_idset, all_flux_gp)).T
    data = data[data[:, 0].argsort()]

    # Prepare header
    if header:
        head_lines = ["# Dataset mapping:"]
        unique_z = sorted(set(all_idset))
        for unique_id in unique_z:
            head_lines.append(f"# {unique_id}: {obj.tra_data_sets[unique_id][-1]}")

        head_lines.append("#" + delimiter.join([
            f"{'BJD':<{width}}",
            f"{'flux':<{width}}",
            f"{'flux_sigma':<{width}}",
            f"{'flux_model':<{width}}",
            f"{'o_c':<{width}}",
            f"{'idset':<{width}}",
            f"{'flux_gp':<{width}}"
        ]))
        head = "\n".join(head_lines)
    else:
        head = ""

    # Save to file
    fmt = []
    for i in range(data.shape[1]):
        if i == 5:
            fmt.append(f'%{width}d')
        else:
            fmt.append(f'%{width}.{precision}f')
    fmt = tuple(fmt)

    np.savetxt(file, data, delimiter=delimiter, fmt=fmt, header=head)

    if print_data:
        for row in data:
            formatted_row = []
            for i, val in enumerate(row):
                if i == 5:
                    formatted_row.append(f"{int(val):{width}d}")
                else:
                    formatted_row.append(f"{val:.{precision}f}")
            print(' '.join(formatted_row))

    return data


def export_transit_model(obj, file='transit_model.dat', delimiter=' ',  print_data=False,  header=True, width=10, precision=6):


    if len(obj.transit_results[3]) == 0:
        return

    transit_results_rich = obj.transit_results[3]

 
    all_BJD = np.array(transit_results_rich[0], dtype=float)
    all_flux = np.array(transit_results_rich[1], dtype=float)

    # Stack and sort by BJD
    data = np.vstack((all_BJD, all_flux)).T
    data = data[data[:, 0].argsort()]

    # Prepare header
    if header:
        head_lines = [""]
        head_lines.append("#" + delimiter.join([
            f"{'BJD':<{width}}",
            f"{'flux':<{width}}"
        ]))
        head = "\n".join(head_lines)
    else:
        head = ""

    # Save to file
    fmt = []
    for i in range(data.shape[1]):
        if i == 5:
            fmt.append(f'%{width}d')
        else:
            fmt.append(f'%{width}.{precision}f')
    fmt = tuple(fmt)

    np.savetxt(file, data, delimiter=delimiter, fmt=fmt, header=head)

    if print_data:
        for row in data:
            formatted_row = []
            for i, val in enumerate(row):
                if i == 5:
                    formatted_row.append(f"{int(val):{width}d}")
                else:
                    formatted_row.append(f"{val:.{precision}f}")
            print(' '.join(formatted_row))

    return data
 






def export_RV_data(obj, file='RV_data.txt', delimiter=' ',  print_data=False,  header = True, width = 10, precision = 3):

    if len(obj.fit_results.model_data[5])==0:
        return
    
    if header == True:
        head = "{}       ,   {},   {},    {},    {},    {},  ".format("BJD", "rvs","rvs_sigma","idset","o_c","model_rvs")
        for i in range(obj.npl):
            head = head + "model pl_{},    ".format(i+1) 
        if obj.doGP == True:
            head = head + "{},   {},   {},   ".format("GP model mu", "GP model var","GP model std")   

    else:
        head = " "
    
    if obj.doGP == True:
        RV_data = np.concatenate((obj.fit_results.model_data, np.array(obj.gp_model_data)))
    else:
        RV_data = obj.fit_results.model_data 
    
    np.savetxt(file, RV_data.T, delimiter=delimiter, fmt='%'+'%s.%sf'%(int(width), int(precision)), header=head)
    
 
    
    return
    


def export_RV_data_old(obj, idset_ts, file="RV_data.txt",  jitter=False, o_c=False,
                   print_data=False, remove_offset = False, width = 10, precision = 3):

    if len(obj.filelist.idset)==0:
        return

    output_file = str(file)
    f = open(output_file, 'w')

    idset_ts = np.array(np.atleast_1d(idset_ts)) #-1


    JD = obj.fit_results.rv_model.jd
    if not o_c:
        rv = obj.fit_results.rv_model.rvs
    else:
        rv = obj.fit_results.rv_model.o_c
    id_set = obj.filelist.idset

    if jitter==True:
        sigma = add_jitter(obj,obj.fit_results.rv_model.rv_err, id_set)
    else:
        sigma =  obj.fit_results.rv_model.rv_err


    for i in range(len(idset_ts)):
        for ii in range(len(JD)):
             if int(id_set[ii]) != int(idset_ts[i]):
                 continue
             else:
                 if not remove_offset:
                     rv[ii] = rv[ii] - float(obj.rvoff[i])

                 f.write('{0:{width}.{precision}f}  {1:{width}.{precision}f}  {2:{width}.{precision}f}  {3:{width}.{precision}f}  \n'.format(float(JD[ii]), float(rv[ii]), float(sigma[ii]), idset_ts[i], width = width, precision = precision )   )

                 if print_data:
                     print('{0:{width}.{precision}f}  {1:{width}.{precision}f}  {2:{width}.{precision}f}  {3:{width}.{precision}f}'.format(float(JD[ii]), float(rv[ii]), float(sigma[ii]), idset_ts[i], width = width, precision = precision )   )

    f.close()
    print('Done!')
    return
    
    
    

def export_RV_model(obj, file='RV_model.txt', delimiter=' ',  print_data=False,  header = True, width = 10, precision = 3):

    if len(obj.fit_results.model_data[5])==0:
        return
    
    if header == True:
        head = "{}       ,   {},    ".format("BJD", "RV model")
        for i in range(obj.npl):
            head = head + "model pl_{},    ".format(i+1) 
        if obj.doGP == True:
            head = head + "{},   {},   {},   ".format("GP model mu", "GP model var","GP model std")   

    else:
        head = " "
    
    if obj.doGP == True:
        RV_data = np.concatenate((obj.fit_results.fit.T, np.array(obj.gp_model_curve)))
    else:
        RV_data = obj.fit_results.fit.T
    
    np.savetxt(file, RV_data.T, delimiter=delimiter, fmt='%'+'%s.%sf'%(int(width), int(precision)), header=head)
    
    
    return    
    
    

def export_RV_model_old(obj, file="RV_model.txt", width = 10, precision = 4,print_data=False):

    if len(obj.fit_results.rv_model.jd)==0:
        return

    #if not os.path.exists(path):
   #     os.makedirs(path)

    output_file = str(file)
    f = open(output_file, 'w')


    JD = obj.fit_results.model_jd

    if obj.doGP == True:
        y_model = obj.fit_results.model + obj.gp_model_curve[0]
    else:
        y_model = obj.fit_results.model

    for i in range(len(JD)):
        f.write('{0:{width}.{precision}f}  {1:{width}.{precision}f} \n'.format(float(JD[i]), float(y_model[i]), width = width, precision = precision) )
    f.close()

    if print_data:
        print('{0:{width}.{precision}f}  {1:{width}.{precision}f} \n'.format(float(JD[i]), float(y_model[i]), width = width, precision = precision) )


    print('Done!')
    return


def export_orbital_evol(obj, file="orb_evol.txt", planet = 1, width = 10, precision = 6):

    k = int(planet-1)


    if len(obj.evol_T[k])==0 or k < 0:
        print("No N-body integrations done?")
        return

    output_file = str(file)
    f = open(output_file, 'w')

    #obj.evol_T_energy
    #obj.evol_energy

    #obj.evol_momentum['lx']
    #obj.evol_momentum['ly']
    #obj.evol_momentum['lz']



    T = obj.evol_T[k]
    a = obj.evol_a[k]
    e = obj.evol_e[k]
    om = obj.evol_p[k]
    M = obj.evol_M[k]

    inc = obj.evol_i[k]
    Om =obj.evol_Om[k]


    for i in range(len(T)):
        f.write('{0:{width}.{precision}f}  {1:{width}.{precision}f}  {2:{width}.{precision}f}  {3:{width}.{precision}f}  {4:{width}.{precision}f}  {5:{width}.{precision}f}  {6:{width}.{precision}f}  \n'.format(float(T[i]), float(a[i]), float(e[i]), float(om[i]),float(M[i]),float(inc[i]),float(Om[i]), width = width, precision = precision) )
    f.close()
    print('Done!')
    return

def check_temp_RV_file(obj):

    for i in range(obj.filelist.ndset):

        if os.path.exists(obj.filelist.files[i].path) and os.access(obj.filelist.files[i].path, os.W_OK):
           # print(obj.filelist.files[i].path)
            continue


        else:
            dirname, basename = os.path.split(obj.filelist.files[i].path)

            if platform.system() == 'Darwin':
                dirname = '/tmp/es_%s'%(randomString(3))
            else:
                dirname = '/tmp/%s'%(randomString(3))


            dirname = '%s/%s'%(dirname,randomString(5))

            #else:
            #    dirname = os.path.abspath(dirname)+'_'+randomString(3)+'/'+randomString(3)
            obj.filelist.files[i].path = os.path.join(dirname, basename)


            os.makedirs(dirname)
            f  = open(obj.filelist.files[i].path, 'wb') # open the file
            for j in range(len(obj.rv_data_sets[i][0])):
                if str(obj.rv_data_sets[i][0][j]).startswith("#"):
                    continue
                text = b"%s  %s  %s \n"%(bytes(str(obj.rv_data_sets[i][0][j]).encode()),bytes(str(obj.rv_data_sets[i][1][j]).encode()),bytes(str(obj.rv_data_sets[i][2][j]).encode()) )
                f.write(text)
            f.close()



def modify_temp_RV_file_old(obj, file_n = 0, add_error = 0, data_to_keep = None):

    if obj.filelist.ndset < file_n +1:
        print("No RV file # %s"%(file_n+1))
        return
    elif not os.path.exists(obj.filelist.files[file_n].path):
        return
    else:
        if add_error < 0:
            sign = -1
        else:
            sign = 1
        new_error = []
        for j in range(len(obj.rv_data_sets[file_n][0])):
            k = obj.rv_data_sets[file_n][2][j]**2 + add_error**2 *sign
            if k < 0:
                print("You seem to subtract %s from the error budget. As a result, the RV uncertainty of one or more elements would be negative. Errors cannot be negative. Please subtract another value"%add_error)
                return
            new_error.append(k)
        f  = open(obj.filelist.files[file_n].path, 'wb') # open the file

        for j in range(len(obj.rv_data_sets[file_n][0])):
            if str(obj.rv_data_sets[file_n][0][j]).startswith("#") or  data_to_keep != None and  j not in data_to_keep:
                continue
            text = b"%s  %s  %s \n"%(bytes(str(obj.rv_data_sets[file_n][0][j]).encode()),bytes(str(obj.rv_data_sets[file_n][1][j]).encode()),bytes(str(np.sqrt(new_error[j])).encode()) )
            f.write(text)
        f.close()


        obj.filelist.read_rvfiles(obj.rvoff)

        return obj


def common_member(a, b):
    a_set = set(a)
    b_set = set(b)

    if (a_set & b_set):
        return True
    else:
        return False



def bin_data(JD,rv,sigma,idset, bin_size = 1.0):


    binning_list = []
    v_l = [1100000,1212234324]
    for i in range(len(JD)):
        b_l = [x for x,z in enumerate(JD) if abs(z - JD[i]) < bin_size]
        if common_member(b_l, v_l):
            continue
        else:
            binning_list.append(b_l)
            v_l = b_l

    mj_all = []
    mr_all = []
    ms_all = []
    mi_all = []

    for x in range(len(binning_list)):

            mj_all.append(np.mean(JD[binning_list[x]]))
            mr_all.append(np.average(rv[binning_list[x]], weights=1./sigma[binning_list[x]]))
            #ms_all.append(np.average(ms/np.sqrt(len(ms)), weights=1./ms) )
            ms_all.append(np.average(sigma[binning_list[x]]) )
            #ms_all.append( np.sqrt( (np.average(ms/np.sqrt(len(ms)), weights=1./ms)**2.0) + (abs(max(mr)-min(mr))**2.0) ) )
            #ms_all.append( np.sqrt( (np.average(ms/np.sqrt(len(ms)), weights=1./ms)**2.0) + np.std(mr)**2.0) ) )
            #print np.median(mr), np.std(mr)
            mi_all.append(np.mean(idset[binning_list[x]]))

    JD2, indices = np.unique(np.asarray(mj_all), return_index=True)

    ind    = np.array(indices)
    mr_all = np.array(mr_all)
    mj_all = np.array(mj_all)
    ms_all = np.array(ms_all)
    mi_all = np.array(mi_all)

    mr_all = mr_all[ind]
    mj_all = mj_all[ind]
    ms_all = ms_all[ind]
    mi_all = mi_all[ind]

    print(len(JD2))

    return  mj_all, mr_all, ms_all, mi_all



def bin_rv_data(obj, file_n = 0, bin_size = 1.0, bin_tf = False):

    if bin_tf == False:
        obj.rv_data_sets[file_n] = dill.copy(obj.rv_data_sets_init[file_n])
        return

    else:

        JD    = np.array(obj.rv_data_sets[file_n][0])
        rv    = np.array(obj.rv_data_sets[file_n][1])
        sigma = np.array(obj.rv_data_sets[file_n][2])
        idset = np.array(obj.rv_data_sets[file_n][3])


        mj_all,mr_all,ms_all,mi_all = bin_data(JD,rv,sigma,idset, bin_size = bin_size)

        obj.rv_data_sets[file_n] = np.array([mj_all,mr_all,ms_all,mi_all])

        #obj.rv_data_sets[file_n][0] = dill.copy(mj_all)
        #obj.rv_data_sets[file_n][1] = dill.copy(mr_all)
        #obj.rv_data_sets[file_n][2] = dill.copy(ms_all)
        #obj.rv_data_sets[file_n][3] = dill.copy(mi_all)

        return obj



def bin_rv_dataOld(obj, file_n = 0, bin_size = 1.0, bin_tf = False):

    if bin_tf == False:
        obj.rv_data_sets[file_n] = dill.copy(obj.rv_data_sets_init[file_n])
        return

    else:

        JD    = np.array(obj.rv_data_sets[file_n][0])
        rv    = np.array(obj.rv_data_sets[file_n][1])
        sigma = np.array(obj.rv_data_sets[file_n][2])
        idset = np.array(obj.rv_data_sets[file_n][3])


        mask = np.zeros(len(JD))
        mj_all = []
        mr_all = []
        ms_all = []
        mi_all = []


        for x in range(len(JD)):

            JD_int = JD.astype(int)

            mask = (JD_int != JD_int[x]).astype(int)

            mj = np.ma.masked_array(JD, mask=mask).compressed()
            mr = np.ma.masked_array(rv, mask=mask).compressed()
            ms = np.ma.masked_array(sigma, mask=mask).compressed()
            mi = np.ma.masked_array(idset, mask=mask).compressed()


            mj_all.append(np.mean(mj))
            mr_all.append(np.average(mr, weights=1./ms))
            #ms_all.append(np.average(ms/np.sqrt(len(ms)), weights=1./ms) )
            ms_all.append(np.average(ms) )
            mi_all.append(np.mean(mi))

            #ms_all.append( np.sqrt( (np.average(ms/np.sqrt(len(ms)), weights=1./ms)**2.0) + (abs(max(mr)-min(mr))**2.0) ) )
            #ms_all.append( np.sqrt( (np.average(ms/np.sqrt(len(ms)), weights=1./ms)**2.0) + np.std(mr)**2.0) ) )
            #print np.median(mr), np.std(mr)

        JD, indices = np.unique(np.asarray(mj_all), return_index=True)

        ind    = np.array(indices)
        mr_all = np.array(mr_all)
        mj_all = np.array(mj_all)
        ms_all = np.array(ms_all)
        mi_all = np.array(mi_all)

        mr_all = mr_all[ind]
        mj_all = mj_all[ind]
        ms_all = ms_all[ind]
        mi_all = mi_all[ind]


        obj.rv_data_sets[file_n] = np.array([mj_all,mr_all,ms_all,mi_all])

        #obj.rv_data_sets[file_n][0] = dill.copy(mj_all)
        #obj.rv_data_sets[file_n][1] = dill.copy(mr_all)
        #obj.rv_data_sets[file_n][2] = dill.copy(ms_all)
        #obj.rv_data_sets[file_n][3] = dill.copy(mi_all)

        return obj



def inject_signal(obj, P = 100, M0 = 0, K = 10):

    obj2 = dill.copy(obj)


    obj2.add_planet(K,P,0.0,0.0,M0,90.0,0.0)
    obj2.fitting(minimize_fortran=True, minimize_loglik=True, outputfiles=[0,1,0], amoeba_starts=0)

    new_RVs = obj.fit_results.rv_model.rvs + obj2.fit_results.model_rvs

    obj2.fit_results.rv_model.rvs = dill.copy(new_RVs)
#    del obj2

    #return [obj.fit_results.rv_model.jd, new_RVs, obj.fit_results.rv_model.rv_err]

    return obj2

def identify_power_peaks(x,y,sig_level=np.array([]), power_level=np.array([]) ):

    per_ind = argrelextrema(y, np.greater)
    per_x   = x[per_ind]
    per_y   = y[per_ind]

    peaks_sort = sorted(range(len(per_y)), key=lambda k: per_y[k], reverse=True)

    per_x   = per_x[peaks_sort]
    per_y   = per_y[peaks_sort]

    peaks_pos = [per_x,per_y]

    ################## text generator #################
    text_peaks = """
"""
    if power_level.size != 0 and sig_level.size != 0:

        text_peaks = text_peaks +"""FAP levels
-----------------------------------
"""
        for ii in range(len(power_level)):
            text_peaks = text_peaks +"""
%.2f per cent = %.4f"""%(power_level[ii]*100.0,sig_level[ii])

    text_peaks = text_peaks + """
----------------------------------------------
The 10 strongest peaks
----------------------------------------------
"""

    if len(per_x)  < 10:
        max_peaks = len(per_x)
    else:
        max_peaks = 10

    for j in range(max_peaks):
        text_peaks = text_peaks +"""
period = %.2f [d], power = %.4f"""%(per_x[j],per_y[j])
        if sig_level.size != 0 and per_y[j] > sig_level[-1]:
            text_peaks = text_peaks +"""  significant"""

    return text_peaks , peaks_pos


def run_detection_rate(obj2,options):

    from tqdm import tqdm

    results = {"P":[],"K":[],"M0":[],'gls':[],'power':[],'peaks':[],'input':options}


    for P in tqdm(options['P']):
        for K in tqdm(options['K']):
            for M0 in tqdm(options['M0']):

                #obj2 = dill.copy(inject_signal(obj, K=K,P=P,M0=M0))

                sinewave = K * np.sin(2 * np.pi * 1/P * obj2.fit_results.rv_model.jd + np.radians(M0))



                if options['gls_o-c']:
                    new_RVs = obj2.fit_results.rv_model.o_c + sinewave
                else:
                    new_RVs = obj2.fit_results.rv_model.rvs + sinewave

                RV_per = gls.Gls((obj2.fit_results.rv_model.jd, new_RVs, obj2.fit_results.rv_model.rv_err), fast=True, verbose=False, norm="ZK",ofac=options['gls_ofac'], fbeg=1/options['gls_Pmax'], fend=1/options['gls_Pmin'])


                text_peaks, pos_peaks = identify_power_peaks(1/RV_per.freq, RV_per.power, power_level = options['power_levels'], sig_level = RV_per.powerLevel(np.array(options['power_levels'])) )

                results["P"].append(P)
                results["K"].append(K)
                results["M0"].append(M0)
                results["gls"].append(dill.copy(RV_per.powerLevel(np.array(options['power_levels']))))
                results["peaks"].append(dill.copy(pos_peaks[0][0:10]))
                results["power"].append(dill.copy(pos_peaks[1][0:10]))

    return results



def modify_temp_RV_file(obj, file_n = 0, add_error = 0, data_to_keep = None):

    if obj.filelist.ndset < file_n +1:
        print("No RV file # %s"%(file_n+1))
        return
    elif not os.path.exists(obj.filelist.files[file_n].path):
        return
    else:
        if add_error < 0:
            sign = -1
        else:
            sign = 1
        new_error = []
        for j in range(len(obj.rv_data_sets[file_n][0])):
            k = obj.rv_data_sets[file_n][2][j]**2 + add_error**2 *sign
            if k < 0:
                print("You seem to subtract %s from the error budget. As a result, the RV uncertainty of one or more elements would be negative. Errors cannot be negative. Please subtract another value"%add_error)
                return
            new_error.append(k)
        f  = open(obj.filelist.files[file_n].path, 'wb') # open the file

        org_data_file = obj.rv_data_sets[file_n]

        for j in range(len(org_data_file[0])):
            if str(org_data_file[0][j]).startswith("#") or  data_to_keep != None and  j not in data_to_keep:
                continue
            text = b"%s  %s  %s \n"%(bytes(str(org_data_file[0][j]).encode()),bytes(str(org_data_file[1][j]).encode()),bytes(str(np.sqrt(new_error[j])).encode()) )
            f.write(text)
        f.close()


        obj.filelist.read_rvfiles(obj.rvoff)

        return obj




### some experimets! ###
def sigma_clip(obj, type = 'RV', sigma_clip = 10, file_n = 0, add_error = 0, remove_mean = False, verbose = True):


    if type == 'RV':

        if sigma_clip == None:
            modify_temp_RV_file(obj, file_n = file_n, add_error = add_error, data_to_keep = None)
            return
        else:
            obj2 = dill.copy(obj)

            modify_temp_RV_file(obj2, file_n = file_n, add_error = add_error, data_to_keep = None)
            #obj2.epoch = obj.epoch
            obj2.fitting(outputfiles=[1,1,1], minimize_fortran=True, minimize_loglik=True,amoeba_starts=0)

            JD_data  = obj2.fit_results.rv_model.jd[obj2.filelist.idset==file_n]
            o_c_data = obj2.fit_results.rv_model.o_c[obj2.filelist.idset==file_n]
            data_ind = obj2.filelist.idset

            c, low, upp = pdf.sigmaclip(o_c_data, sigma_clip, sigma_clip)
            remaining_idx    = [x for x, z in enumerate(o_c_data) if z in c]
            removed_idx    = [x for x, z in enumerate(o_c_data) if z not in c]

            modify_temp_RV_file(obj, file_n = file_n, add_error = add_error, data_to_keep = remaining_idx)

            del obj2

            if verbose:
                print("\n %s clipped epochs:"%type)
                for z in JD_data[removed_idx]:
                    print(z)

            return obj


    if type == 'act':
        if len(obj.act_data_sets[file_n]) == 0:
            print("No act. file # %s"%(file_n))
            return

        org_epoch     = obj.act_data_sets_init[file_n][0]
        org_data      = obj.act_data_sets_init[file_n][1]
        org_data_sig  = obj.act_data_sets_init[file_n][2]

        org_data_mean = org_data - np.mean(org_data)

        if sigma_clip != None:

            c, low, upp = pdf.sigmaclip(org_data_mean, sigma_clip, sigma_clip)
            remaining_idx    = [x for x, z in enumerate(org_data_mean) if z in c]
            removed_idx      = [x for x, z in enumerate(org_data_mean) if z not in c]

            obj.act_data_sets[file_n][1] = np.take(obj.act_data_sets_init[file_n][1], remaining_idx)
            obj.act_data_sets[file_n][0] = np.take(obj.act_data_sets_init[file_n][0], remaining_idx)
            obj.act_data_sets[file_n][2] = np.take(obj.act_data_sets_init[file_n][2], remaining_idx)

            new_org_data      = obj.act_data_sets[file_n][1]
            new_org_data_mean = new_org_data - np.mean(new_org_data)

            if verbose:
                print("\n %s clipped epochs:"%type)
                for z in org_epoch[removed_idx]:
                    print(z)

            if remove_mean == True:
                obj.act_data_sets[file_n][1] = new_org_data_mean

        else:
            if remove_mean == True:
                obj.act_data_sets[file_n][0] = org_epoch
                obj.act_data_sets[file_n][1] = org_data_mean
                obj.act_data_sets[file_n][2] = org_data_sig

            else:
                obj.act_data_sets[file_n][0] = org_epoch
                obj.act_data_sets[file_n][1] = org_data
                obj.act_data_sets[file_n][2] = org_data_sig

        return obj


    if type == 'tra':
        if len(obj.tra_data_sets[file_n]) == 0:
            print("No transit file # %s"%(file_n))
            return

        org_epoch     = obj.tra_data_sets_init[file_n][0]
        org_data      = obj.tra_data_sets_init[file_n][1]
        org_data_sig  = obj.tra_data_sets_init[file_n][2]
        org_data_air  = obj.tra_data_sets_init[file_n][3]
        org_data_o_c  = obj.tra_data_sets_init[file_n][4]

        org_data_mean = org_data_o_c - np.mean(org_data_o_c)

        if sigma_clip != None:

            c, low, upp = pdf.sigmaclip(org_data_mean, sigma_clip, sigma_clip)
            remaining_idx    = [x for x, z in enumerate(org_data_mean) if z in c]
            removed_idx      = [x for x, z in enumerate(org_data_mean) if z not in c]

            obj.tra_data_sets[file_n][3] = np.take(obj.tra_data_sets_init[file_n][3], remaining_idx)
            obj.tra_data_sets[file_n][0] = np.take(obj.tra_data_sets_init[file_n][0], remaining_idx)
            obj.tra_data_sets[file_n][2] = np.take(obj.tra_data_sets_init[file_n][2], remaining_idx)
            obj.tra_data_sets[file_n][1] = np.take(obj.tra_data_sets_init[file_n][1], remaining_idx)
            obj.tra_data_sets[file_n][4] = np.take(obj.tra_data_sets_init[file_n][4], remaining_idx)

            new_org_data      = obj.tra_data_sets[file_n][1]
            new_org_data_mean = new_org_data - np.mean(new_org_data)

            if verbose:
                print("\n %s clipped epochs:"%type)
                for z in org_epoch[removed_idx]:
                    print(z)

            #if remove_mean == True:
           #     obj.tra_data_sets[file_n][1] = new_org_data_mean

        else:
           # if remove_mean == True:
          #      obj.tra_data_sets[file_n][0] = org_epoch
           #     obj.tra_data_sets[file_n][1] = org_data_mean
          #      obj.tra_data_sets[file_n][2] = org_data_sig
          #
          #  else:
            obj.tra_data_sets[file_n][0] = org_epoch
            obj.tra_data_sets[file_n][1] = org_data
            obj.tra_data_sets[file_n][2] = org_data_sig
            obj.tra_data_sets[file_n][3] = org_data_air
            obj.tra_data_sets[file_n][4] = org_data_o_c
        return obj





### some experimets! ###
def sigma_clip_act(JD_,act_,sigma_,idset_, sigma_clip = 10, verbose = True):



    act = act_ - np.mean(act_)


    c, low, upp = pdf.sigmaclip(act, sigma_clip, sigma_clip)
    remaining_idx    = [x for x, z in enumerate(act) if z in c]
    removed_idx      = [x for x, z in enumerate(act) if z not in c]

    act   = np.take(act_, remaining_idx)
    JD    = np.take(JD_, remaining_idx)
    sigma = np.take(sigma_, remaining_idx)
    idset = np.take(idset_, remaining_idx)

    if verbose:
        print("\n activity clipped epochs:")
        for z in JD_[removed_idx]:
            print(z)

    return JD,act,sigma,idset


def sigma_clip_rvs(JD_,rvs_,sigma_,idset_, sigma_clip = 10, verbose = True):



    rvs = rvs_ - np.mean(rvs_)


    c, low, upp = pdf.sigmaclip(rvs, sigma_clip, sigma_clip)
    remaining_idx    = [x for x, z in enumerate(rvs) if z in c]
    removed_idx      = [x for x, z in enumerate(rvs) if z not in c]

    rvs   = np.take(rvs_, remaining_idx)
    JD    = np.take(JD_, remaining_idx)
    sigma = np.take(sigma_, remaining_idx)
    idset = np.take(idset_, remaining_idx)

    if verbose:
        print("\n activity clipped epochs:")
        for z in JD_[removed_idx]:
            print(z)

    return JD,rvs,sigma,idset

def transit_data_norm(obj,  file_n = 0, norm = False, verbose = True):


    if len(obj.tra_data_sets[file_n]) == 0:
        print("No transit file # %s"%(file_n))
        return

    if norm == True:
        obj.tra_data_sets[file_n][0] = obj.tra_data_sets_init[file_n][0]
        obj.tra_data_sets[file_n][1] = obj.tra_data_sets_init[file_n][1]/np.mean(obj.tra_data_sets_init[file_n][1])
        obj.tra_data_sets[file_n][2] = obj.tra_data_sets_init[file_n][2]/np.mean(obj.tra_data_sets_init[file_n][1])
    else:
        obj.tra_data_sets[file_n][0] = dill.copy(obj.tra_data_sets_init[file_n][0])
        obj.tra_data_sets[file_n][1] = dill.copy(obj.tra_data_sets_init[file_n][1])
        obj.tra_data_sets[file_n][2] = dill.copy(obj.tra_data_sets_init[file_n][2])

    return obj



### some experimets! ###
def gen_RV_curve(obj,x=None):
    obj2 = dill.copy(obj)

    if len(x) > 3:
        f  = open('datafiles/RV_curve', 'wb') # open the file

        for j in range(len(x)):
            #print(fit_new.rv_data_sets[i][0][j])
            text = b"%s  %s  %s \n"%(bytes(str(x[j]).encode()),bytes(str(0.0).encode()),bytes(str(1.0).encode()) )
            f.write(text)

        f.close()

        obj2.add_dataset("RV_curve", "datafiles/RV_curve",0.0,0.0)  # the last two entries are initial offset and jitter


    os.system("rm datafiles/RV_curve")
    obj2.fitting(outputfiles=[0,1,0], minimize_fortran=True, minimize_loglik=True, amoeba_starts=0, print_stat=False)
    jd        = obj2.fit_results.rv_model.jd
#    rvs       = obj2.fit_results.rv_model.rvs
    o_c       = obj2.fit_results.rv_model.o_c*(-1)

    return np.array([jd,o_c])



#############################

def file_from_path(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)


def run_command_with_timeout(args, secs, output=False, pipe=False): # set output=True if you need to save the output
    '''
    Run a command and kill if it takes too long.
    '''


   # print(args)
    if not (pipe):
        text=tempfile.TemporaryFile() # because PIPE usually has too low capacity
        if 'win' in sys.platform[0:3]:
            proc = Popen(args, shell=True, creationflags=subprocess.CREATE_NEW_PROCESS_GROUP, stdout=text, stderr=text)
        else:
            proc = Popen(args, shell=True, preexec_fn=os.setsid, stdout=text, stderr=text)
    else:
        if 'win' in sys.platform[0:3]:
            proc = Popen(args, shell=True, creationflags=subprocess.CREATE_NEW_PROCESS_GROUP, stdout=PIPE, stderr=PIPE)
        else:
            proc = Popen(args, shell=True, preexec_fn=os.setsid, stdout=PIPE, stderr=PIPE)
   # print(proc)
    proc_thread = Thread(target=proc.wait)
    proc_thread.start()
    proc_thread.join(secs)
    if proc_thread.is_alive():
        #print (proc.pid)
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        except OSError:
            pass
        #print(args)
        print('Process #{} killed after {} seconds'.format(proc.pid, secs))
        flag = -1
        return '',flag
    if not (pipe):
        text.seek(0)
        string_to_output=text.readlines()
    else:
        text=proc.communicate()[0]
        string_to_output=text.splitlines()
       # text.close()
    for i in range(len(string_to_output)):
        string_to_output[i]=string_to_output[i].decode('utf-8').split()
    if not (pipe):
        text.close()
    flag = 1
    if (output):
        return string_to_output,flag # besides the flag which informs about successful termination we also return all the console output in case we want to save it in a variable
    else:
        return '',flag


def run_command_with_timeout_old(args, secs, output=False, pipe=False): # set output=True if you need to save the output
    proc = Popen(args, shell=True, preexec_fn=os.setsid, stdout=PIPE)
    proc_thread = Thread(target=proc.wait)
    proc_thread.start()
    proc_thread.join(secs)
    text = proc.communicate()[0]
    flag = 1
    if proc_thread.is_alive():
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        except OSError:

            print('Process #{} killed after {} seconds'.format(proc.pid, secs))
            flag = -1
            #text = '0 0 0 0'
            return text.decode('utf-8'),flag
    #return proc, flag , text.decode('utf-8')
    return text.decode('utf-8'),flag





def phase_RV_planet_signal(obj,index):

    if not obj.type_fit.get('RV', False):
        obj.ph_data[index] = []
        obj.ph_model[index] = []
        return

    if obj.npl ==0 or len(obj.fit_results.rv_model.jd) ==0:
        return
    else:

        # Advanced indexing below creates independent output arrays.
        sorted_periods = sorted([k for k, v in obj.P.items() if obj.use_planet[k]], key=obj.P.get) 

        rank = sorted_periods.index(index)

        ############ phase fold fix for sparse model ######
        model_time_phase = np.array( (obj.fit_results.model_jd -obj.fit_results.model_jd[0] + (obj.fit_results.model_jd[0] - obj.epoch) )%obj.P[index] )

        model_shift = obj.P[index] - (obj.fit_results.rv_model.jd[0] - obj.epoch )%obj.P[index]

        model_time_phase = (model_time_phase + model_shift)%obj.P[index]

        sort = sorted(range(len(model_time_phase)), key=lambda k: model_time_phase[k])
        model_time_phase  = model_time_phase[sort]



        phased_model      = obj.fit_results.res['fit'].T[2+rank][sort]

        ############ phase data ######
        data_time_phase = np.array( (obj.fit_results.rv_model.jd  - obj.fit_results.rv_model.jd[0])%obj.P[index] )

        sort = sorted(range(len(data_time_phase)), key=lambda k: data_time_phase[k])
        data_time_phase      = data_time_phase[sort]
        phased_data          = obj.fit_results.model_data[6+rank][sort] + obj.fit_results.rv_model.o_c[sort]  
        
        phased_data_o_c      = obj.fit_results.rv_model.o_c[sort]#  - obj.fit_results.rv_model.rvs[sort]
        phased_data_err      = obj.fit_results.rv_model.rv_err[sort]
        phased_data_idset    = obj.fit_results.idset[sort]


        if obj.doGP == True:
            phased_data = phased_data - obj.gp_model_data[0][sort]
            phased_data_o_c = phased_data_o_c - obj.gp_model_data[0][sort]

        model = [model_time_phase,  phased_model]
        data  = [data_time_phase,  phased_data, phased_data_err, phased_data_idset,phased_data_o_c]



        #####################
        obj.ph_data[index] = data
        obj.ph_model[index] = model

        return data, model



def phase_RV_planet_signal_old(obj,planet):

    if obj.npl ==0 or len(obj.fit_results.rv_model.jd) ==0:
        return
    else:

        copied_obj = dill.copy(obj)

        if(copied_obj.mod_dynamical):
            copied_obj.mod_dynamical = False

        index = planet - 1
        #print(index)
        ############################################
        ############# here is the trick!  ##########
        ############################################
        pp0 =  copied_obj.K[index]  # we define a variable to be the planet amplitude Kj

        #copied_obj.params.planet_params[7*index+0] = 0.0000001 # then we set Kj to be 0, i.e. remove the j-th planet signal
        copied_obj.K[index] = 0.0000001 # then we set Kj to be 0, i.e. remove the j-th planet signal
        copied_obj.fitting(minimize_loglik=True, amoeba_starts=0,
                           outputfiles=[1,1,1],return_flag=False, npoints=int(obj.model_npoints),
                          # model_max=int(max(obj.fit_results.model_jd)-max(copied_obj.fit_results.rv_model.jd)),
                          # model_min=int(copied_obj.epoch -min(obj.fit_results.model_jd)))
                           model_max=int(copied_obj.model_max),
                           model_min=int(copied_obj.model_min))
        # and we create the static Nplanet model for the data and the model curve
        # now this model residuals will contain ONLY the j-th planet signal + the best fit residuals

        #copied_obj.params.planet_params[7*index+0] = pp0 # we restore Kj to its best fit value.
        copied_obj.K[index] = pp0 # we restore Kj to its best fit value.

        ############################################
        #########      trick is over      ##########
        ############################################


        ############ phase fold fix for sparse model ######
        model_time_phase = np.array( (copied_obj.fit_results.model_jd -copied_obj.fit_results.model_jd[0] + (copied_obj.fit_results.model_jd[0] - copied_obj.epoch) )%copied_obj.P[index] )

        model_shift = copied_obj.P[index] - (copied_obj.fit_results.rv_model.jd[0] - copied_obj.epoch )%copied_obj.P[index]

        model_time_phase = (model_time_phase + model_shift)% copied_obj.P[index]

        sort = sorted(range(len(model_time_phase)), key=lambda k: model_time_phase[k])
        model_time_phase  = model_time_phase[sort]

        phased_model      = obj.fit_results.model[sort] - copied_obj.fit_results.model[sort]

        ############ phase data ######
        data_time_phase = np.array( (copied_obj.fit_results.rv_model.jd  - copied_obj.fit_results.rv_model.jd[0])% copied_obj.P[index] )

        sort = sorted(range(len(data_time_phase)), key=lambda k: data_time_phase[k])
        data_time_phase      = data_time_phase[sort]
        phased_data          = copied_obj.fit_results.rv_model.o_c[sort]#  - copied_obj.fit_results.rv_model.rvs[sort]
        phased_data_o_c      = obj.fit_results.rv_model.o_c[sort]#  - copied_obj.fit_results.rv_model.rvs[sort]
        phased_data_err      = copied_obj.fit_results.rv_model.rv_err[sort]
        phased_data_idset    = copied_obj.fit_results.idset[sort]


        if copied_obj.doGP == True:
            phased_data = phased_data - copied_obj.gp_model_data[0][sort]
            phased_data_o_c = phased_data_o_c - copied_obj.gp_model_data[0][sort]

        model = [model_time_phase,  phased_model]
        data  = [data_time_phase,  phased_data, phased_data_err, phased_data_idset,phased_data_o_c]


        del copied_obj

        #####################
        obj.ph_data[planet-1] = data
        obj.ph_model[planet-1] = model

        return data, model




def find_planets(obj,fend=0.75,RV_jitter=False):

    power_levels = np.array([0.1,0.01,0.001])

    # check if RV data is present
    if obj.ndset <= 0:
         return

    # the first one on the data GLS
    if obj.gls.power.max() <= obj.gls.powerLevel(obj.auto_fit_FAP_level):
         return obj

    else:
        if obj.npl !=0:
            #for j in range(obj.npl):

            for j in range(9):
                if not bool(obj.use_planet[j]):
                    continue
                obj.remove_planet(j)

        mean_anomaly_from_gls = np.degrees((((obj.epoch - float(obj.gls.hpstat["T0"]) )% (obj.gls.hpstat["P"]) )/ (obj.gls.hpstat["P"]) ) * 2*np.pi)

        obj.add_planet(obj.gls.hpstat["amp"],obj.gls.hpstat["P"],0.0,0.0,mean_anomaly_from_gls -90.0,90.0,0.0, useK=True,useP=True,usee=obj.auto_fit_allow_ecc,usew=obj.auto_fit_allow_ecc,useM0=True,usei=False, usecap=False)

        obj.fitting(fileinput=False,outputfiles=[1,1,1], doGP=False,   minimize_fortran=True, fortran_kill=3, timeout_sec= 3)

        run_gls_o_c(obj,fend=fend, RV_jitter=True)#,as_main=True)

        

        #now inspect the residuals

        for i in range(1,int(obj.auto_fit_max_pl)):

            if obj.gls_o_c.power.max() <= obj.gls_o_c.powerLevel(obj.auto_fit_FAP_level):
                #for j in range(obj.npl):

                for j in range(9):
                    if not bool(obj.use_planet[j]):
                        continue

                    #obj.use.update_use_planet_params_one_planet(j,True,True,obj.auto_fit_allow_ecc,obj.auto_fit_allow_ecc,True,False,False)
                    obj.P_use[j] = True
                    obj.K_use[j] = True
                    obj.e_use[j] = obj.auto_fit_allow_ecc
                    obj.w_use[j] = obj.auto_fit_allow_ecc
                    obj.M0_use[j] = True
                    obj.i_use[j]  = False
                    obj.Node_use[j] = False


                    obj.fitting(fileinput=False,outputfiles=[1,1,1], doGP=False,   minimize_fortran=True, fortran_kill=10, timeout_sec= 10)
                    obj = run_gls_o_c(obj,fend=fend, RV_jitter=RV_jitter)
                return obj

            else:
                mean_anomaly_from_gls = np.degrees((((obj.epoch - float(obj.gls_o_c.hpstat["T0"]) )% (obj.gls_o_c.hpstat["P"]) )/ (obj.gls_o_c.hpstat["P"]) ) * 2*np.pi)


                obj.add_planet(obj.gls_o_c.hpstat["amp"],obj.gls_o_c.hpstat["P"],0.0,0.0,mean_anomaly_from_gls -90.0,90.0,0.0, useK=True,useP=True,usee=obj.auto_fit_allow_ecc,usew=obj.auto_fit_allow_ecc,useM0=True,usei=False, usecap=False)

                obj.fitting(fileinput=False,outputfiles=[1,1,1], doGP=False,   minimize_fortran=True, fortran_kill=3, timeout_sec= 3)
                run_gls_o_c(obj,fend=fend, RV_jitter=RV_jitter)
            #else:
             #   continue

        #for j in range(obj.npl):

        for j in range(9):
            if not bool(obj.use_planet[j]):
                continue

            #obj.use.update_use_planet_params_one_planet(j,True,True,obj.auto_fit_allow_ecc,obj.auto_fit_allow_ecc,True,False,False)
            obj.P_use[j] = True
            obj.K_use[j] = True
            obj.e_use[j] = obj.auto_fit_allow_ecc
            obj.w_use[j] = obj.auto_fit_allow_ecc
            obj.M0_use[j] = True
            obj.i_use[j]  = False
            obj.Node_use[j] = False

        obj.fitting(fileinput=False,outputfiles=[1,1,1], doGP=False,   minimize_fortran=True, fortran_kill=3, timeout_sec= 3)
        run_gls_o_c(obj,fend=fend,RV_jitter=RV_jitter)
    return obj





def run_gls(obj,fend =1.0,fbeg=10000, ofac=10, RV_jitter=False):


    if RV_jitter:
        error_list = add_jitter(obj, obj.fit_results.rv_model.rv_err, obj.fit_results.idset)
    else:
        error_list = obj.fit_results.rv_model.rv_err

    #fbeg = abs(max(obj.fit_results.rv_model.jd)-min(obj.fit_results.rv_model.jd))  * 2.0
    #omega = 1/ np.logspace(np.log10(fend), np.log10(fbeg), num=int(1000))



    if len(obj.fit_results.rv_model.jd) > 5:
        #RV_per = gls.Gls((obj.fit_results.rv_model.jd, obj.fit_results.rv_model.rvs, obj.fit_results.rv_model.rv_err),
       # fast=True,  verbose=False, norm='ZK',ofac=10, fbeg=omega[-1], fend=omega[0],)

        RV_per = gls.Gls((obj.fit_results.rv_model.jd, obj.fit_results.rv_model.rvs, error_list), fast=True,  verbose=False, norm='ZK',ofac=ofac, fbeg=1/fbeg, fend=1/fend)  

        obj.gls = RV_per
    else:
        return obj

    return obj

 

def run_gls_o_c(obj,fend =1.0,fbeg=10000, as_main = False, ofac=10, RV_jitter=False):

    #fbeg = abs(max(obj.fit_results.rv_model.jd)-min(obj.fit_results.rv_model.jd))  * 2.0
    #omega = 1/ np.logspace(np.log10(fend), np.log10(fbeg), num=int(1000))


    if RV_jitter:
        error_list = add_jitter(obj, obj.fit_results.rv_model.rv_err, obj.fit_results.idset)
    else:
        error_list = obj.fit_results.rv_model.rv_err


    if len(obj.fit_results.rv_model.jd) > 5:
        #RV_per_res = gls.Gls((obj.fit_results.rv_model.jd, obj.fit_results.rv_model.o_c, obj.fit_results.rv_model.rv_err),
        #fast=True,  verbose=False, norm='ZK', ofac=10, fbeg=omega[-1], fend=omega[0],)

        RV_per_res = gls.Gls((obj.fit_results.rv_model.jd, obj.fit_results.rv_model.o_c, error_list), fast=True,  verbose=False, norm='ZK',ofac=ofac, fbeg=1/fbeg, fend=1/fend)  


        if as_main == False:
            obj.gls_o_c = RV_per_res
        elif as_main == True:
            obj.gls = RV_per_res
    else:
        return obj

    return obj


def is_float(n):
    '''
    Given a string n, verify if it expresses a valid float.
    Casting n to string in case an object of type float or similar is given as an argument
    '''
    return re.match(r'^-?\d*(\.\d+)?(E-?\d+)?$', str(n))

# Given a float or string, verify if it expresses an integer. Possible to introduce upper and lower bounds and if the inequalities on either side should be strong or weak .
def is_int(s,bounded=[False,False],bounds=[0,0],equal=[False,False]):
    if is_float(s): # if it is an int, it is certainly float as well
        n=float(s) # we need n as a number, not as a string, for comparisons with bounds later
        is_an_int=float(s).is_integer()
    else:
        is_an_int=False
    # is_an_int now contains an information if s is an int, but without bounds. Let's introduce bounds:
    if(is_an_int): # if it's not an int at all we don't need to check any further
        if(bounded[0]): # if there is a lower bound let's apply it
            if (n<bounds[0] or (not equal[0] and n==bounds[0])):
                is_an_int=False
    if(is_an_int): # if the lower bound returned False we don't need to check any further
        if(bounded[1]): # if there is a lower bound let's apply it
            if (n>bounds[1] or (not equal[1] and n==bounds[1])):
                is_an_int=False
    return is_an_int

# If save_wrong_lines is enabled we will save a string 'wrong_line' instead of this line and save indices at which this occurred, otherwise we will skip this line
def convert_array_to_float(a,save_wrong_lines=False):
    converting_warnings=Warning_log([],'Converting array to float')
    b=[]
    if (save_wrong_lines):
        wrong_indices=[]
    for i in range(len(a)):
        if not is_float(a[i]):
            if not (save_wrong_lines):
                converting_warnings.update_warning_list('Array passed to convert_array_to_float function should only contain floats! Line %d skipped'%(i+1))
            else:
                b.append('wrong_line')
                wrong_indices=np.concatenate((wrong_indices,np.atleast_1d(i)))
        else:
            b.append(float(a[i]))

    converting_warnings.print_warning_log()
    if (save_wrong_lines):
        return np.array(b),wrong_indices
    else:
        return np.array(b)

def convert_array_to_int(a, save_wrong_lines=False):
    converting_warnings=Warning_log([],'Converting array to int')
    b=[]
    if (save_wrong_lines):
        wrong_indices=[]
    for i in range(len(a)):
        if not is_int(a[i]):
            if not (save_wrong_lines):
                converting_warnings.update_warning_list('Array passed to convert_array_to_int function should only contain ints! Line %d skipped'%(i+1))
            else:
                b.append('wrong_line')
                wrong_indices=np.concatenate((wrong_indices,np.atleast_1d(i)))
        else:
            b.append(int(a[i]))

    converting_warnings.print_warning_log()
    if (save_wrong_lines):
        return np.array(b),wrong_indices
    else:
        return np.array(b)

#for convenient reading of the input file
def read_file_as_array_of_arrays(inputfile):
    a=open(inputfile, 'r')
    b=a.readlines() # b as array of strings
    a.close()
    c=[]
    ic=0 # iterator for values in c
    for i in range(len(b)):
        b[i]=np.atleast_1d(b[i].split()) # turn a row of b into an array of arrays
        c.append([]) # need to make a separate array so every element is of correct type
        # convert each string that represents a float into float
        for j in range(0,len(b[i])):
            if (is_float(b[i][j])):
                c[ic].append(float(b[i][j]))
            elif not (b[i][j][-1]==':'): # ignore comments, which can be place by the user as strings which end with a collon, in the comments use underline instead of space or an error will arise
                c[ic].append(b[i][j])
        ic=ic+1
    #c = np.array(c, dtype=float)

    return c


#for convenient reading of the input file  the second is a hack so the mcmc lnL line is skipped! TBFixed
def read_file_as_array_of_arrays_mcmc(inputfile):
    a=open(inputfile, 'r')
    b=a.readlines() # b as array of strings
    a.close()
    c=[]
    ic=0 # iterator for values in c
    for i in range(len(b)):
        b[i]=np.atleast_1d(b[i].split()) # turn a row of b into an array of arrays
        c.append([]) # need to make a separate array so every element is of correct type
        # convert each string that represents a float into float
        for j in range(1,len(b[i])):
            if (is_float(b[i][j])):
                c[ic].append(float(b[i][j]))
            elif not (b[i][j][-1]==':'): # ignore comments, which can be place by the user as strings which end with a collon, in the comments use underline instead of space or an error will arise
                c[ic].append(float(b[i][j]))
        ic=ic+1

    c = np.array(c, dtype=float)
    return c



def verify_array_with_bounds(ar,bounds):
    '''Verify if values of array ar fit withind declared bounds, if too many/too few bounds do as much as can be done'''
    if (len(ar)<=len(bounds)):
        num=len(ar) # number of values to check
    elif (len(ar)>len(bounds)):
        num=len(bounds) # number of values to check
    verification=True # initial value
    for i in range(num):
        # check if some of the values doesn't fit in the bounds, if so return False
        if (ar[i]<bounds[i][0] or ar[i]>bounds[i][1]):
            verification=False
            break

    return verification





def _latex_table_parameters(obj):
    """Shared fitted-parameter selection and labels for both LaTeX exports."""
    planets = [i for i, enabled in enumerate(obj.use_planet) if enabled]
    types = obj.type_fit
    rv, transit = types.get('RV', False), types.get('Transit', False)
    orbital = rv or types.get('TTV', False) or types.get('AST', False)
    rows = []

    def add(name, label, slots=planets, scalar=False, bounds=None):
        rows.append((name, label, slots, scalar, bounds or name + '_bound'))

    if orbital:
        add('K', r'$K$ [m\,s$^{-1}$]')
    add('P', r'$P$ [day]')
    if getattr(obj, 'hkl', False):
        add('e_sinw', r'$e\sin\omega$')
        add('e_cosw', r'$e\cos\omega$')
        if orbital:
            add('lamb', r'$\lambda$ [deg]')
    else:
        add('e', r'$e$')
        add('w', r'$\omega$ [deg]')
        if orbital:
            add('M0', r'$M_0$ [deg]')
    if getattr(obj, 'mod_dynamical', False) or transit or types.get('AST', False):
        add('i', r'$i$ [deg]')
    if getattr(obj, 'mod_dynamical', False) or types.get('AST', False):
        add('Node', r'$\Omega$ [deg]')
    if transit:
        add('t0', r'$t_0$ [day]')
        add('pl_rad', r'$R_{\rm p}/R_\star$')
        add('pl_a', r'$a/R_\star$')
    if hasattr(obj, 'omega_dot_use') and any(obj.omega_dot_use[i] for i in planets):
        add('omega_dot', r'$\dot\omega$ [deg\,yr$^{-1}$]', bounds='omega_dot_bounds')
    for active, prefix, datasets in ((rv, '', 'rv_data_sets'), (transit, 'tra_', 'tra_data_sets')):
        if not active:
            continue
        if rv and not prefix:
            add('rv_lintr', r'RV linear trend [m\,s$^{-1}$\,day$^{-1}$]', [0], True, 'rv_lintr_bounds')
            add('rv_quadtr', r'RV quadratic trend [m\,s$^{-1}$\,day$^{-2}$]', [0], True, 'rv_quadtr_bounds')
        data = getattr(obj, datasets, {})
        slots = [i for i in range(len(data)) if len(data[i])]
        for name, label in ((('rvoff', r'RV offset [m\,s$^{-1}$]'), ('jitt', r'RV jitter [m\,s$^{-1}$]'))
                            if not prefix else (('tra_off', 'Transit offset'), ('tra_jitt', 'Transit jitter'))):
            for i in slots:
                add(name, label + ' %s' % (i + 1), [i], bounds=name + '_bounds')
        if getattr(obj, prefix + 'doGP', False):
            kernel = getattr(obj, prefix + 'gp_kernel', '')
            suffix = {'RotKernel': 'rot', 'SHOKernel': 'sho', 'Matern32': 'mat',
                      'RealTerm': 'drw', 'dSHOKernel': 'double_sho'}.get(kernel)
            if suffix:
                base = prefix + 'GP_' + suffix
                values = getattr(obj, base + '_params', [])
                labels = getattr(obj, base + '_str', [])
                for i in range(len(values)):
                    label = labels[i] if i < len(labels) else base + ' %s' % (i + 1)
                    add(base + '_params', label, [i], bounds=base + '_bounds')
    return planets, rows


def _latex_number(value, precision):
    """Avoid turning small nonzero uncertainties/priors into printed zeroes."""
    value = float(value)
    if not np.isfinite(value):
        return r'\mathrm{--}'
    if value != 0 and abs(value) < .5 * 10. ** (-precision):
        mantissa, exponent = f'{value:.{precision}e}'.split('e')
        return mantissa + r'\times10^{' + str(int(exponent)) + '}'
    return f'{value:.{precision}f}'


def _latex_fit_cell(value, error, precision, asymmetric):
    number = _latex_number(value, precision)
    if error is None:
        return '$' + number + '$'
    errors = np.abs(np.asarray(error, dtype=float).reshape(-1))
    if not errors.size or not np.isfinite(errors).all():
        return '$' + number + '$'
    lower, upper = errors[0], errors[-1]
    if max(lower, upper) == 0:
        return '$' + number + '$'
    use_asymmetric = asymmetric is True or (asymmetric is None and
                                           not np.isclose(lower, upper, rtol=1e-6, atol=0))
    if use_asymmetric:
        return '$' + number + '^{+' + _latex_number(upper, precision) + '}_{-' + _latex_number(lower, precision) + '}$'
    return '$' + number + r'\pm ' + _latex_number(max(lower, upper), precision) + '$'


def _latex_parameter_error(obj, name, slot, scalar):
    if name.startswith('GP_') and name.endswith('_params'):
        packed = getattr(getattr(obj, 'param_errors', None), 'GP_params_errors', None)
        if packed is not None and slot < len(packed):
            return packed[slot]
    error_name = name[:-7] + '_err' if name.endswith('_params') else name + '_err'
    errors = getattr(obj, error_name, None)
    if errors is not None:
        return errors if scalar else errors[slot]
    # Compatibility with older fit objects that only stored packed errors.
    if name in ('K', 'P', 'e', 'w', 'M0', 'i', 'Node'):
        packed = getattr(getattr(obj, 'param_errors', None), 'planet_params_errors', None)
        if packed is not None:
            return packed[7 * slot + ('K', 'P', 'e', 'w', 'M0', 'i', 'Node').index(name)]
    return None


def _latex_render_table(planets, rows, width, prior=False):
    columns = max(1, len(planets))
    headings = ['Planet ' + chr(98 + i) for i in planets] if planets else ['Value']
    text = [r'\begin{table}[ht]', r'\centering', r'\caption{}',
            r'\label{table:priors}' if prior else r'\label{table:parameters}',
            r'\begin{tabular}{l' + 'r' * columns + '}', r'\hline\hline',
            'Parameter & ' + ' & '.join(headings) + r' \\', r'\hline']
    for label, cells in rows:
        # Dataset/global parameters occupy the first value column; pad the
        # others so every row has exactly the advertised number of columns.
        cells = cells + [''] * (columns - len(cells))
        text.append(label.ljust(30) + ' & ' + ' & '.join(cell.rjust(width) for cell in cells) + r' \\')
    text += [r'\hline', r'\end{tabular}']
    if prior:
        text.append(r'\par\small $\mathcal{N}(\mu,\sigma)$ denotes a Gaussian; '
                    r'$\mathcal{J}(a,b)$ a log-uniform prior; all are truncated to the listed bounds.')
    text.append(r'\end{table}')
    return '\n'.join(text) + '\n'


def _latex_save_table(text, path, file_name, return_text):
    if return_text:
        return text
    from pathlib import Path
    (Path(path) / file_name).write_text(text, encoding='utf-8')
    print('Done')


def latex_pl_param_table(obj, width=10, precision=2, asymmetric=None,
                         file_name='test.tex', path='./', return_text=False):
    """Export fit-specific parameters, automatically preserving unequal errors.

    asymmetric=None selects per-cell formatting; False uses the larger error
    with +/-; True explicitly requests upper/lower error notation.
    Zero/missing uncertainties are displayed as values, not invented errors.
    """
    if asymmetric not in (None, True, False):
        raise ValueError('asymmetric must be None, True, or False')
    planets, parameters = _latex_table_parameters(obj)
    rows = []
    for name, label, slots, scalar, bounds in parameters:
        values = getattr(obj, name, None)
        if values is None:
            continue
        cells = [_latex_fit_cell(values if scalar else values[i],
                                 _latex_parameter_error(obj, name, i, scalar), precision, asymmetric)
                 for i in slots]
        rows.append((label, cells))
    orbital = any(obj.type_fit.get(k, False) for k in ('RV', 'TTV', 'AST'))
    if orbital:
        for name, label in (('semimajor', r'$a$ [au]'),
                            ('masses', r'$m$ [$M_{\rm Jup}$]' if getattr(obj, 'mod_dynamical', False)
                             or obj.type_fit.get('AST', False) else r'$m\sin i$ [$M_{\rm Jup}$]')):
            values = getattr(obj, name, None)
            if values is not None and len(values):
                projected = name == 'masses' and not getattr(obj, 'mod_dynamical', False) and not obj.type_fit.get('AST', False)
                rows.append((label, [_latex_fit_cell(
                    values[i] * abs(np.sin(np.radians(obj.i[i]))) if projected else values[i],
                    None, precision, asymmetric) for i in planets]))
        if hasattr(obj, 'epoch'):
            rows.append(('Epoch [day]', [_latex_fit_cell(obj.epoch, None, precision, asymmetric)]))
            if hasattr(obj, 'M0') and hasattr(obj, 'P'):
                rows.append((r'$t_{\rm peri}$ [day]', [
                    _latex_fit_cell(float(obj.epoch) - np.radians(obj.M0[i]) / (2 * np.pi) * obj.P[i],
                                    None, precision, asymmetric) for i in planets]))
    results = getattr(obj, 'fit_results', None)
    if obj.type_fit.get('RV', False) and results is not None:
        for name, label in (('chi2', r'$\chi^2$'), ('reduced_chi2', r'$\chi_\nu^2$'),
                            ('rms', r'RMS [m\,s$^{-1}$]'), ('wrms', r'WRMS [m\,s$^{-1}$]')):
            if hasattr(results, name):
                rows.append((label, [_latex_fit_cell(getattr(results, name), None, precision, asymmetric)]))
        if hasattr(obj, 'loglik'):
            rows.append((r'$\ln\mathcal{L}$', [_latex_fit_cell(obj.loglik, None, precision, asymmetric)]))
        if hasattr(results, 'jd'):
            rows.append((r'$N_{\rm RV}$', [str(len(results.jd))]))
    return _latex_save_table(_latex_render_table(planets, rows, width), path, file_name, return_text)


def latex_prior_table(obj, width=10, precision=2, file_name='prior_table.tex',
                      path='./', return_text=False):
    """Export the same fitted parameters with distributions and hard bounds.

    Gaussian notation uses standard deviation (not variance), including the
    truncation actually used by the nested sampler.
    """
    planets, parameters = _latex_table_parameters(obj)
    rows = []
    for name, label, slots, scalar, bounds_name in parameters:
        bounds = getattr(obj, bounds_name, None)
        if bounds is None:
            continue
        base = name[:-7] if name.endswith('_params') else name
        normal, jeff = getattr(obj, base + '_norm_pr', None), getattr(obj, base + '_jeff_pr', None)
        cells = []
        for i in slots:
            lo, hi = bounds[i]
            limits = '[' + _latex_number(lo, precision) + ', ' + _latex_number(hi, precision) + ']'
            if normal is not None and normal[i][2]:
                mu, sigma = normal[i][:2]
                cell = r'\mathcal{N}(' + _latex_number(mu, precision) + ', ' + _latex_number(sigma, precision) + r')\;\in ' + limits
            elif jeff is not None and jeff[i][2]:
                a, b = jeff[i][:2]
                cell = r'\mathcal{J}(' + _latex_number(a, precision) + ', ' + _latex_number(b, precision) + r')\;\in ' + limits
            else:
                cell = r'\mathcal{U}(' + _latex_number(lo, precision) + ', ' + _latex_number(hi, precision) + ')'
            cells.append('$' + cell + '$')
        rows.append((label, cells))
    return _latex_save_table(_latex_render_table(planets, rows, width, prior=True), path, file_name, return_text)


def _f_test_parameter_count(obj):
    count = getattr(obj.fit_results, 'mfit', None)
    if count is not None:
        return count
    from fit_statistics import component_parameter_counts
    prepared = dill.copy(obj)
    prepared.prepare_for_mcmc(rtg=prepared.rtg)
    return component_parameter_counts(prepared)['RV']


def f_test(obj, obj2 = None, alpha = 0.01, lnL=False):


    chi2 = obj.fit_results.chi2
    ndata = len(obj.fit_results.jd)
    if chi2 is None or not ndata:
        print("Nothing to test. Fit RV data first.")
        return
    par2 = _f_test_parameter_count(obj)
   #     self.value_reduced_chi2.setText("%.4f"%(fit.fit_results.reduced_chi2))
        #self.value_loglik.setText("%.4f"%(fit.fit_results.loglik))
   #     self.value_loglik.setText("%.4f"%(fit.loglik))


    if obj2 is None:
        obj2 = dill.copy(obj)
        obj2.npl = 0
        obj2.use_planet = [0] * len(obj2.use_planet)
        if obj2.fitting(return_flag=True) != 1:
            print("Cannot perform F-test: null model fitting failed.")
            return
    else:
        obj2 = dill.copy(obj2)

        if len(obj.fit_results.jd) != len(obj2.fit_results.jd):
            print("not the same data, test makes no sense")
            return


    chi1 = obj2.fit_results.chi2
    par1 = _f_test_parameter_count(obj2)

    if par2 <= par1:
        print("Nothing to test. The tested model must have more free parameters than the null model.")
        return
    if ndata <= par2 or not np.isfinite(chi2) or chi2 <= 0:
        print("Cannot perform F-test: invalid chi-squared or no residual degrees of freedom.")
        return

    chi2_red = chi2/(ndata - par2)


    if lnL == True:
        F = 2*(obj.loglik - obj2.loglik) # in case \Delta lnL must be tested.
    else:
        if abs(par2-par1) > 0:
            F = ((chi1 - chi2)/(par2-par1))/chi2_red # else standard f-test
        else:
            print("Nothing to test. The Tested model has == or < Npar. than the null model.")
            return

    p_value = pdf.f.sf(F, par2 - par1, ndata - par2, loc=0, scale=1)

    print("""
\chi^2 null model = %s
\chi^2 tested model = %s
lnL null model = %s
lnL tested model = %s
N parametrs null model = %s
N parametrs tested model = %s
F value = %s
p-value = %s
alpha value = %s
"""%(chi1,chi2,obj2.loglik,obj.loglik,par1,par2,F,p_value,alpha))


    if float(p_value) < alpha:
        print("Null hypothesis rejected")
        print("Probability = ", (1.0-float(p_value))*100.0,'%')
    else:
        print("Null hypothesis cannot be rejected")



def plot_gp(obj, curve=False):

    import matplotlib.pyplot as plt

    color="#ff7f0e"
    colors = ['b','g','r']

    x     = obj.fit_results.rv_model.jd
    y     = obj.fit_results.rv_model.o_c
    y_err = obj.fit_results.rv_model.rv_err
    idset = obj.filelist.idset


    if curve==True:
        x_model = np.linspace(min(x), max(x), 5000) #obj.fit_results.model_jd
        mu,var,std = obj.gp_model_curve

    else:
        x_model = x
        mu,var,std = obj.gp_model_data

    #print(mu[0:10])
    #print(y[0:10])

    for i in range(obj.filelist.ndset):
        plt.errorbar(x[idset==i],y[idset==i], yerr=y_err[idset==i], fmt=".",color=colors[i],  capsize=0);

    plt.plot(x_model, mu, color = '0.5' );
    plt.fill_between(x_model ,mu+std,  mu-std, color=color, alpha=0.3, edgecolor="none")



def plot_transit_gp(obj, curve=False):

    import matplotlib.pyplot as plt

    color="#ff7f0e"
    colors = ['b','g','r']

    x     = obj.tra_data_sets[0][0]
    y     = obj.tra_data_sets[0][1]
    y_err = obj.tra_data_sets[0][2]
    #idset = obj.filelist.idset


    if curve==True:
        x_model = np.linspace(min(x), max(x), 5000) #obj.fit_results.model_jd
        mu,var,std = obj.tra_gp_model_curve

    else:
        x_model = x
        mu,var,std = obj.tra_gp_model_data

    #print(mu[0:10])
    #print(y[0:10])

    #for i in range(obj.filelist.ndset):
        #plt.errorbar(x[idset==i],y[idset==i], yerr=y_err[idset==i], fmt=".",color=colors[i],  capsize=0);
    plt.errorbar(x,y, yerr=y_err, fmt=".",color=colors[0],  capsize=0);

    plt.plot(x_model, mu, color = '0.5' );
    plt.fill_between(x_model ,mu+std,  mu-std, color=color, alpha=0.3, edgecolor="none")


 
 
 


####################### mass_semimajor ###########################################
def mass_a_from_Kepler_fit(a, npl, m0, incl=None):
    """Legacy five-element interface; incl defaults to edge-on degrees.

    Pass incl for physical masses. Prefer mass_a_from_jacobi_el when all
    seven orbital elements are available. Outputs are Mjup and Jacobi AU.
    """
    rows = np.asarray(a, dtype=float).reshape(npl, 5)
    inclinations = np.broadcast_to(90.0 if incl is None else incl, (npl,))
    elements = np.column_stack((rows, inclinations, np.zeros(npl)))
    return mass_a_from_jacobi_el(elements, npl, m0)



def check_swift(path='./'):
    """Ensure the Python SWIFT extension is built for the current interpreter."""
    if __package__ == 'RV_mod':
        from swift import is_available
    else:
        from ..swift import is_available
    if not is_available():
        builder = Path(__file__).resolve().parents[2] / 'compile_swift.py'
        subprocess.run([sys.executable, str(builder)], check=True)
        if not is_available():
            raise RuntimeError("SWIFT was built but its Python extension cannot be imported")
    return True


def run_stability(obj, timemax=3000.0, timestep=10, timeout_sec=1000.0,
                  stab_save_dir='run', remove_stab_save_dir=True,
                  integrator='symba', output_step=1.0,
                  use_fixed_steps=False, fixed_steps=1000):
    """Run the SWIFT Python library and populate the stability plot arrays.

    Initial conditions and histories pass directly as arrays to swift.f95.
    Integration uses no files or worker processes; native sessions are serialized.
    Planet results retain the existing order of increasing orbital period.
    """
    from time import perf_counter
    start_time = perf_counter()
    if __package__ == 'RV_mod':
        from swift import integrate
    else:
        from ..swift import integrate
    used_indices = sorted((j for j in range(9) if obj.use_planet[j]),
                          key=lambda j: obj.P[j])
    if len(used_indices) != obj.npl:
        raise ValueError("obj.npl does not match the enabled planets")
    planets = [[obj.masses[j] / 1047.5654817267318, obj.semimajor[j],
                obj.e[j], obj.i[j], obj.w[j], obj.Node[j], obj.M0[j]]
               for j in used_indices]
    save_dir = None
    if not remove_stab_save_dir:
        save_dir = Path(stab_save_dir)
        if not save_dir.is_absolute():
            save_dir = Path.cwd() / 'stability' / save_dir
    result = integrate(
        obj.params.stellar_mass, planets, timemax=timemax,
        timestep=timestep, timeout_sec=timeout_sec, integrator=integrator,
        output_step=output_step, use_fixed_steps=use_fixed_steps,
        fixed_steps=fixed_steps, gr_step=getattr(obj, 'GR_step', 1),
        save_dir=save_dir)
    from integration_period import record_period_source
    record_period_source(obj, used_indices)
    obj.swift_evol = result
    if not result.get('completed', True):
        print(result['status'])
    energy = result['energy']
    obj.evol_T_energy = energy[:, 0] / 365.25
    obj.evol_energy = energy[:, 1]
    for column, axis in enumerate(('lx', 'ly', 'lz'), start=2):
        obj.evol_momentum[axis] = energy[:, column]
    for k, orbit in enumerate(result['orbits']):
        obj.evol_T[k] = orbit[:, 0] / 365.25
        for attribute, column in (('evol_a', 2), ('evol_e', 3),
                                  ('evol_i', 4), ('evol_Om', 5),
                                  ('evol_p', 6), ('evol_M', 7)):
            getattr(obj, attribute)[k] = orbit[:, column]
        obj.evol_Per[k] = a_to_P(obj.evol_a[k], obj.params.stellar_mass)
    print(f"--- SWIFT {integrator} N-body: {perf_counter() - start_time:.5f} seconds ---")
    return obj



def run_stability_rebound(
    obj,
    timemax=3000.0,
    timestep=10.0 / 365.25,
    output_step=1.0,
    integrator="whfast",
    exact_finish_time=0,
    collision_distance=None,
    escape_distance=None,
    move_to_com=True,
    sort_by_period=True,
    safe_mode=0,
    corrector=11,
    verbose=True,
    coordinates='JACOBI',
    use_fixed_steps=False,
    fixed_steps=1000,
):
    """
    Run an N-body stability integration with REBOUND.

    This version follows the Exo-Striker run_stability style:
        input  : obj
        output : obj

    It fills:
        obj.evol_T
        obj.evol_a
        obj.evol_e
        obj.evol_i
        obj.evol_Om
        obj.evol_p
        obj.evol_M
        obj.evol_Per
        obj.evol_T_energy
        obj.evol_energy
        obj.evol_momentum

    It also stores the full REBOUND output dictionary in:
        obj.rebound_evol

    Units:
        length = au
        time   = yr
        mass   = Msun

    Input angles are assumed to be in degrees.
    Output angles are in degrees.
    Periods in obj.evol_Per are converted to days, matching Exo-Striker usage.
    """
    from time import perf_counter
    start_time = perf_counter()

    try:
        import rebound
    except ImportError:
        raise ImportError(
            "REBOUND is not installed. Install it with: pip install rebound"
        )

    if verbose:
        try:
            print("REBOUND version:", rebound.__version__)
        except Exception:
            print("REBOUND version: unknown")

    integrator = str(integrator).lower()

    mjup_to_msun = 1.0 / 1047.5654817267318
    deg_to_rad = np.pi / 180.0
    rad_to_deg = 180.0 / np.pi

    used_indices = [j for j in range(9) if bool(obj.use_planet[j])]

    if sort_by_period:
        used_indices.sort(key=lambda j: obj.P[j])

    npl = len(used_indices)

    if npl == 0:
        obj.rebound_evol = {
            "stable": False,
            "status": "no active planets",
            "used_indices": used_indices,
        }

        if verbose:
            print("REBOUND stability aborted: no active planets.")

        return obj

    timemax = float(timemax)
    timestep = float(timestep)

    if use_fixed_steps:
        output_step = timemax / float(fixed_steps)
    else:
        output_step = float(output_step)

    if timemax <= 0.0:
        raise ValueError("timemax must be positive.")

    if timestep <= 0.0:
        raise ValueError("timestep must be positive.")

    if output_step <= 0.0:
        raise ValueError("output_step must be positive.")

    sim = rebound.Simulation()

    sim.G = 4.0 * np.pi * np.pi
 
    sim.add(m=float(obj.params.stellar_mass))

    for j in used_indices:
        sim.add(
            m=float(obj.masses[j]) * mjup_to_msun,
            a=float(obj.semimajor[j]),
            e=float(obj.e[j]),
            inc=float(obj.i[j]) * deg_to_rad,
            Omega=float(obj.Node[j]) * deg_to_rad,
            omega=float(obj.w[j]) * deg_to_rad,
            M=float(obj.M0[j]) * deg_to_rad
            #primary=sim.particles[0],
            #jacobi_masses = True
        )
 

    if move_to_com:
        sim.move_to_com()

    sim.integrator = integrator
    sim.dt = timestep

    if integrator == "whfast":
        print("coordinates:",sim.integrator.coordinates)
    
        if hasattr(sim, "ri_whfast"):
            try:
                sim.ri_whfast.safe_mode = int(safe_mode)
            except Exception:
                pass

            try:
                sim.ri_whfast.corrector = int(corrector)
            except Exception:
                pass
        else:
            if verbose:
                print("Warning: this REBOUND version has no sim.ri_whfast.")
                print("Running WHFast with default WHFast settings.")

    if verbose:
        print("running REBOUND stability with integrator:", integrator)
        print("number of active planets:", npl)
        print("timemax [yr]:", timemax)
        print("timestep [yr]:", timestep)
        print("output_step [yr]:", output_step)

        
    times = np.arange(0.0, timemax + 0.5 * output_step, output_step)
    nt = len(times)

    evol = {}
    evol["stable"] = True
    evol["status"] = "finished"
    evol["used_indices"] = used_indices
    evol["integrator"] = integrator
    evol["timestep"] = timestep
    evol["output_step"] = output_step
    evol["timemax_requested"] = timemax

    evol["t"] = np.full(nt, np.nan)

    array_keys = [
        "a",
        "e",
        "inc",
        "Omega",
        "omega",
        "pomega",
        "M",
        "lambda",
        "P",
        "x",
        "y",
        "z",
        "vx",
        "vy",
        "vz",
    ]

    for key in array_keys:
        evol[key] = np.full((npl, nt), np.nan)

    evol["energy"] = np.full(nt, np.nan)
    evol["rel_energy_error"] = np.full(nt, np.nan)

    try:
        E0 = sim.energy()
    except Exception:
        E0 = np.nan

    def wrap_360(angle_deg):
        return angle_deg % 360.0

    def check_stop_conditions(sim):
        ps = sim.particles

        if escape_distance is not None:
            r_esc = float(escape_distance)

            for kk in range(1, npl + 1):
                r = np.sqrt(ps[kk].x**2 + ps[kk].y**2 + ps[kk].z**2)

                if r > r_esc:
                    return (
                        False,
                        "escape: particle %d reached r = %.8g au" % (kk, r),
                    )

        if collision_distance is not None:
            r_col = float(collision_distance)

            for kk in range(1, npl + 1):
                for ll in range(kk + 1, npl + 1):
                    dx = ps[kk].x - ps[ll].x
                    dy = ps[kk].y - ps[ll].y
                    dz = ps[kk].z - ps[ll].z

                    dr = np.sqrt(dx * dx + dy * dy + dz * dz)

                    if dr < r_col:
                        return (
                            False,
                            "close encounter: particles %d and %d reached dr = %.8g au"
                            % (kk, ll, dr),
                        )

        return True, "ok"

    last_good_index = -1

    for it, t in enumerate(times):

        try:
            sim.integrate(float(t), exact_finish_time=exact_finish_time)

        except Exception as err:
            evol["stable"] = False
            evol["status"] = (
                "REBOUND integration error at t = %.8g yr: %s"
                % (t, str(err))
            )
            break

        ok, msg = check_stop_conditions(sim)

        if not ok:
            evol["stable"] = False
            evol["status"] = msg + " at t = %.8g yr" % t
            break

        evol["t"][it] = sim.t

        try:
            E = sim.energy()
            evol["energy"][it] = E

            if np.isfinite(E0) and E0 != 0.0:
                evol["rel_energy_error"][it] = (E - E0) / E0
            else:
                evol["rel_energy_error"][it] = np.nan

        except Exception:
            evol["energy"][it] = np.nan
            evol["rel_energy_error"][it] = np.nan

        ps = sim.particles

        for k in range(npl):
            p = ps[k + 1]

            evol["x"][k, it] = p.x
            evol["y"][k, it] = p.y
            evol["z"][k, it] = p.z

            evol["vx"][k, it] = p.vx
            evol["vy"][k, it] = p.vy
            evol["vz"][k, it] = p.vz

            try:
                #orb = p.orbit(primary=ps[0])
                orb = p.orbit()
                
                evol["a"][k, it] = orb.a
                evol["e"][k, it] = orb.e

                evol["inc"][k, it] = wrap_360(orb.inc * rad_to_deg)
                evol["Omega"][k, it] = wrap_360(orb.Omega * rad_to_deg)
                evol["omega"][k, it] = wrap_360(orb.omega * rad_to_deg)
                evol["pomega"][k, it] = wrap_360(orb.pomega * rad_to_deg)
                evol["M"][k, it] = wrap_360(orb.M * rad_to_deg)
                evol["lambda"][k, it] = wrap_360(orb.l * rad_to_deg)

                mcentral = float(obj.params.stellar_mass)
                mplanet = float(obj.masses[used_indices[k]]) * mjup_to_msun

                if orb.a > 0.0:
                    evol["P"][k, it] = np.sqrt(
                        orb.a**3 / (mcentral + mplanet)
                    )
                else:
                    evol["P"][k, it] = np.nan

            except Exception:
                evol["a"][k, it] = np.nan
                evol["e"][k, it] = np.nan
                evol["inc"][k, it] = np.nan
                evol["Omega"][k, it] = np.nan
                evol["omega"][k, it] = np.nan
                evol["pomega"][k, it] = np.nan
                evol["M"][k, it] = np.nan
                evol["lambda"][k, it] = np.nan
                evol["P"][k, it] = np.nan

        last_good_index = it

    if last_good_index >= 0:
        n_keep = last_good_index + 1

        evol["t"] = evol["t"][:n_keep]
        evol["energy"] = evol["energy"][:n_keep]
        evol["rel_energy_error"] = evol["rel_energy_error"][:n_keep]

        for key in array_keys:
            evol[key] = evol[key][:, :n_keep]

    else:
        evol["t"] = np.array([])
        evol["energy"] = np.array([])
        evol["rel_energy_error"] = np.array([])

        for key in array_keys:
            evol[key] = np.empty((npl, 0))

    if len(evol["t"]) > 0:
        evol["timemax_reached"] = float(evol["t"][-1])
    else:
        evol["timemax_reached"] = np.nan

    obj.rebound_evol = evol
    obj.rebound_stable = evol["stable"]
    obj.rebound_status = evol["status"]
    obj.rebound_used_indices = used_indices
    from integration_period import record_period_source
    record_period_source(obj, used_indices)

    if not hasattr(obj, "evol_momentum") or obj.evol_momentum is None:
        obj.evol_momentum = {}

    obj.evol_T_energy = evol["t"]
    obj.evol_energy = evol["rel_energy_error"]

    obj.evol_momentum["lx"] = np.full_like(evol["t"], np.nan, dtype=float)
    obj.evol_momentum["ly"] = np.full_like(evol["t"], np.nan, dtype=float)
    obj.evol_momentum["lz"] = np.full_like(evol["t"], np.nan, dtype=float)

    obj.evol_T = {}
    obj.evol_a = {}
    obj.evol_e = {}
    obj.evol_i = {}
    obj.evol_Om = {}
    obj.evol_p = {}
    obj.evol_M = {}
    obj.evol_Per = {}

    for k in range(npl):
        obj.evol_T[k] = evol["t"]
        obj.evol_a[k] = evol["a"][k]
        obj.evol_e[k] = evol["e"][k]
        obj.evol_i[k] = evol["inc"][k]
        obj.evol_Om[k] = evol["Omega"][k]
        obj.evol_p[k] = evol["omega"][k]
        obj.evol_M[k] = evol["M"][k]

        # REBOUND period is in years because time unit is yr.
        # Convert to days for consistency with Exo-Striker.
        obj.evol_Per[k] = evol["P"][k] * 365.25

    if verbose:
        print("REBOUND stability done.")
        print("status:", evol["status"])
        print("stable:", evol["stable"])
        print("timemax reached [yr]:", evol["timemax_reached"])

    print(f"--- REBOUND {integrator} N-body: {perf_counter() - start_time:.5f} seconds ---")
    return obj



 





def run_stability_arb(obj, timemax=3000.0, timestep=10, timeout_sec=1000.0,
                      stab_save_dir='./', integrator=None, *, backend=None,
                      output_step=1.0, use_fixed_steps=False, fixed_steps=1000,
                      remove_stab_save_dir=True, **rebound_options):
    """Integrate arbitrary planets through the current SWIFT or REBOUND API.

    timestep is in days for either backend; timemax/output_step are years.
    Results populate the usual evol_* plot arrays, ordered by orbital period.
    Fitted planet inputs are untouched. An explicit integrator selects its
    backend; otherwise the selected backend's nbody_* settings supply it.
    """
    from types import SimpleNamespace

    swift_integrators = ('symba', 'mvs', 'mvs_gr')
    rebound_integrators = ('whfast', 'mercurius', 'saba', 'ias15', 'trace', 'janus')
    if integrator is not None:
        integrator = str(integrator).lower()
    if backend is None:
        backend = 'rebound' if integrator in rebound_integrators else 'swift'
    backend = str(backend).lower()
    if backend not in ('swift', 'rebound'):
        raise ValueError("backend must be swift or rebound")
    if integrator is None:
        settings = getattr(obj, 'nbody_' + backend, None)
        integrator = getattr(settings, 'integrator',
                             'symba' if backend == 'swift' else 'whfast')
        integrator = str(integrator).lower()
    allowed = swift_integrators if backend == 'swift' else rebound_integrators
    if integrator not in allowed:
        raise ValueError(f"Unsupported {backend} integrator: {integrator}")
    if backend == 'swift' and rebound_options:
        raise ValueError("REBOUND options require backend='rebound'")
    enabled = [j for j in range(9) if obj.pl_arb_use[j]]
    if not enabled or len(enabled) != obj.npl_arb:
        raise ValueError("obj.npl_arb must match the enabled arbitrary planets and be positive")

    # Adapt arbitrary inputs to the interface used by fitted-orbit runs.
    # Separate output dictionaries preserve previous results on failure.
    model = SimpleNamespace(
        params=SimpleNamespace(stellar_mass=obj.arb_st_mass),
        use_planet=dict(obj.pl_arb_use), npl=obj.npl_arb,
        masses=dict(obj.mass_arb), semimajor=dict(obj.a_arb),
        e=dict(obj.e_arb), i=dict(obj.i_arb), w=dict(obj.w_arb),
        Node=dict(obj.Node_arb), M0=dict(obj.M0_arb),
        P={j: a_to_P(obj.a_arb[j], obj.arb_st_mass) for j in enabled},
        GR_step=getattr(obj, 'GR_step', 1),
        evol_momentum={'lx': [], 'ly': [], 'lz': []})
    outputs = ('evol_T', 'evol_a', 'evol_e', 'evol_i', 'evol_Om',
               'evol_p', 'evol_M', 'evol_Per')
    for attribute in outputs:
        setattr(model, attribute, {})
    options = dict(timemax=timemax, integrator=integrator,
                   output_step=output_step, use_fixed_steps=use_fixed_steps,
                   fixed_steps=fixed_steps)
    if backend == 'swift':
        model = run_stability(model, timestep=timestep, timeout_sec=timeout_sec,
                              stab_save_dir=stab_save_dir,
                              remove_stab_save_dir=remove_stab_save_dir, **options)
    else:
        model = run_stability_rebound(model, timestep=timestep / 365.25,
                                      **options, **rebound_options)
    for attribute in (*outputs, 'evol_T_energy', 'evol_energy', 'evol_momentum'):
        setattr(obj, attribute, getattr(model, attribute))
    if backend == 'rebound':
        for attribute in ('rebound_evol', 'rebound_stable', 'rebound_status',
                          'rebound_used_indices'):
            setattr(obj, attribute, getattr(model, attribute))
    else:
        obj.swift_evol = model.swift_evol
    # These arrays now describe arbitrary planets, not the fitted RV system.
    obj.phase_period_source = None
    return obj


def update_jacobi_masses(obj):
    """Refresh physical derived values for active GUI slots, including old sessions."""
    slots = [i for i, active in enumerate(obj.use_planet) if active]
    masses = np.zeros(len(obj.use_planet))
    axes = np.zeros(len(obj.use_planet))
    rows = [[obj.K[i], obj.P[i], obj.e[i], obj.w[i], obj.M0[i], obj.i[i], obj.Node[i]]
            for i in slots]
    physical, jacobi_axes = mass_a_from_jacobi_el(rows, len(slots), obj.stellar_mass)
    masses[slots], axes[slots] = physical, jacobi_axes
    obj.masses, obj.semimajor = masses.tolist(), axes.tolist()
    # Keep consumers of saved fit results consistent with the GUI and exports.
    for name in ('fit_results', 'rvmod'):
        result = getattr(obj, name, None)
        if result is not None:
            result.mass, result.a = masses.copy(), axes.copy()
    return masses, axes
