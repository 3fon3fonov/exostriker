 #!/usr/bin/python3

__author__ = 'Trifon Trifonov'

import sys, os
#sys.path.insert(0, '../lib')
import numpy as np
import jac2astrocen
import corner_ES as corner

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

matplotlib.rcParams['axes.formatter.useoffset'] = False

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
        fit_new.tra_jitt_use  = {k: False for k in range(20)}
        fit_new.tra_off_use  = {k: False for k in range(20)}
        fit_new.tra_dil_use  = {k: False for k in range(20)}

        fit_new.init_tra_jitter()
        fit_new.init_tra_offset()
        fit_new.init_tra_dilution()
        fit_new.init_tra_lintr()
        fit_new.init_tra_quadtr()


    fit_new.cwd = dill.copy(fit.cwd)

    return fit_new






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


def transit_tperi(per, ecc, om, ma, epoch, primary = True):
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


def ma_from_epoch(per, t_peri, epoch):
    '''
    '''
    ma =  np.degrees(2.0*np.pi*( (epoch-t_peri)/per % 1.))

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

            tr_ind = np.where(np.logical_and(t >= tran_t0-0.07, t <= tran_t0+0.07))
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
    '''Calculates the actual masses and Jacobi semimajor axes of a
       system for using the parameters P, K and e from a Kepler fit.
       The output of the semi-major axis is in AU.
       if mass_type="J" the mass is in Jupiter masses (deafault)
       if mass_type="E" the mass is in Erath masses (deafault)
       else, e.g.,        mass_type="" mass is in Solar units.
    '''
    THIRD = 1.0/3.0
    PI    = 3.14159265358979e0
    TWOPI = 2.0*PI
    GMSUN = 1.32712440018e20
    AU=1.49597892e11

    mass = np.zeros(10)
    ap    = np.zeros(9)
    pl_mass = np.zeros(9)
    mpold = pl_mass

#*******G is set to be unit, and s, m, kg as unit of time, length and mass
#*****  and there is a reason for that! later I might correct for that.
    mtotal = m_s
    f = 5e-6

    #mass[0] = m_s


    for i in range(len(K)):

        T = P[i]*86400.0

        # we need innitial guess for each planet mass
        dm =2
        mass[i+1] = abs(K[i])*(T*m_s**2.0/(TWOPI*GMSUN))**THIRD * np.sqrt(1.0-ecc[i]**2.0)/abs(np.sin(np.radians(incl[i])))
        #mpold[i] = mass[i+1]
        # This is a simple iteration to solve for mp
        mass[0] = m_s
        mpold[i] = 0

        while (dm >= f):

            if i == 0:
                mtotal = m_s
                mass[i+1] = abs(K[i])*(T*(m_s + mpold[i])**2.0/(TWOPI*GMSUN))**THIRD * np.sqrt(1.0-ecc[i]**2.0)/abs(np.sin(np.radians(incl[i])))
            else:
                mtotal = m_s
                for j in range(i):
                    mtotal = mtotal + mass[j+1]
                mass[i+1] = abs(K[i])*(T*(mtotal + mpold[i])**2.0/(TWOPI*GMSUN))**THIRD * np.sqrt(1.0-ecc[i]**2.0)/abs(np.sin(np.radians(incl[i])))

            #dm = (mpold[i] - mass[i+1])
            dm = abs((mass[i+1] - mpold[i])/mass[i+1] )

            mpold[i] =  mass[i+1]


        ap[i] = (GMSUN*(mtotal + mass[i+1])*(T/TWOPI)**2)**THIRD

    for i in range(len(K)):

        ap[i] = ap[i]/AU # to be in AU
        if mass_type=="J":
            pl_mass[i] = mass[i+1]*1047.5654817267318 # to be in Jup. masses
        elif  mass_type=="E":
            pl_mass[i] = mass[i+1]*1047.5654817267318 * 317.82838 # to be in Earth. masses
        else:
            pl_mass[i] = mass[i+1]


        # I have seen that 1 Sol Mass = 1047.92612 Jup. masses???
    return pl_mass,ap

####################### mass_semimajor #########################################
def get_mass_a(obj, mass_type="J"):
    '''Calculates the actual masses and Jacobi semimajor axes of a
       system for using the parameters P, K and e from a Kepler fit.
       The output of the semi-major axis is in AU.
       if mass_type="J" the mass is in Jupiter masses (deafault)
       if mass_type="E" the mass is in Erath masses (deafault)
       else, e.g.,        mass_type="" mass is in Solar units.
    '''
    THIRD = 1.0/3.0
    PI    = 3.14159265358979e0
    TWOPI = 2.0*PI
    GMSUN = 1.32712440018e20
    AU=1.49597892e11

    mass = np.zeros(10)
    ap    = np.zeros(9)
    pl_mass = np.zeros(9)
    mpold = pl_mass

#*******G is set to be unit, and s, m, kg as unit of time, length and mass
#*****  and there is a reason for that! later I might correct for that.
    mtotal = obj.st_mass[0]
    f = 5e-6
    for i in range(obj.npl):

        T = obj.P[i]*86400.0
        mass[0] = obj.st_mass[0]

        # we need innitial guess for each planet mass
        dm = 0
        mass[i+1] = abs(obj.K[i])*(T*(obj.st_mass[0])**2.0/(TWOPI*GMSUN))**THIRD * np.sqrt(1.0-obj.e[i]**2.0)/abs(np.sin(np.radians(obj.i[i])))
        mpold[i] = mass[i+1]
        # This is a simple iteration to solve for mp
        while (dm <= 0):

            if i == 0:
                mtotal = obj.st_mass[0]
                mass[i+1] = abs(obj.K[i])*(T*(obj.st_mass[0] + mpold[i])**2.0/(TWOPI*GMSUN))**THIRD * np.sqrt(1.0-obj.e[i]**2.0)/abs(np.sin(np.radians(obj.i[i])))
                mtotal = obj.st_mass[0]
                for j in range(i):
                    mtotal = mtotal + mass[j+1]
                mass[i+1] = abs(obj.K[i])*(T*(mtotal + mpold[i])**2.0/(TWOPI*GMSUN))**THIRD * np.sqrt(1.0-obj.e[i]**2.0)/abs(np.sin(np.radians(obj.i[i])))

            dm = (mpold[i] - mass[i+1])
            mpold[i] =  mpold[i] + f
           # print mass[i+1], mpold[i]

        ap[i] = (GMSUN*(mtotal + mass[i+1])*(T/TWOPI)**2)**THIRD

#    for i in range(npl+1):
#        mass[i] = mass[i]*GMSUN
    for i in range(obj.npl):

        ap[i] = ap[i]/AU # to be in AU
        if mass_type=="J":
            pl_mass[i] = mass[i+1]*1047.5654817267318 # to be in Jup. masses
        elif  mass_type=="E":
            pl_mass[i] = mass[i+1]*1047.5654817267318 * 317.82838 # to be in Earth. masses
        else:
            pl_mass[i] = mass[i+1]


        # I have seen that 1 Sol Mass = 1047.92612 Jup. masses???
    return pl_mass,ap



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

def add_mcmc_samples(obj,sampler):

    bestfit_labels      = ["median","mean","mode","best_samp","best_gui","none",
                           "mass","use_Me","use_Mj","use_Ms",
                           "semimajor","radius","use_Re","use_Rj","use_Rs","use_lambda", "use_ppm"]
    bestfit_labels_bool = [obj.mcmc_save_median,obj.mcmc_save_means,obj.mcmc_save_mode,
                           obj.mcmc_save_maxlnL,False,False,False,
                           True,False,False,False,False,True,False,False,False,False]


    sampler.lbf             = {k: np.array([obj.e_for_mcmc[k], True]) for k in range(len(obj.e_for_mcmc))}
    for k in range(17):
        sampler.lbf[bestfit_labels[k]] = bestfit_labels_bool[k]

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

    sampler.lbf["cornerplot"] = cornerplot_opt
    sampler.lbf["OrigLabels"] = dill.copy(obj.e_for_mcmc)

    obj.mcmc_sampler=sampler
    obj.sampler_saved=True




def add_ns_samples(obj,sampler):

    bestfit_labels      = ["median","mean","mode","best_samp","best_gui","none",
                           "mass","use_Me","use_Mj","use_Ms",
                           "semimajor","radius","use_Re","use_Rj","use_Rs","use_lambda","use_ppm"]
    bestfit_labels_bool = [obj.ns_save_median,obj.ns_save_means,obj.ns_save_mode,
                           obj.ns_save_maxlnL,False,False,False,
                           True,False,False,False,False,True,False,False,False,False]



    obj.ns_sampler= dill.copy(sampler.results)
    obj.ns_sampler.lbf     = {k: np.array([obj.e_for_mcmc[k], True]) for k in range(len(obj.e_for_mcmc))}

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




def cornerplot(obj, level=(100.0-68.3)/2.0, type_plot = 'mcmc', **kwargs):

    #obj = dill.copy(copied_obj)
    '''Generates a corner plot visualizing the mcmc samples. Optionally samples can be read from a file.'''
    #self.mcmc_sample_file = 'mcmc_samples'+'_%s'%mod
    #self.corner_plot_file = 'cornerplot.png'


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
    #r_s   = np.random.normal(loc=obj.stellar_radius,    scale=obj.stellar_radius_err,    size=len(samples[:,0]))
   # L_s   = np.random.normal(loc=obj.stellar_luminosity,scale=obj.stellar_luminosity_err,size=len(samples[:,0]))
   # vsini = np.random.normal(loc=obj.stellar_vsini,     scale=obj.stellar_vsini_err,     size=len(samples[:,0]))
   #

    ######### define new samples, labels and best-fit params to be refilled again
    ######### with masses, semi-major axes, etc. (to be re-worked).

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


    samp              = [samp[i] for i in index_to_keep]
    samp_labels       = [samp_labels[i] for i in index_to_keep]
    samp_best_fit_par = [samp_best_fit_par[i] for i in index_to_keep]

    letters = ['b','c','d','e','f','g','h'] #... For the planets

    if mod_labels['use_lambda']:

        for i in range(obj.npl):
            let = letters[i]

            if not 'e$_%s$'%let in labels or not '$\omega_%s$'%let in labels:
                continue

            Ma_    = np.hstack(samples[:,[ii for ii, j in enumerate(labels) if j == 'MA$_%s$'%let]])
            omega_ = np.hstack(samples[:,[ii for ii, j in enumerate(labels) if j == '$\omega_%s$'%let]])

            lambda_ = np.array(Ma_ + omega_)%360.0

            samp.append(lambda_)
            samp_labels.append(r' $\lambda_%s$'%let)


            if mod_labels['mean']:
                samp_best_fit_par.append((np.mean(Ma_) + np.mean(omega_))%360.0 )
            elif mod_labels['median']:
                samp_best_fit_par.append((np.median(Ma_) + np.median(omega_))%360.0 )
            elif mod_labels['best_gui']:
                samp_best_fit_par.append((Ma_[np.argmax(ln)] + omega_[np.argmax(ln)])%360.0 )
            else:
                samp_best_fit_par.append((Ma_[np.argmax(ln)] + omega_[np.argmax(ln)])%360.0 )


    if mod_labels['mass'] or mod_labels['semimajor']:

        m_s   = np.random.normal(loc=obj.stellar_mass,      scale=obj.stellar_mass_err,      size=len(ss))
        K,P, ecc, esinw, ecosw, incl,masses,semimajor = [],[],[],[],[],[],[],[]

        for i in range(9):
            if not bool(obj.use_planet[i]):
                continue
            let = letters[i]


            if not 'K$_%s$'%let in labels or not 'P$_%s$'%let in labels:
                continue

            K.append(np.hstack(samples[:,[ii for ii, j in enumerate(labels) if j == 'K$_%s$'%let]]))
            P.append(np.hstack(samples[:,[ii for ii, j in enumerate(labels) if j == 'P$_%s$'%let]]))

            if obj.hkl == True and '$e sin(\omega_%s)$'%let in labels and '$e cos(\omega_%s)$'%let in labels:

                esinw.append(np.hstack(samples[:,[ii for ii, j in enumerate(labels) if j == '$e sin(\omega_%s)$'%let]]))
                ecosw.append(np.hstack(samples[:,[ii for ii, j in enumerate(labels) if j == '$e cos(\omega_%s)$'%let]]))

                ecc.append(np.sqrt(esinw[i]**2 + ecosw[i]**2))
            elif obj.hkl == False and 'e$_%s$'%let in labels:
                ecc.append(np.hstack(samples[:,[ii for ii, j in enumerate(labels) if j == 'e$_%s$'%let]]))

            else:
                ecc.append([0]*len(K[i]))
                print("Warning, no eccentricity samples found for planet %s ! Assuming ecc = 0"%str(i+1))

            if mod_labels['use_Me']:
                M_fact = 317.82838
                mass_lab = r'[M$_\oplus$]'
            elif mod_labels['use_Mj']:
                M_fact = 1
                mass_lab = r'[M$_{\rm Jup.}$]'
            elif mod_labels['use_Ms']:
                M_fact = 1.0/1047.5654817267318
                mass_lab = r'[M$_\odot$]'


            if 'i$_%s$'%let in labels:
                incl.append(np.hstack(samples[:,[ii for ii, j in enumerate(labels) if j == 'i$_%s$'%let]]))
                samp_labels.append(r'm$_%s$ %s'%(let,mass_lab))
            else:
                if obj.copl_incl == True:
                    incl.append(incl[0])
                else:
                    incl.append([obj.i[i]]*len(K[i]))
                if obj.i[i] == 90.0:
                    samp_labels.append(r'm $\sin i_%s$ %s'%(let,mass_lab))
                else:
                    samp_labels.append(r'm$_%s$ %s'%(let,mass_lab))

        K = np.transpose(K)
        P = np.transpose(P)
        ecc = np.transpose(ecc)
        incl = np.transpose(incl)

        for k in range(len(m_s)):

            m_a = get_mass_a_samples(np.array(K[k]),np.array(P[k]), np.array(ecc[k]), np.array(incl[k]), m_s[k], mass_type="J")
            masses.append(m_a[0])
            semimajor.append(m_a[1])

        masses = np.transpose(masses)
        semimajor = np.transpose(semimajor)


        if mod_labels['mass']:
            for i in range(9):
                if not bool(obj.use_planet[i]):
                    continue
                samp.append(np.array(masses[i] * M_fact))

                if mod_labels['mean']:
                    samp_best_fit_par.append( np.mean(masses[i]) * M_fact)
                elif mod_labels['median']:
                    samp_best_fit_par.append( np.median(masses[i]) * M_fact)
                elif mod_labels['best_gui']:
                    samp_best_fit_par.append(obj.masses[i]*M_fact)
                else:
                    samp_best_fit_par.append(masses[i][np.argmax(ln)]*M_fact)


        if mod_labels['semimajor']:

            for i in range(9):
                if not bool(obj.use_planet[i]):
                    continue
                samp.append(np.array(semimajor[i]))
                let = letters[i]
                samp_labels.append(r'a$_%s$ [au]'%let)

                if mod_labels['mean']:
                    samp_best_fit_par.append(np.mean(semimajor[i]))
                elif mod_labels['median']:
                    samp_best_fit_par.append(np.median(semimajor[i]))
                elif mod_labels['best_gui']:
                    samp_best_fit_par.append(obj.semimajor[i])
                else:
                    samp_best_fit_par.append(semimajor[i][np.argmax(ln)])





    if mod_labels['radius']:
        r_s   = np.random.normal(loc=obj.stellar_radius,    scale=obj.stellar_radius_err,    size=len(samples[:,0]))

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
#            r_s   = np.random.normal(loc=obj.stellar_radius,    scale=obj.stellar_radius_err,    size=len(samples[:,0]))


    ################### Transpose is needed for the cornerplot. ###################

    samples_ = np.transpose(samp)
    #labels = samp_labels
    best_fit_par =samp_best_fit_par

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
 


    if verbose:
        print("   ")
        print("   ")
        if mod_labels['mean']:
            print("Means and their 1 sigma errors")
        elif mod_labels['median']:
            print("Median and their 1 sigma errors")
        elif mod_labels['best_gui']:
            print("Best-fit (GUI) and their 1 sigma errors")
        else:
            print("Best-fit (max. samp. lnL) and their 1 sigma errors")


        for i in range(len(samp_best_fit_par)):
            ci = np.percentile(samp[i], [level, 100.0-level])
            if mod_labels['mean']:
                print(samp_labels[i],'=', np.mean(samp[i]), "- %s"%(np.mean(samp[i])-ci[0]), "+ %s"%(ci[1]  - np.mean(samp[i]) ))
            elif mod_labels['median']:
                print(samp_labels[i],'=', np.median(samp[i]), "- %s"%(np.median(samp[i])-ci[0]), "+ %s"%(ci[1]  - np.median(samp[i]) ))
            elif mod_labels['best_gui']:
                print(samp_labels[i],'=', samp_best_fit_par[i], "- %s"%(samp_best_fit_par[i]-ci[0]), "+ %s"%(ci[1]  - samp_best_fit_par[i] ))
            else:
                print(samp_labels[i],'=', samp[i][np.argmax(ln)], "- %s"%(samp[i][np.argmax(ln)]-ci[0]), "+ %s"%(ci[1]  - samp[i][np.argmax(ln)] ))

        print(" ")
        print("Median Absolute Deviation values")

        mad = get_MAD_of_samples(samples_,len(samp_labels))
        for i in range(len(samp_labels)):
            print(samp_labels[i],' MAD =', mad[i])


    range_stab =[]
#    for i in range(len(samp)):
#        range_stab.append([0.0,1.0])

    if mod_labels['none']==True:
        best_fit_par = None
        
        



    if "max. time" in labels:
        samples_stab = []

        samples_ = np.transpose(samples_)
        for i in range(len(samples_)):
            samples_stab.append(samples_[i][np.where(samples[:,-1]> cornerplot_opt["stab_threshold"])])
#        print(max(samples[:,-1])-1000.0,max(samples[:,-1]))
            range_stab.append([min(samples_[i]),max(samples_[i])])
        N_samp = len(np.atleast_1d(samples_stab[0]))



        if verbose:
            print("   ")
            print("   ")
            if mod_labels['mean']:
                print("Means and their 1 sigma errors (stable)")
            elif mod_labels['median']:
                print("Median and their 1 sigma errors (stable)")
            elif mod_labels['best_gui']:
                print("Best-fit (GUI) and their 1 sigma errors (stable)")
            else:
                print("Best-fit (max. samp. lnL) and their 1 sigma errors (stable)")


            for i in range(len(samp_best_fit_par)):
                ci = np.percentile(samples_stab[i], [level, 100.0-level])
                if mod_labels['mean']:
                    print(samp_labels[i],'=', np.mean(samples_stab[i]), "- %s"%(np.mean(samples_stab[i])-ci[0]), "+ %s"%(ci[1]  - np.mean(samples_stab[i]) ))
                elif mod_labels['median']:
                    print(samp_labels[i],'=', np.median(samples_stab[i]), "- %s"%(np.median(samples_stab[i])-ci[0]), "+ %s"%(ci[1]  - np.median(samples_stab[i]) ))
                elif mod_labels['best_gui']:
                    print(samp_labels[i],'=', samp_best_fit_par[i], "- %s"%(samp_best_fit_par[i]-ci[0]), "+ %s"%(ci[1]  - samp_best_fit_par[i] ))
                else:
                    print(samp_labels[i],'=', samples_stab[i][np.argmax(ln)], "- %s"%(samples_stab[i][np.argmax(ln)]-ci[0]), "+ %s"%(ci[1]  - samples_stab[i][np.argmax(ln)] ))

            print(" ")

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
            print("very few, or no stable samples !!!!!")
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

    #axes = fig.axes()
  #  axes.ticklabel_format(style='plain')

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





def phase_RV_planet_signal(obj,planet):

    if obj.npl ==0 or len(obj.fit_results.rv_model.jd) ==0:
        return
    else:

        copied_obj = dill.copy(obj)

        if(copied_obj.mod_dynamical):
            copied_obj.mod_dynamical = False

        index = planet - 1
        #print(index)
 


        ############ phase fold fix for sparse model ######
        model_time_phase = np.array( (copied_obj.fit_results.model_jd -copied_obj.fit_results.model_jd[0] + (copied_obj.fit_results.model_jd[0] - copied_obj.epoch) )%copied_obj.P[index] )

        model_shift = copied_obj.P[index] - (copied_obj.fit_results.rv_model.jd[0] - copied_obj.epoch )%copied_obj.P[index]

        model_time_phase = (model_time_phase + model_shift)% copied_obj.P[index]

        sort = sorted(range(len(model_time_phase)), key=lambda k: model_time_phase[k])
        model_time_phase  = model_time_phase[sort]

        phased_model      = copied_obj.fit_results.res['fit'].T[2+index][sort]

        ############ phase data ######
        data_time_phase = np.array( (copied_obj.fit_results.rv_model.jd  - copied_obj.fit_results.rv_model.jd[0])% copied_obj.P[index] )

        sort = sorted(range(len(data_time_phase)), key=lambda k: data_time_phase[k])
        data_time_phase      = data_time_phase[sort]
        phased_data          = copied_obj.fit_results.model_data[6+index][sort] + copied_obj.fit_results.rv_model.o_c[sort]  
        
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





def find_planets(obj,fend=0.75):

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

        run_gls_o_c(obj,fend=fend)

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
                    obj = run_gls_o_c(obj)
                return obj

            else:
                mean_anomaly_from_gls = np.degrees((((obj.epoch - float(obj.gls_o_c.hpstat["T0"]) )% (obj.gls_o_c.hpstat["P"]) )/ (obj.gls_o_c.hpstat["P"]) ) * 2*np.pi)

                #obj.add_planet(obj.gls_o_c.hpstat["amp"],obj.gls_o_c.hpstat["P"],0.0,0.0,mean_anomaly_from_gls -90.0,90.0,0.0)

                obj.add_planet(obj.gls.hpstat["amp"],obj.gls.hpstat["P"],0.0,0.0,mean_anomaly_from_gls -90.0,90.0,0.0, useK=True,useP=True,usee=obj.auto_fit_allow_ecc,usew=obj.auto_fit_allow_ecc,useM0=True,usei=False, usecap=False)

                obj.fitting(fileinput=False,outputfiles=[1,1,1], doGP=False,   minimize_fortran=True, fortran_kill=3, timeout_sec= 3)
                run_gls_o_c(obj,fend=fend)

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
        run_gls_o_c(obj,fend=fend)
    return obj





def run_gls(obj,fend =1.0,fbeg=10000):


    #fbeg = abs(max(obj.fit_results.rv_model.jd)-min(obj.fit_results.rv_model.jd))  * 2.0
    omega = 1/ np.logspace(np.log10(fend), np.log10(fbeg), num=int(1000))



    if len(obj.fit_results.rv_model.jd) > 5:
        RV_per = gls.Gls((obj.fit_results.rv_model.jd, obj.fit_results.rv_model.rvs, obj.fit_results.rv_model.rv_err),
        fast=True,  verbose=False, norm='ZK',ofac=10, fbeg=omega[-1], fend=omega[0],)

        obj.gls = RV_per
    else:
        return obj

    return obj

def run_gls_o_c(obj,fend =1.0,fbeg=10000, as_main = False):

    #fbeg = abs(max(obj.fit_results.rv_model.jd)-min(obj.fit_results.rv_model.jd))  * 2.0
    omega = 1/ np.logspace(np.log10(fend), np.log10(fbeg), num=int(1000))


    if len(obj.fit_results.rv_model.jd) > 5:
        RV_per_res = gls.Gls((obj.fit_results.rv_model.jd, obj.fit_results.rv_model.o_c, obj.fit_results.rv_model.rv_err),
        fast=True,  verbose=False, norm='ZK', ofac=10, fbeg=omega[-1], fend=omega[0],)

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





def latex_pl_param_table(obj, width = 10, precision = 2, asymmetric = False, file_name='test.tex', path='./', return_text=False):


    if asymmetric != True:

        text = '''
    \\begin{table}[ht]
    % \\begin{adjustwidth}{-4.0cm}{}
    % \\resizebox{0.69\\textheight}{!}
    % {\\begin{minipage}{1.1\\textwidth}

    \centering
    \caption{{}}
    \label{table:}

    \\begin{tabular}{lrrrrrrrr}     % 2 columns

    \hline\hline  \\noalign{\\vskip 0.7mm}
    '''


        text = text + '''Parameter \hspace{0.0 mm}'''
        #for i in range(obj.npl):

        for i in range(9):
            if not bool(obj.use_planet[i]):
                continue
            text = text + '''& Planet %s '''%chr(98+i)
        text = text + '''\\\\
    \hline \\noalign{\\vskip 0.7mm}

        '''
        if obj.type_fit["RV"]== True or obj.type_fit["TTV"]== True:
            text = text + '''{0:{width}s}'''.format("$K$ [m\,s$^{-1}$]", width = 30)
            for i in range(9):
                if not bool(obj.use_planet[i]):
                    continue
                text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.K[i], max(np.abs(obj.param_errors.planet_params_errors[7*i])), width = width, precision = precision)
            text = text + '''\\\\
        '''

        text = text + '''{0:{width}s}'''.format("$P$ [day]", width = 30)
        for i in range(9):
            if not bool(obj.use_planet[i]):
                continue
            text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.P[i], max(np.abs(obj.param_errors.planet_params_errors[7*i +1])), width = width, precision = precision)
        text = text + '''\\\\
        '''

        if obj.hkl == False:
            text = text + '''{0:{width}s}'''.format("$e$  ", width = 30)
            for i in range(9):
                if not bool(obj.use_planet[i]):
                    continue
                text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.e[i], max(np.abs(obj.param_errors.planet_params_errors[7*i +2])), width = width, precision = precision)
            text = text + '''\\\\
        '''
            text = text + '''{0:{width}s}'''.format("$\omega$ [deg]", width = 30)
            for i in range(9):
                if not bool(obj.use_planet[i]):
                    continue
                text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.w[i], max(np.abs(obj.param_errors.planet_params_errors[7*i +3])), width = width, precision = precision)
            text = text + '''\\\\
        '''
        else:
            text = text + '''{0:{width}s}'''.format("$e sin(\omega)$  ", width = 30)
            for i in range(9):
                if not bool(obj.use_planet[i]):
                    continue
                text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.e_sinw[i], max(np.abs(obj.e_sinw_err[i])), width = width, precision = precision)
            text = text + '''\\\\
        '''
            text = text + '''{0:{width}s}'''.format("$e cos(\omega)$", width = 30)
            for i in range(9):
                if not bool(obj.use_planet[i]):
                    continue
                text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.e_cosw[i], max(np.abs(obj.e_cosw_err[i])), width = width, precision = precision)
            text = text + '''\\\\
        '''

        if obj.type_fit["RV"]== True or obj.type_fit["TTV"]== True:

            if obj.hkl == False:
                text = text + '''{0:{width}s}'''.format("$M_{\\rm 0}$ [deg]", width = 30)
                for i in range(9):
                    if not bool(obj.use_planet[i]):
                        continue
                    text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.M0[i], max(np.abs(obj.param_errors.planet_params_errors[7*i +4])), width = width, precision = precision)
                text = text + '''\\\\
        '''
            else:
                text = text + '''{0:{width}s}'''.format("$\lambda$ [deg]", width = 30)
                for i in range(9):
                    if not bool(obj.use_planet[i]):
                        continue
                    text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.lamb[i], max(np.abs(obj.lamb_err[i])), width = width, precision = precision)
                text = text + '''\\\\
        '''
        if obj.mod_dynamical == True:
            text = text + '''{0:{width}s}'''.format("$i$ [deg]", width = 30)
            for i in range(9):
                if not bool(obj.use_planet[i]):
                    continue
                text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.i[i], max(np.abs(obj.param_errors.planet_params_errors[7*i +5])), width = width, precision = precision)
            text = text + '''\\\\
        '''
            text = text + '''{0:{width}s}'''.format("$\Omega$ [deg]", width = 30)
            for i in range(9):
                if not bool(obj.use_planet[i]):
                    continue
                text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.Node[i], max(np.abs(obj.param_errors.planet_params_errors[7*i +6])), width = width, precision = precision)
            text = text + '''\\\\
        '''

        if obj.type_fit["Transit"] == True:
            text = text + '''{0:{width}s}'''.format("$t_{\\rm 0}$ [day]", width = 30)
            for i in range(9):
                if not bool(obj.use_planet[i]):
                    continue
                text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.t0[i], max(np.abs(obj.t0_err[i])), width = width, precision = precision)
            text = text + '''\\\\
        '''
            text = text + '''{0:{width}s}'''.format("Rad. [$R_\odot$]", width = 30)
            for i in range(9):
                if not bool(obj.use_planet[i]):
                    continue
                text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.pl_rad[i], max(np.abs(obj.pl_rad_err[i])), width = width, precision = precision)
            text = text + '''\\\\
        '''
            text = text + '''{0:{width}s}'''.format("$a$ [$R_\odot$]", width = 30)
            for i in range(9):
                if not bool(obj.use_planet[i]):
                    continue
                text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.pl_a[i], max(np.abs(obj.pl_a_err[i])), width = width, precision = precision)
            text = text + '''\\\\
        '''

        text = text + '''{0:{width}s}'''.format("$a$ [au]", width = 30)
        for i in range(9):
            if not bool(obj.use_planet[i]):
                continue
            text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.fit_results.a[i], 0, width = width, precision = precision)
        text = text + '''\\\\
        '''
        text = text + '''{0:{width}s}'''.format("$m \sin i$ [$M_{\\rm jup}$]", width = 30)
        for i in range(9):
            if not bool(obj.use_planet[i]):
                continue
            text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(obj.fit_results.mass[i], 0, width = width, precision = precision)
        text = text + '''\\\\
        '''
        text = text + '''{0:{width}s}'''.format("$t_{\omega}$ [day]", width = 30)
        for i in range(9):
            if not bool(obj.use_planet[i]):
                continue
            text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format((float(obj.epoch) - (np.radians(obj.M0[i])/(2*np.pi))*obj.P[i] ), 0, width = width, precision = precision)
        text = text + '''\\\\
        '''

        if obj.type_fit["RV"]== True:


            text = text + '''{0:{width}s}'''.format("RV lin. trend [m\,s$^{-1}$\,day$^{-1}$]", width = 30)
            text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(float(obj.rv_lintr),float(max(np.abs(obj.rv_lintr_err))) , width = 30, precision = 6)
            text = text + '''\\\\
        '''

            text = text + '''{0:{width}s}'''.format("RV quad. trend [m\,s$^{-1}$\,day$^{-1}$]", width = 30)
            text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(float(obj.rv_quadtr),float(max(np.abs(obj.rv_quadtr_err))) , width = 30, precision = 6)
            text = text + '''\\\\
        '''

            for i in range(obj.ndset):
                text = text + '''{0:{width}s}'''.format("RV$_{\\rm off}$ %s [m\,s$^{-1}$]"%(i+1), width = 30)
                text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(float(obj.rvoff[i]), float(max(np.abs(obj.rvoff_err[i]))), width = width, precision = precision)
                text = text + '''\\\\
        '''
            for i in range(obj.ndset):
                text = text + '''{0:{width}s}'''.format("RV$_{\\rm jit}$ %s [m\,s$^{-1}$]"%(i+1), width = 30)
                text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(float(obj.jitt[i]), float(max(np.abs(obj.jitt_err[i]))), width = width, precision = precision)
                text = text + '''\\\\
        '''

        if obj.type_fit["Transit"]== True:

            for i in range(10):
                if len(obj.tra_data_sets[i]) != 0:
                    text = text + '''{0:{width}s}'''.format("Tran.$_{\\rm off}$ %s"%(i+1), width = 30)
                    text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(float(obj.tra_off[i]), float(max(np.abs(obj.tra_off_err[i]))), width = width, precision = precision)
                    text = text + '''\\\\ \\noalign{\\vskip 0.9mm}
        '''
            for i in range(10):
                if len(obj.tra_data_sets[i]) != 0:
                    text = text + '''{0:{width}s}'''.format("Tran.$_{\\rm jit}$ %s"%(i+1), width = 30)
                    text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(float(obj.tra_jitt[i]), float(max(np.abs(obj.tra_jitt_err[i]))), width = width, precision = precision)
                    text = text + '''\\\\ \\noalign{\\vskip 0.9mm}
        '''
        if obj.doGP == True:

            if obj.gp_kernel == 'RotKernel':
                for i in range(4):
                    text = text + '''{0:{width}s}'''.format("%s"%(obj.GP_rot_str[i]), width = 30)
                    text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(float(obj.GP_rot_params[i]), float(max(np.abs(obj.param_errors.GP_params_errors[i]))), width = width, precision = precision)
                    text = text + '''\\\\ \\noalign{\\vskip 0.9mm}
        '''


            elif obj.gp_kernel == 'SHOKernel':
                for i in range(3):
                    text = text + '''{0:{width}s}'''.format("%s"%(obj.GP_sho_str[i]), width = 30)
                    text = text + '''& {0:{width}.{precision}f} $\pm$ {1:{width}.{precision}f} '''.format(float(obj.GP_sho_params[i]), float(max(np.abs(obj.param_errors.GP_params_errors[i]))), width = width, precision = precision)
                    text = text + '''\\\\ \\noalign{\\vskip 0.9mm}
        '''



        text = text + '''{0:{width}s}'''.format("$\chi^2$", width = 30)
        text = text + '''& {0:{width}.{precision}f} '''.format(float(obj.fit_results.chi2), width = width, precision = precision)
        text = text + '''\\\\
        '''
        text = text + '''{0:{width}s}'''.format("$\chi_{\\nu}^2$", width = 30)
        text = text + '''& {0:{width}.{precision}f} '''.format(float(obj.fit_results.reduced_chi2), width = width, precision = precision)
        text = text + '''\\\\
        '''
        text = text + '''{0:{width}s}'''.format("$rms$ [m\,s$^{-1}$]", width = 30)
        text = text + '''& {0:{width}.{precision}f} '''.format(float(obj.fit_results.rms), width = width, precision = precision)
        text = text + '''\\\\
        '''
        text = text + '''{0:{width}s}'''.format("$wrms$ [m\,s$^{-1}$]", width = 30)
        text = text + '''& {0:{width}.{precision}f} '''.format(float(obj.fit_results.wrms), width = width, precision = precision)
        text = text + '''\\\\
        '''

        text = text + '''{0:{width}s}'''.format("$-\ln\mathcal{L}$", width = 30)
        text = text + '''& {0:{width}.{precision}f} '''.format(float(obj.loglik), width = width, precision = precision)
        text = text + '''\\\\
        '''
        text = text + '''{0:{width}s}'''.format("N$_{\\rm RV}$ data", width = 30)
        text = text + '''& {0:{width}.{precision}f} '''.format(len(obj.fit_results.jd), width = width, precision = 0)
        text = text + '''\\\\
        '''

        text = text + '''{0:{width}s}'''.format("Epoch", width = 30)
        text = text + '''& {0:{width}.{precision}f} '''.format(obj.epoch, width = width, precision = precision)
        text = text + '''\\\\
        '''

        text = text + '''\\\\
    \hline \\noalign{\\vskip 0.7mm}

        '''

        text = text + '''
    \end{tabular}

    % \end{minipage}}
    % \end{adjustwidth}

    %\\tablefoot{\small }

    \end{table}
    '''

    elif asymmetric == True:

        text = '''
    \\begin{table}[ht]
    % \\begin{adjustwidth}{-4.0cm}{}
    % \\resizebox{0.69\\textheight}{!}
    % {\\begin{minipage}{1.1\\textwidth}

    \centering
    \caption{{}}
    \label{table:}

    \\begin{tabular}{lrrrrrrrr}     % 2 columns

    \hline\hline  \\noalign{\\vskip 0.7mm}
    '''


        text = text + '''Parameter \hspace{0.0 mm}'''
        for i in range(9):
            if not bool(obj.use_planet[i]):
                continue
            text = text + '''& Planet %s '''%chr(98+i)
        text = text + '''\\\\
    \hline \\noalign{\\vskip 0.7mm}

        '''

        if obj.type_fit["RV"]== True or obj.type_fit["TTV"]== True:

            text = text + '''{0:{width}s}'''.format("$K$ [m\,s$^{-1}$]", width = 30)
            for i in range(9):
                if not bool(obj.use_planet[i]):
                    continue
                text = text + '''& {0:{width}.{precision}f}$_{{-{1:{width2}.{precision}f}}}^{{+{2:{width2}.{precision}f}}}$ '''.format(obj.K[i], obj.param_errors.planet_params_errors[7*i][0], obj.param_errors.planet_params_errors[7*i][1], width = width, width2 = 0, precision = precision)
            text = text + '''\\\\ \\noalign{\\vskip 0.9mm}
        '''

        text = text + '''{0:{width}s}'''.format("$P$ [day]", width = 30)
        for i in range(9):
            if not bool(obj.use_planet[i]):
                continue
            text = text + '''& {0:{width}.{precision}f}$_{{-{1:{width2}.{precision}f}}}^{{+{2:{width2}.{precision}f}}}$ '''.format(obj.P[i], obj.param_errors.planet_params_errors[7*i +1][0], obj.param_errors.planet_params_errors[7*i +1][1], width = width, width2 = 0, precision = precision)
        text = text + '''\\\\ \\noalign{\\vskip 0.9mm}
        '''

        if obj.hkl == False:
            text = text + '''{0:{width}s}'''.format("$e$ ", width = 30)
            for i in range(9):
                if not bool(obj.use_planet[i]):
                    continue
                text = text + '''& {0:{width}.{precision}f}$_{{-{1:{width2}.{precision}f}}}^{{+{2:{width2}.{precision}f}}}$ '''.format(obj.e[i], obj.param_errors.planet_params_errors[7*i +2][0], obj.param_errors.planet_params_errors[7*i +2][1], width = width, width2 = 0, precision = precision)
            text = text + '''\\\\ \\noalign{\\vskip 0.9mm}
        '''
            text = text + '''{0:{width}s}'''.format("$\omega$ [deg]", width = 30)
            for i in range(9):
                if not bool(obj.use_planet[i]):
                    continue
                text = text + '''& {0:{width}.{precision}f}$_{{-{1:{width2}.{precision}f}}}^{{+{2:{width2}.{precision}f}}}$ '''.format(obj.w[i], obj.param_errors.planet_params_errors[7*i +3][0], obj.param_errors.planet_params_errors[7*i +3][1], width = width, width2 = 0, precision = precision)
            text = text + '''\\\\ \\noalign{\\vskip 0.9mm}
        '''
        else:
            text = text + '''{0:{width}s}'''.format("$e sin(\omega)$  ", width = 30)
            for i in range(9):
                if not bool(obj.use_planet[i]):
                    continue
                text = text + '''& {0:{width}.{precision}f}$_{{-{1:{width2}.{precision}f}}}^{{+{2:{width2}.{precision}f}}}$'''.format(obj.e_sinw[i],  obj.e_sinw_err[i][0], obj.e_sinw_err[i][1], width = width, width2 = 0,  precision = precision)
            text = text + '''\\\\
        '''
            text = text + '''{0:{width}s}'''.format("$e cos(\omega)$", width = 30)
            for i in range(9):
                if not bool(obj.use_planet[i]):
                    continue
                text = text + '''& {0:{width}.{precision}f}$_{{-{1:{width2}.{precision}f}}}^{{+{2:{width2}.{precision}f}}}$'''.format(obj.e_cosw[i],  obj.e_cosw_err[i][0], obj.e_cosw_err[i][1], width = width, width2 = 0, precision = precision)
            text = text + '''\\\\
        '''
        if obj.type_fit["RV"]== True or obj.type_fit["TTV"]== True:

            if obj.hkl == False:
                text = text + '''{0:{width}s}'''.format("$M_{\\rm 0}$ [deg]", width = 30)
            for i in range(9):
                if not bool(obj.use_planet[i]):
                    continue
                    text = text + '''& {0:{width}.{precision}f}$_{{-{1:{width2}.{precision}f}}}^{{+{2:{width2}.{precision}f}}}$ '''.format(obj.M0[i], obj.param_errors.planet_params_errors[7*i +4][0], obj.param_errors.planet_params_errors[7*i +4][1], width = width, width2 = 0, precision = precision)
                text = text + '''\\\\ \\noalign{\\vskip 0.9mm}
        '''
            else:
                text = text + '''{0:{width}s}'''.format("$\lambda$ [deg]", width = 30)
                for i in range(9):
                    if not bool(obj.use_planet[i]):
                        continue
                    text = text + '''& {0:{width}.{precision}f}$_{{-{1:{width2}.{precision}f}}}^{{+{2:{width2}.{precision}f}}}$ '''.format(obj.lamb[i],  obj.lamb_err[i][0], obj.lamb_err[i][1], width = width,  width2 = 0, precision = precision)
                text = text + '''\\\\
        '''
        if obj.mod_dynamical == True:
            text = text + '''{0:{width}s}'''.format("$i$ [deg]", width = 30)
            for i in range(9):
                if not bool(obj.use_planet[i]):
                    continue
                text = text + '''& {0:{width}.{precision}f}$_{{-{1:{width2}.{precision}f}}}^{{+{2:{width2}.{precision}f}}}$ '''.format(obj.i[i], obj.param_errors.planet_params_errors[7*i +5][0], obj.param_errors.planet_params_errors[7*i +5][1], width = width, width2 = 0, precision = precision)
            text = text + '''\\\\ \\noalign{\\vskip 0.9mm}
        '''
            text = text + '''{0:{width}s}'''.format("$\Omega$ [deg]", width = 30)
            for i in range(9):
                if not bool(obj.use_planet[i]):
                    continue
                text = text + '''& {0:{width}.{precision}f}$_{{-{1:{width2}.{precision}f}}}^{{+{2:{width2}.{precision}f}}}$ '''.format(obj.Node[i], obj.param_errors.planet_params_errors[7*i +6][0], obj.param_errors.planet_params_errors[7*i +6][1], width = width, width2 = 0, precision = precision)
            text = text + '''\\\\ \\noalign{\\vskip 0.9mm}
        '''


        if obj.type_fit["Transit"] == True:
            text = text + '''{0:{width}s}'''.format("$t_{\\rm 0}$ [day]", width = 30)
            for i in range(9):
                if not bool(obj.use_planet[i]):
                    continue
                text = text + '''& {0:{width}.{precision}f}$_{{-{1:{width2}.{precision}f}}}^{{+{2:{width2}.{precision}f}}}$ '''.format(obj.t0[i], obj.t0_err[i][0], obj.t0_err[i][1], width = width, width2 = 0, precision = precision)
            text = text + '''\\\\ \\noalign{\\vskip 0.9mm}
        '''
            text = text + '''{0:{width}s}'''.format("Rad. [$R_\odot$]", width = 30)
            for i in range(9):
                if not bool(obj.use_planet[i]):
                    continue
                text = text + '''& {0:{width}.{precision}f}$_{{-{1:{width2}.{precision}f}}}^{{+{2:{width2}.{precision}f}}}$ '''.format(obj.pl_rad[i], obj.pl_rad_err[i][0], obj.pl_rad_err[i][1], width = width, width2 = 0, precision = precision)
            text = text + '''\\\\ \\noalign{\\vskip 0.9mm}
        '''
            text = text + '''{0:{width}s}'''.format("$a$ [$R_\odot$]", width = 30)
            for i in range(9):
                if not bool(obj.use_planet[i]):
                    continue
                text = text + '''& {0:{width}.{precision}f}$_{{-{1:{width2}.{precision}f}}}^{{+{2:{width2}.{precision}f}}}$ '''.format(obj.pl_a[i], obj.pl_a_err[i][0], obj.pl_a_err[i][1], width = width, width2 = 0, precision = precision)
            text = text + '''\\\\ \\noalign{\\vskip 0.9mm}
        '''

        text = text + '''{0:{width}s}'''.format("$a$ [au]", width = 30)
        for i in range(9):
            if not bool(obj.use_planet[i]):
                continue
            text = text + '''& {0:{width}.{precision}f}$_{{-{1:{width2}.{precision}f}}}^{{+{2:{width2}.{precision}f}}}$ '''.format(obj.fit_results.a[i], 0,0, width = width, width2 = 0, precision = precision)
        text = text + '''\\\\ \\noalign{\\vskip 0.9mm}
        '''
        text = text + '''{0:{width}s}'''.format("$m \sin i$ [$M_{\\rm jup}$]", width = 30)
        for i in range(9):
            if not bool(obj.use_planet[i]):
                continue
            text = text + '''& {0:{width}.{precision}f}$_{{-{1:{width2}.{precision}f}}}^{{+{2:{width2}.{precision}f}}}$ '''.format(obj.fit_results.mass[i], 0,0, width = width, width2 = 0, precision = precision)
        text = text + '''\\\\ \\noalign{\\vskip 0.9mm}
        '''
        text = text + '''{0:{width}s}'''.format("$t_{\omega}$ [day]", width = 30)
        for i in range(9):
            if not bool(obj.use_planet[i]):
                continue
            text = text + '''& {0:{width}.{precision}f}$_{{-{1:{width2}.{precision}f}}}^{{+{2:{width2}.{precision}f}}}$ '''.format((float(obj.epoch) - (np.radians(obj.M0[i])/(2*np.pi))*obj.P[i] ), 0,0, width = width, width2 = 0, precision = precision)
        text = text + '''\\\\ \\noalign{\\vskip 0.9mm}
        '''

        if obj.type_fit["RV"]== True:

            text = text + '''{0:{width}s}'''.format("RV lin. trend [m\,s$^{-1}$\,day$^{-1}$]", width = 30)
            text = text + '''& {0:{width}.{precision}f}$_{{-{1:{width2}.{precision}f}}}^{{+{2:{width2}.{precision}f}}}$ '''.format(float(obj.rv_lintr),float(obj.rv_lintr_err[0]),float(obj.rv_lintr_err[1]) , width = width, width2 = 0, precision = precision)
            text = text + '''\\\\
        '''

            text = text + '''{0:{width}s}'''.format("RV quad. trend [m\,s$^{-1}$\,day$^{-1}$]", width = 30)
            text = text + '''& {0:{width}.{precision}f}$_{{-{1:{width2}.{precision}f}}}^{{+{2:{width2}.{precision}f}}}$  '''.format(float(obj.rv_quadtr),float(obj.rv_quadtr_err[0]),float(obj.rv_quadtr_err[1]) , width = width, width2 = 0, precision = precision)
            text = text + '''\\\\
        '''

            for i in range(obj.ndset):
                text = text + '''{0:{width}s}'''.format("RV$_{\\rm off}$ %s [m\,s$^{-1}$]"%(i+1), width = 30)
                text = text + '''& {0:{width}.{precision}f}$_{{-{1:{width2}.{precision}f}}}^{{+{2:{width2}.{precision}f}}}$ '''.format(float(obj.rvoff[i]), obj.rvoff_err[i][0], obj.rvoff_err[i][1], width = width, width2 = 0, precision = precision)
                text = text + '''\\\\ \\noalign{\\vskip 0.9mm}
        '''
            for i in range(obj.ndset):
                text = text + '''{0:{width}s}'''.format("RV$_{\\rm jit}$ %s [m\,s$^{-1}$]"%(i+1), width = 30)
                text = text + '''& {0:{width}.{precision}f}$_{{-{1:{width2}.{precision}f}}}^{{+{2:{width2}.{precision}f}}}$ '''.format(float(obj.jitt[i]), obj.jitt_err[i][0], obj.jitt_err[i][1], width = width, width2 = 0, precision = precision)
                text = text + '''\\\\ \\noalign{\\vskip 0.9mm}
        '''

        if obj.type_fit["Transit"]== True:

            for i in range(20):
                if len(obj.tra_data_sets[i]) != 0:
                    text = text + '''{0:{width}s}'''.format("Tran.$_{\\rm off}$ %s"%(i+1), width = 30)
                    text = text + '''& {0:{width}.{precision}f}$_{{-{1:{width2}.{precision}f}}}^{{+{2:{width2}.{precision}f}}}$ '''.format(float(obj.tra_off[i]), obj.tra_off_err[i][0], obj.tra_off_err[i][1], width = width, width2 = 0, precision = precision)
                    text = text + '''\\\\ \\noalign{\\vskip 0.9mm}
        '''
            for i in range(20):
                if len(obj.tra_data_sets[i]) != 0:
                    text = text + '''{0:{width}s}'''.format("Tran.$_{\\rm jit}$ %s"%(i+1), width = 30)
                    text = text + '''& {0:{width}.{precision}f}$_{{-{1:{width2}.{precision}f}}}^{{+{2:{width2}.{precision}f}}}$ '''.format(float(obj.tra_jitt[i]), obj.tra_jitt_err[i][0], obj.tra_jitt_err[i][1], width = width, width2 = 0, precision = precision)
                    text = text + '''\\\\ \\noalign{\\vskip 0.9mm}
        '''

        if obj.doGP == True:

            if obj.gp_kernel == 'RotKernel':
                for i in range(4):
                    text = text + '''{0:{width}s}'''.format("%s"%(obj.GP_rot_str[i]), width = 30)
                    text = text + '''& {0:{width}.{precision}f}$_{{-{1:{width2}.{precision}f}}}^{{+{2:{width2}.{precision}f}}}$ '''.format(float(obj.GP_rot_params[i]), obj.param_errors.GP_params_errors[i][0], obj.param_errors.GP_params_errors[i][1], width = width, width2 = 0, precision = precision)
                    text = text + '''\\\\ \\noalign{\\vskip 0.9mm}
        '''


            elif obj.gp_kernel == 'SHOKernel':
                for i in range(3):
                    text = text + '''{0:{width}s}'''.format("%s"%(obj.GP_sho_str[i]), width = 30)
                    text = text + '''& {0:{width}.{precision}f}$_{{-{1:{width2}.{precision}f}}}^{{+{2:{width2}.{precision}f}}}$ '''.format(float(obj.GP_sho_params[i]), obj.param_errors.GP_params_errors[i][0], obj.param_errors.GP_params_errors[i][1], width = width, width2 = 0, precision = precision)
                    text = text + '''\\\\ \\noalign{\\vskip 0.9mm}
        '''


        text = text + '''{0:{width}s}'''.format("$\chi^2$", width = 30)
        text = text + '''& {0:{width}.{precision}f} '''.format(float(obj.fit_results.chi2), width = width, precision = precision)
        text = text + '''\\\\
        '''
        text = text + '''{0:{width}s}'''.format("$\chi_{\\nu}^2$", width = 30)
        text = text + '''& {0:{width}.{precision}f} '''.format(float(obj.fit_results.reduced_chi2), width = width, precision = precision)
        text = text + '''\\\\
        '''
        text = text + '''{0:{width}s}'''.format("$rms$ [m\,s$^{-1}$]", width = 30)
        text = text + '''& {0:{width}.{precision}f} '''.format(float(obj.fit_results.rms), width = width, precision = precision)
        text = text + '''\\\\
        '''
        text = text + '''{0:{width}s}'''.format("$wrms$ [m\,s$^{-1}$]", width = 30)
        text = text + '''& {0:{width}.{precision}f} '''.format(float(obj.fit_results.wrms), width = width, precision = precision)
        text = text + '''\\\\
        '''

        text = text + '''{0:{width}s}'''.format("$-\ln\mathcal{L}$", width = 30)
        text = text + '''& {0:{width}.{precision}f} '''.format(float(obj.loglik), width = width, precision = precision)
        text = text + '''\\\\
        '''
        text = text + '''{0:{width}s}'''.format("N$_{\\rm RV}$ data", width = 30)
        text = text + '''& {0:{width}.{precision}f} '''.format(len(obj.fit_results.jd), width = width, precision = 0)
        text = text + '''\\\\
        '''

        text = text + '''{0:{width}s}'''.format("Epoch", width = 30)
        text = text + '''& {0:{width}.{precision}f} '''.format(obj.epoch, width = width, precision = precision)
        text = text + '''\\\\
        '''

        text = text + '''\\\\
    \hline \\noalign{\\vskip 0.7mm}

        '''

        text = text + '''
    \end{tabular}

    % \end{minipage}}
    % \end{adjustwidth}

    %\\tablefoot{\small }

    \end{table}
    '''

    else:
        print("asymmetric must be True or False")
        return


    if return_text == True:
        return text
    else:

        table_file = open("%s/%s"%(path,file_name), 'w')
        table_file.write(text)
        table_file.close()
        print("Done")
        return





def latex_prior_table(obj, width = 10, precision = 2,  file_name='prior_table.tex', path='./', return_text = False):


    text = '''
\\begin{table}[ht]
% \\begin{adjustwidth}{-4.0cm}{}
% \\resizebox{0.69\\textheight}{!}
% {\\begin{minipage}{1.1\\textwidth}

\centering
\caption{{}}
\label{table:}

\\begin{tabular}{lrrrrrrrr}     % 2 columns

\hline\hline  \\noalign{\\vskip 0.7mm}
'''


    text = text + '''Parameter \hspace{0.0 mm}'''
    for i in range(9):
        if not bool(obj.use_planet[i]):
            continue
        text = text + '''& Planet %s '''%chr(98+i)
    text = text + '''\\\\
\hline \\noalign{\\vskip 0.7mm}

    '''
    if obj.type_fit["RV"] == True or obj.type_fit["TTV"] == True :
        text = text + '''{0:{width}s}'''.format("$K$ [m\,s$^{-1}$]", width = 30)
        for i in range(9):
            if not bool(obj.use_planet[i]):
                continue
            if obj.K_norm_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg = "$\mathcal{N}$",obj.K_norm_pr[i][0],obj.K_norm_pr[i][1],"$^2$"
            elif obj.K_jeff_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg = "$\mathcal{J}$",obj.K_jeff_pr[i][0],obj.K_jeff_pr[i][1],""
            else:
                sign,f_arg,s_arg,pow_arg = "$\mathcal{U}$",obj.K_bound[i][0],obj.K_bound[i][1],""

            text = text + '''& {0:s}({1:{width}.{precision}f},{2:{width}.{precision}f}{3:s})'''.format(sign, f_arg,s_arg,pow_arg, width = width, precision = precision)
        text = text + '''\\\\
    '''

    text = text + '''{0:{width}s}'''.format("$P$ [day]", width = 30)
    for i in range(9):
        if not bool(obj.use_planet[i]):
            continue
        if obj.P_norm_pr[i][2]==True:
            sign,f_arg,s_arg,pow_arg  = "$\mathcal{N}$",obj.P_norm_pr[i][0],obj.P_norm_pr[i][1],"$^2$"
        elif obj.P_jeff_pr[i][2]==True:
            sign,f_arg,s_arg,pow_arg  = "$\mathcal{J}$",obj.P_jeff_pr[i][0],obj.P_jeff_pr[i][1],""
        else:
            sign,f_arg,s_arg,pow_arg  = "$\mathcal{U}$",obj.P_bound[i][0],obj.P_bound[i][1],""

        text = text + '''& {0:s}({1:{width}.{precision}f},{2:{width}.{precision}f}{3:s})'''.format(sign, f_arg,s_arg,pow_arg, width = width, precision = precision)
    text = text + '''\\\\
    '''

    if obj.hkl == False:
        text = text + '''{0:{width}s}'''.format("$e$  ", width = 30)
        for i in range(9):
            if not bool(obj.use_planet[i]):
                continue
            if obj.e_norm_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{N}$",obj.e_norm_pr[i][0],obj.e_norm_pr[i][1],"$^2$"
            elif obj.e_jeff_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{J}$",obj.e_jeff_pr[i][0],obj.e_jeff_pr[i][1],""
            else:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{U}$",obj.e_bound[i][0],obj.e_bound[i][1],""

            text = text + '''& {0:s}({1:{width}.{precision}f},{2:{width}.{precision}f}{3:s})'''.format(sign, f_arg,s_arg,pow_arg, width = width, precision = precision)
        text = text + '''\\\\
    '''

        text = text + '''{0:{width}s}'''.format("$\omega$ [deg]", width = 30)
        for i in range(9):
            if not bool(obj.use_planet[i]):
                continue
            if obj.w_norm_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{N}$",obj.w_norm_pr[i][0],obj.w_norm_pr[i][1],"$^2$"
            elif obj.w_jeff_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{J}$",obj.w_jeff_pr[i][0],obj.w_jeff_pr[i][1],""
            else:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{U}$",obj.w_bound[i][0],obj.w_bound[i][1],""

            text = text + '''& {0:s}({1:{width}.{precision}f},{2:{width}.{precision}f}{3:s})'''.format(sign, f_arg,s_arg,pow_arg, width = width, precision = precision)
        text = text + '''\\\\
    '''

        if obj.type_fit["RV"] == True or obj.type_fit["TTV"] == True :
            text = text + '''{0:{width}s}'''.format("$M_{\\rm 0}$ [deg]", width = 30)
            for i in range(9):
                if not bool(obj.use_planet[i]):
                    continue
                if obj.M0_norm_pr[i][2]==True:
                    sign,f_arg,s_arg,pow_arg  = "$\mathcal{N}$",obj.M0_norm_pr[i][0],obj.M0_norm_pr[i][1],"$^2$"
                elif obj.M0_jeff_pr[i][2]==True:
                    sign,f_arg,s_arg,pow_arg  = "$\mathcal{J}$",obj.M0_jeff_pr[i][0],obj.M0_jeff_pr[i][1],""
                else:
                    sign,f_arg,s_arg,pow_arg  = "$\mathcal{U}$",obj.M0_bound[i][0],obj.M0_bound[i][1],""

                text = text + '''& {0:s}({1:{width}.{precision}f},{2:{width}.{precision}f}{3:s})'''.format(sign, f_arg,s_arg,pow_arg, width = width, precision = precision)
            text = text + '''\\\\
    '''

    elif obj.hkl == True:

        text = text + '''{0:{width}s}'''.format("$e\sin(\omega)$  ", width = 30)
        for i in range(9):
            if not bool(obj.use_planet[i]):
                continue
            if obj.e_sinw_norm_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{N}$",obj.e_sinw_norm_pr[i][0],obj.e_sinw_norm_pr[i][1],"$^2$"
            elif obj.e_sinw_jeff_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{J}$",obj.e_sinw_jeff_pr[i][0],obj.e_sinw_jeff_pr[i][1],""
            else:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{U}$",obj.e_sinw_bound[i][0],obj.e_sinw_bound[i][1],""

            text = text + '''& {0:s}({1:{width}.{precision}f},{2:{width}.{precision}f}{3:s})'''.format(sign, f_arg,s_arg,pow_arg, width = width, precision = precision)
        text = text + '''\\\\
    '''

        text = text + '''{0:{width}s}'''.format("$e\cos(\omega)$  ", width = 30)
        for i in range(9):
            if not bool(obj.use_planet[i]):
                continue
            if obj.e_cosw_norm_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{N}$",obj.e_cosw_norm_pr[i][0],obj.e_cosw_norm_pr[i][1],"$^2$"
            elif obj.e_cosw_jeff_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{J}$",obj.e_cosw_jeff_pr[i][0],obj.e_cosw_jeff_pr[i][1],""
            else:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{U}$",obj.e_cosw_bound[i][0],obj.e_cosw_bound[i][1],""

            text = text + '''& {0:s}({1:{width}.{precision}f},{2:{width}.{precision}f}{3:s})'''.format(sign, f_arg,s_arg,pow_arg, width = width, precision = precision)
        text = text + '''\\\\
    '''

        text = text + '''{0:{width}s}'''.format("$\lambda$ [deg]", width = 30)
        for i in range(9):
            if not bool(obj.use_planet[i]):
                continue
            if obj.lamb_norm_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{N}$",obj.lamb_norm_pr[i][0],obj.lamb_norm_pr[i][1],"$^2$"
            elif obj.lamb_jeff_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{J}$",obj.lamb_jeff_pr[i][0],obj.lamb_jeff_pr[i][1],""
            else:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{U}$",obj.lamb_bound[i][0],obj.lamb_bound[i][1],""

            text = text + '''& {0:s}({1:{width}.{precision}f},{2:{width}.{precision}f}{3:s})'''.format(sign, f_arg,s_arg,pow_arg, width = width, precision = precision)
        text = text + '''\\\\
    '''

    if obj.mod_dynamical == True:
        text = text + '''{0:{width}s}'''.format("$i$ [deg]", width = 30)
        for i in range(9):
            if not bool(obj.use_planet[i]):
                continue
            if obj.i_norm_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{N}$",obj.i_norm_pr[i][0],obj.i_norm_pr[i][1],"$^2$"
            elif obj.i_jeff_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{J}$",obj.i_jeff_pr[i][0],obj.i_jeff_pr[i][1],""
            else:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{U}$",obj.i_bound[i][0],obj.i_bound[i][1],""

            text = text + '''& {0:s}({1:{width}.{precision}f},{2:{width}.{precision}f}{3:s})'''.format(sign, f_arg,s_arg,pow_arg, width = width, precision = precision)
        text = text + '''\\\\
    '''

        text = text + '''{0:{width}s}'''.format("$\Omega$ [deg]", width = 30)
        for i in range(9):
            if not bool(obj.use_planet[i]):
                continue
            if obj.Node_norm_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{N}$",obj.Node_norm_pr[i][0],obj.Node_norm_pr[i][1],"$^2$"
            elif obj.Node_jeff_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{J}$",obj.Node_jeff_pr[i][0],obj.Node_jeff_pr[i][1],""
            else:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{U}$",obj.Node_bound[i][0],obj.Node_bound[i][1],""

            text = text + '''& {0:s}({1:{width}.{precision}f},{2:{width}.{precision}f}{3:s})'''.format(sign, f_arg,s_arg,pow_arg, width = width, precision = precision)
        text = text + '''\\\\
    '''

    if obj.type_fit["Transit"] == True:

        text = text + '''{0:{width}s}'''.format("$t_{\\rm 0}$ [day]", width = 30)
        for i in range(9):
            if not bool(obj.use_planet[i]):
                continue
            if obj.i_norm_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{N}$",obj.t0_norm_pr[i][0],obj.t0_norm_pr[i][1],"$^2$"
            elif obj.t0_jeff_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{J}$",obj.t0_jeff_pr[i][0],obj.t0_jeff_pr[i][1],""
            else:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{U}$",obj.t0_bound[i][0],obj.t0_bound[i][1],""

            text = text + '''& {0:s}({1:{width}.{precision}f},{2:{width}.{precision}f}{3:s})'''.format(sign, f_arg,s_arg,pow_arg, width = width, precision = precision)
        text = text + '''\\\\
    '''

        text = text + '''{0:{width}s}'''.format("Rp/$R_\star$", width = 30)
        for i in range(9):
            if not bool(obj.use_planet[i]):
                continue
            if obj.pl_rad_norm_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{N}$",obj.pl_rad_norm_pr[i][0],obj.pl_rad_norm_pr[i][1],"$^2$"
            elif obj.pl_rad_jeff_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{J}$",obj.pl_rad_jeff_pr[i][0],obj.pl_rad_jeff_pr[i][1],""
            else:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{U}$",obj.pl_rad_bound[i][0],obj.pl_rad_bound[i][1],""

            text = text + '''& {0:s}({1:{width}.{precision}f},{2:{width}.{precision}f}{3:s})'''.format(sign, f_arg,s_arg,pow_arg, width = width, precision = precision)
        text = text + '''\\\\
    '''


        text = text + '''{0:{width}s}'''.format("a/$R_\star$", width = 30)
        for i in range(9):
            if not bool(obj.use_planet[i]):
                continue
            if obj.pl_a_norm_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{N}$",obj.pl_a_norm_pr[i][0],obj.pl_a_norm_pr[i][1],"$^2$"
            elif obj.pl_a_jeff_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{J}$",obj.pl_a_jeff_pr[i][0],obj.pl_a_jeff_pr[i][1],""
            else:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{U}$",obj.pl_a_bound[i][0],obj.pl_a_bound[i][1],""

            text = text + '''& {0:s}({1:{width}.{precision}f},{2:{width}.{precision}f}{3:s})'''.format(sign, f_arg,s_arg,pow_arg, width = width, precision = precision)
        text = text + '''\\\\
    '''


    if obj.type_fit["RV"]== True:

        text = text + '''{0:{width}s}'''.format("RV lin. trend [m\,s$^{-1}$\,day$^{-1}$]", width = 30)
        i = 0
        if obj.rv_lintr_norm_pr[i][2]==True:
            sign,f_arg,s_arg,pow_arg  = "$\mathcal{N}$",obj.rv_lintr_norm_pr[i][0],obj.rv_lintr_norm_pr[i][1],"$^2$"
        elif obj.rv_lintr_jeff_pr[i][2]==True:
            sign,f_arg,s_arg,pow_arg  = "$\mathcal{J}$",obj.rv_lintr_jeff_pr[i][0],obj.rv_lintr_jeff_pr[i][1],""
        else:
            sign,f_arg,s_arg,pow_arg  = "$\mathcal{U}$",obj.rv_lintr_bounds[i][0],obj.rv_lintr_bounds[i][1],""

        text = text + '''& {0:s}({1:{width}.{precision}f},{2:{width}.{precision}f}{3:s})'''.format(sign, f_arg,s_arg,pow_arg, width = width, precision = precision)
        text = text + '''\\\\
    '''

        text = text + '''{0:{width}s}'''.format("RV quad. trend [m\,s$^{-1}$\,day$^{-1}$]", width = 30)
        i = 0
        if obj.rv_quadtr_norm_pr[i][2]==True:
            sign,f_arg,s_arg,pow_arg  = "$\mathcal{N}$",obj.rv_quadtr_norm_pr[i][0],obj.rv_quadtr_norm_pr[i][1],"$^2$"
        elif obj.rv_quadtr_jeff_pr[i][2]==True:
            sign,f_arg,s_arg,pow_arg  = "$\mathcal{J}$",obj.rv_quadtr_jeff_pr[i][0],obj.rv_quadtr_jeff_pr[i][1],""
        else:
            sign,f_arg,s_arg,pow_arg  = "$\mathcal{U}$",obj.rv_quadtr_bounds[i][0],obj.rv_quadtr_bounds[i][1],""

        text = text + '''& {0:s}({1:{width}.{precision}f},{2:{width}.{precision}f}{3:s})'''.format(sign, f_arg,s_arg,pow_arg, width = width, precision = precision)
        text = text + '''\\\\
    '''

        for i in range(obj.ndset):
            text = text + '''{0:{width}s}'''.format("RV$_{\\rm off}$ %s [m\,s$^{-1}$]"%(i+1), width = 30)
            if obj.rvoff_norm_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{N}$",obj.rvoff_norm_pr[i][0],obj.rvoff_norm_pr[i][1],"$^2$"
            elif obj.rvoff_jeff_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{J}$",obj.rvoff_jeff_pr[i][0],obj.rvoff_jeff_pr[i][1],""
            else:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{U}$",obj.rvoff_bounds[i][0],obj.rvoff_bounds[i][1],""

            text = text + '''& {0:s}({1:{width}.{precision}f},{2:{width}.{precision}f}{3:s})'''.format(sign, f_arg,s_arg,pow_arg, width = width, precision = precision)
            text = text + '''\\\\
    '''

        for i in range(obj.ndset):
            text = text + '''{0:{width}s}'''.format("RV$_{\\rm jit}$ %s [m\,s$^{-1}$]"%(i+1), width = 30)
            if obj.jitt_norm_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{N}$",obj.jitt_norm_pr[i][0],obj.jitt_norm_pr[i][1],"$^2$"
            elif obj.jitt_jeff_pr[i][2]==True:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{J}$",obj.jitt_jeff_pr[i][0],obj.jitt_jeff_pr[i][1],""
            else:
                sign,f_arg,s_arg,pow_arg  = "$\mathcal{U}$",obj.jitt_bounds[i][0],obj.jitt_bounds[i][1],""

            text = text + '''& {0:s}({1:{width}.{precision}f},{2:{width}.{precision}f}{3:s})'''.format(sign, f_arg,s_arg,pow_arg, width = width, precision = precision)
            text = text + '''\\\\
    '''

        if obj.doGP == True:

            if obj.gp_kernel == 'RotKernel':
                for i in range(4):
                    text = text + '''{0:{width}s}'''.format("%s"%(obj.GP_rot_str[i]), width = 30)
                    if obj.GP_rot_norm_pr[i][2]==True:
                        sign,f_arg,s_arg,pow_arg  = "$\mathcal{N}$",obj.GP_rot_norm_pr[i][0],obj.GP_rot_norm_pr[i][1],"$^2$"
                    elif obj.GP_rot_jeff_pr[i][2]==True:
                        sign,f_arg,s_arg,pow_arg  = "$\mathcal{J}$",obj.GP_rot_jeff_pr[i][0],obj.GP_rot_jeff_pr[i][1],""
                    else:
                        sign,f_arg,s_arg,pow_arg  = "$\mathcal{U}$",obj.GP_rot_bounds[i][0],obj.GP_rot_bounds[i][1],""

                    text = text + '''& {0:s}({1:{width}.{precision}f},{2:{width}.{precision}f}{3:s})'''.format(sign, f_arg,s_arg,pow_arg, width = width, precision = precision)
                    text = text + '''\\\\
    '''

            elif obj.gp_kernel == 'SHOKernel':
                for i in range(3):
                    text = text + '''{0:{width}s}'''.format("%s"%(obj.GP_sho_str[i]), width = 30)
                    if obj.GP_sho_norm_pr[i][2]==True:
                        sign,f_arg,s_arg,pow_arg  = "$\mathcal{N}$",obj.GP_sho_norm_pr[i][0],obj.GP_sho_norm_pr[i][1],"$^2$"
                    elif obj.GP_sho_jeff_pr[i][2]==True:
                        sign,f_arg,s_arg,pow_arg  = "$\mathcal{J}$",obj.GP_sho_jeff_pr[i][0],obj.GP_sho_jeff_pr[i][1],""
                    else:
                        sign,f_arg,s_arg,pow_arg  = "$\mathcal{U}$",obj.GP_sho_bounds[i][0],obj.GP_sho_bounds[i][1],""

                    text = text + '''& {0:s}({1:{width}.{precision}f},{2:{width}.{precision}f}{3:s})'''.format(sign, f_arg,s_arg,pow_arg, width = width, precision = precision)
                    text = text + '''\\\\
    '''

    if obj.type_fit["Transit"]== True:



        for i in range(20):
            if len(obj.tra_data_sets[i]) != 0:
                text = text + '''{0:{width}s}'''.format("Tran.$_{\\rm off}$ %s [m\,s$^{-1}$]"%(i+1), width = 30)
                if obj.tra_off_norm_pr[i][2]==True:
                    sign,f_arg,s_arg,pow_arg  = "$\mathcal{N}$",obj.tra_off_norm_pr[i][0],obj.tra_off_norm_pr[i][1],"$^2$"
                elif obj.tra_off_jeff_pr[i][2]==True:
                    sign,f_arg,s_arg,pow_arg  = "$\mathcal{J}$",obj.tra_off_jeff_pr[i][0],obj.tra_off_jeff_pr[i][1],""
                else:
                    sign,f_arg,s_arg,pow_arg  = "$\mathcal{U}$",obj.tra_off_bounds[i][0],obj.tra_off_bounds[i][1],""

                text = text + '''& {0:s}({1:{width}.{precision}f},{2:{width}.{precision}f}{3:s})'''.format(sign, f_arg,s_arg,pow_arg, width = width, precision = precision)
                text = text + '''\\\\
    '''

        for i in range(20):
            if len(obj.tra_data_sets[i]) != 0:
                text = text + '''{0:{width}s}'''.format("Tran.$_{\\rm jit}$ %s [m\,s$^{-1}$]"%(i+1), width = 30)
                if obj.tra_jitt_norm_pr[i][2]==True:
                    sign,f_arg,s_arg,pow_arg  = "$\mathcal{N}$",obj.tra_jitt_norm_pr[i][0],obj.tra_jitt_norm_pr[i][1],"$^2$"
                elif obj.tra_jitt_jeff_pr[i][2]==True:
                    sign,f_arg,s_arg,pow_arg  = "$\mathcal{J}$",obj.tra_jitt_jeff_pr[i][0],obj.tra_jitt_jeff_pr[i][1],""
                else:
                    sign,f_arg,s_arg,pow_arg  = "$\mathcal{U}$",obj.tra_jitt_bounds[i][0],obj.tra_jitt_bounds[i][1],""

                text = text + '''& {0:s}({1:{width}.{precision}f},{2:{width}.{precision}f}{3:s})'''.format(sign, f_arg,s_arg,pow_arg, width = width, precision = precision)
                text = text + '''\\\\
    '''

        if obj.tra_doGP == True:

            if obj.tra_gp_kernel == 'RotKernel':
                for i in range(4):
                    text = text + '''{0:{width}s}'''.format("%s"%(obj.tra_GP_rot_str[i]), width = 30)
                    if obj.tra_GP_rot_norm_pr[i][2]==True:
                        sign,f_arg,s_arg,pow_arg  = "$\mathcal{N}$",obj.tra_GP_rot_norm_pr[i][0],obj.tra_GP_rot_norm_pr[i][1],"$^2$"
                    elif obj.tra_GP_rot_jeff_pr[i][2]==True:
                        sign,f_arg,s_arg,pow_arg  = "$\mathcal{J}$",obj.tra_GP_rot_jeff_pr[i][0],obj.tra_GP_rot_jeff_pr[i][1],""
                    else:
                        sign,f_arg,s_arg,pow_arg  = "$\mathcal{U}$",obj.tra_GP_rot_bounds[i][0],obj.tra_GP_rot_bounds[i][1],""

                    text = text + '''& {0:s}({1:{width}.{precision}f},{2:{width}.{precision}f}{3:s})'''.format(sign, f_arg,s_arg,pow_arg, width = width, precision = precision)
                    text = text + '''\\\\
    '''

            elif obj.tra_gp_kernel == 'SHOKernel':
                for i in range(3):
                    text = text + '''{0:{width}s}'''.format("%s"%(obj.tra_GP_sho_str[i]), width = 30)
                    if obj.tra_GP_sho_norm_pr[i][2]==True:
                        sign,f_arg,s_arg,pow_arg  = "$\mathcal{N}$",obj.tra_GP_sho_norm_pr[i][0],obj.tra_GP_sho_norm_pr[i][1],"$^2$"
                    elif obj.tra_GP_sho_jeff_pr[i][2]==True:
                        sign,f_arg,s_arg,pow_arg  = "$\mathcal{J}$",obj.tra_GP_sho_jeff_pr[i][0],obj.tra_GP_sho_jeff_pr[i][1],""
                    else:
                        sign,f_arg,s_arg,pow_arg  = "$\mathcal{U}$",obj.tra_GP_sho_bounds[i][0],obj.tra_GP_sho_bounds[i][1],""

                    text = text + '''& {0:s}({1:{width}.{precision}f},{2:{width}.{precision}f}{3:s})'''.format(sign, f_arg,s_arg,pow_arg, width = width, precision = precision)
                    text = text + '''\\\\
    '''

    text = text + '''\\\\
\hline \\noalign{\\vskip 0.7mm}
    '''

    text = text + '''
\end{tabular}

% \end{minipage}}
% \end{adjustwidth}

%\\tablefoot{\small }

\end{table}
'''

    if return_text == True:
        return text
    else:

        table_file = open("%s/%s"%(path,file_name), 'w')
       # table_file = open(file_name, 'w')
        table_file.write(text)
        table_file.close()
        print("Done")
        return


def f_test(obj, obj2 = None, alpha = 0.01, lnL=False):


    chi2 = obj.fit_results.chi2
    ndata = len(obj.fit_results.jd)
    par2 = obj.fit_results.mfit
   #     self.value_reduced_chi2.setText("%.4f"%(fit.fit_results.reduced_chi2))
        #self.value_loglik.setText("%.4f"%(fit.fit_results.loglik))
   #     self.value_loglik.setText("%.4f"%(fit.loglik))


    if obj2 == None:
        obj2 = dill.copy(obj)
        obj2.npl = 0
        obj2.fitting()
    else:
        obj2 = dill.copy(obj2)

        if len(obj.fit_results.jd) != len(obj2.fit_results.jd):
            print("not the same data, test makes no sense")
            return


    chi1 = obj2.fit_results.chi2
    par1 = obj2.fit_results.mfit

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
"""%(chi1,chi2,obj.loglik,obj2.loglik,par1,par2,F,p_value,alpha))


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
def mass_a_from_Kepler_fit(a,npl,m0):
    '''Calculates the actual masses and Jacobi semimajor axes of a
       system for assumed sin(i) using the parameters P, K and e from a Kepler fit
       The output is now in Mjup and AU
    '''
    THIRD = 1.0/3.0
    PI    = 3.14159265358979e0
    TWOPI = 2.0*PI
    GMSUN = 1.32712440018e20
    AU=1.49597892e11
    incl = 90.0
    sini = np.sin(PI*(incl/180.0))
    mass  = np.zeros(npl+1)
    ap    = np.zeros(npl)
    pl_mass = np.zeros(npl)
    mpold = pl_mass

#*******G is set to be unit, and s, m, kg as unit of time, length and mass
#*****  and there is a reason for that! later I might correct for that.
    mtotal = m0
    f = 5e-6
    for i in range(npl):

        T = a[5*i+1]*86400.0
        mass[0] = m0

        # we need innitial guess for each planet mass
        dm = 0
        mass[i+1] = abs(a[5*i])*(T*(m0)**2.0/(TWOPI*GMSUN))**THIRD * np.sqrt(1.0-a[5*i+2]**2.0)/abs(sini)
        mpold[i] = mass[i+1]
        # This is a simple iteration to solve for mp
        while (dm <= 0):

            if i == 0:
                mtotal = m0
                mass[i+1] = abs(a[5*i])*(T*(m0 + mpold[i])**2.0/(TWOPI*GMSUN))**THIRD * np.sqrt(1.0-a[5*i+2]**2.0)/abs(sini)
            else:
                mtotal = m0
                for j in range(i):
                    mtotal = mtotal + mass[j+1]
                mass[i+1] = abs(a[5*i])*(T*(mtotal + mpold[i])**2.0/(TWOPI*GMSUN))**THIRD * np.sqrt(1.0-a[5*i+2]**2.0)/abs(sini)

            dm = (mpold[i] - mass[i+1])
            mpold[i] =  mpold[i] + f
           # print mass[i+1], mpold[i]

        ap[i] = (GMSUN*(mtotal + mass[i+1])*(T/TWOPI)**2)**THIRD

#    for i in range(npl+1):
#        mass[i] = mass[i]*GMSUN
    for i in range(npl):

        ap[i] = ap[i]/AU # to be in AU
        pl_mass[i] = mass[i+1]*1047.5654817267318 # to be in Jup. masses
        # I have seen that 1 Sol Mass = 1047.92612 Jup. masses???
    return pl_mass,ap



def check_swift(path='./'):

    lib_dir = PureWindowsPath(".\\lib")
    root_dir = PureWindowsPath(".\\")


    if not os.path.exists('%s'%os.path.join(Path(lib_dir),'libswift.a')):
        print("Installing the swift N-body lib for a first time!")
        result6, flag6 = run_command_with_timeout('%s'%os.path.join(Path(root_dir),'install_swift.sh'), 600,output=True)
        #print(result6)
        print("Installation DONE!")

def run_stability(obj, timemax=3000.0, timestep=10, timeout_sec=1000.0, stab_save_dir = 'run', remove_stab_save_dir = True, integrator='symba' ):


    stab_dir = PureWindowsPath(".\\stability")
    mvs_dir = PureWindowsPath("..\\mvs")
    mvsgr_dir = PureWindowsPath("..\\mvs_gr")
    symba_dir = PureWindowsPath("..\\symba")
    root_dir = PureWindowsPath("..\\..")
 
    os.chdir('%s'%Path(stab_dir))
    os.system("mkdir %s"%stab_save_dir)
    os.chdir("%s"%Path(stab_save_dir))


    print("running stability with: %s"%integrator)
    ##### crate the param.in file (change only the "t_max" and the "dt" for now) ######
    param_file = open('param.in', 'wb')

    max_time = float(timemax)*365.25 # make it is days

    param_file.write(b"""0.0d0 %s %s
%s %s

F T T T T F
0.0001 50.0 50.0 -1. T
bin.dat
unknown
"""%(bytes(str(max_time).encode()),
 bytes(str(timestep).encode()),
 bytes(str(max_time/1e4).encode()),
 bytes(str(max_time/1e3).encode())  ))

    param_file.close()

    #os.system("cp param.in test_param.in__")


    getin_file = open('geninit_j.in', 'wb')
    getin_file.write(b"""1
%s
%s
1.d0
pl.in
"""%(bytes(str(obj.params.stellar_mass).encode()), bytes(str(obj.npl).encode() ) ))



    for j in range(9):
        if not bool(obj.use_planet[j]):
            continue
        getin_file.write(b'%s \n'%bytes(str(obj.masses[j]/1047.5654817267318).encode()))
        getin_file.write(b'%s %s %s %s %s %s \n'%(bytes(str(obj.semimajor[j]).encode()),
                                                 bytes(str(obj.e[j]).encode()),
                                                 bytes(str(obj.i[j]).encode()),
                                                 bytes(str(obj.w[j]).encode()),
                                                 bytes(str(obj.Node[j]).encode()),
                                                 bytes(str(obj.M0[j]).encode() )) )

    getin_file.close()

    # runnning fortran codes
    result, flag = run_command_with_timeout('%s < geninit_j.in'%os.path.join(Path(mvs_dir),'geninit_j3_in_days'), timeout_sec)

    if integrator=='symba':
        result, flag = run_command_with_timeout('%s << EOF \nparam.in \npl.in \n1e-40 \nEOF'%os.path.join(Path(symba_dir),'swift_symba5_j'), timeout_sec)
    elif integrator=='mvs':
        result, flag = run_command_with_timeout('%s << EOF \nparam.in \npl.in \nEOF'%os.path.join(Path(mvs_dir),'swift_mvs_j'), timeout_sec)
    elif integrator=='mvs_gr':
        result, flag = run_command_with_timeout('%s << EOF \nparam.in \npl.in \n%s \nEOF'%(os.path.join(Path(mvsgr_dir),'swift_mvs_j_GR'),int(obj.GR_step)), timeout_sec)

    if not os.path.exists("energy.out"):
        os.chdir('%s'%Path(root_dir))
        print("something went wrong!!! No output generated.")
        return obj


    obj.evol_T_energy   = np.genfromtxt("energy.out",skip_header=0, unpack=True,skip_footer=1, usecols = [0]) /  365.25
    obj.evol_energy   = np.genfromtxt("energy.out",skip_header=0, unpack=True,skip_footer=1, usecols = [1])
   # obj.evol_momentum = np.genfromtxt("energy.out",skip_header=0, unpack=True,skip_footer=1, usecols = [2])

    obj.evol_momentum['lx'] = np.genfromtxt("energy.out",skip_header=0, unpack=True,skip_footer=1, usecols = [2])
    obj.evol_momentum['ly'] = np.genfromtxt("energy.out",skip_header=0, unpack=True,skip_footer=1, usecols = [3])
    obj.evol_momentum['lz'] = np.genfromtxt("energy.out",skip_header=0, unpack=True,skip_footer=1, usecols = [4])



    for k in range(obj.npl):

        if integrator=='symba':
            result, flag = run_command_with_timeout('%s << EOF \nparam.in \npl.in \n%s \nEOF'%(os.path.join(Path(symba_dir),'follow_symba2'),int(k+2)),timeout_sec)
            if "win" in sys.platform[0:3]:
                result, flag = run_command_with_timeout('ren follow_symba.out pl_%s.out'%(k+1),timeout_sec)
            else:
                result, flag = run_command_with_timeout('mv follow_symba.out pl_%s.out'%(k+1),timeout_sec)

        elif integrator=='mvs' or integrator=='mvs_gr':

            result, flag = run_command_with_timeout('%s << EOF \nparam.in \npl.in \n-%s \nEOF'%(os.path.join(Path(mvs_dir),'follow2'),int(k+2)),timeout_sec)
            if "win" in sys.platform[0:3]: 
                result, flag = run_command_with_timeout('ren follow2.out pl_%s.out'%(k+1),timeout_sec)
            else:
                result, flag = run_command_with_timeout('mv follow2.out pl_%s.out'%(k+1),timeout_sec)

        obj.evol_T[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [0]) /  365.25
        obj.evol_a[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [2])
        obj.evol_e[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [3])
        obj.evol_p[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [6])
        obj.evol_M[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [7])

        obj.evol_i[k]  = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [4])
        obj.evol_Om[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [5])

        obj.evol_Per[k] = a_to_P(obj.evol_a[k],obj.params.stellar_mass)



    os.chdir('..')
    if remove_stab_save_dir == True:

        if "win" in sys.platform[0:3]:
            os.system("rd -r %s"%stab_save_dir)
        else:
            os.system("rm -r %s"%stab_save_dir)
    os.chdir('..')

    print("stability with: %s done!"%integrator)

    return obj





def run_stability_old(obj, timemax=3000.0, timestep=10, timeout_sec=1000.0, stab_save_dir = './run', remove_stab_save_dir = True, integrator='symba' ):


    os.chdir('./stability/')
    os.system("mkdir %s"%stab_save_dir)
    os.chdir("./%s"%stab_save_dir)


    print("running stability with: %s"%integrator)
    ##### crate the param.in file (change only the "t_max" and the "dt" for now) ######
    param_file = open('param.in', 'wb')

    max_time = float(timemax)*365.25 # make it is days

    param_file.write(b"""0.0d0 %s %s
%s %s

F T T T T F
0.0001 50.0 50.0 -1. T
bin.dat
unknown
"""%(bytes(str(max_time).encode()),
 bytes(str(timestep).encode()),
 bytes(str(max_time/1e4).encode()),
 bytes(str(max_time/1e3).encode())  ))

    param_file.close()

    #os.system("cp param.in test_param.in__")


    getin_file = open('geninit_j.in', 'wb')
    getin_file.write(b"""1
%s
%s
1.d0
pl.in
"""%(bytes(str(obj.params.stellar_mass).encode()), bytes(str(obj.npl).encode() ) ))



    for j in range(9):
        if not bool(obj.use_planet[j]):
            continue
        #getin_file.write(b'%s \n'%bytes(str(obj.fit_results.mass[j]/1047.348644).encode()))
        getin_file.write(b'%s \n'%bytes(str(obj.masses[j]/1047.5654817267318).encode()))
#        getin_file.write(b'%s %s %s %s %s %s \n'%(bytes(str(obj.fit_results.a[j]).encode()),
        getin_file.write(b'%s %s %s %s %s %s \n'%(bytes(str(obj.semimajor[j]).encode()),
                                                 bytes(str(obj.e[j]).encode()),
                                                 bytes(str(obj.i[j]).encode()),
                                                 bytes(str(obj.w[j]).encode()),
                                                 bytes(str(obj.Node[j]).encode()),
                                                 bytes(str(obj.M0[j]).encode() )) )

    getin_file.close()

    # runnning fortran codes
    result, flag = run_command_with_timeout('../mvs/geninit_j3_in_days < geninit_j.in', timeout_sec)

    if integrator=='symba':
        result, flag = run_command_with_timeout('../symba/swift_symba5_j << EOF \nparam.in \npl.in \n1e-40 \nEOF', timeout_sec)
    elif integrator=='mvs':
        result, flag = run_command_with_timeout('../mvs/swift_mvs_j << EOF \nparam.in \npl.in \nEOF', timeout_sec)
    elif integrator=='mvs_gr':
        result, flag = run_command_with_timeout('../mvs_gr/swift_mvs_j_GR << EOF \nparam.in \npl.in \n%s \nEOF'%int(obj.GR_step), timeout_sec)

    #print('./swift_mvs_j_GR << EOF \nparam.in \npl.in \n%s \nEOF'%obj.GR_step)

    if not os.path.exists("energy.out"):
        os.chdir('../../')
        print("something went wrong!!! No output generated.")
        return obj


    obj.evol_T_energy   = np.genfromtxt("energy.out",skip_header=0, unpack=True,skip_footer=1, usecols = [0]) /  365.25
    obj.evol_energy   = np.genfromtxt("energy.out",skip_header=0, unpack=True,skip_footer=1, usecols = [1])
   # obj.evol_momentum = np.genfromtxt("energy.out",skip_header=0, unpack=True,skip_footer=1, usecols = [2])

    obj.evol_momentum['lx'] = np.genfromtxt("energy.out",skip_header=0, unpack=True,skip_footer=1, usecols = [2])
    obj.evol_momentum['ly'] = np.genfromtxt("energy.out",skip_header=0, unpack=True,skip_footer=1, usecols = [3])
    obj.evol_momentum['lz'] = np.genfromtxt("energy.out",skip_header=0, unpack=True,skip_footer=1, usecols = [4])



    for k in range(obj.npl):

        if integrator=='symba':
            result, flag = run_command_with_timeout('../symba/follow_symba2 << EOF \nparam.in \npl.in \n%s \nEOF'%(k+2),timeout_sec)
            result, flag = run_command_with_timeout('mv follow_symba.out pl_%s.out'%(k+1),timeout_sec)
        elif integrator=='mvs' or integrator=='mvs_gr':
            result, flag = run_command_with_timeout('../mvs/follow2 << EOF \nparam.in \npl.in \n-%s \nEOF'%(k+2),timeout_sec)
            result, flag = run_command_with_timeout('mv follow2.out pl_%s.out'%(k+1),timeout_sec)

        obj.evol_T[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [0]) /  365.25
        obj.evol_a[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [2])
        obj.evol_e[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [3])
        obj.evol_p[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [6])
        obj.evol_M[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [7])

        obj.evol_i[k]  = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [4])
        obj.evol_Om[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [5])

        obj.evol_Per[k] = a_to_P(obj.evol_a[k],obj.params.stellar_mass)



    os.chdir('../')
    if remove_stab_save_dir == True:
        os.system("rm -r %s"%stab_save_dir)
    os.chdir('../')

    print("stability with: %s done!"%integrator)

    return obj









def run_stability_arb(obj, timemax=3000.0, timestep=10, timeout_sec=1000.0, stab_save_dir = './', integrator='symba'):

#if not os.path.exists(directory):
#    os.makedirs(directory)



    if integrator=='symba':
        os.chdir('./stability/symba/')
    elif integrator=='mvs':
        os.chdir('./stability/mvs/')
    elif integrator=='mvs_gr':
        os.chdir('./stability/mvs_gr/')

    print("running stability with: %s"%integrator)
    ##### crate the param.in file (change only the "t_max" and the "dt" for now) ######
    param_file = open('param.in', 'wb')

    max_time = float(timemax)*365.25 # make it is days

    param_file.write(b"""0.0d0 %s %s
%s %s

F T T T T F
0.0001 50.0 50.0 -1. T
bin.dat
unknown
"""%(bytes(str(max_time).encode()),
 bytes(str(timestep).encode()),
 bytes(str(max_time/1e4).encode()),
 bytes(str(max_time/1e3).encode())  ))

    param_file.close()

    #os.system("cp param.in test_param.in__")


    getin_file = open('geninit_j.in', 'wb')
    getin_file.write(b"""1
%s
%s
1.d0
pl.in
"""%(bytes(str(obj.arb_st_mass).encode()), bytes(str(obj.npl_arb).encode() ) ))



    for j in range(9):
        if obj.pl_arb_use[j] == True:
            getin_file.write(b'%s \n'%bytes(str(obj.mass_arb[j]/1047.5654817267318).encode()))
            getin_file.write(b'%s %s %s %s %s %s \n'%(bytes(str(obj.a_arb[j]).encode()),
                                                 bytes(str(obj.e_arb[j]).encode()),
                                                 bytes(str(obj.i_arb[j]).encode()),
                                                 bytes(str(obj.w_arb[j]).encode()),
                                                 bytes(str(obj.Node_arb[j]).encode()),
                                                 bytes(str(obj.M0_arb[j]).encode() )) )
        else:
            continue

#
    getin_file.close()

    # runnning fortran codes
    result, flag = run_command_with_timeout('./geninit_j3_in_days < geninit_j.in', timeout_sec)

    if integrator=='symba':
        result, flag = run_command_with_timeout('./swift_symba5_j << EOF \nparam.in \npl.in \n1e-40 \nEOF', timeout_sec)
    elif integrator=='mvs':
        result, flag = run_command_with_timeout('./swift_mvs_j << EOF \nparam.in \npl.in \nEOF', timeout_sec)
    elif integrator=='mvs_gr':
        result, flag = run_command_with_timeout('./swift_mvs_j_GR << EOF \nparam.in \npl.in \n%s \nEOF'%int(obj.GR_step), timeout_sec)


    if not os.path.exists("energy.out"):
        os.chdir('../../')
        print("something went wrong!!! No output generated.")
        return obj

    obj.evol_T_energy   = np.genfromtxt("energy.out",skip_header=0, unpack=True,skip_footer=1, usecols = [0])  /  365.25
    obj.evol_energy   = np.genfromtxt("energy.out",skip_header=0, unpack=True,skip_footer=1, usecols = [1])
#    obj.evol_momentum = np.genfromtxt("energy.out",skip_header=0, unpack=True,skip_footer=1, usecols = [2])

    obj.evol_momentum['lx'] = np.genfromtxt("energy.out",skip_header=0, unpack=True,skip_footer=1, usecols = [2])
    obj.evol_momentum['ly'] = np.genfromtxt("energy.out",skip_header=0, unpack=True,skip_footer=1, usecols = [3])
    obj.evol_momentum['lz'] = np.genfromtxt("energy.out",skip_header=0, unpack=True,skip_footer=1, usecols = [4])


    for k in range(obj.npl_arb):

        if integrator=='symba':
            result, flag = run_command_with_timeout('./follow_symba2 << EOF \nparam.in \npl.in \n%s \nEOF'%(k+2),timeout_sec)
            result, flag = run_command_with_timeout('mv follow_symba.out pl_%s.out'%(k+1),timeout_sec)
        elif integrator=='mvs' or integrator=='mvs_gr':
            result, flag = run_command_with_timeout('./follow2 << EOF \nparam.in \npl.in \n-%s \nEOF'%(k+2),timeout_sec)
            result, flag = run_command_with_timeout('mv follow2.out pl_%s.out'%(k+1),timeout_sec)

        obj.evol_T[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [0]) /  365.25
        obj.evol_a[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [2])
        obj.evol_e[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [3])
        obj.evol_p[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [6])
        obj.evol_M[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [7])

        obj.evol_i[k]  = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [4])
        obj.evol_Om[k] = np.genfromtxt("pl_%s.out"%(k+1),skip_header=0, unpack=True,skip_footer=1, usecols = [5])

        obj.evol_Per[k] = a_to_P(obj.evol_a[k],obj.params.stellar_mass)


    try:
        os.system('rm *.out *.dat *.in')
        #os.system('mv *.out *.dat *.in last_run')
    except OSError:
        pass

    os.chdir('../../')

    print("stability with: %s done!"%integrator)

    return obj



def run_copl_fit_stab(obj, incl_max=90.0, incl_min=90.0, incl_step = 1.0, save_output=True, output_file="./copl_incl.txt", fit_bf = False,
 timemax=3000.0, timestep=10, timeout_sec=1000.0, stab_save_dir = './run', remove_stab_save_dir = True, integrator='symba',a_threshold =10, e_max =0.9):

    """So far only RVs can be fitted!!!"""

    incl_fit = dill.copy(obj)
    incl_fit.mod_dynamical=True

    if save_output == True:
        f = open(output_file,"w")


    incl_range = np.arange(incl_max,incl_min,-incl_step)

    for incl in incl_range:

        for i in range(incl_fit.npl):
            incl_fit.params.planet_params[7*i+5] = incl
            incl_fit.use.use_planet_params[7*i+5] = False

            if fit_bf:
                incl_fit.use.update_use_planet_params_one_planet(i,True,True,True,True,True,False,False)
            else:
                incl_fit.use.update_use_planet_params_one_planet(i,False,False,False,False,False,False,False)


        incl_fit.fitting(outputfiles=[1,1,1], doGP=False,  minimize_fortran=True, minimize_loglik=False, amoeba_starts=0, print_stat=False)
        incl_fit.fitting(outputfiles=[1,1,1], doGP=False,  minimize_fortran=True, minimize_loglik=True, amoeba_starts=0, print_stat=False)
        incl_fit.fitting(outputfiles=[1,1,1], doGP=False,  minimize_fortran=True, minimize_loglik=True, amoeba_starts=10, print_stat=False)

        run_stability(incl_fit, timemax=timemax, timestep=timestep, timeout_sec=timeout_sec, stab_save_dir = stab_save_dir, remove_stab_save_dir = remove_stab_save_dir, integrator=integrator)

        for i in range(incl_fit.npl):
            export_orbital_evol(incl_fit, file='planet_%s_%s.txt'%(i,incl), planet = i+1, width = 10, precision = 6)


        stab_amd = int(get_AMD_stab(incl_fit))


        stab = 1
        for i in range(incl_fit.npl):
            if max(incl_fit.evol_e[i]) > e_max:
                stab = 0


        print("%s   %s "%(incl,incl_fit.loglik))

        if save_output == True:
            f.write("%s"%incl_fit.loglik)
            for i in range(incl_fit.npl):
                for z in range(7):
                    f.write("%s  " %(incl_fit.params.planet_params[7*i+z]))
            f.write("%s %s\n"%(stab,stab_amd))

    if save_output == True:
        f.close()

    return obj





