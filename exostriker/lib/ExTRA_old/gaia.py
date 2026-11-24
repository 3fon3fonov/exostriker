import numpy as np
from astropy.time import Time

from .astrometry import *
from .hipparcos import *
from .useful import *

__all__ = ["angle_trafo", "gaia_init", "gaia_JD", "res_to_orbit_gaia", "decode_gaia_flag", "decode_transit_proc_flags"]

#### Psi=pi/2 - theta #psi=HIP2 theta=GAIADR4

def angle_trafo(theta):
    psi=-theta+np.pi/2

    psi=((psi + 180) % 360) - 180
    return psi


def gaia_init(data,standardepoch="2017.5"):
    """transforms astrometric gaia dr4 data into HIP format
    Input:
    ----------
    data    format(5,N)    [time,measurement,err,parallaxfactor,angle]
    ----------
    Keys:
    ----------
    standardepoch   str     the epoch of the input gaia measurement, typically always 2017.5
    ----------
    Output:
    ----------
    data    format(7,N)    [cos,sin,parallaxfactor,cos*t,sin*t,measurement,err]
    """

    GAIA_EPOCH=Time(standardepoch, format='jyear',scale="tcb")

    

    relative_time=Time(data[0], format='jd', scale='tcb').jyear-GAIA_EPOCH.jyear
    gaia_angle=angle_trafo(np.radians(data[4]))
    #gaia_angle=angle_trafo(data[4])

    t_gaia=data[0] #gets returned seperatly, not important for fitting
    
    A2=relative_time
    A3=np.cos(gaia_angle)
    A4=np.sin(gaia_angle)
    A5=data[3]
    A6=A3*A2
    A7=A4*A2
    A8=data[1]
    A9=data[2]
    transformed=[A3,A4,A5,A6,A7,A8,A9]


    return transformed,t_gaia

def gaia_JD(gaia_ad,format="jd",Sepoch=None):
    if Sepoch==None:
        Sepoch=J2017()

    Sepoch_jyear=Time(Sepoch,format="jd").jyear

    A3,A4,A5,A6,A7,A8,A9=gaia_ad
    frac=A7/A4
    epoch=frac+Sepoch_jyear
    if format=="jd":
        JD=Sepoch+(epoch-Sepoch_jyear)*365.25 #JD for standard epoch J2017.5
        return JD
    if format=="jyear":
        return epoch
    if format=="relative":
        return frac


def res_to_orbit_gaia(residuals,gaia_ad,orbitfit):
    
     
    #gaia

    t=gaia_JD(gaia_ad)
    res=gaia_ad.copy()
    res[-2]=residuals

    

    res_2D=hip_2d(res)

    hip_x=res_2D[0]
    hip_x_err=res_2D[1]
    hip_y=res_2D[2]
    hip_y_err=res_2D[3]

    #print(res_2D)
    orb_x,orb_y=orbit_total(orbitfit,t)

    res_orb_x=hip_x+orb_x
    res_orb_y=hip_y+orb_y

    return res_orb_x,hip_x_err,res_orb_y,hip_y_err


def decode_gaia_flag(value):
    value = int(value)
    # --- IPD result (bits 0–3) ---
    ipd_codes = {
        0: "SUCCESS",
        1: "ILLEGAL_AMPLITUDE",
        2: "ILLEGAL_COORDINATE",
        3: "DEPRECATED_ILLEGAL_OUTPUT",
        4: "MAX_ITER_EXCEEDED",
        5: "FAILED",
        6: "PENDING",
        7: "NO_LSF_PSF",
        8: "NOT_ENOUGH_SAMPLES",
        9: "ODD_LSF_PSF_SUCCESS",
        14: "NOT_ATTEMPTED",
        15: "CRASHED",
    }

    ipd_result_code = value & 0x000F
    ipd_result = ipd_codes.get(ipd_result_code, f"UNKNOWN ({ipd_result_code})")

    # --- Flags (bits 4–15) ---
    flag_defs = {
        4:  "non_nominal_ipd",
        5:  "ipd_not_available",
        6:  "cosmetic_issue",
        7:  "cosmic_removed",
        8:  "saturation_removed",
        9:  "window_part_discarded",
        10: "no_window",
        11: "bad_initial_centroid",
        12: "odd_background",
        13: "fallback_bias_mitigation",
        14: "non_target_removed",
        15: "psf_parameters_clamped",
    }

    active_flags = [
        name for bit, name in flag_defs.items()
        if value & (1 << bit)
    ]
    flag_numbers = np.array([
        bit for bit, name in flag_defs.items()
        if value & (1 << bit)
    ])

    return {
        "value": value,
        "ipd_result_code": ipd_result_code,
        "ipd_result": ipd_result,
        "flags": active_flags,
        "flag_numbers": flag_numbers
    }

def decode_transit_proc_flags(value):
    """
    Decode Gaia epoch_astrometry transit_proc_flags.

    Parameters
    ----------
    value : int
        Integer transit_proc_flags value.

    Returns
    -------
    dict
        Dictionary containing the decoded fields.
    """
    value = int(value)

    # ----- Bits 0-3 : IPD strategy -----
    ipd_strategy = {
        0: "Normal IPD (no ambiguity detected or tested)",
        1: "Multiple Source IPD: earlier object selected (or higher μ for large AC separation)",
        2: "Multiple Source IPD: later object selected (or lower μ for large AC separation)",
        3: "Multiple Source IPD: disturbing sources masked during IPD",
        4: "DIPD solution (window model template + combined AF solution)",
        5: "TukeyBiWeight IPD initialisation (internal IDU-AttAc only)",
    }

    strategy_code = value & 0x000F
    strategy = ipd_strategy.get(
        strategy_code,
        f"Reserved/unknown ({strategy_code})"
    )

    # ----- Bits 4-6 : Bias strategy -----
    bias_strategy = {
        0: "UNDEFINED / failure",
        1: "PreScan bias only",
        2: "PreScan + common NU baseline",
        3: "Full PEM NU",
    }

    bias_code = (value >> 4) & 0x7
    bias = bias_strategy.get(
        bias_code,
        f"Reserved ({bias_code})"
    )

    # ----- Individual flags -----
    flags = []

    if value & 0x0080:
        flags.append("Default colour used")

    if value & 0x0100:
        flags.append("Epoch (transit-level) colour used")
    else:
        flags.append("Mean (source-level) colour used")

    if value & 0x0200:
        flags.append("Single-parameter colour")
    else:
        flags.append("Multi-parameter colour")

    if value & 0x0400:
        flags.append("Colour missing or dubious")

    if value & 0x0800:
        flags.append("Dead pixel columns affect some windows")

    if value & 0x1000:
        flags.append("Anomalously low samples present")

    if value & 0x2000:
        flags.append("Anomalously high samples present")

    if value & 0x4000:
        flags.append("Problems in initial centroid")

    if value & 0x8000:
        flags.append("Attitude issues affecting some windows")

    return {
        "value": value,
        "ipd_strategy_code": strategy_code,
        "ipd_strategy": strategy,
        "bias_strategy_code": bias_code,
        "bias_strategy": bias,
        "flags": flags,
    }