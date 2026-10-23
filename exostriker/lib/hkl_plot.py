"""Derived Keplerian posterior columns from h=e sin(w), k=e cos(w), lambda."""
import numpy as np


def derived_hkl_columns(fit, samples, labels, options):
    if not getattr(fit, 'hkl', False):
        return []
    result=[]
    labels=list(labels)
    for slot, active in enumerate(fit.use_planet):
        if not active:
            continue
        letter='bcdefghij'[slot]
        def values(root):
            name=getattr(fit,root+'_str')[slot]
            return np.asarray(samples[:,labels.index(name)],dtype=float) if name in labels else np.full(len(samples),getattr(fit,root)[slot],dtype=float)
        h,k,longitude=values('e_sinw'),values('e_cosw'),values('lamb')
        eccentricity=np.hypot(h,k)
        omega=np.degrees(np.arctan2(h,k)) % 360.
        # A circular orbit has no periastron. For plotting M, adopt the
        # omega=0 reference convention, hence M=lambda, while omega stays NaN.
        anomaly=(longitude-np.where(eccentricity>0,omega,0.)) % 360.
        omega=np.where(eccentricity>0,omega,np.nan)
        for option,label,draws in [('derive_e',f'e$_{letter}$',eccentricity),
                                  ('derive_omega',rf'$\omega_{letter}$ [deg]',omega),
                                  ('derive_M0',f'MA$_{letter}$ [deg]',anomaly)]:
            if options.get(option,False) and label not in labels:
                result.append((option,slot,label,draws))
    return result
