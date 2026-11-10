"""Independent TTVFast + batman local-chord transit backup.

No native transit extension or REBOUND dependency. Event positions and speeds
are converted to the RVmod physical length convention before occultation.
The local straight-line path is an approximation, not a full photodynamical
trajectory. Overlapping planet curves are multiplied, as in standard batman use.
"""
import numpy as np
import batman
import ttvfast

GM_SUN=1.32712440018e20
AU=1.49597892e11
RSUN=6.957e8
DAY=86400.
G_DENSITY=6.67408e-11


def physical_radius(stellar_mass, stellar_radius=1., rho_star=None):
    if not np.isfinite(stellar_mass) or stellar_mass<=0:
        raise ValueError('Stellar mass must be finite and positive.')
    if rho_star is None:
        radius=float(stellar_radius)*RSUN/AU
    else:
        if not np.isfinite(rho_star) or rho_star<=0:
            raise ValueError('Stellar density must be finite and positive.')
        radius=(3*stellar_mass*GM_SUN/(4*np.pi*G_DENSITY*rho_star))**(1/3)/AU
    if not np.isfinite(radius) or radius<=0:
        raise ValueError('Stellar radius must be finite and positive.')
    return radius


def predict_events(planets, masses, stellar_mass, epoch, end, step):
    """Return (number, t0 JD, rsky AU, vsky AU/day), original planet order.

    TTVFast evolves relative times to avoid secular rounding of absolute JD.
    Its public wrapper's G sets its length scale, which is normalized here to
    the GM_SUN/AU/DAY convention used by RVmod. No orbit/phase rescaling occurs.
    """
    rows=np.asarray(planets,float).reshape(-1,7);masses=np.asarray(masses,float)
    if masses.shape!=(len(rows),):raise ValueError('One physical mass per planet is required.')
    if np.any(~np.isfinite(rows)) or np.any(~np.isfinite(masses)) or np.any(masses<0):
        raise ValueError('Require finite orbital elements and nonnegative masses.')
    if np.any(rows[:,1]<=0) or np.any(rows[:,2]<0) or np.any(rows[:,2]>=1):
        raise ValueError('Require positive periods and eccentricities in [0,1).')
    if not np.isfinite(stellar_mass) or stellar_mass<=0:raise ValueError('Stellar mass must be positive.')
    if not np.isfinite(epoch) or not np.isfinite(end) or end<=epoch:
        raise ValueError('Require a finite TTVFast end time greater than epoch.')
    if not np.isfinite(step) or step<=0:raise ValueError('TTVFast step must be positive.')
    if not len(rows):return []
    order=np.argsort(rows[:,1],kind='stable')
    bodies=[ttvfast.models.Planet(masses[j],rows[j,1],rows[j,2],rows[j,5],rows[j,6],rows[j,3],rows[j,4]) for j in order]
    source_g=float(ttvfast.models.planets_to_params(stellar_mass,bodies)[0])
    length_scale=(GM_SUN*DAY**2/AU**3/source_g)**(1/3)
    output=np.asarray(ttvfast.ttvfast(bodies,float(stellar_mass),0.,float(step),float(end-epoch),input_flag=0)['positions']).T
    if output.size==0:return [np.empty((0,4)) for _ in rows]
    output=output[np.isfinite(output[:,2])&(output[:,2]>=0)]
    events=[None]*len(rows)
    for k,j in enumerate(order):
        selected=output[output[:,0]==k]
        table=selected[:,1:5].copy()
        table[:,1]+=epoch;table[:,2:4]*=length_scale
        events[j]=table
    return events


def occultation_from_separation(distance, radius_ratio, limb_dark, coefficients):
    """Batman's existing occultation laws, with externally supplied separation."""
    z=np.ascontiguousarray(distance,float)
    if z.ndim!=1 or np.any(~np.isfinite(z)) or np.any(z<0):
        raise ValueError('Projected separations must be finite, nonnegative and one-dimensional.')
    if not np.isfinite(radius_ratio) or radius_ratio<0:
        raise ValueError('Radius ratio must be finite and nonnegative.')
    if not len(z) or radius_ratio==0:return np.ones(len(z))
    p=batman.TransitParams();p.t0=0.;p.per=3.;p.a=10.;p.inc=90.;p.ecc=0.;p.w=90.
    p.rp=float(radius_ratio);p.limb_dark=str(limb_dark);p.u=list(coefficients)
    model=batman.TransitModel(p,np.zeros(len(z)),max_err=.001)
    # Do not use batman's fixed Keplerian trajectory. light_curve keeps ds when
    # the orbital parameters have not changed since initialization.
    model.ds=z
    return model.light_curve(p)


def event_light_curve(times, events, radius_ratio, stellar_radius_au, limb_dark, coefficients):
    """Instantaneous isolated-planet flux using t0, rsky AND vsky.

    z(t)^2 = b^2 + [v_perp*(t-t0)]^2. Event windows follow the actual
    first/fourth-contact estimate, with no hard-coded 0.3-day truncation.
    """
    if not np.isfinite(radius_ratio) or radius_ratio<0:
        raise ValueError('Radius ratio must be finite and nonnegative.')
    times=np.asarray(times,float)
    if times.ndim!=1 or np.any(~np.isfinite(times)):raise ValueError('Times must be finite and one-dimensional.')
    if not np.isfinite(stellar_radius_au) or stellar_radius_au<=0:
        raise ValueError('Stellar radius must be positive.')
    flux=np.ones(len(times))
    for _,t0,rsky,vsky in np.asarray(events,float).reshape(-1,4):
        if not np.all(np.isfinite([t0,rsky,vsky])) or rsky<0 or vsky<=0:
            raise ValueError('Invalid TTVFast event geometry.')
        b=rsky/stellar_radius_au;speed=vsky/stellar_radius_au
        if b>=1+radius_ratio:continue
        half=np.sqrt((1+radius_ratio)**2-b*b)/speed
        dt=times-t0;mask=abs(dt)<=half
        if not np.any(mask):continue
        separation=np.hypot(b,speed*dt[mask])
        flux[mask]*=occultation_from_separation(separation,radius_ratio,limb_dark,coefficients)
    return flux
