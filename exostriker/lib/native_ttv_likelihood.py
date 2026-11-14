"""Independent timing-only TTV likelihood and legacy GUI result adapter.

No TTVFast, batman, SciPy optimizer, or Python event search is used. The single
native entry integrates, finds conjunctions, matches event numbers and scores
all enabled timing datasets. Depth/duration likelihoods remain a legacy option.
"""
from dataclasses import dataclass, replace
import numpy as np
from native_transit_model import _native, _LOCK, _mass_axes


@dataclass(frozen=True)
class NativeTTVLikelihoodResult:
    dataset: np.ndarray
    planet: np.ndarray                 # original zero-based input row
    number: np.ndarray                 # one-based event since search start
    observed: np.ndarray
    uncertainty: np.ndarray
    data_model: np.ndarray
    residual: np.ndarray
    point_loglik: np.ndarray
    by_dataset: np.ndarray             # ten original dataset slots
    loglik: float
    events: dict                      # original planet row -> [N, JD, AU, AU/day]
    legacy_model: list                 # the existing five-section GUI structure
    planet_slots: np.ndarray           # maps input row -> public GUI slot


def _integers(values, shape, name, lower, upper):
    raw=np.broadcast_to(np.asarray(values,dtype=float),shape)
    if not np.all(np.isfinite(raw)) or np.any(raw!=np.floor(raw)) or np.any(raw<lower) or np.any(raw>upper):
        raise ValueError(f'{name} must contain integers in {lower}..{upper}.')
    return raw.astype(np.int32)


def evaluate_ttv_likelihood(planets, number, observed, uncertainty, *, planet=0, dataset=0,
                            stellar_mass=1., epoch=0., start=None, end=None,
                            step=.1, tolerance=1e-14, jacobi_convention='ttvfast',
                            samples_per_orbit=512, max_samples=1000000,
                            max_events=None, max_steps=1000000):
    """Score up to 10 timing datasets and 9 planets, using one native call.

    Planet rows are K,P,e,omega,M0,inclination,Omega, with M0 at epoch.
    Event 1 is the first foreground sky-separation minimum inside [start,end],
    independently for each planet, even if its impact parameter misses the star.
    Observed uncertainties and times are days; output t0 is absolute JD. The
    event table returns physical AU and AU/day, not TTVFast's slightly different
    G-dependent length convention. No transit offset, jitter, LD or LC input is
    needed for timing-only data. Dataset/planet/sample ordering is preserved.
    """
    rows=np.asarray(planets,float)
    if rows.ndim!=2 or rows.shape[1]!=7 or not 1<=len(rows)<=9:
        raise ValueError('Supply 1..9 rows of K,P,e,omega,M0,inclination,Omega.')
    if not np.all(np.isfinite(rows)) or np.any(rows[:,1]<=0) or np.any(rows[:,2]<0) or np.any(rows[:,2]>=1):
        raise ValueError('Require finite elements, positive periods and eccentricities in [0,1).')
    observed=np.asarray(observed,float);uncertainty=np.asarray(uncertainty,float)
    if observed.ndim!=1 or not len(observed) or uncertainty.shape!=observed.shape:
        raise ValueError('Supply nonempty, matching one-dimensional times and errors.')
    if not np.all(np.isfinite(observed)) or not np.all(np.isfinite(uncertainty)) or np.any(uncertainty<=0):
        raise ValueError('TTV times must be finite and uncertainties strictly positive.')
    labels=_integers(dataset,observed.shape,'dataset',0,9)
    bodies=_integers(planet,observed.shape,'planet',0,len(rows)-1)
    numbers=_integers(number,observed.shape,'transit number',1,np.iinfo(np.int32).max)
    for d in np.unique(labels):
        if len(np.unique(bodies[labels==d]))!=1:
            raise ValueError('Each timing dataset must refer to exactly one planet.')
    if start is None:start=epoch
    if end is None:raise ValueError('Supply the end of the TTV event-search interval.')
    if not np.all(np.isfinite([stellar_mass,epoch,start,end,step,tolerance])) or stellar_mass<=0 or end<=start or step<=0:
        raise ValueError('Require positive stellar mass/step and finite increasing search bounds.')
    if not 1e-15<=tolerance<=1e-3:raise ValueError('tolerance must lie in 1e-15..1e-3.')
    for name,value,low,high in [('samples_per_orbit',samples_per_orbit,32,1000000),
                              ('max_samples',max_samples,4,2000000),('max_steps',max_steps,1,np.iinfo(np.int32).max)]:
        if not np.isscalar(value) or not np.isfinite(value) or value!=int(value) or not low<=value<=high:
            raise ValueError(f'{name} must be an integer in {low}..{high}.')
    if jacobi_convention not in ('standard','ttvfast'):
        raise ValueError('jacobi_convention must be standard or ttvfast.')
    estimate=2*(end-start)/np.min(rows[:,1])+32
    if max_events is None:
        if not np.isfinite(estimate) or estimate>200000:
            raise ValueError('Event search is too long; narrow the interval.')
        max_events=max(32,int(np.ceil(estimate)))
    if not np.isscalar(max_events) or not np.isfinite(max_events) or max_events!=int(max_events) or not 1<=max_events<=200000:
        raise ValueError('max_events must be an integer in 1..200000.')
    masses,_=_mass_axes(rows,stellar_mass)
    order=np.argsort(rows[:,1],kind='stable');inverse=np.argsort(order)
    with _LOCK:
        output=_native().ttv_likelihood(rows[order,1:],masses[order],float(stellar_mass),float(epoch),
            float(start),float(end),float(step),float(tolerance),int(max_steps),int(jacobi_convention=='ttvfast'),
            int(samples_per_orbit),int(max_samples),int(max_events),labels+1,inverse[bodies]+1,numbers,
            observed,uncertainty)
    counts,table,model,residual,points,by_dataset,loglik,status=output
    if status:
        messages={1:'Invalid native TTV geometry.',2:'TTV Kepler solver failed.',
                  4:'TTV integration exceeded max_steps.',5:'TTV integration failed.',
                  7:'Invalid native TTV inputs.',8:'A requested transit number is outside the computed event sequence.',
                  9:'TTV event-search sample/memory limit exceeded.',10:'TTV event output capacity exceeded.'}
        raise ValueError(messages.get(status,f'Native TTV evaluation failed (status {status}).'))
    events={int(original):table[:counts[j],j,:].copy() for j,original in enumerate(order)}
    calc_data={d:[] for d in range(10)};calc_model={d:[] for d in range(10)}
    short=[[],[]];full=[[],[],[],[],[],[]]
    for d in range(10):
        mask=labels==d
        if not np.any(mask):continue
        body=int(bodies[mask][0]);event=events[body];selected=event[numbers[mask]-1]
        calc_data[d]=[numbers[mask].copy(),model[mask].copy(),uncertainty[mask].copy(),
                      selected[:,2].copy(),selected[:,3].copy(),[],[]]
        calc_model[d]=[event[:,0].astype(int),event[:,1].copy(),event[:,2].copy(),event[:,3].copy(),
                       np.ones(len(event)),np.zeros(len(event))]
        short=[numbers[mask].copy(),model[mask].copy()];full=calc_model[d]
    legacy=[float(loglik),short,full,calc_data,calc_model]
    return NativeTTVLikelihoodResult(labels.copy(),bodies.copy(),numbers.copy(),observed.copy(),uncertainty.copy(),
        model,residual,points,by_dataset,float(loglik),events,legacy,np.arange(len(rows)))


def from_session_ttv_likelihood(session, *, store=True, publish=False, **options):
    """Read enabled legacy TTV slots, preserving their one-based numbering.

    The TTV-only legacy path places input elements at ttv_times[0], which is
    therefore the default phase epoch here too. Explicit epoch/start/end options
    can override that convention. GUI planet assignments retain public slots.
    """
    targets=[session.fit_results]
    if getattr(session,'rvmod',None) is not None and session.rvmod is not session.fit_results:
        targets.append(session.rvmod)
    if store:
        for target in targets:
            target.native_ttv_likelihood=None;target.native_ttv_loglik=None
            target.native_ttv_results=None;target.native_ttv_likelihood_status='not calculated'
    if publish:session.ttv_results=None
    try:
        slots=list(np.flatnonzero(session.use_planet))
        if not slots:raise ValueError('No active planets for native TTV evaluation.')
        from mass_parameterization import mass_model_session
        model_session = mass_model_session(session)
        rows=np.array([[getattr(model_session,name)[i] for name in ('K','P','e','w','M0','i','Node')] for i in slots])
        if getattr(session,'hkl',False):
            h=np.array([session.e_sinw[i] for i in slots]);k=np.array([session.e_cosw[i] for i in slots])
            rows[:,2]=np.hypot(h,k);rows[:,3]=np.degrees(np.arctan2(h,k))
            rows[:,4]=(np.array([session.lamb[i] for i in slots])-rows[:,3])%360.
        if getattr(session,'gr_flag',False):raise ValueError('Native TTV evaluation does not support GR.')
        if any(session.omega_dot[i]!=0 for i in slots):raise ValueError('Native TTV evaluation does not support precession.')
        mixed=getattr(session,'mixed_fit',{0:False})
        if bool(np.asarray(mixed[0]).item()):raise ValueError('Native TTV evaluation does not support mixed dynamical subsets.')
        # The legacy dedicated TTV model integrates all planets regardless of
        # the RV Keplerian/dynamical toggle. Keep that behavior in this adapter.
        numbers=[];observed=[];errors=[];datasets=[];planets=[]
        for d in range(10):
            data=session.ttv_data_sets[d]
            if not len(data) or not data[4]:continue
            if len(data)>11 and (data[-2] or data[-1]):
                raise ValueError('This native TTV entry is timing-only; use legacy for depth/duration likelihoods.')
            raw=float(data[3])
            if not np.isfinite(raw) or raw!=int(raw) or int(raw)-1 not in slots:
                raise ValueError(f'TTV dataset {d} refers to an inactive or invalid planet slot.')
            t=np.asarray(data[1],float)
            if not len(t):raise ValueError(f'Enabled TTV dataset {d} is empty.')
            if np.shape(data[0])!=t.shape or np.shape(data[2])!=t.shape:
                raise ValueError(f'TTV dataset {d} numbers/times/errors have different shapes.')
            numbers.append(np.asarray(data[0]));observed.append(t);errors.append(np.asarray(data[2]))
            datasets.append(np.full(len(t),d));planets.append(np.full(len(t),slots.index(int(raw)-1)))
        if not observed:raise ValueError('No enabled TTV observations are loaded.')
        convention=getattr(session,'native_jacobi_convention','auto')
        if convention=='auto':convention='ttvfast'
        defaults=dict(stellar_mass=session.params.stellar_mass,epoch=session.ttv_times[0],
                      start=session.ttv_times[0],end=session.ttv_times[2],step=session.ttv_times[1],
                      tolerance=getattr(session,'native_ttv_tolerance',1e-14),
                      jacobi_convention=convention)
        defaults.update(options)
        result=evaluate_ttv_likelihood(rows,np.concatenate(numbers),np.concatenate(observed),np.concatenate(errors),
            dataset=np.concatenate(datasets),planet=np.concatenate(planets),**defaults)
        result=replace(result,planet_slots=np.asarray(slots,dtype=int))
    except Exception:
        if store:
            for target in targets:target.native_ttv_likelihood_status='failed'
        raise
    if store:
        for target in targets:
            target.native_ttv_likelihood=result;target.native_ttv_loglik=result.loglik
            target.native_ttv_results=result.legacy_model;target.native_ttv_likelihood_status='ready'
    if publish:session.ttv_results=result.legacy_model
    return result
