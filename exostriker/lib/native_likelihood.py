"""Native Gaussian RV + transit likelihood, with no optimization or GP fallback."""
from dataclasses import dataclass
import numpy as np
from native_transit_model import TransitRVModel, evaluate, from_session, _native, _LOCK


@dataclass(frozen=True)
class NativeLikelihoodResult(TransitRVModel):
    kind: np.ndarray                 # 0 RV, 1 transit
    observed: np.ndarray
    uncertainty: np.ndarray
    covariate: np.ndarray
    data_model: np.ndarray           # prediction in the observed-data frame
    residual: np.ndarray
    variance: np.ndarray
    point_loglik: np.ndarray
    by_dataset: np.ndarray           # [zero-based dataset slot, kind]
    loglik: float
    rv_loglik: float
    transit_loglik: float
    chi2: np.ndarray                 # [RV, transit]


def _table(value, shape, name, fill=0.):
    if value is None:
        return np.full(shape,fill,dtype=float)
    if np.shape(value)!=shape:
        raise ValueError(f'{name} must have shape {shape}.')
    out=np.asarray(value,float)
    if not np.all(np.isfinite(out)):
        raise ValueError(f'{name} must be finite.')
    return out


def evaluate_likelihood(planets, times, observed, uncertainty, *, kind='transit', dataset=0,
                        covariate=None, baseline=None, rv_offsets=None, rv_jitters=None,
                        rv_trend=None, **model_options):
    """One Fortran call returns physical models and Gaussian RV/transit lnL.

    kind is 'rv'/'transit' or integer 0/1, scalar or per observation. Dataset
    slots are 0..59 independently for RV and transit. Baseline is (60,3):
    linear coefficient, quadratic coefficient, dilution. Covariate is elapsed
    time from the dataset reference or raw airmass (caller supplies it).
    Transit offsets/jitters and exposure/LD options use evaluate's existing API.
    Empty datasets are rejected; zero-planet baseline-only models are allowed.
    """
    times=np.asarray(times,float)
    if times.ndim!=1 or not len(times):
        raise ValueError('Supply a nonempty one-dimensional observation-time array.')
    if np.shape(observed)!=times.shape or np.shape(uncertainty)!=times.shape:
        raise ValueError('Observed values and uncertainties must match times.')
    observed=np.asarray(observed,float);uncertainty=np.asarray(uncertainty,float)
    if not np.all(np.isfinite(observed)) or not np.all(np.isfinite(uncertainty)) or np.any(uncertainty<0):
        raise ValueError('Observations must be finite and uncertainties nonnegative.')
    raw=np.broadcast_to(np.asarray(kind),times.shape)
    if raw.dtype.kind in 'USO':
        if not np.all(np.isin(raw,['rv','transit'])):
            raise ValueError("kind must be 'rv'/'transit' or 0/1.")
        kinds=(raw=='transit').astype(np.int32)
    else:
        if not np.all(np.isin(raw,[0,1])):
            raise ValueError('kind must contain only 0 (RV) or 1 (transit).')
        kinds=raw.astype(np.int32)
    cov=_table(covariate,times.shape,'covariate')
    base=_table(baseline,(60,3),'baseline')
    if baseline is None:base[:,2]=1.
    rv_sets=np.column_stack((_table(rv_offsets,(60,),'rv_offsets'),_table(rv_jitters,(60,),'rv_jitters')))
    trend=_table(rv_trend,(2,),'rv_trend')
    if np.any(rv_sets[:,1]<0) or np.any(base[:,2]<0):
        raise ValueError('RV jitter and dilution must be nonnegative.')
    payload=(observed,uncertainty,kinds,cov,base,rv_sets,trend)
    return evaluate(planets,times,dataset=dataset,_observations=payload,**model_options)


def _evaluate_prepared(planets,times,labels,rows,axis,masses,stellar_mass,stellar_radius,epoch,
                       laws,sets,exposure,samples,nquad,model,step,tolerance,max_steps,phase,
                       transit_input,convention,payload):
    observed,uncertainty,kinds,cov,base,rv_sets,trend=payload
    if not np.all(np.isfinite(exposure)) or np.any(exposure<0):
        raise ValueError('Exposure durations must be finite and nonnegative.')
    # RV measurements are instantaneous; only transit records are exposure-averaged.
    counts=np.where(kinds==1,samples[labels],1)
    owners=np.repeat(np.arange(len(times)),counts)
    sample_index=np.arange(len(owners))-np.repeat(np.cumsum(counts)-counts,counts)
    fractions=np.zeros(len(owners));mask=counts[owners]>1
    fractions[mask]=sample_index[mask]/(counts[owners[mask]]-1)-.5
    elapsed=times-float(epoch)
    expanded=elapsed[owners]+exposure[labels[owners]]*fractions
    order=np.argsort(expanded,kind='stable')
    planet_order=np.argsort(planets[:,1],kind='stable') if model=='nbody' else np.arange(len(planets))
    physical_mass=np.zeros(9);physical_mass[:len(planets)]=masses[planet_order]
    t0=np.zeros(9);t0[:len(planets)]=transit_input[planet_order]
    with _LOCK:
        output=_native().rv_transit_likelihood(expanded[order],(owners[order]+1).astype(np.int32),elapsed,
            labels+1,kinds+1,observed,uncertainty,cov,rows[planet_order],physical_mass,
            float(stellar_mass),float(stellar_radius),float(epoch),int(model=='nbody'),
            int(convention=='ttvfast'),int(phase=='t0'),t0,laws,sets,rv_sets,base,trend,nquad,
            float(step or 1.),float(tolerance),max_steps)
    (rv,flux,data_model,residual,variance,point_loglik,by_dataset,components,loglik,chi2,
     radii,rv_parts,mean_anomaly,transit_epoch,status)=output
    if status:
        messages={1:'Invalid physical model or limb-darkening profile.',2:'Kepler solver did not converge.',
                  4:'N-body integration exceeded max_steps.',5:'N-body integration failed.',
                  6:'Invalid likelihood data, variance, or nonfinite model/score.'}
        raise ValueError(messages.get(status,f'Native likelihood failed (status {status}).'))
    restore=np.argsort(planet_order)
    jitter=np.where(kinds==0,rv_sets[labels,1],sets[labels,1])
    return NativeLikelihoodResult(times.copy(),labels.copy(),rv,flux,flux-sets[labels,0],jitter,
        axis,radii[:len(planets)][restore],rv_parts[:,restore],mean_anomaly[:len(planets)][restore],
        transit_epoch[:len(planets)][restore] if model=='keplerian' else None,
        kinds.copy(),observed.copy(),uncertainty.copy(),cov.copy(),data_model,residual,variance,
        point_loglik,by_dataset,float(loglik),float(components[0]),float(components[1]),chi2)


def from_session_likelihood(session, *, include_rv=None, include_transits=None, phase='M0',
                            store=True, **options):
    """Reevaluate loaded RV/LC observations; preserve legacy scores and parameters.

    Defaults follow type_fit RV/Transit flags. No GP, separate TTV/TDV data,
    astrometry, priors or stability penalty is included in this Gaussian score.
    GP-enabled selected components raise rather than silently dropping a GP.
    """
    targets=[session.fit_results]
    if getattr(session,'rvmod',None) is not None and session.rvmod is not session.fit_results:
        targets.append(session.rvmod)
    if store:
        for target in targets:
            target.native_likelihood=None;target.native_loglik=None
            target.native_rv_loglik=None;target.native_transit_loglik=None
            target.native_likelihood_status='not calculated'
    try:
        fit_types=getattr(session,'type_fit',{})
        include_rv=fit_types.get('RV',True) if include_rv is None else bool(include_rv)
        include_transits=fit_types.get('Transit',True) if include_transits is None else bool(include_transits)
        rtg=getattr(session,'rtg',[False]*4)
        if include_rv and len(rtg)>1 and rtg[1] and any(len(d) for d in getattr(session,'rv_data_sets',{}).values()):
            raise ValueError('Native Gaussian likelihood does not support RV GPs yet.')
        if include_transits and len(rtg)>3 and rtg[3] and any(
                len(data) and data[9] for data in session.tra_data_sets.values()):
            raise ValueError('Native Gaussian likelihood does not support transit GPs yet.')
        if include_transits and getattr(session,'get_TTVs',[False])[0]:
            raise ValueError('Native physical likelihood does not support independent per-event transit-time shifts yet.')
        times=[];values=[];errors=[];labels=[];kinds=[];covariates=[]
        rv_slots=[]
        for kind,name,selected in ((0,'rv_data_sets',include_rv),(1,'tra_data_sets',include_transits)):
            if not selected:continue
            for d,data in getattr(session,name,{}).items():
                if not len(data):continue
                t=np.asarray(data[0],float)
                if not len(t):continue
                times.append(t);values.append(np.asarray(data[1],float));errors.append(np.asarray(data[2],float))
                labels.append(np.full(len(t),d));kinds.append(np.full(len(t),kind))
                if kind==1:
                    covariates.append(np.asarray(data[3],float) if len(data)>10 and data[10] else t-t[0])
                else:
                    covariates.append(np.zeros(len(t)));rv_slots.append(d)
        if not times:raise ValueError('No observations are loaded for the selected likelihood components.')
        base=np.zeros((60,3));base[:,2]=1.
        for d in range(60):
            data=getattr(session,'tra_data_sets',{}).get(d,[])
            default_dilution=data[8] if len(data)>8 else 1.
            base[d]=[getattr(session,'tra_lintr',{}).get(d,0.),getattr(session,'tra_quadtr',{}).get(d,0.),
                     getattr(session,'tra_dil',{}).get(d,default_dilution)]
        offsets=np.array([getattr(session,'rvoff',{}).get(d,0.) for d in range(60)])
        jitter=np.array([getattr(session,'jitt',{}).get(d,0.) for d in range(60)])
        if getattr(session,'jit_flag',0) and rv_slots:
            jitter[:]=jitter[min(rv_slots)]
        extra=dict(baseline=base,rv_offsets=offsets,rv_jitters=jitter,
                   rv_trend=[getattr(session,'rv_lintr',0.),getattr(session,'rv_quadtr',0.)])
        extra.update(options)
        def evaluator(planets,times,**model_options):
            return evaluate_likelihood(planets,times,np.concatenate(values),np.concatenate(errors),
                kind=np.concatenate(kinds),covariate=np.concatenate(covariates),**model_options)
        result=from_session(session,np.concatenate(times),dataset=np.concatenate(labels),phase=phase,
                            _evaluator=evaluator,**extra)
    except Exception:
        if store:
            for target in targets:target.native_likelihood_status='failed'
        raise
    if store:
        for target in targets:
            target.native_likelihood=result;target.native_loglik=result.loglik
            target.native_rv_loglik=result.rv_loglik;target.native_transit_loglik=result.transit_loglik
            target.native_likelihood_status='ready'
    return result
