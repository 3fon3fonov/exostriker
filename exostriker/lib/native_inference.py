"""Bridge the established inference vectors/results to the joint native model.

The sampler API and parameter layout remain unchanged. Context contains only
picklable model inputs; each trial creates fresh parameter dictionaries.
"""
from types import SimpleNamespace
import copy
import numpy as np
from native_likelihood import from_session_likelihood
from native_transit_model import from_session, _mass_axes


def configure(session, opt):
    from activity_gp import configure as configure_activity_gp
    configure_activity_gp(session, opt)
    backend = getattr(session, 'transit_likelihood_backend', 'auto')
    if backend not in ('auto', 'native', 'legacy'):
        raise ValueError('transit_likelihood_backend must be auto, native or legacy.')
    reason = None
    if not (session.type_fit.get('Transit') or session.type_fit.get('TTV')): reason = 'transit/timing likelihood is not selected'
    elif session.rtg[1] or session.rtg[3]: reason = 'Gaussian processes'
    elif opt['AST']: reason = 'astrometry'
    elif opt['AMD_stab'] or opt['Nbody_stab']: reason = 'stability penalties'
    elif opt['gr_flag']: reason = 'general relativity'
    elif session.get_TTVs[0]: reason = 'independent per-event transit times'
    elif session.mod_dynamical and bool(np.asarray(session.mixed_fit[0]).item()): reason = 'mixed dynamical subsets'
    elif any(session.omega_dot[i] != 0 or session.omega_dot_use[i] for i in opt['planet_slots']): reason = 'apsidal precession'
    elif any(session.K[i] < 0 for i in opt['planet_slots']): reason = 'negative RV amplitudes'
    if reason is None:
        from native_joint_fit import unsupported_reason
        reason = unsupported_reason(session)
    if backend == 'native' and reason:
        raise ValueError('Native inference does not support ' + reason + '.')
    enabled = backend != 'legacy' and reason is None
    session.likelihood_backend_used = 'native' if enabled else 'legacy'
    session.likelihood_backend_reason = '' if enabled else (reason or 'explicit legacy selection')
    opt['native_context'] = None
    opt['native_joint_problem'] = None
    if not enabled: return
    names = ('type_fit', 'rtg', 'mod_dynamical', 'epoch', 'stellar_radius', 'rho',
             'gr_flag', 'use_planet', 'ld_m', 'ld_u_lin', 'ld_u_quad', 'ld_u_nonlin',
             'tra_off', 'tra_jitt', 'tra_lintr', 'tra_quadtr', 'tra_dil', 'rvoff', 'jitt',
             'rv_lintr', 'rv_quadtr', 'K', 'P', 'e', 'w', 'M0', 'i', 'Node', 't0',
             'pl_a', 'pl_rad', 'e_sinw', 'e_cosw', 'lamb', 'hkl', 'omega_dot',
             'copl_incl', 'jit_flag', 'get_TTVs')
    values = {name: copy.deepcopy(getattr(session, name)) for name in names}
    values.update(rv_data_sets=session.rv_data_sets, tra_data_sets=session.tra_data_sets,
                  native_jacobi_convention=getattr(session, 'native_jacobi_convention', 'auto'),
                  params=SimpleNamespace(stellar_mass=session.params.stellar_mass),
                  fit_results=None, rvmod=None)
    context=SimpleNamespace(**values)
    opt['native_context'] = context
    from native_joint_fit import prepare
    opt['native_joint_problem'] = prepare(copy.deepcopy(session),free_indices=[],bounds=[])
    opt['native_joint_phase']='t0' if opt['native_joint_problem'].config[6] else 'M0'
    context.epoch=float(opt['native_joint_problem'].settings[0])


def trial_session(par, npl, ndset, tr_model, opt):
    s = copy.copy(opt['native_context'])
    slots = opt['planet_slots']; datasets = opt['dataset_slots']
    # Never mutate the context shared by sampler calls or the GUI session.
    for name in ('K','P','e','w','M0','i','Node','e_sinw','e_cosw','lamb','t0','pl_a','pl_rad',
                 'rvoff','jitt','tra_off','tra_jitt','tra_lintr','tra_quadtr','tra_dil',
                 'ld_u_lin','ld_u_quad','ld_u_nonlin'):
        setattr(s, name, copy.deepcopy(getattr(s, name)))
    rows = np.asarray(par[2*ndset:2*ndset+7*npl]).reshape(npl,7)
    names = ('K','P','e_sinw','e_cosw','lamb','i','Node') if s.hkl else ('K','P','e','w','M0','i','Node')
    for j, slot in enumerate(slots):
        for k, name in enumerate(names): getattr(s,name)[slot] = rows[j,k]
    for j, slot in enumerate(datasets):
        s.rvoff[slot] = par[j]; s.jitt[slot] = abs(par[ndset+j])
    start = 2*ndset+7*npl
    s.rv_lintr, s.rv_quadtr = par[start:start+2]
    start += 2
    for j, slot in enumerate(slots):
        s.t0[slot], s.pl_a[slot], s.pl_rad[slot] = par[start+3*j:start+3*j+3]
    start += 3*npl
    loaded = [d for d in range(60) if len(s.tra_data_sets[d])]; count = len(loaded)
    for j,d in enumerate(loaded):
        s.tra_off[d] = par[start+j]; s.tra_jitt[d] = abs(par[start+count+j])
        s.tra_lintr[d] = par[start+2*count+j]; s.tra_quadtr[d] = par[start+3*count+j]
        law = tr_model[0][d]
        sizes = {'uniform':0,'linear':1,'quadratic':2,'nonlinear':4}
        if law not in sizes: raise ValueError('Unsupported limb-darkening law: '+str(law))
        if sizes[law]:
            index = start+4*count+npl+int(tr_model[3][d])
            getattr(s, {'linear':'ld_u_lin','quadratic':'ld_u_quad','nonlinear':'ld_u_nonlin'}[law])[d] = list(par[index:index+sizes[law]])
    for d,index in opt.get('dilution_indices',{}).items(): s.tra_dil[d] = par[index]
    # The stellar mass is the final entry of the established parameter vector.
    s.params = SimpleNamespace(stellar_mass=float(par[-1]))
    if s.mod_dynamical and s.copl_incl and len(slots):
        for slot in slots: s.i[slot] = s.i[slots[0]]
    return s


def evaluate_loglik(par, npl, ndset, tr_model, opt):
    """Scalar callback used by optimizers and both sampling drivers."""
    if opt.get('native_joint_problem') is not None:
        return opt['native_joint_problem'].loglik(par)
    return evaluate_trial(par,npl,ndset,tr_model,opt).loglik


def evaluate_trial(par, npl, ndset, tr_model, opt):
    s = trial_session(par,npl,ndset,tr_model,opt)
    return from_session_likelihood(s, phase='M0' if s.mod_dynamical else 't0', store=False,
                                   step=opt['dt']/86400., tolerance=opt['eps'])


def transit_results(s, result, dense, dense_times, dense_labels):
    """Keep all five legacy transit-result sections and all 60 dataset slots."""
    columns = [[[] for _ in range(60)] for _ in range(7)]
    physical, baseline, shifted, velocities = {}, {}, {}, {}
    for d in range(60):
        mask = (result.kind==1)&(result.dataset==d)
        if not np.any(mask): continue
        t=result.times[mask]; off=s.tra_off[d]
        model=result.data_model[mask]+off
        x=result.covariate[mask]
        trend=s.tra_lintr[d]*x+s.tra_quadtr[d]*x*x
        vals=(t,result.observed[mask]+off,np.sqrt(result.variance[mask]),model,
              result.residual[mask],result.residual[mask],trend if s.tra_data_sets[d][10] else model)
        for col,value in zip(columns,vals): col[d]=value.copy()
        physical[d]=result.flux[mask].copy();baseline[d]=model.copy()
        shifted[d]=result.data_model[mask].copy();velocities[d]=result.rv[mask].copy()
    joined=[np.concatenate([v for v in col if len(v)]) if any(len(v) for v in col) else np.empty(0) for col in columns]
    joined[6]=joined[3].copy()  # Legacy no-GP all-data convention.
    if dense is None:
        rt=np.empty(0);rf=np.empty(0)
    else:
        rt=dense_times;rf=dense.flux.copy()
        for d in np.unique(dense_labels):
            mask=dense_labels==d; data=s.tra_data_sets[d];t=np.asarray(data[0])
            if data[10]:
                order=np.argsort(t,kind='stable');x=np.interp(rt[mask],t[order],np.asarray(data[3])[order])
            else:x=rt[mask]-t[0]
            dil=s.tra_dil[d]
            rf[mask]=(rf[mask]+s.tra_lintr[d]*x+s.tra_quadtr[d]*x*x)*dil+1-dil
        order=np.argsort(rt);rt=rt[order];rf=rf[order]
    residual=joined[4];sigma=joined[2];chi=float(result.chi2[1]);n=len(residual)
    stats=[chi,chi/n if n else np.nan,float(np.sqrt(np.mean(residual**2))) if n else np.nan,
           float(np.sqrt(np.average(residual**2,weights=1/sigma**2))) if n else np.nan,n]
    output=np.empty(5,dtype=object)
    sep=np.empty((7,60),dtype=object)
    for i,col in enumerate(columns):
        for d,value in enumerate(col):sep[i,d]=value
    all_data=np.empty(7,dtype=object)
    for i,value in enumerate(joined):all_data[i]=value
    rich=np.empty(2,dtype=object);rich[:]=[rt,rf]
    output[:]=[result.transit_loglik,sep,all_data,rich,stats]
    return output,(physical,baseline,shifted,velocities)


def publish_results(obj, par, opt, tr_model, errors, free_indices):
    """Evaluate the final trial and adapt it to the established GUI containers."""
    from RV_mod.rvmod import Rvfit, RvModelWrap, PlanetParamsWrap, PlanetParamsErrorsWrap, ParamsWrap
    from native_model_results import model_time_grid
    npl=len(opt['planet_slots']);ndset=len(opt['dataset_slots'])
    s=trial_session(par,npl,ndset,tr_model,opt)
    phase=opt.get('native_joint_phase','M0' if s.mod_dynamical or not s.type_fit.get('Transit',False) else 't0')
    numerical=dict(step=opt['dt']/86400.,tolerance=opt['eps'])
    numerical.update(opt.get('native_joint_model_options',{}))
    joint=opt.get('native_joint_result')
    shared_propagation=joint is not None and joint.diagnostics['physical_evaluations'].get('shared',0)>0
    if joint is not None:
        from native_joint_fit import likelihood_from_joint
        result=likelihood_from_joint(s,joint)
    else:
        result=from_session_likelihood(s,phase=phase,store=False,**numerical)
    # Dense output is requested only after inference, never per sampler trial.
    grid=model_time_grid(obj,obj.epoch)
    transit_grid=model_time_grid(obj,obj.epoch,component='transit')
    loaded=[d for d in range(60) if s.type_fit.get('Transit',False) and len(s.tra_data_sets[d])]
    rich_times=[];rich_labels=[]
    for d in loaded:
        t=np.asarray(s.tra_data_sets[d][0]);factor=int(obj.tra_model_fact)
        if factor<1:raise ValueError('tra_model_fact must be positive.')
        rich=np.linspace(t[:-1]-.3,t[1:]+.3,factor+1,axis=1)[:,:-1].ravel()
        rich=np.append(rich,t[-1]);rich_times.append(rich);rich_labels.append(np.full(len(rich),d))
    rich_t=np.concatenate(rich_times) if rich_times else np.empty(0);rich_d=np.concatenate(rich_labels) if rich_labels else np.empty(0,dtype=int)
    # Evaluate separate RV/LC display grids and rich LC samples together.
    shared=bool(getattr(obj,'calculate_transit_model',True))
    grid_slots=(loaded if shared else loaded[:1]) or [0]
    same_grid=np.array_equal(grid,transit_grid)
    shared_t=np.tile(transit_grid,len(grid_slots));shared_d=np.repeat(grid_slots,len(transit_grid))
    transit_start=0
    if not same_grid:
        transit_start=len(grid)
        shared_t=np.r_[grid,shared_t];shared_d=np.r_[np.full(len(grid),grid_slots[0]),shared_d]
    times=np.r_[shared_t,rich_t];labels=np.r_[shared_d,rich_d]
    if npl:
        dense=from_session(s,times,dataset=labels,phase=phase,_shared_nbody=shared_propagation,**numerical)
    else:
        from native_transit_model import TransitRVModel
        dense=TransitRVModel(times,labels,np.zeros(len(times)),np.ones(len(times)),
            np.ones(len(times)),np.zeros(len(times)),np.empty(0),np.empty(0),
            np.empty((len(times),0)),np.empty(0),None)
    from dataclasses import replace
    cut=len(shared_t)
    rich_model=replace(dense,times=dense.times[cut:],dataset=dense.dataset[cut:],
                       flux=dense.flux[cut:])
    tr,curves=transit_results(s,result,rich_model,rich_t,rich_d)
    obj.transit_results=tr
    r=Rvfit(); slots=opt['planet_slots'];ds=opt['dataset_slots']
    r.npl=npl;r.ndset=len(ds);r.npl_pos=np.asarray(obj.use_planet).copy()
    r.epoch=s.epoch;r.stellar_mass=s.params.stellar_mass;r.coplar_inc=s.copl_incl
    r.model_type='dyn' if s.mod_dynamical else 'kep'
    mask=result.kind==0;index=np.flatnonzero(mask)
    index=index[np.argsort(result.times[index])]
    r.jd=result.times[index].copy();r.idset=result.dataset[index].copy()
    off=np.array([s.rvoff[d] for d in r.idset]);elapsed=r.jd-s.epoch
    r.rvs=result.observed[index]-off;r.rv_err=result.uncertainty[index].copy()
    r.o_c=result.residual[index].copy();r.model_rvs=result.rv[index].copy()
    r.rv_model=RvModelWrap(r.jd,r.rvs,r.rv_err,r.o_c)
    r.ndata=r.Ndata=len(index) if len(index) else tr[4][4]
    # Rvfit statistics describe the RV component in a joint fit; LC-only
    # nuisance parameters must not consume its degrees of freedom.
    rv_end=2*ndset+7*npl+2
    r.mfit=sum(i<rv_end or i==len(par)-1 for i in free_indices) if len(index) else len(free_indices)
    r.dof=r.Ndata-r.mfit
    r.loglik=result.rv_loglik if len(index) else result.transit_loglik
    r.chi2=float(result.chi2[0]) if len(index) else tr[4][0]
    residual=r.o_c if len(index) else np.asarray(tr[2][4],float)
    variance=result.variance[index] if len(index) else np.asarray(tr[2][2],float)**2
    r.rms=float(np.sqrt(np.mean(residual**2))) if len(residual) else np.nan
    weights_variance=r.rv_err**2 if len(index) and np.all(r.rv_err>0) else variance
    r.wrms=float(np.sqrt(np.average(residual**2,weights=1/weights_variance))) if len(residual) else np.nan
    r.reduced_chi2=r.chi2/r.dof if r.dof>0 else np.nan
    r.rms2=r.rms;r.reduced_chi22=r.reduced_chi2
    r.to_sort=np.argsort([s.P[i] for i in slots],kind='stable')
    # Preserve the established h/k/lambda representation in RV parameter fields.
    names=('K','P','e_sinw','e_cosw','lamb','i','Node') if s.hkl else ('K','P','e','w','M0','i','Node')
    fields=('K','P','e','w','M0','incl','cap0m')
    for k,(name,field) in enumerate(zip(names,fields)):
        values=np.zeros((9,2))
        for j,i in enumerate(slots):
            values[i]=[getattr(s,name)[i],max(errors[2*ndset+7*j+k])]
        setattr(r,field,values)
    rows=np.array([[s.K[i],s.P[i],s.e[i],s.w[i],s.M0[i],s.i[i],s.Node[i]] for i in slots]).reshape(-1,7)
    if s.hkl and npl:
        h=np.array([s.e_sinw[i] for i in slots]);k=np.array([s.e_cosw[i] for i in slots])
        rows[:,2]=np.hypot(h,k);rows[:,3]=np.degrees(np.arctan2(h,k))
    rows[:,4]=result.mean_anomaly
    masses,axes=_mass_axes(rows,s.params.stellar_mass)
    if s.mod_dynamical:
        from native_transit_model import jacobi_axis_factor, session_jacobi_convention
        axes *= jacobi_axis_factor(rows,masses,s.params.stellar_mass,session_jacobi_convention(s))
    r.mass=np.zeros(9);r.a=np.zeros(9)
    r.mass[list(slots)]=masses*1047.5654817267318;r.a[list(slots)]=axes
    for field in ('t0','RR','aR','a_new','m','tw','h','k','lamb'):
        setattr(r,field,np.zeros((9,2)))
    for j,i in enumerate(slots):
        r.t0[i,0]=s.t0[i] if result.transit_epoch is None else result.transit_epoch[j]
        r.RR[i,0]=s.pl_rad[i];r.aR[i,0]=result.scaled_axis[j]
        r.h[i,0]=rows[j,2]*np.sin(np.radians(rows[j,3]));r.k[i,0]=rows[j,2]*np.cos(np.radians(rows[j,3]))
        r.lamb[i,0]=(rows[j,3]+rows[j,4])%360.
        r.planet_params[7*i:7*i+7]=[getattr(s,name)[i] for name in names]
        r.M0[i,0]=(rows[j,3]+rows[j,4])%360. if s.hkl else rows[j,4]
        r.planet_params[7*i+4]=r.M0[i,0]
        r.planet_params_errors[7*i:7*i+7]=errors[2*ndset+7*j:2*ndset+7*j+7]
        obj.M0[i]=rows[j,4];obj.lamb[i]=(rows[j,3]+rows[j,4])%360.
        if result.transit_epoch is not None: obj.t0[i]=result.transit_epoch[j]
        obj.i[i]=s.i[i]
        obj.masses[i]=r.mass[i];obj.semimajor[i]=r.a[i]
    capacity=max(20,max(ds,default=-1)+1)
    r.offsets=np.zeros(capacity);r.jitters=np.zeros(capacity)
    r.offset_errors=np.zeros((capacity,2));r.jitter_errors=np.zeros((capacity,2))
    for j,d in enumerate(ds):
        r.offsets[d]=s.rvoff[d];r.jitters[d]=abs(s.jitt[ds[0] if s.jit_flag else d])
        r.offset_errors[d]=errors[j];r.jitter_errors[d]=errors[ndset+j]
    extent=max(ds,default=-1)+1
    r.V0=np.column_stack((r.offsets[:extent],np.max(r.offset_errors[:extent],axis=1)))
    r.jitt=np.column_stack((r.jitters[:extent],np.max(r.jitter_errors[:extent],axis=1)))
    start=2*ndset+7*npl
    r.rv_lintr=s.rv_lintr;r.rv_quadtr=s.rv_quadtr
    r.rv_lintr_err=np.asarray(errors[start]);r.rv_quadtr_err=np.asarray(errors[start+1])
    r.lin_trend_out=np.array([r.rv_lintr,max(r.rv_lintr_err)])
    r.quad_trend_out=np.array([r.rv_quadtr,max(r.rv_quadtr_err)])
    r.params=PlanetParamsWrap(r.offsets,r.jitters,r.planet_params,r.rv_lintr,r.stellar_mass)
    r.params_errors=PlanetParamsErrorsWrap(r.offset_errors,r.jitter_errors,r.planet_params_errors,r.rv_lintr_err,np.asarray(errors[-1]))
    r.stat=ParamsWrap(r.params,r.params_errors,dof=r.dof)
    r.stat_array_saved=True;r.ind_orig=r.idset+1;r.ind_used=np.asarray(ds)+1
    r.model_jd=grid.copy();delta=grid-s.epoch
    r.model=dense.rv[:len(grid)]+s.rv_lintr*delta+s.rv_quadtr*delta**2
    parts=dense.rv_components[:len(grid),r.to_sort]
    r.fit=np.column_stack((grid,r.model,parts))
    obs_parts=result.rv_components[index][:,r.to_sort]
    compact_ids=np.array([ds.index(int(d))+1 for d in r.idset],dtype=int)
    data=np.column_stack((r.jd,r.rvs,r.rv_err,compact_ids,r.o_c,r.model_rvs,obs_parts))
    r.model_data=data.T
    r.res_dict={name:getattr(r,name) for name in ('jd','rvs','rv_err','idset','o_c','model_rvs','loglik',
                'reduced_chi2','chi2','rms','K','P','e','w','M0','incl','cap0m','t0','aR','RR',
                'a_new','m','tw','h','k','lamb','V0','jitt','lin_trend_out','quad_trend_out',
                'Ndata','mfit','rms2','reduced_chi22','epoch','fit')}
    active=np.asarray(slots,dtype=int)[r.to_sort]
    for field in fields+('t0','aR','RR','a_new','m','tw','h','k','lamb'):
        r.res_dict[field]=getattr(r,field)[active].copy()
    r.res_dict.update(model_data=data,idset=compact_ids-1,V0=r.V0[list(ds)].copy(),
                      jitt=r.jitt[list(ds)].copy(),w_dot=np.zeros((npl,2)),
                      jup_mass=(masses*1047.5654817267318)[r.to_sort],jacobi_major_ax=axes[r.to_sort])
    r.res=r.res_dict.copy()
    r.native_likelihood=result;r.native_loglik=result.loglik
    r.native_rv_loglik=result.rv_loglik;r.native_transit_loglik=result.transit_loglik
    r.native_likelihood_status='ready'
    physical,baseline,shifted,velocities=curves
    r.transit_data_jd={d:np.asarray(s.tra_data_sets[d][0]).copy() for d in loaded}
    r.transit_data_models=physical;r.transit_data_models_with_baseline=baseline
    r.transit_data_models_with_offset=shifted;r.transit_data_rv=velocities
    r.transit_data_phase=phase;r.transit_data_mean_anomaly=result.mean_anomaly.copy()
    r.transit_data_t0=result.transit_epoch
    if loaded:
        r.transit_data_model_dataset=loaded[0];r.transit_data_model=physical[loaded[0]]
    if shared and loaded:
        r.transit_model_jd=r.model_jd if same_grid else transit_grid.copy()
        r.transit_model_dataset=loaded[0]
        r.transit_models={d:dense.flux[transit_start+j*len(transit_grid):transit_start+(j+1)*len(transit_grid)].copy() for j,d in enumerate(loaded)}
        r.transit_models_with_offset={d:r.transit_models[d]-s.tra_off[d] for d in loaded}
        r.transit_model=r.transit_models[loaded[0]]
        r.transit_model_with_offset=r.transit_models_with_offset[loaded[0]]
        r.transit_model_status='ready'
    obj.fit_results=r;obj.rvmod=copy.deepcopy(r)
    obj.params=r.params;obj.param_errors=r.params_errors
    obj.loglik=result.loglik;obj.chi2=r.chi2;obj.reduced_chi2=r.reduced_chi2;obj.rms=r.rms
    obj.flag=1;obj.stat_saved=True;obj.fit_performed=True;obj.model_saved=True
    if getattr(obj,'calculate_transit_events',False):
        event_grid=model_time_grid(obj,obj.epoch,component='all')
        obj.update_transit_events(start=event_grid[0],end=event_grid[-1],phase=phase,**numerical)
    return result
