"""Deterministic, shared RV/activity magnetic cycle, independent of celerite."""
import copy
import numpy as np

GLOBAL = ('period', 'phase', 'rv_amplitude')


def parameter(value, bounds, fit=True):
    """Create a value, fit flag, bounds and priors for one cycle parameter."""

    return dict(
        value=float(value),
        fit=fit,
        bounds=list(bounds),
        normal=[float(value), max(abs(value) * .1, 1.), False],
        jeffreys=[max(bounds[0], 1e-12), bounds[1], False],
        error=[0., 0.],
    )


def ensure(s):
    """Add missing cycle settings without replacing saved session choices."""

    from activity_gp import ensure as ensure_activity
    ensure_activity(s)

    # Only read RV times when a reference epoch still needs to be chosen.
    needs_times = not hasattr(s, 'magnetic_cycle') or s.magnetic_cycle['reference_epoch'] == 0.
    times = []

    if needs_times:
        for data in s.rv_data_sets.values():
            if len(data):
                times.append(np.asarray(data[0], float))

    time = np.concatenate(times) if times else np.array([])

    if not hasattr(s, 'magnetic_cycle'):
        span = float(np.ptp(time)) if len(time) else 0.
        scatter = []

        for data in s.rv_data_sets.values():
            if len(data):
                scatter.append(np.asarray(data[1], float) - np.median(data[1]))

        if scatter:
            amplitude = float(np.std(np.concatenate(scatter)))
        else:
            amplitude = 1.

        amplitude = max(amplitude, 1.)

        s.magnetic_cycle = dict(
            enabled=False,
            reference_epoch=0.,
            period=parameter(max(3650., 2 * span), [1., 1e7]),
            phase=parameter(0., [-360., 360.]),
            rv_amplitude=parameter(amplitude, [0., 1e6]),
            indicators={},
        )

    cycle = s.magnetic_cycle

    if cycle['reference_epoch'] == 0. and len(time):
        cycle['reference_epoch'] = float(np.min(time))

    # New files receive their own response settings. Existing files keep theirs.
    for slot, data in s.act_data_sets.items():
        if len(data) and slot not in cycle['indicators']:
            amplitude = max(float(np.std(data[1])) * np.sqrt(2), 1e-12)
            cycle['indicators'][slot] = dict(
                enabled=False,
                amplitude=parameter(amplitude, [-1e6, 1e6]),
            )

    for indicator in cycle['indicators'].values():
        indicator.setdefault('phase_offset', parameter(0., [-180., 180.], fit=False))

    return cycle


def active(s):
    """Check whether the RV magnetic-cycle component is enabled."""

    return bool(getattr(s, 'magnetic_cycle', {}).get('enabled', False) and s.type_fit.get('RV'))


def selected(s):
    """Return the loaded activity slots included in the magnetic cycle."""

    cycle = ensure(s)

    if not active(s):
        return []

    slots = []
    for slot, indicator in cycle['indicators'].items():
        if indicator['enabled'] and len(s.act_data_sets.get(slot, [])):
            slots.append(slot)

    return slots


def setting(s, slot, name):
    """Find a cycle parameter or its shared activity offset/jitter."""

    if slot is None:
        return s.magnetic_cycle[name]
    if name in ('amplitude', 'phase_offset'):
        return s.magnetic_cycle['indicators'][slot][name]
    return s.activity_gp[slot][name]


def pack(s, par, flags, labels, bounds, normal, jeffreys):
    """Append cycle parameters and priors to the inference vector."""

    ensure(s)
    s.magnetic_cycle_indices = {}
    if not active(s):
        return
    if not any((len(parameter_settings) for parameter_settings in s.rv_data_sets.values())):
        raise ValueError('Magnetic cycle requires at least one RV dataset.')
    if not np.isfinite(s.magnetic_cycle['reference_epoch']):
        raise ValueError('Magnetic cycle reference epoch must be finite.')
    entries = [(None, name) for name in GLOBAL]
    for slot in selected(s):
        data = np.asarray([s.act_data_sets[slot][k] for k in range(3)], float)
        if not np.all(np.isfinite(data)) or np.any(data[2] <= 0):
            raise ValueError(f'Magnetic cycle activity {slot + 1} needs finite data and positive errors.')
        entries.extend(((slot, name) for name in ('amplitude', 'phase_offset')))
        for name in ('offset', 'jitter'):
            if (slot, name) not in s.activity_gp_indices:
                entries.append((slot, name))
    units = {
        'period': 'd',
        'phase': 'deg',
        'phase_offset': 'deg',
        'rv_amplitude': 'm/s',
        'amplitude': 'U',
        'offset': 'U',
        'jitter': 'U',
    }
    for (slot, name) in entries:
        parameter_settings = setting(s, slot, name)
        value = float(parameter_settings['value'])
        low, high = map(float, parameter_settings['bounds'])
        normal_prior, jeffreys_prior = (list(parameter_settings['normal']), list(parameter_settings['jeffreys']))
        label = f"Magnetic cycle {('RV' if slot is None else 'activity ' + str(slot + 1))} {name} [{units[name]}]"
        if normal_prior[2] and jeffreys_prior[2]:
            raise ValueError(label + ': choose one prior family.')
        if jeffreys_prior[2]:
            if jeffreys_prior[0] <= 0 or jeffreys_prior[1] <= jeffreys_prior[0]:
                raise ValueError(label + ': Jeffreys support must be strictly positive and ordered.')
            low, high = (max(low, jeffreys_prior[0]), min(high, jeffreys_prior[1]))
        if (not np.all(np.isfinite([value, low, high])) or low >= high or not low <= value <= high
                or (name == 'period' and (value <= 0 or low <= 0))
                or (name in ('jitter', 'rv_amplitude') and low < 0)
                or (normal_prior[2] and (not np.all(np.isfinite(normal_prior[:2])) or normal_prior[1] <= 0))):
            raise ValueError(label + ': invalid starting value, bounds or prior.')
        s.magnetic_cycle_indices[slot, name] = len(par)
        par.append(value)
        flags.append(bool(parameter_settings['fit']))
        labels.append(label)
        bounds.append([low, high])
        normal.append(normal_prior)
        jeffreys.append(jeffreys_prior)


def configure(s, opt):
    """Copy the cycle inputs needed by the likelihood evaluation."""

    opt['magnetic_cycle'] = None

    if not active(s):
        return

    indices = copy.deepcopy(s.magnetic_cycle_indices)

    for slot in selected(s):
        for name in ('offset', 'jitter'):
            if (slot, name) not in indices:
                indices[(slot, name)] = s.activity_gp_indices[(slot, name)]

    data = {}
    for slot in selected(s):
        data[slot] = tuple(
            np.array(s.act_data_sets[slot][column], float, copy=True)
            for column in range(3)
        )

    if opt.get('activity_gp'):
        gp_slots = set(opt['activity_gp']['data'])
    else:
        gp_slots = set()

    context = dict(
        indices=indices,
        reference_epoch=s.magnetic_cycle['reference_epoch'],
        data=data,
        gp_slots=gp_slots,
    )

    opt['magnetic_cycle'] = context

    if opt.get('activity_gp'):
        opt['activity_gp']['magnetic_cycle'] = context


def mean(context, par, time, slot=None):
    """Evaluate the RV sine or one activity sine with its phase offset."""

    t = np.asarray(time, float)

    if not context or (slot is not None and slot not in context['data']):
        return np.zeros_like(t)

    def value_at(key):
        return float(par[context['indices'][key]])

    period = value_at((None, 'period'))

    if period <= 0:
        raise ValueError('Magnetic cycle period must be positive.')

    phase = np.deg2rad(value_at((None, 'phase')))

    if slot is not None:
        phase += np.deg2rad(value_at((slot, 'phase_offset')))

    angle = 2 * np.pi * np.remainder((t - context['reference_epoch']) / period, 1.) + phase

    if slot is None:
        amplitude = value_at((None, 'rv_amplitude'))
    else:
        amplitude = value_at((slot, 'amplitude'))

    return amplitude * np.sin(angle)


def independent_activity(context, par):
    """Likelihood/chi² for cycle-only indicators; GP-overlapping data occur once."""
    loglik = chi2 = 0.
    if context:
        for (slot, (t, y, err)) in context['data'].items():
            if slot in context['gp_slots']:
                continue
            offset = par[context['indices'][slot, 'offset']]
            jitter = par[context['indices'][slot, 'jitter']]
            if jitter < 0:
                return (-np.inf, np.inf)
            variance = err ** 2 + jitter ** 2
            residual = y - offset - mean(context, par, t, slot)
            q = float(np.sum(residual ** 2 / variance))
            chi2 += q
            loglik -= .5 * (q + float(np.log(2 * np.pi * variance).sum()))
    return (loglik, chi2)


def publish_parameters(s, par, errors, opt):
    """Copy fitted cycle values and uncertainties back to the session."""

    if not opt.get('magnetic_cycle'):
        return
    for ((slot, name), index) in s.magnetic_cycle_indices.items():
        v = setting(s, slot, name)
        v['value'] = float(par[index])
        v['error'] = list(map(float, errors[index]))


def apply_rv_model(s):
    """Add the deterministic mean to fresh orbital results, once per result."""
    if not active(s):
        return
    result = s.fit_results
    if getattr(result, '_magnetic_cycle_applied', False):
        return
    s.prepare_for_mcmc()
    opt = {'activity_gp': None}
    configure(s, opt)
    data = mean(opt['magnetic_cycle'], s.parameters, result.jd)
    curve = mean(opt['magnetic_cycle'], s.parameters, result.model_jd)
    result.o_c = np.asarray(result.o_c) - data
    result.rv_model.o_c = result.o_c
    result.model = np.asarray(result.model) + curve
    result.model_rvs = np.asarray(result.model_rvs) + data
    if getattr(result, 'model_data', np.empty((0,))).ndim == 2 and len(result.model_data) > 5:
        result.model_data = result.model_data.copy()
        result.model_data[4] = result.o_c
        result.model_data[5] = result.model_rvs
    if getattr(result, 'fit', np.empty((0,))).ndim == 2 and len(result.fit) == len(curve):
        result.fit = result.fit.copy()
        result.fit[:, 1] = result.model
    for name in ('res', 'res_dict'):
        mapping = getattr(result, name, None)
        if isinstance(mapping, dict):
            mapping.update(o_c=result.o_c.copy(),model_rvs=result.model_rvs.copy(),
                           model_data=result.model_data.T.copy(),fit=result.fit.copy())
    result._magnetic_cycle_applied = True
    s.magnetic_cycle_model_data = data
    s.magnetic_cycle_model_curve = curve


def publish_activity(s):
    """Publish cycle-only curves; GP curves already include the cycle mean."""
    if not active(s):
        return
    s.prepare_for_mcmc()
    opt = {'activity_gp': dict(data={i: None for i in s.activity_gp if s.rtg[1] and s.activity_gp[i]['enabled'] and len(s.act_data_sets.get(i, []))})}
    configure(s, opt)
    context = opt['magnetic_cycle']
    if not s.rtg[1] or not hasattr(s, 'activity_gp_models'):
        s.activity_gp_models = {}
    for (slot, (t, _, _)) in context['data'].items():
        if slot in context['gp_slots']:
            continue
        grid = np.linspace(float(t.min()), float(t.max()), max(200, min(1000, len(t))))
        offset = s.activity_gp[slot]['offset']['value']
        s.activity_gp_models[slot] = dict(
            time=t.copy(),
            mean=offset + mean(context, s.parameters, t, slot),
            variance=np.zeros_like(t),
            curve_time=grid,
            curve_mean=offset + mean(context, s.parameters, grid, slot),
        )


def orbital_refresh(s, **options):
    """Refresh only the orbital backbone while publishing a Python fit."""
    old = getattr(s, '_magnetic_cycle_orbital_refresh', False)
    s._magnetic_cycle_orbital_refresh = True
    try:
        return s.fitting(**options)
    finally:
        s._magnetic_cycle_orbital_refresh = old


def publish_summary(s):
    """Update the joint RV/activity statistics without counting data twice."""

    s.prepare_for_mcmc()
    from activity_gp import configure as configure_gp, selected as gp_selected
    opt = {}
    configure_gp(s, opt)
    configure(s, opt)
    independent, independent_chi2 = independent_activity(opt['magnetic_cycle'], s.parameters)
    result = s.fit_results
    errors = np.array([np.hypot(e, s.jitt[int(i)]) for (e, i) in zip(result.rv_err, result.idset)])
    residual = np.asarray(result.o_c)
    if s.rtg[1] and gp_selected(s):
        loglik, chi2 = (s.activity_gp_loglik, s.activity_gp_chi2)
    elif s.rtg[1]:
        loglik = s.gps.log_likelihood(residual)
        chi2 = float(residual @ np.ravel(s.gps.apply_inverse(residual)))
    else:
        chi2 = float(np.sum((residual / errors) ** 2))
        loglik = -.5 * (chi2 + float(np.log(2 * np.pi * errors ** 2).sum()))
    count = len(residual) + sum((len(s.act_data_sets[i][0]) for i in set(selected(s)) | (set(gp_selected(s)) if s.rtg[1] else set())))
    s.magnetic_cycle_statistics = dict(
        lnL=float(loglik + independent),
        chi2=float(chi2 + independent_chi2),
        N=count,
    )
    after_gp = residual - s.gp_model_data[0] if s.rtg[1] else residual
    result.rms = float(np.sqrt(np.mean(after_gp ** 2)))
    result.wrms = float(np.sqrt(np.average(after_gp ** 2, weights=1 / errors ** 2)))
    result.chi2 = chi2 + independent_chi2
    dof = count - len(s.par_for_mcmc)
    result.reduced_chi2 = result.chi2 / dof if dof > 0 else np.nan
    result.stat.dof = dof
