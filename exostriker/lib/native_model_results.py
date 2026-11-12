"""Publish RV and transit plotting curves without changing fit scores."""
import numpy as np


def clear_transit_models(result, reason='not calculated'):
    """Remove cached transit predictions and record why they are absent."""

    result.transit_model = None
    result.transit_model_jd = None
    result.transit_model_with_offset = None
    result.transit_model_dataset = None
    result.transit_models = {}
    result.transit_models_with_offset = {}
    result.transit_model_status = reason
    result.transit_events = {}
    result.transit_events_status = 'not calculated'
    result.ttv_model = np.empty((0, 5))
    clear_transit_data_models(result)


def shared_model_reason(session):
    """Retain the established RV plotter for unsupported physical models."""
    if not getattr(session, 'calculate_transit_model', True):
        return 'disabled'
    if session.gr_flag:
        return 'GR is not supported by the shared transit model'
    slots = np.flatnonzero(session.use_planet)
    if any((session.omega_dot[i] != 0 for i in slots)):
        return 'apsidal precession is not supported by the shared transit model'
    if any((session.K[i] < 0 for i in slots)):
        return 'negative K is not supported by the shared transit model'
    if session.mod_dynamical:
        mixed = getattr(session, 'mixed_fit', {0: False})
        if bool(np.asarray(mixed[0]).item()):
            return 'mixed dynamical subsets are not supported by the shared transit model'
    return None


def model_time_grid(session, epoch, *, component='rv'):
    """Use observation limits and the RV plot controls, independent of epoch.

    RV curves prefer RV dates; transit curves prefer photometry dates. The
    all-data range remains available for optional event searches.
    """
    if component not in ('rv', 'transit', 'all'):
        raise ValueError('Plotting component must be rv, transit or all.')
    ranges = {}
    for (name, column) in (('rv_data_sets', 0), ('tra_data_sets', 0), ('ttv_data_sets', 1)):
        bounds = []
        for data in getattr(session, name, {}).values():
            if not len(data):
                continue
            if name == 'ttv_data_sets' and (not bool(data[4])):
                continue
            times = np.asarray(data[column], dtype=float)
            times = times[np.isfinite(times)]
            if times.size:
                bounds.extend((np.min(times), np.max(times)))
        ranges[name] = bounds
    rv, transit, timing = (ranges[name] for name in ('rv_data_sets', 'tra_data_sets', 'ttv_data_sets'))
    if component == 'rv':
        bounds = rv or transit + timing
    elif component == 'transit':
        bounds = transit or timing or rv
    else:
        bounds = [float(epoch)] + rv + transit + timing
    if not bounds:
        bounds = [float(epoch)]
    npoints = session.model_npoints
    if not np.isfinite(npoints) or npoints != int(npoints) or npoints < 2:
        raise ValueError('model_npoints must be an integer >= 2.')
    margins = np.array([session.model_min, session.model_max], dtype=float)
    if not np.all(np.isfinite(margins)) or np.any(margins < 0):
        raise ValueError('Plotting margins model_min/model_max must be finite and nonnegative.')
    lo, hi = (min(bounds) - margins[0], max(bounds) + margins[1])
    if lo == hi:
        periods = [session.P[i] for i in np.flatnonzero(session.use_planet)]
        hi = lo + (min(periods) if periods else 1.)
    return np.linspace(lo, hi, int(npoints))


def update_model_curves(session, **options):
    """One physical calculation for separate RV/transit plotting ranges.

    No likelihoods, observation arrays, or parameter errors are changed.
    The primary curve uses the first loaded transit dataset, or slot 0 if none.
    """
    targets = [session.fit_results]
    if getattr(session, 'rvmod', None) is not None and session.rvmod is not session.fit_results:
        targets.append(session.rvmod)
    reason = shared_model_reason(session)
    for target in targets:
        clear_transit_models(target, reason or 'not calculated')
    if reason:
        if getattr(session, 'calculate_transit_events', False):
            from native_transit_events import update_transit_events
            update_transit_events(session, **options)
        return None
    from native_transit_model import from_session
    from model_epoch import resolve_epoch
    default_epoch = resolve_epoch(session)
    epoch = float(options.get('epoch', default_epoch))
    times = model_time_grid(session, epoch)
    transit_times = model_time_grid(session, epoch, component='transit')
    slots = [d for d in range(60) if len(session.tra_data_sets[d])] or [0]
    # Repeated timestamps with different passbands share the same propagated
    # N-body state. Only the occultation integral changes between datasets.
    same_grid = np.array_equal(times, transit_times)
    query_times = np.tile(transit_times, len(slots))
    query_slots = np.repeat(slots, len(transit_times))
    transit_start = 0
    if not same_grid:
        transit_start = len(times)
        query_times = np.r_[times, query_times]
        query_slots = np.r_[np.full(len(times), slots[0]), query_slots]
    result = from_session(session, query_times, dataset=query_slots, **options)
    flux = result.flux[transit_start:].reshape(len(slots), len(transit_times))
    shifted_flux = result.model_flux[transit_start:].reshape(len(slots), len(transit_times))
    elapsed = times - epoch
    rv = result.rv[:len(times)] + session.rv_lintr * elapsed + session.rv_quadtr * elapsed ** 2
    # Existing phased RV plotting uses component columns in increasing-period
    # order, even though session/native adapter outputs follow enabled-slot order.
    active = list(np.flatnonzero(session.use_planet))
    order = np.argsort([session.P[i] for i in active], kind='stable')
    components = result.rv_components[:len(times), order]
    fit_array = np.column_stack((times, rv, components))
    for target in targets:
        target.model_jd = times.copy()
        target.model = rv.copy()
        target.fit = fit_array.copy()
        from magnetic_cycle import active as cycle_active, configure as configure_cycle, mean as cycle_mean
        if cycle_active(session) and getattr(target, '_magnetic_cycle_applied', False):
            session.prepare_for_mcmc()
            cycle_options = {'activity_gp': None}
            configure_cycle(session, cycle_options)
            cycle = cycle_mean(cycle_options['magnetic_cycle'], session.parameters, times)
            target.model += cycle
            target.fit[:, 1] = target.model
            session.magnetic_cycle_model_curve = cycle
        for field in ('res', 'res_dict'):
            mapping = getattr(target, field, None)
            if isinstance(mapping, dict):
                mapping['fit'] = target.fit.copy()
        target.transit_model_jd = target.model_jd if same_grid else transit_times.copy()
        target.transit_model_dataset = slots[0]
        target.transit_models = {slot: flux[i].copy() for (i, slot) in enumerate(slots)}
        target.transit_models_with_offset = {slot: shifted_flux[i].copy() for (i, slot) in enumerate(slots)}
        target.transit_model = target.transit_models[slots[0]]
        target.transit_model_with_offset = target.transit_models_with_offset[slots[0]]
        target.transit_model_status = 'ready'
    session.model_saved = True
    from native_transit_events import update_transit_events
    event_options = {k: v for (k, v) in options.items() if k in ('epoch', 'model', 'phase', 'step', 'tolerance', 'max_steps', 'rho_star', 'stellar_mass', 'stellar_radius', 'radius_mode', 'jacobi_convention')}
    event_times = model_time_grid(session, epoch, component='all')
    update_transit_events(session, start=event_times[0], end=event_times[-1], **event_options)
    return session.fit_results.transit_model


def clear_transit_data_models(result):
    """Remove predictions cached at the observed transit times."""

    result.transit_data_jd = {}
    result.transit_data_models = {}
    result.transit_data_models_with_baseline = {}
    result.transit_data_models_with_offset = {}
    result.transit_data_rv = {}
    result.transit_data_model = None
    result.transit_data_model_dataset = None
    result.transit_data_phase = None
    result.transit_data_mean_anomaly = None
    result.transit_data_t0 = None


def update_transit_data_models(session, *, dataset=None, phase='M0', **options):
    """Evaluate only actual transit timestamps, leaving dense models untouched.

    The `with_baseline` curves reproduce the old GUI's baseline/dilution frame;
    physical `transit_data_models` remain normalized, uncontaminated flux.
    """
    targets = [session.fit_results]
    if getattr(session, 'rvmod', None) is not None and session.rvmod is not session.fit_results:
        targets.append(session.rvmod)
    for target in targets:
        clear_transit_data_models(target)
    if not getattr(session, 'calculate_transit_model', True):
        return None
    if dataset is None:
        slots = [i for i in range(60) if len(session.tra_data_sets[i])]
    else:
        if not np.isscalar(dataset) or not np.isfinite(dataset) or dataset != int(dataset) or (not 0 <= dataset < 60):
            raise ValueError('dataset must be a single slot in 0..59, or None for all loaded slots.')
        slots = [int(dataset)]
    if not slots or any((not len(session.tra_data_sets[i]) for i in slots)):
        raise ValueError('Load a transit dataset before generating data-time models.')
    times = [np.asarray(session.tra_data_sets[i][0], dtype=float) for i in slots]
    if any((t.ndim != 1 or not t.size for t in times)):
        raise ValueError('Transit datasets must contain nonempty one-dimensional time arrays.')
    from native_transit_model import from_session
    labels = np.concatenate([np.full(len(t), i) for (i, t) in zip(slots, times)])
    result = from_session(session, np.concatenate(times), dataset=labels, phase=phase, **options)
    lengths = np.cumsum([0] + [len(t) for t in times])
    physical, baselined, shifted, velocities = ({}, {}, {}, {})
    for (k, i) in enumerate(slots):
        sl = slice(lengths[k], lengths[k + 1])
        physical[i] = result.flux[sl].copy()
        data = session.tra_data_sets[i]
        # Match the legacy data-time baseline exactly: raw airmass or elapsed
        # time from the first sample (not the integration epoch).
        covariate = np.asarray(
            data[3],
            dtype=float,
        ) if len(data) > 10 and data[10] else times[k] - times[k][0]
        linear = getattr(session, 'tra_lintr', {}).get(i, 0.)
        quadratic = getattr(session, 'tra_quadtr', {}).get(i, 0.)
        dilution = getattr(session, 'tra_dil', {}).get(i, data[8] if len(data) > 8 else 1.)
        baselined[i] = (physical[i] + linear * covariate + quadratic * covariate ** 2) * dilution + 1. - dilution
        shifted[i] = baselined[i] - session.tra_off[i]
        velocities[i] = result.rv[sl].copy()
    for target in targets:
        target.transit_data_jd = {i: t.copy() for (i, t) in zip(slots, times)}
        target.transit_data_models = physical.copy()
        target.transit_data_models_with_baseline = baselined.copy()
        target.transit_data_models_with_offset = shifted.copy()
        target.transit_data_rv = velocities.copy()
        target.transit_data_model_dataset = slots[0]
        target.transit_data_model = target.transit_data_models[slots[0]]
        target.transit_data_phase = phase
        target.transit_data_mean_anomaly = result.mean_anomaly.copy()
        target.transit_data_t0 = None if result.transit_epoch is None else result.transit_epoch.copy()
    return session.fit_results.transit_data_models
