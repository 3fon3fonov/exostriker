"""Optional transit-event diagnostics from native Cartesian trajectories.

Times/durations are days. Event searches use relative times and Hermite
interpolation of positions AND velocities, independently of the plot grid.
"""
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq, minimize_scalar

from native_transit_model import evaluate, _native, _mass_axes, _LOCK, LD_LAWS


def evaluate_events(planets, times, *, dataset=0, samples_per_orbit=512,
                    max_samples=2000000, event_curve_points=0, **options):
    """Find foreground separation minima; report T14, T23 and isolated depths.

    ``times`` supplies two inclusive search bounds. Padded trajectory samples
    allow complete contacts for events whose midpoint lies at either boundary.
    Increase samples_per_orbit to check interpolation convergence. O-C uses
    the initial nominal conjunction and initial osculating period, NOT a fitted
    linear ephemeris. TDV/depth variations are relative to the first event.
    """
    planets = np.asarray(planets, float).reshape(-1, 7)
    bounds = np.asarray(times, float)
    if bounds.shape != (2,) or not np.all(np.isfinite(bounds)) or bounds[1] <= bounds[0]:
        raise ValueError('Supply finite start/end times with end > start.')
    if int(samples_per_orbit) != samples_per_orbit or samples_per_orbit < 32:
        raise ValueError('samples_per_orbit must be an integer >=32.')
    if int(event_curve_points) != event_curve_points or (event_curve_points != 0 and event_curve_points < 3):
        raise ValueError('event_curve_points must be zero or an integer >=3.')
    if not np.isfinite(max_samples) or max_samples != int(max_samples) or max_samples < 2:
        raise ValueError('max_samples must be an integer >=2.')
    labels = np.atleast_1d(dataset)
    if labels.size == 0 or not np.all(np.isfinite(labels)) or np.any(labels != labels.astype(int)) or np.any((labels < 0) | (labels > 59)):
        raise ValueError('Dataset slots must be integers in 0..59.')
    labels = np.unique(labels.astype(int))
    # Reuse all existing parameter/LD/physical consistency checks and phase conversion.
    epoch = float(options.get('epoch', 0.))
    if options.get('rho_star') is not None:
        from native_transit_model import density_geometry
        _, options['stellar_radius'] = density_geometry(planets[:, 1], options['rho_star'], options.get('stellar_mass', 1.))
    checked = evaluate(planets, [epoch], dataset=int(labels[0]), **options)
    mode = options.get('model', 'keplerian')
    radii = np.broadcast_to(np.asarray(options['radius_ratio'], float), (len(planets),))
    if not len(planets):
        return {}
    elements = planets[:, 1:].copy()
    elements[:, 3] = checked.mean_anomaly
    masses = _mass_axes(planets, options.get('stellar_mass', 1.))[0] if mode == 'nbody' else np.zeros(len(planets))
    order = np.argsort(elements[:, 0], kind='stable')
    pmin = np.min(elements[:, 0])
    dynamical_scale = np.min(elements[:, 0]*(1-elements[:, 1])**1.5/np.sqrt(1+elements[:, 1]))
    cadence = dynamical_scale/samples_per_orbit
    lo, hi = bounds-epoch
    padding = np.max(elements[:, 0])/2
    count = int(np.ceil((hi-lo+2*padding)/cadence))+1
    if count > max_samples:
        raise ValueError(f'Event search needs {count} trajectory samples; narrow the interval or raise max_samples.')
    grid = np.linspace(lo-padding, hi+padding, count)
    with _LOCK:
        backend = (_native().rv_transit_geometry_ttvfast
                   if mode == 'nbody' and options.get('jacobi_convention', 'standard') == 'ttvfast'
                   else _native().rv_transit_geometry)
        geometry, status = backend(
            grid, elements[order], masses[order], checked.scaled_axis[order],
            float(options.get('stellar_mass', 1.)), float(options['stellar_radius']),
            int(mode == 'nbody'), float(options.get('step') or pmin/50),
            float(options.get('tolerance', 1.e-12)), int(options.get('max_steps', 1000000)))
    if status:
        raise RuntimeError(f'Native transit event trajectory failed (status {status}).')
    geometry = geometry[:, np.argsort(order), :]
    coefficients = np.zeros((60, 4)) if options.get('coefficients') is None else np.asarray(options['coefficients'], float)
    laws = options.get('limb_dark')
    if laws is None:
        laws = ['uniform']*60
    if isinstance(laws, str):
        laws = [laws]*60
    q = int(options.get('quadrature_order', 128))
    output = {}
    for j in range(len(planets)):
        trajectory = CubicHermiteSpline(grid, geometry[:, j, :3], geometry[:, j, 3:])
        def separation(t):
            r = trajectory(t)
            return np.linalg.norm(r[..., :2], axis=-1)
        def derivative(t):
            r = trajectory(t); v = trajectory(t, 1)
            return np.sum(r[..., :2]*v[..., :2], axis=-1)
        g = derivative(grid)
        candidates = np.flatnonzero((g[:-1] <= 0) & (g[1:] > 0))
        events = []
        for k in candidates:
            mid = brentq(derivative, grid[k], grid[k+1], xtol=1.e-12)
            impact = float(separation(mid))
            if not (lo-1.e-10 <= mid <= hi+1.e-10) or trajectory(mid)[2] <= 0 or impact >= 1+radii[j]:
                continue
            def contact(limit, direction):
                if impact >= limit:
                    return np.nan
                indices = range(k+1, len(grid)) if direction > 0 else range(k, -1, -1)
                for index in indices:
                    edge = grid[index]
                    if direction*(edge-mid) <= 0:
                        continue
                    if separation(edge) >= limit:
                        return brentq(lambda t: separation(t)-limit, min(mid, edge), max(mid, edge), xtol=1.e-12)
                raise RuntimeError('Transit contact is outside the padded trajectory; narrow the search or inspect the orbit.')
            t1 = contact(1+radii[j], -1); t4 = contact(1+radii[j], 1)
            t2 = contact(abs(1-radii[j]), -1); t3 = contact(abs(1-radii[j]), 1)
            events.append((mid, t1, t2, t3, t4, impact, float(np.linalg.norm(trajectory(mid, 1)[:2]))))
        values = np.asarray(events, float).reshape(-1, 7)
        relative_mid = values[:, 0]
        # Nominal conjunction reference calculated with the existing Fortran phase convention.
        e, w = elements[j, 1], np.deg2rad(elements[j, 2])
        f = np.pi/2-w
        ea = 2*np.arctan2(np.sqrt(1-e)*np.sin(f/2), np.sqrt(1+e)*np.cos(f/2))
        reference = ((ea-e*np.sin(ea)-np.deg2rad(elements[j, 3])) % (2*np.pi))*elements[j, 0]/(2*np.pi)
        number = np.rint((relative_mid-reference)/elements[j, 0]).astype(int)
        duration = values[:, 4]-values[:, 1]
        event = dict(epoch=number, t0=relative_mid+epoch, t1=values[:, 1]+epoch,
                     t2=values[:, 2]+epoch, t3=values[:, 3]+epoch, t4=values[:, 4]+epoch,
                     duration=duration, duration_full=values[:, 3]-values[:, 2],
                     impact_parameter=values[:, 5], sky_speed=values[:, 6],
                     linear_t0=epoch+reference+number*elements[j, 0],
                     ttv=relative_mid-reference-number*elements[j, 0],
                     tdv=duration-duration[0] if len(duration) else duration.copy(),
                     reference_t0=epoch+reference, reference_period=elements[j, 0],
                     depth={}, depth_at_t0={}, depth_variation={}, depth_time={}, models={})
        for d in labels:
            def flux(t):
                with _LOCK:
                    result, status = _native().rv_transit_event_flux(
                        np.atleast_1d(separation(t)), float(radii[j]), LD_LAWS[laws[d]], coefficients[d], q)
                if status:
                    raise RuntimeError('Invalid limb-darkening profile for event depth.')
                return result
            depths = []; peaks = []; curves = []
            for mid, t1, _, _, t4, _, _ in events:
                # A coarse search also handles limb-brightened profiles where
                # the deepest point need not coincide with geometric mid-transit.
                probes = np.linspace(t1, t4, 65); sampled = flux(probes)
                m = int(np.argmin(sampled)); left = probes[max(0, m-1)]; right = probes[min(64, m+1)]
                center = (left+right)/2
                minimum = minimize_scalar(lambda dt: float(flux(center+dt)[0]),
                                          bounds=(left-center, right-center), method='bounded',
                                          options={'xatol': 1.e-12})
                choices = [(float(sampled[m]), probes[m]), (float(minimum.fun), center+minimum.x), (float(flux(mid)[0]), mid)]
                value, peak = min(choices)
                depths.append(1-value); peaks.append(peak+epoch)
                if event_curve_points:
                    ts = np.linspace(t1, t4, event_curve_points)
                    curves.append(dict(times=ts+epoch, flux=flux(ts)))
            event['depth'][int(d)] = np.asarray(depths)
            event['depth_time'][int(d)] = np.asarray(peaks)
            event['depth_at_t0'][int(d)] = 1-flux(relative_mid) if len(events) else np.empty(0)
            event['depth_variation'][int(d)] = np.asarray(depths)-depths[0] if depths else np.empty(0)
            event['models'][int(d)] = curves
        output[j] = event
    return output


def update_transit_events(session, *, start=None, end=None, dataset=None, **options):
    """Publish optional event tables, keyed by original zero-based planet slot."""
    from native_model_results import model_time_grid
    from native_transit_model import from_session
    targets = [session.fit_results]
    if getattr(session, 'rvmod', None) is not None and session.rvmod is not session.fit_results:
        targets.append(session.rvmod)
    for target in targets:
        target.transit_events = {}
        target.transit_events_status = 'disabled'
        target.ttv_model = np.empty((0, 5))
    if not getattr(session, 'calculate_transit_events', False):
        return None
    if start is None or end is None:
        bounds = model_time_grid(session, float(options.get('epoch', session.epoch)), component='all')
        start = bounds[0] if start is None else start
        end = bounds[-1] if end is None else end
    if dataset is None:
        dataset = [d for d in range(60) if len(session.tra_data_sets[d])] or [0]
    for target in targets:
        target.transit_events_status = 'not calculated'
    try:
        events = from_session(session, [start, end], dataset=dataset, _evaluator=evaluate_events, **options)
    except Exception:
        for target in targets:
            target.transit_events_status = 'failed'
        raise
    slots = np.flatnonzero(session.use_planet)
    mapped = {int(slots[j]): table for j, table in events.items()}
    rows = [np.column_stack((np.full(len(t['t0']), slot+1), t['epoch'], t['t0'],
                             t['impact_parameter'], t['sky_speed'])) for slot, t in mapped.items()]
    table = np.concatenate(rows) if rows else np.empty((0, 5))
    table = table[np.argsort(table[:, 2], kind='stable')]
    for target in targets:
        target.transit_events = mapped
        target.ttv_model = table.copy()
        target.transit_events_status = 'ready'
    return mapped
