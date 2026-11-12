"""Dataset-free predictions from the RVmod Fortran models.

No session, fake data, fitted offsets or optimizer calls are required.
"""
from collections.abc import Mapping
from dataclasses import dataclass, replace

import numpy as np

from native_transit_model import evaluate, density_geometry, _native, _LOCK, _mass_axes, LD_LAWS


def _validate_rv_orbits(planets, stellar_mass):
    """Validate physical inputs before the iterative mass conversion."""
    if not np.isscalar(stellar_mass) or not np.isfinite(stellar_mass) or stellar_mass <= 0:
        raise ValueError('stellar_mass must be finite and positive.')
    if (np.any(planets[:, 0] < 0) or np.any(planets[:, 1] <= 0)
            or np.any(planets[:, 2] < 0) or np.any(planets[:, 2] >= 1)
            or np.any(planets[:, 5] < 0) or np.any(planets[:, 5] > 180)):
        raise ValueError('Require K>=0, P>0, 0<=e<1 and inclination in 0..180 degrees.')


def _standard_parameters(standard):
    """Pack a named five-parameter astrometric solution in ExTRA order."""
    if isinstance(standard, Mapping):
        fields = ('RA', 'Dec', 'parallax', 'mu_alpha_star', 'mu_delta')
        missing = set(fields)-standard.keys()
        unknown = standard.keys()-set(fields)
        if missing:
            raise ValueError(f'Astrometry standard: missing parameters: {", ".join(sorted(missing))}.')
        if unknown:
            raise ValueError(f'Astrometry standard: unknown parameters: {", ".join(sorted(unknown))}.')
        standard = [standard[field] for field in fields]
    values = np.asarray(standard, float)
    if values.shape != (5,) or not np.all(np.isfinite(values)):
        raise ValueError('Astrometry standard must contain five finite scalar values.')
    return values


def time_grid(times=None, *, start=None, end=None, points=1000, before=0., after=0.):
    """Return explicit times or an inclusive grid padded by days on either side.

    With explicit times and no grid options, preserve order and duplicates.
    With start/end, points controls plot density independently of integration.
    """
    if not np.all(np.isfinite([before, after])) or before < 0 or after < 0:
        raise ValueError('before and after must be finite nonnegative days.')
    if times is not None:
        values = np.asarray(times, float)
        if values.ndim != 1 or not values.size or not np.all(np.isfinite(values)):
            raise ValueError('times must be a nonempty finite one-dimensional array.')
        if start is None and end is None and before == 0 and after == 0:
            return values.copy()
        start = values.min() if start is None else start
        end = values.max() if end is None else end
    if start is None or end is None:
        raise ValueError('Supply times or both start and end; no data bounds are assumed.')
    if not np.all(np.isfinite([start, end])) or end < start:
        raise ValueError('Require finite end >= start.')
    if not np.isscalar(points) or not np.isfinite(points) or points != int(points) or points < 2:
        raise ValueError('points must be an integer >=2.')
    return np.linspace(start-before, end+after, int(points))


def _evaluate_model(planets, times=None, *, start=None, end=None, points=1000, before=0., after=0.,
          epoch=0., model='keplerian', timing=False, max_events=10000,
          samples_per_orbit=512, max_samples=2000000, **options):
    """Return RV, flux and optional conjunction/TTV tables without observations.

    planets: up to nine [K(m/s), P(days), e, omega(deg), M0(deg), i(deg),
    node(deg)] rows. Their count defines the planet count. epoch and all times
    use the same day scale. Other options follow native_transit_model.evaluate.
    RV-only calls default to zero planet radii and physically derived axes.
    Timing reports foreground conjunctions (including non-transiting orbits),
    with impact distance in AU; it never claims every conjunction is a transit.
    """
    grid = time_grid(times, start=start, end=end, points=points, before=before, after=after)
    planets = np.asarray(planets, float)
    if planets.size == 0:
        planets = planets.reshape(0, 7)
    if planets.ndim != 2 or planets.shape[1] != 7 or len(planets) > 9:
        raise ValueError('Expected up to nine [K,P,e,omega,M0,i,node] rows.')
    if not np.all(np.isfinite(planets)) or not np.isfinite(epoch):
        raise ValueError('Planet parameters and epoch must be finite.')
    options = options.copy()
    options.setdefault('radius_ratio', 0.)
    stellar_mass = options.get('stellar_mass', 1.)
    _validate_rv_orbits(planets, stellar_mass)
    if options.get('rho_star') is not None:
        _, options['stellar_radius'] = density_geometry(planets[:, 1], options['rho_star'], stellar_mass)
    if model == 'keplerian' and options.get('scaled_axis') is None and options.get('rho_star') is None:
        _, axes = _mass_axes(planets, stellar_mass)
        options['scaled_axis'] = axes*1.49597892e11/(options.get('stellar_radius', 1.)*6.957e8)
    prediction = evaluate(planets, grid, epoch=epoch, model=model, **options)
    events = {}
    if timing and len(planets):
        if grid.max() <= grid.min():
            raise ValueError('Timing search requires a nonzero time span.')
        for name, value, minimum in [('max_events', max_events, 1),
                                     ('samples_per_orbit', samples_per_orbit, 32),
                                     ('max_samples', max_samples, 4)]:
            if not np.isscalar(value) or not np.isfinite(value) or value != int(value) or value < minimum:
                raise ValueError(f'{name} must be an integer >={minimum}.')
        elements = planets[:, 1:].copy()
        elements[:, 3] = prediction.mean_anomaly
        order = np.argsort(elements[:, 0], kind='stable')
        masses = _mass_axes(planets, stellar_mass)[0] if model == 'nbody' else np.zeros(len(planets))
        step = options.get('step') or np.min(planets[:, 1])/50
        with _LOCK:
            counts, table, status = _native().timing_model(
                elements[order], masses[order], float(stellar_mass), float(epoch),
                float(grid.min()), float(grid.max()), int(model == 'nbody'),
                prediction.scaled_axis[order]*options.get('stellar_radius', 1.),
                float(step), float(options.get('tolerance', 1.e-12)),
                int(options.get('max_steps', 1000000)),
                int(options.get('jacobi_convention', 'standard') == 'ttvfast'),
                int(samples_per_orbit), int(max_samples), int(max_events))
        if status:
            raise ValueError(f'Native timing model failed (status {status}); check search limits and parameters.')
        for k, slot in enumerate(order):
            values = table[:counts[k], k, :].copy()
            # Reference the first conjunction in this interval, without fitting an ephemeris.
            linear = values[0, 1]+np.arange(len(values))*planets[slot, 1] if len(values) else np.empty(0)
            events[int(slot)] = dict(time=values[:, 1], number=values[:, 0].astype(int),
                                     impact_au=values[:, 2], speed_au_per_day=values[:, 3],
                                     linear_time=linear, ttv=values[:, 1]-linear)
    return dict(times=prediction.times, rv=prediction.rv, flux=prediction.flux,
                rv_components=prediction.rv_components, events=events, prediction=prediction)


# Preserve the initial array-based API for existing callers.
model = _evaluate_model


def astrometric_model(orbits, times=None, *, standard, epoch=2457389., observer=None,
                      scan_angle=0., model='keplerian', stellar_mass=1., step=None,
                      tolerance=1.e-12, max_steps=1000000, jacobi_convention='standard',
                      start=None, end=None, points=1000, before=0., after=0.):
    """Return ExTRA-compatible east/north/along-scan predictions [mas].

    Keplerian orbits: [P(day), e, omega(rad), i(rad), node(rad), Tperi(JD),
    stellar semimajor axis(mas)], exactly ExTRA.orbit_total's convention.
    N-body orbits: RV rows [K,P,e,omega,M0,i,node], angles in degrees;
    masses and stellar reflex motion are derived from K and stellar_mass.
    standard: [RA(deg), Dec(deg), parallax(mas), mu_alpha_star, mu_delta(mas/yr)].
    observer: (N,3) barycentric AU, allowing Gaia's own positions. When omitted,
    ExTRA supplies Earth positions for nonzero parallax. Scan angle is radians,
    north through east. RA parallax follows ExTRA's existing convention.
    """
    grid = time_grid(times, start=start, end=end, points=points, before=before, after=after)
    standard = _standard_parameters(standard)
    rows = np.asarray(orbits, float)
    if rows.size == 0:
        rows = rows.reshape(0, 7)
    if rows.ndim != 2 or rows.shape[1] != 7 or len(rows) > 9 or not np.all(np.isfinite(rows)):
        raise ValueError('orbits must contain up to nine finite seven-parameter rows.')
    if standard.shape != (5,) or not np.all(np.isfinite(standard)) or not np.isfinite(epoch):
        raise ValueError('Supply five finite standard parameters and a finite epoch.')
    if observer is None:
        if standard[2] != 0:
            from ExTRA.astrometry import earth_position
            observer = earth_position(grid).T
        else:
            observer = np.zeros((len(grid), 3))
    observer = np.asarray(observer, float)
    angles = np.broadcast_to(np.asarray(scan_angle, float), grid.shape).copy()
    if observer.shape != (len(grid), 3) or not np.all(np.isfinite(observer)) or not np.all(np.isfinite(angles)):
        raise ValueError('observer must have finite shape (N,3); scan angles must be finite.')
    if model not in ('keplerian', 'nbody'):
        raise ValueError("model must be 'keplerian' or 'nbody'.")
    with _LOCK:
        if model == 'nbody' and len(rows):
            _validate_rv_orbits(rows, stellar_mass)
            if jacobi_convention not in ('standard', 'ttvfast'):
                raise ValueError('Unknown Jacobi convention.')
            step = np.min(rows[:, 1])/50 if step is None else step
            if not np.isfinite(step) or step <= 0 or not np.isfinite(tolerance) or not 1.e-15 <= tolerance <= 1.e-3:
                raise ValueError('Invalid integration step or tolerance.')
            if not np.isfinite(max_steps) or int(max_steps) != max_steps or max_steps < 1:
                raise ValueError('max_steps must be a positive integer.')
            masses, axes = _mass_axes(rows, stellar_mass)
            order = np.argsort(grid, kind='stable')
            bodies = np.argsort(rows[:, 1], kind='stable')
            east, north, along, parts, status = _native().astrometry_nbody(
                grid[order], rows[bodies, 1:], masses[bodies], axes[bodies], float(stellar_mass),
                standard, observer[order], float(epoch), angles[order], float(step), float(tolerance),
                int(max_steps), int(jacobi_convention == 'ttvfast'))
            inverse = np.argsort(order)
            east, north, along = east[inverse], north[inverse], along[inverse]
            parts = parts[inverse][:, np.argsort(bodies)]
        else:
            east, north, along, parts, status = _native().astrometry_model(
                grid, rows, standard, observer, float(epoch), angles)
    if status:
        raise ValueError(f'Native astrometry model failed (status {status}).')
    return dict(times=grid, east=east, north=north, along_scan=along, components=parts)


class _Fields:
    """Allow both result.RV.BJD and result['RV']['BJD'] access."""

    def __getitem__(self, name):
        return getattr(self, name)


@dataclass(frozen=True)
class RVModel(_Fields):
    """Sampled stellar RV and named planetary contributions in m/s."""

    BJD: np.ndarray
    model: np.ndarray
    components: dict
    data: object = None
    uncertainty: object = None
    noise: object = None
    requested_BJD: object = None


@dataclass(frozen=True)
class TransitModel(_Fields):
    """Sampled combined normalized transit flux."""

    BJD: np.ndarray
    model: np.ndarray
    data: object = None
    uncertainty: object = None
    noise: object = None
    requested_BJD: object = None


@dataclass(frozen=True)
class TransitEvents(_Fields):
    """Actual foreground transits of one planet; all durations are days."""

    N_transit: np.ndarray
    BJD: np.ndarray
    BJD_minus_mean: np.ndarray
    Mean_BJD: float
    Depth: np.ndarray
    Duration: np.ndarray
    Duration_full: np.ndarray
    T1: np.ndarray
    T2: np.ndarray
    T3: np.ndarray
    T4: np.ndarray
    Linear_BJD: np.ndarray
    O_C: np.ndarray
    Impact_parameter: np.ndarray
    Sky_speed: np.ndarray
    Depth_at_midtransit: np.ndarray
    Depth_BJD: np.ndarray
    TDV: np.ndarray
    Depth_variation: np.ndarray
    data: object = None
    uncertainty: object = None
    noise: object = None


@dataclass(frozen=True)
class AstrometryModel(_Fields):
    """Sky offsets and named stellar-reflex contributions in mas."""

    BJD: np.ndarray
    East: np.ndarray
    North: np.ndarray
    along_scan: np.ndarray
    components: dict
    data: object = None
    uncertainty: object = None
    noise: object = None
    requested_BJD: object = None


@dataclass(frozen=True)
class ModelResult(_Fields):
    """Group dataset-free observable curves and per-planet transit tables."""

    RV: object
    Transit: object
    TTV: dict
    Astrometry: object = None
    settings: object = None


def _planet_parameters(planets, options):
    """Pack named planets in caller order, retaining names for result tables."""
    fields = ('K', 'P', 'e', 'omega', 'M0', 'inclination', 'node')
    extra = ('radius_ratio', 'scaled_axis')
    if isinstance(planets, Mapping):
        if 'K' in planets or 'P' in planets:
            named = [('planet_1', planets)]
        else:
            named = list(planets.items())
    else:
        values = list(planets)
        if not values or not isinstance(values[0], Mapping):
            rows = np.asarray(values, float)
            if not values:
                rows = rows.reshape(0, 7)
            return tuple(f'planet_{i+1}' for i in range(len(rows))), rows, options.copy()
        named = [(f'planet_{i+1}', value) for i, value in enumerate(values)]
    rows = []
    packed_options = options.copy()
    for name, parameters in named:
        if not isinstance(name, str) or not name:
            raise ValueError('Planet names must be nonempty strings.')
        if not isinstance(parameters, Mapping):
            raise ValueError(f'Planet {name}: supply a parameter dictionary.')
        missing = set(fields)-parameters.keys()
        unknown = parameters.keys()-set(fields+extra)
        if missing:
            raise ValueError(f'Planet {name}: missing parameters: {", ".join(sorted(missing))}.')
        if unknown:
            raise ValueError(f'Planet {name}: unknown parameters: {", ".join(sorted(unknown))}.')
        if any(np.ndim(parameters[field]) != 0 for field in parameters):
            raise ValueError(f'Planet {name}: each parameter must be a scalar.')
        rows.append([parameters[field] for field in fields])
    for field in extra:
        present = [field in parameters for _, parameters in named]
        if any(present):
            if field in packed_options:
                raise ValueError(f'Supply {field} in planet dictionaries or as an option, not both.')
            if field == 'scaled_axis' and not all(present):
                raise ValueError('Supply scaled_axis for every planet or omit it for every planet.')
            packed_options[field] = [parameters.get(field, 0.) for _, parameters in named]
    return tuple(name for name, _ in named), np.asarray(rows, float).reshape(-1, 7), packed_options


def _transit_tables(planets, names, prediction, epoch, mode, options, profile,
                    samples_per_orbit, max_samples, max_events):
    """Read complete transit diagnostics from the native event model."""
    for name, value, minimum in [('samples_per_orbit', samples_per_orbit, 32),
                                 ('max_samples', max_samples, 4), ('max_events', max_events, 1)]:
        if (not np.isscalar(value) or not np.isfinite(value) or value != int(value)
                or value < minimum or value > np.iinfo(np.int32).max):
            raise ValueError(f'{name} must be a 32-bit integer >={minimum}.')
    if not np.isscalar(profile) or not np.isfinite(profile) or int(profile) != profile or not 0 <= profile < 60:
        raise ValueError('event_dataset must be an integer in 0..59.')
    profile = int(profile)
    stellar_mass = options.get('stellar_mass', 1.)
    stellar_radius = options.get('stellar_radius', 1.)
    if options.get('rho_star') is not None:
        _, stellar_radius = density_geometry(planets[:, 1], options['rho_star'], stellar_mass)
    radii = np.broadcast_to(np.asarray(options.get('radius_ratio', 0.), float), (len(planets),))
    elements = planets[:, 1:].copy()
    elements[:, 3] = prediction.mean_anomaly
    order = np.argsort(planets[:, 1], kind='stable')
    masses = _mass_axes(planets, stellar_mass)[0] if mode == 'nbody' else np.zeros(len(planets))
    laws = options.get('limb_dark')
    law = LD_LAWS[laws[profile]] if laws is not None else 0
    coeff = np.zeros(4)
    if options.get('coefficients') is not None:
        coeff = np.asarray(options['coefficients'], float)[profile]
    counts = np.zeros(len(planets), int)
    table = np.empty((0, len(planets), 12))
    if len(planets) and np.any(radii > 0):
        if prediction.times.max() <= prediction.times.min():
            raise ValueError('Transit events require a nonzero time span; use transit_events=False for one time.')
        with _LOCK:
            counts, table, status = _native().transit_event_model(
                elements[order], masses[order], prediction.scaled_axis[order], radii[order],
                float(stellar_mass), float(stellar_radius), float(epoch),
                float(prediction.times.min()), float(prediction.times.max()), int(mode == 'nbody'),
                float(options.get('step') or np.min(planets[:, 1])/50),
                float(options.get('tolerance', 1.e-12)), int(options.get('max_steps', 1000000)),
                int(options.get('jacobi_convention', 'standard') == 'ttvfast'),
                int(samples_per_orbit), int(max_samples), int(max_events), law, coeff,
                int(options.get('quadrature_order', 128)))
        if status:
            errors = {9: 'Transit search exceeds max_samples or its memory limit.',
                      10: 'Transit count exceeds max_events.',
                      11: 'Transit contacts lie outside the padded search interval.'}
            raise ValueError(errors.get(status, f'Native transit event model failed (status {status}).'))
    output = {}
    for k, slot in enumerate(order):
        values = table[:counts[k], k, :]
        relative = values[:, 0]
        mean = float(np.mean(relative)) if len(relative) else np.nan
        # O-C uses the first returned event and initial period; retain skipped orbits.
        cycles = np.rint((relative-relative[0])/planets[slot, 1]) if len(relative) else np.empty(0)
        linear = relative[0]+cycles*planets[slot, 1] if len(relative) else np.empty(0)
        output[names[slot]] = TransitEvents(
            N_transit=np.arange(1, len(relative)+1), BJD=relative+epoch,
            BJD_minus_mean=relative-mean, Mean_BJD=mean+epoch,
            Depth=values[:, 9], Duration=values[:, 5], Duration_full=values[:, 6],
            T1=values[:, 1]+epoch, T2=values[:, 2]+epoch,
            T3=values[:, 3]+epoch, T4=values[:, 4]+epoch,
            Linear_BJD=linear+epoch, O_C=relative-linear,
            Impact_parameter=values[:, 7], Sky_speed=values[:, 8],
            Depth_at_midtransit=values[:, 10], Depth_BJD=values[:, 11]+epoch,
            TDV=values[:, 5]-values[0, 5] if len(relative) else np.empty(0),
            Depth_variation=values[:, 9]-values[0, 9] if len(relative) else np.empty(0))
    # Keep caller order even when the integrator sorts planets by period.
    return {name: output[name] for name in names}



def GetModelDefaults():
    """Return fresh editable model dictionaries; only RV is enabled initially.

    None sampling values inherit GetModel's shared sampling arguments.
    Noise sigma uses m/s (RV), relative flux (Transit), days (TTV), or mas
    (Astrometry along-scan). Time jitter sigma is always in days.
    """
    settings = {}
    for name in ('RV', 'Transit', 'TTV', 'Astrometry'):
        entry = dict(enabled=name == 'RV', BJD=None, start=None, end=None,
                     points=None, before=None, after=None, parameters={},
                     noise=dict(sigma=0., seed=None))
        if name != 'TTV':
            entry['time_jitter'] = dict(sigma=0., seed=None)
        settings[name] = entry
    return settings


def _seed_sequence(seed, name):
    """Reject ambiguous seeds while preserving caller/global random state."""
    if seed is not None:
        if isinstance(seed, (bool, np.bool_)) or not isinstance(seed, (int, np.integer)) or seed < 0:
            raise ValueError(f'{name} must be a nonnegative integer or None.')
        seed = int(seed)
    return np.random.SeedSequence(seed)


def _model_settings(models, transit_events):
    """Merge documented defaults and reject misspelled configuration keys."""
    settings = GetModelDefaults()
    if models is not None:
        if not isinstance(models, Mapping):
            raise ValueError('models must be a dictionary keyed by RV, Transit, TTV, Astrometry.')
        unknown = models.keys()-settings.keys()
        if unknown:
            raise ValueError(f'Unknown models: {", ".join(sorted(unknown))}.')
        for name, supplied in models.items():
            if isinstance(supplied, (bool, np.bool_)):
                settings[name]['enabled'] = bool(supplied)
                continue
            if not isinstance(supplied, Mapping):
                raise ValueError(f'models[{name}] must be a boolean or settings dictionary.')
            unknown = supplied.keys()-settings[name].keys()
            if unknown:
                raise ValueError(f'Model {name}: unknown settings: {", ".join(sorted(unknown))}.')
            for key, value in supplied.items():
                if key in ('noise', 'time_jitter'):
                    if not isinstance(value, Mapping) or value.keys()-{'sigma', 'seed'}:
                        raise ValueError(f'{name}.{key} must contain only sigma and seed.')
                    settings[name][key].update(value)
                elif key == 'parameters':
                    if not isinstance(value, Mapping):
                        raise ValueError(f'{name}.parameters must be a dictionary.')
                    settings[name][key] = dict(value)
                else:
                    settings[name][key] = value
            if not isinstance(settings[name]['enabled'], (bool, np.bool_)):
                raise ValueError(f'{name}.enabled must be True or False.')
    if transit_events is not None:
        if not isinstance(transit_events, (bool, np.bool_)):
            raise ValueError('transit_events must be True or False.')
        if models is not None and 'TTV' in models and settings['TTV']['enabled'] != transit_events:
            raise ValueError('Conflicting TTV.enabled and transit_events switches.')
        settings['TTV']['enabled'] = bool(transit_events)
    return settings


def _sigma_array(sigma, shape, name):
    """Validate scalar or per-sample Gaussian standard deviations."""
    try:
        supplied = np.asarray(sigma, float)
    except (TypeError, ValueError) as error:
        raise ValueError(f'{name} must contain numeric standard deviations.') from error
    if not np.all(np.isfinite(supplied)) or np.any(supplied < 0):
        raise ValueError(f'{name} must contain finite nonnegative standard deviations.')
    try:
        values = np.broadcast_to(supplied, shape).copy()
    except ValueError as error:
        raise ValueError(f'{name} must be scalar or match the sample shape {shape}.') from error
    return values


def _gaussian_noise(values, settings, stream, name):
    """Return synthetic measurements and errors without changing the model."""
    sigma = _sigma_array(settings['sigma'], values.shape, f'{name}.noise.sigma')
    sequence = stream
    if settings['seed'] is not None:
        sequence = _seed_sequence(settings['seed'], f'{name}.noise.seed')
    rng = np.random.default_rng(sequence)
    noise = rng.normal(0., sigma)
    data = values+noise
    if not np.all(np.isfinite(noise)) or not np.all(np.isfinite(data)):
        raise ValueError(f'{name}: Gaussian simulation produced nonfinite data.')
    return data, sigma, noise


def _sample_times(settings, shared, stream, name):
    """Resolve an observable's explicit BJD or grid, then optional time scatter."""
    values = {key: settings[key] if settings[key] is not None else shared[key]
              for key in ('BJD', 'start', 'end', 'points', 'before', 'after')}
    if settings['BJD'] is None and (settings['start'] is not None or settings['end'] is not None):
        values['BJD'] = None
    if values['BJD'] is not None:
        if settings['start'] is not None or settings['end'] is not None:
            raise ValueError(f'{name}: use explicit BJD or start/end, not both.')
        requested = np.asarray(values['BJD'], float)
        if requested.ndim == 0:
            requested = requested.reshape(1)
        if requested.ndim != 1 or not requested.size or not np.all(np.isfinite(requested)):
            raise ValueError(f'{name}.BJD must be a finite timestamp or nonempty 1-D array.')
        requested = requested.copy()
    else:
        requested = time_grid(start=values['start'], end=values['end'], points=values['points'],
                              before=values['before'], after=values['after'])
    actual = requested.copy()
    if name != 'TTV':
        jitter = settings['time_jitter']
        sigma = _sigma_array(jitter['sigma'], requested.shape, f'{name}.time_jitter.sigma')
        sequence = stream
        if jitter['seed'] is not None:
            sequence = _seed_sequence(jitter['seed'], f'{name}.time_jitter.seed')
        actual += np.random.default_rng(sequence).normal(0., sigma)
        if not np.all(np.isfinite(actual)):
            raise ValueError(f'{name}: timestamp jitter produced nonfinite BJD.')
    return requested, actual


def _profile_options(options):
    """Accept one limb profile or the existing 60-profile native arrays."""
    prepared = options.copy()
    law = prepared.get('limb_dark')
    if isinstance(law, str):
        if law not in LD_LAWS:
            raise ValueError(f'Unknown limb-darkening law: {law}.')
        prepared['limb_dark'] = [law]*60
    coefficients = prepared.get('coefficients')
    if coefficients is not None:
        values = np.asarray(coefficients, float)
        if values.ndim == 1 and len(values) <= 4:
            expanded = np.zeros((60, 4))
            expanded[:, :len(values)] = values
            prepared['coefficients'] = expanded
    return prepared


def _observable_options(common, specific, rv_only=False):
    """Resolve model parameters; RV-only predictions have no occulting disks."""
    options = common.copy()
    options.update(specific)
    if rv_only:
        options['radius_ratio'] = 0.
        for name in ('limb_dark', 'coefficients', 'offsets', 'jitters'):
            options.pop(name, None)
    return _profile_options(options)


def _astrometry_curve(rows, names, grid, epoch, mode, options, supplied, mean_anomaly):
    """Derive named stellar reflex orbits and evaluate the astrometric model."""
    settings = dict(supplied)
    settings.setdefault('epoch', epoch)
    settings.setdefault('model', mode)
    settings.setdefault('stellar_mass', options.get('stellar_mass', 1.))
    if 'standard' not in settings:
        raise ValueError('Enable Astrometry with parameters.standard: [RA, Dec, parallax, mu_alpha_star, mu_delta].')
    standard = _standard_parameters(settings['standard'])
    settings['standard'] = standard
    if 'orbits' not in settings:
        if settings['model'] == 'nbody':
            orbit_rows = rows.copy()
            orbit_rows[:, 4] = mean_anomaly
            settings['orbits'] = orbit_rows
        else:
            masses, axes = _mass_axes(rows, settings['stellar_mass'])
            order = np.argsort(rows[:, 1], kind='stable')
            total = np.empty(len(rows))
            total[order] = settings['stellar_mass']+np.cumsum(masses[order])
            amplitude = axes*masses/total*standard[2]
            # ExTRA describes the stellar orbit: reverse companion omega.
            settings['orbits'] = np.column_stack((
                rows[:, 1], rows[:, 2], np.deg2rad(rows[:, 3]+180.),
                np.deg2rad(rows[:, 5]), np.deg2rad(rows[:, 6]),
                epoch-mean_anomaly/360.*rows[:, 1], amplitude))
    result = astrometric_model(times=grid, **settings)
    count = result['components'].shape[1]
    astro_names = names if count == len(names) else tuple(f'planet_{i+1}' for i in range(count))
    components = {name: result['components'][:, i, :] for i, name in enumerate(astro_names)}
    return AstrometryModel(result['times'], result['east'], result['north'],
                           result['along_scan'], components)


def GetModel(planets, times=None, *, BJD=None, models=None, seed=None,
             start=None, end=None, points=1000, before=0., after=0., epoch=0.,
             model='keplerian', transit_events=None, event_dataset=0,
             samples_per_orbit=512, max_samples=2000000, max_events=10000,
             astrometry=None, **options):
    """Predict selected observables and optional Gaussian synthetic measurements.

    planets is a named dictionary of K, P, e, omega, M0, inclination, node and
    optional radius_ratio/scaled_axis. K is m/s, P is days, angles are degrees.
    models is an editable dictionary from GetModelDefaults(): RV defaults on;
    Transit, TTV and Astrometry default off. Entries accept enabled, BJD or grid
    settings, parameters and noise={sigma, seed}. Sampled curves also accept
    time_jitter={sigma, seed}, with sigma in days. Top-level seed initializes
    independent streams, so toggling an observable does not alter other draws.

    Explicit BJD (or the times alias) accepts one epoch or an irregular array.
    Each enabled observable can override it with its own BJD. Arrays retain
    order and duplicates. model contains noiseless predictions; data contains
    predictions plus Gaussian noise, with uncertainty and noise returned too.
    TTV.data simulates measured mid-transit BJD; Astrometry.data simulates
    along-scan offsets. Disabled sampled structures are None; disabled TTV is
    {}. Supplied astrometry parameters do not enable its switch automatically.
    """
    if BJD is not None and times is not None:
        raise ValueError('Supply BJD or times, not both.')
    if (BJD is not None or times is not None) and (start is not None or end is not None):
        raise ValueError('Supply explicit BJD/times or shared start/end, not both.')
    settings = _model_settings(models, transit_events)
    streams = _seed_sequence(seed, 'seed').spawn(8)
    shared = dict(BJD=BJD if BJD is not None else times, start=start, end=end,
                  points=points, before=before, after=after)
    names, rows, options = _planet_parameters(planets, options)
    if rows.ndim != 2 or rows.shape[1] != 7 or len(rows) > 9:
        raise ValueError('Expected up to nine seven-parameter planets.')
    if not np.all(np.isfinite(rows)) or not np.isfinite(epoch):
        raise ValueError('Planet parameters and epoch must be finite.')
    _validate_rv_orbits(rows, options.get('stellar_mass', 1.))
    if model not in ('keplerian', 'nbody'):
        raise ValueError("model must be 'keplerian' or 'nbody'.")
    rv = None
    transit = None
    tables = {}
    astro = None
    for index, name in enumerate(('RV', 'Transit', 'TTV', 'Astrometry')):
        config = settings[name]
        if not config['enabled']:
            continue
        requested, grid = _sample_times(config, shared, streams[2*index], name)
        specific = config['parameters']
        if name in ('RV', 'Transit'):
            prepared = _observable_options(options, specific, rv_only=name == 'RV')
            curves = _evaluate_model(rows, grid, epoch=epoch, model=model, **prepared)
            prediction = curves['prediction']
            values = prediction.rv if name == 'RV' else prediction.flux
            data, uncertainty, noise = _gaussian_noise(values, config['noise'], streams[2*index+1], name)
            if name == 'RV':
                components = {planet: prediction.rv_components[:, i] for i, planet in enumerate(names)}
                rv = RVModel(grid, values, components, data, uncertainty, noise, requested)
            else:
                transit = TransitModel(grid, values, data, uncertainty, noise, requested)
        elif name == 'TTV':
            # Only prepare the orbital state at epoch; no RV/flux plot is sampled.
            specific = dict(settings['Transit']['parameters'], **specific)
            prepared = _observable_options(options, specific)
            prediction = _evaluate_model(rows, [epoch], epoch=epoch, model=model, **prepared)['prediction']
            prediction = replace(prediction, times=grid)
            tables = _transit_tables(rows, names, prediction, epoch, model, prepared, event_dataset,
                                     samples_per_orbit, max_samples, max_events)
            sequence = streams[2*index+1]
            if config['noise']['seed'] is not None:
                sequence = _seed_sequence(config['noise']['seed'], 'TTV.noise.seed')
            planet_streams = sequence.spawn(len(names))
            for slot, planet in enumerate(names):
                event = tables[planet]
                per_planet = dict(config['noise'], seed=None)
                sigma = per_planet['sigma']
                if isinstance(sigma, Mapping):
                    unknown = sigma.keys()-set(names)
                    if unknown:
                        raise ValueError(f'TTV.noise.sigma has unknown planets: {", ".join(sorted(unknown))}.')
                    per_planet['sigma'] = sigma.get(planet, 0.)
                data, uncertainty, noise = _gaussian_noise(event.BJD, per_planet, planet_streams[slot], 'TTV')
                tables[planet] = replace(event, data=data, uncertainty=uncertainty, noise=noise)
        else:
            supplied = dict(astrometry or {})
            supplied.update(specific)
            mean_anomaly = np.mod(rows[:, 4], 360.)
            if options.get('phase') == 't0':
                prepared = _observable_options(options, {}, rv_only=True)
                prediction = _evaluate_model(rows, [epoch], epoch=epoch, model=model, **prepared)['prediction']
                mean_anomaly = prediction.mean_anomaly
            astro = _astrometry_curve(rows, names, grid, epoch, model, options, supplied, mean_anomaly)
            data, uncertainty, noise = _gaussian_noise(astro.along_scan, config['noise'],
                                                       streams[2*index+1], name)
            astro = replace(astro, data=data, uncertainty=uncertainty, noise=noise, requested_BJD=requested)
    return ModelResult(rv, transit, tables, astro, settings)
