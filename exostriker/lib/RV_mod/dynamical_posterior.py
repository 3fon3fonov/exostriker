"""Stability-conditioned posteriors, independent of Qt and command-line code.

Provided by Trifon Trifonov with the aid of CODEX; use with caution.
"""
from dataclasses import dataclass, asdict
from pathlib import Path
import csv
import json
import numpy as np
from .jacobi import mass_a_from_jacobi_el

SWIFT_INTEGRATORS = ('symba', 'mvs', 'mvs_gr')
REBOUND_INTEGRATORS = ('whfast', 'mercurius', 'saba', 'ias15', 'trace', 'janus')
INTEGRATORS = SWIFT_INTEGRATORS + REBOUND_INTEGRATORS


@dataclass(frozen=True)
class IntegrationOptions:
    years: float = 1000.0
    timestep_days: float = 1.0
    output_years: float = 1.0
    timeout_seconds: float = 60.0
    integrator: str = 'symba'
    max_fractional_a_change: float = 0.1
    max_eccentricity: float = 0.95
    derive_e: bool = False
    derive_omega: bool = False
    derive_M0: bool = False
    extended_diagnostics: bool = False
    nv2016: bool = False
    nv2016_points: int = 256
    pair: tuple = ()  # 1-based GUI planet slots, inner then outer
    mmr: tuple = ()   # outer/inner period ratio p:q, p > q

    def __post_init__(self):
        if self.nv2016 and (len(self.pair) != 2 or len(self.mmr) != 2 or self.mmr[0]-self.mmr[1] != 1):
            raise ValueError('NV2016 requires --pair INNER OUTER and first-order --mmr k k-1')
        if isinstance(self.nv2016_points, bool) or int(self.nv2016_points) != self.nv2016_points or self.nv2016_points < 2:
            raise ValueError('nv2016_points must be an integer >= 2')
        if self.pair and (len(self.pair) != 2 or any(int(v) != v or not 1 <= v <= 9 for v in self.pair) or self.pair[0] == self.pair[1]):
            raise ValueError('pair must contain two distinct 1-based planet slots')
        if self.mmr and (not self.pair or len(self.mmr) != 2 or any(int(v) != v or v < 1 for v in self.mmr) or self.mmr[0] <= self.mmr[1]):
            raise ValueError('MMR requires a pair and integers p > q > 0')
        for key in ('years', 'timestep_days', 'output_years', 'timeout_seconds',
                    'max_fractional_a_change', 'max_eccentricity'):
            if not np.isfinite(getattr(self, key)) or getattr(self, key) <= 0:
                raise ValueError(f'{key} must be finite and positive')
        if self.max_eccentricity >= 1 or self.integrator not in INTEGRATORS:
            raise ValueError('Require eccentricity < 1 and a supported SWIFT or REBOUND integrator')


def extract_posterior(fit, source='auto'):
    """Return samples, normalized weights, labels. Prefer preserved raw NS rows.

    Never combine raw NS weights/logL with equal-weight resampled rows.
    """
    if source == 'auto':
        available = [s for s in ('ns', 'mcmc')
                     if np.size(getattr(getattr(fit, s + '_sampler', None), 'samples', []))]
        if len(available) != 1:
            raise ValueError('Select source="ns" or "mcmc" explicitly when both/neither exist')
        source = available[0]
    if source not in ('ns', 'mcmc'):
        raise ValueError('source must be auto, ns or mcmc')
    sampler = getattr(fit, source + '_sampler')
    labels = list(sampler.lbf['OrigLabels'])
    raw = source == 'ns' and hasattr(sampler, 'samples_ns')
    samples = np.asarray(sampler.samples_ns if raw else sampler.samples, dtype=float)
    if samples.ndim != 2 or not len(samples) or samples.shape[1] != len(labels):
        raise ValueError('Missing posterior or incompatible original labels')
    weights = np.ones(len(samples))
    if source == 'ns' and (raw or not getattr(fit, 'nest_weighted', False)):
        logwt = np.asarray(sampler.logwt, dtype=float)
        if logwt.shape != (len(samples),):
            raise ValueError('NS weights do not match raw sample rows')
        weights = np.exp(logwt - np.max(logwt))
    if not np.isfinite(weights).all() or weights.sum() <= 0:
        raise ValueError('Invalid posterior weights')
    return samples, weights / weights.sum(), labels


def orbital_draw(fit, sample, labels, *, validate_mapping_only=False):
    """Apply free orbital columns; retain fixed elements and ignore nuisance terms.

    Uses saved packing flags, not hard-coded sample columns. Refuses stale
    mappings. The session must retain the planet configuration used for sampling.
    Returns period-sorted Jacobi elements and their original GUI slots.
    """
    flags = np.asarray(fit.f_for_mcmc, dtype=int)
    if len(flags) != len(sample) or list(fit.e_for_mcmc) != list(labels):
        raise ValueError('Sampling metadata changed since sampling; use the original session')
    # Older sessions predate explicit sampling slots. Index the stored mask,
    # which may be a list or a dictionary; iterating a dict yields its keys.
    active_slots = tuple(i for i in range(len(fit.use_planet)) if fit.use_planet[i])
    slots = tuple(getattr(fit, 'mcmc_planet_slots', active_slots))
    ndset = len(getattr(fit, 'mcmc_dataset_slots', range(fit.ndset)))
    if set(slots) != set(active_slots) or not slots:
        raise ValueError('Active planets changed since sampling')
    roots = ['K', 'P', 'e_sinw', 'e_cosw', 'lamb', 'i', 'Node'] if fit.hkl else ['K', 'P', 'e', 'w', 'M0', 'i', 'Node']
    elements = np.array([[getattr(fit, r)[s] for r in roots] for s in slots], dtype=float)
    stellar_mass = float(fit.stellar_mass)
    for col, index in enumerate(flags):
        orbit_index = index - 2 * ndset
        if 0 <= orbit_index < 7 * len(slots):
            planet, param = divmod(orbit_index, 7)
            expected = getattr(fit, roots[param] + '_str')[slots[planet]]
            if labels[col] != expected:
                raise ValueError('Orbital labels disagree with saved parameter packing')
            elements[planet, param] = sample[col]
        elif labels[col] == getattr(fit, 'stellar_mass_str', 'Stellar mass'):
            stellar_mass = float(sample[col])
    if validate_mapping_only:
        return
    if fit.hkl:
        h, k, longitude = elements[:, 2:5].T.copy()
        elements[:, 2] = np.hypot(h, k)
        elements[:, 3] = np.where(elements[:, 2] > 0, np.degrees(np.arctan2(h, k)) % 360, 0.)
        elements[:, 4] = (longitude - elements[:, 3]) % 360
    order = np.argsort(elements[:, 1], kind='stable')
    elements = elements[order]
    if getattr(fit, 'copl_incl', False):
        elements[:, 5] = elements[0, 5]
    masses, axes = mass_a_from_jacobi_el(elements, len(elements), stellar_mass)
    planets = np.column_stack((masses / 1047.5654817267318, axes,
                               elements[:, 2], elements[:, 5], elements[:, 3],
                               elements[:, 6], elements[:, 4]))
    return stellar_mass, planets, np.asarray(slots)[order]


def classify_history(result, planets, options):
    """Classify recorded histories; incomplete/error runs never count as stable."""
    summaries, failure_time, reason = [], options.years, ''
    complete = result.get('completed', True)
    for initial, orbit in zip(planets, result['orbits']):
        orbit = np.asarray(orbit)
        if not len(orbit):
            complete = False
            summaries.append({})
            continue
        t, a, e = orbit[:, 0] / 365.25, orbit[:, 2], orbit[:, 3]
        bad = (~np.isfinite(orbit).all(axis=1) | (a <= 0) | (e < 0)
               | (e > options.max_eccentricity)
               | (np.abs(a / initial[1] - 1) > options.max_fractional_a_change))
        if bad.any() and (not reason or t[np.flatnonzero(bad)[0]] <= failure_time):
            failure_time = float(t[np.flatnonzero(bad)[0]])
            reason = 'orbital_threshold'
        complete &= bool(t[-1] >= options.years * (1 - 1e-8))
        summaries.append(dict(a_initial=float(initial[1]), e_initial=float(initial[2]),
                              a_min=float(np.nanmin(a)), a_max=float(np.nanmax(a)),
                              e_min=float(np.nanmin(e)), e_max=float(np.nanmax(e)),
                              e_mean=float(np.nanmean(e)), end_year=float(t[-1])))
    complete &= len(result['orbits']) == len(planets)
    status = 'unstable' if reason else ('stable' if complete else 'incomplete')
    return dict(status=status, failure_year=failure_time if reason else None,
                planets=summaries)



def angular_diagnostics(degrees):
    """Circular mean, shortest covering-arc half-width, concentration and excursion.

    Half-width is descriptive, not a test proving libration. An unresolved
    circulating angle can alias into a narrow arc at a coarse output cadence.
    """
    x = np.asarray(degrees, dtype=float) % 360
    x = x[np.isfinite(x)]
    if not len(x):
        return dict(mean=np.nan, semiamp=np.nan, resultant=np.nan, excursion=np.nan)
    vector = np.mean(np.exp(1j*np.radians(x)))
    ordered = np.sort(x)
    gaps = np.diff(np.r_[ordered, ordered[0]+360])
    excursion = np.ptp(np.degrees(np.unwrap(np.radians(x))))
    return dict(mean=float(np.degrees(np.angle(vector)) % 360) if abs(vector)>1e-12 else np.nan,
                semiamp=float((360-np.max(gaps))/2), resultant=float(abs(vector)),
                excursion=float(excursion))


def _valid_histories(result, planets, failure_year=None):
    """Common valid samples before failure; shared by all diagnostics."""
    orbits = [np.asarray(o) for o in result['orbits']]
    if len(orbits) != len(planets) or any(not len(o) for o in orbits):
        return np.array([]), []
    times = orbits[0][:, 0]
    for orbit in orbits[1:]:
        times = np.intersect1d(times, orbit[:, 0])
    if failure_year is not None:
        times = times[times < failure_year*365.25]
    histories = [o[np.searchsorted(o[:, 0], times)] for o in orbits]
    if not len(times):
        return np.array([]), []
    good = np.ones(len(times), dtype=bool)
    for o in histories:
        good &= np.isfinite(o).all(axis=1) & (o[:, 2] > 0) & (o[:, 3] >= 0) & (o[:, 3] < 1)
    # Stop at the first invalid state, rather than joining disjoint intervals.
    end = int(np.flatnonzero(~good)[0]) if (~good).any() else len(times)
    histories = [o[:end] for o in histories]
    return times[:end], histories


def describe_dynamics(result, planets, slots, stellar_mass, options, failure_year=None):
    """Describe the common valid pre-failure history; no hard-coded planet count."""
    from itertools import combinations
    from .jacobi import AU, GMSUN
    times, histories = _valid_histories(result, planets, failure_year)
    end = len(times)
    if not end:
        return {}
    values = {'diagnostic_points': end, 'diagnostic_years': float((times[end-1]-times[0])/365.25)}
    def linear(prefix, x):
        values[prefix+'_mean'] = float(np.mean(x))
        values[prefix+'_semiamp'] = float(np.ptp(x)/2)
    def angular(prefix, x):
        for key, value in angular_diagnostics(x).items():
            values[prefix+'_'+key] = value
    periods, longitudes, varpis = [], [], []
    enclosed = float(stellar_mass)
    for orbit, initial, slot in zip(histories, planets, slots):
        enclosed += initial[0]
        period = 2*np.pi*np.sqrt((orbit[:, 2]*AU)**3/(GMSUN*enclosed))/86400
        periods.append(period)
        varpi = orbit[:, 6]+orbit[:, 5]
        longitude = orbit[:, 7]+varpi
        varpis.append(varpi)
        longitudes.append(longitude)
        prefix = f'pl{slot+1}'
        for name, x in [('a', orbit[:, 2]), ('e', orbit[:, 3]), ('i', orbit[:, 4]), ('P', period)]:
            linear(prefix+'_'+name, x)
        for name, x in [('omega', orbit[:, 6]), ('Omega', orbit[:, 5]),
                        ('M', orbit[:, 7]), ('varpi', varpi), ('lambda', longitude)]:
            angular(prefix+'_'+name, x)
    pairs = list(combinations(range(len(slots)), 2))
    for inner, outer in pairs:
        prefix = f'pl{slots[inner]+1}_pl{slots[outer]+1}'
        linear(prefix+'_period_ratio', periods[outer]/periods[inner])
        ii, io = np.radians(histories[inner][:, 4]), np.radians(histories[outer][:, 4])
        node_delta = np.radians(histories[outer][:, 5]-histories[inner][:, 5])
        cos_mutual = np.cos(ii)*np.cos(io)+np.sin(ii)*np.sin(io)*np.cos(node_delta)
        linear(prefix+'_mutual_i', np.degrees(np.arccos(np.clip(cos_mutual, -1, 1))))
        angular(prefix+'_delta_omega', histories[outer][:, 6]-histories[inner][:, 6])
        angular(prefix+'_delta_varpi', varpis[outer]-varpis[inner])
    if options.mmr:
        inner, outer = [list(slots).index(v-1) for v in options.pair]
        p, q = options.mmr
        prefix = f'pl{options.pair[0]}_pl{options.pair[1]}_MMR{p}_{q}'
        for k in range(p-q+1):
            theta = q*longitudes[inner]-p*longitudes[outer]+(p-q-k)*varpis[inner]+k*varpis[outer]
            angular(prefix+f'_theta{k}', theta)
    if not options.extended_diagnostics:
        # Match make_stab.py's compact result list. Keep only the selected
        # pair (default first two in period order), not every possible pair.
        pair = options.pair or tuple(int(s)+1 for s in slots[:2])
        pair_prefix = f'pl{pair[0]}_pl{pair[1]}' if len(pair) == 2 else ''
        allowed = {f'pl{slot+1}_e_{stat}' for slot in slots for stat in ('mean', 'semiamp')}
        if pair_prefix:
            allowed.add(pair_prefix+'_period_ratio_mean')
            allowed.update(pair_prefix+'_delta_omega_'+stat for stat in ('mean', 'semiamp'))
        if options.mmr:
            p, q = options.mmr
            allowed.update(pair_prefix+f'_MMR{p}_{q}_theta{k}_'+stat
                           for k in range(p-q+1) for stat in ('mean', 'semiamp'))
        values = {key: value for key, value in values.items() if key in allowed}
    return values

def nv2016_diagnostics(result, planets, slots, stellar_mass, options, failure_year=None):
    from .resonance import normal_form, poincare_elements, DELTA_CRITICAL, NV2016_REFERENCE
    times, histories = _valid_histories(result, planets, failure_year)
    if not len(times):
        raise ValueError('No common valid orbital samples for NV2016 mapping')
    pair = [list(slots).index(v-1) for v in options.pair]
    elements, mutual = poincare_elements(histories, np.asarray(planets)[:, 0], stellar_mass, pair)
    mapped = normal_form(np.asarray(planets)[pair, 0]/stellar_mass, *elements, k=options.mmr[0])
    summary = {name+'_initial': float(mapped[name][0]) for name in
               ('Psi', 'Psi_cos_psi', 'Psi_sin_psi')}
    summary['psi_initial_deg'] = float(np.degrees(mapped['psi'][0]) % 360)
    summary['delta_initial'] = mapped['delta']
    trapezoid = getattr(np, 'trapezoid', None) or np.trapz
    for name in ('delta_series', 'Psi', 'Psi_cos_psi', 'Psi_sin_psi'):
        values = mapped[name]
        label = 'delta' if name == 'delta_series' else name
        summary[label+'_mean'] = float(trapezoid(values, times)/(times[-1]-times[0])) if len(times) > 1 else float(values[0])
        summary[label+'_semiamp'] = float(np.ptp(values)/2)
    summary['hamiltonian_semiamp'] = float(np.ptp(mapped['hamiltonian'])/2)
    summary['Phi2_scaled_semiamp'] = float(np.ptp(mapped['Phi2_scaled'])/2)
    summary['mutual_i_max_deg'] = float(np.max(mutual))
    summary['eccentricity_max'] = float(np.max(elements[1]))
    warnings = []
    if summary['eccentricity_max'] > .1:
        warnings.append('Low-eccentricity expansion: e > 0.1; higher-order terms can be important.')
    if summary['mutual_i_max_deg'] > 5:
        warnings.append('Planar approximation: mutual inclination exceeds 5 degrees.')
    if abs(elements[0][0, 0]/elements[0][0, 1]/mapped['alpha_res']-1) > .1:
        warnings.append('Pair is far from nominal commensurability; near-resonance expansion may not apply.')
    if len(planets) > 2:
        warnings.append('Isolated-pair normal form: perturbations from other planets are omitted.')
    selection = np.unique(np.linspace(0, len(times)-1, min(len(times), options.nv2016_points), dtype=int))
    history = {name: mapped[name][selection].tolist() for name in
               ('Psi', 'psi', 'Psi_cos_psi', 'Psi_sin_psi', 'delta_series', 'hamiltonian')}
    history['time_years'] = (times[selection]/365.25).tolist()
    return dict(summary=summary, history=history, warnings=warnings,
                reference=NV2016_REFERENCE, pair=list(options.pair), mmr=list(options.mmr),
                f1=mapped['f1'], f2=mapped['f2'], alpha_res=mapped['alpha_res'],
                delta_critical=DELTA_CRITICAL,
                coordinates='heliocentric canonical Poincare; pair initial angular-momentum plane',
                reference_axes='first common integration sample; coefficients at nominal resonance',
                full_history_points=len(times))


def _integrate_rebound(mass, planets, options):
    """Return Jacobi histories in the same days/degrees format as SWIFT."""
    from time import monotonic
    try:
        import rebound
    except ImportError as exc:
        raise RuntimeError('REBOUND is not installed; install it with: pip install rebound') from exc

    deadline = monotonic() + options.timeout_seconds
    sim = rebound.Simulation()
    sim.G = 4 * np.pi**2
    sim.integrator = options.integrator
    sim.dt = options.timestep_days / 365.25
    sim.add(m=mass)
    for m, a, e, inc, omega, node, anomaly in planets:
        sim.add(m=m, a=a, e=e, inc=np.radians(inc), omega=np.radians(omega),
                Omega=np.radians(node), M=np.radians(anomaly))
    sim.move_to_com()
    timed_out = False

    def heartbeat(pointer):
        nonlocal timed_out
        if monotonic() >= deadline:
            timed_out = True
            pointer.contents.stop()

    try:
        sim.heartbeat = heartbeat
    except AttributeError as exc:
        # REBOUND 5.1.x's attribute guard rejects the callback storage name.
        if "'_hb'" not in str(exc):
            raise
        from rebound.simulation import AFF
        heartbeat_callback = AFF(heartbeat)  # Keep alive until integration ends.
        sim._heartbeat = heartbeat_callback
    histories = [[] for _ in planets]
    times = np.r_[np.arange(0., options.years, options.output_years), options.years]
    for target in times:
        if target > sim.t:
            try:
                sim.integrate(float(target), exact_finish_time=0)
            except rebound.GenericError as exc:
                raise RuntimeError(str(exc)) from exc
        if timed_out or monotonic() >= deadline:
            raise TimeoutError('REBOUND integration exceeded timeout_seconds')
        for index, orbit in enumerate(sim.orbits()):
            histories[index].append([
                sim.t * 365.25, index + 1, orbit.a, orbit.e,
                np.degrees(orbit.inc), np.degrees(orbit.Omega) % 360,
                np.degrees(orbit.omega) % 360, np.degrees(orbit.M) % 360,
                orbit.a * (1-orbit.e), orbit.a * (1+orbit.e), planets[index, 0]/mass])
    return dict(orbits=[np.asarray(h) for h in histories], completed=True)


def _integrate_orbits(mass, planets, slots, options, gr_step):
    """Process worker: only small initial-condition arrays cross the boundary."""
    try:
        if options.integrator in REBOUND_INTEGRATORS:
            result = _integrate_rebound(mass, planets, options)
        else:
            from ..swift import integrate
            result = integrate(mass, planets, timemax=options.years,
                           timestep=options.timestep_days, output_step=options.output_years,
                           timeout_sec=options.timeout_seconds, integrator=options.integrator,
                           gr_step=gr_step)
        record = classify_history(result, planets, options)
        record['planet_slots'] = slots.tolist()
        record['dynamics'] = describe_dynamics(result, planets, slots, mass, options, record['failure_year'])
        if options.nv2016:
            try:
                resonance = nv2016_diagnostics(result, planets, slots, mass, options, record['failure_year'])
                record['resonance_nv2016'] = resonance
                prefix = f'pl{options.pair[0]}_pl{options.pair[1]}_NV2016_'
                record['dynamics'].update({prefix+key: value for key, value in resonance['summary'].items()
                                          if key.endswith(('_initial', '_initial_deg', '_mean'))})
            except ValueError as exc:
                # Normal-form applicability does not change the N-body stability verdict.
                record['resonance_nv2016'] = {'error': str(exc)}
        return record
    except TimeoutError as exc:
        return dict(status='timeout', error=str(exc))
    except (ValueError, RuntimeError) as exc:
        return dict(status='error', error=str(exc))


def integrate_draw(fit, sample, labels, options):
    """One draw; no session mutations or fitting. Suitable for a GUI worker."""
    try:
        mass, planets, slots = orbital_draw(fit, sample, labels)
    except ValueError as exc:
        return dict(status='error', error=str(exc))
    return _integrate_orbits(mass, planets, slots, options, getattr(fit, 'GR_step', 1))


def integrate_posterior(fit, options=None, source='auto', count=None, seed=42, progress=None, workers=1):
    """Integrate all rows or weighted draws (with replacement) reproducibly.

    workers selects independent spawned CPU processes (default 1).
    progress(done, total, record) runs in the caller. Returns plain arrays/dicts for
    later GUI integration. No changes are made to the input fit or its sampler.
    """
    if isinstance(workers, bool) or int(workers) != workers or workers < 1:
        raise ValueError('workers must be a positive integer')
    workers = int(workers)
    options = options or IntegrationOptions()
    samples, weights, labels = extract_posterior(fit, source)
    if options.pair and any(not fit.use_planet[v-1] for v in options.pair):
        raise ValueError('The selected pair must contain enabled planets')
    if options.pair and fit.P[options.pair[0]-1] >= fit.P[options.pair[1]-1]:
        raise ValueError('Specify the inner planet first in --pair')
    if count is None:
        indices = np.arange(len(samples))
        selected_weights = weights.copy()
    else:
        if count < 1 or int(count) != count:
            raise ValueError('count must be a positive integer or None')
        indices = np.random.default_rng(seed).choice(len(samples), int(count), p=weights)
        selected_weights = np.full(len(indices), 1 / len(indices))
    # Validate metadata once before launching potentially expensive integrations.
    orbital_draw(fit, samples[indices[0]], labels, validate_mapping_only=True)
    cache = {}
    unique, multiplicity = np.unique(indices, return_counts=True)
    counts = dict(zip(unique.tolist(), multiplicity.tolist()))
    completed = 0

    def report(index, record):
        nonlocal completed
        cache[int(index)] = record
        completed += counts[int(index)]
        if progress:
            progress(completed, len(indices), dict(record, sample_index=int(index)))

    if workers == 1:
        for index in unique:
            report(index, integrate_draw(fit, samples[index], labels, options))
    else:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        import multiprocessing
        # Spawn avoids inheriting Qt threads and SWIFT locks/native state.
        with ProcessPoolExecutor(max_workers=min(workers, len(unique)),
                                 mp_context=multiprocessing.get_context('spawn')) as pool:
            futures = {}
            for index in unique:
                try:
                    mass, planets, slots = orbital_draw(fit, samples[index], labels)
                except ValueError as exc:
                    report(index, dict(status='error', error=str(exc)))
                    continue
                future = pool.submit(_integrate_orbits, mass, planets, slots,
                                     options, getattr(fit, 'GR_step', 1))
                futures[future] = int(index)
            for future in as_completed(futures):
                # Worker crashes propagate; never silently classify them as unstable.
                report(futures[future], future.result())
    records = [dict(cache[int(index)], sample_index=int(index)) for index in indices]
    stable = np.array([r['status'] == 'stable' for r in records])
    from hkl_plot import derived_hkl_columns
    derived = derived_hkl_columns(fit, samples[indices], labels, asdict(options))
    return dict(samples=samples[indices], weights=selected_weights, labels=labels,
                indices=indices, stable=stable, records=records, options=asdict(options),
                source=source, seed=seed, workers=workers,
                derived_labels=[column[2] for column in derived],
                derived_samples=np.column_stack([column[3] for column in derived]) if derived else np.empty((len(indices), 0)))



def map_nv2016_posterior(fit, options, source='auto', count=None, seed=42):
    """Epoch coordinates only, without rerunning integrations.

    Integrated sessions retain exact original samples/weights/statuses in
    dynamical_posterior; use those instead of the GUI's augmented sampler.
    """
    if not options.nv2016:
        raise ValueError('Enable nv2016 when mapping resonance coordinates')
    saved = getattr(fit, 'dynamical_posterior', None)
    if saved:
        if source != 'auto' and saved['source'] not in ('auto', source):
            raise ValueError('Requested source disagrees with the integrated posterior')
        samples, weights, labels = saved['samples'], saved['weights'], saved['labels']
        sample_indices = saved['indices']
        statuses = [record['status'] for record in saved['records']]
        source = saved['source']
    else:
        samples, weights, labels = extract_posterior(fit, source)
        sample_indices = np.arange(len(samples))
        statuses = ['not_integrated']*len(samples)
    if count is None:
        selection = np.arange(len(samples))
        selected_weights = np.array(weights, copy=True)
    else:
        if isinstance(count, bool) or int(count) != count or count < 1:
            raise ValueError('count must be a positive integer or None')
        selection = np.random.default_rng(seed).choice(len(samples), int(count), p=weights)
        selected_weights = np.full(len(selection), 1/len(selection))
    orbital_draw(fit, samples[selection[0]], labels, validate_mapping_only=True)
    cache = {}
    for index in np.unique(selection):
        try:
            mass, planets, slots = orbital_draw(fit, samples[index], labels)
            if any(v-1 not in slots for v in options.pair):
                raise ValueError('Selected pair must contain active planets')
            if list(slots).index(options.pair[0]-1) >= list(slots).index(options.pair[1]-1):
                raise ValueError('Specify the inner planet first in --pair')
            histories = []
            for planet in planets:
                history = np.zeros((1, 11))
                history[0, 2:8] = planet[[1, 2, 3, 5, 4, 6]]
                histories.append(history)
            mapped = nv2016_diagnostics({'orbits': histories}, planets, slots, mass, options)
            mapped['summary'] = {key:value for key,value in mapped['summary'].items()
                                 if not key.endswith(('_mean', '_semiamp'))}
            mapped['warnings'].append('Epoch-only mapping: no new integration or time averages.')
        except ValueError as exc:
            mapped = {'error': str(exc)}
        cache[int(index)] = mapped
    records = [dict(status=statuses[index], sample_index=int(sample_indices[index]),
                    resonance_nv2016=cache[int(index)]) for index in selection]
    return dict(records=records, indices=np.asarray(sample_indices)[selection], weights=selected_weights,
                stable=np.array([record['status'] == 'stable' for record in records]),
                options=asdict(options), seed=seed, source=source)


def diagnostic_label(key):
    """Plot labels follow P$_b$, K$_b$, e$_b$ (slot 1=b, 2=c, ...)."""
    import re
    nv = re.fullmatch(r'pl(\d+)_pl(\d+)_NV2016_(.+)', key)
    if nv:
        first, second, quantity = nv.groups()
        pair = chr(96+int(first)+1)+chr(96+int(second)+1)
        names = {'delta_initial': r'$\delta_0$', 'delta_mean': r'$\langle\delta\rangle$',
                 'Psi_initial': r'$\Psi_0$', 'Psi_mean': r'$\langle\Psi\rangle$',
                 'psi_initial_deg': r'$\psi_0$ [deg]',
                 'Psi_cos_psi_initial': r'$(\Psi\cos\psi)_0$',
                 'Psi_sin_psi_initial': r'$(\Psi\sin\psi)_0$',
                 'Psi_cos_psi_mean': r'$\langle\Psi\cos\psi\rangle$',
                 'Psi_sin_psi_mean': r'$\langle\Psi\sin\psi\rangle$'}
        return 'NV2016 '+pair+' '+names.get(quantity, quantity)
    match = re.fullmatch(r'pl(\d+)(?:_pl(\d+))?_(.+)_(mean|semiamp|resultant|excursion)', key)
    if not match:
        return {'diagnostic_points': 'Recorded points',
                'diagnostic_years': 'Recorded duration [yr]'}.get(key, key)
    first, second, quantity, statistic = match.groups()
    inner = chr(ord('a') + int(first))
    outer = chr(ord('a') + int(second)) if second else ''
    pair = inner + outer
    prefix = {'mean': 'mean ', 'semiamp': 'ampl. ',
              'resultant': 'concentration ', 'excursion': 'excursion '}[statistic]
    if quantity == 'period_ratio':
        return prefix + f'P$_{outer}$/P$_{inner}$'
    resonance = re.fullmatch(r'MMR(\d+)_(\d+)_theta(\d+)', quantity)
    if resonance:
        p, q, angle = resonance.groups()
        symbol = rf'\theta_{{{angle},{pair}}}'
        suffix = f' ({p}:{q})'
    else:
        symbols = {'a': 'a', 'e': 'e', 'i': 'i', 'P': 'P',
                   'omega': r'\omega', 'Omega': r'\Omega', 'M': 'M',
                   'varpi': r'\varpi', 'lambda': r'\lambda',
                   'delta_omega': r'\Delta\omega', 'delta_varpi': r'\Delta\varpi',
                   'mutual_i': r'i_{\mathrm{mut}}'}
        symbol = symbols.get(quantity, quantity)
        # Simple parameters match the existing K/P/e label style exactly.
        if quantity == 'mutual_i':
            body = r'$i_{\mathrm{mut},' + pair + '}$'
        elif quantity in ('a', 'e', 'i', 'P', 'M'):
            body = f'{symbol}$_{pair}$'
        else:
            body = '$' + symbol + '_{' + pair + '}$'
        units = ' [AU]' if quantity == 'a' else ' [d]' if quantity == 'P' else ''
        if quantity not in ('a', 'e', 'P'):
            units = ' [deg]'
        return prefix + body + ('' if statistic == 'resultant' else units)
    return prefix + '$' + symbol + '$' + suffix + ('' if statistic == 'resultant' else ' [deg]')


def dynamical_columns(result):
    keys = sorted({key for r in result['records'] for key in r.get('dynamics', {})})
    labels = list(result.get('derived_labels', [])) + [diagnostic_label(key) for key in keys]
    dynamics = np.array([[r.get('dynamics', {}).get(key, np.nan) for key in keys]
                             for r in result['records']], dtype=float).reshape(len(result['records']), len(keys))
    derived = np.asarray(result.get('derived_samples', np.empty((len(result['records']), 0))))
    return labels, np.column_stack((derived, dynamics))

def export_posterior(result, directory):
    """Write aligned posterior rows, stable subset, diagnostics, and interval table."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    samples, weights, stable = (result[k] for k in ('samples', 'weights', 'stable'))
    diagnostic_labels, diagnostics = dynamical_columns(result)
    samples = np.column_stack((samples, diagnostics))
    labels = list(result['labels']) + diagnostic_labels
    for name, mask in [('posterior', np.ones(len(samples), dtype=bool)), ('stable_posterior', stable)]:
        with (directory / (name + '.csv')).open('w', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(['sample_index', 'weight', 'status'] + labels)
            norm = weights[mask].sum()
            for j in np.flatnonzero(mask):
                writer.writerow([result['indices'][j], weights[j] / norm,
                                 result['records'][j]['status']] + samples[j].tolist())
    np.savez_compressed(directory / 'posterior.npz', samples=samples, weights=weights,
                        stable=stable, indices=result['indices'], labels=labels)
    metadata = {k: result[k] for k in ('records', 'options', 'source', 'seed', 'workers')}
    metadata['stable_weight_fraction'] = float(weights[stable].sum())
    metadata['credit'] = 'Provided by Trifon Trifonov with the aid of CODEX; use with caution.'
    # Full mapped trajectories are stored once per source sample in resonance.npz.
    metadata['records'] = [dict(record, resonance_nv2016={k:v for k,v in record['resonance_nv2016'].items() if k != 'history'})
                           if 'resonance_nv2016' in record else record for record in result['records']]
    (directory / 'diagnostics.json').write_text(json.dumps(metadata, indent=2))
    with (directory / 'summary.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['population', 'parameter', 'q15.85', 'median', 'q84.15'])
        for name, mask in [('original', np.ones(len(samples), dtype=bool)), ('stable', stable)]:
            if not mask.any():
                continue
            for col, label in enumerate(labels):
                x, w = samples[mask, col], weights[mask]
                order = np.argsort(x)
                q = np.interp([.1585, .5, .8415], (np.cumsum(w[order]) - .5*w[order])/w.sum(), x[order])
                writer.writerow([name, label] + q.tolist())
    if result['options'].get('nv2016', False):
        from .resonance import export_resonance
        export_resonance(result, directory)
    return directory


def save_integrated_session(fit, result, filename):
    """Save a GUI-readable session with integration columns, like make_stab.py.

    GUI samplers are equal-weight, so weighted rows are systematically resampled
    with all derived columns and likelihoods using the SAME indices. The exact
    tested rows/weights/statuses remain in saved.dynamical_posterior.
    """
    import dill
    filename = Path(filename)
    if filename.exists():
        raise FileExistsError(filename)
    saved = dill.copy(fit)
    saved.dynamical_posterior = dill.copy(result)
    source = result['source']
    if source == 'auto':
        available = [s for s in ('ns', 'mcmc')
                     if np.size(getattr(getattr(fit, s + '_sampler', None), 'samples', []))]
        if len(available) != 1:
            raise ValueError('Cannot determine sampled source')
        source = available[0]
    sampler = getattr(saved, source + '_sampler')
    original = getattr(fit, source + '_sampler')
    n = len(result['samples'])
    weights = np.asarray(result['weights'], dtype=float)
    if np.allclose(weights, np.full(n, 1/n), rtol=1e-12, atol=0):
        selection = np.arange(n)
    else:
        cumulative = np.cumsum(weights / weights.sum())
        cumulative[-1] = 1
        positions = (np.arange(n) + np.random.default_rng(result['seed']).random()) / n
        selection = np.searchsorted(cumulative, positions)
    # Align likelihoods with the raw/equal-weight population used by extraction.
    raw = source == 'ns' and hasattr(original, 'samples_ns')
    likelihoods = np.asarray(getattr(original, 'logl_ns' if raw else
        ('logl' if source == 'ns' else 'lnL'), []), dtype=float).reshape(-1)
    population = np.asarray(original.samples_ns if raw else original.samples)
    aligned = likelihoods.shape == (len(population),)
    # Old equal-weight NS sessions may have retained the raw logL vector.
    if source == 'ns' and not raw and getattr(fit, 'nest_weighted', False):
        aligned = False
    logl = likelihoods[result['indices']] if aligned else np.full(n, np.nan)
    derived_labels, diagnostics = dynamical_columns(result)
    columns = list(diagnostics.T)
    # Unknown/failed runs get zero (never accepted by the legacy time overlay).
    lifetime = [result['options']['years'] if r['status'] == 'stable' else
                max(0., float(np.nextafter(r.get('failure_year') or 0., 0.))) if r['status'] == 'unstable' else 0.
                for r in result['records']]
    columns.append(lifetime)
    labels = list(result['labels']) + derived_labels + ['max. time']
    samples = np.column_stack([result['samples'], np.asarray(columns).T])[selection]
    sampler.samples = samples
    setattr(sampler, 'logl' if source == 'ns' else 'lnL', logl[selection])
    sampler.lbf['OrigLabels'] = np.asarray(labels)
    for j, label in enumerate(labels):
        sampler.lbf[j] = np.array([label, True])
    sampler.lbf['cornerplot']['stab_threshold'] = float(np.nextafter(result['options']['years'], 0))
    # Full tested population remains visible; max. time enables the stable overlay.
    stats = getattr(saved, 'nest_stat' if source == 'ns' else 'mcmc_stat')
    for name, func in [('mean', np.nanmean), ('median', np.nanmedian)]:
        finite_columns = np.isfinite(samples).any(axis=0)
        stats[name] = np.full(samples.shape[1], np.nan)
        stats[name][finite_columns] = func(samples[:, finite_columns], axis=0)
    best = int(np.nanargmax(logl[selection])) if np.isfinite(logl[selection]).any() else 0
    stats['best'] = samples[best].copy()
    # Use median for the optional mode marker rather than stale-width statistics.
    stats['mode'] = stats['median'].copy()
    sampler.lbf['mode'] = False
    sampler.lbf['best_gui'] = False
    if not aligned:
        sampler.lbf['best_samp'] = False
        sampler.lbf['median'] = True
    if source == 'ns':
        # Keep raw evidence arrays intact in the original sampler backup, not
        # falsely associated with the augmented posterior's new columns.
        saved.dynamical_posterior_original_sampler = dill.copy(original)
        sampler.samples_ns = samples.copy()
        sampler.logl_ns = logl[selection].copy()
        sampler.logwt = np.full(n, -np.log(n))
        sampler.weights = np.full(n, 1/n)
        saved.nest_weighted = True
    saved.dynamical_posterior_session = True
    saved.dynamical_posterior_credit = 'Provided by Trifon Trifonov with the aid of CODEX; use with caution.'
    with filename.open('xb') as stream:
        dill.dump(saved, stream)
    return filename
