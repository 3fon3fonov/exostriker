"""Read-only session adapters and pointwise posterior bands for GetModel.

Instrument offsets, jitters and GP parameters are nuisance parameters here;
models describe physical orbital signals. Weighted NS rows retain their weights.
"""
from copy import deepcopy
from dataclasses import dataclass, replace
from math import erf, sqrt

import numpy as np

from model_only import GetModel, _model_settings, _standard_parameters
from model_epoch import resolve_epoch
from native_transit_model import LD_COUNTS, LD_LAWS, session_jacobi_convention


_GP_COUNTS = {'RotKernel': 4, 'SHOKernel': 3, 'Matern32': 3, 'RealTerm': 2, 'dSHOKernel': 6}


class SessionModelAdapter:
    """Reuse one session's parameter mapping without modifying its fit or sampler.

    sample, when supplied, is one row in the saved f_for_mcmc column order.
    Sparse planet slots retain names planet_1 ... planet_9. Both MCMC and NS
    use this order. Repacking or sorting a fit after sampling invalidates that
    correspondence; preserve the original sampling session.
    """

    def __init__(self, session, *, radius_mode=None, phase=None):
        self.session = session
        self.slots = tuple(int(slot) for slot in np.flatnonzero(session.use_planet))
        self.datasets = tuple(i for i, row in session.rv_data_sets.items() if len(row))
        self.transits = tuple(i for i, row in session.tra_data_sets.items() if len(row))
        self.mode = 'nbody' if session.mod_dynamical else 'keplerian'
        self.epoch = resolve_epoch(session)
        self.phase = phase or ('t0' if self.mode == 'keplerian' and session.type_fit.get('Transit', False) else 'M0')
        self.radius_mode = radius_mode or ('density' if session.rho else ('scaled' if self.transits and self.mode == 'keplerian' else 'stellar'))
        if session.gr_flag or any(session.omega_dot[i] != 0 for i in self.slots):
            raise ValueError('GetModel session adapters currently support Newtonian orbits without apsidal precession.')
        if bool(np.asarray(session.mixed_fit[0]).item()):
            raise ValueError('Mixed dynamical subsets cannot be represented by one GetModel dynamics setting.')
        if session.get_TTVs[0]:
            raise ValueError('Independent fitted per-event time shifts require a custom posterior evaluator.')
        self.orbits = np.array([[getattr(session, name)[slot] for name in
                                ('mp' if getattr(session, 'mass_parameterization', False) else 'K', 'P', 'e', 'w', 'M0', 'i', 'Node')]
                               for slot in self.slots], float).reshape(-1, 7)
        if session.hkl:
            for j, slot in enumerate(self.slots):
                self.orbits[j, 2:5] = [session.e_sinw[slot], session.e_cosw[slot], session.lamb[slot]]
        self.geometry = np.array([[session.t0[slot], session.pl_a[slot], session.pl_rad[slot]]
                                  for slot in self.slots], float).reshape(-1, 3)
        self.mass = float(session.params.stellar_mass)
        self.radius = float(session.stellar_radius)
        self.trends = np.array([session.rv_lintr, session.rv_quadtr], float)
        self.laws = list(session.ld_m)
        self.coefficients = np.zeros((60, 4))
        tables = {'linear': 'ld_u_lin', 'quadratic': 'ld_u_quad', 'nonlinear': 'ld_u_nonlin'}
        for slot, law in enumerate(self.laws):
            if law == 'uniform':
                continue
            if law not in tables:
                raise ValueError(f'Unsupported session limb profile: {law}.')
            count = LD_COUNTS[LD_LAWS[law]]
            self.coefficients[slot, :count] = getattr(session, tables[law])[slot]
        self.parameter_map = {}
        self.ignored_indices = set()
        self._layout()

    def _layout(self):
        """Map the saved packed vector, including GP and limb-group block sizes."""
        session = self.session
        nd = len(self.datasets)
        npl = len(self.slots)
        self.ignored_indices.update(range(2*nd))
        for j in range(npl):
            for column in range(7):
                self.parameter_map[2*nd+7*j+column] = ('orbit', j, column)
        position = 2*nd+7*npl
        self.parameter_map[position] = ('trend', 0)
        self.parameter_map[position+1] = ('trend', 1)
        position += 2
        gp = _GP_COUNTS[session.gp_kernel] if session.rtg[1] else 0
        self.ignored_indices.update(range(position, position+gp))
        position += gp
        for j in range(npl):
            for column in range(3):
                self.parameter_map[position+3*j+column] = ('geometry', j, column)
        position += 3*npl
        count = len(self.transits)
        self.ignored_indices.update(range(position, position+2*count))
        position += 2*count
        gp = _GP_COUNTS[session.tra_gp_kernel] if session.rtg[3] else 0
        self.ignored_indices.update(range(position, position+gp))
        position += gp
        # Non-orbital transit trends and dilution are not silently mapped to orbits.
        position += 2*count+npl
        for slot in self.transits:
            if session.ld_gr[slot] != slot:
                continue
            size = LD_COUNTS[LD_LAWS[self.laws[slot]]]
            for column in range(size):
                self.parameter_map[position] = ('limb', slot, column)
                position += 1
        packed = np.asarray(getattr(session, 'parameters', []), float)
        if packed.size:
            for column in range(5):
                self.parameter_map[len(packed)-6+column] = ('astrometry', column)
            self.parameter_map[len(packed)-1] = ('mass',)

    def inputs(self, sample=None):
        """Return named planets, GetModel keywords and RV/sky offset terms."""
        session = self.session
        orbits = self.orbits.copy()
        geometry = self.geometry.copy()
        coefficients = self.coefficients.copy()
        trends = self.trends.copy()
        mass = self.mass
        correction = np.zeros(5)
        if sample is not None:
            stored_slots = tuple(getattr(session, 'mcmc_planet_slots', self.slots))
            stored_datasets = tuple(getattr(session, 'mcmc_dataset_slots', self.datasets))
            if stored_slots != self.slots or stored_datasets != self.datasets:
                raise ValueError('Saved posterior slot order differs from the session. Use its original sampling session.')
            flags = np.asarray(session.f_for_mcmc)
            values = np.asarray(sample, float)
            if flags.ndim != 1 or values.shape != flags.shape or not np.all(np.isfinite(values)):
                raise ValueError('A posterior row must match the saved f_for_mcmc column order.')
            if np.any(flags != flags.astype(int)) or len(set(flags)) != len(flags):
                raise ValueError('Invalid saved free-parameter indices.')
            packed = np.asarray(session.parameters, float)
            if packed.ndim != 1 or not len(packed) or np.any(flags < 0) or np.any(flags >= len(packed)):
                raise ValueError('The sampling session must retain parameters and f_for_mcmc.')
            for index, value in zip(flags.astype(int), values):
                target = self.parameter_map.get(index)
                if target is None:
                    if index in self.ignored_indices:
                        continue
                    raise ValueError(f'Posterior parameter at packed index {index} is unsupported by the orbital adapter; use a custom evaluator.')
                if target[0] == 'orbit':
                    orbits[target[1], target[2]] = value
                elif target[0] == 'geometry':
                    geometry[target[1], target[2]] = value
                elif target[0] == 'trend':
                    trends[target[1]] = value
                elif target[0] == 'mass':
                    mass = value
                elif target[0] == 'astrometry':
                    field = ('ast_alpha', 'ast_delta', 'ast_pi', 'ast_mu_alpha', 'ast_mu_delta')[target[1]]
                    correction[target[1]] = value-getattr(session, field)[0]
                else:
                    for slot in range(60):
                        if session.ld_gr[slot] == target[1]:
                            coefficients[slot, target[2]] = value
        if session.hkl and len(orbits):
            omega = np.rad2deg(np.arctan2(orbits[:, 2], orbits[:, 3]))
            eccentricity = np.hypot(orbits[:, 2], orbits[:, 3])
            orbits[:, 4] = (orbits[:, 4]-omega) % 360.
            orbits[:, 2] = eccentricity
            orbits[:, 3] = omega
        if self.mode == 'nbody' and session.copl_incl and len(orbits):
            orbits[:, 5] = orbits[0, 5]
        if getattr(session, 'mass_parameterization', False):
            from mass_parameterization import masses_to_k
            orbits[:, 0] = masses_to_k(orbits, orbits[:, 0].copy(), mass, self.mode == 'nbody')
        planets = {}
        fields = ('K', 'P', 'e', 'omega', 'M0', 'inclination', 'node')
        for j, slot in enumerate(self.slots):
            planet = dict(zip(fields, orbits[j]))
            planet['radius_ratio'] = geometry[j, 2]
            if self.radius_mode == 'scaled':
                planet['scaled_axis'] = geometry[j, 1]
            planets[f'planet_{slot+1}'] = planet
        options = dict(epoch=self.epoch, model=self.mode, stellar_mass=mass, stellar_radius=self.radius,
                       limb_dark=self.laws.copy(), coefficients=coefficients, phase=self.phase)
        if self.radius_mode == 'density':
            options['rho_star'] = geometry[0, 1] if len(geometry) else session.stellar_rho
        elif self.radius_mode not in ('scaled', 'stellar'):
            raise ValueError('radius_mode must be density, scaled or stellar.')
        if self.mode == 'nbody':
            options['jacobi_convention'] = session_jacobi_convention(session)
        if self.phase == 't0':
            options['t0'] = geometry[:, 0]
        return planets, options, trends, correction

    def __call__(self, sample=None, **kwargs):
        """Evaluate one session/sample, including deterministic RV trends."""
        planets, options, trends, correction = self.inputs(sample)
        options.update(deepcopy(kwargs))
        selected = _model_settings(options.get('models'), options.get('transit_events'))
        if selected['Transit']['enabled'] or selected['TTV']['enabled']:
            for slot in self.transits:
                dilution = getattr(self.session, 'tra_dil', {}).get(slot, 1.)
                if dilution != 1. or self.session.tra_lintr[slot] != 0. or self.session.tra_quadtr[slot] != 0.:
                    raise ValueError('Session dilution/non-orbital transit trends require a custom model evaluator.')
        if np.any(correction):
            models = selected
            options['models'] = models
            astro = models['Astrometry']
            if astro['enabled']:
                parameters = astro.setdefault('parameters', {})
                if 'standard' not in parameters and options.get('astrometry'):
                    parameters.update(options['astrometry'])
                if 'standard' not in parameters:
                    raise ValueError('Supply the absolute best-fit astrometric standard solution for posterior propagation.')
                standard = _standard_parameters(parameters['standard']).copy()
                standard[0] += correction[0]/(np.cos(np.deg2rad(standard[1]))*3.6e6)
                standard[1] += correction[1]/3.6e6
                standard[2:] += correction[2:]
                parameters['standard'] = standard
        result = GetModel(planets, **options)
        if result.RV is not None and np.any(trends):
            elapsed = result.RV.BJD-options['epoch']
            trend = trends[0]*elapsed+trends[1]*elapsed**2
            components = dict(result.RV.components, trend=trend)
            rv = replace(result.RV, model=result.RV.model+trend, data=result.RV.data+trend, components=components)
            result = replace(result, RV=rv)
        if result.Astrometry is not None and np.any(correction[:2]):
            astro = result.Astrometry
            configuration = result.settings['Astrometry']['parameters']
            angle = np.broadcast_to(configuration.get('scan_angle', 0.), astro.BJD.shape)
            shift = correction[0]*np.sin(angle)+correction[1]*np.cos(angle)
            astro = replace(astro, East=astro.East+correction[0], North=astro.North+correction[1],
                            along_scan=astro.along_scan+shift, data=astro.data+shift)
            result = replace(result, Astrometry=astro)
        return result


def GetModelFromSession(session, sample=None, *, radius_mode=None, phase=None, **kwargs):
    """Evaluate a loaded .ses fit or one of its saved MCMC/NS parameter rows."""
    return SessionModelAdapter(session, radius_mode=radius_mode, phase=phase)(sample, **kwargs)


@dataclass(frozen=True)
class PosteriorSamples:
    values: np.ndarray
    weights: np.ndarray
    source: str


def _sampler_array(sampler, name):
    """Read stored arrays without invoking dill-restored sampler methods.

    Some legacy dynesty Results classes were pickled by value. Their Python
    bytecode is not portable between interpreter versions, but their dictionary
    entries remain usable. Bypass custom __getattr__ and __getitem__ methods.
    """
    if isinstance(sampler, dict):
        return dict.get(sampler, name, [])
    try:
        attributes = object.__getattribute__(sampler, '__dict__')
    except AttributeError:
        return []
    return attributes.get(name, [])


def PosteriorSamplesFromSession(session, source='mcmc'):
    """Read aligned retained MCMC rows or raw NS rows and normalized weights."""
    if source.lower() == 'mcmc':
        sampler = session.mcmc_sampler
        values = np.asarray(_sampler_array(sampler, 'samples'), float)
        weights = np.ones(len(values))
    elif source.lower() in ('ns', 'nested'):
        sampler = session.ns_sampler
        values = np.asarray(_sampler_array(sampler, 'samples_ns'), float)
        logwt = np.asarray(_sampler_array(sampler, 'logwt'), float)
        if values.ndim != 2 or logwt.shape != (len(values),):
            raise ValueError('NS propagation requires aligned raw samples_ns and logwt, not resampled samples.')
        if not np.all(np.isfinite(logwt) | np.isneginf(logwt)) or not np.any(np.isfinite(logwt)):
            raise ValueError('Nested-sampling log weights must contain finite posterior mass.')
        weights = np.exp(logwt-np.max(logwt))
    else:
        raise ValueError('source must be mcmc or ns.')
    if values.ndim != 2 or not len(values) or not np.all(np.isfinite(values)):
        raise ValueError(f'No finite retained {source} sample matrix was found.')
    if values.shape[1] != len(session.f_for_mcmc):
        raise ValueError('Posterior columns differ from this session free-parameter mapping.')
    return PosteriorSamples(values.copy(), weights/np.sum(weights), source)


@dataclass(frozen=True)
class CredibleBand:
    """Pointwise equal-tail credible intervals; coverage reports transit survival."""

    BJD: np.ndarray
    best: np.ndarray
    median: np.ndarray
    intervals: dict
    draws: np.ndarray
    coverage: np.ndarray


@dataclass(frozen=True)
class ModelEnsemble:
    reference: object
    RV: object
    Transit: object
    TTV: dict
    Astrometry: dict
    indices: np.ndarray
    weights: np.ndarray
    effective_sample_size: float
    rejected: tuple
    retained_mass: float


def _band(times, best, draws, weights, levels, minimum_coverage):
    """Compute weighted inverse-CDF quantiles without losing NS row weights."""
    draws = np.asarray(draws, float)
    median = np.full(len(best), np.nan)
    intervals = {level: np.full((2, len(best)), np.nan) for level in levels}
    coverage = np.zeros(len(best))
    for column in range(len(best)):
        keep = np.isfinite(draws[:, column])
        coverage[column] = np.sum(weights[keep])
        if coverage[column]+1.e-14 < minimum_coverage or not np.any(keep):
            continue
        values = draws[keep, column]
        local_weights = weights[keep]
        order = np.argsort(values, kind='stable')
        cumulative = np.cumsum(local_weights[order])/np.sum(local_weights)
        cumulative[-1] = 1.
        def quantile(probability):
            index = min(np.searchsorted(cumulative, probability), len(order)-1)
            return values[order[index]]
        median[column] = quantile(.5)
        for level in levels:
            tail = .5*(1.-erf(level/sqrt(2.)))
            intervals[level][:, column] = [quantile(tail), quantile(1.-tail)]
    return CredibleBand(np.asarray(times).copy(), np.asarray(best).copy(), median, intervals, draws, coverage)


def GetModelEnsemble(evaluator, samples, *, reference, weights=None, max_samples=None,
                     seed=42, levels=(1, 2, 3), minimum_coverage=.95, on_invalid='raise',
                     event_periods=None):
    """Propagate parameter rows through any callback returning a GetModel result.

    evaluator(row) must keep BJD grids fixed and disable simulation noise/jitter.
    reference is the best-fit result, not a noisy prediction or marginal median.
    No rows are capped unless max_samples is explicitly supplied. A cap resamples
    posterior mass with replacement and returns equal weights and original indices.
    Bands describe latent model uncertainty, not measurement-noise prediction.
    Event times are matched by orbital cycle, never by array position.
    """
    samples = np.asarray(samples, float)
    if samples.ndim != 2 or not len(samples) or not np.all(np.isfinite(samples)):
        raise ValueError('samples must be a nonempty finite two-dimensional matrix.')
    weights = np.ones(len(samples)) if weights is None else np.asarray(weights, float)
    if weights.shape != (len(samples),) or not np.all(np.isfinite(weights)) or np.any(weights < 0) or not np.any(weights > 0):
        raise ValueError('weights must be finite, nonnegative and aligned with sample rows.')
    weights = weights/np.max(weights)
    weights = weights/weights.sum()
    if not levels or any(not np.isfinite(level) or level <= 0 for level in levels):
        raise ValueError('Credible-band levels must be positive sigma equivalents.')
    if not 0 < minimum_coverage <= 1 or on_invalid not in ('raise', 'skip'):
        raise ValueError('Use minimum_coverage in (0,1] and on_invalid=raise or skip.')
    indices = np.flatnonzero(weights > 0)
    weights = weights[indices]
    if max_samples is not None:
        if isinstance(max_samples, bool) or int(max_samples) != max_samples or max_samples < 1:
            raise ValueError('max_samples must be a positive integer or None.')
        if len(indices) > max_samples:
            indices = np.random.default_rng(seed).choice(indices, size=int(max_samples), p=weights)
            weights = np.full(len(indices), 1./len(indices))
    results = []
    retained = []
    rejected = []
    for position, index in enumerate(indices):
        try:
            result = evaluator(samples[index])
            for name in ('RV', 'Transit', 'Astrometry'):
                base = getattr(reference, name)
                current = getattr(result, name)
                if (base is None) != (current is None):
                    raise ValueError('Posterior draws must retain the same enabled observables.')
                if base is not None and not np.array_equal(base.BJD, current.BJD):
                    raise ValueError('Posterior BJD grids differ; disable timestamp jitter.')
                if current is not None and np.any(current.noise):
                    raise ValueError('Disable simulated measurement noise for latent-model bands.')
            if set(result.TTV) != set(reference.TTV):
                raise ValueError('Posterior planet names differ from the reference.')
            if any(np.any(event.noise) for event in result.TTV.values()):
                raise ValueError('Disable simulated TTV measurement noise for latent-model bands.')
            results.append(result)
            retained.append(position)
        except (ValueError, RuntimeError) as error:
            if on_invalid == 'raise':
                raise ValueError(f'Posterior row {index}: {error}') from error
            rejected.append((int(index), str(error)))
    if not results:
        raise ValueError('No valid posterior models remain.')
    indices = indices[retained]
    weights = weights[retained]
    retained_mass = float(weights.sum())
    weights /= weights.sum()
    def curve(name, field):
        base = getattr(reference, name)
        if base is None:
            return None
        draws = [getattr(getattr(result, name), field) for result in results]
        return _band(base.BJD, getattr(base, field), draws, weights, levels, minimum_coverage)
    astro = {}
    if reference.Astrometry is not None:
        astro = {field: curve('Astrometry', field) for field in ('East', 'North', 'along_scan')}
    events = {}
    for name, base in reference.TTV.items():
        if not len(base.BJD):
            events[name] = {}
            continue
        if event_periods is None or name not in event_periods:
            raise ValueError('Supply event_periods[name] to align posterior transits by orbital cycle.')
        period = float(event_periods[name])
        if not np.isfinite(period) or period <= 0:
            raise ValueError('Event alignment periods must be finite and positive.')
        cycles = np.rint((base.BJD-base.BJD[0])/period).astype(int)
        if len(set(cycles)) != len(cycles):
            raise ValueError('Reference transits cannot be uniquely matched by the supplied period.')
        aligned = {field: [] for field in ('BJD', 'Depth', 'Duration', 'Duration_full', 'O_C')}
        for result in results:
            event = result.TTV[name]
            draw_cycles = np.rint((event.BJD-base.BJD[0])/period).astype(int)
            if len(set(draw_cycles)) != len(draw_cycles):
                raise ValueError('Posterior transits cannot be uniquely aligned by orbital cycle.')
            positions = {cycle: index for index, cycle in enumerate(draw_cycles)}
            for field in aligned:
                values = np.full(len(cycles), np.nan)
                for column, cycle in enumerate(cycles):
                    if cycle in positions:
                        index = positions[cycle]
                        # One common reference ephemeris for all O-C draws.
                        values[column] = event.BJD[index]-base.Linear_BJD[column] if field == 'O_C' else getattr(event, field)[index]
                aligned[field].append(values)
        events[name] = {field: _band(base.BJD, getattr(base, field), values, weights, levels, minimum_coverage)
                        for field, values in aligned.items()}
    combined_weights = np.bincount(indices, weights=weights, minlength=len(samples))
    effective_size = 1./np.sum(combined_weights**2)
    return ModelEnsemble(reference, curve('RV', 'model'), curve('Transit', 'model'), events, astro,
                         indices, weights, effective_size, tuple(rejected), retained_mass)


def GetPosteriorModels(session, *, source='mcmc', samples=None, weights=None,
                       max_samples=None, seed=42, levels=(1, 2, 3), minimum_coverage=.95,
                       on_invalid='raise', radius_mode=None, phase=None, **kwargs):
    """Propagate a session's MCMC/NS rows into reusable GetModel credible bands.

    Leave simulation noise and timestamp jitter disabled. Instrument/GP nuisance
    columns are marginalized through their correlated orbital sample rows, not
    added as deterministic signals. No fitting or sampler execution occurs here.
    """
    adapter = SessionModelAdapter(session, radius_mode=radius_mode, phase=phase)
    if samples is None:
        posterior = PosteriorSamplesFromSession(session, source)
        samples = posterior.values
        if weights is None:
            weights = posterior.weights
    reference = adapter(**kwargs)
    planets, _, _, _ = adapter.inputs()
    periods = {name: parameters['P'] for name, parameters in planets.items()}
    return GetModelEnsemble(lambda row: adapter(row, **kwargs), samples, reference=reference,
                            weights=weights, max_samples=max_samples, seed=seed, levels=levels,
                            minimum_coverage=minimum_coverage, on_invalid=on_invalid,
                            event_periods=periods)
