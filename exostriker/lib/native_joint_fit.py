"""Native bounded LM/Simplex fits of normalized Gaussian RV/LC/TTV/TDV lnL.

Python packs inputs and publishes the final result. Trial models, mass conversion,
likelihood, optimization, and covariance are evaluated in Fortran, without callbacks.
"""
from dataclasses import dataclass
import copy
import time
import numpy as np
from native_transit_model import _native, _LOCK, LD_LAWS, LD_COUNTS, session_jacobi_convention


@dataclass
class JointFitResult:
    """Hold the native solution, model samples and fit diagnostics."""

    parameters: np.ndarray
    models: np.ndarray
    variances: np.ndarray
    loglik: float
    components: np.ndarray  # RV, LC, timing, duration, normalized mid-transit flux
    covariance: np.ndarray
    errors: np.ndarray
    diagnostics: dict
    free_indices: np.ndarray
    meta: np.ndarray
    data: np.ndarray
    event_counts: np.ndarray
    event_table: np.ndarray
    rv: np.ndarray
    flux: np.ndarray
    rv_components: np.ndarray
    mean_anomaly: np.ndarray


@dataclass
class JointFitProblem:
    """Hold the packed data and parameter mapping for a native fit."""

    initial: np.ndarray
    free_indices: np.ndarray
    bounds: np.ndarray
    scales: np.ndarray
    packed: np.ndarray
    mapping: np.ndarray
    meta: np.ndarray
    data: np.ndarray
    sample_time: np.ndarray
    sample_owner: np.ndarray
    laws: np.ndarray
    config: np.ndarray
    settings: np.ndarray
    npl: int
    ignored: tuple = ()
    priors: np.ndarray = None
    parameter_labels: dict = None

    def loglik(self, parameters):
        """Sampler evaluation: only a scalar crosses the Fortran boundary."""

        parameters = np.asarray(parameters, float)
        if parameters.shape != self.initial.shape or not np.all(np.isfinite(parameters)):
            return -np.inf

        with _LOCK:
            value, status = _native().native_joint_loglik(
                parameters,
                self.packed,
                self.mapping,
                self.meta,
                self.data,
                self.sample_time,
                self.sample_owner,
                self.laws,
                self.config,
                self.settings,
                int(self.config[7]),
                self.npl,
            )
        return float(value) if status == 0 and np.isfinite(value) else -np.inf

    def solve(self, method='LM', *, maxiter=500, tolerance=1e-7, optimizer_options=None):
        """Run the chosen native optimizer and unpack its diagnostics."""

        if self.npl == 0 and method == 'LM':
            from dataclasses import replace
            # Preserve the RV LM convention: optimize means with fixed jitter.
            jitter_targets = np.r_[np.arange(5, 124, 2), np.arange(125, 664, 9)]
            jitter_indices = self.mapping[jitter_targets] - 1
            keep = ~np.isin(self.free_indices, jitter_indices[jitter_indices >= 0])
            if not np.all(keep):
                reduced = replace(
                    self,
                    free_indices=self.free_indices[keep],
                    bounds=self.bounds[keep],
                    scales=self.scales[keep],
                    ignored=self.ignored + tuple(self.free_indices[~keep]),
                )
                return reduced.solve(
                    method,
                    maxiter=maxiter,
                    tolerance=tolerance,
                    optimizer_options=optimizer_options,
                )

        methods = {'evaluate': 0, 'LM': 1, 'Simplex': 2, 'Nelder-Mead': 2}
        if method not in methods:
            raise ValueError('Select LM, Simplex, or evaluate.')
        if not np.isfinite(tolerance) or tolerance <= 0:
            raise ValueError('tolerance must be positive.')
        if int(maxiter) != maxiter or maxiter < 0:
            raise ValueError('maxiter must be a nonnegative integer.')

        settings = self.settings.copy()
        settings[6] = tolerance
        free = self.free_indices
        bounds = self.bounds
        scales = self.scales
        if not len(free):
            # f2py cannot marshal zero-dimensional work matrices on all NumPy builds.
            free = np.array([len(self.initial) - 1], np.int32)
            bounds = np.array([[self.initial[-1], self.initial[-1]]])
            scales = np.ones(1)
        with _LOCK:
            call_options = {'optimizer_options': optimizer_options} if optimizer_options is not None else {}
            if self.priors is not None and np.any(self.priors[:, [2, 5]]):
                call_options['prior_options'] = self.priors
            out = _native().native_joint_fit(
                self.initial,
                free + 1,
                bounds[:, 0],
                bounds[:, 1],
                scales,
                self.packed,
                self.mapping,
                self.meta,
                self.data,
                self.sample_time,
                self.sample_owner,
                self.laws,
                self.config,
                settings,
                methods[method] if len(self.free_indices) else 0,
                int(maxiter),
                int(self.config[7]),
                self.npl,
                **call_options,
            )
        par, mu, var, ll, comp, cov, err, diag, event_counts, event_table, physical_rv, physical_flux, parts, mean_anomaly = out
        if not len(self.free_indices):
            cov = cov[:0, :0]
            err = err[:0]
        if diag[0] in (2, 3):
            if diag[0] == 3 and diag[4] == 4:
                observed = self.data[self.meta[:, 0] <= 2, 0] + self.settings[0]
                span = f'observation times {observed.min():.12g} to {observed.max():.12g}; ' if len(observed) else ''
                timing = f'timing search {self.settings[3]:.12g} to {self.settings[4]:.12g}; ' if np.any(self.meta[:, 0] > 2) else ''
                raise ValueError(
                    'Native N-body model exceeded the 1,000,000 integration-step limit '
                    f'before optimization (model status 4). Epoch {self.settings[0]:.12g}; '
                    f'{span}{timing}requested integration step {self.settings[1]:.12g} days. '
                    'Check the orbital epoch, full versus reduced Julian dates, time coverage '
                    'and integration step. This limit is separate from optimizer iterations/evaluations.')
            raise ValueError(f'Native fit input/model failed (status {diag[0]}, model status {diag[4]}).')
        diagnostics = dict(
            backend='fortran',
            method=method,
            converged=bool(diag[0] == 0),
            fit_components=dict(zip(('timing', 'duration', 'depth'), map(bool, self.config[9:12]))),
            iterations=int(diag[1]),
            evaluations=int(diag[2]),
            physical_evaluations=dict(rv_transit=int(diag[5]), timing=int(diag[6]), shared=int(diag[7])),
            integration_steps=int(diag[8]),
            covariance_status='available' if len(self.free_indices) and diag[3] == 0 and (method != 'evaluate') else 'unavailable',
            covariance_method='observed Hessian' if methods[method] == 2 and (not self.config[0]) else 'Gaussian Fisher',
            covariance_reason='optimizer has not converged' if diag[0] != 0 else 'active bounds, singular/invalid curvature, or no error calculation' if diag[3] != 0 or method == 'evaluate' else '',
            ignored_parameters=self.ignored,
            objective='log posterior' if 'prior_options' in call_options else 'normalized Gaussian lnL; no priors',
        )
        diagnostics['covariance_code'] = int(diag[3])
        if method == 'evaluate' or not len(self.free_indices):
            diagnostics['covariance_reason'] = 'Model evaluation does not estimate parameter errors.'
        elif diag[0] == 0 and diag[3] != 0:
            reasons = {
                2: 'Model derivatives could not be evaluated for covariance.',
                3: 'A central covariance step would cross a parameter bound.',
                4: 'A fitted parameter is at a bound; symmetric local errors are unavailable.',
                5: 'The posterior information matrix is singular or not positive definite.',
            }
            diagnostics['covariance_reason'] = reasons.get(
                int(diag[3]),
                diagnostics['covariance_reason'],
            )
        if len(self.free_indices) and diag[3] != 0 and (method != 'evaluate'):
            q = (par[self.free_indices] - self.initial[self.free_indices]) / self.scales
            margin = np.maximum(
                2e-3,
                np.abs(q) * 2e-4,
            ) if methods[method] == 2 and (not self.config[0]) else np.full(len(q), 1e-7)
            near = ((par[self.free_indices] - self.bounds[:, 0]) / self.scales <= margin) | ((self.bounds[:, 1] - par[self.free_indices]) / self.scales <= margin)
            affected = []
            for k in np.flatnonzero(near):
                index = int(self.free_indices[k])
                label = (self.parameter_labels or {}).get(index, f'parameter[{index}]')
                affected.append(f'{label} = {par[index]:.12g} (bounds [{self.bounds[k, 0]:.12g}, {self.bounds[k, 1]:.12g}])')
            diagnostics['covariance_bound_parameters'] = affected
            if affected:
                diagnostics['covariance_reason'] += ' At or near bounds: ' + '; '.join(affected) + '.'
            physical = []
            for dataset in np.unique(self.meta[self.meta[:, 0] == 2, 1]) - 1:
                if self.laws[dataset] != 2:
                    continue
                offset = 10 * self.npl + 124 + 9 * dataset + 2
                indices = self.mapping[offset:offset + 2] - 1
                if np.any(indices < 0) or not np.any(np.isin(indices, self.free_indices)):
                    continue
                total = float(np.sum(par[indices]))
                if 1.0 - 0.001 <= total <= 1.0 + 1e-12:
                    physical.append(f'Transit dataset {dataset+1}: quadratic limb-darkening '
                        f'u1 + u2 = {total:.9g}, near the physical limit of 1 (zero limb intensity).')
            diagnostics['covariance_physical_constraints'] = physical
            if physical:
                diagnostics['covariance_reason'] += ' ' + ' '.join(physical)
        logprior = 0.
        if 'prior_options' in call_options:
            import RV_mod as rv
            indices = np.flatnonzero(np.any(self.priors[:, [2, 5]], axis=1))
            logprior = float(rv.lnprior(par[indices], np.tile([-np.inf, np.inf], (len(indices), 1)), [self.priors[indices, :3], self.priors[indices, 3:]]))
        diagnostics.update(loglik=float(ll), logprior=logprior, logposterior=float(ll) + logprior)
        if 'prior_options' in call_options:
            diagnostics['covariance_method'] += ' + prior Hessian'
        period_indices = self.mapping[1:10 * self.npl:10] - 1
        restore = np.argsort(np.argsort(par[period_indices], kind='stable'), kind='stable')
        return JointFitResult(par,mu,var,float(ll),comp,cov,err,diagnostics,self.free_indices.copy(),
                              self.meta.copy(),self.data.copy(),event_counts[restore],event_table[:,restore,:],physical_rv,physical_flux,parts[:,restore],mean_anomaly[restore])



def has_enabled_priors(s):
    """Whether an enabled fitted coordinate has a Gaussian/Jeffreys prior."""
    s.prepare_for_mcmc(rtg=s.rtg)
    return any((bool(row[2]) for rows in (s.nr_pr_for_mcmc, s.jeff_pr_for_mcmc) for row in rows))


def unsupported_reason(s):
    """Explain why this session needs a different likelihood backend."""

    from magnetic_cycle import active as cycle_active
    if cycle_active(s):
        return 'deterministic magnetic cycles (Python likelihood)'
    if not any(s.use_planet) and s.type_fit.get('TTV'):
        return 'timing fits without planets'
    if s.type_fit.get('AST', False):
        return 'astrometry'
    if s.rtg[1] or s.rtg[3]:
        return 'Gaussian processes'
    if s.gr_flag:
        return 'general relativity'
    if s.optim_AMD_stab or s.optim_Nbody_stab:
        return 'stability penalties'
    if s.get_TTVs[0]:
        return 'independent per-event transit-time shifts'
    if bool(np.asarray(s.mixed_fit[0]).item()):
        return 'mixed dynamical subsets'
    if any((s.omega_dot[i] != 0 or s.omega_dot_use[i] for i in np.flatnonzero(s.use_planet))):
        return 'apsidal precession'
    return None


def prepare(
    s,
    *,
    free_indices=None,
    bounds=None,
    exposure_time=0.,
    supersample_factor=1,
    phase=None,
    use_priors=True,
):
    """Pack the established full inference vector; sparse GUI slots stay intact."""
    reason = unsupported_reason(s)
    if reason:
        raise ValueError('Native joint fitting does not support ' + reason + '. Select SciPy for this fit.')
    import RV_mod as rv
    rv.migrate_ttv_datasets(s)
    from RV_mod.functions import ttv_component_selection
    selected_components = ttv_component_selection(getattr(s, 'ttv_fit_components', None))
    s.prepare_for_mcmc(rtg=s.rtg)
    slots = list(s.mcmc_planet_slots)
    ds = list(s.mcmc_dataset_slots)
    n = len(slots)
    nd = len(ds)
    initial = np.array(s.parameters, float)
    nf = len(initial)
    pv = np.zeros(10 * n + 684)
    mapping = np.zeros(len(pv), np.int32)
    scales = np.maximum(abs(initial) * .02, 1e-3)
    def link(target,index,value=None,scale=None):
        """Link an inference parameter to its native model position."""

        mapping[target] = index + 1
        pv[target] = initial[index] if value is None else value
        if scale is not None:
            scales[index] = scale
    for (j, slot) in enumerate(slots):
        start = 2 * nd + 7 * j
        end = 2 * nd + 7 * n + 2 + 3 * j
        for k in range(7):
            link(
                10 * j + k,
                start + k,
                scale=[max(abs(s.K[slot]) * .02, 1.), s.P[slot] * .002, .01, .5, .5, .2, .5][k],
            )
        for k in range(3):
            link(
                10 * j + 7 + k,
                end + k,
                scale=[s.P[slot] * .002, max(abs(s.pl_a[slot]) * .02, .1), .005][k],
            )
    off = 10 * n
    pv[off:off + 4] = [s.params.stellar_mass, s.stellar_radius, s.rv_lintr, s.rv_quadtr]
    link(off, nf - 1, scale=max(s.params.stellar_mass * .02, 1e-3))
    for k in range(2):
        link(off + 2 + k, 2 * nd + 7 * n + k)
    for (j, d) in enumerate(ds):
        link(off + 4 + 2 * d, j, scale=max(np.std(s.rv_data_sets[d][1]) * .1, 1.))
        link(off + 5 + 2 * d, nd + j, scale=max(np.mean(s.rv_data_sets[d][2]), 1.))
    loaded = [d for d in range(60) if len(s.tra_data_sets[d])]
    count = len(loaded)
    start = 2 * nd + 10 * n + 2
    laws = np.zeros(60, np.int32)
    for d in range(60):
        target = off + 124 + 9 * d
        law = s.ld_m[d]
        laws[d] = LD_LAWS[law]
        coeff = np.zeros(4)
        if law != 'uniform':
            coeff[:LD_COUNTS[laws[d]]] = getattr(
                s,
                dict(linear='ld_u_lin', quadratic='ld_u_quad', nonlinear='ld_u_nonlin')[law],
            )[d]
        dil = getattr(s, 'tra_dil', {}).get(d, s.tra_data_sets[d][8] if d in loaded else 1.)
        pv[target:target + 9] = [
            s.tra_off[d],
            s.tra_jitt[d],
            *coeff,
            s.tra_lintr[d],
            s.tra_quadtr[d],
            dil,
        ]
        if d in loaded:
            j = loaded.index(d)
            for (k, src) in enumerate([start + j, start + count + j]):
                link(target + k, src, scale=.0001)
            for k in range(LD_COUNTS[laws[d]]):
                link(target + 2 + k, start + 4 * count + n + int(s.ld_gr_ind[d]) + k, scale=.02)
            link(target + 6, start + 2 * count + j, scale=1e-4)
            link(target + 7, start + 3 * count + j, scale=1e-4)
            if d in s.tra_dil_parameter_indices:
                link(target + 8, s.tra_dil_parameter_indices[d], scale=.02)
    types = s.type_fit
    haslc = types.get('Transit', False) and bool(loaded)
    hasrv = types.get('RV', False) and bool(ds)
    epoch = float(s.epoch if haslc or hasrv else s.ttv_times[0])
    if hasrv or haslc:
        from model_epoch import resolve_epoch
        epoch = resolve_epoch(s, hasrv, haslc)
    if not np.isfinite(epoch):
        raise ValueError('The orbital epoch must be finite.')
    metas = []
    observations = []
    for (kind, sets, enabled) in [(1, s.rv_data_sets, hasrv), (2, s.tra_data_sets, haslc)]:
        if not enabled:
            continue
        for (d, row) in sets.items():
            if not len(row):
                continue
            t = np.asarray(row[0], float)
            y = np.asarray(row[1], float)
            err = np.asarray(row[2], float)
            if not len(t):
                continue
            if t.shape != y.shape or t.shape != err.shape:
                raise ValueError('Dataset arrays have different lengths.')
            x = np.asarray(row[3], float) if kind == 2 and len(row) > 10 and row[10] else t - t[0]
            metas.extend([[kind, d + 1, 1, 1]] * len(t))
            observations.extend(np.column_stack((t - epoch, y, err, x)))
    nr = len(metas)
    if types.get('TTV', False):
        for (d, row) in s.ttv_data_sets.items():
            if not len(row) or not row[4]:
                continue
            slot = int(row[3]) - 1
            if row[3] != int(row[3]) or slot not in slots:
                raise ValueError('Invalid TTV planet assignment.')
            j = slots.index(slot)
            target = off + 664 + 2 * d
            pv[target:target + 2] = [value if np.isfinite(value) else 0. for value in (row[10], row[11])]
            # With only TDV data, optional shared density and planet radius are
            # enabled through the same physical GUI parameters as light curves.
            if s.rho and s.pl_a_use[slots[0]]:
                link(target, 2 * nd + 7 * n + 3, scale=max(s.pl_a[slots[0]] * .02, .1))
            if s.pl_rad_use[slot]:
                link(target + 1, 2 * nd + 7 * n + 4 + 3 * j, scale=.005)
            nn = np.asarray(row[0])
            t = np.asarray(row[1])
            err = np.asarray(row[2])
            if not len(t) or nn.shape != t.shape or err.shape != t.shape or np.any(nn != nn.astype(int)) or np.any(nn < 1):
                raise ValueError('Invalid TTV times, errors, or one-based event numbers.')
            kinds = [(3, t, err if selected_components[0] else np.ones(len(t)))]
            if row[13] and selected_components[1]:
                kinds.append((4, np.asarray(row[8]), np.asarray(row[9])))
            if row[12] and selected_components[2]:
                kinds.append((5, np.asarray(row[6]), np.asarray(row[7])))
            for (kind, y, e) in kinds:
                if y.shape != t.shape or e.shape != t.shape:
                    raise ValueError('Invalid depth/duration array lengths.')
                metas.extend([[kind, d + 1, j + 1, int(num)] for num in nn])
                observations.extend(np.column_stack((t, y, e, np.zeros(len(t)))))
    if not metas:
        raise ValueError('No observations for the selected fit components.')
    meta = np.asarray(metas, np.int32)
    data = np.asarray(observations, float)
    if not np.all(np.isfinite(data)) or np.any(data[:, 2] <= 0):
        raise ValueError('Require finite data and strictly positive errors.')
    exposure = np.broadcast_to(np.asarray(exposure_time, float), (60,))
    samples = np.broadcast_to(np.asarray(supersample_factor, float), (60,))
    if not np.all(np.isfinite(exposure)) or np.any(exposure < 0) or (not np.all(np.isfinite(samples))) or np.any(samples != samples.astype(int)) or np.any(samples < 1) or np.any(samples > 1001):
        raise ValueError('Invalid exposure duration or supersampling factor.')
    expanded = []
    owners = []
    for i in range(nr):
        d = meta[i, 1] - 1
        cnt = int(samples[d]) if meta[i, 0] == 2 else 1
        expanded.extend(data[i, 0] + exposure[d] * (np.linspace(-.5, .5, cnt) if cnt > 1 else np.zeros(1)))
        owners.extend([i + 1] * cnt)
    order = np.argsort(expanded, kind='stable')
    # f2py accepts a one-entry unused buffer for timing-only fits.
    expanded = np.asarray(expanded)[order] if nr else np.array([0.])
    owners = np.asarray(owners, np.int32)[order] if nr else np.array([1], np.int32)
    dyn = bool(s.mod_dynamical)
    phase = getattr(s, 'native_joint_phase', 'auto') if phase is None else phase
    if phase == 'auto':
        phase = 't0' if haslc and (not dyn) else 'M0'
    if phase not in ('M0', 't0') or (phase == 't0' and (dyn or not haslc)):
        raise ValueError('Use M0 for N-body or timing-only fits; t0 is available for Keplerian light curves.')
    phase = phase == 't0'
    propagation = getattr(s, 'native_joint_propagation', 'shared')
    if propagation not in ('shared', 'separate'):
        raise ValueError('native_joint_propagation must be shared or separate.')
    cfg = np.array(
        [
            dyn,
            session_jacobi_convention(s) == 'ttvfast',
            s.hkl,
            s.copl_incl,
            s.jit_flag,
            s.rho,
            phase,
            min(
                200000,
                max(
                    64,
                    int(2 * (s.ttv_times[2] - s.ttv_times[0]) / min((s.P[i] for i in slots), default=1.) + 32),
                ),
            ),
            min(ds, default=0) + 1,
            *selected_components,
            propagation == 'shared',
        ],
        np.int32,
    )
    if not any((kind <= 2 or selected_components[kind - 3] for kind in meta[:, 0])):
        raise ValueError('Select at least one fitted timing, duration, depth, RV, or transit component.')
    settings = np.array([epoch, s.time_step_model, s.dyn_model_accuracy * 1e-13, s.ttv_times[0], s.ttv_times[2], 0., 1e-7, 0.])
    if settings[1] <= 0 or not 1e-15 <= settings[2] <= 0.001:
        raise ValueError('Invalid native integrator settings.')
    if free_indices is None:
        selected = list(map(int, s.f_for_mcmc))
        btable = {i: b for (i, b) in zip(selected, s.b_for_mcmc)}
        # prepare_for_mcmc historically enables these only for LC fits.
        if types.get('TTV') and any((len(r) and r[4] and (r[12] and selected_components[2] or (r[13] and selected_components[1])) for r in s.ttv_data_sets.values())):
            for (j, slot) in enumerate(slots):
                for k,flag,bnd in [
                    (1,s.pl_a_use[slot] and s.rho and j==0,s.stellar_rho_bounds),
                    (2,s.pl_rad_use[slot],s.pl_rad_bound[slot]),
                ]:
                    index = 2 * nd + 7 * n + 2 + 3 * j + k
                    if flag and index not in selected:
                        selected.append(index)
                        btable[index] = bnd
        full_bounds = np.asarray([btable[i] for i in selected], float).reshape(-1, 2)
    else:
        selected = list(map(int, free_indices))
        full_bounds = np.asarray(bounds, float).reshape(-1, 2)
        if len(selected) != len(full_bounds):
            raise ValueError('Provide bounds for every free index.')
    if any((i < 0 or i >= nf for i in selected)) or len(set(selected)) != len(selected):
        raise ValueError('Invalid free parameter indices.')
    inactive = set()
    for j in range(n):
        if phase:
            inactive.add(2 * nd + 7 * j + 4)
        else:
            inactive.add(2 * nd + 7 * n + 2 + 3 * j)
        if dyn and (not s.rho):
            inactive.add(2 * nd + 7 * n + 3 + 3 * j)
        if s.rho and j > 0:
            inactive.add(2 * nd + 7 * n + 3 + 3 * j)
        if dyn and s.copl_incl and (j > 0):
            inactive.add(2 * nd + 7 * j + 5)
    if s.jit_flag:
        inactive.update((nd + j for j in range(1, nd)))
    if not hasrv:
        inactive.update([*range(2 * nd), 2 * nd + 7 * n, 2 * nd + 7 * n + 1])
    if not haslc:
        inactive.update(mapping[off + 124:off + 664][mapping[off + 124:off + 664] > 0] - 1)
        # LD grouping can pack fewer entries than this upper bound; explicitly
        # preserve anything mapped into timing metadata, then remove unused maps.
    if n == 0:
        # A flat physical model cannot constrain stellar mass or limb darkening.
        inactive.add(nf - 1)
        for d in loaded:
            inactive.update(mapping[off + 124 + 9 * d + 2:off + 124 + 9 * d + 6] - 1)
            if not s.tra_lintr_use[d] and (not s.tra_quadtr_use[d]) and (s.tra_lintr[d] == 0) and (s.tra_quadtr[d] == 0):
                inactive.add(mapping[off + 124 + 9 * d + 8] - 1)
    used = set(mapping[mapping > 0] - 1)
    # Remove only genuinely inactive coordinates, preserving shared TTV entries.
    if not haslc:
        timing_used = set(mapping[off + 664:][mapping[off + 664:] > 0] - 1)
        inactive -= timing_used
    keep = [k for (k, i) in enumerate(selected) if i in used and i not in inactive and (full_bounds[k, 0] < full_bounds[k, 1])]
    ignored = tuple((i for (k, i) in enumerate(selected) if k not in keep))
    selected = np.asarray(selected, np.int32)[keep]
    full_bounds = full_bounds[keep]
    if np.any(selected < 0) or np.any(selected >= nf) or len(set(selected)) != len(selected):
        raise ValueError('Invalid free parameter indices.')
    # Jitters are positive parameters in native minimization.
    jitter_indices = {nd + j for j in range(nd)} | {start + count + j for j in range(count)}
    for (k, i) in enumerate(selected):
        if i in jitter_indices:
            full_bounds[k, 0] = max(0., full_bounds[k, 0])
    if s.rho and n:
        # The Stellar hard bounds remain authoritative even for explicit API bounds.
        density_index = 2 * nd + 7 * n + 3
        for (k, index) in enumerate(selected):
            if index == density_index:
                full_bounds[k, 0] = max(full_bounds[k, 0], s.stellar_rho_bounds[0])
                full_bounds[k, 1] = min(full_bounds[k, 1], s.stellar_rho_bounds[1])
                if full_bounds[k, 0] >= full_bounds[k, 1]:
                    raise ValueError('Stellar density bounds have no overlap with the requested optimizer bounds.')
    outside = (initial[selected] < full_bounds[:, 0]) | (initial[selected] > full_bounds[:, 1])
    if np.any(outside):
        names = dict(zip(map(int, s.f_for_mcmc), s.e_for_mcmc))
        details = []
        for k in np.flatnonzero(outside):
            index = int(selected[k])
            label = 'stellar density rho_star [kg/m^3]' if s.rho and index == 2 * nd + 7 * n + 3 else str(names.get(index, f'parameter[{index}]'))
            details.append(f'{label}: {initial[index]:.12g} outside [{full_bounds[k, 0]:.12g}, {full_bounds[k, 1]:.12g}]')
        hint = ' Set the starting stellar density within these limits, or adjust the bounds under Limits and Priors > Stellar.' if s.rho else ' Set the starting values within their bounds, or adjust the bounds.'
        raise ValueError('Initial parameters are outside their bounds: ' + '; '.join(details) + '.' + hint)
    if not np.all(np.isfinite(initial)) or np.any(np.isnan(full_bounds)):
        raise ValueError('Nonfinite parameters or invalid bounds.')
    priors = np.zeros((nf, 6))
    if use_priors:
        prior_table = {int(i): np.r_[normal, jeff] for (i, normal, jeff) in zip(s.f_for_mcmc, s.nr_pr_for_mcmc, s.jeff_pr_for_mcmc)}
        if s.rho and n:
            prior_table[2 * nd + 7 * n + 3] = np.r_[s.stellar_rho_norm_pr, s.stellar_rho_jeff_pr]
        for index in selected:
            row = prior_table.get(int(index), np.zeros(6)).copy()
            if not row[2]:
                row[:3] = 0.
            if not row[5]:
                row[3:] = 0.
            if not np.all(np.isfinite(row)) or (row[2] and row[1] <= 0):
                raise ValueError('Enabled Gaussian priors require a finite mean and a positive finite standard deviation.')
            if row[5] and (row[3] <= 0 or row[4] <= row[3] or initial[index] <= 0):
                raise ValueError('Enabled Jeffreys priors require positive parameter values and 0 < lower < upper.')
            priors[index] = row
    labels = dict(zip(map(int, s.f_for_mcmc), map(str, s.e_for_mcmc)))
    return JointFitProblem(
        initial,
        selected,
        full_bounds,
        scales[selected],
        pv,
        mapping,
        meta,
        data,
        expanded,
        owners,
        laws,
        cfg,
        settings,
        n,
        ignored,
        priors,
        labels,
    )


def timing_results(s, problem, result):
    """Final event tables in the established five-section GUI format."""
    from native_transit_model import _mass_axes
    from RV_mod.functions import transit_duration_depth_from_rsky_vsky
    slots = list(s.mcmc_planet_slots)
    rows = np.array([[getattr(s, name)[i] for name in ('K', 'P', 'e', 'w', 'M0', 'i', 'Node')] for i in slots])
    # Published M0/e/w have already been synchronized from h/k/lambda or t0.
    if s.mod_dynamical and s.copl_incl:
        rows[:, 5] = rows[0, 5]
    masses, axes = _mass_axes(rows, s.params.stellar_mass)
    order = np.argsort(rows[:, 1], kind='stable')
    inverse = np.argsort(order)
    rs = s.stellar_radius
    if s.rho:
        rho = s.pl_a[slots[0]]
        rs = (3 * s.params.stellar_mass * 1.32712440018e20 / (4 * np.pi * 6.67408e-11 * rho)) ** (1 / 3) / 6.957e8
    if s.type_fit.get('Transit'):
        if s.rho:
            scaled = (6.67408e-11 * s.pl_a[slots[0]] * (rows[:, 1] * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
        else:
            scaled = np.array([s.pl_a[i] for i in slots])
        input_axes = scaled * rs
    else:
        input_axes = axes * 1.49597892e11 / 6.957e8
    mask = result.meta[:, 0] == 3
    labels = result.meta[mask, 1]
    planet = result.meta[mask, 2] - 1
    numbers = result.meta[mask, 3]
    if not len(numbers):
        return None
    if getattr(result, 'event_table', None) is not None:
        counts = result.event_counts[order]
        table = result.event_table[:, order, :]
    else:
        # Compatibility with results saved before event tables were returned.
        with _LOCK:
            counts, table, *rest = _native().joint_timing_likelihood(
                rows[order, 1:],
                masses[order],
                s.params.stellar_mass,
                problem.settings[0],
                problem.settings[3],
                problem.settings[4],
                int(s.mod_dynamical),
                input_axes[order],
                problem.settings[1],
                min(problem.settings[2], 1e-14),
                1000000,
                int(problem.config[1]),
                512,
                1000000,
                int(problem.config[7]),
                labels,
                inverse[planet] + 1,
                numbers,
                result.data[mask, 0],
                np.ones(len(numbers)),
            )
        if rest[-1]:
            raise ValueError(f'Final native timing table failed: {rest[-1]}')
    cd = {d: [] for d in range(10)}
    cm = {d: [] for d in range(10)}
    short = [[], []]
    full = [[], [], [], [], [], []]
    for d in range(10):
        selected = labels == d + 1
        if not np.any(selected):
            continue
        row = s.ttv_data_sets[d]
        j = int(planet[selected][0])
        jj = inverse[j]
        event = table[:counts[jj], jj]
        ix = numbers[selected] - 1
        rho = row[10]
        rr = row[11]
        if s.rho and s.pl_a_use[slots[0]]:
            rho = s.pl_a[slots[0]]
        if s.pl_rad_use[slots[j]]:
            rr = s.pl_rad[slots[j]]
        if s.type_fit.get('Transit'):
            rho = 3 * s.params.stellar_mass * 1.9884098709677423e30 / (4 * np.pi * (rs * 6.957e8) ** 3)
            rr = s.pl_rad[slots[j]]
        # Legacy depth is normalized flux, and duration is chord T14 in days.
        if rho > 0:
            duration, _, loss, _, _ = transit_duration_depth_from_rsky_vsky(
                event[:, 2],
                event[:, 3],
                rho,
                s.params.stellar_mass,
                rr,
            )
            depth = 1 - np.asarray(loss)
        else:
            duration = np.zeros(len(event))
            depth = np.ones(len(event))
        model = result.models[mask][selected]
        cd[d] = [
            numbers[selected].copy(),
            model.copy(),
            np.asarray(row[2], dtype=float).copy(),
            event[ix, 2].copy(),
            event[ix, 3].copy(),
            depth[ix].copy(),
            np.asarray(duration)[ix].copy(),
        ]
        cm[d] = [
            event[:, 0].astype(int),
            event[:, 1].copy(),
            event[:, 2].copy(),
            event[:, 3].copy(),
            depth,
            np.asarray(duration),
        ]
        short = [numbers[selected].copy(), model.copy()]
        full = cm[d]
    return [float(sum(result.components[2:])), short, full, cd, cm]


def likelihood_from_joint(session, result):
    """Adapt exact native best-trial models, without a second propagation."""
    from native_likelihood import NativeLikelihoodResult
    from native_transit_model import _mass_axes, jacobi_axis_factor
    mask = result.meta[:, 0] <= 2
    kind = result.meta[mask, 0] - 1
    labels = result.meta[mask, 1] - 1
    data = result.data[mask]
    slots = list(np.flatnonzero(session.use_planet))
    names = ('K', 'P', 'e', 'w', 'M0', 'i', 'Node')
    rows = np.array([[getattr(session, name)[i] for name in names] for i in slots]).reshape(-1, 7)
    if session.hkl:
        h = np.array([session.e_sinw[i] for i in slots])
        k = np.array([session.e_cosw[i] for i in slots])
        rows[:, 2] = np.hypot(h, k)
        rows[:, 3] = np.degrees(np.arctan2(h, k))
    if slots and session.mod_dynamical and session.copl_incl:
        rows[:, 5] = rows[0, 5]
    rows[:, 4] = result.mean_anomaly
    mass, axes = _mass_axes(rows, session.params.stellar_mass)
    radius = session.stellar_radius
    if slots and session.rho:
        radius = (3 * session.params.stellar_mass * 1.32712440018e20 / (4 * np.pi * 6.67408e-11 * session.pl_a[slots[0]])) ** (1 / 3) / 6.957e8
    if session.mod_dynamical:
        axes *= jacobi_axis_factor(
            rows,
            mass,
            session.params.stellar_mass,
            session_jacobi_convention(session),
        )
        axis = axes * 1.49597892e11 / (radius * 6.957e8)
    elif slots and session.rho:
        axis = (6.67408e-11 * session.pl_a[slots[0]] * (rows[:, 1] * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
    else:
        axis = np.array([session.pl_a[i] for i in slots], float)
    model = result.models[mask]
    variance = result.variances[mask]
    residual = data[:, 1] - model
    point = -.5 * (residual ** 2 / variance + np.log(variance) + np.log(2 * np.pi))
    by_dataset = np.zeros((60, 2))
    np.add.at(by_dataset, (labels, kind), point)
    chi2 = np.bincount(kind, weights=residual ** 2 / variance, minlength=2)
    shared_slot = min((d for (d, v) in session.rv_data_sets.items() if len(v)), default=0)
    jitter = np.array([abs(session.jitt[shared_slot if session.jit_flag else int(d)]) if k == 0 else abs(session.tra_jitt[int(d)]) for (d, k) in zip(labels, kind)])
    offsets = np.array([session.tra_off[int(d)] for d in labels])
    return NativeLikelihoodResult(data[:,0]+session.epoch,labels,result.rv[mask],result.flux[mask],
        result.flux[mask]-offsets,jitter,axis,np.array([session.pl_rad[i] for i in slots])*radius,
        result.rv_components[mask],result.mean_anomaly.copy(),None,kind,data[:,1],data[:,2],data[:,3],
        model,residual,variance,point,by_dataset,float(sum(result.components[:2])),
        float(result.components[0]),float(result.components[1]),chi2)


def publish_results(s, par, opt, tr_model, errors, flags):
    """Copy a native solution into the established session result fields."""

    from native_inference import publish_results as publish_photometry
    result = opt['native_joint_result']
    problem = opt['native_joint_problem']
    # Resolve the epoch before constructing dense plotting grids as well.
    if not s.type_fit.get('RV') and s.type_fit.get('Transit') and (
            s.epoch==0. or s.epoch==getattr(s,'_transit_default_epoch',None)):
        s._transit_default_epoch = float(problem.settings[0])
    s.epoch = problem.settings[0]
    publish_photometry(s, par, opt, tr_model, errors, flags)
    s.native_joint_phase_used = opt['native_joint_phase']
    if s.type_fit.get('TTV'):
        s.ttv_results = timing_results(s, problem, result)
        s.ttv_model_backend = 'native_joint'
    r = s.fit_results
    r.native_joint_result = result
    r.native_fit_diagnostics = result.diagnostics
    r.native_covariance = result.covariance.copy()
    r.native_parameter_errors = result.errors.copy()
    r.native_free_indices = result.free_indices.copy()
    r.native_component_loglik = result.components.copy()
    r.native_loglik = result.loglik
    r.native_ttv_loglik = float(sum(result.components[2:]))
    if s.type_fit.get('TTV'):
        r.native_ttv_results = s.ttv_results
    s.loglik = result.loglik
    active = np.array([kind <= 2 or problem.config[kind + 6] for kind in result.meta[:, 0]], dtype=bool)
    r.native_joint_chi2 = float(np.sum(((result.data[:, 1] - result.models) ** 2 / result.variances)[active]))
    r.native_joint_ndata = int(sum(active))
    r.native_joint_dof = int(sum(active)) - len(result.free_indices)
    if not s.type_fit.get('RV'):
        r.chi2 = r.native_joint_chi2
        r.Ndata = r.ndata = int(sum(active))
        r.mfit = len(result.free_indices)
        r.dof = r.Ndata - r.mfit
        r.reduced_chi2 = r.chi2 / r.dof if r.dof > 0 else np.nan
    # In joint RV fits retain the legacy RV component's loglik field.
    r.loglik = result.components[0] if s.type_fit.get('RV') else result.loglik
    if not s.type_fit.get('RV') and (not s.type_fit.get('Transit')):
        mask = result.meta[:, 0] == 3
        residual = result.data[mask, 1] - result.models[mask]
        r.rms = float(np.sqrt(np.mean(residual ** 2)))
        r.wrms = float(np.sqrt(np.average(residual ** 2, weights=1 / result.variances[mask])))
    for name in ('Ndata', 'mfit', 'chi2', 'loglik', 'reduced_chi2', 'rms'):
        r.res_dict[name] = getattr(r, name)
        r.res[name] = getattr(r, name)
    s.chi2 = r.chi2
    s.reduced_chi2 = r.reduced_chi2
    s.rms = r.rms
    s.native_fit_diagnostics = result.diagnostics
    s.likelihood_backend_used = 'native'
    s.optimizer_backend_used = 'external' if opt.get('native_external_errors') else 'fortran'
    s.rvmod = copy.deepcopy(r)


def run(s, method=None, *, maxiter=None, tolerance=1e-7, free_indices=None, bounds=None, **options):
    """Fit and publish; None uses configured Fortran LM/Simplex stages.

    Enabled Gaussian/Jeffreys priors give a bounded MAP fit; use_priors=False
    explicitly requests maximum likelihood. Reported loglik remains likelihood.
    Non-convergence is reported in diagnostics; unavailable covariance is NaN.
    GP/astrometry and other unsupported physics raise, never silently disappear.
    """
    import RV_mod as rv
    start = time.perf_counter()
    from fortran_optimizer_options import ensure_options, option_vector
    configured = ensure_options(s)
    problem = prepare(s, free_indices=free_indices, bounds=bounds, **options)
    if getattr(s, 'init_fit', False):
        # Validate the configured starting bounds before removing the free
        # coordinates for a model-only evaluation.
        if s.rho and problem.npl:
            index = 2 * len(s.mcmc_dataset_slots) + 7 * problem.npl + 3
            value = problem.initial[index]
            low, high = s.stellar_rho_bounds
            if not low <= value <= high:
                raise ValueError(f'Stellar density rho_star [kg/m^3]: {value:.12g} outside '
                    f'[{low:.12g}, {high:.12g}]. Set the starting stellar density within these limits, '
                    'or adjust the bounds under Limits and Priors > Stellar.')
        problem.free_indices = np.empty(0, dtype=np.int32)
        problem.bounds = np.empty((0, 2))
        problem.scales = np.empty(0)
        stages = ['evaluate']
    elif method is not None:
        stages = [method]
    else:
        stages = []
        for name,count in [
            (s.Fortran_min_use_1,s.Fortran_min_N_use_1),
            (s.Fortran_min_use_2,s.Fortran_min_N_use_2),
        ]:
            if name not in (
                'LM',
                'Simplex',
                'Nelder-Mead',
            ) or isinstance(count,(bool,np.bool_)) or not isinstance(count,(int,np.integer)) or count<0:
                raise ValueError('Select LM/Simplex and nonnegative integer Fortran stage counts.')
            stages.extend([name] * int(count))
        if not stages:
            stages = ['evaluate']
    history = []
    for stage in stages:
        method_key = 'lm' if stage == 'LM' else 'simplex'
        settings = configured[method_key]
        limit = maxiter if maxiter is not None else settings['maxiter'] or getattr(
            s,
            'native_joint_maxiter',
            5000,
        )
        for restart in range(settings['restarts'] or 1 if method_key == 'simplex' and stage != 'evaluate' else 1):
            result = problem.solve(
                stage,
                maxiter=limit,
                tolerance=tolerance,
                optimizer_options=option_vector(s, method_key),
            )
            result.diagnostics['optimizer_options'] = settings.copy()
            history.append(result.diagnostics.copy())
            problem.initial = result.parameters.copy()
    result.diagnostics.update(stages=history, elapsed_seconds=time.perf_counter() - start)
    s.bound_error = False
    s.bound_error_msg = ''
    opt = dict(
        rho=s.rho,
        native_joint_phase='t0' if problem.config[6] else 'M0',
        native_joint_result=result,
        native_joint_problem=problem,
        native_context=copy.copy(s),
        planet_slots=s.mcmc_planet_slots,
        dataset_slots=s.mcmc_dataset_slots,
        dilution_indices=getattr(s, 'tra_dil_parameter_indices', {}),
        get_TTVs=s.get_TTVs,
        dt=problem.settings[1] * 86400.,
        eps=problem.settings[2],
        native_joint_model_options={k: options[k] for k in ('exposure_time', 'supersample_factor') if k in options},
    )
    # Dedicated TTV-only elements retain the legacy ttv_times[0] reference epoch.
    opt['native_context'].epoch = problem.settings[0]
    tr_model = np.array([s.ld_m, s.ld_u, s.ld_gr, s.ld_gr_ind], dtype=object)
    flags = result.free_indices
    label_map = dict(zip(map(int, s.f_for_mcmc), s.e_for_mcmc))
    nd = len(s.mcmc_dataset_slots)
    for (j, slot) in enumerate(s.mcmc_dataset_slots):
        label_map[j] = s.rvoff_str[slot]
        label_map[nd + j] = s.jitt_str[slot]
    names = (
        'K',
        'P',
        'e_sinw',
        'e_cosw',
        'lamb',
        'i',
        'Node',
    ) if s.hkl else ('K', 'P', 'e', 'w', 'M0', 'i', 'Node')
    for (j, slot) in enumerate(s.mcmc_planet_slots):
        for (k, name) in enumerate(names):
            label_map[2 * nd + 7 * j + k] = getattr(s, name + '_str')[slot]
    for (j, slot) in enumerate(s.mcmc_planet_slots):
        for (k, name) in enumerate(('t0_str', 'pl_a_str', 'pl_rad_str')):
            label_map[2 * len(s.mcmc_dataset_slots) + 7 * problem.npl + 2 + 3 * j + k] = s.stellar_rho_str if k == 1 and s.rho else getattr(
                s,
                name,
            )[slot]
    labels = np.array([label_map.get(int(i), f'parameter[{i}]') for i in flags])
    err = np.column_stack((result.errors, result.errors))
    rv.return_results(s,result.parameters[flags],labels,result.parameters.copy(),flags,problem.npl,
        [0]*len(s.mcmc_dataset_slots),s.tra_data_sets,tr_model,s.tr_params,s.epoch,s.params.stellar_mass,
        problem.bounds,[],[],[],s.rtg,s.mixed_fit,err,'native',opt)
    s.parameters = result.parameters.copy()
    s.logprior = result.diagnostics['logprior']
    s.logposterior = result.diagnostics['logposterior']
    if result.diagnostics['objective'] == 'log posterior':
        print(f'Log prior: {s.logprior:.12g}; log posterior: {s.logposterior:.12g}')
    if not result.diagnostics['converged']:
        print('Native optimizer did not converge; local covariance errors are unavailable.')
        print('Covariance details:', result.diagnostics['covariance_reason'])
    elif result.diagnostics['covariance_status'] != 'available' and stages[-1] != 'evaluate':
        print('Native covariance unavailable:', result.diagnostics['covariance_reason'])
    print('lnL components (RV, LC, TTV, duration, depth):', result.components)
    print(f'--- Fortran joint fit: {result.diagnostics["elapsed_seconds"]:.5f} seconds; '
          f'lnL={result.loglik:.12g}; converged={result.diagnostics["converged"]} ---')
    return s
