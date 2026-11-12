"""Joint RV/activity GP, sharing the existing celerite kernel in time.

Dense reference backend: all series, measurement errors and cross-covariances
enter one normalized Gaussian likelihood. RV = G + b_RV G'; activity series
= offset + a G + b G'. G has the existing RV kernel's amplitude convention.
"""
import copy
import numpy as np
from scipy.linalg import cho_factor, cho_solve

FIELDS = ('response', 'derivative', 'offset', 'jitter')
DIFFERENTIABLE = ('SHOKernel', 'dSHOKernel', 'Matern32')


def ensure(session):
    """Add settings for loaded activity files while keeping existing choices."""

    if not hasattr(session, 'activity_gp'):
        session.activity_gp = {}
    if not hasattr(session, 'activity_gp_rv_derivative'):
        session.activity_gp_rv_derivative = dict(
            value=0.,
            fit=False,
            bounds=[-1e5, 1e5],
            normal=[0., 1., False],
            jeffreys=[.001, 1e5, False],
            error=[0., 0.],
        )
    for (slot, data) in getattr(session, 'act_data_sets', {}).items():
        if not len(data) or slot in session.activity_gp:
            continue
        y = np.asarray(data[1], float)
        rv = [np.asarray(row[1], float) for row in getattr(session, 'rv_data_sets', {}).values() if len(row)]
        scale = max(float(np.std(np.concatenate(rv))) if rv else 1., 1.)
        amplitude = max(float(np.std(y)), 1e-12) / scale
        settings = dict(enabled=False)
        for name, value, fit, bounds in [('response', amplitude, True, [-1e6, 1e6]),
                                        ('derivative', 0., False, [-1e6, 1e6]),
                                        ('offset', float(np.median(y)), True, [-1e12, 1e12]),
                                        ('jitter', 0., False, [0., 1e6])]:
            settings[name] = dict(
                value=value,
                fit=fit,
                bounds=bounds,
                normal=[value, max(abs(value), 1.), False],
                jeffreys=[1e-12, 1e6, False],
                error=[0., 0.],
            )
        session.activity_gp[slot] = settings
    return session.activity_gp


def selected(session):
    """Return the loaded activity slots selected for the joint GP."""

    return [slot for slot, setting in ensure(session).items()
            if setting['enabled'] and len(session.act_data_sets.get(slot, []))]


def validate(session):
    """Check selected data and the kernel support for derivative responses."""

    slots = selected(session)
    if not slots or not session.rtg[1]:
        return []
    if not session.type_fit.get('RV'):
        raise ValueError('Joint activity GP requires RV data and Use RV GP to be enabled.')
    if not any((len(row) for row in session.rv_data_sets.values())):
        raise ValueError('Joint activity GP requires at least one RV dataset.')
    settings = [session.activity_gp_rv_derivative] + [session.activity_gp[i]['derivative'] for i in slots]
    if session.gp_kernel not in DIFFERENTIABLE and any((s['fit'] or s['value'] != 0 for s in settings)):
        raise ValueError('GP derivative responses require SHO, double SHO, or Matérn 3/2. '
                         'Set derivative values to zero and uncheck Fit for the selected kernel.')
    for i in slots:
        data = np.asarray([session.act_data_sets[i][k] for k in range(3)], float)
        if data.ndim != 2 or data.shape[0] != 3 or (not data.shape[1]) or (not np.all(np.isfinite(data))) or np.any(data[2] <= 0):
            raise ValueError(f'Activity dataset {i + 1} needs finite times/values and positive uncertainties.')
    return slots


def pack(session, par, flags, labels, bounds, normal, jeffreys):
    """Insert before the fixed six-entry astrometry/stellar-mass tail."""
    slots = validate(session)
    session.activity_gp_indices = {}
    if not slots:
        return
    entries = [(None, 'rv_derivative', session.activity_gp_rv_derivative)]
    entries += [(slot, name, session.activity_gp[slot][name]) for slot in slots for name in FIELDS]
    for (slot, name, setting) in entries:
        session.activity_gp_indices[slot, name] = len(par)
        value = float(setting['value'])
        low, high = map(float, setting['bounds'])
        nr = list(setting['normal'])
        jp = list(setting['jeffreys'])
        if nr[2] and jp[2]:
            raise ValueError('Choose only one prior family for each activity GP parameter.')
        if jp[2]:
            if jp[0] <= 0 or jp[1] <= jp[0]:
                raise ValueError('Activity GP Jeffreys priors require 0 < lower < upper.')
            low, high = (max(low, jp[0]), min(high, jp[1]))
        if (not np.all(np.isfinite([value, low, high])) or low >= high or not low <= value <= high
                or (name == 'jitter' and low < 0) or (nr[2] and (not np.all(np.isfinite(nr[:2])) or nr[1] <= 0))):
            raise ValueError(f"Activity GP {(slot + 1 if slot is not None else 'RV')} {name}: invalid value, bounds, or prior.")
        par.append(value)
        flags.append(bool(setting['fit']))
        label = 'RV GP derivative [d]' if slot is None else f'Activity {slot + 1} GP {name}'
        labels.append(label)
        bounds.append([low, high])
        normal.append(nr)
        jeffreys.append(jp)


def configure(session, opt):
    """Copy the selected activity inputs for optimization or sampling."""

    slots = validate(session)
    opt['activity_gp'] = None
    if slots:
        opt['activity_gp'] = dict(
            indices=copy.deepcopy(session.activity_gp_indices),
            data={i: tuple((np.array(session.act_data_sets[i][k], float, copy=True) for k in range(3))) for i in slots},
            kernel=session.gp_kernel,
        )


def publish_parameters(session, par, errors, opt):
    """Store fitted activity responses, nuisance parameters and errors."""

    context = opt.get('activity_gp')
    if not context:
        return
    for ((slot, name), index) in context['indices'].items():
        setting = session.activity_gp_rv_derivative if slot is None else session.activity_gp[slot][name]
        setting['value'] = float(par[index])
        setting['error'] = list(map(float, errors[index]))


def kernel_matrices(kernel, lag):
    """k(t-t'), d k/dt and d² k/dt² from celerite's exponential terms."""
    ar, cr, ac, bc, cc, dc = kernel.coefficients
    r = np.abs(lag)
    k = np.zeros_like(r)
    first = k.copy()
    second = k.copy()
    for (a, c) in zip(ar, cr):
        e = a * np.exp(-c * r)
        k += e
        first -= c * e
        second += c * c * e
    for (a, b, c, d) in zip(ac, bc, cc, dc):
        e = np.exp(-c * r)
        co = np.cos(d * r)
        si = np.sin(d * r)
        k += e * (a * co + b * si)
        first += e * ((-c * a + d * b) * co + (-c * b - d * a) * si)
        second += e * (((c * c - d * d) * a - 2 * c * d * b) * co + ((c * c - d * d) * b + 2 * c * d * a) * si)
    return (k, np.sign(lag) * first, second)


def covariance(kernel, x, a, b, y=None, c=None, d=None):
    """Build the covariance between direct and derivative GP responses."""

    if y is None:
        y, c, d = (x, a, b)
    a, b, c, d = map(np.asarray, (a, b, c, d))
    lag = np.asarray(x)[:, None] - np.asarray(y)[None, :]
    if not np.any(b) and (not np.any(d)):
        return a[:, None] * c * kernel.get_value(lag)
    k, kp, kpp = kernel_matrices(kernel, lag)
    return a[:, None] * c * k + (b[:, None] * c - a[:, None] * d) * kp - b[:, None] * d * kpp


class JointGP:
    """Condition RV residuals and selected activity data on one time GP."""

    def __init__(self, kernel, time, residual, errors, context, par):
        self.kernel = kernel
        self.nrv = len(time)
        indices = context['indices']

        def parameter_value(slot, name):
            return float(par[indices[slot, name]])

        # RVs define the amplitude convention for the shared latent process.
        self.rv_b = parameter_value(None, 'rv_derivative')
        times = [np.asarray(time)]
        values = [np.asarray(residual)]
        noise = [np.asarray(errors) ** 2]
        responses = [np.ones(self.nrv)]
        derivative_responses = [np.full(self.nrv, self.rv_b)]
        self.slices = {}
        self.offsets = {}
        total = self.nrv

        # Append the selected indicators after the RV observations.
        for (slot, (t, y, err)) in context['data'].items():
            jitter = parameter_value(slot, 'jitter')
            if jitter < 0:
                raise ValueError('Activity GP jitter must be nonnegative.')
            self.slices[slot] = slice(total, total + len(t))
            total += len(t)
            self.offsets[slot] = parameter_value(slot, 'offset')
            from magnetic_cycle import mean as cycle_mean
            baseline = cycle_mean(context.get('magnetic_cycle'), par, t, slot)
            times.append(t)
            values.append(y - self.offsets[slot] - baseline)
            noise.append(err ** 2 + jitter ** 2)
            responses.append(np.full(len(t), parameter_value(slot, 'response')))
            derivative_responses.append(np.full(len(t), parameter_value(slot, 'derivative')))

        self.t, self.y, self.noise, self.a, self.b = map(np.concatenate, (times, values, noise, responses, derivative_responses))
        if not np.all(np.isfinite(np.r_[self.t, self.y, self.noise, self.a, self.b])) or np.any(self.noise <= 0):
            raise ValueError('Joint activity GP inputs must be finite with positive noise variance.')

        # Add measurement noise before the Cholesky likelihood calculation.
        K = covariance(kernel, self.t, self.a, self.b)
        K[np.diag_indices_from(K)] += self.noise
        if not np.all(np.isfinite(K)):
            raise ValueError('Joint activity GP covariance is not finite.')

        self.factor = cho_factor(K, lower=True, check_finite=False)
        self.alpha = cho_solve(self.factor, self.y, check_finite=False)
        self.chi2 = float(self.y @ self.alpha)
        self.loglik = float(
            -.5 * (self.chi2
                   + 2 * np.log(self.factor[0].diagonal()).sum()
                   + len(self.t) * np.log(2 * np.pi))
        )
        if not np.isfinite(self.loglik):
            raise ValueError('Joint activity GP likelihood is not finite.')

    def predict(self, time, response=1., derivative=None, *, return_var=True):
        """Predict a chosen response at the requested times."""

        t = np.asarray(time, float)
        a = np.full(len(t), response)
        b = np.full(len(t), self.rv_b if derivative is None else derivative)
        cross = covariance(self.kernel, t, a, b, self.t, self.a, self.b)
        mean = cross @ self.alpha
        if not return_var:
            return mean
        k0, _, kpp = kernel_matrices(self.kernel, np.zeros(1))
        prior = a * a * k0[0] - b * b * kpp[0]
        var = np.maximum(
            0.,
            prior - np.sum(cross * cho_solve(self.factor, cross.T, check_finite=False).T, axis=1),
        )
        return (mean, var)


def evaluate(kernel, time, residual, errors, context, par):
    """Return the joint log likelihood, rejecting invalid numerical trials."""

    try:
        return JointGP(kernel, time, residual, errors, context, par).loglik
    except (np.linalg.LinAlgError, ValueError, FloatingPointError, OverflowError):
        return -np.inf


def publish_model(session, get_lnl):
    """Use the same joint conditioning for RV data, RV curves and activity plots."""
    session.prepare_for_mcmc()
    opt = {}
    configure(session, opt)
    from magnetic_cycle import configure as configure_cycle, mean as cycle_mean
    configure_cycle(session, opt)
    result = session.fit_results
    time = np.asarray(result.rv_model.jd)
    y = np.asarray(result.rv_model.o_c)
    errors = np.array([np.hypot(err, session.jitt[int(slot)]) for (err, slot) in zip(result.rv_err, result.idset)])
    joint = JointGP(session.gps.kernel, time, y, errors, opt['activity_gp'], session.parameters)
    mean, var = joint.predict(time)
    session.gp_model_data = [mean, var, np.sqrt(var)]
    mean_grid, var_grid = joint.predict(result.model_jd)
    session.gp_model_curve = [mean_grid, var_grid, np.sqrt(var_grid)]
    session.activity_gp_models = {}
    for (slot, (t, _, _)) in opt['activity_gp']['data'].items():
        setting = session.activity_gp[slot]
        m, v = joint.predict(t, setting['response']['value'], setting['derivative']['value'])
        grid = np.linspace(float(np.min(t)), float(np.max(t)), max(200, min(1000, len(t))))
        curve = joint.predict(
            grid,
            setting['response']['value'],
            setting['derivative']['value'],
            return_var=False,
        )
        session.activity_gp_models[slot] = dict(
            time=t.copy(),
            mean=m + setting['offset']['value'] + cycle_mean(
                opt.get('magnetic_cycle'),
                session.parameters,
                t,
                slot,
            ),
            variance=v,
            curve_time=grid,
            curve_mean=curve + setting['offset']['value'] + cycle_mean(
                opt.get('magnetic_cycle'),
                session.parameters,
                grid,
                slot,
            ),
        )
    session.activity_gp_loglik = joint.loglik
    session.activity_gp_ndata = len(joint.y) - joint.nrv
    session.activity_gp_chi2 = joint.chi2
    if get_lnl:
        residual = y - mean
        result.chi2 = float(np.sum((residual / errors) ** 2))
        result.rms = float(np.sqrt(np.mean(residual ** 2)))
        result.wrms = float(np.sqrt(np.average(residual ** 2, weights=1 / errors ** 2)))
        dof = len(joint.y) - len(session.par_for_mcmc)
        result.reduced_chi2 = joint.chi2 / dof if dof > 0 else np.nan
        result.stat.dof = dof
        session.loglik = joint.loglik
