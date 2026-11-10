"""Opt-in Keplerian/N-body RV/transit model; no fitting or GUI state changes."""
from dataclasses import dataclass
import hashlib
import importlib.util
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import sysconfig
import tempfile
import threading

import numpy as np

_LOCK = threading.RLock()
_LIBRARY = None
_JACOBI = None
LD_LAWS = {'uniform': 0, 'linear': 1, 'quadratic': 2, 'nonlinear': 3,
           'squareroot': 4, 'logarithmic': 5, 'power2': 6}
LD_COUNTS = (0, 1, 2, 4, 2, 2, 2)


@dataclass(frozen=True)
class TransitRVModel:
    times: np.ndarray
    dataset: np.ndarray
    rv: np.ndarray
    flux: np.ndarray
    model_flux: np.ndarray
    jitter: np.ndarray
    scaled_axis: np.ndarray
    planet_radius: np.ndarray
    rv_components: np.ndarray
    mean_anomaly: np.ndarray
    transit_epoch: object


def _native():
    global _LIBRARY
    if _LIBRARY is not None:
        return _LIBRARY
    source = Path(__file__).resolve().parent / 'RV_mod' / 'rvmod'
    debug = os.environ.get('EXOSTRIKER_TRANSIT_DEBUG', '0') == '1'
    digest = hashlib.sha256((np.__version__ + str(debug) + 'joint-fit-api-1').encode())
    for path in [source / 'rvmod_for.f95', source / 'joint_fit.f95', source / 'shared_nbody.f95', *sorted(source.glob('*.inc'))]:
        digest.update(path.read_bytes())
    module_name = 'rvtransit_' + digest.hexdigest()[:16]
    cache = Path(os.environ.get('EXOSTRIKER_TRANSIT_CACHE_DIR',
                                str(Path.home() / '.cache/exostriker/transit-model')))
    cache.mkdir(parents=True, exist_ok=True)
    binary = cache / (module_name + sysconfig.get_config_var('EXT_SUFFIX'))
    if not binary.exists():
        if not shutil.which('gfortran') or not shutil.which('meson') or not shutil.which('ninja'):
            raise RuntimeError('Native transit models require gfortran, meson and ninja on PATH.')
        with tempfile.TemporaryDirectory(prefix='rv-transit-model-') as directory:
            directory = Path(directory)
            for path in [source / 'rvmod_for.f95', source / 'joint_fit.f95', source / 'shared_nbody.f95', *source.glob('*.inc')]:
                shutil.copy2(path, directory / path.name)
            env = os.environ.copy()
            options = '-O0 -g -fcheck=all -fbacktrace' if debug else '-O3'
            env['FFLAGS'] = (env.get('FFLAGS', '') + ' ' + options
                            + ' -std=legacy -fmax-stack-var-size=32768')
            include = '-I' + str(directory)
            env['FFLAGS'] += ' ' + (subprocess.list2cmdline([include])
                                    if sys.platform == 'win32' else shlex.quote(include))
            command = [sys.executable, '-m', 'numpy.f2py', '-c', '--backend', 'meson',
                       '-m', module_name, 'rvmod_for.f95', 'joint_fit.f95', 'shared_nbody.f95', 'only:', 'rv_transit_model', 'rv_transit_nbody', 'rv_transit_model_parts', 'rv_transit_nbody_parts', 'rv_transit_phased_model', 'rv_transit_geometry', 'rv_transit_event_flux', 'rv_transit_density', 'rv_transit_nbody_ttvfast_parts', 'rv_transit_geometry_ttvfast', 'rv_transit_likelihood', 'ttv_likelihood', 'native_joint_fit', 'native_joint_loglik', 'joint_timing_likelihood', 'joint_shared_nbody', ':']
            result = subprocess.run(command, cwd=directory, env=env, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, text=True, timeout=180)
            if result.returncode:
                (cache / 'build-error.log').write_text(result.stdout)
                raise RuntimeError('Cannot build native transit model:\n' + result.stdout[-6000:])
            # Publish atomically, including when separate GUI workers build at once.
            with tempfile.NamedTemporaryFile(dir=cache, delete=False) as target:
                staged = Path(target.name)
                with (directory / binary.name).open('rb') as compiled:
                    shutil.copyfileobj(compiled, target)
            try:
                staged.replace(binary)
            except PermissionError:
                # Windows can keep a concurrently loaded identical DLL open.
                staged.unlink(missing_ok=True)
                if not binary.exists():
                    raise
    spec = importlib.util.spec_from_file_location(module_name, binary)
    library = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(library)
    _LIBRARY = library
    return library


def evaluate(planets, times, *, stellar_radius=1., scaled_axis=None, radius_ratio,
             epoch=0., dataset=0, offsets=None, jitters=None, limb_dark=None,
             coefficients=None, exposure_time=0., supersample_factor=1,
             quadrature_order=128, model='keplerian', stellar_mass=1.,
             step=None, tolerance=1.e-12, max_steps=1000000, phase='M0', t0=None, rho_star=None, jacobi_convention='standard', _observations=None, _shared_nbody=False):
    """Return RV [m/s], normalized flux, flux-offset and per-sample jitter.

    planets rows: K [m/s], P [days], e, omega [deg], M(epoch) [deg], i [deg],
    node [deg]. Dataset labels are zero-based GUI slots 0..59, never renumbered.
    offsets/jitters are length-60 arrays; coefficients have shape (60,4).
    stellar_radius is in solar radii; scaled_axis = a/Rstar and radius_ratio =
    Rp/Rstar. Their physical radii are returned in solar radii. model='nbody'
    uses RV-derived masses and physical Jacobi axes with one stellar radius;
    scaled_axis can be omitted and, if supplied, must match those derived axes.
    For Keplerian models phase='t0' takes per-planet t0 and derives M0 in Fortran
    for both RV and flux. phase='M0' returns the next nominal conjunction epoch.
    N-body models require M0 and return no analytic transit_epoch.
    rho_star [kg/m^3] replaces scaled_axis; Fortran derives all Keplerian axes
    or one physical N-body stellar radius from stellar_mass and density.
    """
    planets = np.asarray(planets, dtype=float)
    if planets.size == 0:
        planets = planets.reshape(0, 7)
    if planets.ndim != 2 or planets.shape[1] != 7 or len(planets) > 9:
        raise ValueError('Expected up to 9 rows of [K,P,e,omega,M,i,node].')
    times = np.asarray(times, dtype=float)
    if times.ndim != 1:
        raise ValueError('times must be a one-dimensional array.')
    labels = np.broadcast_to(np.asarray(dataset, dtype=float), times.shape)
    if not np.all(np.isfinite(labels)) or np.any(labels != np.floor(labels)):
        raise ValueError('Dataset slots must be integers in 0..59.')
    if np.any(labels < 0) or np.any(labels >= 60):
        raise ValueError('Dataset slots must be integers in 0..59.')
    labels = labels.astype(np.int32)
    if jacobi_convention not in ('standard', 'ttvfast'):
        raise ValueError("jacobi_convention must be 'standard' or 'ttvfast'.")
    if model not in ('keplerian', 'nbody'):
        raise ValueError("model must be 'keplerian' or 'nbody'.")
    if phase not in ('M0', 't0'):
        raise ValueError("phase must be 'M0' or 't0'.")
    if model == 'nbody' and phase != 'M0':
        raise ValueError("N-body models require phase='M0' at the integration epoch; transit times follow the trajectory.")
    if phase == 't0' and t0 is None:
        raise ValueError("Supply t0 for phase='t0'.")
    transit_input = np.zeros(len(planets)) if t0 is None else np.broadcast_to(
        np.asarray(t0, dtype=float), (len(planets),)).copy()
    if rho_star is not None:
        if scaled_axis is not None:
            raise ValueError('Supply rho_star or scaled_axis, not both.')
        density_axis, stellar_radius = density_geometry(planets[:, 1], rho_star, stellar_mass)
        if model == 'keplerian':
            scaled_axis = density_axis
    if model == 'nbody':
        if not np.isfinite(stellar_radius) or stellar_radius <= 0:
            raise ValueError('stellar_radius must be finite and positive.')
        masses, axes_au = _mass_axes(planets, stellar_mass)
        axis = axes_au * jacobi_axis_factor(planets, masses, stellar_mass, jacobi_convention) * 1.49597892e11 / (stellar_radius * 6.957e8)
        if scaled_axis is not None:
            supplied = np.broadcast_to(np.asarray(scaled_axis, dtype=float), axis.shape)
            if not np.allclose(supplied, axis, rtol=1.e-8, atol=1.e-10):
                raise ValueError('N-body scaled_axis must match RV-derived axes and stellar_radius; omit it to derive them.')
    else:
        if scaled_axis is None:
            raise ValueError('scaled_axis is required for the Keplerian model.')
        axis = np.broadcast_to(np.asarray(scaled_axis, dtype=float), (len(planets),)).copy()
    radii = np.broadcast_to(np.asarray(radius_ratio, dtype=float), (len(planets),))
    rows = np.column_stack((planets, axis, radii))
    sets = np.zeros((60, 6))
    for value, col in ((offsets, 0), (jitters, 1)):
        if value is not None:
            if np.shape(value) != (60,):
                raise ValueError('offsets and jitters must have exactly 60 entries.')
            sets[:, col] = value
    laws = np.zeros(60, dtype=np.int32)
    if limb_dark is not None:
        if len(limb_dark) != 60 or any(law not in LD_LAWS for law in limb_dark):
            raise ValueError('limb_dark must contain 60 supported law names.')
        laws[:] = [LD_LAWS[law] for law in limb_dark]
    if coefficients is not None:
        if np.shape(coefficients) != (60, 4):
            raise ValueError('coefficients must have shape (60,4), zero padded.')
        sets[:, 2:] = coefficients
    exposure = np.broadcast_to(np.asarray(exposure_time, dtype=float), (60,)).copy()
    samples = np.broadcast_to(np.asarray(supersample_factor, dtype=float), (60,))
    if (not np.all(np.isfinite(samples)) or np.any(samples != np.floor(samples))
            or np.any(samples < 1) or np.any(samples > 1001)):
        raise ValueError('supersample_factor must contain integers in 1..1001.')
    if not np.isscalar(quadrature_order) or not np.isfinite(quadrature_order):
        raise ValueError('quadrature_order must be an integer in 16..1024.')
    if quadrature_order != int(quadrature_order) or not 16 <= quadrature_order <= 1024:
        raise ValueError('quadrature_order must be an integer in 16..1024.')
    if model == 'nbody':
        if step is None:
            step = np.min(planets[:, 1])/50 if len(planets) else 1.
        if not np.isfinite(step) or step <= 0 or not np.isfinite(tolerance) or not 1e-15 <= tolerance <= 1e-3:
            raise ValueError('step must be positive and tolerance must be in 1e-15..1e-3.')
        if not np.isfinite(max_steps) or max_steps != int(max_steps) or not 1 <= max_steps <= np.iinfo(np.int32).max:
            raise ValueError('max_steps must be a positive 32-bit integer.')
        if not np.all(np.isfinite(exposure)) or np.any(exposure < 0):
            raise ValueError('Exposure durations must be finite and nonnegative.')
    if _shared_nbody:
        if model!='nbody' or _observations is not None:
            raise ValueError('Shared dense output requires an N-body model-only request.')
        return _run_shared_model(planets,masses,radii,times,labels,axis,stellar_mass,stellar_radius,
                                 epoch,laws,sets,exposure,samples.astype(np.int32),step,tolerance,int(max_steps),jacobi_convention)
    if _observations is not None:
        from native_likelihood import _evaluate_prepared
        return _evaluate_prepared(planets, times, labels, rows, axis,
            masses if model == 'nbody' else np.zeros(len(planets)), stellar_mass,
            stellar_radius, epoch, laws, sets, exposure, samples.astype(np.int32),
            int(quadrature_order), model, step, tolerance, int(max_steps),
            phase, transit_input, jacobi_convention, _observations)
    with _LOCK:
        # Some f2py versions cannot wrap zero-sized observation buffers. Still
        # validate all native parameters, then discard the dummy prediction.
        native_times = times if times.size else np.array([float(epoch)])
        native_labels = labels+1 if times.size else np.array([1], dtype=np.int32)
        if model == 'keplerian' and not len(planets):
            # Some f2py versions reject the empty phase-vector buffer.
            mean_anomaly = np.empty(0)
            transit_epoch = np.empty(0)
            rv, flux, model_flux, jitter, planet_radius, components, status = _native().rv_transit_model_parts(
                native_times, native_labels, rows, float(stellar_radius), float(epoch), laws, sets,
                exposure, samples.astype(np.int32), int(quadrature_order))
        elif model == 'keplerian':
            (rv, flux, model_flux, jitter, planet_radius, components,
             mean_anomaly, transit_epoch, status) = _native().rv_transit_phased_model(
                native_times, native_labels, rows, float(stellar_radius), float(epoch), laws, sets,
                exposure, samples.astype(np.int32), int(quadrature_order), int(phase == 't0'), transit_input)
        else:
            mean_anomaly = np.mod(planets[:, 4], 360.)
            transit_epoch = None  # N-body transit epochs must be measured on the trajectory.
            rv, flux, model_flux, jitter, planet_radius, components, status = _run_nbody(
                planets, masses, radii, native_times, native_labels, stellar_mass, stellar_radius,
                epoch, laws, sets, exposure, samples.astype(np.int32), int(quadrature_order),
                step, tolerance, int(max_steps), jacobi_convention)
    if status:
        errors = {1: 'Invalid transit model input/profile or orbit intersects the star.',
                  2: 'Transit model Kepler solver did not converge.',
                  4: 'N-body integration exceeded max_steps.',
                  5: 'N-body integration failed to advance or produced nonfinite output.'}
        raise ValueError(errors.get(status, 'Native transit model failed: ' + str(status)))
    return TransitRVModel(times.copy(), labels.copy(), rv[:len(times)], flux[:len(times)],
                          model_flux[:len(times)], jitter[:len(times)], axis, planet_radius, components[:len(times)],
                          mean_anomaly, transit_epoch)


def density_geometry(period, rho_star, stellar_mass=1.):
    """Native shared-density conversion: kg/m^3, days, solar stellar mass.

    Keplerian a/Rstar uses the GUI negligible-planet-mass convention. N-body
    uses the returned physical radius and full Jacobi masses instead.
    """
    period = np.atleast_1d(np.asarray(period, float))
    if period.size == 0:
        # Avoid zero-sized f2py input buffers.
        axis, radius = density_geometry([1.], rho_star, stellar_mass)
        return axis[:0], radius
    if np.ndim(rho_star) or np.ndim(stellar_mass):
        raise ValueError('rho_star and stellar_mass must be shared scalars.')
    with _LOCK:
        axis, radius, status = _native().rv_transit_density(period, float(rho_star), float(stellar_mass))
    if status:
        raise ValueError('Stellar density, stellar mass and periods must be finite and positive.')
    return axis, radius


def _mass_axes(planets, stellar_mass):
    """Load the existing RV conversion without importing the GUI/session package."""
    global _JACOBI
    if _JACOBI is None:
        path = Path(__file__).resolve().parent / 'RV_mod' / 'jacobi.py'
        spec = importlib.util.spec_from_file_location('_exostriker_transit_jacobi', path)
        _JACOBI = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(_JACOBI)
    if np.any(planets[:, 0] < 0) or np.any(planets[:, 5] < 0) or np.any(planets[:, 5] > 180):
        raise ValueError('Require nonnegative K and inclinations in 0..180 degrees.')
    masses, axes = _JACOBI.mass_a_from_jacobi_el(planets, len(planets), float(stellar_mass))
    return masses/1047.5654817267318, axes


def jacobi_axis_factor(planets, masses, stellar_mass, convention):
    """TTVFast vs standard initial axis ratio, retaining original slot order."""
    factor = np.ones(len(planets))
    if convention == 'ttvfast':
        order = np.argsort(planets[:, 1], kind='stable')
        inner = stellar_mass + np.r_[0., np.cumsum(masses[order])[:-1]]
        if len(order):
            factor[order] = (stellar_mass/inner)**(1/3)
    return factor


def session_jacobi_convention(session):
    """Old transit/TTV likelihoods use TTVFast; RV-only uses standard RVmod."""
    convention = getattr(session, 'native_jacobi_convention', 'auto')
    if convention == 'auto':
        kinds = getattr(session, 'type_fit', {})
        convention = 'ttvfast' if kinds.get('Transit', False) or kinds.get('TTV', False) else 'standard'
    return convention


def _run_nbody(planets, masses, radii, times, labels, stellar_mass, stellar_radius,
               epoch, laws, sets, exposure, samples, nquad, step, tolerance, max_steps, jacobi_convention):
    if not len(planets):
        # Empty element/mass buffers are rejected by some f2py wrappers. The
        # zero-planet native model has exactly the same baseline in both modes.
        return _native().rv_transit_model_parts(
            times, labels, np.empty((0, 9)), float(stellar_radius), float(epoch),
            laws, sets, exposure, samples, nquad)
    # Expand exposure samples in relative time before sorting, avoiding loss of
    # subsecond precision at large absolute epochs. Restore caller ordering later.
    counts = samples[labels-1]
    owners = np.repeat(np.arange(len(times)), counts)
    starts = np.repeat(np.cumsum(counts)-counts, counts)
    sample_index = np.arange(len(owners))-starts
    fractions = np.zeros(len(owners))
    mask = counts[owners] > 1
    fractions[mask] = sample_index[mask]/(counts[owners[mask]]-1) - .5
    expanded = (times-epoch)[owners] + exposure[labels[owners]-1]*fractions
    order = np.argsort(expanded, kind='stable')
    planet_order = np.argsort(planets[:, 1], kind='stable')
    backend = (_native().rv_transit_nbody_ttvfast_parts if jacobi_convention == 'ttvfast'
               else _native().rv_transit_nbody_parts)
    rv, flux, model_flux, jitter, physical_radii, components, status = backend(
        expanded[order], labels[owners[order]], planets[planet_order, 1:], masses[planet_order],
        radii[planet_order], float(stellar_mass), float(stellar_radius), laws, sets,
        nquad, float(step), float(tolerance), max_steps)
    # Accumulate each exposure in its original dataset/observation slot.
    def average(values):
        return np.bincount(owners[order], weights=values, minlength=len(times))/counts
    averaged_parts = np.zeros((len(times), len(planets)))
    np.add.at(averaged_parts, owners[order], components)
    averaged_parts /= counts[:, None]
    return (average(rv), average(flux), average(model_flux), average(jitter),
            physical_radii[np.argsort(planet_order)],
            averaged_parts[:, np.argsort(planet_order)], status)


def from_session(session, times=None, *, dataset=0, radius_mode=None, model=None, phase='M0', **options):
    """Read GUI parameters without modifying session state or existing models.

    model defaults to the session's dynamical setting. radius_mode defaults to
    'density' for GUI rho mode, otherwise 'stellar' for N-body and 'scaled'
    for Keplerian. 'stellar' derives a/Rstar from the
    RV-derived Jacobi axes and stellar_radius. phase='t0' uses GUI t0 for
    Keplerian predictions and derives M0 inside Fortran for both RV and flux.
    Omitting times evaluates only the selected dataset's observation times.
    """
    if times is None:
        if not np.isscalar(dataset) or not np.isfinite(dataset) or dataset != int(dataset) or not 0 <= dataset < 60:
            raise ValueError('Select one dataset slot in 0..59 when times is omitted.')
        data = session.tra_data_sets[int(dataset)]
        if not len(data):
            raise ValueError('The selected transit dataset is empty.')
        times = data[0]
    if session.gr_flag:
        raise ValueError('GR is not supported by this Newtonian transit evaluator.')
    slots = list(np.flatnonzero(session.use_planet))
    if model is None:
        model = 'nbody' if session.mod_dynamical else 'keplerian'
    if model not in ('keplerian', 'nbody'):
        raise ValueError("model must be 'keplerian' or 'nbody'.")
    if model == 'nbody':
        mixed = getattr(session, 'mixed_fit', {0: False})
        if bool(np.asarray(mixed[0]).item()):
            raise ValueError('Mixed dynamical subsets are not supported by this transit evaluator.')
    if radius_mode is None:
        radius_mode = 'density' if getattr(session, 'rho', False) or options.get('rho_star') is not None else ('stellar' if model == 'nbody' else 'scaled')
    if any(session.omega_dot[i] != 0 for i in slots):
        raise ValueError('Precession is not supported by this transit evaluator.')
    rows = np.array([[getattr(session, name)[i] for name in
                      ('K', 'P', 'e', 'w', 'M0', 'i', 'Node')] for i in slots]).reshape(-1, 7)
    if getattr(session, 'hkl', False) and len(rows):
        h = np.array([session.e_sinw[i] for i in slots], dtype=float)
        k = np.array([session.e_cosw[i] for i in slots], dtype=float)
        omega = np.degrees(np.arctan2(h, k))
        rows[:, 2] = np.hypot(h, k)
        rows[:, 3] = omega
        rows[:, 4] = (np.array([session.lamb[i] for i in slots])-omega) % 360.
    if model == 'nbody' and getattr(session, 'copl_incl', False) and len(rows):
        rows[:, 5] = rows[0, 5]
    stellar_radius = options.get('stellar_radius', session.stellar_radius)
    stellar_mass = options.get('stellar_mass', session.params.stellar_mass)
    if radius_mode != 'density' and (not np.isfinite(stellar_radius) or stellar_radius <= 0):
        raise ValueError('stellar_radius must be finite and positive.')
    if radius_mode == 'density':
        axis = None
        # Legacy transit fitting stores the shared density in the first active
        # planet's pl_a entry; stellar_rho is a separate stellar-properties field.
        options.setdefault('rho_star', session.pl_a[slots[0]] if slots else session.stellar_rho)
    elif radius_mode == 'scaled':
        if getattr(session, 'rho', False):
            raise ValueError("GUI pl_a contains density in rho mode; use radius_mode='stellar' for this evaluator.")
        axis = [session.pl_a[i] for i in slots]
    elif radius_mode == 'stellar':
        # Same conversion already used by the RV session; no independent t0/aR fit.
        masses, axes_au = _mass_axes(rows, stellar_mass)
        convention = options.get('jacobi_convention', session_jacobi_convention(session))
        if model == 'nbody':
            axes_au *= jacobi_axis_factor(rows, masses, stellar_mass, convention)
        axis = axes_au * 1.49597892e11 / (stellar_radius * 6.957e8)
    else:
        raise ValueError("radius_mode must be 'scaled', 'stellar' or 'density'.")
    coeff = np.zeros((60, 4))
    tables = {'linear': session.ld_u_lin, 'quadratic': session.ld_u_quad,
              'nonlinear': session.ld_u_nonlin}
    for d, law in enumerate(session.ld_m):
        if law == 'uniform':
            continue
        if law not in tables:
            raise ValueError('Unsupported GUI limb-darkening law: ' + str(law))
        values = tables[law][d]
        coeff[d, :LD_COUNTS[LD_LAWS[law]]] = values
    from model_epoch import resolve_epoch
    epoch=resolve_epoch(session)
    defaults = dict(model=model, phase=phase, stellar_mass=stellar_mass,
                    stellar_radius=stellar_radius, scaled_axis=axis,
                    radius_ratio=[session.pl_rad[i] for i in slots], epoch=epoch,
                    offsets=[session.tra_off[d] for d in range(60)],
                    jitters=[session.tra_jitt[d] for d in range(60)],
                    limb_dark=list(session.ld_m), coefficients=coeff)
    if model == 'nbody':
        defaults['jacobi_convention'] = session_jacobi_convention(session)
    if phase == 't0':
        defaults['t0'] = [session.t0[i] for i in slots]
    defaults.update(options)
    evaluator = defaults.pop('_evaluator', evaluate)
    return evaluator(rows, times, dataset=dataset, **defaults)


def _run_shared_model(planets,masses,radii,times,labels,axis,stellar_mass,stellar_radius,
                      epoch,laws,sets,exposure,samples,step,tolerance,max_steps,convention):
    """Internal final-plot adapter using the same dense trajectory as the fit."""
    original_size=len(times)
    if not original_size:times=np.array([epoch]);labels=np.array([0],np.int32)
    counts=samples[labels]
    owners=np.repeat(np.arange(len(times)),counts)
    index=np.arange(len(owners))-np.repeat(np.cumsum(counts)-counts,counts)
    fraction=np.zeros(len(owners));mask=counts[owners]>1
    fraction[mask]=index[mask]/(counts[owners[mask]]-1)-.5
    expanded=(times-epoch)[owners]+exposure[labels[owners]]*fraction
    ordered=np.argsort(expanded,kind='stable');bodies=np.argsort(planets[:,1],kind='stable')
    with _LOCK:
        rv,flux,event_count,event_table,steps,parts,status=_native().joint_shared_nbody(
            expanded[ordered],np.full(len(owners),2,np.int32),labels[owners[ordered]]+1,
            planets[bodies,1:],masses[bodies],radii[bodies],stellar_mass,stellar_radius,
            epoch,epoch,epoch,laws,sets,step,min(tolerance,1e-14),max_steps,int(convention=='ttvfast'),1)
    if status:raise ValueError(f'Shared plotting trajectory failed (status {status}).')
    avg_rv=np.zeros(len(times));avg_flux=np.zeros(len(times));avg_parts=np.zeros((len(times),len(planets)))
    np.add.at(avg_rv,owners[ordered],rv);np.add.at(avg_flux,owners[ordered],flux)
    np.add.at(avg_parts,owners[ordered],parts[:,np.argsort(bodies)])
    avg_rv/=counts;avg_flux/=counts;avg_parts/=counts[:,None]
    size=original_size
    return TransitRVModel(times[:size].copy(),labels[:size].copy(),avg_rv[:size],avg_flux[:size],
        avg_flux[:size]-sets[labels[:size],0],sets[labels[:size],1],axis,radii*stellar_radius,
        avg_parts[:size],np.mod(planets[:,4],360.),None)
