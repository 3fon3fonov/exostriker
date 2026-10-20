"""Shared plotting and posterior exports for examples provided by Trifon Trifonov with the aid of CODEX; use with caution."""
from pathlib import Path
import csv
import os
import sys
import warnings

# Resolve paths before importing Exo-Striker (older versions change directory).
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'exostriker' / 'datafiles'
EXAMPLES = ROOT / 'examples'
sys.path.insert(0, str(ROOT))
_START_DIRECTORY = Path.cwd()
import exostriker.lib.RV_mod as rv
os.chdir(_START_DIRECTORY)
import dill
import numpy as np
import matplotlib.pyplot as plt


def uniform_bound(fit, key, index, limits, plural=False):
    """Set a finite uniform prior, explicitly disabling Gaussian/log priors."""
    getattr(fit, key + ('_bounds' if plural else '_bound'))[index] = np.array(limits, dtype=float)
    for suffix in ('_norm_pr', '_jeff_pr'):
        prior = np.array(getattr(fit, key + suffix)[index], copy=True)
        prior[2] = False
        getattr(fit, key + suffix)[index] = prior


def check_priors(fit):
    fit.prepare_for_mcmc(rtg=fit.rtg)
    bounds = np.asarray(fit.b_for_mcmc, dtype=float)
    start = np.asarray(fit.par_for_mcmc, dtype=float)
    if not len(start) or not np.isfinite(bounds).all() or np.any(bounds[:, 0] >= bounds[:, 1]):
        raise ValueError('Every free parameter needs finite, increasing bounds.')
    if np.any(start < bounds[:, 0]) or np.any(start > bounds[:, 1]):
        raise ValueError('A starting value lies outside its prior bounds; inspect the fit.')
    print('\nFree parameters and uniform prior bounds:')
    for label, value, (low, high) in zip(fit.e_for_mcmc, start, bounds):
        print(f'{label:30s} {value:14.6g}  [{low:g}, {high:g}]')
    return bounds


def fit_rvs(fit, dynamical=False, starts=10):
    fit.update_mod_dynamical(dynamical)
    fit.fitting(minimize_fortran=True, minimize_loglik=True, outputfiles=[1, 1, 1],
                amoeba_starts=starts, fortran_kill=30, timeout_sec=600,
                eps=1000, dt=2, npoints=2000, model_max=0, model_min=0)
    if fit.flag != 1 or not np.isfinite(fit.loglik):
        raise RuntimeError('RV fit failed; inspect the starting parameters.')
    # The native likelihood uses jitter squared; report its positive branch.
    for slot in range(fit.ndset):
        fit.jitt[slot] = abs(fit.jitt[slot])
    print(f'{"Dynamical" if dynamical else "Keplerian"} lnL = {fit.loglik:.4f}')
    for i in range(fit.npl):
        print(f'Planet {i+1}: P={fit.P[i]:.6f} d, K={fit.K[i]:.4f} m/s, e={fit.e[i]:.5f}')
    return fit


def save_session(fit, directory, target):
    directory = Path(directory).expanduser().resolve()
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f'{target}.ses'
    with path.open('wb') as stream:
        dill.dump(fit, stream)
    print('Session:', path)
    return path


def plot_rvs(fit, directory, target):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    model = fit.fit_results.rv_model
    times = np.asarray(model.jd)
    residual = np.asarray(model.o_c)
    ids = np.asarray(fit.fit_results.idset, dtype=int)
    # Reconstruct offset-corrected RVs from the fitted planetary model + residual.
    values = np.asarray(fit.fit_results.model_data[5]) + residual
    errors = np.sqrt(np.asarray(model.rv_err)**2 + np.array([fit.jitt[i]**2 for i in ids]))
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(10, 6),
                             gridspec_kw={'height_ratios': [3, 1]}, layout='constrained')
    axes[0].plot(fit.fit_results.model_jd, fit.fit_results.model, 'k-', lw=1)
    for index in np.unique(ids):
        chosen = ids == index
        for ax, y in zip(axes, (values, residual)):
            ax.errorbar(times[chosen], y[chosen], yerr=errors[chosen], fmt='.',
                        label=f'Dataset {index+1}', alpha=.8)
    axes[0].legend(); axes[0].set_ylabel('RV [m/s]'); axes[0].set_title(fit.name)
    axes[1].axhline(0, color='k', lw=.8)
    axes[1].set(xlabel='JD [day]', ylabel='O − C [m/s]')
    fig.savefig(directory / f'{target}_rv.png', dpi=160)
    figures = [fig]
    for index in range(fit.npl):
        data, curve = rv.phase_RV_planet_signal(fit, index)  # zero-based planet slot
        fig, ax = plt.subplots(figsize=(8, 4), layout='constrained')
        ax.plot(curve[0], curve[1], 'k-', lw=1)
        dataset_ids = np.asarray(data[3], dtype=int)
        for slot in np.unique(dataset_ids):
            chosen = dataset_ids == slot
            err = np.sqrt(np.asarray(data[2])[chosen]**2 + fit.jitt[slot]**2)
            ax.errorbar(np.asarray(data[0])[chosen], np.asarray(data[1])[chosen],
                        yerr=err, fmt='.', label=f'Dataset {slot+1}')
        ax.set(xlabel='Orbital phase [day]', ylabel='RV [m/s]', title=f'{fit.name}, planet {index+1}')
        ax.legend(); fig.savefig(directory / f'{target}_phase_{index+1}.png', dpi=160)
        figures.append(fig)
    return figures


def plot_transit(fit, directory, target):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    times, flux, errors, model, residual = [np.asarray(v, dtype=float) for v in fit.transit_results[2][:5]]
    rich_time, rich_flux = fit.transit_results[3]
    origin = int(times.min())
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(9, 6),
                             gridspec_kw={'height_ratios': [3, 1]}, layout='constrained')
    axes[0].errorbar(times-origin, flux, yerr=errors, fmt='.', alpha=.5)
    axes[0].plot(np.asarray(rich_time)-origin, rich_flux, 'k-', lw=1.5)
    axes[0].set(ylabel='Relative flux', title=fit.name)
    axes[1].errorbar(times-origin, residual, yerr=errors, fmt='.', alpha=.5)
    axes[1].axhline(0, color='k', lw=.8)
    axes[1].set(xlabel=f'BJD − {origin} [day]', ylabel='O − C')
    fig.savefig(directory / f'{target}_transit.png', dpi=160)
    phase = ((times-fit.t0[0] + fit.P[0]/2) % fit.P[0]) - fit.P[0]/2
    phase_figure, phase_axes = plt.subplots(2, 1, sharex=True, figsize=(8, 6),
        gridspec_kw={'height_ratios': [3, 1]}, layout='constrained')
    phase_axes[0].errorbar(phase, flux, yerr=errors, fmt='.', alpha=.2)
    order = np.argsort(phase)
    phase_axes[0].plot(phase[order], model[order], 'k-', lw=1.5)
    phase_axes[0].set(ylabel='Relative flux', title=fit.name)
    phase_axes[1].errorbar(phase, residual, yerr=errors, fmt='.', alpha=.2)
    phase_axes[1].axhline(0, color='k', lw=.8)
    phase_axes[1].set(xlabel='Time from mid-transit [day]', ylabel='O − C', xlim=(-.2, .2))
    phase_figure.savefig(directory / f'{target}_transit_phase.png', dpi=160)
    return [fig, phase_figure]


def configure_sampling(fit, threads=1):
    fit.Dynamic_nest = True
    fit.ns_samp_method = 'rwalk'
    fit.ns_samp_bound = 'multi'
    fit.live_points_fact = 10
    fit.stop_crit = .1
    fit.ns_pfrac = 1.0
    fit.ns_use_stop = True
    fit.ns_threads = threads
    fit.ns_progress = True
    fit.ns_maxiter = {0: False, 1: 1_000_000}
    fit.ns_maxcall = {0: False, 1: 10_000_000}
    fit.nest_weighted = True
    fit.ns_save_sampler = True
    fit.ns_fileoutput = False
    fit.ns_save_maxlnL = True
    fit.ns_save_median = fit.ns_save_means = fit.ns_save_mode = False
    fit.nest_percentile_level = 68.3


def posterior_arrays(result, seed=42):
    """Resample raw rows and their lnL together, retaining raw weights for diagnostics."""
    raw = np.asarray(result.samples_ns, dtype=float)
    logl = np.asarray(result.logl_ns, dtype=float)
    logwt = np.asarray(result.logwt, dtype=float)
    if raw.ndim != 2 or not len(raw) or len(raw) != len(logl) or len(raw) != len(logwt):
        raise ValueError('Raw nested-sampling arrays are missing or misaligned.')
    weights = np.exp(logwt - np.max(logwt))
    weights /= weights.sum()
    positions = (np.random.default_rng(seed).random() + np.arange(len(raw))) / len(raw)
    cumulative = np.cumsum(weights); cumulative[-1] = 1.
    chosen = np.searchsorted(cumulative, positions)
    return raw, logl, weights, raw[chosen], logl[chosen]


def export_posterior(fit, directory, target):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    result = fit.ns_sampler
    raw, logl, weights, posterior, post_logl = posterior_arrays(result)
    labels = list(result.lbf['OrigLabels'])
    if len(labels) != raw.shape[1]:
        raise ValueError('Parameter labels do not match the posterior width.')
    # Include a unique index: duplicate human-readable labels remain unambiguous.
    columns = [f'p{i+1}: {label}' for i, label in enumerate(labels)]
    def write_csv(path, header, rows):
        with path.open('w', newline='') as stream:
            writer = csv.writer(stream); writer.writerow(header); writer.writerows(rows)
    write_csv(directory / f'{target}_posterior.csv', columns + ['logL'],
              np.column_stack((posterior, post_logl)))
    write_csv(directory / f'{target}_posterior_raw.csv', columns + ['logL', 'logwt', 'weight'],
              np.column_stack((raw, logl, result.logwt, weights)))
    tail = (100 - fit.nest_percentile_level) / 2
    quantiles = np.percentile(posterior, [tail, 50, 100-tail], axis=0)
    best = raw[np.argmax(logl)]
    ess = 1 / np.sum(weights**2)
    with (directory / f'{target}_posterior_summary.csv').open('w', newline='') as stream:
        stream.write(f'# Examples provided by Trifon Trifonov with the aid of CODEX; use with caution.\n# ESS = {ess:.3f}\n')
        stream.write(f'# lnZ = {result.logz[-1]:.6f} +/- {result.logzerr[-1]:.6f}\n')
        stream.write(f'# credible interval = {fit.nest_percentile_level}%\n')
        writer = csv.writer(stream)
        writer.writerow(['parameter', 'label', 'median', 'err_minus', 'err_plus', 'best_lnL'])
        for i, label in enumerate(labels):
            low, median, high = quantiles[:, i]
            writer.writerow([f'p{i+1}', label, median, median-low, high-median, best[i]])
    print(f'lnZ = {result.logz[-1]:.4f} ± {result.logzerr[-1]:.4f}; ESS = {ess:.1f}')
    print(f'niter = {result.niter}; ncall = {np.sum(result.ncall)}; efficiency = {result.eff:.2f}%')
    if ess < 200:
        warnings.warn('Low effective sample size: do not interpret these intervals as converged.')
    if np.any(np.ptp(posterior, axis=0) == 0):
        warnings.warn('Posterior has unresolved columns; skipping the corner plot. Run longer.')
        return posterior
    # Use the project's bundled corner implementation; no external LaTeX required.
    import exostriker.lib.corner_ES as corner
    figure = corner.corner(posterior, labels=labels, quantiles=[.16, .5, .84], show_titles=True)
    figure.savefig(directory / f'{target}_corner.png', dpi=160, bbox_inches='tight')
    figure.savefig(directory / f'{target}_corner.pdf', bbox_inches='tight')
    return posterior


def finish(fit, output, target, nested=False, threads=1):
    if nested:
        check_priors(fit)
        configure_sampling(fit, threads)
        fit = rv.run_nestsamp(fit)
        if fit.bound_error:
            raise RuntimeError(fit.bound_error_msg)
        export_posterior(fit, output, target)
    save_session(fit, output, target)
    return fit


def command_line():
    import argparse
    parser = argparse.ArgumentParser(description='Exo-Striker example provided by Trifon Trifonov with the aid of CODEX; use with caution.')
    parser.add_argument('--nested', action='store_true', help='Also run nested sampling (can take hours).')
    parser.add_argument('--threads', type=int, default=1, help='Nested-sampling processes; default: 1.')
    parser.add_argument('--output', type=Path, help='Output directory; default: examples/results/TARGET.')
    parser.add_argument('--show', action='store_true', help='Show plot windows after saving.')
    args = parser.parse_args()
    if args.threads < 1:
        parser.error('--threads must be positive')
    return args
