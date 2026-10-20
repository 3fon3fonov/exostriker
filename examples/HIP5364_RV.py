#!/usr/bin/env python3
"""HIP 5364 (Eta Ceti): RV fitting and nested sampling.

Example provided by Trifon Trifonov with the aid of CODEX; use with caution. Follows the nu Oph notebook workflow.
Run: python examples/HIP5364_RV.py [--nested] [--threads 1]
"""
from example_tools import (DATA, EXAMPLES, rv, np, dill, plt, uniform_bound,
                           check_priors, fit_rvs, plot_rvs, save_session,
                           finish, command_line)

TARGET = 'HIP5364'


def build_fit():
    fit = rv.signal_fit(name='HIP 5364 / Eta Ceti — two RV planets')
    fit.init_fit = False
    fit.stellar_mass = fit.params.stellar_mass = 1.7  # illustrative fixed mass [Msun]
    for i, filename in enumerate(('hip5364.vels', 'hip5364_crires.vels')):
        data = np.loadtxt(DATA / filename)
        fit.add_dataset(filename, str(DATA / filename), float(np.median(data[:, 1])),
                        0., usejitter=False)
    # K [m/s], P [d], e, omega [deg], M0 [deg], i [deg], Omega [deg].
    # Starting guesses follow the original Eta Ceti example in this repository.
    fit.add_planet(50., 400., .10, 200., 230., 90., 0., usew=True)
    fit.add_planet(50., 770., .10, 170., 170., 90., 0., usew=True)
    set_priors(fit)
    return fit


def set_priors(fit):
    # Fixed tutorial priors; do not tighten these around posterior estimates.
    bounds = {
        0: {'K': (1., 150.), 'P': (300., 500.), 'e': (0., .6),
            'w': (0., 360.), 'M0': (0., 360.)},
        1: {'K': (1., 150.), 'P': (600., 1000.), 'e': (0., .6),
            'w': (0., 360.), 'M0': (0., 360.)},
    }
    for i, parameters in bounds.items():
        for key, limits in parameters.items():
            uniform_bound(fit, key, i, limits)
    for i in range(2):
        uniform_bound(fit, 'rvoff', i, (-2000., 2000.), plural=True)
        uniform_bound(fit, 'jitt', i, (0., 100.), plural=True)


def optimize(fit, dynamical=True):
    # First stabilize the orbital solution at fixed e and zero jitter.
    fit_rvs(fit, starts=10)
    for i in range(2):
        fit.e_use[i] = True
        fit.jitt[i] = 5.  # nonzero start so Simplex can explore the noise scale
        fit.jitt_use[i] = True
    fit_rvs(fit, starts=10)
    keplerian = dill.copy(fit)
    if dynamical:
        fit_rvs(fit, dynamical=True, starts=5)
        print(f'Delta lnL (dynamical − Keplerian) = {fit.loglik-keplerian.loglik:.4f}')
    check_priors(fit)
    return fit, keplerian


def main():
    args = command_line()
    output = (args.output or EXAMPLES / 'results' / TARGET).expanduser().resolve()
    fit, keplerian = optimize(build_fit())
    save_session(keplerian, output, TARGET + '_keplerian')
    plot_rvs(keplerian, output, TARGET + '_keplerian')
    plot_rvs(fit, output, TARGET + '_dynamical')
    finish(fit, output, TARGET, nested=args.nested, threads=args.threads)
    if args.show:
        plt.show()
    else:
        plt.close('all')


if __name__ == '__main__':
    main()
