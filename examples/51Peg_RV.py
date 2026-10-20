#!/usr/bin/env python3
"""51 Peg: a circular one-planet RV model with two instruments.

Example provided by Trifon Trifonov with the aid of CODEX; use with caution. Follows the nu Oph notebook workflow.
Run: python examples/51Peg_RV.py [--nested] [--threads 1]
"""
from example_tools import (DATA, EXAMPLES, rv, np, plt, uniform_bound,
                           check_priors, fit_rvs, plot_rvs, finish, command_line)

TARGET = '51Peg'


def build_fit():
    fit = rv.signal_fit(name='51 Peg — circular RV orbit')
    fit.init_fit = False
    fit.stellar_mass = fit.params.stellar_mass = 1.0  # illustrative fixed mass [Msun]
    for filename in ('51Peg_ELODIE.dat', '51Peg_HARPS.dat'):
        data = np.loadtxt(DATA / filename)
        # The two files have very different velocity zero points, in m/s.
        fit.add_dataset(filename, str(DATA / filename), float(np.median(data[:, 1])),
                        3., usejitter=False)
    fit.add_planet(55., 4.230785, 0., 0., 180., 90., 0., usee=False, usew=False)
    set_priors(fit)
    return fit


def set_priors(fit):
    # Circular orbit: e and omega stay fixed, avoiding an undefined omega at e=0.
    for key, limits in {'K': (1., 100.), 'P': (4.22, 4.24), 'M0': (0., 360.)}.items():
        uniform_bound(fit, key, 0, limits)
    for i, limits in enumerate(((-34000., -32500.), (-500., 500.))):
        uniform_bound(fit, 'rvoff', i, limits, plural=True)
        uniform_bound(fit, 'jitt', i, (0., 50.), plural=True)


def optimize(fit):
    # Refine the known short period with LM before likelihood Simplex.
    # A 5% period perturbation spans many aliases over this long baseline.
    fit.fitting(minimize_fortran=True, minimize_loglik=False,
                outputfiles=[1, 1, 1], amoeba_starts=5, fortran_kill=30)
    for i in range(2):
        fit.jitt_use[i] = True
    fit_rvs(fit, starts=10)
    check_priors(fit)
    return fit


def main():
    args = command_line()
    output = (args.output or EXAMPLES / 'results' / TARGET).expanduser().resolve()
    fit = optimize(build_fit())
    plot_rvs(fit, output, TARGET)
    finish(fit, output, TARGET, nested=args.nested, threads=args.threads)
    if args.show:
        plt.show()
    else:
        plt.close('all')


if __name__ == '__main__':
    main()
