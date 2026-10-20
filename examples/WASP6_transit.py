#!/usr/bin/env python3
"""WASP-6: fit the bundled detrended transit light curve.

Example provided by Trifon Trifonov with the aid of CODEX; use with caution. Follows the nu Oph notebook workflow.
Run: python examples/WASP6_transit.py [--nested] [--threads 1]
"""
from example_tools import (DATA, EXAMPLES, rv, np, plt, uniform_bound,
                           check_priors, plot_transit, finish, command_line)

TARGET = 'WASP6'


def build_fit():
    fit = rv.signal_fit(name='WASP-6 — transit-only example')
    fit.type_fit = {'RV': False, 'Transit': True, 'TTV': False, 'AST': False}
    fit.rtg = [False, False, True, False]  # RV, RV GP, transit, transit GP
    fit.init_fit = False
    fit.stellar_mass = fit.params.stellar_mass = .9  # illustrative; fixed
    fit.add_transit_dataset('WASP-6 detrended', str(DATA / 'WASP_6_lc_clean.tran'))
    # Only columns 1–3 are used for time/flux/error; the extra columns in the
    # clean file record its detrending provenance. No airmass term is fitted.
    fit.tra_data_sets[0][3] = np.zeros_like(fit.tra_data_sets[0][0])
    fit.tra_data_sets_init[0][3] = np.zeros_like(fit.tra_data_sets[0][0])
    fit.add_planet(K=0., P=3.361, e=0., w=90., M0=0., i=88.5,
                   useK=False, useP=True, usee=False, usew=False,
                   useM0=False, usei=True)
    # The detrended file spans 27 days and multiple transits. Fit the period,
    # conjunction time, a/Rstar, radius ratio and inclination together.
    fit.t0[0] = 2458354.03
    fit.epoch = fit.t0[0]
    fit.pl_a[0] = 11.
    fit.pl_rad[0] = .145
    fit.t0_use[0] = fit.pl_a_use[0] = fit.pl_rad_use[0] = True
    fit.tra_off[0] = 0.
    fit.tra_jitt[0] = .0005
    fit.tra_jitt_use[0] = False
    fit.ld_m[0] = 'quadratic'
    fit.ld_u[0] = fit.ld_u_quad[0] = [.12, .35]
    fit.ld_u_quad_use[0] = [False, False]
    fit.tra_lintr_use[0] = fit.tra_quadtr_use[0] = False
    fit.tra_model_fact = 10  # integer oversampling factor for NumPy linspace
    fit.SciPy_min_use_1 = 'Nelder-Mead'
    fit.SciPy_min_N_use_1 = 2
    fit.SciPy_min_N_use_2 = 0
    fit.Simplex_opt.update(disp=False, maxiter=10000, maxfev=20000, xatol=1e-7, fatol=1e-6)
    fit.Simplex_custom_steps = True
    set_priors(fit)
    return fit


def set_priors(fit):
    for key, limits in {'P': (3.35, 3.37), 't0': (2458353.85, 2458354.2), 'pl_a': (5., 20.),
                        'pl_rad': (.05, .25), 'i': (80., 90.)}.items():
        uniform_bound(fit, key, 0, limits)
    uniform_bound(fit, 'tra_off', 0, (-.02, .02), plural=True)
    uniform_bound(fit, 'tra_jitt', 0, (0., .01), plural=True)


def optimize(fit):
    rv.run_SciPyOp(fit)
    if fit.bound_error or not np.isfinite(fit.loglik):
        raise RuntimeError(fit.bound_error_msg or 'Nonfinite transit likelihood')
    fit.tra_jitt_use[0] = True
    rv.run_SciPyOp(fit)
    if fit.bound_error or not np.isfinite(fit.loglik):
        raise RuntimeError(fit.bound_error_msg or 'Nonfinite transit likelihood')
    print(f'Transit lnL = {fit.loglik:.4f}')
    for key in ('P', 't0', 'pl_a', 'pl_rad', 'i'):
        print(key, getattr(fit, key)[0])
    print('Optimizer diagnostics:', fit.scipy_fit_diagnostics)
    check_priors(fit)
    return fit


def main():
    args = command_line()
    output = (args.output or EXAMPLES / 'results' / TARGET).expanduser().resolve()
    fit = optimize(build_fit())
    plot_transit(fit, output, TARGET)
    finish(fit, output, TARGET, nested=args.nested, threads=args.threads)
    if args.show:
        plt.show()
    else:
        plt.close('all')


if __name__ == '__main__':
    main()
