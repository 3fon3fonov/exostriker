"""GUI-free calculations executed against a private session in a child process."""


def optimize(fit, *, ff=20, m_ln=True, doGP=False, gp_kernel_id=-1,
             minimize_fortran=True, discovery=False, reset_errors=True, **options):
    import dill
    import RV_mod as rv
    old_errors = dill.copy(fit.param_errors.GP_params_errors)
    options.update(doGP=doGP, kernel_id=gp_kernel_id, minimize_fortran=minimize_fortran)
    if discovery:
        fit.fitting(minimize_loglik=False, amoeba_starts=ff, **options)
    elif minimize_fortran and not doGP and ff > 0 and not fit.init_fit:
        fit.fitting_sequence(**options)
    else:
        fit.fitting(minimize_loglik=True if ff == 0 or not minimize_fortran or doGP else m_ln,
                    amoeba_starts=0 if not minimize_fortran else ff, **options)
    if getattr(fit, 'bound_error', False):
        return fit
    if fit.doGP:
        rv.get_RV_gps_model(fit, get_lnl=fit.get_GP_lnl)
    for index in range(9):
        if fit.type_fit.get('RV', False) and fit.use_planet[index]:
            try:
                rv.phase_RV_planet_signal(fit, index)
            except Exception:
                rv.phase_RV_planet_signal_old(fit, index)
    if not reset_errors and fit.init_fit:
        fit.param_errors.GP_params_errors = old_errors
    return fit


def calculate(fit, operation, options):
    import RV_mod as rv
    print(f'Running {operation} calculation in an isolated process…', flush=True)
    if operation == 'optimize':
        return optimize(fit, **options)
    if operation == 'scipy':
        fit.init_fit = options['ff'] == 0
        rv.run_SciPyOp(fit)
        return fit
    if operation == 'nested':
        return rv.run_nestsamp(fit, **options)
    if operation == 'mcmc':
        return rv.run_mcmc(fit, **options)
    if operation == 'stability':
        return rv.run_stability(fit, **options)
    if operation == 'stability_rebound':
        return rv.run_stability_rebound(fit, **options)
    if operation == 'stability_arb':
        return rv.run_stability_arb(fit, **options)
    if operation == 'corner':
        rv.cornerplot(fit, **options)
        return fit
    if operation == 'mlp':
        import mlp
        if len(fit.fit_results.rv_model.jd) <= 5:
            return fit
        residual = options.pop('residual')
        gp = options.pop('gp')
        if gp and not fit.doGP:
            raise ValueError('No GP model available for residual MLP')
        model = fit.fit_results.rv_model
        values = model.o_c if residual or gp else model.rvs
        if gp:
            values = values - fit.gp_model_data[0]
        datasets = [(model.jd[fit.fit_results.idset == i],
                     values[fit.fit_results.idset == i],
                     model.rv_err[fit.fit_results.idset == i]) for i in range(fit.ndset)]
        fit.mlp = mlp.Gls(datasets, fast=True, verbose=False, nojit=True, norm='dlnL', **options)
        return fit
    if operation == 'tls':
        import numpy as np
        from transitleastsquares_ES import transitleastsquares_ES
        residual = options.pop('residual')
        slots = [i for i in range(20) if len(fit.tra_data_sets[i])]
        values = np.concatenate([fit.transit_results[1][4 if residual else 1][i]
                                 for i in slots])
        if residual:
            values = values + 1.0
        times = np.concatenate([fit.transit_results[1][0][i] for i in slots])
        fit.tls = transitleastsquares_ES(times, values).power(**options)
        return fit
    raise ValueError(f'Unknown calculation: {operation}')
