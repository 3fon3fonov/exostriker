"""Local Gaussian uncertainties for a scalar negative-log-probability objective.

This is independent of the optimizer. It never substitutes a pseudo-inverse for
an unidentified covariance or presents a boundary-truncated Gaussian as valid.
"""
import sys
import numpy as np


def numerical_uncertainties(objective, parameters, bounds, scales=None, *, converged=True, progress=False):
    """Return checked covariance and symmetric 1-sigma errors, or a diagnostic.

    Central differences at two step sizes must agree in the Hessian's whitened
    coordinates. Derivatives use scaled coordinates, then covariance is returned
    in the original parameter units. Up to four progressively smaller stencils
    are tried; the usual successful case costs 4*n*n + 1 objective evaluations.
    Set progress=True to show evaluation count, rate, elapsed time and ETA.
    """
    x = np.asarray(parameters, dtype=float).copy()
    n = len(x)
    result = dict(status='unavailable', message='', nfev=0, covariance=None,
                  correlation=None, errors=np.full(n, np.nan), hessian=None)

    def unavailable(message, indices=None):
        result['message'] = message
        if indices is not None:
            result['parameter_indices'] = np.flatnonzero(indices).tolist()
        return result

    if not n:
        result.update(status='skipped', message='No free parameters.')
        return result
    if not converged:
        return unavailable('The final optimizer stage did not converge; refine the fit first.')
    b = np.asarray(bounds, dtype=float)
    scale = np.maximum(np.abs(x), 1.) if scales is None else np.asarray(scales, dtype=float).copy()
    if (x.ndim != 1 or b.shape != (n, 2) or scale.shape != x.shape
            or not np.isfinite(x).all() or not np.isfinite(scale).all()
            or np.any(scale <= 0) or np.isnan(b).any() or np.any(b[:, 0] >= b[:, 1])):
        return unavailable('Invalid parameter values, scales or bounds.')
    width = b[:, 1] - b[:, 0]
    scale = np.minimum(scale, width)
    distance = np.minimum(x - b[:, 0], b[:, 1] - x)
    if np.any(distance <= 1e-7 * scale):
        return unavailable('A free parameter is at or too close to a bound.', distance <= 1e-7 * scale)
    bar = None
    if progress:
        from tqdm.std import tqdm
        bar = tqdm(total=4*n*n + 1, desc='SciPy Hessian (pass 1/4)', unit='eval',
                   file=sys.stdout, mininterval=.5, miniters=1, ncols=110, ascii=True)
    steps = np.minimum(np.finfo(float).eps**.25 * scale, .2 * distance)
    minimum = 16 * np.abs(np.spacing(np.maximum(np.abs(x), 1.)))

    def evaluate(offset):
        point = x + offset
        result['nfev'] += 1
        try:
            value = float(objective(point.copy()))
        finally:
            if bar is not None:
                bar.update(1)
        if not np.isfinite(value):
            raise ValueError('The objective is non-finite near the best fit.')
        return value

    def derivatives(h):
        f0 = base
        hz = h / scale
        hessian, gradient = np.zeros((n, n)), np.zeros(n)
        for i in range(n):
            ei = np.zeros(n)
            ei[i] = h[i]
            plus, minus = evaluate(ei) - f0, evaluate(-ei) - f0
            hessian[i, i] = (plus + minus) / hz[i]**2
            gradient[i] = (plus - minus) / (2 * hz[i])
            for j in range(i):
                ej = np.zeros(n)
                ej[j] = h[j]
                mixed = ((evaluate(ei + ej) - f0) - (evaluate(ei - ej) - f0)
                         - (evaluate(-ei + ej) - f0) + (evaluate(-ei - ej) - f0))
                hessian[i, j] = hessian[j, i] = mixed / (4 * hz[i] * hz[j])
        return hessian, gradient

    reason = 'Numerical curvature could not be resolved.'
    try:
        base = evaluate(np.zeros(n))
        for attempt in range(4):
            if bar is not None and attempt:
                bar.total = 1 + (attempt + 1) * 4*n*n
                bar.set_description(f'SciPy Hessian (pass {attempt + 1}/4; retry)')
            if np.any(steps / 2 <= minimum):
                return unavailable('Derivative steps are too small to resolve safely.')
            coarse, _ = derivatives(steps)
            hessian, gradient = derivatives(steps / 2)
            if not np.isfinite(hessian).all() or np.any(np.diag(hessian) <= 0):
                reason = 'The curvature is not positive definite; a parameter may be unconstrained.'
                steps /= 4
                continue
            d = np.sqrt(np.diag(hessian))
            normalized = hessian / np.outer(d, d)
            eigenvalues = np.linalg.eigvalsh(normalized)
            if eigenvalues[0] <= 1e-10 or eigenvalues[-1] / eigenvalues[0] > 1e10:
                reason = 'The curvature is singular or poorly conditioned; parameters are degenerate.'
                steps /= 4
                continue
            chol = np.linalg.cholesky(normalized)
            difference = (hessian - coarse) / np.outer(d, d)
            whitened = np.linalg.solve(chol, difference)
            whitened = np.linalg.solve(chol, whitened.T).T
            if np.linalg.norm(whitened, ord=2) > .05:
                reason = 'Curvature changes with derivative step size; uncertainties are unreliable.'
                steps /= 4
                continue
            covariance_z = np.linalg.solve(normalized, np.eye(n)) / np.outer(d, d)
            if float(gradient @ covariance_z @ gradient) > .01:
                return unavailable('The solution is not sufficiently stationary; refine the fit first.')
            covariance = covariance_z * np.outer(scale, scale)
            errors = np.sqrt(np.diag(covariance))
            if not np.isfinite(covariance).all() or not np.isfinite(errors).all():
                return unavailable('The covariance contains non-finite values.')
            if np.any(errors >= distance):
                return unavailable('A 1-sigma interval reaches a parameter bound; use profiles or sampling.', errors >= distance)
            result.update(status='ok', message='Approximate symmetric 1-sigma errors from local curvature.',
                          covariance=covariance, correlation=covariance / np.outer(errors, errors),
                          errors=errors, hessian=hessian / np.outer(scale, scale))
            return result
    except Exception as exc:
        return unavailable(f'Uncertainty calculation failed: {exc}')
    finally:
        if bar is not None:
            bar.close()
    return unavailable(reason)
