"""Initialize a local MCMC ensemble inside hard bounds with finite posteriors."""
import numpy as np


def initialize_walkers(center, bounds, nwalkers, scale, log_prob, labels=None,
                       rng=None, max_attempts=10):
    center = np.asarray(center, dtype=float)
    bounds = np.asarray(bounds, dtype=float)
    if center.ndim != 1 or not len(center) or bounds.shape != (len(center), 2):
        raise ValueError('MCMC requires at least one fitted parameter and matching bounds')
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError('MCMC starting spread must be finite and positive')
    if not np.isfinite(center).all() or np.isnan(bounds).any():
        raise ValueError('MCMC starting parameters and bounds must not contain NaN')
    labels = list(labels) if labels is not None else list(range(len(center)))
    bad = (bounds[:, 0] >= bounds[:, 1]) | (center < bounds[:, 0]) | (center > bounds[:, 1])
    if bad.any():
        raise ValueError(f'MCMC starting value/bounds are invalid for {labels[np.flatnonzero(bad)[0]]}')
    # A local interval on both sides; at a boundary only its interior is used.
    # Maintain representable spread for large parameters such as BJD epochs.
    width = np.maximum(scale, 32 * np.abs(np.spacing(center)))
    lower = np.maximum(center-width, np.nextafter(bounds[:, 0], bounds[:, 1]))
    upper = np.minimum(center+width, np.nextafter(bounds[:, 1], bounds[:, 0]))
    if np.any(lower >= upper):
        raise ValueError('MCMC bounds/spread leave no representable interval for walkers')
    rng = np.random if rng is None else rng
    positions = np.empty((nwalkers, len(center)))
    probabilities = np.full(nwalkers, -np.inf)
    for _ in range(max_attempts):
        missing = np.flatnonzero(~np.isfinite(probabilities))
        positions[missing] = rng.uniform(lower, upper, size=(len(missing), len(center)))
        probabilities[missing] = log_prob(positions[missing])
        if np.isfinite(probabilities).all():
            return positions, probabilities
    raise ValueError('Cannot initialize MCMC: some walkers have non-finite log probability. '
                     'Check the starting model, fitted-parameter bounds and priors; '
                     'try a valid fitted solution or a smaller starting spread.')
