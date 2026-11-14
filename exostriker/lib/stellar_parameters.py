"""Independent stellar coordinates and propagation to derived quantities."""
import copy
import numpy as np

SOLAR_MASS = 1.9884098709677423e30
SOLAR_RADIUS = 6.957e8
DENSITY_FACTOR = 3 * SOLAR_MASS / (4 * np.pi * SOLAR_RADIUS ** 3)


def ensure_stellar_settings(session):
    """Migrate sessions while retaining the established density/aR modes."""
    from stellar_density import ensure_density_settings
    ensure_density_settings(session)
    defaults = {
        'stellar_mass_radius': False,
        'planet_radius_stellar_source': 'auto',
        'stellar_radius_parameter_index': None,
        'stellar_mass_use': False,
        'stellar_radius_use': False,
        'stellar_rho_use': False,
        'stellar_mass_err': [0., 0.],
        'stellar_radius_err': [0., 0.],
        'stellar_rho_err': [0., 0.],
    }
    for name, value in defaults.items():
        if not hasattr(session, name):
            setattr(session, name, copy.deepcopy(value))

    # Older sessions stored symmetric stellar errors as a scalar.
    for root in ('stellar_mass', 'stellar_radius', 'stellar_rho'):
        name = root + '_err'
        uncertainty = np.asarray(getattr(session, name), dtype=float)
        if uncertainty.size == 1:
            setattr(session, name, np.repeat(uncertainty.reshape(-1), 2))


def radius_mode(session):
    return bool(getattr(session, 'stellar_mass_radius', False))


def fitted_stellar(session):
    """Whether stellar coordinates require the joint inference objective."""
    return bool(getattr(session, 'stellar_mass_use', False) or
                getattr(session, 'stellar_radius_use', False) or radius_mode(session))


def sync_stellar(session):
    """Derive the dependent stellar coordinate using SI density units."""
    ensure_stellar_settings(session)
    mass = float(session.params.stellar_mass)
    session.stellar_mass = mass
    if radius_mode(session):
        if mass <= 0 or session.stellar_radius <= 0:
            raise ValueError('Stellar mass and radius must be positive.')
        session.stellar_rho = DENSITY_FACTOR * mass / session.stellar_radius ** 3
    elif session.rho:
        slots = np.flatnonzero(session.use_planet)
        if len(slots):
            session.stellar_rho = float(session.pl_a[int(slots[0])])
        if mass <= 0 or session.stellar_rho <= 0:
            raise ValueError('Stellar mass and density must be positive.')
        session.stellar_radius = (DENSITY_FACTOR * mass / session.stellar_rho) ** (1 / 3)


def set_radius_mode(session, enabled):
    """Select physical mass/radius geometry, retaining the individual aR route."""
    ensure_stellar_settings(session)
    enabled = bool(enabled)
    if enabled:
        from stellar_density import set_density_mode
        if session.rho:
            sync_stellar(session)
            set_density_mode(session, False)
        session.stellar_rho_use = False
    session.stellar_mass_radius = enabled
    sync_stellar(session)


def validate_stellar_prior(session, root):
    """Validate a stellar coordinate before packing bounds and priors."""
    bounds = getattr(session, root + '_bounds')
    normal = getattr(session, root + '_norm_pr')
    jeffreys = getattr(session, root + '_jeff_pr')
    if not np.all(np.isfinite(bounds)) or not 0 < bounds[0] < bounds[1]:
        raise ValueError(root + ' bounds must be positive and ordered.')
    if normal[2] and jeffreys[2]:
        raise ValueError('Select only one prior family for ' + root + '.')
    if normal[2] and (not np.all(np.isfinite(normal[:2])) or normal[1] <= 0):
        raise ValueError(root + ' Gaussian prior requires a positive finite sigma.')
    if jeffreys[2] and (not np.all(np.isfinite(jeffreys[:2])) or not 0 < jeffreys[0] < jeffreys[1]):
        raise ValueError(root + ' Jeffreys prior requires positive ordered limits.')
    return [max(bounds[0], jeffreys[0]), min(bounds[1], jeffreys[1])] if jeffreys[2] else bounds


def publish_stellar(session, par, errors, covariance=None, free_indices=()):
    """Publish fitted values and propagate full local covariance to rho or R."""
    ensure_stellar_settings(session)
    session.params.stellar_mass = float(par[-1])
    session.stellar_mass_err = np.asarray(errors[-1], float).copy()
    index = session.stellar_radius_parameter_index
    if radius_mode(session) and index is not None:
        session.stellar_radius = float(par[index])
        if index in free_indices:
            session.stellar_radius_err = np.asarray(errors[index], float).copy()
    sync_stellar(session)
    if not (radius_mode(session) or session.rho):
        return
    mass = session.stellar_mass
    if radius_mode(session):
        dependent = 'stellar_rho_err'
        indices = [len(par) - 1, index]
        gradient = [session.stellar_rho / mass, -3 * session.stellar_rho / session.stellar_radius]
    else:
        slots = np.flatnonzero(session.use_planet)
        if not len(slots):
            return
        density_index = 2 * len(session.mcmc_dataset_slots) + 7 * len(slots) + 3
        dependent = 'stellar_radius_err'
        indices = [len(par) - 1, density_index]
        gradient = [session.stellar_radius / (3 * mass), -session.stellar_radius / (3 * session.stellar_rho)]
    free = list(map(int, free_indices))
    weights = np.zeros(len(free))
    for index, derivative in zip(indices, gradient):
        if index in free:
            weights[free.index(index)] = derivative
    active = weights != 0
    if not np.any(active):
        error = 0.
    elif covariance is None:
        error = np.nan
    else:
        block = np.asarray(covariance)[np.ix_(active, active)]
        variance = weights[active] @ block @ weights[active]
        error = np.sqrt(variance) if np.isfinite(variance) and variance >= 0 else np.nan
    setattr(session, dependent, np.array([error, error]))


def stellar_posterior(session, samples, flags, parameters):
    """Derive M, R and rho for each posterior sample, preserving correlations."""
    samples = np.asarray(samples, float)
    flags = list(map(int, flags))
    def coordinate(index):
        if index in flags:
            return samples[:, flags.index(index)]
        return np.full(len(samples), parameters[index])
    mass = coordinate(len(parameters) - 1)
    if radius_mode(session):
        radius = coordinate(session.stellar_radius_parameter_index)
        density = DENSITY_FACTOR * mass / radius ** 3
    elif session.rho:
        slots = np.flatnonzero(session.use_planet)
        if not len(slots):
            return {}
        index = 2 * len(session.mcmc_dataset_slots) + 7 * len(slots) + 3
        density = coordinate(index)
        radius = (DENSITY_FACTOR * mass / density) ** (1 / 3)
    else:
        radius = np.full(len(samples), session.stellar_radius)
        density = DENSITY_FACTOR * mass / radius ** 3
    return dict(stellar_mass=mass, stellar_radius=radius, stellar_rho=density)


def publish_stellar_posterior(session, samples, level, use_mad=False):
    """Report dependent stellar errors from transformed posterior draws."""
    if radius_mode(session):
        root = 'stellar_rho'
    elif session.rho:
        root = 'stellar_radius'
    else:
        return
    values = samples.get(root)
    if values is None or not len(values):
        return
    if use_mad:
        error = np.median(np.abs(values - np.median(values)))
        errors = [error, error]
    else:
        lower, upper = np.percentile(values, [level, 100. - level])
        center = getattr(session, root)
        errors = [center - lower, upper - center]
    setattr(session, root + '_err', np.asarray(errors, float))
