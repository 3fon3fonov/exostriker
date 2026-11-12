"""Small plotting/fixture helpers for the GetModel synthetic-model notebook."""
from pathlib import Path
from types import SimpleNamespace
import csv

import numpy as np

EPOCH = 2459000.
SPAN = 3200.
STANDARD = dict(RA=120., Dec=-20., parallax=40., mu_alpha_star=40., mu_delta=-30.)


def write_csv(path, columns):
    """Write named one-dimensional NumPy columns without a dataframe dependency."""
    values = [np.asarray(value) for value in columns.values()]
    if not values or any(value.ndim != 1 for value in values):
        raise ValueError('CSV columns must be nonempty named one-dimensional arrays.')
    if len({len(value) for value in values}) != 1:
        raise ValueError('CSV columns must have matching lengths.')
    path = Path(path)
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(zip(*values))
    return path


def giant_planets():
    """Two nearly edge-on, approximately one-Mjup planets in Jacobi orbits."""
    planets = {}
    interior_mass = 1.
    for name, mass, period, eccentricity, omega, anomaly in (
            ('b', 1., 320., .08, 40., 20.),
            ('c', 1.1, 680., .12, 110., 300.)):
        companion_mass = mass/1047.5654817267318
        inclination = 89.98
        velocity = (2*np.pi/(period*86400.))**(1/3)
        velocity *= 1.32712440018e20**(1/3)*companion_mass
        velocity *= np.sin(np.deg2rad(inclination))/(interior_mass+companion_mass)**(2/3)
        velocity /= np.sqrt(1-eccentricity**2)
        planets[name] = dict(K=velocity, P=period, e=eccentricity, omega=omega,
                             M0=anomaly, inclination=inclination, node=0., radius_ratio=.10045)
        interior_mass += companion_mass
    return planets


def earth_observer(times):
    """Offline Astropy built-in Earth ephemeris, barycentric equatorial AU."""
    from astropy.coordinates import get_body_barycentric, solar_system_ephemeris
    from astropy.time import Time
    from astropy import units as u
    with solar_system_ephemeris.set('builtin'):
        return get_body_barycentric('earth', Time(times, format='jd', scale='tdb')).xyz.to_value(u.au).T


def plot_band(ax, band, *, origin=EPOCH, scale=1., label='Best fit', color='black'):
    """Plot nested grey 3/2/1-sigma pointwise bands around the best-fit curve."""
    x = band.BJD-origin
    for level, shade in ((3, '.90'), (2, '.75'), (1, '.55')):
        lower, upper = band.intervals[level]
        ax.fill_between(x, lower*scale, upper*scale, color=shade, label=f'{level}σ credible band')
    ax.plot(x, band.best*scale, color=color, lw=1.2, label=label)
    return ax


def build_demo_session(directory):
    """Fit synthetic RVs, retain transit geometry, and write a GUI-readable .ses.

    Eccentricities, omega, inclination and node are fixed tutorial assumptions.
    K, P, M0 and the instrumental offset are optimized with native LM.
    """
    import dill
    from exostriker.lib import RV_mod as rv
    from exostriker.lib.model_only import GetModel
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    planets = giant_planets()
    times = EPOCH+np.sort(np.random.default_rng(81).uniform(0., SPAN, 180))
    simulation = GetModel(planets, BJD=times, epoch=EPOCH, seed=42,
                          models={'RV': {'noise': {'sigma': 1.5}}})
    data_path = directory/'GetModel_giant_planets.vels'
    np.savetxt(data_path, np.column_stack((times, simulation.RV.data+12., simulation.RV.uncertainty)),
               fmt='%.12f', header='BJD RV_m_per_s uncertainty_m_per_s; synthetic, seed=42; offset=12 m/s')
    fit = rv.signal_fit(name='GetModel synthetic two-giant-planet tutorial')
    fit.init_fit = False
    fit.epoch = EPOCH
    fit.stellar_mass = fit.params.stellar_mass = 1.
    fit.stellar_radius = 1.
    fit.add_dataset('Synthetic RV', str(data_path), 12., 0., useoffset=True, usejitter=False)
    for index, planet in enumerate(planets.values()):
        fit.add_planet(planet['K']*1.01, planet['P']*1.0001, planet['e'], planet['omega'],
                       planet['M0']+.1, planet['inclination'], planet['node'],
                       useK=True, useP=True, usee=False, usew=False, useM0=True, usei=False, usecap=False)
        fit.pl_rad[index] = planet['radius_ratio']
        fit.K_bound[index] = np.array([1., 80.])
        fit.P_bound[index] = np.array([planet['P']*.8, planet['P']*1.2])
        fit.M0_bound[index] = np.array([0., 360.])
    fit.rvoff_bounds[0] = np.array([-100., 100.])
    fit.fit_native(method='LM', maxiter=300, tolerance=1.e-8, use_priors=False)
    result = fit.fit_results.native_joint_result
    if not result.diagnostics['converged'] or not np.all(np.isfinite(result.covariance)):
        raise RuntimeError('Synthetic fixture fit did not converge with finite covariance.')
    fit.getmodel_demo = dict(synthetic=True, seed=42, standard=STANDARD.copy(),
                             fixed_parameters='e, omega, inclination, node, stellar mass, radii')
    # Keep fixture paths portable; the notebook repairs this relative path after loading.
    for entry in fit.filelist.files:
        entry.path = data_path.name
    session_path = directory/'GetModel_giant_planets.ses'
    with session_path.open('wb') as stream:
        dill.dump(fit, stream)
    return fit, session_path


def load_demo_session(path):
    """Load the included tutorial session and repair only its local data path."""
    import dill
    path = Path(path).expanduser().resolve()
    with path.open('rb') as stream:
        fit = dill.load(stream)
    sampler = getattr(fit, 'ns_sampler', None)
    if isinstance(sampler, dict):
        # A legacy Results.__getattr__ may contain bytecode from another Python
        # version. Keep its stored fields without executing restored methods.
        fit.ns_sampler = SimpleNamespace(**dict(dict.items(sampler)))
    fit.cwd = str(Path(__file__).resolve().parents[1]/'exostriker')
    for entry in fit.filelist.files:
        candidate = path.parent/Path(entry.path).name
        if candidate.exists():
            entry.path = str(candidate)
    return fit
