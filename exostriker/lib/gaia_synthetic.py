"""Synthetic absolute Gaia along-scan measurements in Exo-Striker's format."""
from dataclasses import dataclass
import inspect
from pathlib import Path

import numpy as np
from astropy.time import Time


HEADER_FIELDS = ('RAdeg', 'DEdeg', 'Plx', 'pm_RA', 'pm_DE',
                 'e_RA', 'e_DE', 'e_Plx', 'e_pmRA', 'e_pmDE')
REFERENCE_EPOCH = 2017.5


def _source_id(value):
    """Validate an integer identifier without rounding large Gaia IDs."""
    if isinstance(value, (bool, float, np.floating)):
        raise ValueError('source_id must be a positive integer or digit string.')
    text = str(value)
    if not text.isascii() or not text.isdigit() or int(text) <= 0:
        raise ValueError('source_id must be a positive integer or digit string.')
    return int(text)


def _query_catalogue(source_id):
    """Retrieve the DR3 five-parameter solution and formal uncertainties."""
    source_id = _source_id(source_id)
    try:
        from astroquery.gaia import Gaia
    except ImportError as error:
        raise ImportError('Install astroquery or supply catalogue for offline generation.') from error

    columns = ('ra', 'dec', 'parallax', 'pmra', 'pmdec', 'ra_error',
               'dec_error', 'parallax_error', 'pmra_error', 'pmdec_error',
               'ref_epoch')
    query = ('SELECT ' + ', '.join(columns)
             + f' FROM gaiadr3.gaia_source WHERE source_id = {source_id}')
    result = Gaia.launch_job(query).get_results()
    if len(result) != 1:
        raise ValueError(f'Gaia DR3 source {source_id} was not found.')
    catalogue = {}
    for field, column in zip((*HEADER_FIELDS, 'ref_epoch'), columns):
        value = result[column][0]
        if np.ma.is_masked(value):
            raise ValueError(f'Gaia DR3 source {source_id} has no {column}.')
        catalogue[field] = float(value)
    return catalogue


@dataclass(frozen=True)
class GaiaDR4Simulation:
    """Eight-column data plus noiseless model, orbit, noise and provenance."""

    data: np.ndarray
    model: np.ndarray
    orbital_al: np.ndarray
    baseline_al: np.ndarray
    noise: np.ndarray
    header: dict
    source_id: object
    catalogue_epoch: float

    def write(self, path, *, overwrite=False):
        """Write absolute centroids and the stellar header used by the importer."""
        path = Path(path)
        with path.open('w' if overwrite else 'x') as stream:
            stream.write('# Synthetic Gaia DR4 nominal scanning-law astrometry\n')
            stream.write(f'# Source {self.source_id}; reference epoch J{REFERENCE_EPOCH} TCB\n')
            stream.write('# Absolute AL centroids include parallax and proper motion\n')
            stream.write('# One synthetic observation per nominal field-of-view transit\n')
            stream.write('# Independent Gaussian noise; no Gaia instrument noise model\n')
            stream.write('# Catalogue errors retained without covariance propagation\n')
            stream.write(f'# Original catalogue epoch J{self.catalogue_epoch}\n#\n#\n')
            fields = (*HEADER_FIELDS, 'RefEpoch')
            stream.write('# ' + ' '.join(fields) + '\n')
            stream.write('# ' + ' '.join(str(self.header[key]) for key in fields) + '\n#\n')
            stream.write('# transit_id ccd t_JD_TCB centroid_AL_mas sigma_AL_mas '
                         'parallax_factor_AL scan_angle_deg flag\n')
            np.savetxt(stream, self.data, fmt=['%d', '%d', '%.12f', '%.15g',
                                             '%.15g', '%.15g', '%.15g', '%d'])
        return path


def GetGaiaDR4(planets, source_id=None, *, epoch=None, stellar_mass=None, sigma_al=0.1,
               seed=None, add_noise=True, catalogue=None, scan_data=None,
               output=None, overwrite=False, inclinations=None, Omega=None,
               omega_convention='stellar'):
    """Create absolute along-scan astrometry from GetModel orbital parameters.

    planets uses GetModel's named dictionaries: K [m/s], P [days], e,
    omega, M0, inclination and node [degrees]. M0 refers to epoch [TCB JD];
    stellar_mass is in solar masses. Uses Keplerian stellar reflex orbits.
    sigma_al is a positive scalar or one mas uncertainty per nominal transit.
    Gaussian noise is optional and reproducible with seed.

    source_id queries Gaia DR3 through astroquery. Alternatively catalogue is
    a mapping with HEADER_FIELDS and ref_epoch [Julian year, normally 2016].
    Missing formal errors in a supplied catalogue default to zero. Position
    is propagated linearly to J2017.5; perspective acceleration and covariance
    propagation are omitted, matching the importer's five-parameter model.

    gaiascanlaw supplies the nominal DR4 schedule and AL parallax factors.
    For offline use scan_data is (TCB JD, scan radians north through east,
    AL parallax factor), with equally sized 1-D arrays. One row represents a
    field-of-view transit, not nine independent CCD measurements. IDs are
    synthetic sequential IDs; CCD and flag are zero.

    Returns GaiaDR4Simulation; output optionally writes an importable .d file.
    Neither queries nor optional imports occur until this function is called.

    planets can also be a fit object. Its active planets, epoch and stellar
    mass are copied through SessionModelAdapter. inclinations and Omega
    override angles in degrees: scalar for all planets, a sequence in active
    planet order, or a mapping keyed by planet name (e.g. planet_1).
    The fit is never modified. Dynamical fits require an explicit Keplerian
    orbital configuration instead of silently discarding their dynamics.
    omega_convention='stellar' matches Exo-Striker's Gaia likelihood, which
    uses the fitted omega directly. 'companion' adds 180 degrees, matching
    generic GetModel astrometry's companion-to-stellar conversion.
    """
    from model_only import _planet_parameters, _validate_rv_orbits

    if hasattr(planets, 'use_planet'):
        from model_posterior import SessionModelAdapter

        adapter = SessionModelAdapter(planets, phase='M0')
        if adapter.mode != 'keplerian':
            raise ValueError('Gaia DR4 synthesis currently requires Keplerian orbits; '
                             'supply an explicit orbital configuration for a dynamical fit.')
        planets, fit_options, trends, correction = adapter.inputs()
        if epoch is not None and float(epoch) != fit_options['epoch']:
            raise ValueError('Use the fit orbital epoch; changing epoch requires rephasing M0.')
        epoch = fit_options['epoch']
        if stellar_mass is None:
            stellar_mass = fit_options['stellar_mass']
    if epoch is None:
        raise ValueError('Supply epoch (TCB JD) for an explicit orbital configuration.')
    if stellar_mass is None:
        stellar_mass = 1.
    names, rows, unused_options = _planet_parameters(planets, {})
    rows = rows.copy()
    for values, column, label in ((inclinations, 5, 'inclinations'), (Omega, 6, 'Omega')):
        if values is None:
            continue
        if isinstance(values, dict):
            unknown = values.keys() - set(names)
            if unknown:
                raise ValueError(f'{label}: unknown planet names: ' + ', '.join(sorted(unknown)))
            for index, name in enumerate(names):
                if name in values:
                    rows[index, column] = float(values[name])
        else:
            try:
                rows[:, column] = np.broadcast_to(np.asarray(values, dtype=float), (len(rows),))
            except ValueError as error:
                raise ValueError(f'{label}: supply a scalar or one value per active planet.') from error
    fields = ('K', 'P', 'e', 'omega', 'M0', 'inclination', 'node')
    planets = {name: dict(zip(fields, row)) for name, row in zip(names, rows)}
    _validate_rv_orbits(rows, stellar_mass)
    if not np.all(np.isfinite(rows)):
        raise ValueError('Orbital parameters must be finite.')
    if np.any((rows[:, 0] > 0.) & (np.abs(np.sin(np.deg2rad(rows[:, 5]))) < 1.e-12)):
        raise ValueError('Nonzero K requires inclination strictly between 0 and 180 degrees.')
    if omega_convention not in ('stellar', 'companion'):
        raise ValueError("omega_convention must be 'stellar' or 'companion'.")

    if source_id is not None:
        source_id = _source_id(source_id)
    if catalogue is None:
        if source_id is None:
            raise ValueError('Supply source_id or catalogue.')
        catalogue = _query_catalogue(source_id)
    else:
        catalogue = dict(catalogue)

    required = (*HEADER_FIELDS[:5], 'ref_epoch')
    missing = set(required) - catalogue.keys()
    if missing:
        raise ValueError('Missing catalogue fields: ' + ', '.join(sorted(missing)))
    header = {key: float(catalogue.get(key, 0.)) for key in HEADER_FIELDS}
    catalogue_epoch = float(catalogue['ref_epoch'])
    if not np.all(np.isfinite([*header.values(), catalogue_epoch, epoch])):
        raise ValueError('Catalogue parameters and epoch must be finite.')
    if not 0 <= header['RAdeg'] < 360 or not -90 < header['DEdeg'] < 90:
        raise ValueError('Require RA in [0,360) and Dec strictly between -90 and 90 degrees.')
    if header['Plx'] <= 0 or any(header[key] < 0 for key in HEADER_FIELDS[5:]):
        raise ValueError('Require positive parallax and nonnegative catalogue errors.')

    if scan_data is None:
        try:
            import gaiascanlaw as gs
        except ImportError as error:
            raise ImportError('Install gaiascanlaw or supply scan_data for offline generation.') from error

        if 'parallax_factor' not in inspect.signature(gs.scanlaw).parameters:
            raise ImportError(
                'This gaiascanlaw build does not supply AL parallax factors. '
                'The PyPI 0.2.0 build differs from the GitHub source. '
                'Install the upstream build with the current interpreter: '
                'python -m pip install --upgrade --force-reinstall --no-deps '
                'git+https://github.com/zpenoyre/gaiascanlaw.git '
                '(then restart the Python/Jupyter kernel), or supply scan_data '
                'containing TCB JD, scan angles in radians and AL parallax factors.')

        years, angles, factors, across = gs.scanlaw(
            header['RAdeg'], header['DEdeg'], tend=gs.tdr4,
            obstype='astrometry', parallax_factor=True)
        times = Time(years, format='decimalyear', scale='tcb').jd
    else:
        times, angles, factors = scan_data
    times = np.asarray(times, dtype=float)
    angles = np.asarray(angles, dtype=float)
    factors = np.asarray(factors, dtype=float)
    if (times.ndim != 1 or len(times) == 0 or angles.shape != times.shape
            or factors.shape != times.shape
            or not np.all(np.isfinite([times, angles, factors]))):
        raise ValueError('Scan data must contain three nonempty finite 1-D arrays of equal length.')
    uncertainty = np.broadcast_to(np.asarray(sigma_al, dtype=float), times.shape).copy()
    if not np.all(np.isfinite(uncertainty)) or np.any(uncertainty <= 0):
        raise ValueError('sigma_al must contain finite positive uncertainties in mas.')

    # Record the catalogue reference epoch explicitly. The importer retains
    # its internal J2017.5 time basis and uses RefEpoch for baseline removal.
    header['RefEpoch'] = REFERENCE_EPOCH
    elapsed = REFERENCE_EPOCH - catalogue_epoch
    header['RAdeg'] += header['pm_RA'] * elapsed / (
        3600000. * np.cos(np.deg2rad(header['DEdeg'])))
    header['RAdeg'] %= 360.
    header['DEdeg'] += header['pm_DE'] * elapsed / 3600000.
    if not -90 < header['DEdeg'] < 90:
        raise ValueError('Propagated catalogue position crosses a celestial pole.')

    # Use the same omega and RV-to-angular-axis formula as ast_loglik_gaia /
    # a_from_K_in_mas, including its AU constant. A 180-degree reversal here
    # would generate a signal opposite to the one fitted by Exo-Striker.
    from model_only import GetModel
    standard = [header['RAdeg'], header['DEdeg'], header['Plx'], 0., 0.]
    amplitudes = np.zeros(len(rows))
    nonzero = rows[:, 0] > 0.
    amplitudes[nonzero] = (
        rows[nonzero, 0] * rows[nonzero, 1] * 86400.
        * np.sqrt(1. - rows[nonzero, 2] ** 2) * header['Plx']
        / (2. * np.pi * np.sin(np.deg2rad(rows[nonzero, 5])) * 1.496e11))
    omega_shift = 180. if omega_convention == 'companion' else 0.
    orbits = np.column_stack((
        rows[:, 1], rows[:, 2], np.deg2rad(rows[:, 3] + omega_shift),
        np.deg2rad(rows[:, 5]), np.deg2rad(rows[:, 6]),
        epoch - rows[:, 4] / 360. * rows[:, 1], amplitudes))
    orbit = GetModel(
        planets, BJD=times, epoch=epoch, stellar_mass=stellar_mass,
        models={'RV': False, 'Astrometry': {
            'enabled': True, 'parameters': {
                'standard': standard, 'observer': np.zeros((len(times), 3)),
                'scan_angle': angles, 'orbits': orbits}}}).Astrometry.along_scan

    relative_years = Time(times, format='jd', scale='tcb').jyear - REFERENCE_EPOCH
    scan_east = np.cos(np.pi / 2. - angles)
    scan_north = np.sin(np.pi / 2. - angles)
    baseline = (factors * header['Plx'] + relative_years * (
        scan_east * header['pm_RA'] + scan_north * header['pm_DE']))
    model = baseline + orbit
    # Follow GetModel's seed validation and leave global RNG state untouched.
    if seed is not None and (isinstance(seed, (bool, np.bool_))
                            or not isinstance(seed, (int, np.integer)) or seed < 0):
        raise ValueError('seed must be a nonnegative integer or None.')
    noise = np.zeros(len(times))
    if add_noise:
        noise = np.random.default_rng(seed).normal(0., uncertainty)
    data = np.column_stack((np.arange(1, len(times) + 1), np.zeros(len(times)),
                            times, model + noise, uncertainty, factors,
                            np.rad2deg(angles), np.zeros(len(times))))
    result = GaiaDR4Simulation(data, model, orbit, baseline, noise, header,
                              source_id, catalogue_epoch)
    if output is not None:
        result.write(output, overwrite=overwrite)
    return result
