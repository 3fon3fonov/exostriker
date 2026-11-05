"""Nesvorny & Vokrouhlicky (2016), ApJ 823, 72, equations 11--24.

First-order, planar, low-eccentricity normal form. G=1, star mass=1 units;
Psi, psi and delta are dimensionless and independent of the chosen units.
Coefficients f1,f2 are evaluated at nominal resonance, as in Table 1 and
published Figures 1--2. This is an approximate model, not a stability test.
"""
from functools import lru_cache
import numpy as np

NV2016_REFERENCE = 'https://doi.org/10.3847/0004-637X/823/2/72'
DELTA_CRITICAL = (27/32)**(1/3)


@lru_cache(maxsize=32)
def resonant_coefficients(k):
    """Equation 11 including the 2:1 indirect term -1/sqrt(alpha)."""
    from scipy.integrate import quad
    if isinstance(k, bool) or int(k) != k or k < 2:
        raise ValueError('Require first-order k:(k-1), integer k >= 2')
    alpha = ((k-1)/k)**(2/3)
    def laplace(m):
        def denominator(x):
            return 1-2*alpha*np.cos(x)+alpha*alpha
        b = 2/np.pi*quad(lambda x: np.cos(m*x)/np.sqrt(denominator(x)),
                         0, np.pi, epsabs=1.e-12, epsrel=1.e-12)[0]
        db = 2/np.pi*quad(lambda x: np.cos(m*x)*(np.cos(x)-alpha)/denominator(x)**1.5,
                          0, np.pi, epsabs=1.e-12, epsrel=1.e-12)[0]
        return b, db
    b, db = laplace(k)
    c, dc = laplace(k-1)
    return -k*b-alpha*db/2, (k-.5)*c+alpha*dc/2-(alpha**-.5 if k == 2 else 0.)


def normal_form(masses, axes, eccentricities, mean_longitudes, varpis, k=2):
    """Map heliocentric canonical planar elements to the NV2016 variables.

    masses: planet/star ratios (inner, outer). Other arrays: (time, 2),
    axes in any common length unit; angles in radians. Reference axes and
    constants are frozen at the first sample. delta_series measures drift
    of the approximate invariants; delta is its initial value.
    """
    mass = np.asarray(masses, dtype=float)
    a, e, longitude, varpi = [np.atleast_2d(np.asarray(x, dtype=float))
                              for x in (axes, eccentricities, mean_longitudes, varpis)]
    if mass.shape != (2,) or not len(a) or a.shape[1:] != (2,) or any(x.shape != a.shape for x in (e, longitude, varpi)):
        raise ValueError('Require two planet masses and matching (time, 2) element arrays')
    if (not np.all(np.isfinite(mass)) or np.any(mass <= 0)
            or not all(np.all(np.isfinite(x)) for x in (a, e, longitude, varpi))
            or np.any(a <= 0) or np.any(e < 0) or np.any(e >= 1)
            or np.any(a[:, 0] >= a[:, 1])):
        raise ValueError('Require finite, bound, radially ordered orbits and positive masses')
    f = np.array(resonant_coefficients(k))
    total = 1+mass
    reduced = mass/total
    action = reduced*np.sqrt(total*a)
    reference = action[0]
    n = np.sqrt(total/a[0]**3)
    # Stable evaluation of 1-sqrt(1-e^2), including circular orbits.
    gamma = action*e**2/(1+np.sqrt(1-e**2))
    A, B = f/np.sqrt(reference)
    norm = np.hypot(A, B)
    sigma = k*longitude[:, 1, None]-(k-1)*longitude[:, 0, None]-varpi
    z = np.sqrt(2*gamma)*np.exp(1j*sigma)
    resonant = (A*z[:, 0]+B*z[:, 1])/norm
    free = (B*z[:, 0]-A*z[:, 1])/norm
    phi1, phi2 = abs(resonant)**2/2, abs(free)**2/2
    nu = 1.5*((k-1)**2*n[0]/reference[0]+k*k*n[1]/reference[1])
    C = np.prod(mass)*norm/a[0, 1]
    scale = (nu/C)**(2/3)
    # Eq. 15 after substituting K1=Lambda1+(k-1)*sum(Gamma),
    # K2=Lambda2-k*sum(Gamma). Avoid subtracting large Kepler terms.
    detuning = k*n[1]-(k-1)*n[0]
    changes = action-reference
    ns_over_2nu = detuning/(2*nu)+gamma.sum(axis=1)-3/(2*nu)*(
        k*n[1]*changes[:, 1]/reference[1]-(k-1)*n[0]*changes[:, 0]/reference[0])
    delta_series = scale*(ns_over_2nu-phi2)
    Psi = scale*phi1
    # Compute Cartesian combinations without dividing by a vanishing amplitude.
    magnitude = abs(resonant)
    unit = np.divide(resonant, magnitude, out=np.zeros_like(resonant), where=magnitude > 0)
    psi = np.where(magnitude > 0, np.angle(resonant), np.nan)
    x, y = Psi*unit.real, Psi*unit.imag
    delta = float(delta_series[0])
    return dict(Psi=Psi, psi=psi, delta=delta, delta_series=delta_series,
                Psi_cos_psi=x, Psi_sin_psi=y, Phi2_scaled=scale*phi2,
                hamiltonian=-(Psi-delta)**2-np.sqrt(2*Psi)*unit.real,
                f1=float(f[0]), f2=float(f[1]), alpha_res=((k-1)/k)**(2/3))


def _jacobi_cartesian(history, central_mass):
    """Vectorized elliptic Jacobi elements -> position/velocity (G=1)."""
    a, e = history[:, 2], history[:, 3]
    inc, node, omega, mean = np.radians(history[:, 4:8]).T
    mean = (mean+np.pi) % (2*np.pi)-np.pi
    # Bracketed Newton iterations converge also for large elliptic eccentricity.
    E = np.where(e < .8, mean, np.where(mean >= 0, np.pi, -np.pi))
    lower, upper = np.full_like(E, -np.pi), np.full_like(E, np.pi)
    for _ in range(64):
        residual = E-e*np.sin(E)-mean
        lower = np.where(residual < 0, E, lower)
        upper = np.where(residual > 0, E, upper)
        candidate = E-residual/(1-e*np.cos(E))
        next_E = np.where((candidate < lower) | (candidate > upper), (lower+upper)/2, candidate)
        if np.max(abs(next_E-E)) < 1.e-13:
            E = next_E
            break
        E = next_E
    cosO, sinO, cosw, sinw, cosi, sini = np.cos(node), np.sin(node), np.cos(omega), np.sin(omega), np.cos(inc), np.sin(inc)
    p = np.column_stack((cosO*cosw-sinO*sinw*cosi, sinO*cosw+cosO*sinw*cosi, sinw*sini))
    q = np.column_stack((-cosO*sinw-sinO*cosw*cosi, -sinO*sinw+cosO*cosw*cosi, cosw*sini))
    root = np.sqrt(1-e*e)
    pos = a[:, None]*((np.cos(E)-e)[:, None]*p+(root*np.sin(E))[:, None]*q)
    factor = np.sqrt(central_mass/a)/(1-e*np.cos(E))
    vel = factor[:, None]*(-np.sin(E)[:, None]*p+(root*np.cos(E))[:, None]*q)
    return pos, vel


def poincare_elements(histories, masses, stellar_mass, pair):
    """SWIFT Jacobi histories -> heliocentric canonical elements for a pair.

    Convert all planets to astrocentric positions and barycentric momenta,
    then use v_canonical = p/reduced_mass (NV2016 eqs. 7--12). Orient
    longitudes in the pair's initial angular-momentum plane, making results
    independent of a common sky-plane rotation. Finite mutual inclinations
    are only a planar approximation; report them explicitly.
    """
    masses = np.asarray(masses, dtype=float)/stellar_mass
    positions, velocities = [], []
    center_r = np.zeros((len(histories[0]), 3))
    center_v = center_r.copy()
    enclosed = 1.
    for history, mass in zip(histories, masses):
        r, v = _jacobi_cartesian(np.asarray(history), enclosed+mass)
        r, v = r+center_r, v+center_v
        positions.append(r)
        velocities.append(v)
        center_r += mass*r/(enclosed+mass)-mass*center_r/(enclosed+mass)
        center_v += mass*v/(enclosed+mass)-mass*center_v/(enclosed+mass)
        enclosed += mass
    bary_v = sum(m*v for m, v in zip(masses, velocities))/(1+masses.sum())
    h = [np.cross(positions[j], velocities[j]-bary_v) for j in pair]
    normal = sum(masses[j]*moment[0] for j, moment in zip(pair, h))
    if not np.isfinite(normal).all() or np.linalg.norm(normal) == 0:
        raise ValueError('Cannot define the pair angular-momentum plane')
    normal /= np.linalg.norm(normal)
    direction = np.eye(3)[np.argmin(abs(normal))]
    xaxis = direction-normal*np.dot(direction, normal)
    xaxis /= np.linalg.norm(xaxis)
    yaxis = np.cross(normal, xaxis)
    arrays = [[] for _ in range(4)]
    for j in pair:
        r, v = positions[j], (1+masses[j])*(velocities[j]-bary_v)
        gm = 1+masses[j]
        radius = np.linalg.norm(r, axis=1)
        axis = 1/(2/radius-np.sum(v*v, axis=1)/gm)
        if not np.all(np.isfinite(axis)) or np.any(axis <= 0):
            raise ValueError('Unbound heliocentric canonical orbit; NV2016 model unavailable')
        angular = np.cross(r, v)
        eccentricity_vector = np.cross(v, angular)/gm-r/radius[:, None]
        eccentricity = np.linalg.norm(eccentricity_vector, axis=1)
        # Standard Omega+omega measured in the fixed pair reference plane.
        # Use the nonsingular longitude directly when the orbit lies in it.
        hunit = angular/np.linalg.norm(angular, axis=1)[:, None]
        nodes = np.cross(normal, angular)
        has_node = np.linalg.norm(nodes, axis=1) > 1.e-10*np.linalg.norm(angular, axis=1)
        Omega = np.arctan2(nodes@yaxis, nodes@xaxis)
        omega = np.arctan2(np.sum(np.cross(nodes, eccentricity_vector)*hunit, axis=1),
                           np.sum(nodes*eccentricity_vector, axis=1))
        varpi = np.where(has_node, Omega+omega,
                         np.arctan2(eccentricity_vector@yaxis, eccentricity_vector@xaxis))
        # e*sin(E), e*cos(E); avoid division by e at circularity.
        esin = np.sum(r*v, axis=1)/np.sqrt(gm*axis)
        ecos = 1-radius/axis
        mean = np.arctan2(esin, ecos)-esin
        longitude = mean+varpi
        argument_latitude = np.arctan2(np.sum(np.cross(nodes, r)*hunit, axis=1), np.sum(nodes*r, axis=1))
        circular_longitude = np.where(has_node, Omega+argument_latitude, np.arctan2(r@yaxis, r@xaxis))
        longitude = np.where(eccentricity > 1.e-12, longitude, circular_longitude)
        for target, value in zip(arrays, (axis, eccentricity, longitude, varpi)):
            target.append(value)
    mutual = np.degrees(np.arccos(np.clip(np.sum(h[0]*h[1], axis=1)/
                         (np.linalg.norm(h[0], axis=1)*np.linalg.norm(h[1], axis=1)), -1, 1)))
    return [np.column_stack(x) for x in arrays], mutual


def structure_curves(delta):
    """Equilibria and separatrix intersections in (delta, Psi*cos(psi))."""
    delta = np.asarray(delta, dtype=float)
    curves = {name: np.full(delta.shape, np.nan) for name in
              ('stable1', 'stable2', 'unstable', 'separatrix_lower', 'separatrix_upper')}
    for index, d in np.ndenumerate(delta):
        roots = np.roots([1., 0., -2*d, 1.])  # canonical X=sqrt(2*Psi)*cos(psi)
        roots = sorted(float(r.real) for r in roots if abs(r.imag) < 1.e-7)
        negative = roots[0]
        curves['stable1'][index] = -negative**2/2
        if d > DELTA_CRITICAL and len(roots) == 3:
            small, saddle = roots[1:]
            curves['stable2'][index], curves['unstable'][index] = small**2/2, saddle**2/2
            offset = np.sqrt(max(0., 4*d-2*saddle*saddle))
            left, right = -saddle-offset, -saddle+offset
            curves['separatrix_lower'][index] = -left**2/2
            curves['separatrix_upper'][index] = np.sign(right)*right**2/2
    return curves


def export_resonance(result, directory):
    """Export weighted draw coordinates, unique mapped traces, and structure figure.

    Structure figure uses epoch coordinates (not products of separately averaged
    Psi and psi). Histories use fixed initial scaling and retain delta drift.
    """
    import csv
    from pathlib import Path
    directory = Path(directory)
    fields = ['delta_initial', 'Psi_initial', 'psi_initial_deg',
              'Psi_cos_psi_initial', 'Psi_sin_psi_initial', 'delta_mean',
              'delta_semiamp', 'Psi_mean', 'Psi_semiamp', 'Psi_cos_psi_mean',
              'Psi_sin_psi_mean', 'hamiltonian_semiamp', 'Phi2_scaled_semiamp',
              'mutual_i_max_deg', 'eccentricity_max']
    valid_rows = []
    unique = {}
    with (directory/'resonance.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['sample_index', 'weight', 'status', *fields, 'model_notes'])
        for row, record in enumerate(result['records']):
            resonance = record.get('resonance_nv2016', {})
            summary = resonance.get('summary', {})
            notes = '; '.join(resonance.get('warnings', [])) or resonance.get('error', '')
            writer.writerow([int(result['indices'][row]), result['weights'][row], record['status'],
                             *(summary.get(field, np.nan) for field in fields), notes])
            if summary:
                valid_rows.append(row)
                unique.setdefault(int(result['indices'][row]), resonance)
    names = ('time_years', 'Psi', 'psi', 'Psi_cos_psi', 'Psi_sin_psi', 'delta_series', 'hamiltonian')
    lengths = [len(r['history']['time_years']) for r in unique.values()]
    offsets = np.r_[0, np.cumsum(lengths)].astype(int)
    arrays = {name: np.concatenate([r['history'][name] for r in unique.values()]) if unique else np.array([])
              for name in names}
    rows = np.array([[record.get('resonance_nv2016', {}).get('summary', {}).get(field, np.nan)
                      for field in fields] for record in result['records']], dtype=float)
    np.savez_compressed(directory/'resonance.npz', summary=rows, summary_columns=fields,
        indices=result['indices'], weights=result['weights'], stable=result['stable'],
        trace_sample_indices=np.array(list(unique), dtype=int), trace_offsets=offsets,
        pair=result['options']['pair'], mmr=result['options']['mmr'],
        delta_critical=DELTA_CRITICAL, reference=NV2016_REFERENCE, **arrays)
    if not valid_rows:
        return
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    delta, x = rows[:, 0], rows[:, 3]
    good = np.isfinite(delta) & np.isfinite(x) & (result['weights'] > 0)
    if not good.any():
        return
    dmin, dmax = min(-1., delta[good].min()-.3), max(5., delta[good].max()+.3)
    ylim = max(5., np.max(abs(x[good]))*1.15)
    d = np.linspace(dmin, dmax, 1200)
    curves = structure_curves(d)
    fig, ax = plt.subplots(figsize=(6.2, 6.2))
    ax.fill_between(d, curves['separatrix_lower'], curves['separatrix_upper'], color='.9', label='Resonant island (normal form)')
    for name, style, label in [('stable1', '-', 'Stable equilibrium'), ('stable2', '-.', 'Second stable equilibrium'),
                               ('unstable', '--', 'Unstable equilibrium'), ('separatrix_lower', '-', None),
                               ('separatrix_upper', '-', None)]:
        ax.plot(d, curves[name], style, color='black', lw=1, label=label)
    ax.plot(d[d > DELTA_CRITICAL], -d[d > DELTA_CRITICAL], ':', color='black', lw=.8, label=r'$-\delta$ approximation')
    # Resampling for visualization only; exact weights remain in CSV/NPZ.
    ids = np.flatnonzero(good)
    w = result['weights'][ids]/result['weights'][ids].sum()
    selected = np.random.default_rng(result['seed']).choice(ids, min(3000, max(1000, len(ids))), p=w)
    for stable, color, label in [(False, '#b55c30', 'Other / not integrated'), (True, '#007f82', 'N-body stable draws')]:
        chosen = selected[np.asarray(result['stable'])[selected] == stable]
        if len(chosen):
            ax.scatter(delta[chosen], x[chosen], s=8, alpha=.18, c=color, linewidths=0, label=label, rasterized=True)
    ax.axvline(DELTA_CRITICAL, color='.6', lw=.7, zorder=0)
    ax.set(xlim=(dmin, dmax), ylim=(-ylim, ylim), xlabel=r'$\delta$', ylabel=r'$\Psi\cos\psi$',
           title=f'NV2016 {result["options"]["mmr"][0]}:{result["options"]["mmr"][1]} — epoch posterior coordinates')
    ax.tick_params(direction='in', top=True, right=True)
    ax.legend(fontsize=7, loc='best')
    fig.tight_layout()
    fig.savefig(directory/'resonance_structure.pdf')
    fig.savefig(directory/'resonance_structure.png', dpi=180)
    plt.close(fig)
