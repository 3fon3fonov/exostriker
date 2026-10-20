"""Physical masses and Jacobi axes using the native RV backend convention.

K = (2*pi*G/P)**(1/3) * m*sin(i)/(M_inner+m)**(2/3)/sqrt(1-e**2).
M_inner includes the star and planets with shorter periods. Units and mass
conversion match MA_J_cop_fin; semi-major axes follow Kepler's third law.
"""
import numpy as np

GMSUN = 1.32712440018e20
AU = 1.49597892e11
TWOPI = 2.0 * np.pi
THIRD = 1.0 / 3.0


def mass_a_from_jacobi_el(pl_par, npl, m0):
    """
    Computes physical planet masses and Jacobi semimajor axes from RV
    Keplerian parameters interpreted as Jacobi elements.

    Input
    -----
    pl_par : list
        List of planet parameter lists:

            pl_par[i][0] = K       [m/s]
            pl_par[i][1] = P       [days]
            pl_par[i][2] = e
            pl_par[i][3] = omega   [deg]
            pl_par[i][4] = M       [deg]
            pl_par[i][5] = incl    [deg]
            pl_par[i][6] = Omega   [deg]

    npl : int
        Number of planets.

    m0 : float
        Stellar mass in solar masses.

    Output
    ------
    pl_mass : ndarray
        Planet masses in Jupiter masses.

    ap : ndarray
        Jacobi semimajor axes in AU.
    """

    pl_mass = np.zeros(npl)
    ap = np.zeros(npl)

    if not np.isfinite(m0) or m0 <= 0:
        raise ValueError("Stellar mass must be finite and positive.")
    if npl != len(pl_par):
        raise ValueError("Planet count does not match the orbital elements.")
    if npl == 0:
        return pl_mass, ap
    elements = np.asarray(pl_par, dtype=float)
    if elements.shape != (npl, 7) or not np.all(np.isfinite(elements)):
        raise ValueError("Expected finite [K, P, e, omega, M, inclination, node] rows.")
    # Use the same period ordering as the native Jacobi backend, then restore
    # the caller's order (which may contain sparse GUI planet slots).
    order = np.argsort(elements[:, 1], kind='stable')

    # cumulative Jacobi inner mass, in solar masses
    m_inner = float(m0)

    MSUN_TO_MJUP = 1047.5654817267318

    for i in order:

        K = abs(float(pl_par[i][0]))          # m/s
        P_days = float(pl_par[i][1])         # days
        e = float(pl_par[i][2])
        inc = float(pl_par[i][5])            # deg

        P = P_days * 86400.0                 # seconds

        if P <= 0.0:
            raise ValueError("Planet %d: period must be positive." % (i+1))

        if e < 0.0 or e >= 1.0:
            raise ValueError("Planet %d: eccentricity must satisfy 0 <= e < 1." % (i+1))

        sini = np.sin(np.radians(inc))

        if np.abs(sini) <= np.finfo(float).eps:
            raise ValueError("Planet %d: sin(incl) is zero. RV mass is undefined." % (i+1))

        kbar = K * np.sqrt(1.0 - e*e) / np.abs(sini)

        # Jacobi RV mass equation:
        #
        # kbar = (2 pi GMSUN / P)^(1/3)
        #        * mp / (m_inner + mp)^(2/3)
        #
        # Here mp and m_inner are both in solar masses.
        coeff = (TWOPI * GMSUN / P)**THIRD

        def f(mp):
            return coeff * mp / (m_inner + mp)**(2.0/3.0) - kbar

        if K == 0.0:
            ap[i] = (GMSUN * m_inner * (P / TWOPI)**2.0)**THIRD / AU
            continue

        # Initial low-mass estimate
        mp_guess = kbar * m_inner**(2.0/3.0) / coeff

        lo = 0.0
        hi = max(2.0 * mp_guess, 1.0e-12)

        while f(hi) < 0.0:
            hi *= 2.0
            if not np.isfinite(hi):
                raise ValueError("RV mass cannot be bracketed for planet %d." % (i+1))

        # Robust bisection
        for j in range(100):
            mid = 0.5 * (lo + hi)
            if f(mid) > 0.0:
                hi = mid
            else:
                lo = mid

        mp = 0.5 * (lo + hi)

        # Jacobi two-body total mass for planet i
        mtot_jacobi = m_inner + mp

        # Jacobi semimajor axis
        ap[i] = (GMSUN * mtot_jacobi * (P / TWOPI)**2.0)**THIRD / AU

        # Output planet mass in Jupiter masses
        pl_mass[i] = mp * MSUN_TO_MJUP

        # Add this planet to the cumulative inner Jacobi mass
        m_inner += mp

    return pl_mass, ap

