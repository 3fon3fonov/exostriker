"""
Provides vector astrometry functions.
"""
#Thank you Anthony Brown from Leiden University for letting me copy this code to safe some work



import numpy as np
from .useful import *
__all__ = [
    "spherical_to_cartesian",
    "cartesian_to_spherical",
    "normal_triad",
    "elementary_rotation_matrix",
    "phase_space_to_astrometry",
    "astrometry_to_phase_space",
    "EpochPropagation"
]


def spherical_to_cartesian(r, phi, theta):
    """
    Convert spherical to Cartesian coordinates. The input can be scalars or
    1-dimensional numpy arrays. Note that the angle coordinates follow the astronomical
    convention of using elevation (declination, latitude) rather than its complement
    (pi/2-elevation), where the latter is commonly used in the mathematical treatment of
    spherical coordinates.

    Parameters
    ----------
    r : float or float array
        Length of input Cartesian vector.
    phi : float or float array
        Longitude-like angle (e.g., right ascension, ecliptic longitude) in radians
    theta : float or float array
        Latitide-like angle (e.g., declination, ecliptic latitude) in radians

    Returns
    -------
    x, y, z : float or float array
        The Cartesian vector components x, y, z
    """
    ctheta = np.cos(theta)
    x = r * np.cos(phi) * ctheta
    y = r * np.sin(phi) * ctheta
    z = r * np.sin(theta)
    return x, y, z


def cartesian_to_spherical(x, y, z):
    r"""
    Convert Cartesian to spherical coordinates. The input can be scalars or
    1-dimensional numpy arrays. Note that the angle coordinates follow the astronomical
    convention of using elevation (declination, latitude) rather than its complement
    (pi/2-elevation), which is commonly used in the mathematical treatment of spherical
    coordinates.

    Parameters
    ----------
    x : float or float array
        Cartesian vector component along the X-axis
    y : float or float array
        Cartesian vector component along the Y-axis
    z : float or float array
        Cartesian vector component along the Z-axis

    Returns
    -------
    r, phi, theta : float or float array
        The spherical coordinates r=np.sqrt(x*x+y*y+z*z), longitude phi, latitude theta.

    Raises
    ------
    ValueError
        When one of the input (x,y,z) coordinates correspond to r=0.

    Notes
    -----
        The longitude angle is between 0 and :math:`+2\pi`.
    """
    rCylSq = x * x + y * y
    r = np.sqrt(rCylSq + z * z)
    if np.any(r == 0.0):
        raise ValueError("Error: one or more of the input points is at distance zero.")
    phi = np.arctan2(y, x)
    phi = np.where(phi < 0.0, phi + 2 * np.pi, phi)
    phi=np.float64(phi)
    return r, phi, np.arctan2(z, np.sqrt(rCylSq))


def normal_triad(phi, theta):
    """
    Calculate the so-called normal triad [p, q, r] associated with a spherical
    coordinate system. The three vectors are:

    p - The unit tangent vector in the direction of increasing longitudinal angle phi.
    q - The unit tangent vector in the direction of increasing latitudinal angle theta.
    r - The unit vector toward the point (phi, theta).

    Parameters
    ----------
    phi : float or array
        Longitude-like angle (e.g., right ascension, ecliptic longitude) in radians
    theta : float or array
        Latitide-like angle (e.g., declination, ecliptic latitude) in radians

    Returns
    -------
    p, q, r : array
        The normal triad as the vectors p, q, r as (N,3) arrays
    """
    sphi = np.sin(phi)
    stheta = np.sin(theta)
    cphi = np.cos(phi)
    ctheta = np.cos(theta)
    p = np.array([-sphi, cphi, np.zeros_like(phi)])
    q = np.array([-stheta * cphi, -stheta * sphi, ctheta])
    r = np.array([ctheta * cphi, ctheta * sphi, stheta])
    return p, q, r


def elementary_rotation_matrix(axis, rotationAngle):
    """
    Construct an elementary rotation matrix describing a rotation around the x, y, or
    z-axis.

    Parameters
    ----------
    axis : str
        Axis around which to rotate ("x", "X", "y", "Y", "z", or "Z")
    rotationAngle : float
        The rotation angle in radians

    Returns
    -------
    rmat : array
        The rotation matrix

    Raises
    ------
    ValueError
        If an unsupported rotation axis string is supplied.

    Examples
    --------
    >>> rotmat = elementaryRotationMatrix("y", np.pi/6.0)
    """
    if axis.upper() == "X":
        return np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, np.cos(rotationAngle), np.sin(rotationAngle)],
                [0.0, -np.sin(rotationAngle), np.cos(rotationAngle)],
            ]
        )
    elif axis.upper() == "Y":
        return np.array(
            [
                [np.cos(rotationAngle), 0.0, -np.sin(rotationAngle)],
                [0.0, 1.0, 0.0],
                [np.sin(rotationAngle), 0.0, np.cos(rotationAngle)],
            ]
        )
    elif axis.upper() == "Z":
        return np.array(
            [
                [np.cos(rotationAngle), np.sin(rotationAngle), 0.0],
                [-np.sin(rotationAngle), np.cos(rotationAngle), 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
    else:
        raise ValueError("Unknown rotation axis " + axis + "!")


def phase_space_to_astrometry(x, y, z, vx, vy, vz):
    r"""
    From the given phase space coordinates calculate the astrometric observables,
    including the radial velocity, which here is seen as the sixth astrometric
    parameter. The phase space coordinates are assumed to represent barycentric (i.e.
    centred on the Sun) positions and velocities.

    Parameters
    ----------
    x : float or float array
        The x component of the barycentric position vector (in pc or kpc).
    y : float or float array
        The y component of the barycentric position vector (in pc or kpc).
    z : float or float array
        The z component of the barycentric position vector (in pc or kpc).
    vx : float or float array
        The x component of the barycentric velocity vector (in km/s).
    vy : float or float array
        The y component of the barycentric velocity vector (in km/s).
    vz : float or float array
        The z component of the barycentric velocity vector (in km/s).

    Returns
    -------
    phi : float or float array
        The longitude-like angle of the position of the source (radians).
    theta : float or float array
        The latitude-like angle of the position of the source (radians).
    parallax : float or float array
        The parallax of the source (in mas or muas, see notes)
    muphistar : float or float array
        The proper motion in the longitude-like angle, multiplied by cos(theta) (mas/yr
        or muas/yr, see notes)
    mutheta : float or float array
        The proper motion in the latitude-like angle (mas/yr or muas/yr, see notes)
    vrad : float or float array
        The radial velocity (km/s)

    Notes
    -----
    This function has no mechanism to deal with units. The velocity units are always
    assumed to be km/s, and the code is set up such that for positions in pc, the return
    units for the astrometry are radians, milliarcsec, milliarcsec/year and km/s. For
    positions in kpc the return units are: radians, microarcsec, microarcsec/year, and
    km/s.

    The doppler factor :math:`k=1/(1-v_\mathrm{rad}/c)` is not used in the calculations.
    This is not a problem for sources moving at typical velocities of Galactic stars.

    """
    u, phi, theta = cartesian_to_spherical(x, y, z)
    parallax = au_mas_parsec / u
    p, q, r = normal_triad(phi, theta)
    velocitiesArray = np.array([vx, vy, vz])
    if np.isscalar(u):
        muphistar = np.dot(p, velocitiesArray) * parallax / au_km_year_per_sec
        mutheta = np.dot(q, velocitiesArray) * parallax / au_km_year_per_sec
        vrad = np.dot(r, velocitiesArray)
    else:
        muphistar = np.sum(p * velocitiesArray, axis=0) * parallax / au_km_year_per_sec
        mutheta = np.sum(q * velocitiesArray, axis=0) * parallax / au_km_year_per_sec
        vrad = np.sum(r * velocitiesArray, axis=0)

    return phi, theta, parallax, muphistar, mutheta, vrad


def astrometry_to_phase_space(phi, theta, parallax, muphistar, mutheta, vrad):
    r"""
    From the input astrometric parameters calculate the phase space coordinates. The
    output phase space coordinates represent barycentric (i.e. centred on the Sun)
    positions and velocities.

    Parameters
    ----------
    phi : float or float array
        The longitude-like angle of the position of the source (radians).
    theta : float or float array
        The latitude-like angle of the position of the source (radians).
    parallax : float or float array
        The parallax of the source (in mas or muas, see notes)
    muphistar : float or float array
        The proper motion in the longitude-like angle, multiplied by cos(theta) (mas/yr
        or muas/yr, see notes)
    mutheta : float or float array
        The proper motion in the latitude-like angle (mas/yr or muas/yr, see notes)
    vrad : float or float array
        The radial velocity (km/s)

    Returns
    -------
    x : float or float array
        The x component of the barycentric position vector (in pc or kpc).
    y : float or float array
        The y component of the barycentric position vector (in pc or kpc).
    z : float or float array
        The z component of the barycentric position vector (in pc or kpc).
    vx : float or float array
        The x component of the barycentric velocity vector (in km/s).
    vy : float or float array
        The y component of the barycentric velocity vector (in km/s).
    vz : float or float array
        The z component of the barycentric velocity vector (in km/s).

    Raises
    ------
    ValueError
        If any of the input parallaxes is non-positive.

    Notes
    -----
    This function has no mechanism to deal with units. The code is set up such that for
    input astrometry with parallaxes and proper motions in mas and mas/yr, and radial
    velocities in km/s, the phase space coordinates are in pc and km/s. For input
    astrometry with parallaxes and proper motions in muas and muas/yr, and radial
    velocities in km/s, the phase space coordinates are in kpc and km/s. Only positive
    parallaxes are accepted, an exception is thrown if this condition is not met.

    The doppler factor :math:`k=1/(1-v_mathrm{rad}/c)` is not used in the calculations.
    This is not a problem for sources moving at typical velocities of Galactic stars.

    **This function should not be used when the parallaxes have relative uncertainties
    larger than about 20 per cent** (see http://arxiv.org/abs/1507.02105 for example).
    For astrometric data with relatively large parallax errors you should consider doing
    your analysis in the data space and use forward modelling of some kind.
    """
    if np.any(parallax <= 0.0):
        raise ValueError("One or more of the input parallaxes is non-positive")
    x, y, z = spherical_to_cartesian(au_mas_parsec / parallax, phi, theta)
    p, q, r = normal_triad(phi, theta)
    transverseMotionArray = np.array(
        [
            muphistar * au_km_year_per_sec / parallax,
            mutheta * au_km_year_per_sec / parallax,
            vrad,
        ]
    )
    if np.isscalar(parallax):
        velocityArray = np.dot(np.transpose(np.array([p, q, r])), transverseMotionArray)
        vx = velocityArray[0]
        vy = velocityArray[1]
        vz = velocityArray[2]
    else:
        vx = np.zeros_like(parallax)
        vy = np.zeros_like(parallax)
        vz = np.zeros_like(parallax)
        velocityArray = np.einsum(
            "kij,ji->ki", np.stack((p, q, r), axis=-1), transverseMotionArray
        )
        vx = velocityArray[0]
        vy = velocityArray[1]
        vz = velocityArray[2]
    return x, y, z, vx, vy, vz




class EpochPropagation:
    """
    Provides methods for transforming the astrometry and radial velocity of a given
    source to a different epoch.

    The formulae for rigorous epoch transformation can be found on the `Gaia
    documentation pages
    <https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu3ast/sec_cu3ast_intro/ssec_cu3ast_intro_tansforms.html>`_.

    Attributes
    ----------
    mastorad : float
        Numerical factor to convert milliarcseconds ro radians.
    """

    def __init__(self):
        """
        Class constructor/initializer.
        """
        self.mastorad = np.pi / (180 * 3600 * 1000)

    def propagate_astrometry(
        self, phi, theta, parallax, muphistar, mutheta, vrad, t0, t1
    ):
        """
        Propagate the astrometric parameters of a source from the reference epoch t0 to
        the new epoch t1.

        Parameters
        ----------
        phi : float
            Longitude at reference epoch (radians).
        theta : float
            Latitude at reference epoch (radians).
        parallax : float
            Parallax at the reference epoch (mas).
        muphistar : float
            Proper motion in longitude (including np.cos(latitude) term) at reference
            epoch (mas/yr).
        mutheta : float
            Proper motion in latitude at reference epoch (mas/yr).
        vrad : float
            Radial velocity at reference epoch (km/s). Can be set to 0 km/s if not known.
        t0 : float
            Reference epoch (Julian years).
        t1 : float
            New epoch (Julian years).

        Returns
        -------
        phi1, theta1, parallax1, muphistar1, mutheta1, murad1 : float or array
            Astrometric parameters, including the "radial proper motion" (NOT the radial
            velocity), at the new epoch.
        """

        t = t1 - t0
        p0, q0, r0 = normal_triad(phi, theta)

        # Convert input data to units of radians and Julian year. Use ICRS coordinate
        # names internally to avoid errors in translating the formulae to code.
        pmra0 = muphistar * self.mastorad
        pmdec0 = mutheta * self.mastorad
        pmr0 = vrad * parallax / au_km_year_per_sec * self.mastorad
        pmtot0sqr = (muphistar**2 + mutheta**2) * self.mastorad**2

        # Proper motion vector
        pmvec0 = pmra0 * p0 + pmdec0 * q0

        f = (1 + 2 * pmr0 * t + (pmtot0sqr + pmr0**2) * t**2) ** (-0.5)
        u = (r0 * (1 + pmr0 * t) + pmvec0 * t) * f

        _, phi1, theta1 = cartesian_to_spherical(u[0], u[1], u[2])
        parallax1 = parallax * f
        pmr1 = (pmr0 + (pmtot0sqr + pmr0**2) * t) * f**2
        pmvec1 = (pmvec0 * (1 + pmr0 * t) - r0 * pmtot0sqr * t) * f**3
        p1, q1, r1 = normal_triad(phi1, theta1)
        muphistar1 = np.sum(p1 * pmvec1 / self.mastorad, axis=0)
        mutheta1 = np.sum(q1 * pmvec1 / self.mastorad, axis=0)
        murad1 = pmr1 / self.mastorad

        return phi1, theta1, parallax1, muphistar1, mutheta1, murad1

    def propagate_pos(self, phi, theta, parallax, muphistar, mutheta, vrad, t0, t1):
        """
        Propagate the position of a source from the reference epoch t0 to the new epoch
        t1.

        Parameters
        ----------
        phi : float
            Longitude at reference epoch (radians).
        theta : float
            Latitude at reference epoch (radians).
        parallax : float
            Parallax at the reference epoch (mas).
        muphistar : float
            Proper motion in longitude (including np.cos(latitude) term) at reference
            epoch (mas/yr).
        mutheta : float
            Proper motion in latitude at reference epoch (mas/yr).
        vrad : float
            Radial velocity at reference epoch (km/s). Can be set to 0 km/s if not known.
        t0 : float
            Reference epoch (Julian years).
        t1 : float
            New epoch (Julian years).

        Returns
        -------
        phi, theta : float
            Coordinates the new epoch (in radians)
        """
        (
            phi1,
            theta1,
            _,
            _,
            _,
            _,
        ) = self.propagate_astrometry(
            phi, theta, parallax, muphistar, mutheta, vrad, t0, t1
        )
        return phi1, theta1

    def propagate_astrometry_and_covariance_matrix(self, a0, c0, t0, t1):
        """
        Propagate the astrometric parameters iand radial proper motion, as well as the
        corresponding covariance matrix, from epoch t0 to epoch t1.

        Code based on the Hipparcos Fortran implementation by Lennart Lindegren.

        Parameters
        ----------
        a0 : array_like
            6-element vector: (phi, theta, parallax, muphistar, mutheta, vrad) in units
            of (radians, radians, mas, mas/yr, mas/yr, km/s). Shape of a should be (6,)
            or (6,N), with N the number of sources for which the astrometric parameters
            are provided. The value of vrad can be set to 0 km/s if the radial velocity is not know. An appropriate uncertainty should then be provided (see below).
        c0 : array_like
            Covariance matrix stored in a 6x6 element array. This can be constructed
            from the columns listed in the Gaia catalogue. The units are [mas^2,
            mas^2/yr, mas^2/yr^2] for the various elements. Note that the elements in
            the 6th row and column should be:
            c[6,i] = c[i,6] = c[i,3] * vrad / auKmYearPerSec
            for i = 1,..,5 and
            c[6,6] = c[3,3] * (vrad^2+vrad_error^2) / auKmYearPerSec^2 +
            (parallax*vrad_error/auKmYearPerSec)^2
            The shape of c0 should be (6,6) or (N,6,6). If the radial velocity is not know the uncertainty on the radial velocity (vrad_error) should be set to the velocity dispersion of the population the source is drawn from.
            To construct the covariance matrix the convenience function :py:meth:`pygaia.utils.construct_covariance_matrix` is provided.
        t0 : float
            Reference epoch (Julian years).
        t1 : float
            New epoch (Julian years).

        Returns
        -------
        a1, c1 : array
            Astrometric parameters, including the "radial proper motion" (NOT the radial
            velocity), and the covariance matrix at the new epoch as a 2D matrix with
            the new variances on the diagonal and the covariance in the off-diagonal
            elements. Matrix shapes are (6, N) and (N, 6, 6).
        """

        zero, one, two, three = 0, 1, 2, 3
        tau = t1 - t0

        # Calculate the normal triad [p0 q0 r0] at t0
        p0, q0, r0 = normal_triad(a0[0], a0[1])

        # Convert to internal units (radians, Julian year)
        par0 = a0[2] * self.mastorad
        pma0 = a0[3] * self.mastorad
        pmd0 = a0[4] * self.mastorad
        pmr0 = a0[5] * a0[2] / au_km_year_per_sec * self.mastorad

        # Proper motion vector
        pmvec0 = pma0 * p0 + pmd0 * q0

        # Auxiliary quantities
        tau2 = tau * tau
        pm02 = pma0**2 + pmd0**2
        w = one + pmr0 * tau
        f2 = one / (one + two * pmr0 * tau + (pm02 + pmr0**2) * tau2)
        f = np.sqrt(f2)
        f3 = f2 * f
        f4 = f2 * f2

        # Position vector and parallax at t1
        u = (r0 * w + pmvec0 * tau) * f
        _, ra, dec = cartesian_to_spherical(u[0], u[1], u[2])
        par = par0 * f

        # Proper motion vector and radial proper motion at t1
        pmvec = (pmvec0 * (one + pmr0 * tau) - r0 * pm02 * tau) * f3
        pmr = (pmr0 + (pm02 + pmr0**2) * tau) * f2

        # Normal triad at t1
        p, q, r = normal_triad(ra, dec)

        # Convert parameters at t1 to external units (mas, Julian year)
        pma = np.sum(p * pmvec, axis=0)
        pmd = np.sum(q * pmvec, axis=0)

        a = np.zeros_like(a0)
        a[0] = ra
        a[1] = dec
        a[2] = par / self.mastorad
        a[3] = pma / self.mastorad
        a[4] = pmd / self.mastorad
        a[5] = pmr / self.mastorad

        # Auxiliary quantities for the partial derivatives

        pmz = pmvec0 * f - three * pmvec * w
        pp0 = np.sum(p * p0, axis=0)
        pq0 = np.sum(p * q0, axis=0)
        pr0 = np.sum(p * r0, axis=0)
        qp0 = np.sum(q * p0, axis=0)
        qq0 = np.sum(q * q0, axis=0)
        qr0 = np.sum(q * r0, axis=0)
        ppmz = np.sum(p * pmz, axis=0)
        qpmz = np.sum(q * pmz, axis=0)

        jacobian = np.zeros_like(c0)
        if c0.ndim == 2:
            jacobian = jacobian[np.newaxis, :, :]

        # Partial derivatives
        jacobian[:, 0, 0] = pp0 * w * f - pr0 * pma0 * tau * f
        jacobian[:, 0, 1] = pq0 * w * f - pr0 * pmd0 * tau * f
        jacobian[:, 0, 2] = zero
        jacobian[:, 0, 3] = pp0 * tau * f
        jacobian[:, 0, 4] = pq0 * tau * f
        jacobian[:, 0, 5] = -pma * tau2

        jacobian[:, 1, 0] = qp0 * w * f - qr0 * pma0 * tau * f
        jacobian[:, 1, 1] = qq0 * w * f - qr0 * pmd0 * tau * f
        jacobian[:, 1, 2] = zero
        jacobian[:, 1, 3] = qp0 * tau * f
        jacobian[:, 1, 4] = qq0 * tau * f
        jacobian[:, 1, 5] = -pmd * tau2

        jacobian[:, 2, 0] = zero
        jacobian[:, 2, 1] = zero
        jacobian[:, 2, 2] = f
        jacobian[:, 2, 3] = -par * pma0 * tau2 * f2
        jacobian[:, 2, 4] = -par * pmd0 * tau2 * f2
        jacobian[:, 2, 5] = -par * w * tau * f2

        jacobian[:, 3, 0] = -pp0 * pm02 * tau * f3 - pr0 * pma0 * w * f3
        jacobian[:, 3, 1] = -pq0 * pm02 * tau * f3 - pr0 * pmd0 * w * f3
        jacobian[:, 3, 2] = zero
        jacobian[:, 3, 3] = (
            pp0 * w * f3 - two * pr0 * pma0 * tau * f3 - three * pma * pma0 * tau2 * f2
        )
        jacobian[:, 3, 4] = (
            pq0 * w * f3 - two * pr0 * pmd0 * tau * f3 - three * pma * pmd0 * tau2 * f2
        )
        jacobian[:, 3, 5] = ppmz * tau * f2

        jacobian[:, 4, 0] = -qp0 * pm02 * tau * f3 - qr0 * pma0 * w * f3
        jacobian[:, 4, 1] = -qq0 * pm02 * tau * f3 - qr0 * pmd0 * w * f3
        jacobian[:, 4, 2] = zero
        jacobian[:, 4, 3] = (
            qp0 * w * f3 - two * qr0 * pma0 * tau * f3 - three * pmd * pma0 * tau2 * f2
        )
        jacobian[:, 4, 4] = (
            qq0 * w * f3 - two * qr0 * pmd0 * tau * f3 - three * pmd * pmd0 * tau2 * f2
        )
        jacobian[:, 4, 5] = qpmz * tau * f2

        jacobian[:, 5, 0] = zero
        jacobian[:, 5, 1] = zero
        jacobian[:, 5, 2] = zero
        jacobian[:, 5, 3] = two * pma0 * w * tau * f4
        jacobian[:, 5, 4] = two * pmd0 * w * tau * f4
        jacobian[:, 5, 5] = (w**2 - pm02 * tau2) * f4

        jacobian_transposed = np.zeros_like(jacobian)
        for i in range(jacobian.shape[0]):
            jacobian_transposed[i] = jacobian[i].T

        if c0.ndim == 2:
            c = np.matmul(
                jacobian, np.matmul(c0[np.newaxis, :, :], jacobian_transposed)
            )
        else:
            c = np.matmul(jacobian, np.matmul(c0, jacobian_transposed))

        return a, np.squeeze(c)
