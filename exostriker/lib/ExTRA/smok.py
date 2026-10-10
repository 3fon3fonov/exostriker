"""Scaled Modelling of Kinematics (SMOK).

SMOK coordinates are a local Cartesian representation of a scaled barycentric
position. Angles are in radians and vectors use the first axis for Cartesian
components, matching :func:`ExTRA.vectorastrometry.normal_triad`.
"""

import numpy as np

from .vectorastrometry import cartesian_to_spherical, normal_triad
from .useful import au_km_year_per_sec,pc_in_km,julian_year_seconds

__all__ = [
	"sss_to_smok",
	"combine_sss_inverse_covariance",
	"smok_to_sss",
	"smok_coordinates",
	"smok_vector",
	"change_comparison_point",
	"propagate_smok_to_spherical",
	"change_comparison_point_6par",
	"smok_at_new_comparison_point"
]


def combine_sss_inverse_covariance(sss_1, covariance_1, sss_2, covariance_2):
	"""Combine two aligned SSS solutions by inverse-covariance weighting.

	The two solutions must describe the same five parameters in the same
	coordinate system and at the same epoch. Propagate and transform HIP
	before calling this function if it is not already aligned with Gaia.

	Parameters
	----------
	sss_1, sss_2 : array_like, shape (5,)
		``(alpha, delta, parallax, mu_alpha_star, mu_delta)`` vectors.
	covariance_1, covariance_2 : array_like, shape (5, 5)
		Covariance matrices in the same parameter ordering and units as the
		corresponding SSS vectors.

	Returns
	-------
	common_sss, common_covariance : ndarray
		The inverse-covariance weighted solution and its covariance matrix.
	"""
	sss_1 = np.asarray(sss_1, dtype=float)
	sss_2 = np.asarray(sss_2, dtype=float)
	covariance_1 = np.asarray(covariance_1, dtype=float)
	covariance_2 = np.asarray(covariance_2, dtype=float)

	weight_1 = np.linalg.pinv(covariance_1)
	weight_2 = np.linalg.pinv(covariance_2)
	common_covariance = np.linalg.pinv(weight_1 + weight_2)
	common_sss = common_covariance @ (weight_1 @ sss_1 + weight_2 @ sss_2)
	return common_sss, common_covariance


def sss_to_smok(
	sss, alpha_c, delta_c, vrad=0.0, parallax_c=None, angles_in_degrees=True
):
	"""Convert a five-parameter SSS into a complete SMOK state.

	Parameters
	----------
	sss : array_like
		Five standard astrometric parameters ``(alpha, delta, parallax,
		mu_alpha_star, mu_delta)``. Parallax is in mas and proper motions
		are in mas/year. Angles are in degrees by default.
	alpha_c, delta_c : float
		Fixed comparison-point right ascension and declination.
	vrad : float, optional
		Radial velocity in km/s. The default value of zero omits perspective
		motion.
	parallax_c : float, optional
		Parallax in mas defining the SMOK scale. Defaults to the SSS parallax.
	angles_in_degrees : bool, optional
		Interpret SSS and comparison-point angles as degrees when true.

	Returns
	-------
	 ndarray
		The six-component state ``(a, d, r, adot, ddot, rdot)``. The first
		three values are dimensionless and the derivatives are per Julian year.
	"""
	sss = np.asarray(sss)
	
	alpha, delta = sss[:2]
	parallax = sss[2]
	if angles_in_degrees:
		alpha, delta, alpha_c, delta_c = np.radians(
			[alpha, delta, alpha_c, delta_c]
		)
	if parallax_c is None:
		parallax_c = parallax

	position = np.array(
		[
			np.cos(delta) * np.cos(alpha),
			np.cos(delta) * np.sin(alpha),
			np.sin(delta),
		]
	)
	p_star, q_star, radial_star = normal_triad(alpha, delta)
	position_scale = parallax_c / parallax
	mu_radial = vrad * parallax * np.pi / (180 * 3600 * 1000) / au_km_year_per_sec
	motion = position_scale * (
		sss[3] * np.pi / (180 * 3600 * 1000) * p_star
		+ sss[4] * np.pi / (180 * 3600 * 1000) * q_star
		+ mu_radial * radial_star
	)
	return np.concatenate((smok_coordinates(position, alpha_c, delta_c), smok_coordinates(motion, alpha_c, delta_c)))


def smok_to_sss(
	a, d, r, adot, ddot, rdot, alpha_c, delta_c, parallax_c,
	angles_in_degrees=True
):
	"""Convert a complete SMOK state to five standard astrometric parameters.

	Parameters
	----------
	a, d, r : float
		SMOK position coordinates.
adot, ddot, rdot : float
		SMOK coordinate derivatives per Julian year.
alpha_c, delta_c : float
		Comparison-point coordinates. These use degrees when
		``angles_in_degrees`` is true, otherwise radians.
parallax_c : float
		Parallax in mas defining the SMOK scale.

	Returns
	-------
	 ndarray
		``(alpha, delta, parallax, mu_alpha_star, mu_delta)``. The angles use
		degrees by default, parallax is in mas, and proper motions are in
		mas/year.
	"""
	if angles_in_degrees:
		alpha_c, delta_c = np.radians([alpha_c, delta_c])

	position = np.asarray(smok_vector(a, d, r, alpha_c, delta_c), dtype=float)
	motion = np.asarray(smok_vector(adot, ddot, rdot, alpha_c, delta_c), dtype=float)
	position_length = np.linalg.norm(position)

	unit_position = position / position_length
	_, alpha, delta = cartesian_to_spherical(*unit_position)
	p_star, q_star, _ = normal_triad(alpha, delta)
	rad_to_mas = 180 * 3600 * 1000 / np.pi

	distance_c_km = pc_in_km * 1000.0 / parallax_c
	#rv in kms
	radial_velocity = (
		np.dot(unit_position, motion) * distance_c_km / julian_year_seconds
	)  

	if angles_in_degrees:
		alpha, delta = np.degrees([alpha, delta])

	return np.array([
		alpha,
		delta,
		parallax_c / position_length,
		rad_to_mas * np.dot(p_star, motion) / position_length,
		rad_to_mas * np.dot(q_star, motion) / position_length,
		radial_velocity
	])
def smok_coordinates(vector, alpha_c, delta_c):
	"""Project scaled Cartesian vectors onto a SMOK normal triad.

	Parameters
	----------
	vector : array_like
		Cartesian vector with shape ``(3,)`` or ``(3, N)``.
	alpha_c, delta_c : float
		Right ascension and declination of the fixed comparison point, in
		radians. This low-level helper does not convert angle units.

	Returns
	-------
	a, d, r : ndarray or float
		SMOK coordinates, corresponding to the components along the local
		east, north, and radial axes.
	"""
	vector = np.asarray(vector)
	if vector.ndim == 0 or vector.shape[0] != 3:
		raise ValueError("vector must have shape (3,) or (3, N)")

	p, q, radial = normal_triad(alpha_c, delta_c)
	component_shape = (3,) + (1,) * (vector.ndim - 1)
	p = np.reshape(p, component_shape)
	q = np.reshape(q, component_shape)
	radial = np.reshape(radial, component_shape)
	return (
		np.sum(vector * p, axis=0),
		np.sum(vector * q, axis=0),
		np.sum(vector * radial, axis=0),
	)


def smok_vector(a, d, r, alpha_c, delta_c):
	"""Reconstruct scaled Cartesian vectors from SMOK coordinates.

	The comparison-point angles must be in radians. The returned vector uses
	the first axis for Cartesian components and can have shape ``(3,)`` or
	``(3, N)``.
	"""
	p, q, radial = normal_triad(alpha_c, delta_c)
	component_shape = (3,) + (1,) * np.ndim(a)
	p = np.reshape(p, component_shape)
	q = np.reshape(q, component_shape)
	radial = np.reshape(radial, component_shape)
	return p * a + q * d + radial * r


def change_comparison_point(a, d, r, alpha_c, delta_c, new_alpha_c, new_delta_c):
	"""Express SMOK coordinates relative to a new comparison point.

	All comparison-point angles must be in radians. This changes only the
	coordinate representation, not the underlying scaled vector.
	"""
	vector = smok_vector(a, d, r, alpha_c, delta_c)
	return smok_coordinates(vector, new_alpha_c, new_delta_c)

def change_comparison_point_6par(
	a, d, r, adot, ddot, rdot,
	alpha_c, delta_c, new_alpha_c, new_delta_c,
):
	"""Express a full SMOK state (position + velocity) relative to a new
	comparison point.

	All comparison-point angles must be in radians. This changes only the
	coordinate representation, not the underlying scaled vectors.

	Returns
	-------
	ndarray
		(a, d, r, adot, ddot, rdot) in the new comparison-point frame.
	"""
	a2, d2, r2 = change_comparison_point(a, d, r, alpha_c, delta_c, new_alpha_c, new_delta_c)
	adot2, ddot2, rdot2 = change_comparison_point(
		adot, ddot, rdot, alpha_c, delta_c, new_alpha_c, new_delta_c
	)
	return np.array([a2, d2, r2, adot2, ddot2, rdot2])

def smok_at_new_comparison_point(
	a, d, r, adot, ddot, rdot,
	alpha_c, delta_c, parallax_c,
	new_alpha_c, new_delta_c,
	angles_in_degrees=True,
):
	"""Re-express a SMOK state at a new comparison point, and return the
	standard astrometric parameters there.

	Parameters
	----------
	a, d, r, adot, ddot, rdot : float
		SMOK state relative to the old comparison point (alpha_c, delta_c).
	alpha_c, delta_c : float
		Old comparison-point coordinates.
	parallax_c : float
		Parallax in mas defining the SMOK scale (unchanged by a comparison-
		point switch, since it's a pure rotation of the same vector).
	new_alpha_c, new_delta_c : float
		New comparison-point coordinates, same angle convention as
		alpha_c/delta_c (degrees if angles_in_degrees, else radians).
	angles_in_degrees : bool
		If True, all angle inputs/outputs are degrees; internally converted
		to radians for the rotation.

	Returns
	-------
	ndarray
		(alpha, delta, parallax, mu_alpha_star, mu_delta, v_r) at the star's
		propagated/rotated position, expressed relative to the new
		comparison point. Angles in degrees by default, v_r in m/s.
	"""
	if angles_in_degrees:
		alpha_c_rad, delta_c_rad, new_alpha_c_rad, new_delta_c_rad = np.radians(
			[alpha_c, delta_c, new_alpha_c, new_delta_c]
		)
	else:
		alpha_c_rad, delta_c_rad = alpha_c, delta_c
		new_alpha_c_rad, new_delta_c_rad = new_alpha_c, new_delta_c

	a2, d2, r2, adot2, ddot2, rdot2 = change_comparison_point_6par(
		a, d, r, adot, ddot, rdot,
		alpha_c_rad, delta_c_rad, new_alpha_c_rad, new_delta_c_rad,
	)

	return smok_to_sss(
		a2, d2, r2, adot2, ddot2, rdot2,
		new_alpha_c_rad, new_delta_c_rad, parallax_c,
		angles_in_degrees=False,  # already converted above
	)


def propagate_smok_to_spherical(
	a, d, r, adot, ddot, rdot,
	alpha_c, delta_c, parallax_c,
	tau, angles_in_degrees=True,
):
	"""Propagate a SMOK state to one or more epochs and return sky position.

	Parameters
	----------
	a, d, r, adot, ddot, rdot : float
		SMOK position and velocity, relative to (alpha_c, delta_c).
	alpha_c, delta_c : float
		Comparison-point coordinates (degrees if angles_in_degrees, else
		radians).
	parallax_c : float
		Parallax in mas defining the SMOK scale.
	tau : float or array_like
		Propagation interval(s), in Julian years.

	Returns
	-------
	ndarray
		Shape (3,) if tau is scalar, or (3, N) if tau has N elements:
		(alpha, delta, parallax) at each epoch. Degrees by default.
	"""
	scalar_input = np.isscalar(tau) or np.ndim(tau) == 0

	if angles_in_degrees:
		alpha_c_rad, delta_c_rad = np.radians([alpha_c, delta_c])
	else:
		alpha_c_rad, delta_c_rad = alpha_c, delta_c

	tau = np.atleast_1d(np.asarray(tau, dtype=float))

	p_c, q_c, r_c = normal_triad(alpha_c_rad, delta_c_rad)

	a_p = a + adot * tau
	d_p = d + ddot * tau
	r_p = r + rdot * tau

	position = (
		a_p[np.newaxis, :] * p_c[:, np.newaxis]
		+ d_p[np.newaxis, :] * q_c[:, np.newaxis]
		+ r_p[np.newaxis, :] * r_c[:, np.newaxis]
	)

	position_length = np.linalg.norm(position, axis=0)
	unit_position = position / position_length

	_, alpha, delta = cartesian_to_spherical(*unit_position)
	parallax = parallax_c / position_length

	if angles_in_degrees:
		alpha, delta = np.degrees([alpha, delta])

	result = np.array([alpha, delta, parallax])
	return result[:, 0] if scalar_input else result