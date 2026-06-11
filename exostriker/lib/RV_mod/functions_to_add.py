def calculate_ttv_period(P_in, P_out):
    """
    Calculate the TTV period based on the periods of two near-resonant planets.
    
    Parameters:
    P_in (float): Orbital period of the inner planet.
    P_out (float): Orbital period of the outer planet.
    
    Returns:
    float: TTV period (P_TTV)
    """
    # Estimate the resonance order (j) based on the ratio of periods
    j = round(P_out / P_in)

 
    # Calculate the TTV period
    P_TTV = abs((j - 1) / P_in - j / P_out)**-1
    
    return P_TTV
    
    
import numpy as np


def possible_perturber_periods(P_planet, P_super, planet_is="inner",
                               j_min=2, j_max=10, tol=None):
    """
    Estimate possible perturber periods from a known planet period and
    an observed TTV super-period, assuming first-order j:(j-1) resonances.

    Parameters
    ----------
    P_planet : float
        Orbital period of the known planet.

    P_super : float
        Observed TTV super-period.

    planet_is : str
        Either "inner" or "outer".
        Use "inner" if the known planet is the inner planet.
        Use "outer" if the known planet is the outer planet.

    j_min : int
        Minimum resonance index j.
        j=2 corresponds to 2:1, j=3 to 3:2, etc.

    j_max : int
        Maximum resonance index j to test.

    tol : float or None
        Optional maximum fractional distance from exact resonance.
        Example: tol=0.05 keeps only solutions within 5 percent of exact resonance.

    Returns
    -------
    list of dict
        Possible perturber periods and resonance information.
    """

    if P_planet <= 0:
        raise ValueError("P_planet must be positive.")

    if P_super <= 0:
        raise ValueError("P_super must be positive.")

    if planet_is not in ["inner", "outer"]:
        raise ValueError("planet_is must be either 'inner' or 'outer'.")

    s = 1.0 / P_super
    solutions = []

    for j in range(j_min, j_max + 1):

        if planet_is == "inner":
            P_inner = P_planet

            # From:
            # s = abs((j - 1)/P_inner - j/P_outer)
            #
            # Therefore:
            # j/P_outer = (j - 1)/P_inner +/- s

            for sign in [+1, -1]:
                denom = (j - 1) / P_inner + sign * s

                if denom <= 0:
                    continue

                P_outer = j / denom

                if P_outer <= P_inner:
                    continue

                period_ratio = P_outer / P_inner
                resonance_ratio = j / (j - 1)
                frac_offset = (period_ratio / resonance_ratio) - 1.0

                if tol is not None and abs(frac_offset) > tol:
                    continue

                solutions.append({
                    "known_planet": "inner",
                    "perturber": "outer",
                    "resonance": f"{j}:{j-1}",
                    "P_known": P_inner,
                    "P_perturber": P_outer,
                    "period_ratio": period_ratio,
                    "exact_resonance_ratio": resonance_ratio,
                    "fractional_offset": frac_offset,
                    "wide_or_narrow": "wide" if frac_offset > 0 else "narrow"
                })

        elif planet_is == "outer":
            P_outer = P_planet

            # From:
            # s = abs((j - 1)/P_inner - j/P_outer)
            #
            # Therefore:
            # (j - 1)/P_inner = j/P_outer +/- s

            for sign in [+1, -1]:
                denom = j / P_outer + sign * s

                if denom <= 0:
                    continue

                P_inner = (j - 1) / denom

                if P_inner >= P_outer:
                    continue

                period_ratio = P_outer / P_inner
                resonance_ratio = j / (j - 1)
                frac_offset = (period_ratio / resonance_ratio) - 1.0

                if tol is not None and abs(frac_offset) > tol:
                    continue

                solutions.append({
                    "known_planet": "outer",
                    "perturber": "inner",
                    "resonance": f"{j}:{j-1}",
                    "P_known": P_outer,
                    "P_perturber": P_inner,
                    "period_ratio": period_ratio,
                    "exact_resonance_ratio": resonance_ratio,
                    "fractional_offset": frac_offset,
                    "wide_or_narrow": "wide" if frac_offset > 0 else "narrow"
                })

    return solutions    
    
def a_from_K_in_mas(parallax,P,e,i,K):
    """
    Calculate the semi-major axis 'a' given K, parallax, period P, eccentricity e, and inclination i.

    Parameters:
    K : float
        Radial velocity amplitude (m/s)
    parallax : float
        Parallax in miliarcseconds
    P : float
        Orbital period in days
    e : float
        Eccentricity
    i : float
        Inclination in degerees

    Returns:
    a : float
        Semi-major axis in mas
"""    
    i = np.radians(i)

    a = (K * P *86400.0 * np.sqrt(1.0-e**2)* parallax) / (2*np.pi*np.sin(i)*1.496e11)

    return a    
