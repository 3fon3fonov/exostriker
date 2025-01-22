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
