from __future__ import division, print_function
import numpy
from os import path
import transitleastsquares_ES.tls_constants as tls_constants


def catalog_info_KIC(KIC_ID):
    """Takes KIC_ID, returns stellar information from online catalog using Vizier"""

    if type(KIC_ID) is not int:
        raise TypeError('KIC_ID ID must be of type "int"')
    try:
        from astroquery.vizier import Vizier
    except:
        raise ImportError("Package astroquery required but failed to import")

    columns = ["Teff", "log(g)", "Rad", "E_Rad", "e_Rad", "Mass", "E_Mass", "e_Mass"]
    catalog = "J/ApJS/229/30/catalog"
    result = (
        Vizier(columns=columns)
        .query_constraints(KIC=KIC_ID, catalog=catalog)[0]
        .as_array()
    )
    Teff = result[0][0]
    logg = result[0][1]
    radius = result[0][2]
    radius_max = result[0][3]
    radius_min = result[0][4]
    mass = result[0][5]
    mass_max = result[0][6]
    mass_min = result[0][7]
    return Teff, logg, radius, radius_min, radius_max, mass, mass_min, mass_max


def catalog_info_EPIC(EPIC_ID):
    """Takes EPIC_ID, returns stellar information from online catalog using Vizier"""

    try:
        from astroquery.vizier import Vizier
    except:
        raise ImportError("Package astroquery required but failed to import")
    if type(EPIC_ID) is not int:
        raise TypeError('EPIC_ID ID must be of type "int"')
    if (EPIC_ID < 201000001) or (EPIC_ID > 251813738):
        raise TypeError("EPIC_ID ID must be in range 201000001 to 251813738")
    columns = ["Teff", "logg", "Rad", "E_Rad", "e_Rad", "Mass", "E_Mass", "e_Mass"]
    catalog = "IV/34/epic"
    result = (
        Vizier(columns=columns)
        .query_constraints(ID=EPIC_ID, catalog=catalog)[0]
        .as_array()
    )
    Teff = result[0][0]
    logg = result[0][1]
    radius = result[0][2]
    radius_max = result[0][3]
    radius_min = result[0][4]
    mass = result[0][5]
    mass_max = result[0][6]
    mass_min = result[0][7]
    return Teff, logg, radius, radius_min, radius_max, mass, mass_min, mass_max


def catalog_info_TIC(TIC_ID):
    """Takes TIC_ID, returns stellar information from online catalog using Vizier"""
    if type(TIC_ID) is not int:
        raise TypeError('TIC_ID ID must be of type "int"')
    try:
        from astroquery.mast import Catalogs
    except:
        raise ImportError("Package astroquery required but failed to import")

    result = Catalogs.query_criteria(catalog="Tic", ID=TIC_ID).as_array()
    Teff = result[0][64]
    logg = result[0][66]
    radius = result[0][70]
    radius_max = result[0][71]
    radius_min = result[0][71]
    mass = result[0][72]
    mass_max = result[0][73]
    mass_min = result[0][73]
    return Teff, logg, radius, radius_min, radius_max, mass, mass_min, mass_max


def catalog_info(EPIC_ID=None, TIC_ID=None, KIC_ID=None):
    """Takes EPIC ID, returns limb darkening parameters u (linear) and
        a,b (quadratic), and stellar parameters. Values are pulled for minimum
        absolute deviation between given/catalog Teff and logg. Data are from:
        - K2 Ecliptic Plane Input Catalog, Huber+ 2016, 2016ApJS..224....2H
        - New limb-darkening coefficients, Claret+ 2012, 2013,
          2012A&A...546A..14C, 2013A&A...552A..16C"""

    if (EPIC_ID is None) and (TIC_ID is None) and (KIC_ID is None):
        raise ValueError("No ID was given")
    if (EPIC_ID is not None) and (TIC_ID is not None):
        raise ValueError("Only one ID allowed")
    if (EPIC_ID is not None) and (KIC_ID is not None):
        raise ValueError("Only one ID allowed")
    if (TIC_ID is not None) and (KIC_ID is not None):
        raise ValueError("Only one ID allowed")

    # KOI CASE (Kepler K1)
    if KIC_ID is not None:
        Teff, logg, radius, radius_min, radius_max, mass, mass_min, mass_max = catalog_info_KIC(
            KIC_ID
        )

    # EPIC CASE (Kepler K2)
    if EPIC_ID is not None:
        Teff, logg, radius, radius_min, radius_max, mass, mass_min, mass_max = catalog_info_EPIC(
            EPIC_ID
        )

    # TESS CASE
    if TIC_ID is not None:
        Teff, logg, radius, radius_min, radius_max, mass, mass_min, mass_max = catalog_info_TIC(
            TIC_ID
        )
        ld = numpy.genfromtxt(
            path.join(tls_constants.resources_dir, "ld_claret_tess.csv"),
            skip_header=1,
            delimiter=",",
            dtype="f8, int32, f8, f8",
            names=["logg", "Teff", "a", "b"],
        )
    else:  # Limb darkening is the same for K1 (KIC) and K2 (EPIC)
        ld = numpy.genfromtxt(
            path.join(tls_constants.resources_dir, "JAA546A14limb1-4.csv"),
            skip_header=1,
            delimiter=",",
            dtype="f8, int32, f8, f8, f8",
            names=["logg", "Teff", "u", "a", "b"],
        )

    if logg is None:
        logg = 4
        warnings.warn("No logg in catalog. Proceeding with logg=4")

    if Teff is None:
        Teff = 6000
        warnings.warn("No Teff in catalog. Proceeding with Teff=6000")

    """From here on, K2 and TESS catalogs work the same:
        - Take Teff from star catalog and find nearest entry in LD catalog
        - Same for logg, but only for the Teff values returned before
        - Return stellar parameters and best-match LD
    """
    nearest_Teff = ld["Teff"][(numpy.abs(ld["Teff"] - Teff)).argmin()]
    idx_all_Teffs = numpy.where(ld["Teff"] == nearest_Teff)
    relevant_lds = numpy.copy(ld[idx_all_Teffs])
    idx_nearest = numpy.abs(relevant_lds["logg"] - logg).argmin()
    a = relevant_lds["a"][idx_nearest]
    b = relevant_lds["b"][idx_nearest]

    mass = numpy.array(mass)
    mass_min = numpy.array(mass_min)
    mass_max = numpy.array(mass_max)
    radius = numpy.array(radius)
    radius_min = numpy.array(radius_min)
    radius_max = numpy.array(radius_max)

    if mass == 0.0:
        mass = numpy.nan
    if mass_min == 0.0:
        mass_min = numpy.nan
    if mass_max == 0.0:
        mass_max = numpy.nan
    if radius == 0.0:
        radius = numpy.nan
    if radius_min == 0.0:
        radius_min = numpy.nan
    if radius_max == 0.0:
        radius_max = numpy.nan

    return ((a, b), mass, mass_min, mass_max, radius, radius_min, radius_max)
