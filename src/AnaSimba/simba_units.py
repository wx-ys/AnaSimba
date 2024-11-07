'''
units and descriptions of parameters
'''

from pynbody import units
import numpy as np


# Define common units used in simba simulations, run with gizmo
UnitLength = units.kpc / units.h
UnitMass = 1e10 * units.Msol / units.h
UnitVel = units.km / units.s
UnitTime = (0.978 * units.Gyr / units.h)
UnitMassdTime = UnitMass / UnitTime
UnitInternalEnergy = UnitVel**2
UnitDensity = UnitMass/(UnitLength**3)


UnitComvingLength = units.a * UnitLength
UnitComvingDensity = UnitDensity/ (units.a**3)
UnitComvingVel = UnitVel * units.a**(1,2)

UnitPressure = (UnitMass / UnitLength) * (units.km / units.s / units.kpc) ** 2
UnitNo = units.no_unit

# Define parameters that will not be converted in physical_units()
# global NotneedtransGCPa
NotneedtransGCPa = [ ]


def halo_pa_name(
    field: str,
) -> str:
    """
    This function modifies the name of a halo parameter to a custom-defined name.

    Parameters:
    -----------
    field : str
        The standard name of the halo parameter as it is typically used in the data.

    Returns:
    --------
    str
        The custom name corresponding to the input parameter. If the input parameter
        does not match any predefined names, the function returns the input parameter
        name unchanged.

    Notes:
    ------
    The function uses a dictionary `Matchfield` to map standard halo parameter names to
    their custom names. Currently, `Matchfield` is empty. If the input `field` is found
    in `Matchfield`, the corresponding custom name is returned. Otherwise, the function
    returns the original `field` name.
    """
    Matchfield = {}
    if field in Matchfield:
        return Matchfield[field]
    else:
        return field


def galaxy_pa_name(
    field: str,
) -> str:
    """
    This function modifies the name of a galaxy parameter to a custom-defined name.

    Parameters:
    -----------
    field : str
        The standard name of the halo parameter as it is typically used in the data.

    Returns:
    --------
    str
        The custom name corresponding to the input parameter. If the input parameter
        does not match any predefined names, the function returns the input parameter
        name unchanged.

    Notes:
    ------
    The function uses a dictionary `Matchfield` to map standard galaxy parameter names to
    their custom names. Currently, `Matchfield` is empty. If the input `field` is found
    in `Matchfield`, the corresponding custom name is returned. Otherwise, the function
    returns the original `field` name.
    """
    Matchfield = {}
    if field in Matchfield:
        return Matchfield[field]
    else:
        return field


def snapshot_pa_name(
    field: str,
) -> str:
    """
    This function modifies the name of a particle parameter to a custom-defined name.

    Parameters:
    -----------
    field : str
        The standard name of the parameter as it is typically used in the data.

    Returns:
    --------
    str
        The custom name corresponding to the input parameter. If the input parameter
        does not match any predefined names, the function returns the input parameter
        name unchanged.

    Notes:
    ------
    The function uses a dictionary `Matchfield` to map standard parameter names to
    their custom names. If the input `field` is found in `Matchfield`, the corresponding
    custom name is returned. Otherwise, the function returns the original `field` name.
    """
    Matchfield = {
        'Coordinates': 'pos',
        'Density': 'rho',
        'ParticleIDs': 'iord',
        'Potential': 'pot',
        'Masses': 'mass',
        'Velocities': 'vel',
        'StellarFormationTime': 'aform',
        'Metallicity': 'metals',
        'InternalEnergy': 'u',
        'StarFormationRate': 'sfr',
    }
    if field in Matchfield:
        return Matchfield[field]
    else:
        return field


def snapshot_units(
    field: str,
) -> units.Unit:
    """
    This function provides the unit corresponding to a given particle parameter.

    Parameters:
    -----------
    field : str
        The name of the particle parameter for which the unit is requested.

    Returns:
    --------
    unit
        The unit associated with the input parameter. If the parameter is not
        defined in `Matchfieldunits`, the function will raise a KeyError.

    """
    Matchfieldunits = {
        'Coordinates': UnitComvingLength,
        'Velocities': units.km * units.a ** (1, 2) / units.s,
        'ParticleIDs': UnitNo,
        'Masses': (UnitMass),
        'InternalEnergy': (UnitVel) ** 2,
        'Density': (UnitMass) / (UnitComvingLength) ** 3,
        'SmoothingLength': UnitComvingLength,
        'ElectronAbundance': UnitNo,
        'StarFormationRate': units.Msol / units.yr,
        'Metallicity': UnitNo,
        'StellarFormationTime': UnitNo,
        'BH_Mass': UnitMass,
        'BH_Mdot': UnitMassdTime,
        'BH_Mass_AlphaDisk': UnitMass,
        'BH_AccretionLength': UnitComvingLength,
        'BH_NProgs': UnitNo,
        'HaloID': UnitNo,
        'ID_Generations': UnitNo,
        'Potential': (UnitVel) ** 2 / units.a,
        'Dust_Masses': UnitMass,
        'Dust_Metallicity': UnitNo,
        
        
        'FractionH2': UnitNo,
        'GrackleHI': UnitNo,
        'GrackleHII': UnitNo,
        'GrackleHM': UnitNo,
        'GrackleHeI': UnitNo,
        'GrackleHeII': UnitNo,
        'GrackleHeIII': UnitNo,
        'NWindLaunches': UnitNo,
        'NeutralHydrogenAbundance': UnitNo,
        'Sigma': UnitNo,
        'AGS-Softening': UnitNo,
        
        
    
    }
    if field in Matchfieldunits:
        return Matchfieldunits[field]
    else:
        print(f"Parameter '{field}' not found in Matchfieldunits. use unitno")
        return UnitNo


def groupcat_units(
    field: str,
) -> units.Unit:
    """
    This function provides the unit corresponding to a given halo or galaxy parameter.

    Parameters:
    -----------
    field : str
        The name of the halo or galaxy parameter for which the unit is requested.

    Returns:
    --------
    unit
        The unit associated with the input parameter. If the parameter is not
        defined in `Matchfieldunits`, the function will raise a KeyError.
    """
    Matchfieldunits = {
        ### halo properties

    }
    if field in Matchfieldunits:
        return Matchfieldunits[field]
    else:
        raise KeyError(f"Parameter '{field}' not found in Matchfieldunits.")


