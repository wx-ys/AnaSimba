'''
Some useful tools

potential: cal_potential, cal_acceleration, Function.
galaxy profile: single: profile(), all: Profile_1D(). Class.
...
'''

from typing import List
import multiprocessing as mp
import re
import math

import numpy as np
import h5py
from tqdm import tqdm
from pynbody import units, filt
from pynbody.array import SimArray
from pynbody.analysis.profile import Profile 

from .Anatools import Orbit
from .pytreegrav import PotentialTarget, AccelTarget
from .simba_snapshot import Basehalo


def cal_potential(sim, targetpos):
    """
    Calculates the gravitational potential at target positions.

    Parameters:
    -----------
    sim : object
        The simulation data object containing particle positions and masses.
    targetpos : array-like
        The positions where the gravitational potential needs to be calculated.

    Returns:
    --------
    phi : SimArray
        The gravitational potential at the target positions.
    """

    try:
        eps = sim.properties.get('eps', 0)
    except:
        eps = 0
    if eps == 0:
        print('Calculate the gravity without softening length')
    pot = PotentialTarget(
        targetpos,
        sim['pos'].view(np.ndarray),
        sim['mass'].view(np.ndarray),
        np.repeat(eps, len(targetpos)).view(np.ndarray),
    )
    phi = SimArray(pot, units.G * sim['mass'].units / sim['pos'].units)
    phi.sim = sim
    return phi


def cal_acceleration(sim, targetpos):
    """
    Calculates the gravitational acceleration at specified target positions.

    Parameters:
    -----------
    sim : object
        The simulation data object containing particle positions and masses.
    targetpos : array-like
        The positions where the gravitational acceleration needs to be calculated.

    Returns:
    --------
    acc : SimArray
        The gravitational acceleration at the target positions.
    """
    try:
        eps = sim.properties.get('eps', 0)
    except:
        eps = 0
    if eps == 0:
        print('Calculate the gravity without softening length')
    accelr = AccelTarget(
        targetpos,
        sim['pos'].view(np.ndarray),
        sim['mass'].view(np.ndarray),
        np.repeat(eps, len(targetpos)).view(np.ndarray),
    )
    acc = SimArray(
        accelr, units.G * sim['mass'].units / sim['pos'].units / sim['pos'].units
    )
    acc.sim = sim
    return acc

class profile(Profile):
    
    def _calculate_x(self, sim):
        if self.zmax:
            return SimArray(np.abs(sim['z']), sim['z'].units)
        else:
            return ((sim['pos'][:, 0:self.ndim] ** 2).sum(axis=1)) ** (1, 2)
    
    def __init__(self, sim, rmin = 0.1, rmax = 30, 
                 nbins=100, ndim=2, type='lin', weight_by='mass', calc_x=None, **kwargs):
        
        zmax = kwargs.get('zmax', None)
        self.zmax = zmax
        if isinstance(zmax, str):
            zmax = units.Unit(zmax)
        
        if self.zmax:
            if isinstance(rmin, str):
                rmin = units.Unit(rmin)
            if isinstance(rmax, str):
                rmax = units.Unit(rmax)
            self.rmin = rmin
            self.rmax = rmax

            assert ndim in [2, 3]
            if ndim == 3:
                sub_sim = sim[
                    filt.Disc(rmax, zmax) & ~filt.Disc(rmin, zmax)]
            else:
                sub_sim = sim[(filt.BandPass('x', rmin, rmax) |
                            filt.BandPass('x', -rmax, -rmin)) &
                            filt.BandPass('z', -zmax, zmax)]

            Profile.__init__(
                self, sub_sim, nbins=nbins, weight_by=weight_by, 
                ndim=ndim, type=type, **kwargs)
        else:
            Profile.__init__(
                self, sim, rmin=rmin, rmax=rmax, nbins=nbins, weight_by=weight_by, 
                ndim=ndim, type=type, **kwargs)
    
        
    def _setup_bins(self):
        Profile._setup_bins(self)
        if self.zmax:
            dr = self.rmax - self.rmin

            if self.ndim == 2:
                self._binsize = (
                    self['bin_edges'][1:] - self['bin_edges'][:-1]) * dr
            else:
                area = SimArray(
                    np.pi * (self.rmax ** 2 - self.rmin ** 2), "kpc^2")
                self._binsize = (
                    self['bin_edges'][1:] - self['bin_edges'][:-1]) * area
    def _get_profile(self, name):
        """Return the profile of a given kind"""
        x = name.split(",")
        find = re.search(r'_\d+',name)
        if name in self._profiles:
            return self._profiles[name]

        elif x[0] in Profile._profile_registry:
            args = x[1:]
            self._profiles[name] = Profile._profile_registry[x[0]](self, *args)
            try:
                self._profiles[name].sim = self.sim
            except AttributeError:
                pass
            return self._profiles[name]

        elif name in list(self.sim.keys()) or name in self.sim.all_keys():
            self._profiles[name] = self._auto_profile(name)
            self._profiles[name].sim = self.sim
            return self._profiles[name]

        elif name[-5:] == "_disp" and (name[:-5] in list(self.sim.keys()) or name[:-5] in self.sim.all_keys()):
            self._profiles[name] = self._auto_profile(
                name[:-5], dispersion=True)
            self._profiles[name].sim = self.sim
            return self._profiles[name]

        elif name[-4:] == "_rms" and (name[:-4] in list(self.sim.keys()) or name[:-4] in self.sim.all_keys()):
            self._profiles[name] = self._auto_profile(name[:-4], rms=True)
            self._profiles[name].sim = self.sim
            return self._profiles[name]

        elif name[-4:] == "_med" and (name[:-4] in list(self.sim.keys()) or name[:-4] in self.sim.all_keys()):
            self._profiles[name] = self._auto_profile(name[:-4], median=True)
            self._profiles[name].sim = self.sim
            return self._profiles[name]
        
        elif name[-4:] == "_sum" and (name[:-4] in list(self.sim.keys()) or name[:-4] in self.sim.all_keys()):
            self._profiles[name] = self._auto_profile(name[:-4], sum=True)
            self._profiles[name].sim = self.sim
            return self._profiles[name]

        elif name[0:2] == "d_" and (name[2:] in list(self.keys()) or name[2:] in self.derivable_keys() or name[2:] in self.sim.all_keys()):
            #            if np.diff(self['dr']).all() < 1e-13 :
            self._profiles[name] = np.gradient(self[name[2:]], self['dr'][0])
            self._profiles[name] = self._profiles[name] / self['dr'].units
            return self._profiles[name]
            # else :
            #    raise RuntimeError, "Derivatives only possible for profiles of fixed bin width."
        elif find and (name[:find.start()] in list(self.sim.keys()) or name[:find.start()] in self.sim.all_keys()):
            self._profiles[name] = self._auto_profile(name[:find.start()], q = float(name[find.start()+1:]))
            self._profiles[name].sim = self.sim
            return self._profiles[name]
            
        else:
            raise KeyError(name + " is not a valid profile")

    def _auto_profile(self, name, dispersion=False, rms=False, median=False,sum=False, q=None ):
        result = np.zeros(self.nbins)

        # force derivation of array if necessary:
        self.sim[name]

        for i in range(self.nbins):
            subs = self.sim[self.binind[i]]
            name_array = subs[name].view(np.ndarray)
            mass_array = subs[self._weight_by].view(np.ndarray)

            if dispersion:
                sq_mean = (name_array ** 2 * mass_array).sum() / \
                    self['weight_fn'][i]
                mean_sq = (
                    (name_array * mass_array).sum() / self['weight_fn'][i]) ** 2
                try:
                    result[i] = math.sqrt(sq_mean - mean_sq)
                except ValueError:
                    # sq_mean<mean_sq occasionally from numerical roundoff
                    result[i] = 0

            elif rms:
                result[i] = np.sqrt(
                    (name_array ** 2 * mass_array).sum() / self['weight_fn'][i])
            elif sum:
                result[i] = name_array.sum()
            elif median:
                if len(subs) == 0:
                    result[i] = np.nan
                else:
                    sorted_name = sorted(name_array)
                    result[i] = sorted_name[int(np.floor(0.5 * len(subs)))]
            elif q:
                if len(subs) == 0:
                    result[i] = np.nan
                else:
                    sorted_name = sorted(name_array)
                    weight_array = mass_array[np.argsort(name_array)]
                    cumw = np.cumsum(weight_array) / np.sum(weight_array)
                    imin = min(
                            np.arange(len(sorted_name)), key=lambda x: abs(cumw[x] - q/100))
                    inc = q/100 - cumw[imin]
                    lowval = sorted_name[imin]
                    if inc > 0:
                        nextval = sorted_name[imin + 1]
                    else:
                        if imin == 0:
                            nextval = lowval
                        else:
                            nextval = sorted_name[imin - 1]

                    result[i] = lowval + inc * (nextval - lowval)
                    #result[i] = sorted_name[cumw*100>q].min()+sorted_name[cumw*100<q].max()
            else:
                result[i] = (name_array * mass_array).sum() / self['weight_fn'][i]

        result = result.view(SimArray)
        result.units = self.sim[name].units
        result.sim = self.sim
        return result
    
@Profile.profile_property
def v_circ(p, grav_sim=None):
    """Circular velocity, i.e. rotation curve. Calculated by computing the gravity
    in the midplane - can be expensive"""
    # print("Profile v_circ -- this routine assumes the disk is in the x-y plane")
    grav_sim = grav_sim or p.sim
    cal_2 = np.sqrt(2) / 2
    basearray = np.array(
        [
            (1, 0, 0),
            (0, 1, 0),
            (-1, 0, 0),
            (0, -1, 0),
            (cal_2, cal_2, 0),
            (-cal_2, cal_2, 0),
            (cal_2, -cal_2, 0),
            (-cal_2, -cal_2, 0),
        ]
    )
    R = p['rbins'].in_units('kpc').copy()
    POS = np.array([(0, 0, 0)])
    for j in R:
        binsr = basearray * j
        POS = np.concatenate((POS, binsr), axis=0)
    POS = SimArray(POS, R.units)
    ac = cal_acceleration(grav_sim, POS)
    ac.convert_units('kpc Gyr**-2')
    POS.convert_units('kpc')
    velall = np.diag(np.dot(ac - ac[0], -POS.T))
    if 'units' in dir(velall):
        velall.units = units.kpc**2 / units.Gyr**2
    else:
        velall = SimArray(velall, units.kpc**2 / units.Gyr**2)
    velTrue = np.zeros(len(R))
    for i in range(len(R)):
        velTrue[i] = np.mean(velall[i + 1 : 8 * (i + 1) + 1])
    velTrue[velTrue < 0] = 0
    velTrue = np.sqrt(velTrue)
    velTrue = SimArray(velTrue, units.kpc / units.Gyr)
    velTrue.convert_units('km s**-1')
    velTrue.sim = grav_sim.ancestor
    return velTrue
@Profile.profile_property
def pot(p, grav_sim=None):
    grav_sim = grav_sim or p.sim
    cal_2 = np.sqrt(2) / 2
    basearray = np.array(
        [
            (1, 0, 0),
            (0, 1, 0),
            (-1, 0, 0),
            (0, -1, 0),
            (cal_2, cal_2, 0),
            (-cal_2, cal_2, 0),
            (cal_2, -cal_2, 0),
            (-cal_2, -cal_2, 0),
        ]
    )
    R = p['rbins'].in_units('kpc').copy()
    POS = np.array([(0, 0, 0)])
    for j in R:
        binsr = basearray * j
        POS = np.concatenate((POS, binsr), axis=0)
    POS = SimArray(POS, R.units)
    po = cal_potential(grav_sim, POS)
    po.convert_units('km**2 s**-2')
    poall = np.zeros(len(R))
    for i in range(len(R)):
        poall[i] = np.mean(po[i + 1 : 8 * (i + 1) + 1])

    poall = SimArray(poall, po.units)
    poall.sim = grav_sim.ancestor
    return poall

@Profile.profile_property
def omega(p):
    """Circular frequency Omega = v_circ/radius (see Binney & Tremaine Sect. 3.2)"""
    prof = p['v_circ'] / p['rbins']
    prof.convert_units('km s**-1 kpc**-1')
    return prof

class Profile_1D:
    _properties={}
    def __init__(
        self, sim, rmin=0.1, rmax=100.0, zmax = 5.,nbins=100, type='lin', **kwargs
    ):
        """
        Initializes the profile object for different types of particles in the simulation.

        Parameters:
        -----------
        sim : object
            The simulation data object containing particles of different types (e.g., stars, gas, dark matter).
        rmin : float, optional
            The minimum radius for the profile (default is 0.1).
        rmax : float, optional
            The maximum radius for the profile (default is 100.0).
        zmax : float, optional
            maximum height to consider (default is 5.0).
        nbins : int, optional
            The number of bins to use in the profile (default is 100).
        type : str, optional
            The type of profile ('lin' for linear or other types as needed, default is 'lin').

        **kwargs : additional keyword arguments
            Additional parameters to pass to the Profile initialization.
            
        Usage: str like 'A-B-C'
                A: the parameter key,  d_A, derivatives, A_disp, A_med, A_rms, A_30 ...
                B: family, 'star', 'gas', 'dm', 'all'
                C: direction and dims, 'z', 'Z', 'r', 'R'; 'z' vertical and 3 dims, 'Z' 2dims ... 
            examples : 'vr-star-R'
        """
        print(
            "Profile_1D -- assumes it's already at the center, and the disk is in the x-y plane"
        )
        print("If not, please use face_on()")
        self.__P={'all':{}, 'star':{}, 'gas':{}, 'dm':{}}
        self.__P['all']['r']=profile(sim, rmin=rmin, rmax=rmax, nbins=nbins,ndim=3, type=type, **kwargs)
        self.__P['all']['R']=profile(sim, rmin=rmin, rmax=rmax, nbins=nbins, ndim=2, type=type, **kwargs)
        self.__P['all']['Z']=profile(sim, rmin=rmin, rmax=rmax, nbins=nbins, ndim=2, type=type, zmax = zmax, **kwargs)
        self.__P['all']['z']=profile(sim, rmin=rmin, rmax=rmax, nbins=nbins, ndim=3, type=type, zmax = zmax,**kwargs)
        
        self.__P['star']['r']=profile(sim.s, rmin=rmin, rmax=rmax, nbins=nbins,ndim=3, type=type, **kwargs)
        self.__P['star']['R']=profile(sim.s, rmin=rmin, rmax=rmax, nbins=nbins, ndim=2, type=type, **kwargs)
        self.__P['star']['Z']=profile(sim.s, rmin=rmin, rmax=rmax, nbins=nbins, ndim=2, type=type, zmax = zmax, **kwargs)
        self.__P['star']['z']=profile(sim, rmin=rmin, rmax=rmax, nbins=nbins, ndim=3, type=type, zmax = zmax,**kwargs)
        try:
            self.__P['gas']['r']=profile(sim.g, rmin=rmin, rmax=rmax, nbins=nbins,ndim=3, type=type, **kwargs)
        except:
            print('No gas r')
        try:
            self.__P['gas']['R']=profile(sim.g, rmin=rmin, rmax=rmax, nbins=nbins, ndim=2, type=type, **kwargs)
        except:
            print('No gas R')
        try:
            self.__P['gas']['Z']=profile(sim.g, rmin=rmin, rmax=rmax, nbins=nbins, ndim=2, type=type, zmax = zmax, **kwargs)
        except:
            print('No gas Z')
        try:
            self.__P['gas']['z']=profile(sim.g, rmin=rmin, rmax=rmax, nbins=nbins, ndim=3, type=type, zmax = zmax,**kwargs)
        except:
            print('No gas z')
        
        self.__P['dm']['r']=profile(sim.dm, rmin=rmin, rmax=rmax, nbins=nbins, ndim=3, type=type, **kwargs)
        self.__P['dm']['R']=profile(sim.dm, rmin=rmin, rmax=rmax, nbins=nbins, ndim=2, type=type, **kwargs)
        self.__P['dm']['Z']=profile(sim.dm, rmin=rmin, rmax=rmax, nbins=nbins, ndim=2, type=type, zmax = zmax, **kwargs)
        self.__P['dm']['z']=profile(sim.dm, rmin=rmin, rmax=rmax, nbins=nbins, ndim=3, type=type, zmax = zmax,**kwargs)

    
    def _util_fa(self, ks):
        if set(['star', 's', 'Star']) & set(ks):
            return 'star'
        if set(['gas', 'g', 'Gas']) & set(ks):
            return 'gas'
        if set(['dm', 'DM']) & set(ks):
            return 'dm'
        if set(['all', 'ALL']) & set(ks):
            return 'all'
        return 'all'
    
    def _util_pr(self, ks):
        if set(['r']) & set(ks):
            return 'r'
        if set(['z']) & set(ks):
            return 'z'
        if set(['R']) & set(ks):
            return 'R'
        if set(['Z']) & set(ks):
            return 'Z'
        return 'R'   

    def __getitem__(self, key):

        if isinstance(key, str):
            ks = key.split('-')
            if len(ks) > 1:
                return self.__P[self._util_fa(ks)][self._util_pr(ks)][ks[0]]
            else:
                if key in self._properties:
                    return self._properties[key](self)
                else:
                    return self.__P['all']['R'][key]
        else:
            print('Type error, should input a str')
            return
    @staticmethod
    def profile_property(fn):
        Profile_1D._properties[fn.__name__] = fn
        return fn
    
@Profile_1D.profile_property    
def Qgas(self):
    '''
    Toomre-Q for gas
    '''
    return (
        self['kappa-all-R']
        * self['vrxy_disp-gas-R']
        / (np.pi * self['density-gas-R'] * units.G)
    ).in_units("")
    
@Profile_1D.profile_property  
def Qstar(self):
    '''
    Toomre-Q parameter
    '''
    return (
        self['kappa-all-R']
        * self['vrxy_disp-star-R']
        / (3.36 * self['density-star-R'] * units.G)
    ).in_units("")
    
@Profile_1D.profile_property  
def Qs(self):
    '''
    Toomre-Q parameter
    '''
    return (
        self['kappa-all-R']
        * self['vrxy_disp-star-R']
        / (np.pi * self['density-star-R'] * units.G)
    ).in_units("")
    
@Profile_1D.profile_property  
def Q2ws(self):
    '''
    Toomre Q of two component. Wang & Silk (1994)
    '''
    Qs = self['Qs']
    Qg = self['Qgas']
    return (Qs * Qg) / (Qs + Qg)

@Profile_1D.profile_property  
def Q2thin(self):
    '''
    The effective Q of two component thin disk. Romeo & Wiegert (2011) eq. 6.
    '''
    w = (
        2
        * self['vrxy_disp-star-R']
        * self['vrxy_disp-gas-R']
        / ((self['vrxy_disp-star-R']) ** 2 + self['vrxy_disp-gas-R'] ** 2)
    ).in_units("")
    Qs = self['Qs']
    Qg = self['Qgas']
    q = [Qs * Qg / (Qs + w * Qg)]
    return [
        (
            Qs[i] * Qg[i] / (Qs[i] + w[i] * Qg[i])
            if Qs[i] > Qg[i]
            else Qs[i] * Qg[i] / (w[i] * Qs[i] + Qg[i])
        )
        for i in range(len(w))
    ]
    
@Profile_1D.profile_property  
def Q2thick(self):
    '''
    The effective Q of two component thick disk. Romeo & Wiegert (2011) eq. 9.
    '''
    w = (
        2
        * self['vrxy_disp-star-R']
        * self['vrxy_disp-gas-R']
        / ((self['vrxy_disp-star-R']) ** 2 + self['vrxy_disp-gas-R'] ** 2)
    ).in_units("")
    Ts = 0.8 + 0.7 * (self['vz_disp-star-R'] / self['vrxy_disp-star-R']).in_units(
        ""
    )
    Tg = 0.8 + 0.7 * (self['vz_disp-gas-R'] / self['vrxy_disp-gas-R']).in_units("")
    Qs = self['Qs']
    Qg = self['Qgas']
    Qs = Qs * Ts
    Qg = Qg * Tg
    return [
        (
            Qs[i] * Qg[i] / (Qs[i] + w[i] * Qg[i])
            if Qs[i] > Qg[i]
            else Qs[i] * Qg[i] / (w[i] * Qs[i] + Qg[i])
        )
        for i in range(len(w))
    ]









