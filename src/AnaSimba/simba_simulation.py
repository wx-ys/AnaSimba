'''
Load Simba data and process it.
'''

from functools import reduce
from copy import deepcopy


from pynbody.snapshot import SimSnap, new
from pynbody.family import get_family
from pynbody.array import SimArray
from pynbody import filt,units
from pynbody.simdict import SimDict
import h5py

import AnaSimba.simba_python as simba
from AnaSimba.simba_units import snapshot_pa_name,snapshot_units
from AnaSimba.simba_snapshot import Basehalo,read_Snap_properties



class Snapshot(SimSnap):

    def __init__(
        self,
        BasePath: str,
        Snap: int,
    ):
        SimSnap.__init__(self)
        self._num_particles = 0
        self._filename = "<created>"
        self._create_arrays(["pos", "vel"], 3)
        self._create_arrays(["mass"], 1)
        self._family_slice[get_family('dm')] = slice(0, 0)
        self._family_slice[get_family('star')] = slice(0, 0)
        self._family_slice[get_family('gas')] = slice(0, 0)
        self._family_slice[get_family('bh')] = slice(0, 0)
        self._decorate()
        self.__set_Snapshot_property(BasePath, Snap)
        self.properties['filedir'] = BasePath
        self.properties['filepath'] = BasePath
        self._filename = self.properties['run']
        
        self.properties['baseunits'] = [
            units.Unit(x) for x in ('kpc', 'km s^-1', 'Msol')
        ]
        self.properties['staunit'] = [

        ]
        for i in self.properties:
            if isinstance(self.properties[i], SimArray):
                self.properties[i].sim = self
        self._filename = self._filename + '_' + 'snapshot' + str(Snap)

        self.__set_load_particle()
        self._canloadPT = True
        self.__PT_loaded = {'halo': set(), 'galaxy': set()}

        self.__GC_loaded = {'halo': set(), 'galaxy': set()}

        self.__pos = SimArray([0.0, 0.0, 0.0], units.kpc)
        self.__pos.sim = self
        self.__vel = SimArray([0.0, 0.0, 0.0], units.km / units.s)
        self.__vel.sim = self
        self.__phi = SimArray([0.0, 0.0, 0.0], units.km**2 / units.s**2)
        self.__phi.sim = self
        self.__acc = SimArray([0.0, 0.0, 0.0], units.km / units.s**2)
        self.__acc.sim = self
        
        __star_pa = simba.loadable_snapshotfields(BasePath,Snap,'star')
        __gas_pa = simba.loadable_snapshotfields(BasePath,Snap,'gas')
        __bh_pa = simba.loadable_snapshotfields(BasePath,Snap,'bh')
        __dm_pa = simba.loadable_snapshotfields(BasePath,Snap,'dm')
        
        __halo_pa = simba.loadable_catalogfields(BasePath,Snap,'halo_data')
        __galaxy_pa = simba.loadable_catalogfields(BasePath,Snap,'galaxy_data')



    def load_particle(
        self, galaxyID = -1, haloID = -1, decorate=True,**kwargs
    ) -> SimSnap:
        catalog = simba.loadSingle(self.properties['filedir'],self.snapshot,galaxyID = galaxyID,haloID=haloID)
        groupType = 'Galaxy' if galaxyID>=0 else 'Halo'
        
        order = kwargs.get('order', self.load_particle_para['particle_field'])
        f = new(
            dm=int(catalog['ndm']),
            star=int(catalog['nstar']),
            gas=int(catalog['ngas']),
            bh=int(catalog['nbh']),
            order=order,
        )

        for party in self.load_particle_para['particle_field'].split(","):
            if len(f[get_family(party)]) > 0:
                if len(self.load_particle_para[party + '_fields']) > 0:
                    self.load_particle_para[party + '_fields'] = list(
                        set(
                            self.load_particle_para[party + '_fields']
                            + self.load_particle_para['Basefields']
                        )
                    )
                else:
                    self.load_particle_para[party + '_fields'] = list.copy(
                        self.load_particle_para['Basefields']
                    )

                if party == 'dm':
                    loaddata = simba.loadSubset(
                        self.properties['filedir'],
                        self.snapshot,
                        party,
                        self.load_particle_para[party + '_fields'],
                        galaxyID,
                        haloID
                    )
                    for i in self.load_particle_para[party + '_fields']:
                        f.dm[snapshot_pa_name(i)] = SimArray(
                            loaddata[i], snapshot_units(i)
                        )

                if party == 'star':
                    loaddata = simba.loadSubset(
                        self.properties['filedir'],
                        self.snapshot,
                        party,
                        self.load_particle_para[party + '_fields'],
                        galaxyID,
                        haloID
                    )
                    for i in self.load_particle_para[party + '_fields']:
                        f.s[snapshot_pa_name(i)] = SimArray(
                            loaddata[i], snapshot_units(i)
                        )

                if party == 'gas':
                    loaddata = simba.loadSubset(
                        self.properties['filedir'],
                        self.snapshot,
                        party,
                        self.load_particle_para[party + '_fields'],
                        galaxyID,
                        haloID
                    )
                    for i in self.load_particle_para[party + '_fields']:
                        f.g[snapshot_pa_name(i)] = SimArray(
                            loaddata[i], snapshot_units(i)
                        )

                if party == 'bh':
                    loaddata = simba.loadSubset(
                        self.properties['filedir'],
                        self.snapshot,
                        party,
                        self.load_particle_para[party + '_fields'],
                        galaxyID,
                        haloID
                    )
                    for i in self.load_particle_para[party + '_fields']:
                        f.bh[snapshot_pa_name(i)] = SimArray(
                            loaddata[i], snapshot_units(i)
                        )
        f.properties = deepcopy(self.properties)
        for i in f.properties:
            if isinstance(f.properties[i], SimArray):
                f.properties[i].sim = f
        f._filename = self.filename + '_' + groupType + '_' + str(max(galaxyID,haloID))
        if decorate:
            return Basehalo(f)
        
        return f
    
    def __repr__(self):
        return "<Snapshot \"" + self.filename + "\" len=" + str(len(self)) + ">"

    def __getattr__(self, name):
        try:
            return super().__getattr__(name)
        except:
            pass

        try:
            return self.properties[name]
        except:
            pass

        raise AttributeError(
            "%r object has no attribute %r" % (type(self).__name__, name)
        )
    
    
    def __set_load_particle(self):
        pa = {}
        pa['particle_field'] = 'dm,star,gas,bh'
        pa['Basefields'] = ['Coordinates', 'Velocities', 'Masses', 'ParticleIDs']
        pa['star_fields'] = []
        pa['gas_fields'] = []
        pa['dm_fields'] = []
        pa['bh_fields'] = []
        self.load_particle_para = pa
        
    def __set_Snapshot_property(self, BasePath: str, Snap: int):
        """
        Init properties (Simdict) from Base Path and Snap.
        """

        SnapshotHeader = simba.loadHeader(BasePath, Snap)
        self.properties = SimDict()
        for i in SnapshotHeader['parameters']:
            self.properties[i] = SnapshotHeader['parameters'][i]
        for i in SnapshotHeader['simulation_attributes']:
            self.properties[i] = SnapshotHeader['simulation_attributes'][i]
        self.properties['a'] = SnapshotHeader['simulation_attributes']['scale_factor']
        self.properties['z'] = (1 / self.properties['a']) - 1       # Redshift
        self.properties['omegaB0'] = SnapshotHeader['simulation_attributes']['omega_baryon'] 
        self.properties['read_Snap_properties'] = SnapshotHeader['parameters']
        for i in self.properties:
            if 'sim' in dir(self.properties[i]):
                self.properties[i].sim = self
        self.properties['filedir'] = BasePath
        self.properties['Snapshot'] = Snap
        self.properties['run'] = SnapshotHeader['simulation_attributes']['basename'].split('_')[1]
        
    @property
    def snapshot(self):
        return self.properties['Snapshot']
        