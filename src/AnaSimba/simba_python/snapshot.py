""" 
Simba Simulation.
catalog.py: File I/O related to the snapshots. 

author: Shuai Lu
"""


import h5py
import six

from .util import *
from .catalog import loadSingle,loadObjects



def snapPath(basePath,snapNum):
    snap_BasePath = search_files(f'{basePath}snapshots/*{snapNum:03d}.hdf5')
    if snap_BasePath:
        if len(snap_BasePath)>1:
            print(f'Multifiles found,{snap_BasePath}, use the {snap_BasePath[0]}')
        snap_BasePath = snap_BasePath[0]
        return snap_BasePath
    else:
        print('No snapshot file')
        return
    
def loadable_snapshotfields(basePath, snapNum, partType) -> list:
    partNum = partTypeNum(partType)
    partname = "PartType" + str(partNum)
    
    with h5py.File(snapPath(basePath,snapNum),'r') as f:
        fields = list(f[f"{partname}"].keys())
    return fields
    
def loadSubset(basePath, snapNum, partType, fields=['Coordinates','Masses','Velocities','ParticleIDs'], galaxyID = -1, haloID = -1) ->dict:
    result = {}
    
    if (haloID < 0 and galaxyID < 0) or (haloID >= 0 and galaxyID >= 0):
        raise Exception("Must specify either haloID or galaxyID (and not both).")
    grouptype = "galaxy_data" if galaxyID >= 0 else "halo_data"
    
    partNum = partTypeNum(partType)
    partname = "PartType" + str(partNum)
    realname = NumTopartname(partNum)
    pa_start = realname + 'list_start' if len(realname) < 3 else realname[:1]+'list_start'
    pa_end = realname + 'list_end' if len(realname) < 3 else realname[:1]+'list_end'
    
    catalog = loadSingle(basePath,snapNum,galaxyID = galaxyID,haloID=haloID)
    result['count'] = catalog[f'n{realname}']
    if result['count']==0:
        return result
    start = catalog[pa_start]
    end = catalog[pa_end]
    select = slice(start,end)
    if galaxyID >= 0:
        partID = loadObjects(basePath,snapNum,grouptype='global_lists',fields=f'galaxy_{pa_end[:-8]}list')
        select = (partID[f'galaxy_{pa_end[:-8]}list'] == galaxyID)
    
    with h5py.File(snapPath(basePath,snapNum),'r') as f:
        if not fields:
            fields = list(f[f'{partname}'].keys())
        if isinstance(fields, six.string_types):
            fields = [fields]
        for i in fields:
            result[i] = f[f'{partname}/{i}'][select]
    
    return result
    

def loadGalaxy(basePath, snapNum, id, partType, fields=['Coordinates','Masses','Velocities','ParticleIDs']) -> dict:
    
    return loadSubset(basePath,snapNum,partType,fields=fields,galaxyID=id)
    

def loadHalo(basePath, snapNum, id, partType, fields=['Coordinates','Masses','Velocities','ParticleIDs']) -> dict:
    return loadSubset(basePath,snapNum,partType,fields=fields,haloID=id)
