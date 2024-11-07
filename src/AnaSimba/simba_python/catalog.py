""" 
Simba Simulation.
catalog.py: File I/O related to the catalogs. 

author: Shuai Lu
"""

import h5py
import six

from .util import *

def catlogPath(basePath,snapNum):
    '''
    file path to catalog file
    '''
    catalog_BasePath = search_files(f'{basePath}catalogs/*{snapNum:03d}.hdf5')
    if catalog_BasePath:
        if len(catalog_BasePath)>1:
            print(f'Multifiles found,{catalog_BasePath}, use the {catalog_BasePath[0]}')
        catalog_BasePath = catalog_BasePath[0]
        return catalog_BasePath
    else:
        print('No catalog file')
        return
    
    
def loadGalaxies(basePath, snapNum, fields=None, dicts_fields=None,lists_fields = None):

    return loadObjects(basePath, snapNum, grouptype = 'galaxy_data', fields=fields,
                       dicts_fields=dicts_fields , lists_fields =lists_fields )
            

def loadHalos(basePath, snapNum, fields=None,dicts_fields=None,lists_fields = None):
    
    return loadObjects(basePath, snapNum, grouptype = 'halo_data', fields=fields,
                       dicts_fields=dicts_fields , lists_fields =lists_fields )

def loadable_catalogfields(basePath, snapNum, grouptype = 'galaxy_data'):
    result = {}
    with h5py.File(catlogPath(basePath,snapNum),'r') as f:
        result['fields'] = list(f[f'{grouptype}'].keys())
        result['lists_fields'] = list(f[f'{grouptype}/lists'].keys())
        result['dicts_fields'] = list(f[f'{grouptype}/dicts'].keys())
    return result

def loadObjects(basePath, snapNum, grouptype = 'galaxy_data', fields=None , dicts_fields=None,lists_fields = None, id=None):
    result = {}
    with h5py.File(catlogPath(basePath,snapNum),'r') as f:
        if not fields:
            fields = list(f[f'{grouptype}'].keys())
        if isinstance(fields, six.string_types):
            fields = [fields]
        if dicts_fields and 'dicts' not in fields:
            fields.append('dicts')
        if lists_fields and 'lists' not in fields:
            fields.append('lists')
            
        if isinstance(dicts_fields, six.string_types):
            dicts_fields = [dicts_fields]
        if isinstance(lists_fields, six.string_types):
            lists_fields = [lists_fields]
        for i in fields:
            if i in ['dicts','lists']:
                if i =='dicts' and dicts_fields:
                    result['dicts']={}
                    for j in dicts_fields:
                        if id:
                            result['dicts'][j] = f[f'{grouptype}/{i}/{j}'][tuple([id])]
                        else:
                            result['dicts'][j] = f[f'{grouptype}/{i}/{j}'][...]
                if i =='lists' and lists_fields:
                    result['lists']={}
                    for j in lists_fields:
                        if id:
                            result['lists'][j] = f[f'{grouptype}/{i}/{j}'][tuple([id])]
                        else:
                            result['lists'][j] = f[f'{grouptype}/{i}/{j}'][...]
                continue
            if id:
                result[i] = f[f'{grouptype}/{i}'][tuple([id])]
            else:
                result[i] = f[f'{grouptype}/{i}'][...]
    return result
    
    

def loadHeader(basePath, snapNum):
    header = {}
    with h5py.File(catlogPath(basePath,snapNum),'r') as f:
        simulation_attributes = {}
        for i,j in f['simulation_attributes'].attrs.items():
            simulation_attributes[i] = j
        parameters = {}
        for i,j in f['simulation_attributes/parameters'].attrs.items():
            parameters[i] = j
        units = {}
        for i,j in f['simulation_attributes/units'].attrs.items():
            units[i] = j
        header['simulation_attributes'] = simulation_attributes
        header['parameters'] = parameters
        header['units'] = units
    
    return header
    


def loadSingle(basePath :str, snapNum: int, haloID:int =-1, galaxyID:int =-1,dicts_fields=None,lists_fields = None):
    if (haloID < 0 and galaxyID < 0) or (haloID >= 0 and galaxyID >= 0):
        raise Exception("Must specify either haloID or galaxyID (and not both).")

    grouptype = "galaxy_data" if galaxyID >= 0 else "halo_data"
    searchID = galaxyID if galaxyID >= 0 else haloID
    return loadObjects(basePath=basePath,snapNum=snapNum,grouptype=grouptype, id = searchID,
                       dicts_fields = dicts_fields, lists_fields=lists_fields )