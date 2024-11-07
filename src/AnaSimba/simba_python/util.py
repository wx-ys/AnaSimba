import glob


def search_files(pattern):
    files = glob.glob(pattern)
    return files


def partTypeNum(partType) -> int:
    """ Mapping between common names and numeric particle types. """
    if str(partType).isdigit():
        return int(partType)
        
    if str(partType).lower() in ['gas','cells','g',0]:
        return 0
    if str(partType).lower() in ['dm','darkmatter',1]:
        return 1
    if str(partType).lower() in ['dmlowres',2]:
        return 2 # only zoom simulations, not present in full periodic boxes
    if str(partType).lower() in ['tracer','tracers','tracermc','trmc',3]:
        return 3
    if str(partType).lower() in ['star','stars','stellar','s',4]:
        return 4 # only those with GFM_StellarFormationTime>0
    if str(partType).lower() in ['wind',4]:
        return 4 # only those with GFM_StellarFormationTime<0
    if str(partType).lower() in ['bh','bhs','blackhole','blackholes',5]:
        return 5
    
    raise Exception("Unknown particle type name.")


def NumTopartname(partNum: int) -> str:
    if partNum==0:
        return 'gas'
    if partNum==1:
        return 'dm'
    if partNum==4:
        return 'star'
    if partNum==5:
        return 'bh'