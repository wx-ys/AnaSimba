## Introduction
AnaSimba is a python package for processing and analyzing the cosmological simulation [Simba](http://simba.roe.ac.uk/).

## Installation

```
git clone https://github.com/wx-ys/AnaSimba.git
```
Install this in editable mode.
```
cd AnaSimba
pip install -e .
```
AnaSimba uses the following python packages:

* numpy, scipy
* pynbody
* h5py
* tqdm
* six
* numba

## Usage


```python
from AnaSimba import simba_simulation 
BasePath = 'filepath'       
snap=151  #snapshot

Snapshot=simba_simulation.Snapshot(BasePath,snap) # use the data of snapshot151

# load a single galaxy
sub = snapshot.load_particle(galaxyID=30,decorate=True)
sub.physical_units() #in physical unit
sub.face_on(alignwith='star',rmax=8) # face on, based on stellar angular momentum.

```

See [example](example) for more details,



## License

[MIT](LICENSE) © Shuai Lu

## Acknowledgments
* [pytreegrav](https://github.com/mikegrudic/pytreegrav)
* [pynbody](https://github.com/pynbody/pynbody)