import numpy as np
try:
    from . import tools as tbx
except ValueError:
    import tools as tbx
import zipfile

datdir = 'bathy/'
fnm = 'g1230485'


def pull_data():
    tbx.retrieve(
        'https://www.ocean.washington.edu/data/pugetsound/datasets/psdem2005/rasters/tiles/g1230485/g1230485.zip',
        datdir + fnm + '.zip'
    )


def process_data():
    with zipfile.ZipFile(datdir + fnm + '.zip') as zf:
        with zf.open(fnm + '/' + fnm + '.asc') as fl:

            # Get the header information.
            params = {}
            for ln in fl:
                vals = ln.split()
                params[vals[0]] = int(vals[1])
                if vals[0] == 'NODATA_value':
                    break

            # Read the data
            dat = np.loadtxt(fl, dtype=np.float32)

    np.savez_compressed(datdir + fnm + '.npz', elev=dat, **params)
