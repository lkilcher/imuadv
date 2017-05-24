import numpy as np
try:
    from . import tools as tbx
except ValueError:
    import tools as tbx
import zipfile

archive_name = 'g1230485'
fnm = 'bathy-data'


def pull():
    print("Retrieving bathy data...")
    tbx.retrieve(
        'https://www.ocean.washington.edu/data/pugetsound/datasets/psdem2005/rasters/tiles/g1230485/g1230485.zip',
        tbx.datdir + fnm + '.zip',
        hash='6dcc6bddac16284a'
    )


def process():
    print("Processing bathy data...")
    outfile = tbx.datdir + fnm + '.npz'
    if tbx.checkhash(outfile, '37400bc9b7829187'):
        print("   hash check passed; skipping processing.")
        return
    with zipfile.ZipFile(tbx.datdir + fnm + '.zip') as zf:
        with zf.open(archive_name + '/' + archive_name + '.asc') as fl:

            # Get the header information.
            params = {}
            for ln in fl:
                vals = ln.split()
                params[vals[0]] = int(vals[1])
                if vals[0] == 'NODATA_value':
                    break

            # Read the data
            dat = np.loadtxt(fl, dtype=np.float32)

    np.savez_compressed(outfile, elev=dat, **params)
