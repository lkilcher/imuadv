from .. import tools as tbx
from scipy.io import loadmat
import dolfyn.adv.api as avm
from os.path import isfile

datdir = tbx.datdir + 'smb_may2015/'


def load(tag, bindat=True):
    """Load ADP-BT motion-corrected StableMoor data.

    Possible values of `tag` are:
    'SMN-BT'      : Load SM nose ADP BT file.

    'SMN-5s'      : 5s filtered SM nose
    'SMN-10s'     : 10s filtered SM nose
    'SMN-10s'     : 10s filtered SM nose
    'SMN-30s'     : 30s filtered SM nose
    'SMN-unfilt'  : unfiltered SM nose data

    """
    for fl in [tag,
               tag + '.h5',
               datdir + tag,
               datdir + tag + '.h5']:
        if isfile(fl):
            return avm.load(fl)
    if tag == 'SMN-BT':
        # This file is from Sam Harding.
        dat = loadmat(datdir + 'BT_IMU_sync_data.mat')
        dat['t_IMU'] = dat['t_IMU'][0] - 366
        return dat
    dset, filt_tag = tag.split('-')
    if bindat:
        prefix = 'bindat'
    else:
        prefix = 'mcdat'
    return avm.load(datdir + '{}_{}_filt-{}.h5'.format(dset, prefix, filt_tag))
