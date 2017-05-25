import numpy as np
import dolfyn.adv.api as avm
import scipy.signal as sig
#import ttm.sm2015 as data_api
import os
from .base import datdir, load

filt_freqs = {
    'unfilt': 0.0,
    '5s': 1. / 5,
    '10s': 1. / 10,
    '30s': 1. / 30
}

doppler_noise = [2e-5, 2e-5, 0]

eps_freqs = np.array([[.1, 3],
                      [.1, 3],
                      [.1, 3], ])


if __name__ == '__main__':

    try:
        os.mkdir(datdir)
    except:
        pass

    # 4800 points is 5min at 16hz
    binner = avm.TurbBinner(4800, 16)

    bt = load('SMN-BT')

    bt_var = 'ubt'
    bt_var = 'ubt2'

    # Load the 'raw' (unprocessed) data that corresponds to the file above.
    dat = data_api.load('SM_Nose', coordsys='raw')

    dat.props['noise'] = doppler_noise

    # The offset-index between S. Harding's bt file, and the ADV data (`dat`)
    # S. Harding already corrected offsets, so t_IMU is his estimate
    # of ADV time.
    offind = np.abs(bt['t_IMU'][0] - dat.mpltime).argmin()

    # Perform basic motion correction based on IMU data:
    for filt_tag, filt_freq in filt_freqs.iteritems():
        datmc = dat.copy()
        avm.motion.correct_motion(datmc, accel_filtfreq=filt_freq)

        # Crop the data so that it is consistent with the ADP data.
        datmc = datmc.subset(slice(offind, (offind + len(bt['t_IMU']))))

        # Calculate 'rotational' motion of ADP head.
        motcalc = avm.motion.CalcMotion(datmc)
        ur_adp = motcalc.calc_urot(np.array([-1, 0, 0]), to_earth=False)

        datmc.add_data('ubt', bt['BT_INST_interp'], 'orient')
        datmc.add_data('ubt2', bt['BT_INST_interp'] + ur_adp, 'orient')
        datmc.add_data('ur_adp', ur_adp, 'orient')
        bt_vars = {'ubt', 'ubt2'}
        datmc.rotate_vars.update(bt_vars)
        # The ADP data is in the instrument frame, so rotate it into the earth frame.
        avm.rotate.inst2earth(datmc, rotate_vars=bt_vars, force=True)
        if filt_freq > 0:
            filt = sig.butter(2, filt_freq / (datmc.fs / 2))
            datmc.ubt = sig.filtfilt(filt[0], filt[1], datmc.ubt)
            datmc.ubt2 = sig.filtfilt(filt[0], filt[1], datmc.ubt2)
        if filt_freq > 0:
            datmc._u += datmc[bt_var]  # Add the bt to the vel
        datnow = datmc.copy()
        avm.rotate.inst2earth(datmc, reverse=True)

        datmc.save(datdir + 'SMN_mcdat_filt-{}.h5'.format(filt_tag))

        avm.rotate.earth2principal(datnow)
        datbd = binner(datnow)
        datbd.add_data('Spec_uraw', binner.psd(datnow.uraw), 'spec')
        datbd.add_data('Spec_uacc', binner.psd(datnow.uacc), 'spec')
        datbd.add_data('Spec_urot', binner.psd(datnow.urot), 'spec')
        datbd.add_data('Spec_ubt', binner.psd(datnow.ubt), 'spec')
        datbd.add_data('Spec_ubt2', binner.psd(datnow.ubt2), 'spec')
        if filt_freq > 0:
            umot = datnow.uacc + datnow.urot + datnow[bt_var]
        else:
            umot = datnow.uacc + datnow.urot
        datbd.add_data('Spec_umot', binner.psd(umot), 'spec')

        datbd.add_data('Cspec_umot', binner.calc_vel_cpsd(umot), 'spec')
        datbd.add_data('Cspec_uraw', binner.calc_vel_cpsd(datnow.uraw), 'spec')
        datbd.add_data('Cspec_u', binner.calc_vel_cpsd(datnow._u), 'spec')

        epstmp = np.zeros_like(datbd.u)
        Ntmp = 0
        for idx, frq_rng in enumerate(eps_freqs):
            if frq_rng is None:
                continue
            om_rng = frq_rng * 2 * np.pi
            N = ((om_rng[0] < datbd.omega) & (datbd.omega < om_rng[1])).sum()
            epstmp += binner.calc_epsilon_LT83(datbd.Spec[idx] - doppler_noise[idx],
                                               datbd.omega,
                                               np.abs(datbd.U),
                                               om_rng) * N
            Ntmp += N
        epstmp /= Ntmp
        # epstmp[np.abs(datbd.U) < 0.2] = np.NaN
        datbd.add_data('epsilon', epstmp, 'main')

        datbd.save(datdir + 'SMN_bindat_filt-{}.h5'.format(filt_tag))
