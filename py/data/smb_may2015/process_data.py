import numpy as np
import dolfyn.adv.api as avm
import scipy.signal as sig
from .base import datdir, load, tbx

filt_freqs = {
    'unfilt': 0.0,
    '5s': 1. / 5,
    '10s': 1. / 10,
    '30s': 1. / 30
}

fname = 'SMN'

doppler_noise = [2e-5, 2e-5, 0]

eps_freqs = np.array([[.1, 3],
                      [.1, 3],
                      [.1, 3], ])


def pre_process_adv():
    print("Pre-processing SMB-Nose mode vec file.")
    outfile = datdir + fname + '.h5'
    if tbx.checkhash(outfile, '62aa7b3370411cbe'):
        print("   hash check passed; skipping pre-processing.")
        return
    rr = avm.read_nortek(datdir + fname + '.VEC')
    rr = rr.subset(slice(*rr.props['inds_range']))
    rr.pitch, rr.roll, rr.heading = avm.rotate.orient2euler(rr)
    rr.roll -= 180
    rr.roll[rr.roll < -180] += 360
    # Clean
    avm.clean.cleanFill(rr.u, rr.u > 0.5)
    avm.clean.GN2002(rr.u)
    avm.clean.GN2002(rr.v)
    avm.clean.GN2002(rr.w)
    # Remove extraneous data
    rr.pop_data('AnaIn')
    rr.pop_data('AnaIn1')
    rr.pop_data('AnaIn2LSB')
    rr.pop_data('AnaIn2MSB')
    rr.pop_data('orientation_down')
    rr.pop_data('status')
    rr.pop_data('error')
    rr.pop_data('Count')
    rr.groups.pop('#extra')
    rr.save(outfile)


def merge_adv_bt(filt_freqs=filt_freqs):
    print("Processing SMB-Nose mode data (merging bottom-track).")

    # 4800 points is 5min at 16hz
    binner = avm.TurbBinner(4800, 16)

    bt = load('SMN-BT')

    bt_var = 'velbt'
    bt_var = 'velbt2'

    # Load the 'raw' (unprocessed) data that corresponds to the file above.
    dat = avm.load(datdir + fname + '.h5')

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
        vel_adp = motcalc.calc_velrot(np.array([-1, 0, 0]), to_earth=False)

        datmc.add_data('velbt', bt['BT_INST_interp'], 'orient')
        datmc.add_data('velbt2', bt['BT_INST_interp'] + vel_adp, 'orient')
        datmc.add_data('vel_adp', vel_adp, 'orient')
        bt_vars = {'velbt', 'velbt2'}
        datmc.rotate_vars.update(bt_vars)
        # The ADP data is in the instrument frame, so rotate it into the earth frame.
        avm.rotate.inst2earth(datmc, rotate_vars=bt_vars, force=True)
        if filt_freq > 0:
            filt = sig.butter(2, filt_freq / (datmc.fs / 2))
            datmc.velbt = sig.filtfilt(filt[0], filt[1], datmc.velbt)
            datmc.velbt2 = sig.filtfilt(filt[0], filt[1], datmc.velbt2)
        if filt_freq > 0:
            datmc.vel += datmc[bt_var]  # Add the bt to the vel
        datnow = datmc.copy()
        avm.rotate.inst2earth(datmc, reverse=True)

        datmc.save(datdir + 'SMN_mcdat_filt-{}.h5'.format(filt_tag))

        avm.rotate.earth2principal(datnow)
        datbd = binner(datnow)
        datbd.add_data('Spec_velraw', binner.psd(datnow.velraw), 'spec')
        datbd.add_data('Spec_velacc', binner.psd(datnow.velacc), 'spec')
        datbd.add_data('Spec_velrot', binner.psd(datnow.velrot), 'spec')
        datbd.add_data('Spec_velbt', binner.psd(datnow.velbt), 'spec')
        datbd.add_data('Spec_velbt2', binner.psd(datnow.velbt2), 'spec')
        if filt_freq > 0:
            umot = datnow.velacc + datnow.velrot + datnow[bt_var]
        else:
            umot = datnow.velacc + datnow.velrot
        datbd.add_data('Spec_velmot', binner.psd(umot), 'spec')

        datbd.add_data('Cspec_velmot', binner.calc_vel_cpsd(umot), 'spec')
        datbd.add_data('Cspec_velraw', binner.calc_vel_cpsd(datnow.velraw), 'spec')
        datbd.add_data('Cspec_vel', binner.calc_vel_cpsd(datnow.vel), 'spec')

        epstmp = np.zeros_like(datbd.vel)
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
