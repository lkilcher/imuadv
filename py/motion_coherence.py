import numpy as np
import ptools as pt
from scipy.io import loadmat
import scipy.signal as sig
import matplotlib.pyplot as plt
import ttm.sm2015 as data_api
import dolfyn.adv.api as avm

flag = {}
flag['save fig'] = True
#flag['bt_basic_time'] = True
#flag['bt_filt_time'] = True
#flag['nofilt spec'] = True
flag['show cohere'] = True
#flag['phase'] = True

filt_freqs = {'unfilt': 0.0,
              '5s': 1. / 5,
              # '10s': 1. / 10,
              #'30s': 1. / 30,
              #'60s': 1. / 60,
              '5m': 1. / (5 * 60),
}

pii = np.pi * 2

# 4800 points is 5min at 16hz
binner = avm.TurbBinner(4800, 16)

if 'bt' not in vars():
    # This file is from Sam Harding.
    bt = loadmat(data_api.package_root + 'ADCP/BT_IMU_sync_data.mat')
    bt['t_IMU'] = bt['t_IMU'][0] - 366

bt_var = 'ubt'
bt_var = 'ubt2'

if 'dat' not in vars():
    # Load the 'raw' (unprocessed) data that corresponds to the file above.
    dat = data_api.load('SM_Nose', coordsys='raw')

    offind = np.abs(bt['t_IMU'][0] - dat.mpltime).argmin()

    # Perform basic motion correction based on IMU data:
    bindat_filt = {}
    dat_filt = {}
    for filt_tag in np.sort(filt_freqs.keys()):
        filt_freq = filt_freqs[filt_tag]
        datmc = dat.copy()
        avm.motion.correct_motion(datmc, accel_filtfreq=filt_freq)

        # Crop the data so that it is consistent with the ADP data.
        datmc = datmc.subset(slice(offind, (offind + len(bt['t_IMU']))))

        # Calculate 'rotational' motion of ADP head.
        motcalc = avm.motion.CalcMotion(datmc)
        ur_adp = motcalc.calc_urot(np.array([-1, 0, 0]), to_earth=False)

        datmc.add_data('ubt', bt['BT_INST_interp'], 'orient')
        datmc.add_data('ubt2', bt['BT_INST_interp'] + ur_adp, 'orient')
        bt_vars = {'ubt', 'ubt2'}
        datmc.rotate_vars.update(bt_vars)
        # The ADP data is in the instrument frame, so rotate it the earth frame.
        avm.rotate.inst2earth(datmc, rotate_vars=bt_vars, force=True)
        if filt_freq > 0:
            filt = sig.butter(2, filt_freq / (datmc.fs / 2))
            datmc.ubt = sig.filtfilt(filt[0], filt[1], datmc.ubt)
            datmc.ubt2 = sig.filtfilt(filt[0], filt[1], datmc.ubt2)
        if filt_freq > 0:
            datmc._u += datmc[bt_var]  # Add the bt to the vel
        else:
            # The principal angle for 'unfilt' is garbage:
            datmc.props['principal_angle'] = dat_filt['5s'].props['principal_angle']
        avm.rotate.earth2principal(datmc)
        dat_filt[filt_tag] = datmc

        datnow = datmc
        datbd = binner(datnow)
        datbd.add_data('Spec_uraw', binner.psd(datnow.uraw), 'Spec')
        datbd.add_data('Spec_uacc', binner.psd(datnow.uacc), 'Spec')
        datbd.add_data('Spec_urot', binner.psd(datnow.urot), 'Spec')
        datbd.add_data('Spec_ubt', binner.psd(datnow.ubt), 'Spec')
        if filt_freq > 0:
            datbd.add_data('Spec_umot',
                           binner.psd(datnow.uacc +
                                      datnow.urot +
                                      datnow[bt_var]),
                           'Spec')
        else:
            datbd.add_data('Spec_umot',
                           binner.psd(datnow.uacc +
                                      datnow.urot),
                           'Spec')
        bindat_filt[filt_tag] = datbd

    dat = dat.subset(slice(offind, (offind + len(bt['t_IMU']))))


def offset_scale(dat, offset, scale, ):
    return (dat - offset) * scale


def within(dat, minval, maxval):
    return (minval < dat) & (dat < maxval)

dnow_name = 'unfilt'
dnow_name = '5m'
dnow = dat_filt[dnow_name]
dnow.umot = dnow.urot + dnow.uacc
if dnow_name != 'unfilt':
    dnow.ubt = dat_filt['unfilt'].ubt
    dnow.ubt2 = dat_filt['unfilt'].ubt2
if dnow_name == 'unfilt':
    inds = within(np.abs(bindat_filt['5s'].U), 1.0, 1.5)
else:
    inds = within(np.abs(bindat_filt[dnow_name].U), 1.0, 1.5)


#n_fft = binner.n_bin / 4
n_fft = binner.n_bin

ubt = dnow.ubt
#ubt = dnow.ubt2


if flag.get('show cohere', False):

    cpsd = binner.cpsd(dnow.umot, ubt, n_fft=n_fft)
    sp_bt = np.abs(binner.cpsd(ubt, ubt, n_fft=n_fft))
    sp_mot = np.abs(binner.cpsd(dnow.umot, dnow.umot, n_fft=n_fft))

    cpsdr = binner.cpsd(dnow.urot, ubt, n_fft=n_fft)
    sp_rot = np.abs(binner.cpsd(dnow.urot, dnow.urot, n_fft=n_fft))

    cpsda = binner.cpsd(dnow.uacc, ubt, n_fft=n_fft)
    sp_acc = binner.cpsd(dnow.uacc, dnow.uacc, n_fft=n_fft)
    cohfreq = binner.calc_omega() / pii

    # coh = ((np.abs(cpsd[:, inds]) ** 2).mean(1) /
    #        (sp_mot[:, inds].mean(1) * sp_bt[:, inds].mean(1)))
    coh = ((np.abs(cpsd[:, inds].mean(1)) ** 2) /
           (sp_mot[:, inds].mean(1) * sp_bt[:, inds].mean(1)))
    cohr = ((np.abs(cpsdr[:, inds].mean(1)) ** 2) /
            (sp_rot[:, inds].mean(1) * sp_bt[:, inds].mean(1)))
    coha = ((np.abs(cpsda[:, inds].mean(1)) ** 2) /
            (sp_acc[:, inds].mean(1) * sp_bt[:, inds].mean(1)))

    with pt.style['twocol']():

        fig, axs = pt.newfig(300, 3, 1, figsize=8,
                             sharex=True, sharey=True,
                             right=0.75, left=0.1,
                             top=0.96, bottom=0.06,
                             hspace=0.08)

        for iax, ax in enumerate(axs):
            ax.semilogx(cohfreq, coh[iax],
                        'k', label='$u_{mot}$', lw=2, zorder=5)
            ax.semilogx(cohfreq, coha[iax],
                        'b', label='$u_{acc}$')
            ax.semilogx(cohfreq, cohr[iax],
                        'g', label='$u_{rot}$')
            ax.text(0.02, 0.98, ['u', 'v', 'w'][iax] + '-component',
                    transform=ax.transAxes,
                    ha='left', va='top')

        axs[0].set_title('Coherence with $u_{bt}$ (%s coord sys)' %
                         dnow.props['coord_sys'], size='large')
        axs[0].legend(loc='upper left', bbox_to_anchor=[1.02, 1])
        ax.set_ylim([0, 1])
        ax.set_xlim([1e-3, 1])
        ax.set_xlabel('freq [hz]')
        if flag.get('save fig', False):
            fig.savefig(pt.figdir + 'BT_IMU_Coherence01.pdf')


    with pt.style['onecol']():

        fig, ax = pt.newfig(200, 1, 1, figsize=3,
                            sharex=True, sharey=True,
                            right=0.95, left=0.19,
                            top=0.86, bottom=0.16,
                            hspace=0.08)

        ax.semilogx(cohfreq, coh[0],
                    'k', label='$u$',
                    lw=2, zorder=5, )
        ax.semilogx(cohfreq, coh[1],
                    'b', label='$v$',
                    lw=2, zorder=3, )
        ax.semilogx(cohfreq, coh[2],
                    'r', label='$w$',
                    lw=2, zorder=1)

        ax.axhline(np.sqrt(6. / (4 * inds.sum())), linestyle=':', color='k')
        ax.legend(loc='upper left')
        ax.set_ylim([0, 1])
        ax.set_xlim([1e-3, 1])
        ax.set_xlabel('$f\ \mathrm{[Hz]}$')
        ax.set_title('Coherence between IMU motion\nand Doppler profiler bottom track')
        if flag.get('save fig', False):
            fig.savefig(pt.figdir + 'BT_IMU_Coherence02.pdf')


if flag.get('phase', False):

    ph = binner.phase_angle(dnow.umot, ubt, n_fft=n_fft)[:, inds].mean(1)
    phr = binner.phase_angle(dnow.urot, ubt, n_fft=n_fft)[:, inds].mean(1)
    pha = binner.phase_angle(dnow.uacc, ubt, n_fft=n_fft)[:, inds].mean(1)

    fig = plt.newfig(301, 3, 1, figsize=8,
                     sharex=True, sharey=True,
                     right=0.75, left=0.1,
                     top=0.96, bottom=0.06,
                     hspace=0.08, )

    for iax, ax in enumerate(axs):
        ax.semilogx(cohfreq, np.angle(ph[iax]),
                    'k', label='$u_{mot}$', lw=2, zorder=5)
        ax.semilogx(cohfreq, np.angle(pha[iax]),
                    'b', label='$u_{acc}$', lw=1, zorder=1)
        ax.semilogx(cohfreq, np.angle(phr[iax]),
                    'g', label='$u_{rot}$', lw=1, zorder=0)
        # Need to confirm how to calculate angle uncertainty...
        ax.fill_between(cohfreq,
                        (np.angle(ph[iax])) - (1 - np.abs(ph[iax])) * np.pi,
                        (np.angle(ph[iax])) + (1 - np.abs(ph[iax])) * np.pi,
                        color='0.7', zorder=-5)
        ax.text(0.02, 0.98, ['u', 'v', 'w'][iax] + '-component',
                transform=ax.transAxes,
                ha='left', va='top')
        ax.yaxis.grid(True)

    axs[0].set_title('Phase shift with $u_{bt}$ (%s coord sys)' %
                     dnow.props['coord_sys'], size='large')
    axs[0].legend(loc='upper left', bbox_to_anchor=[1.02, 1])
    ax.set_yticks(np.arange(-np.pi, np.pi + 1, np.pi / 2))
    ax.set_yticklabels(['$-\pi$', '$-\pi/2$', '$0$', '$\pi/2$', '$\pi$'])
    ax.set_ylim([-np.pi, np.pi])
    ax.set_xlim([1e-3, 1])
    ax.set_xlabel('freq [hz]')
    if flag.get('save fig', False):
        fig.savefig(pt.figdir + 'BT_IMU_Phase01.pdf')

if flag.get('nofilt spec', False):

    bnow = bindat_filt['unfilt'].copy()
    bnow = bindat_filt['5m'].copy()
    bnow['Spec_ubt'] = bindat_filt['unfilt']['Spec_ubt']

    fig = pt.newfig(400, 3, 1, figsize=8,
                    sharex=True, sharey=True,
                    right=0.75,
                    left=0.13,
                    top=0.96,
                    bottom=0.06,
                    hspace=0.08, )

    for iax, ax in enumerate(axs):
        # ax.loglog(bnow.freq,
        #           (bnow.Spec[iax][inds].mean(0) - 2e-5) * pii,
        #           'g', label='$u$', linewidth=2, zorder=10)
        ax.loglog(bnow.freq,
                  (bnow.Spec_uraw[iax][inds].mean(0) - 2e-5) * pii,
                  'g', label='$u_{raw}$')
        ax.loglog(bnow.freq,
                  bnow.Spec_umot[iax][inds].mean(0) * pii,
                  'k', label='$u_{mot}$', zorder=8, linewidth=1.5)
        ax.loglog(bnow.freq,
                  bnow.Spec_urot[iax][inds].mean(0) * pii,
                  'm', label='$u_{rot}$')
        ax.loglog(bnow.freq,
                  bnow.Spec_uacc[iax][inds].mean(0) * pii,
                  'b', label='$u_{acc}$')
        ax.loglog(bnow.freq,
                  bnow.Spec_ubt[iax][inds].mean(0) * pii,
                  'r', label='$u_{bt}$')
        ax.text(0.98, 0.98, ['u', 'v', 'w'][iax],
                ha='right', va='top', transform=ax.transAxes)

    ax.set_ylabel('$\mathrm{[m^2/s^2/hz]}$')
    ax.set_ylim([1e-5, 1e1])
    ax.set_xlabel('freq [hz]')
    axs[0].legend(loc='upper left', bbox_to_anchor=[1.02, 1])
    axs[0].set_title('Velocity Spectra (%s coord sys)' %
                     dnow.props['coord_sys'], size='large')

    if flag.get('save fig', False):
        fig.savefig(pt.figdir + 'NoMC_Spectra01.pdf')
