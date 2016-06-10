from scipy.io import loadmat
import scipy.signal as sig
import ttm.sm2015 as data_api
import numpy as np
import dolfyn.adv.api as avm
import ptools as pt

flag = {}
#flag['bt_basic_time'] = True
#flag['bt_filt_time'] = True
#flag['bt_filt_spec'] = True
#flag['all spec'] = True
flag['multi spec'] = True

filt_freqs = {
    #'unfilt': 0.0,
    '5s': 1. / 5,
    #'10s': 1. / 10,
    #'30s': 1. / 30
}

pii = 2 * np.pi

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
        dat_filt[filt_tag] = datmc

        avm.rotate.earth2principal(datnow)
        datbd = binner(datnow)
        datbd.Spec_uraw = binner.psd(datnow.uraw)
        datbd.Spec_uacc = binner.psd(datnow.uacc)
        datbd.Spec_urot = binner.psd(datnow.urot)
        datbd.Spec_ubt = binner.psd(datnow.ubt)
        datbd.Spec_ubt2 = binner.psd(datnow.ubt2)
        if filt_freq > 0:
            datbd.Spec_umot = binner.psd(datnow.uacc +
                                         datnow.urot +
                                         datnow[bt_var])
        else:
            datbd.Spec_umot = binner.psd(datnow.uacc +
                                         datnow.urot)
        bindat_filt[filt_tag] = datbd

    dat = dat.subset(slice(offind, (offind + len(bt['t_IMU']))))

doppler_noise = [2e-5, 2e-5, 0]

# The mean velocity is totally contaminated by motion correction, so
# lets fix it here
if 'unfilt' in dat_filt:
    bindat_filt['unfilt']['_u'] = bindat_filt['5s']['_u']
    dat_filt['unfilt']['_u'] = dat_filt['5s']['_u']


def offset_scale(dat, offset, scale, ):
    return (dat - offset) * scale


def within(dat, minval, maxval):
    return (minval < dat) & (dat < maxval)


if flag.get('bt_basic_time', False):

    datmc = dat_filt['30s']

    fig, ax = pt.newfig(301, 1, 1, )
    inds = slice(int(1.2e4), int(1.4e4))
    toff = dat.mpltime[inds.start]

    ax.plot(offset_scale(datmc.mpltime[inds], toff, 24 * 36000),
            datmc.uacc[0][inds],
            label='$u_{acc}$', color='r')
    ax.plot(offset_scale(datmc.mpltime[inds], toff, 24 * 36000),
            datmc.urot[0][inds],
            label='$u_{rot}$', color='g')
    ax.plot(offset_scale(datmc.mpltime[inds], toff, 24 * 36000),
            datmc.uacc[0][inds] + datmc.urot[0][inds],
            label='$u_{mot}$', color='m')
    ax.plot(offset_scale(bt['t_IMU'][inds], toff, 24 * 36000),
            bt['IMU_INST'][0][inds],
            label='$u_{mot2}$', color='y')
    ax.plot(offset_scale(bt['t_IMU'][inds], toff, 24 * 36000),
            bt['BT_INST_interp'][0][inds],
            label='$u_{bt}$', color='b')
    ax.legend()
    ax.set_xlabel('Time [sec]')

    ax.set_ylabel('[m/s]')
    fig.savefig(pt.figdir + 'BT_basic_plot01.pdf')
    with file('fig/BT_basic_plot01.caption', 'w') as fl:
        fl.write("""This shows the basic variables

        This confirms that the datasets overlap and that things
        basically look correct.  Note here that the u_bt data is
        VERY different from the u_mot data. This suggests that
        u_acc is wrong somehow, and provides the impetus for doing
        this work.
        """)
    #ax.set_xlim(datmc.mpltime[inds][[0, -1]])
    #ax.set_ylim([-0.05, 0.05])

if flag.get('bt_filt_time', False):

    for ifilt, (filt_tag, filt_freq) in enumerate(filt_freqs.iteritems()):
        datnow = dat_filt[filt_tag]

        fig, ax = pt.newfig(320 + ifilt, 1, 1)
        inds = slice(int(1.2e4), int(1.4e4))
        toff = datnow.mpltime[inds.start]

        ax.plot(offset_scale(datnow.mpltime[inds], toff, 24 * 36000),
                datnow.uacc[0][inds] + datnow.urot[0][inds],
                label='$u_{mot}$', color='y')
        ax.plot(offset_scale(datnow.mpltime[inds], toff, 24 * 36000),
                ur_adp[0][inds],
                label='$u_{adp:rot}$', color='r')
        ax.plot(offset_scale(bt['t_IMU'][inds], toff, 24 * 36000),
                bt['BT_INST_interp'][0][inds],
                label='$u_{bt}$', color='b')
        ax.plot(offset_scale(bt['t_IMU'][inds], toff, 24 * 36000),
                datnow.ubt[0][inds],
                label="$u_{bt}'$", color='b', linestyle=':')
        ax.plot(offset_scale(bt['t_IMU'][inds], toff, 24 * 36000),
                datnow.ubt2[0][inds] + datnow.uacc[0][inds] + datnow.urot[0][inds],
                label='$u_{mot2}$', color='k', linestyle='-')
        ax.legend()
        ax.set_xlabel('Time [sec]')
        ax.axhline(0, linestyle=':', color='k')

        ax.set_ylabel('[m/s]')
        ax.set_title('{} filter'.format(filt_tag))
        fig.savefig(pt.figdir + 'BT_time_filt{}.pdf'.format(filt_tag))

if flag.get('bt_filt_spec', False):

    line = {'x': np.array([1e-5, 100])}
    line['y'] = 2e-4 * line['x'] ** (-5. / 3) * pii

    with pt.twocol():
        velrange = [1.2, 1.5]
        #velrange = [0.5, 1.0]
        for ifilt, (filt_tag, filt_freq) in enumerate(filt_freqs.iteritems()):
            datbd = bindat_filt[filt_tag]

            fig, axs = pt.newfig(330 + ifilt, 1, 3, figsize=2.4,
                                 right=0.98,
                                 bottom=0.17,
                                 sharex=True, sharey=True)

            inds = within(np.abs(datbd.u), *velrange)

            for iax, ax in enumerate(axs):
                ax.loglog(datbd.freq,
                          (datbd.Spec[iax][inds].mean(0) - doppler_noise[iax]) * pii,
                          'b', label=pt.latex['ue'].spec, linewidth=1.5, zorder=10)
                ax.loglog(datbd.freq,
                          (datbd.Spec_uraw[iax][inds].mean(0) - doppler_noise[iax]) * pii,
                          'k', label=pt.latex['umeas'].spec)
                ax.loglog(datbd.freq,
                          datbd.Spec_umot[iax][inds].mean(0) * pii,
                          'r', label=pt.latex['uhead'].spec, zorder=8, )
                # ax.loglog(datbd.freq,
                #           datbd.Spec_urot[iax][inds].mean(0) * pii,
                #           'm', label='$u_{rot}$')
                # ax.loglog(datbd.freq,
                #           datbd.Spec_uacc[iax][inds].mean(0) * pii,
                #           'b', label='$u_{acc}$')
                # ax.loglog(datbd.freq,
                #           datbd.Spec_ubt[iax][inds].mean(0) * pii,
                #           'r', label='$u_{bt}$')
                ax.plot(line['x'], line['y'], 'k--')
                ax.axvline(filt_freq, linestyle=':', color='k')
                ax.set_xlabel('$f\ \mathrm{[Hz]}$')

            ax.set_xlim([1e-3, 2])
            ax.set_ylim([1e-4, 1])
            # axs[-1].legend(bbox_to_anchor=[1.02, 1], loc='upper left')
            axs[0].legend(loc='lower left',
                          prop=dict(size='small'))
            axs[0].set_ylabel('$\mathrm{[m^2s^{-2}/Hz]}$')
            axs[0].set_xlabel('$f\ \mathrm{[Hz]}$')

            fig.savefig(pt.figdir + 'SM_spec_filt{}.pdf'.format(filt_tag))


if flag.get('all spec'):

    line = {'x': np.array([1e-5, 100])}
    line['y'] = 2e-4 * line['x'] ** (-5. / 3) * pii

    velrange = [1.2, 1.5]
    #velrange = [0.5, 1.0]
    for ifilt, (filt_tag, filt_freq) in enumerate(filt_freqs.iteritems()):
        datbd = bindat_filt[filt_tag]

        fig, ax = pt.newfig(340 + ifilt, 1, 1, figsize=[5.5, 4.5],
                            gridspec_kw=dict(right=0.8,
                                             top=0.95,
                                             left=0.15,
                                             bottom=0.15,
                                             hspace=0.08),
                            sharex=True, sharey=True)

        inds = within(np.abs(datbd.u), *velrange)

        for iax, kwd in enumerate([dict(color='b', label='u'),
                                   dict(color='g', label='v'),
                                   dict(color='r', label='w')]):
            ax.loglog(datbd.freq,
                      (datbd.Spec[iax][inds].mean(0) - doppler_noise[iax]) * pii,
                      linewidth=2, zorder=10, **kwd)
            # ax.loglog(datbd.freq,
            #           (datbd.Spec_uraw[iax][inds].mean(0) - 2e-5) * pii,
            #           'y', label='$u_{raw}$')
            # ax.loglog(datbd.freq,
            #           datbd.Spec_umot[iax][inds].mean(0) * pii,
            #           'k', label='$u_{mot}$', zorder=8, linewidth=1.5)
            # ax.loglog(datbd.freq,
            #           datbd.Spec_urot[iax][inds].mean(0) * pii,
            #           'm', label='$u_{rot}$')
            # ax.loglog(datbd.freq,
            #           datbd.Spec_uacc[iax][inds].mean(0) * pii,
            #           'b', label='$u_{acc}$')
            # ax.loglog(datbd.freq,
            #           datbd.Spec_ubt[iax][inds].mean(0) * pii,
            #           'r', label='$u_{bt}$')
            ax.plot(line['x'], line['y'], 'k--')
            ax.axvline(filt_freq, linestyle=':', color='k')

        ax.set_xlim([1e-3, 2])
        ax.set_ylim([1e-5, 1])
        # if filt_tag == 'unfilt':
        #     axs[0].set_title('unfiltered'.format(filt_tag))
        # else:
        #     axs[0].set_title('{} filter'.format(filt_tag))
        ax.legend(bbox_to_anchor=[1.02, 1], loc='upper left')
        ax.set_xlabel('$f\ \mathrm{[hz]}$')
        ax.set_ylabel('$\mathrm{[m^2s^{-2}/Hz]}$')
        fig.savefig(pt.figdir + 'VelSpec_filt{}.pdf'.format(filt_tag))


if flag.get('multi spec'):

    vard = dict(
        Spec_umot=dict(color='r', lw=1.5, zorder=1,
                       label=pt.latex['uhead'].spec,
                       noise=np.zeros(3),
        ),
        Spec_uraw=dict(color='k', zorder=2,
                       label=pt.latex['umeas'].spec,
                       noise=[1.5e-4, 1.5e-4, 1.5e-5, ],

        ),
        Spec=dict(color='b', lw=1.5, zorder=3,
                  label=pt.latex['ue'].spec,
                  noise=[1.5e-4, 1.5e-4, 1.5e-5, ],
        ),
    )

    dat = bindat_filt['5s']

    with pt.style['twocol']():

        velranges = [(0, 0.5),
                     (1, 1.5),
                     (2, 2.5)]

        fig, axs = pt.newfig(101, 3, len(velranges),
                             figsize=5,
                             right=0.86, bottom=0.1,
                             sharex=True, sharey=True)

        for icol in range(axs.shape[1]):
            vr = velranges[icol]
            umag = np.abs(dat.u)
            inds = (vr[0] < umag) & (umag < vr[1])
            axs[-1, icol].set_xlabel('$f\ \mathrm{[Hz]}$')
            if vr[0] == 0:
                axs[0, icol].set_title(r"$ |\bar{u}| < %0.1f$" % vr[1],
                                       fontsize='medium')
            else:
                axs[0, icol].set_title(r"$%0.1f < |\bar{u}| < %0.1f$" % vr,
                                       fontsize='medium')
            axs[0, icol].text(.9, .9, 'N={}'.format(inds.sum()),
                              ha='right', va='top', fontsize='medium',
                              transform=axs[0, icol].transAxes)
            for irow in range(axs.shape[0]):
                # The col-row loop
                ax = axs[irow, icol]
                for fctr in [1, 1e-2, 1e-4, 1e-6, 1e-8]:
                    ax.loglog(*pt.powline(factor=fctr), linewidth=0.6,
                              linestyle=':', zorder=-6, color='k')
                for v in ['Spec', 'Spec_umot', 'Spec_uraw', ]:
                    # The col-row-var loop
                    kwd = vard[v].copy()
                    n = kwd.pop('noise')[irow]
                    ax.loglog(dat.freq, dat[v][irow, inds].mean(0) * pii - n,
                              **kwd)
        for irow in range(axs.shape[0]):
            # The row-only loop
            axs[irow, 0].set_ylabel('$\mathrm{[m^2s^{-2}/Hz]}$')
            axs[irow, -1].text(1.04, 0.05, '$%s$' % (pt.vel_comps[irow]),
                               ha='left', va='bottom', fontsize='x-large',
                               transform=axs[irow, -1].transAxes, clip_on=False)
        axs[0, -1].legend(loc='upper left', bbox_to_anchor=[1.02, 1.0],
                          handlelength=1.4, handletextpad=0.4,
                          prop=dict(size='medium'))
        ax.set_ylim((1e-4, 1))
        ax.set_xlim((1e-3, 5))

        fig.savefig(pt.figdir + 'SpecFig02_SMnose.pdf')
