import numpy as np
import ptools as pt
import data.smb_may2015 as smdat
import kaimal

flag = {}
#flag['bt_basic_time'] = True
#flag['bt_filt_time'] = True
#flag['bt_filt_spec'] = True
#flag['all spec'] = True
flag['multi spec'] = True
flag['save figs'] = True
flag['add Kaimal'] = True

pii = 2 * np.pi

filt_freqs = {
    #'unfilt': 0.0,
    '5s': 1. / 5,
    #'10s': 1. / 10,
    #'30s': 1. / 30
}

if 'dat_filt' not in vars():
    dat_filt = {}
    bindat_filt = {}

for ftag in filt_freqs:
    if ftag not in dat_filt:
        dat_filt[ftag] = smdat.load('SMN-{}'.format(ftag), bindat=False)
    if ftag not in bindat_filt:
        bindat_filt[ftag] = smdat.load('SMN-{}'.format(ftag), bindat=True)

if 'bt' not in vars():
    bt = smdat.load('SMN-BT')

dnow = bindat_filt['5s']
filtfreq = filt_freqs['5s']

# The mean velocity is totally contaminated by motion correction, so
# lets fix it here
if 'unfilt' in dat_filt:
    bindat_filt['unfilt']['vel'] = bindat_filt['5s']['vel']
    dat_filt['unfilt']['vel'] = dat_filt['5s']['vel']

specnd = kaimal.Kaimal(np.logspace(-3, 4, 1000))
z_adv = 10


def offset_scale(dat, offset, scale, ):
    return (dat - offset) * scale

if flag.get('bt_basic_time', False):

    datmc = dat_filt['30s']

    fig, ax = pt.newfig(301, 1, 1, )
    inds = slice(int(1.2e4), int(1.4e4))
    toff = datmc.mpltime[inds.start]

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
    if flag.get('save figs'):
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
        if flag.get('save figs'):
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

            inds = pt.within(np.abs(datbd.u), *velrange)

            for iax, ax in enumerate(axs):
                ax.loglog(datbd.freq,
                          (datbd.Spec[iax][inds].mean(0) - doppler_noise[iax]) * pii,
                          'b', label=pt.latex['ue'].spec, linewidth=1.5, zorder=10)
                ax.loglog(datbd.freq,
                          (datbd.Spec_velraw[iax][inds].mean(0) - doppler_noise[iax]) * pii,
                          'k', label=pt.latex['umeas'].spec)
                ax.loglog(datbd.freq,
                          datbd.Spec_velmot[iax][inds].mean(0) * pii,
                          'r', label=pt.latex['uhead'].spec, zorder=8, )
                # ax.loglog(datbd.freq,
                #           datbd.Spec_velrot[iax][inds].mean(0) * pii,
                #           'm', label='$u_{rot}$')
                # ax.loglog(datbd.freq,
                #           datbd.Spec_velacc[iax][inds].mean(0) * pii,
                #           'b', label='$u_{acc}$')
                # ax.loglog(datbd.freq,
                #           datbd.Spec_velbt[iax][inds].mean(0) * pii,
                #           'r', label='$u_{bt}$')
                ax.plot(line['x'], line['y'], 'k--')
                ax.axvline(filt_freq, linestyle=':', color='k')
                ax.set_xlabel('$f\ \mathrm{[Hz]}$')

            ax.set_xlim([1e-3, 2])
            ax.set_ylim([1e-4, 1])
            # axs[-1].legend(bbox_to_anchor=[1.02, 1], loc='upper left')
            axs[0].legend(loc='lower left',
                          prop=dict(size='small'))
            axs[0].set_ylabel('$\mathrm{[m^2\,s^{-2}\,Hz^{-1}]}$')
            axs[0].set_xlabel('$f\ \mathrm{[Hz]}$')

            if flag.get('save figs'):
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

        inds = pt.within(np.abs(datbd.u), *velrange)

        for iax, kwd in enumerate([dict(color='b', label='u'),
                                   dict(color='g', label='v'),
                                   dict(color='r', label='w')]):
            ax.loglog(datbd.freq,
                      (datbd.Spec[iax][inds].mean(0) - doppler_noise[iax]) * pii,
                      linewidth=2, zorder=10, **kwd)
            # ax.loglog(datbd.freq,
            #           (datbd.Spec_velraw[iax][inds].mean(0) - 2e-5) * pii,
            #           'y', label='$u_{raw}$')
            # ax.loglog(datbd.freq,
            #           datbd.Spec_velmot[iax][inds].mean(0) * pii,
            #           'k', label='$u_{mot}$', zorder=8, linewidth=1.5)
            # ax.loglog(datbd.freq,
            #           datbd.Spec_velrot[iax][inds].mean(0) * pii,
            #           'm', label='$u_{rot}$')
            # ax.loglog(datbd.freq,
            #           datbd.Spec_velacc[iax][inds].mean(0) * pii,
            #           'b', label='$u_{acc}$')
            # ax.loglog(datbd.freq,
            #           datbd.Spec_velbt[iax][inds].mean(0) * pii,
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
        ax.set_ylabel('$\mathrm{[m^2\,s^{-2}\,Hz^{-1}]}$')
        if flag.get('save figs'):
            fig.savefig(pt.figdir + 'VelSpec_filt{}.pdf'.format(filt_tag))


if flag.get('multi spec'):

    vard = dict(
        Spec_velmot=dict(color='r', lw=1.5, zorder=1,
                       label=pt.latex['uhead'].spec,
                       noise=np.zeros(3),
        ),
        Spec_velraw=dict(color='k', zorder=2,
                       label=pt.latex['umeas'].spec,
                       noise=[1.5e-4, 1.5e-4, 1.5e-5, ],

        ),
        Spec=dict(color='b', lw=1.5, zorder=3,
                  label=pt.latex['ue'].spec,
                  noise=[1.5e-4, 1.5e-4, 1.5e-5, ],
        ),
    )

    with pt.style['twocol']():

        velranges = [(0, 0.5),
                     (1, 1.5),
                     (2, 2.5)]

        fig, axs = pt.newfig(101, 3, len(velranges),
                             figsize=5.4,
                             right=0.86,
                             bottom=0.08, top=0.96,
                             hspace=0.14,
                             sharex=True, sharey=True)

        for icol in range(axs.shape[1]):
            vr = velranges[icol]
            umag = np.abs(dnow.u)
            inds = (vr[0] < umag) & (umag < vr[1])
            ustar2 = (dnow.stress[1:] ** 2).sum(0)[inds].mean() ** 0.5
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
                ax.axvline(filtfreq, linewidth=0.6,
                           linestyle=':', zorder=-6, color='r')
                for fctr in [1, 1e-2, 1e-4, 1e-6, 1e-8]:
                    ax.loglog(*pt.powline(factor=fctr), linewidth=0.6,
                              linestyle=':', zorder=-6, color='k')
                for v in ['Spec', 'Spec_velmot', 'Spec_velraw', ]:
                    # The col-row-var loop
                    kwd = vard[v].copy()
                    n = kwd.pop('noise')[irow]
                    ax.loglog(dnow.freq, dnow[v][irow, inds].mean(0) * pii - n,
                              **kwd)
                if flag['add Kaimal'] and specnd[irow] is not None:
                    umag = np.abs(dnow.U[inds]).mean()
                    f0 = umag / z_adv
                    #print ustar2, umag, f0, ustar2 / f0
                    ax.plot(specnd.freq * f0, specnd[irow] * ustar2 / f0, 'c-',
                            label='Kaimal', zorder=5)

        for irow in range(axs.shape[0]):
            # The row-only loop
            axs[irow, 0].set_ylabel('$\mathrm{[m^2\,s^{-2}\,Hz^{-1}]}$')
            axs[irow, -1].text(1.05, 0.05, '$S\{%s\}$' % (pt.vel_comps[irow]),
                               ha='left', va='bottom', fontsize='large',
                               transform=axs[irow, -1].transAxes, clip_on=False)
        axs[0, -1].legend(loc='upper left', bbox_to_anchor=[1.05, 1.0],
                          handlelength=1.4, handletextpad=0.4, borderaxespad=0,
                          prop=dict(size='medium'))
        ax.set_ylim((1e-4, 1))
        ax.set_xlim((1e-3, 5))

        if flag.get('save figs'):
            fig.savefig(pt.figdir + 'SpecFig02_SMnose.pdf')
