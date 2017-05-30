import data.ttm_june2014 as j14
import ptools as pt
import numpy as np
plt = pt.plt
import kaimal


flag = {}

flag['multi spec'] = True
#flag['multi spec color-vel'] = True
#flag['turb time01'] = True
#flag['epsVu01'] = True
#flag['epsVu02'] = True
#flag['multi spec norm'] = True
#flag['multi spec norm2'] = True
flag['add Kaimal'] = True

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

filtfreq = 0.03

# These come from 11 * sin(20), and 11 * sin(20)*tan(20) for w
lscale = [3.8, 3.8, 1.4]

if 'dat' not in vars():
    dat = j14.load('ttm02b-top', 'pax',
                   bin=True)

specnd = kaimal.Kaimal(np.logspace(-3, 4, 1000))
z_adv = 10

if flag.get('multi spec'):
    with pt.style['twocol']():

        velranges = [(0, 0.5),
                     (1, 1.5),
                     (2, 2.5)]

        fig, axs = pt.newfig(2101, 3, len(velranges),
                             figsize=5.4,
                             right=0.86,
                             bottom=0.08, top=0.96,
                             hspace=0.14,
                             sharex=True, sharey=True)

        for icol in range(axs.shape[1]):
            vr = velranges[icol]
            umag = np.abs(dat.u)
            inds = (vr[0] < umag) & (umag < vr[1])
            ustar2 = (dat.stress[1:] ** 2).sum(0)[inds].mean() ** 0.5
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
                    dnow = dat[v][irow, inds].mean(0) * pt.pii - n
                    if v == 'Spec_velmot':
                        _itmp = dat.freq < filtfreq
                        ax.loglog(dat.freq[~_itmp], dnow[~_itmp], **kwd)
                        kwd['linestyle'] = '--'
                        kwd.pop('label')
                        ax.loglog(dat.freq[_itmp], dnow[_itmp], **kwd)
                    else:
                        ax.loglog(dat.freq, dnow, **kwd)
                if flag['add Kaimal'] and specnd[irow] is not None:
                    f0 = np.abs(dat.U[inds]).mean() / z_adv
                    ax.plot(specnd.freq * f0, specnd[irow] * ustar2 / f0, 'c-',
                            label='Kaimal', zorder=5)
        for irow in range(axs.shape[0]):
            # The col-only loop
            axs[irow, 0].set_ylabel('$\mathrm{[m^2\,s^{-2}\,Hz^{-1}]}$')
            axs[irow, -1].text(1.05, 0.05, '$S\{%s\}$' % (pt.vel_comps[irow]),
                               ha='left', va='bottom', fontsize='large',
                               transform=axs[irow, -1].transAxes, clip_on=False)
            #                    ha='left', va='bottom', fontsize='medium',
            #                    transform=axs[irow, -1].transAxes, clip_on=False)

        axs[0, -1].legend(loc='upper left', bbox_to_anchor=[1.05, 1.0],
                          handlelength=1.4, handletextpad=0.4, borderaxespad=0,
                          prop=dict(size='medium'))
        ax.set_ylim((1e-4, 1))
        ax.set_xlim((1e-3, 5))

        fig.savefig(pt.figdir + 'SpecFig02_TTM02B-top.pdf')

if flag.get('multi spec norm'):

    with pt.style['twocol']():

        velranges = [(0, 0.5),
                     (1, 1.5),
                     (2, 2.5)]

        fig, axs = pt.newfig(2110, 3, len(velranges),
                             figsize=5,
                             right=0.86, bottom=0.1,
                             sharex=True, sharey=True)
        for icol in range(axs.shape[1]):
            vr = velranges[icol]
            umag = np.abs(dat.U)
            #umag = dat.u
            inds = (vr[0] < umag) & (umag < vr[1])
            U = umag[inds].mean()
            ustar2 = (dat.stress[1:, inds] ** 2).sum(0).mean() ** 0.5
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
                for v in ['Spec', 'Spec_velmot', 'Spec_velraw', ]:
                    # The col-row-var loop
                    kwd = vard[v].copy()
                    n = kwd.pop('noise')[irow]
                    spec = dat[v][irow, inds].mean(0) * pt.pii - n
                    f0 = U / z_adv
                    spec *= f0 / ustar2
                    freq = dat.freq / f0
                    ax.loglog(freq, spec, **kwd)
        for irow in range(axs.shape[0]):
            # The col-only loop
            axs[irow, 0].set_ylabel('$\mathrm{[m^2\,s^{-2}\,Hz^{-1}]}$')
            axs[irow, -1].text(1.04, 0.05, '$%s$' % (pt.vel_comps[irow]),
                               ha='left', va='bottom', fontsize='x-large',
                               transform=axs[irow, -1].transAxes, clip_on=False)
            #                    ha='left', va='bottom', fontsize='medium',
            #                    transform=axs[irow, -1].transAxes, clip_on=False)
            for icol in range(axs.shape[1]):
                ax = axs[irow, icol]
                if specnd[irow] is not None:
                    ax.plot(specnd.freq, specnd[irow], 'c-')
        axs[0, -1].legend(loc='upper left', bbox_to_anchor=[1.02, 1.0],
                          handlelength=1.4, handletextpad=0.4,
                          prop=dict(size='medium'))
        ax.set_ylim((1e-2, 100))
        ax.set_xlim((1e-3, 5))

        fig.savefig(pt.figdir + 'SpecFig03_TTM02B-top.pdf')

if flag.get('multi spec norm2'):
    with pt.style['twocol']():

        velranges = [(0, 0.5),
                     (1, 1.5),
                     (2, 2.5)]

        fig, axs = pt.newfig(2111, 3, len(velranges),
                             figsize=5,
                             right=0.86, bottom=0.1,
                             sharex=True, sharey=True)

        def interpavg(x, xp, fp, navg=3):
            fout = np.interp(x, xp, fp)
            xr = np.empty(len(x) + 1)
            xr[1:-1] = (x[1:] + x[:-1]) / 2
            xr[0] = x[0] - np.diff(x[:2]) / 2
            xr[-1] = x[-1] + np.diff(x[-2:]) / 2
            # c = np.zeros(len(x), dtype=np.int)
            for idx in range(len(x)):
                inds = (xr[idx] <= xp) & (xp <= xr[idx + 1])
                # c[idx] = inds.sum()
                if inds.sum() >= navg:
                    fout[idx] = fp[inds].mean()
            return fout

        ustar2 = np.sqrt((dat.stress[:2] ** 2).mean(0))
        U = np.abs(dat.U)
        noise = np.array(vard['Spec']['noise'])[:, None, None]
        specnd = ((dat.Spec * pt.pii - noise) *
                  U[None, :, None] / (z_adv * ustar2[None, :, None]))
        freqnd = dat.freq[None, None, :] / (U[None, :, None] / z_adv)
        fnd = np.logspace(-3, 0, 1000)
        snd = np.empty(list(specnd.shape[:2]) + [len(fnd)], dtype=np.float32)
        for i0 in range(snd.shape[0]):
            for i1 in range(snd.shape[1]):
                snd[i0, i1, :] = interpavg(fnd, freqnd[0, i1], specnd[i0, i1])

        for icol in range(axs.shape[1]):
            vr = velranges[icol]
            umag = np.abs(dat.u)
            inds = (vr[0] < umag) & (umag < vr[1])
            axs[-1, icol].set_xlabel('$fz/U$')
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
                ax.loglog(fnd, snd[irow, inds].mean(0), 'b-')
        for irow in range(axs.shape[0]):
            # The col-only loop
            axs[irow, 0].set_ylabel('$\mathrm{[m^2\,s^{-2}\,Hz^{-1}]}$')
            axs[irow, -1].text(1.04, 0.05, '$%s$' % (pt.vel_comps[irow]),
                               ha='left', va='bottom', fontsize='x-large',
                               transform=axs[irow, -1].transAxes, clip_on=False)
            #                    ha='left', va='bottom', fontsize='medium',
            #                    transform=axs[irow, -1].transAxes, clip_on=False)

        axs[0, -1].legend(loc='upper left', bbox_to_anchor=[1.02, 1.0],
                          handlelength=1.4, handletextpad=0.4,
                          prop=dict(size='medium'))
        ax.set_ylim((1e-4, 1))
        ax.set_xlim((1e-3, 5))

        # fig.savefig(pt.figdir + 'SpecFig04_TTM02B-top.pdf')

if flag.get('multi spec color-vel'):
    with pt.style['onecol']():

        velranges = [(0, 0.5),
                     (1, 1.5),
                     (2, 2.5)]

        fig, axs = pt.newfig(2301, 3, 1,
                             figsize=5,
                             right=0.7, bottom=0.1,
                             sharex=True, sharey=True, squeeze=False)

        cmap = plt.get_cmap('jet')

        for icol in range(len(velranges)):
            vr = velranges[icol]
            umag = np.abs(dat.u)
            inds = (vr[0] < umag) & (umag < vr[1])
            axs[-1, 0].set_xlabel('$f\ \mathrm{[Hz]}$')
            if vr[0] == 0:
                label = r"$ |\bar{u}| < %0.1f$" % vr[1]
            else:
                label = r"$%0.1f < |\bar{u}| < %0.1f$" % vr
            # axs[0, icol].text(.9, .9, 'N={}'.format(inds.sum()),
            #                   ha='right', va='top', fontsize='medium',
            #                   transform=axs[0, icol].transAxes)
            for irow in range(axs.shape[0]):
                # The col-row loop
                ax = axs[irow, 0]
                for fctr in [1, 1e-2, 1e-4, 1e-6, 1e-8]:
                    ax.loglog(*pt.powline(factor=fctr), linewidth=0.6,
                              linestyle=':', zorder=-6, color='k')
                for v in ['Spec']:
                    # The col-row-var loop
                    kwd = vard[v].copy()
                    kwd['label'] = label
                    n = kwd.pop('noise')[irow]
                    kwd['color'] = cmap(float(icol) / (len(velranges) - 1))
                    ax.loglog(dat.freq, dat[v][irow, inds].mean(0) * pt.pii - n,
                              **kwd)
        for irow in range(axs.shape[0]):
            # The col-only loop
            axs[irow, 0].set_ylabel('$\mathrm{[m^2\,s^{-2}\,Hz^{-1}]}$')
            axs[irow, -1].text(1.04, 0.05, '$%s$' % (pt.vel_comps[irow]),
                               ha='left', va='bottom', fontsize='x-large',
                               transform=axs[irow, -1].transAxes, clip_on=False)
            #                    ha='left', va='bottom', fontsize='medium',
            #                    transform=axs[irow, -1].transAxes, clip_on=False)

        axs[0, -1].legend(loc='upper left', bbox_to_anchor=[1.02, 1.0],
                          handlelength=1.4, handletextpad=0.4,
                          prop=dict(size='medium'))
        ax.set_ylim((1e-4, 1))
        ax.set_xlim((1e-3, 5))

        fig.savefig(pt.figdir + 'SpecFig04_TTM02B-top.pdf')

