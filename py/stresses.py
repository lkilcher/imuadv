import ttm.June2014 as j14
import make_SMdata_wBT as smdat
import ptools as pt
import dolfyn.adv.api as avm
import numpy as np
from scipy.integrate import cumtrapz
import kaimal
plt = pt.plt

flag = {}
# flag['amp/phase'] = True
# flag['real/imag'] = True
#flag['multi-real'] = True
#flag['multi-real-vp'] = True
#flag['multi-real-vp1'] = True
flag['uw-real'] = True
#flag['multi-ogive'] = True
flag['save fig'] = True

binners = dict(ttm=avm.TurbBinner(9600, 32),
               sm=avm.TurbBinner(4800, 16))
pii = 2 * np.pi
kappa = 0.4

if 'bdat' not in vars():
    bdat = {}
    bdat['ttm'] = j14.load('ttm02b-top', 'pax', bin=True)
    bdat['ttmb'] = j14.load('ttm02b-bottom', 'pax', bin=True)
    bdat['sm'] = smdat.load('SMN-5s', bindat=True)
    bdat['ttt'] = avm.load('/Users/lkilcher/data/pnnl/TTT_Vector_Feb2011_pax_b5m.h5')
    for nnm, onm in [('Cspec_vel', 'Cspec_u'),
                     ('Cspec_velrot', 'Cspec_urot'),
                     ('Cspec_velraw', 'Cspec_uraw'),
                     ('Cspec_velacc', 'Cspec_uacc'),
                     ('Cspec_velmot', 'Cspec_umot'), ]:
        if onm in bdat['sm']:
            bdat['sm'].add_data(nnm,
                                bdat['sm'].pop_data(onm),
                                'spec')
    bdat['ttm'].dudz = (bdat['ttm'].u - bdat['ttmb'].u) / 0.5
    bdat['ttm'].dvdz = (bdat['ttm'].v - bdat['ttmb'].v) / 0.5
    bdat['ttm'].ustar2 = np.sqrt(bdat['ttm'].upwp_ ** 2 +
                                 bdat['ttm'].vpwp_ ** 2)[:, None]
    bdat['ttm'].S2 = bdat['ttm'].dudz ** 2 + bdat['ttm'].dvdz ** 2
    bdat['ttm'].z = 10
    bdat['sm'].z = 10
    bdat['ttt'].z = 4.6

pairs = [(0, 1), (0, 2), (1, 2)]

do_data = [
    'ttm',
    # 'sm',
    #'ttt',
]

for idat, dat_nm in enumerate(do_data):
    bd = bdat[dat_nm]
    if flag.get('amp/phase'):

        pt.twocol()
        fig, axs = pt.newfig(1000 * idat + 101,
                             3, 2, figsize=5,
                             sharex=True, sharey='col',
                             right=0.8, bottom=0.08)

        velrng = (1, 1.5)

        print('\n{: >4} Stresses:   spec      mean     <spec:uncorr>'.format(dat_nm.upper()))
        inds = (velrng[0] < bd['u']) & (bd['u'] < velrng[1])
        f = bd.freq
        for irow, axrow in enumerate(axs):
            dtmp = bd.Cspec_vel[irow][inds].mean(0) * pii
            axrow[0].loglog(f, np.abs(dtmp), 'b-', label=pt.latex['ue'])
            axrow[1].semilogx(f, np.angle(dtmp), 'b-', label=pt.latex['ue'])
            corr = np.trapz(dtmp, f)
            dtmp = bd.Cspec_velmot[irow][inds].mean(0) * pii
            axrow[0].loglog(f, np.abs(dtmp), 'r-', label=pt.latex['uhead'])
            axrow[1].semilogx(f, np.angle(dtmp), 'r-', label=pt.latex['uhead'])
            dtmp = bd.Cspec_velraw[irow][inds].mean(0) * pii
            axrow[0].loglog(f, np.abs(dtmp), 'k-', zorder=-5, label=pt.latex['umeas'])
            axrow[1].semilogx(f, np.angle(dtmp), 'k-', zorder=-5, label=pt.latex['umeas'])
            uncorr = np.trapz(dtmp, f)
            print(r"       {}{}     {: 8.3f}  {: 8.3f}  {: 8.3f}"
                  .format(pt.vel_comps[pairs[irow][0]],
                          pt.vel_comps[pairs[irow][1]],
                          bd.stress[irow][inds].mean(0) * 1000,
                          corr.real * 1000,
                          uncorr.real * 1000))
            ax = axrow[0]
            ax.text(0.03, 0.97,
                    r'$A\{%s,%s\}$' % (pt.vel_comps[pairs[irow][0]],
                                       pt.vel_comps[pairs[irow][1]]),
                    ha='left', va='top',
                    transform=ax.transAxes)
            ax = axrow[1]
            ax.text(0.03, 0.97,
                    r'$\phi\{%s,%s\}$' % (pt.vel_comps[pairs[irow][0]],
                                          pt.vel_comps[pairs[irow][1]]),
                    ha='left', va='top',
                    transform=ax.transAxes)

        axs[0, -1].legend(loc='upper left', bbox_to_anchor=[1.1, 1])

        axs[0, 0].set_ylim([1e-6, 1])
        axs[0, 1].set_ylim([-np.pi, np.pi])
        for ax in axs[:, 1]:
            ax.yaxis.set_ticks_position('right')
        axs[0, 0].set_xlim(1e-3, 5)
        for ax in axs[-1, :]:
            ax.set_xlabel(r'$f\ \mathrm{[Hz]}$')
        for irow in range(axs.shape[0]):
            axs[irow, 0].set_ylabel('$\mathrm{m^2s^{-2}/Hz}$')

        if flag.get('save fig'):
            fig.savefig(pt.figdir + 'StressSpec_{}_01.pdf'.format(dat_nm.upper()))

    if flag.get('real/imag'):

        pt.twocol()
        fig, axs = pt.newfig(1000 * idat + 102,
                             3, 2, figsize=5,
                             sharex=True, sharey=True,
                             right=0.8, bottom=0.08)

        inds = (velrng[0] < bd['u']) & (bd['u'] < velrng[1])
        f = bd.freq
        for irow, axrow in enumerate(axs):
            dtmp = bd.Cspec_vel[irow][inds].mean(0) * pii
            axrow[0].semilogx(f, dtmp.real, 'b-', label=pt.latex['ue'])
            axrow[1].semilogx(f, dtmp.imag, 'b-', label=pt.latex['ue'])
            dtmp = bd.Cspec_velmot[irow][inds].mean(0) * pii
            axrow[0].semilogx(f, dtmp.real, 'r-', label=pt.latex['uhead'])
            axrow[1].semilogx(f, dtmp.imag, 'r-', label=pt.latex['uhead'])
            dtmp = bd.Cspec_velraw[irow][inds].mean(0) * pii
            axrow[0].semilogx(f, dtmp.real, 'k-', zorder=-5, label=pt.latex['umeas'])
            axrow[1].semilogx(f, dtmp.imag, 'k-', zorder=-5, label=pt.latex['umeas'])
            ax = axrow[0]
            ax.text(0.03, 0.97,
                    r'$C_R\{%s,%s\}$' % (pt.vel_comps[pairs[irow][0]],
                                         pt.vel_comps[pairs[irow][1]]),
                    ha='left', va='top',
                    transform=ax.transAxes)
            ax = axrow[1]
            ax.text(0.03, 0.97,
                    r'$C_I\{%s,%s\}$' % (pt.vel_comps[pairs[irow][0]],
                                         pt.vel_comps[pairs[irow][1]]),
                    ha='left', va='top',
                    transform=ax.transAxes)

        axs[0, -1].legend(loc='upper left', bbox_to_anchor=[1.1, 1])

        axs[0, 0].set_ylim([-0.2, 0.2])
        # axs[0, 1].set_ylim([-np.pi, np.pi])
        axs[0, 0].set_xlim(1e-3, 5)
        for ax in axs[-1, :]:
            ax.set_xlabel(r'$f\ \mathrm{[Hz]}$')
        for irow in range(axs.shape[0]):
            axs[irow, 0].set_ylabel('$\mathrm{m^2s^{-2}/Hz}$')

        if flag.get('save fig'):
            fig.savefig(pt.figdir + 'StressSpec_{}_02.pdf'.format(dat_nm.upper()))

    if flag.get('multi-real'):

        velranges = [(0, 0.5),
                     (1, 1.5),
                     (2, 2.5)]

        pt.twocol()
        fig, axs = pt.newfig(1000 * idat + 103,
                             3, len(velranges),
                             figsize=5,
                             sharex=True, sharey=True,
                             right=0.8, left=0.11,
                             bottom=0.08)

        f = bd.freq
        for icol, vr in enumerate(velranges):
            axcol = axs[:, icol]
            inds = (vr[0] < bd['u']) & (bd['u'] < vr[1])
            axs[0, icol].text(.9, .9, 'N={}'.format(inds.sum()),
                              ha='right', va='top', fontsize='medium',
                              transform=axs[0, icol].transAxes)
            axs[0, icol].set_title(r"$%0.1f < \bar{u} < %0.1f$" % vr,
                                   fontsize='medium')
            for irow, ax in enumerate(axcol):
                # r"$\overline{%s'%s'}=$%0.0e" % (pt.vel_comps[pairs[irow][0]],
                #                                 pt.vel_comps[pairs[irow][1]],
                #                                 np.trapz(dtmp, f), ),
                dtmp = bd.Cspec_vel[irow][inds] * pii
                # rng = np.empty((2, dtmp.shape[1]))
                # for ifrq in range(dtmp.shape[1]):
                #     rng[0, ifrq], mid, rng[1, ifrq] = pt.boot(dtmp[:, ifrq])
                rng = dtmp.std(0)[None, :] * np.array([[-1], [1]]) + dtmp.mean(0)
                ax.semilogx(f, dtmp.real.mean(0), 'b-', label=pt.latex['ue'].cspec)
                ax.fill_between(f, rng[0], rng[1], zorder=-12,
                                facecolor=[0, 0, 1, 0.3], edgecolor='none')
                # ax.semilogx(f, dtmp.T * pii, '-', color='0.7', zorder=-10)
                ax.text(.9, .1,
                        r"%2.2g" % (np.trapz(dtmp.mean(0), f) * 10000, ),
                        ha='right', va='bottom',
                        transform=ax.transAxes)
                dtmp = bd.Cspec_velmot[irow][inds].mean(0) * pii
                ax.semilogx(f, dtmp.real, 'r-', label=pt.latex['uhead'].cspec, zorder=-2)
                dtmp = bd.Cspec_velraw[irow][inds].mean(0) * pii
                ax.semilogx(f, dtmp.real, 'k-', zorder=-5, label=pt.latex['umeas'].cspec)
                ax.axhline(0, color='k', linestyle=':', zorder=-5, lw=1)

        axs[0, -1].legend(loc='upper left', bbox_to_anchor=[1.1, 1],
                          handletextpad=0.2, handlelength=2)

        axs[0, 0].set_ylim([-0.2, 0.2])
        # axs[0, 1].set_ylim([-np.pi, np.pi])
        axs[0, 0].set_xlim(1e-3, 5)
        for icol in range(axs.shape[1]):
            axs[-1, icol].set_xlabel(r'$f\ \mathrm{[Hz]}$')
        for irow in range(axs.shape[0]):
            axs[irow, 0].set_ylabel('$\mathrm{m^2s^{-2}/Hz}$')
            axs[irow, -1].text(1.03, 0.03,
                               r'$C\{%s,%s\}$' % (pt.vel_comps[pairs[irow][0]],
                                                  pt.vel_comps[pairs[irow][1]]),
                               ha='left', va='bottom',
                               transform=axs[irow, -1].transAxes)

        if flag.get('save fig'):
            fig.savefig(pt.figdir + 'StressSpec_{}_03.pdf'.format(dat_nm.upper()))

    if flag.get('multi-real-vp'):

        velranges = [(0, 0.5),
                     (1, 1.5),
                     (2, 2.5)]

        pt.twocol()
        fig, axs = pt.newfig(1000 * idat + 113,
                             3, len(velranges),
                             figsize=5,
                             sharex=True, sharey=True,
                             right=0.8, left=0.11,
                             bottom=0.08)

        f = bd.freq
        for icol, vr in enumerate(velranges):
            axcol = axs[:, icol]
            inds = (vr[0] < bd['u']) & (bd['u'] < vr[1])
            axs[0, icol].text(.9, .9, 'N={}'.format(inds.sum()),
                              ha='right', va='top', fontsize='medium',
                              transform=axs[0, icol].transAxes)
            axs[0, icol].set_title(r"$%0.1f < \bar{u} < %0.1f$" % vr,
                                   fontsize='medium')
            for irow, ax in enumerate(axcol):
                # r"$\overline{%s'%s'}=$%0.0e" % (pt.vel_comps[pairs[irow][0]],
                #                                 pt.vel_comps[pairs[irow][1]],
                #                                 np.trapz(dtmp, f), ),
                dtmp = bd.Cspec_vel[irow][inds].real * pii * f * 1000
                # rng = np.empty((2, dtmp.shape[1]))
                # for ifrq in range(dtmp.shape[1]):
                #     rng[0, ifrq], mid, rng[1, ifrq] = pt.boot(dtmp[:, ifrq])
                ax.semilogx(f, dtmp.mean(0), 'b-',
                            label=pt.latex['ue'].cspec_vp)
                rng = (dtmp.std(0)[None, :] * np.array([[-1], [1]]) + dtmp.mean(0))
                ax.fill_between(f, rng[0], rng[1], zorder=-12,
                                facecolor=[0, 0, 1, 0.3], edgecolor='none')
                # ax.semilogx(f, dtmp.T * pii, '-', color='0.7', zorder=-10)
                ax.text(.94, .06,
                        r"%2.1f" % (np.trapz(dtmp.mean(0) / f, f), ),
                        ha='right', va='bottom',
                        transform=ax.transAxes)
                dtmp = bd.Cspec_velmot[irow][inds].real * pii * f * 1000
                ax.semilogx(f, dtmp.mean(0), 'r-',
                            label=pt.latex['uhead'].cspec_vp, zorder=-2)
                dtmp = bd.Cspec_velraw[irow][inds].real * pii * f * 1000
                ax.semilogx(f, dtmp.mean(0), 'k-', zorder=-5,
                            label=pt.latex['umeas'].cspec_vp)
                ax.axhline(0, color='k', linestyle=':', zorder=-5, lw=1)

            Uhor = np.abs(bd.U[inds]).mean()
            ustar2 = np.sqrt(bd.upwp_[inds]**2 + bd.vpwp_[inds]**2).mean()
            kml = kaimal.Kaimal(f * bd.z / Uhor)
            # axs[1, icol].plot(f, kml.Suw(), 'm-', lw=5)

        axs[0, -1].legend(loc='upper left', bbox_to_anchor=[1.05, 1],
                          handletextpad=0.2, handlelength=1.4)

        axs[0, 0].set_ylim([-5, 5])
        # axs[0, 1].set_ylim([-np.pi, np.pi])
        axs[0, 0].set_xlim(1e-3, 5)
        for icol in range(axs.shape[1]):
            axs[-1, icol].set_xlabel(r'$f\ \mathrm{[Hz]}$')
        for irow in range(axs.shape[0]):
            # axs[irow, 0].set_ylabel('$f \cdot C\{%s,%s\} '
            #                         '[10^{-3}\mathrm{m^2s^{-2}}]$' %
            #                         (pt.vel_comps[pairs[irow][0]],
            #                          pt.vel_comps[pairs[irow][1]]))
            axs[irow, 0].set_ylabel('$[10^{-3}\mathrm{m^2s^{-2}}]$')
            # axs[irow, -1].text(1.06, 0.03,
            #                    r'$f \cdot C\{%s,%s\}$' % (
            #                        pt.vel_comps[pairs[irow][0]],
            #                        pt.vel_comps[pairs[irow][1]]),
            #                    ha='left', va='bottom',
            #                    transform=axs[irow, -1].transAxes)
            axs[irow, 0].text(.06, .9,
                              r'$f \cdot C\{%s,%s\}$' % (
                                  pt.vel_comps[pairs[irow][0]],
                                  pt.vel_comps[pairs[irow][1]]),
                              ha='left', va='top',
                              transform=axs[irow, 0].transAxes)

        if flag.get('save fig'):
            fig.savefig(pt.figdir + 'StressSpec_{}_03vp.pdf'
                        .format(dat_nm.upper()))

    if flag.get('multi-real-vp1'):

        velranges = [(2, 2.5)]

        #fctr = 1000
        fctr = 100
        with pt.style['onecol']():
            fig, axs = pt.newfig(1000 * idat + 103,
                                 3, len(velranges),
                                 figsize=5.6,
                                 squeeze=False,
                                 sharex=True, sharey=True,
                                 right=0.7, left=0.15,
                                 hspace=0.15,
                                 bottom=0.08)

            f = bd.freq
            for icol, vr in enumerate(velranges):
                axcol = axs[:, icol]
                inds = (vr[0] < bd['u']) & (bd['u'] < vr[1])
                ustar2 = np.sqrt(bd.upwp_[inds]**2 + bd.vpwp_[inds]**2).mean()
                f0 = np.abs(bd.U[inds]).mean() / bd.z
                kml = kaimal.Kaimal(bd.freq / f0)

                axs[0, icol].text(.9, .9, 'N={}'.format(inds.sum()),
                                  ha='right', va='top', fontsize='medium',
                                  transform=axs[0, icol].transAxes)
                axs[0, icol].set_title(r"$%0.1f < \bar{u} < %0.1f$" % vr,
                                       fontsize='medium')
                for irow, ax in enumerate(axcol):
                    dtmp = bd.Cspec_vel[irow][inds].real * pii * f * fctr
                    ax.semilogx(f, dtmp.mean(0), 'b-',
                                label=pt.latex['ue'].cspec_vp)
                    rng = (dtmp.std(0)[None, :] * np.array([[-1], [1]]) + dtmp.mean(0))
                    ax.fill_between(f, rng[0], rng[1], zorder=-12,
                                    facecolor=[0, 0, 1, 0.3], edgecolor='none')
                    dtmp = bd.Cspec_velmot[irow][inds].real * pii * f * fctr
                    itmp = f < 0.0333
                    ax.semilogx(f[itmp], dtmp.mean(0)[itmp], 'r--',
                                zorder=-2)
                    ax.semilogx(f[~itmp], dtmp.mean(0)[~itmp], 'r-',
                                label=pt.latex['uhead'].cspec_vp, zorder=-2)
                    dtmp = bd.Cspec_velraw[irow][inds].real * pii * f * fctr
                    ax.semilogx(f, dtmp.mean(0), 'k-', zorder=-5,
                                label=pt.latex['umeas'].cspec_vp)
                    ax.axhline(0, color='k', linestyle=':', zorder=-5, lw=1)
                    ax.axvline(0.03, color='r', ls=':', zorder=-10)
                    # if irow == 1:
                    #     ax.plot(bd.freq,
                    #             -bd.freq * kml.Suw() * ustar2 / f0 * fctr,
                    #             color='c', label='Kaimal')

            axs[0, -1].legend(loc='upper left', bbox_to_anchor=[1.05, 1],
                              borderaxespad=0, prop=dict(size='medium'),
                              handletextpad=0.2, handlelength=1.3)

            axs[0, 0].set_ylim([-3, 3])
            axs[0, 0].set_xlim(1e-3, 5)
            for icol in range(axs.shape[1]):
                axs[-1, icol].set_xlabel(r'$f\ \mathrm{[Hz]}$')
            for irow in range(axs.shape[0]):
                axs[irow, 0].set_ylabel('$[10^{-2}\,\mathrm{m^2s^{-2}}]$')
                axs[irow, -1].text(1.05, .1,
                                   r'$f \cdot C\{%s,%s\}$' % (
                                       pt.vel_comps[pairs[irow][0]],
                                       pt.vel_comps[pairs[irow][1]]),
                                   ha='left', va='bottom', size='large',
                                   transform=axs[irow, 0].transAxes)

            if flag.get('save fig'):
                fig.savefig(pt.figdir + 'StressSpec_{}_04vp.pdf'
                            .format(dat_nm.upper()))

    if flag.get('uw-real'):

        case = 'ebb'
        case = 'flood'
        case = 'both'
        add_ttt = False
        if case == 'ebb':
            inds = bd.u > 1
            ttl = r'Ebb ($\bar{u} > 1\,\mathrm{ms^{-1}}$)'
        elif case == 'flood':
            inds = bd.u < -1
            ttl = r'Flood ($\bar{u} < -1\,\mathrm{ms^{-1}}$)'
        elif case == 'both':
            inds = np.abs(bd.u) > 1
            ttl = r'$|\bar{u}| > 1\,\mathrm{ms^{-1}}$'
        else:
            raise Exception("Invalid case.")

        # plt.figure(65)
        # plt.clf()
        # plt.hist(np.log10(bd.u[inds] ** 2 / bd.dudz[inds] ** 2))
        # factor = (bd.u[inds] ** 2 / bd.dudz[inds] ** 2).mean()
        # utmp = np.array([1, 3])
        # # plt.loglog(bd.u[inds] ** 2, bd.dudz[inds] ** 2, 'k.')
        # # plt.loglog(utmp ** 2, utmp ** 2 / factor, 'k--')

        m73f = np.array([1e-3, 10])
        m73a = 0.1 * m73f ** (-7. / 3)
        df = 2e-2
        fstar_bins = np.arange(df, 10, df)
        kml = kaimal.Kaimal(fstar_bins)
        #fstar_bins = np.logspace(-2, 1, 50)
        Cdata, fstar = kaimal.nd_cospec(bd, inds)
        Cs = kaimal.bin_cospec(Cdata, fstar, fstar_bins)
        Cs2 = Cs * fstar_bins
        if add_ttt:
            #z_ttt = 4.6
            #z_ttt = 10  # This gives better agreement, but why?
            # Is this b/c of shear?
            tmp = bdat['ttt']
            tmp_inds = np.abs(bdat['ttt'].u) > 1
            Cdttt, fsttt = kaimal.ndcs(tmp, tmp_inds)
            Cs_ttt = kaimal.bin_cospec(Cdttt, fsttt, fstar_bins)
            Cs2_ttt = Cs_ttt * fstar_bins

        newfig_kws = dict(figsize=5,
                          sharex=True, sharey=False,
                          right=0.95, left=0.2,
                          bottom=0.09, top=0.94,
                          hspace=0.13)
        pt_kws = dict(color='0.4', ms=3,
                      zorder=-20, alpha=0.3,
                      rasterized=True, mec='none')

        with pt.style['onecol']():

            fig, axs = pt.newfig(1000 * idat + 203,
                                 2, 1, **newfig_kws)

            ax = axs[0]
            # ax.loglog(fstar_bins, np.abs(Cs), 'b-', lw=1,
            #           label='Avg.')
            ax.loglog(fstar_bins, Cs, 'b-', lw=2,
                      label='Avg.',)
            ax.loglog(fstar.flatten(), Cdata.flatten(),
                      '.', label='single',
                      **pt_kws)
            if add_ttt:
                ax.loglog(fstar_bins, Cs_ttt, 'm-', lw=1,
                          label='TTT', zorder=-10)
            ax.loglog(m73f, m73a, 'r--', lw=2,
                      label='$\hat{f}^{-7/3}$')
            ax.loglog(fstar_bins, kml.Suw(), 'k-', lw=2, zorder=-2,
                      label='Kaimal')
            ax.set_xlim([1e-2, 10])
            ax.set_ylim([1e-4, 30])
            # z = 10
            ax.set_ylabel('$\hat{C}\{u,w\}$')
            ax.legend(loc='lower left',
                      prop=dict(size='small'), numpoints=1)
            # ax.set_title(ttl)

            ax = axs[1]
            ax.semilogx(fstar_bins, Cs2, 'b-', lw=2,
                        label='Avg.')
            ax.semilogx(fstar.flatten(),
                        fstar.flatten() * Cdata.flatten(),
                        '.', label='single', **pt_kws)
            ax.semilogx(fstar_bins, fstar_bins * kml.Suw(), 'k-',
                        label='Kaimal',
                        lw=2, zorder=-2)
            ax.set_xlim([1e-2, 10])
            ax.set_ylim([-1, 1])
            ax.set_ylabel('$\hat{f}\cdot \hat{C}\{u,w\}$')

            axs[-1].set_xlabel('$\hat{f}$')
            axs[0].set_title(ttl)

            fig.savefig(pt.figdir + 'CoSpecND02_{}-{}.pdf'.format(
                dat_nm.upper(), case),
                dpi=200)

    if flag.get('multi-ogive'):

        velranges = [(0, 0.5),
                     (1, 1.5),
                     (2, 2.5)]

        fig, axs = pt.newfig(1000 * idat + 503,
                             3, len(velranges),
                             figsize=5,
                             sharex=True, sharey=True,
                             right=0.8, left=0.11,
                             bottom=0.08)

        f = bd.freq
        for icol, vr in enumerate(velranges):
            axcol = axs[:, icol]
            inds = (vr[0] < bd['u']) & (bd['u'] < vr[1])
            axs[0, icol].text(.9, .9, 'N={}'.format(inds.sum()),
                              ha='right', va='top', fontsize='medium',
                              transform=axs[0, icol].transAxes)
            axs[0, icol].set_title(r"$%0.1f < \bar{u} < %0.1f$" % vr,
                                   fontsize='medium')
            for irow, ax in enumerate(axcol):
                ax.text(.9, .1,
                        r"%2.2g" % (np.trapz(bd.Cspec_vel[irow][inds].real.mean(0) *
                                             pii, f) * 10000, ),
                        ha='right', va='bottom',
                        transform=ax.transAxes)
                # r"$\overline{%s'%s'}=$%0.0e" % (pt.vel_comps[pairs[irow][0]],
                #                                 pt.vel_comps[pairs[irow][1]],
                #                                 np.trapz(dtmp, f), ),
                dtmp = cumtrapz(bd.Cspec_vel[irow][inds] * pii, f, initial=0).real
                ax.semilogx(f, dtmp.mean(0), 'b-', label=pt.latex['ue'])
                ax.semilogx(f, dtmp.T, '-', color='0.8', zorder=-10)
                dtmp = cumtrapz(bd.Cspec_velmot[irow][inds] * pii, f, initial=0).real
                ax.semilogx(f, dtmp.mean(0), 'r-', label=pt.latex['uhead'], zorder=-2)
                dtmp = cumtrapz(bd.Cspec_velraw[irow][inds] * pii, f, initial=0).real
                ax.semilogx(f, dtmp.mean(0), 'k-', zorder=-5, label=pt.latex['umeas'])
                ax.axhline(0, color='k', linestyle=':', zorder=-5, lw=1)
        axs[0, -1].legend(loc='upper left', bbox_to_anchor=[1.1, 1])

        axs[0, 0].set_ylim([-0.02, 0.02])
        # axs[0, 1].set_ylim([-np.pi, np.pi])
        axs[0, 0].set_xlim(1e-3, 5)
        for icol in range(axs.shape[1]):
            axs[-1, icol].set_xlabel(r'$f\ \mathrm{[Hz]}$')
        for irow in range(axs.shape[0]):
            axs[irow, 0].set_ylabel('$\mathrm{m^2s^{-2}/Hz}$')
            axs[irow, -1].text(1.03, 0.03,
                               r'$C_R\{%s,%s\}$' % (pt.vel_comps[pairs[irow][0]],
                                                    pt.vel_comps[pairs[irow][1]]),
                               ha='left', va='bottom',
                               transform=axs[irow, -1].transAxes)
            
        if flag.get('save fig'):
            fig.savefig(pt.figdir + 'StressSpec_{}_04.pdf'.format(dat_nm.upper()))

