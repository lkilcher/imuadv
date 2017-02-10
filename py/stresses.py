import ttm.June2014 as j14
import make_SMdata_wBT as smdat
import ptools as pt
import dolfyn.adv.api as avm
import numpy as np
from scipy.integrate import cumtrapz
plt = pt.plt

flag = {}
# flag['amp/phase'] = True
# flag['real/imag'] = True
#flag['multi-real'] = True
#flag['multi-real-vp'] = True
flag['uw-real'] = True
#flag['save fig'] = True
#flag['multi-ogive'] = True

binners = dict(ttm=avm.TurbBinner(9600, 32),
               sm=avm.TurbBinner(4800, 16))
pii = 2 * np.pi

if 'bdat' not in vars():
    bdat = {}
    bdat['ttm'] = j14.load('ttm02b-top', 'pax', bin=True)
    bdat['sm'] = smdat.load('SMN-5s', bindat=True)
    for nnm, onm in [('Cspec_vel', 'Cspec_u'),
                     ('Cspec_velrot', 'Cspec_urot'),
                     ('Cspec_velraw', 'Cspec_uraw'),
                     ('Cspec_velacc', 'Cspec_uacc'),
                     ('Cspec_velmot', 'Cspec_umot'), ]:
        if onm in bdat['sm']:
            bdat['sm'].add_data(nnm,
                                bdat['sm'].pop_data(onm),
                                'spec')

pt.twocol()
pairs = [(0, 1), (0, 2), (1, 2)]

do_data = [
    'ttm',
    # 'sm', 
]

for idat, dat_nm in enumerate(do_data):
    bd = bdat[dat_nm]
    if flag.get('amp/phase'):

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
            # U = bd['u'][inds].mean() + 1j * bd['']
            # kaimal = 
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
            z = 10
            fstar = f * z / Uhor
            kaimal = 14 / (1 + 9.6 * fstar) ** (2.4)
            kaimal *= -1 * z * ustar2 / Uhor
            axs[1, icol].plot(f, kaimal, 'm-', lw=5)

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

    if flag.get('uw-real'):

        newfig_kws = dict(figsize=3,
                          sharex=True, sharey=True,
                          right=0.95, left=0.2,
                          bottom=0.15)
        pt_kws = dict(color='0.4', ms=3,
                      zorder=-5, alpha=0.3,
                      rasterized=True, mec='none')

        with pt.style['onecol']():

            m73f = np.array([1e-3, 10])
            m73a = 0.1 * m73f ** (-7. / 3)
            df = 2e-2
            fstar_bins = np.arange(df, 10, df)
            f = bd.freq
            Cs = np.empty_like(fstar_bins)
            Cs2 = np.empty_like(fstar_bins)
            kaimal = 14 / (1 + 9.6 * fstar_bins) ** (2.4)
            inds = np.abs(bd.u) > 1
            #inds = bd.u < -1
            Uhor = np.abs(bd.U[inds])
            fstar = f * z / Uhor[:, None]
            ustar2 = np.sqrt(bd.upwp_[inds]**2 + bd.vpwp_[inds]**2)
            Cdata = (-np.sign(bd['u'][inds][:, None]) *
                     bd['Cspec_vel'][1][inds].real *
                     Uhor[:, None] / (ustar2[:, None] * z))
            Cdata2 = Cdata * fstar
            for idx, fb in enumerate(fstar_bins):
                itmp = (fb - (df / 2) < fstar) & (fstar < fb + (df / 2))
                Cs[idx] = Cdata[itmp].mean()
                Cs2[idx] = Cdata2[itmp].mean()

            fig, ax = pt.newfig(1000 * idat + 203,
                                1, 1, **newfig_kws)

            ax.loglog(fstar_bins, Cs, 'b-', lw=2,
                      label='Avg.')
            ax.loglog(m73f, m73a, 'r--', lw=2,
                      label='$\hat{f}^{-7/3}$')
            ax.loglog(fstar_bins, kaimal, 'k-', lw=2, zorder=-2,
                      label='Kaimal')
            ax.loglog(fstar.flatten(), Cdata.flatten(),
                      '.', label='single', **pt_kws)
            ax.set_xlim([1e-2, 10])
            ax.set_ylim([1e-4, 10])
            # z = 10
            ax.set_ylabel('$\hat{C}\{u,w\}$')
            ax.set_xlabel('$\hat{f}$')
            ax.legend(loc='lower left',
                      prop=dict(size='small'), numpoints=1)

            fig.savefig(pt.figdir + 'CoSpecND_Log01.pdf', dpi=200)

            fig, ax = pt.newfig(1000 * idat + 213,
                                1, 1, **newfig_kws)

            ax.semilogx(fstar_bins, fstar_bins * kaimal, 'k-',
                        label='Kaimal',
                        lw=2, zorder=-2)
            ax.semilogx(fstar_bins, Cs2, 'b-', lw=2,
                        label='Avg.')
            #ax.plot(m73f, m73a, 'r--', lw=2)
            ax.semilogx(fstar.flatten(), fstar.flatten() * Cdata.flatten(),
                        '.', label='single', **pt_kws)
            ax.set_xlim([1e-2, 10])
            ax.set_ylim([-1, 1])
            ax.set_ylabel('$\hat{f}\cdot \hat{C}\{u,w\}$')
            ax.set_xlabel('$\hat{f}$')
            ax.legend(loc='lower left',
                      prop=dict(size='small'), numpoints=1)
            # z = 10
            # kaimal = 14 / (1 + 9.6 * fstar) ** (2.4)
            # kaimal *= -1 * z * ustar2 / Uhor
            # axs[1, icol].plot(f, kaimal, 'm-', lw=5)
            # fig.savefig('NonDim_CoSpecLog01.pdf')
            fig.savefig(pt.figdir + 'CoSpecND_Lin01.pdf', dpi=200)
            
        

            
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

