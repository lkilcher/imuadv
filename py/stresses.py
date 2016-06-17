import ttm.June2014 as j14
import ttm.sm2015 as sm15
import ptools as pt
import dolfyn.adv.api as avm
import numpy as np
plt = pt.plt

flag = {}
# flag['amp/phase'] = True
# flag['real/imag'] = True
flag['multi-real'] = True
flag['save fig'] = True

binner = avm.TurbBinner(4800, 16)
pii = 2 * np.pi

if 'rdat' not in vars():
    rdat = {}
    rdat['ttm'] = j14.load('ttm02b-top', 'pax', )
    rdat['sm'] = sm15.load('SM_Nose', 'pax')

pairs = [(0, 1), (0, 2), (1, 2)]

if 'bdat' not in vars():
    bdat = {}
for nm in rdat:
    if nm not in bdat:
        rd = rdat[nm]
        bd = bdat[nm] = binner(rd)

        bd.add_data('Cspec_u', np.empty_like(bd.Spec, dtype=np.complex64), 'spec')
        bd.add_data('Cspec_umot', np.empty_like(bd.Spec, dtype=np.complex64), 'spec')
        bd.add_data('Cspec_uraw', np.empty_like(bd.Spec, dtype=np.complex64), 'spec')
        umot = rd.uacc + rd.urot
        for ip, ipair in enumerate(pairs):
            bd.Cspec_u[ip] = binner.cpsd(rd._u[ipair[0]], rd._u[ipair[1]],
                                         n_fft=binner.n_fft)
            bd.Cspec_umot[ip] = binner.cpsd(umot[ipair[0]], umot[ipair[1]],
                                            n_fft=binner.n_fft)
            bd.Cspec_uraw[ip] = binner.cpsd(rd.uraw[ipair[0]], rd.uraw[ipair[1]],
                                            n_fft=binner.n_fft)

pt.twocol()


for idat, dat_nm in enumerate(['ttm', 'sm']):
    bd = bdat[dat_nm]

    if flag.get('amp/phase'):

        fig, axs = pt.newfig(1000 * idat + 101,
                             3, 2, figsize=5,
                             sharex=True, sharey='col',
                             right=0.8, bottom=0.08)

        velrng = (1, 1.5)

        print('\n{: >4} Stresses:   spec      mean     <spec:uncorr>'.format(dat_nm.upper()))
        inds = (velrng[0] < np.abs(bd['u'])) & np.abs(bd['u'] < velrng[1])
        f = bd.freq
        for irow, axrow in enumerate(axs):
            dtmp = bd.Cspec_u[irow][inds].mean(0) * pii
            axrow[0].loglog(f, np.abs(dtmp), 'b-', label=pt.latex['ue'])
            axrow[1].semilogx(f, np.angle(dtmp), 'b-', label=pt.latex['ue'])
            corr = np.trapz(dtmp, f)
            dtmp = bd.Cspec_umot[irow][inds].mean(0) * pii
            axrow[0].loglog(f, np.abs(dtmp), 'r-', label=pt.latex['uhead'])
            axrow[1].semilogx(f, np.angle(dtmp), 'r-', label=pt.latex['uhead'])
            dtmp = bd.Cspec_uraw[irow][inds].mean(0) * pii
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

        inds = (velrng[0] < np.abs(bd['u'])) & np.abs(bd['u'] < velrng[1])
        f = bd.freq
        for irow, axrow in enumerate(axs):
            dtmp = bd.Cspec_u[irow][inds].mean(0) * pii
            axrow[0].semilogx(f, dtmp.real, 'b-', label=pt.latex['ue'])
            axrow[1].semilogx(f, dtmp.imag, 'b-', label=pt.latex['ue'])
            dtmp = bd.Cspec_umot[irow][inds].mean(0) * pii
            axrow[0].semilogx(f, dtmp.real, 'r-', label=pt.latex['uhead'])
            axrow[1].semilogx(f, dtmp.imag, 'r-', label=pt.latex['uhead'])
            dtmp = bd.Cspec_uraw[irow][inds].mean(0) * pii
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
            inds = (vr[0] < np.abs(bd['u'])) & np.abs(bd['u'] < vr[1])
            axs[0, icol].text(.9, .9, 'N={}'.format(inds.sum()),
                              ha='right', va='top', fontsize='medium',
                              transform=axs[0, icol].transAxes)
            if vr[0] == 0:
                axs[0, icol].set_title(r"$ |\bar{u}| < %0.1f$" % vr[1],
                                       fontsize='medium')
            else:
                axs[0, icol].set_title(r"$%0.1f < |\bar{u}| < %0.1f$" % vr,
                                       fontsize='medium')
            for irow, ax in enumerate(axcol):
                dtmp = bd.Cspec_u[irow][inds].mean(0) * pii
                ax.text(.9, .1,
                        r"%2.2g" % (np.trapz(dtmp, f) * 10000, ),
                        ha='right', va='bottom',
                        transform=ax.transAxes)
                # r"$\overline{%s'%s'}=$%0.0e" % (pt.vel_comps[pairs[irow][0]],
                #                                 pt.vel_comps[pairs[irow][1]],
                #                                 np.trapz(dtmp, f), ),
                ax.semilogx(f, dtmp.real, 'b-', label=pt.latex['ue'])
                dtmp = bd.Cspec_umot[irow][inds].mean(0) * pii
                ax.semilogx(f, dtmp.real, 'r-', label=pt.latex['uhead'], zorder=-2)
                dtmp = bd.Cspec_uraw[irow][inds].mean(0) * pii
                ax.semilogx(f, dtmp.real, 'k-', zorder=-5, label=pt.latex['umeas'])
                ax.axhline(0, color='k', linestyle=':', zorder=-5, lw=1)

        axs[0, -1].legend(loc='upper left', bbox_to_anchor=[1.1, 1])

        axs[0, 0].set_ylim([-0.2, 0.2])
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
            fig.savefig(pt.figdir + 'StressSpec_{}_03.pdf'.format(dat_nm.upper()))
