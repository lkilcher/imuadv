import ttm.June2014 as j14
import ptools as pt
import numpy as np
import dolfyn.adv.api as avm
plt = pt.plt

eps_freqs = np.array([[.3, 1],
                      [.3, 1],
                      [.3, 3], ])

flag = {}

#flag['multi spec'] = True
#flag['turb time01'] = True
#flag['epsVu01'] = True
flag['epsVu02'] = True

pii = 2 * np.pi

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

if 'dat' not in vars():
    dat = j14.load('ttm02b-top', 'pax',
                   bin=True)

    binner = avm.TurbBinner(dat.props['n_fft'], dat.fs)
    epstmp = np.zeros_like(dat.u)
    Ntmp = 0
    noise = vard['Spec']['noise']
    for idx, frq_rng in enumerate(eps_freqs):
        if frq_rng is None:
            continue
        om_rng = frq_rng * pii
        N = ((om_rng[0] < dat.omega) & (dat.omega < om_rng[1])).sum()
        sptmp = dat.Spec[idx] - noise[idx] / pii
        sptmp[sptmp < 0] = 0
        tmp = binner.calc_epsilon_LT83(sptmp,
                                       dat.omega,
                                       np.abs(dat.U),
                                       om_rng)
        epstmp += tmp * N
        Ntmp += N
    epstmp /= Ntmp
    # epstmp[np.abs(dat.U) < 0.2] = np.NaN
    dat.add_data('epsilon', epstmp, 'main')

if flag.get('multi spec'):
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
            # The col-only loop
            axs[irow, 0].set_ylabel('$\mathrm{[m^2s^{-2}/Hz]}$')
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

        fig.savefig(pt.figdir + 'SpecFig02_TTM02B-top.pdf')


if flag.get('turb time01'):

    with pt.style['twocol']():

        lgnd_kws = dict(loc='ul', bbox_to_anchor=(1.02, 1), )

        t = (dat.mpltime - (dat.mpltime[0] // 1)) * 24
        fig, axs = pt.newfig(201, 4, 1,
                             figsize=6,
                             hspace=0.17,
                             left=0.12, right=0.81,
                             bottom=0.08, top=0.97,
                             sharex=True, sharey=False)

        ax = axs[0]
        ax.plot(t, dat.u, 'b-', label=r'$\bar{u}$')
        ax.plot(t, dat.v, 'g-', label=r'$\bar{v}$')
        ax.plot(t, dat.w, 'r-', label=r'$\bar{w}$')
        ax.axhline(0, color='k', linestyle=':')
        ax.set_ylabel('$\mathrm{[m/s]}$')
        ax.legend(**lgnd_kws)

        ax = axs[1]
        ax.semilogy(t, dat.tke, 'k-', label=r"tke")
        ax.semilogy(t, dat.upup_, 'b-', label=r"$\overline{u^2}$")
        ax.semilogy(t, dat.vpvp_, 'r-', label=r"$\overline{v^2}$")
        ax.semilogy(t, dat.wpwp_, 'g-', label=r"$\overline{w^2}$")
        ax.set_ylabel(r'$\mathrm{tke}\ \mathrm{[m^2/s^2]}$')
        ax.legend(**lgnd_kws)

        ax = axs[2]
        ax.plot(t, dat.upvp_, 'b-', label=r"$\overline{uv}$")
        ax.plot(t, dat.upwp_, 'r-', label=r"$\overline{uw}$")
        ax.plot(t, dat.vpwp_, 'g-', label=r"$\overline{vw}$")
        ax.set_ylabel(r'$\mathrm{stress}\ \mathrm{[m^2/s^2]}$')
        ax.legend(**lgnd_kws)

        ax = axs[3]
        ax.semilogy(t, dat.epsilon, 'k.')
        ax.set_ylabel('$\epsilon\ \mathrm{[W/kg]}$')

        ax.set_xlim([t[0], t[-1]])
        ax.set_xlabel('Time [Hours since 00:00 June 18, 2014]')

        for iax, ax in enumerate(axs):
            ax.text(0.02, 0.93, '({})'.format('ABCDEFG'[iax]),
                    transform=ax.transAxes,
                    va='top', ha='left')
            ax.fill_between(t, 1, 0, where=dat.u > 0.2,
                            edgecolor='none', facecolor='0.9',
                            zorder=-10,
                            transform=ax.get_xaxis_transform())
            ax.fill_between(t, 1, 0, where=dat.u < -0.2,
                            edgecolor='none', facecolor='0.95',
                            zorder=-10,
                            transform=ax.get_xaxis_transform())

        axs[2].set_ylim([-0.02, 0.08])
        axs[3].set_ylim([1e-6, 1e-2])

        fig.savefig(pt.figdir + 'TurbTime_TTM_01.pdf')

if flag.get('epsVu01'):

    with pt.style['onecol']():

        fig, ax = pt.newfig(202, 1, 1,
                            figsize=3,
                            left=0.2, right=0.95,
                            bottom=0.15,)

        inds = dat.u > 0
        ax.loglog(np.abs(dat.U[inds]), dat.epsilon[inds], 'r.')
        ax.loglog(np.abs(dat.U[~inds]), dat.epsilon[~inds], 'k.')
        for fctr in np.logspace(-2, 1, 4):
            ax.loglog(np.array([1e-3, 100]) ** (1. / 3) * fctr, np.array([1e-6, 0.1]), 'k:')

        ax.set_xlim([1e-2, 1e1])
        ax.set_ylim([1e-6, 1e-2])

        ax.set_ylabel(r'$\epsilon\ \mathrm{[W/kg]}$')
        ax.set_xlabel(r'$\bar{U}$')

        fig.savefig(pt.figdir + 'EpsVU_TTM_01.pdf')

if flag.get('epsVu02'):

    with pt.style['twocol']():

        fig, axs = pt.newfig(203, 1, 2,
                             sharex=True, sharey=True,
                             figsize=3,
                             left=0.1, right=0.85,
                             bottom=0.15, top=0.92)

        u_rngs = np.arange(0.2, 3., 0.2)
        Utmp = np.array([1e-5, 1e5])
        for iax, ax in enumerate(axs):
            if iax == 0:
                inds = (dat.u > 0.6) & ~np.isnan(dat.epsilon)
                ax.set_title('Ebb')
                angoff = dat.principal_angle - np.pi
            else:
                inds = (dat.u < -0.6) & ~np.isnan(dat.epsilon)
                ax.set_title('Flood')
                angoff = dat.principal_angle
            U = np.abs(dat.U[inds])
            eps = dat.epsilon[inds]
            l_ratio = np.exp(np.log(eps / U ** 3).mean())
            clr = -(np.angle(dat.U[inds]) - angoff) * 180 / np.pi
            crange = [-30, 30]
            # clr = dat.mpltime[inds] - dat.mpltime[0]
            # crange = [0, 1]
            scat = ax.scatter(U, eps, c=clr,
                              s=8, marker='o', linewidths=0,
                              cmap='coolwarm',
                              vmin=crange[0], vmax=crange[1])
            cax = fig.add_axes([0.88, 0.15, 0.02, 0.77])
            cbar = plt.colorbar(scat, cax=cax)
            cbar.set_label(r'$\Delta\theta\ \mathrm{[degrees]}$')
            cbar.set_ticks(np.arange(-30, 31, 10))
            # ax.scatter(U, eps, c='0.7',
            #            s=8, marker='o', linewidths=0, )
            for ur in zip(u_rngs[:-1], u_rngs[1:]):
                i2 = pt.within(U, *ur)
                if i2.sum() <= 5:
                    continue
                ax.plot(np.mean(ur), eps[i2].mean(), '.',
                        color=[0, 1, 0],
                        ms=7, zorder=2)
                ax.plot(np.mean(ur), eps[i2].mean(), '.',
                        color='k',
                        ms=9, zorder=1)
                unc = pt.boot(eps[i2])
                # uncvals = (np.nanmean(dat.epsilon[i2]) +
                #            np.array([-1, 1]) * np.nanstd(dat.epsilon[i2]))
                # uncvals[uncvals < 0] = 1e-12
                ax.plot(np.mean(ur) * np.array([1, 1]),
                        [unc[0], unc[-1]], '-', lw=1,
                        color=[0, 1, 0], zorder=2)
                ax.plot(np.mean(ur) * np.array([1, 1]), 
                        [unc[0], unc[-1]], '-', lw=1.8,
                        color='k', zorder= 1)
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlabel(r'$\bar{U}$')
            ticks = np.arange(0.4, 3, 0.4)
            ax.xaxis.set_ticks(ticks)
            ax.xaxis.set_ticklabels(['{:0.1f}'.format(tk) for tk in ticks])
            ax.xaxis.set_ticks(np.arange(0.2, 3, 0.2), minor=True)
            ax.xaxis.grid(True, 'minor')
            ax.plot(Utmp, l_ratio * (Utmp ** 3), 'k-', zorder=-5)
            # print eps.max()
            ax.annotate(r'$\epsilon = \bar{U}^3 \cdot$ %0.1e m$^{-1}$' % (l_ratio),
                        (0.7, l_ratio * 0.7 ** 3),
                        (0.04, 0.04), textcoords='axes fraction',
                        bbox=dict(boxstyle="round", fc='0.87', ec='none'),
                        arrowprops=dict(arrowstyle="simple",
                                        fc='0.8', ec="none",
                                        connectionstyle="arc3,rad=-0.3",
                        ),
            )
            # ax.text(0.96, 0.04,
            #         r'$\epsilon = \bar{U}^3 \cdot$ %0.1e m$^{-1}$' % (l_ratio),
            #         color='r',
            #         transform=ax.transAxes, ha='right', va='bottom', )
        # axs[0].plot(Utmp, 6e-5 * (Utmp ** 3), 'r-')
        # axs[1].plot(Utmp, 2e-5 * (Utmp ** 3), 'r-')
        ax.set_xlim([0.6, 2.4])
        ax.set_ylim([1e-6, 3e-3])
        axs[0].set_ylabel(r'$\epsilon\ \mathrm{[W/kg]}$')
        fig.savefig(pt.figdir + 'EpsVU_TTM_02.pdf')




