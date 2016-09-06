import ttm.June2014 as j14
import spec_SM01 as sm
import ptools as pt
import numpy as np
plt = pt.plt

flag = {}

#flag['multi spec'] = True
#flag['epsVu01'] = True
#flag['epsVu02'] = True
flag['epsVu both'] = True
flag['save figs'] = True

if 'dat' not in vars():
    dat = {}
if 'ttm' not in dat:
    dat['ttm'] = j14.load('ttm02b-top', 'pax',
                          bin=True)
if 'sm' not in dat:
    dat['sm'] = sm.dnow


def gt(val):
    return lambda x: x > val


def lt(val):
    return lambda x: x < val


for idat, (datky, dnow) in enumerate(dat.iteritems()):

    if flag.get('epsVu01'):

        with pt.style['onecol']():

            fig, ax = pt.newfig(221 + idat, 1, 1,
                                figsize=3,
                                left=0.2, right=0.95,
                                bottom=0.15,)

            inds = dnow.u > 0
            ax.loglog(np.abs(dnow.U[inds]), dnow.epsilon[inds], 'r.')
            ax.loglog(np.abs(dnow.U[~inds]), dnow.epsilon[~inds], 'k.')
            for fctr in np.logspace(-2, 1, 4):
                ax.loglog(np.array([1e-3, 100]) ** (1. / 3) * fctr, np.array([1e-6, 0.1]), 'k:')

            ax.set_xlim([1e-2, 1e1])
            ax.set_ylim([1e-6, 1e-2])

            ax.set_ylabel(r'$\epsilon\ \mathrm{[W/kg]}$')
            ax.set_xlabel(r'$\bar{U}$')

            if flag.get('save figs'):
                fig.savefig(pt.figdir + 'EpsVU_{}_01.pdf'.format(datky.upper()))

    if flag.get('epsVu02'):

        with pt.style['twocol']():

            fig, axs = pt.newfig(231 + idat, 1, 2,
                                 sharex=True, sharey=True,
                                 figsize=3,
                                 left=0.1, right=0.85,
                                 bottom=0.15, top=0.92)

            u_rngs = np.arange(0.2, 3., 0.2)
            Utmp = np.array([1e-5, 1e5])
            for iax, ax in enumerate(axs):
                inds = [gt(0.6), lt(-0.6)][iax](dnow.u) & ~np.isnan(dnow.epsilon)
                angoff = dnow.principal_angle - np.pi * [1, 0][iax]
                ax.set_title(['Ebb', 'Flood'][iax])
                U = np.abs(dnow.U[inds])
                eps = dnow.epsilon[inds]
                l_ratio = np.exp(np.log(eps / U ** 3).mean())
                clr = -(np.angle(dnow.U[inds]) - angoff) * 180 / np.pi
                crange = [-30, 30]
                # clr = dnow.mpltime[inds] - dnow.mpltime[0]
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
                    # uncvals = (np.nanmean(dnow.epsilon[i2]) +
                    #            np.array([-1, 1]) * np.nanstd(dnow.epsilon[i2]))
                    # uncvals[uncvals < 0] = 1e-12
                    ax.plot(np.mean(ur) * np.array([1, 1]),
                            [unc[0], unc[-1]], '-', lw=1,
                            color=[0, 1, 0], zorder=2)
                    ax.plot(np.mean(ur) * np.array([1, 1]),
                            [unc[0], unc[-1]], '-', lw=1.8,
                            color='k', zorder=1)
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
            if flag.get('save figs'):
                fig.savefig(pt.figdir + 'EpsVU_{}_02.pdf'.format(datky.upper()))


if flag.get('epsVu both'):

    with pt.style['twocol']():

        fig, axs = pt.newfig(231 + idat, 1, 2,
                             sharex=True, sharey=True,
                             figsize=3,
                             left=0.1, right=0.85,
                             bottom=0.15, top=0.92)

        dnow = dict(u=np.concatenate((dat['sm'].u, dat['ttm'].u)),
                    U=np.concatenate((np.abs(dat['sm'].U), np.abs(dat['ttm'].U))),
                    eps=np.concatenate((dat['sm'].epsilon, dat['ttm'].epsilon)),
                    )

        u_rngs = np.arange(0.2, 3., 0.2)
        Utmp = np.array([1e-5, 1e5])
        for iax, ax in enumerate(axs):
            for ky, d in dat.iteritems():
                inds = [gt(0.6), lt(-0.6)][iax](d.u) & ~np.isnan(d.epsilon)
                if ky == 'sm':
                    marker = 'o'
                    label = 'SM'
                elif ky == 'ttm':
                    marker = 'D'
                    label = 'TTM'
                scat = ax.plot(np.abs(d.U[inds]), d.epsilon[inds], marker, c='k',
                               ms=3, label=label)
            ax.set_title(['Ebb', 'Flood'][iax])
            inds = [gt(0.6), lt(-0.6)][iax](dnow['u']) & ~np.isnan(dnow['eps'])
            U = dnow['U'][inds]
            eps = dnow['eps'][inds]
            l_ratio = np.exp(np.log(eps / U ** 3).mean())

            for ur in zip(u_rngs[:-1], u_rngs[1:]):
                i2 = pt.within(U, *ur)
                clr = [0, 1, 0]
                if i2.sum() < 10:
                    continue
                ax.plot(np.mean(ur), eps[i2].mean(), '.',
                        color=clr,
                        ms=7, zorder=2)
                ax.plot(np.mean(ur), eps[i2].mean(), '.',
                        color='k',
                        ms=9, zorder=1)
                unc = pt.boot(eps[i2])
                ax.plot(np.mean(ur) * np.array([1, 1]),
                        [unc[0], unc[-1]], '-', lw=1,
                        color=clr, zorder=2)
                ax.plot(np.mean(ur) * np.array([1, 1]),
                        [unc[0], unc[-1]], '-', lw=1.8,
                        color='k', zorder=1)
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlabel(r'$\bar{U}$')
            ticks = np.arange(0.4, 3, 0.4)
            ax.xaxis.set_ticks(ticks)
            ax.xaxis.set_ticklabels(['{:0.1f}'.format(tk) for tk in ticks])
            ax.xaxis.set_ticks(np.arange(0.2, 3, 0.2), minor=True)
            ax.xaxis.grid(True, 'minor')
            ax.plot(Utmp, l_ratio * (Utmp ** 3), 'b-', zorder=-5)
            # print eps.max()
            ax.annotate(r'$\epsilon = \bar{U}^3 \cdot$ %0.1e m$^{-1}$' % (l_ratio),
                        (0.7, l_ratio * 0.7 ** 3),
                        (0.04, 0.04), textcoords='axes fraction',
                        bbox=dict(boxstyle="round", fc=[0, 0, 1, 0.3], ec='none'),
                        arrowprops=dict(arrowstyle="simple",
                                        fc=[0, 0, 1, 0.3], ec="none",
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
        axs[1].legend(loc='upper left', bbox_to_anchor=[1.02, 1.0],
                      numpoints=1, handlelength=1, handletextpad=0.2)
        if flag.get('save figs'):
            fig.savefig(pt.figdir + 'EpsVU_03.pdf')

