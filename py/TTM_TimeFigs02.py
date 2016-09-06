import ptools as pt
import ttm.June2014 as j14
import matplotlib.pyplot as plt
import numpy as np
import gis

flg = {}
flg['turb time01'] = True
#flg['turb time02'] = True
flg['epsVprod01'] = True
#flg['epsVprod02'] = True
#flg['save figs'] = True

lgnd_kws = dict(loc='upper left', bbox_to_anchor=(1.02, 1), )


if 'dat' not in vars():
    dat = j14.load('ttm02b-top', 'pax',
                   bin=True)
    dat2 = j14.load('ttm02b-bottom', 'pax',
                    bin=True)
    dat3 = j14.load('ttm01b-top', 'pax',
                    bin=True).subset(slice(3, -2))
    dat4 = j14.load('ttm01b-bottom', 'pax',
                    bin=True).subset(slice(3, -2))

if flg.get('turb time01'):

    with pt.style['twocol']():


        t = (dat.mpltime - (dat.mpltime[0] // 1)) * 24
        fig, axs = pt.newfig(201, 4, 1,
                             figsize=6,
                             hspace=0.17,
                             left=0.12, right=0.81,
                             bottom=0.08, top=0.97,
                             sharex=True, sharey=False)

        ax = axs[0]
        ax.plot(t, np.sqrt(dat.u ** 2 + dat.v ** 2), 'k-', label=r'$\bar{U}$')
        ax.plot(t, dat.u, 'b-', label=r'$\bar{u}$')
        ax.plot(t, dat.v, 'g-', label=r'$\bar{v}$')
        ax.plot(t, dat.w, 'r-', label=r'$\bar{w}$')
        ax.axhline(0, color='k', linestyle=':')
        ax.set_ylabel('$\mathrm{[m/s]}$')
        ax.text(12.7, .97, 'ebb 1',
                va='top', ha='center',
                transform=ax.get_xaxis_transform())
        ax.text(18.9, .97, 'flood',
                va='top', ha='center',
                transform=ax.get_xaxis_transform())
        ax.text(26.2, 0.97, 'ebb 2',
                va='top', ha='center',
                transform=ax.get_xaxis_transform())
        ax.legend(**lgnd_kws)

        ax = axs[1]
        ax.semilogy(t, dat.tke, 'k-', label=r"tke")
        ax.semilogy(t, dat.upup_, 'b-', label=r"$\overline{u^2}$")
        ax.semilogy(t, dat.vpvp_, 'g-', label=r"$\overline{v^2}$")
        ax.semilogy(t, dat.wpwp_, 'r-', label=r"$\overline{w^2}$")
        ax.set_ylabel(r'$\mathrm{[m^2/s^2]}$')
        ax.legend(**lgnd_kws)

        ax = axs[2]
        ax.plot(t, dat.upvp_, '-', c=[0, 0.8, 0.8], label=r"$\overline{uv}$")
        ax.plot(t, dat.upwp_, '-', c=[0.8, 0, 0.8], label=r"$\overline{uw}$")
        ax.plot(t, dat.vpwp_, '-', c=[0.8, 0.8, 0], label=r"$\overline{vw}$")
        ax.set_ylabel(r'$\mathrm{[m^2/s^2]}$')
        ax.legend(**lgnd_kws)

        ax = axs[3]
        ax.semilogy(t, dat.epsilon, 'ko', ms=3, zorder=4,
                    mfc='k', mec='none',
                    label='$\epsilon$')
        ax.set_ylabel('$\mathrm{[W/kg]}$')

        prod = -(dat.u - dat2.u) / 0.5 * dat.upwp_
        inds = prod > 0
        ax.semilogy(t[inds], prod[inds], 'ro', ms=3, zorder=5,
                    mfc='r', mec='none', alpha=.6,
                    label='$P_{uz}$')
        ax.semilogy(t[~inds], -prod[~inds], 'ro', ms=2.5, zorder=6,
                    mfc='none', mec='r', alpha=0.6,
                    label='$-P_{uz}$')

        ax.set_ylim([1e-6, 1e-2])
        #ax.set_ylabel('$\partial u/\partial z\ \mathrm{[s^{-1}]}$')
        lkw = lgnd_kws.copy()
        lkw['numpoints'] = 1
        lkw['handletextpad'] = 0.2
        ax.legend(**lkw)

        ax.set_xlim([t[0], t[-1]])
        ax.set_xlabel('Time [Hours since 00:00 June 18, 2014]')

        for iax, ax in enumerate(axs):
            ax.text(0.02, 0.93, '({})'.format('ABCDEFG'[iax]),
                    transform=ax.transAxes,
                    va='top', ha='left')
            ax.fill_between(t, 1, 0, where=dat.u > 1.0,
                            edgecolor='none', facecolor='0.9',
                            zorder=-10,
                            transform=ax.get_xaxis_transform())
            ax.fill_between(t, 1, 0, where=dat.u < -1.0,
                            edgecolor='none', facecolor='0.95',
                            zorder=-10,
                            transform=ax.get_xaxis_transform())

        axs[2].set_ylim([-0.02, 0.08])
        axs[3].set_ylim([1e-6, 1e-2])

        if flg.get('save figs'):
            fig.savefig(pt.figdir + 'TurbTime_TTM_01.pdf')


if flg.get('turb time02'):

    with pt.style['twocol']():

        t = (dat.mpltime - (dat.mpltime[0] // 1)) * 24
        fig, axs = pt.newfig(202, 4, 1,
                             figsize=6,
                             hspace=0.17,
                             left=0.12, right=0.81,
                             bottom=0.08, top=0.97,
                             sharex=True, sharey=False)

        ax = axs[0]
        ax.semilogy(t, dat.tke, 'k-', label=r"tke")
        ax.semilogy(t, dat.upup_, 'b-', label=r"$\overline{u^2}$")
        ax.semilogy(t, dat.vpvp_, 'r-', label=r"$\overline{v^2}$")
        ax.semilogy(t, dat.wpwp_, 'g-', label=r"$\overline{w^2}$")
        ax.set_ylabel(r'$\mathrm{tke}\ \mathrm{[m^2/s^2]}$')
        ax.legend(**lgnd_kws)

        ax = axs[1]
        ax.plot(t, dat.upvp_, 'b-', label=r"$\overline{uv}$")
        ax.plot(t, dat.upwp_, 'r-', label=r"$\overline{uw}$")
        ax.plot(t, dat.vpwp_, 'g-', label=r"$\overline{vw}$")
        ax.plot(t, dat2.upvp_, 'b.', label=r"$\overline{uv}$")
        ax.plot(t, dat2.upwp_, 'r.', label=r"$\overline{uw}$")
        ax.plot(t, dat2.vpwp_, 'g.', label=r"$\overline{vw}$")
        ax.set_ylabel(r'$\mathrm{stress}\ \mathrm{[m^2/s^2]}$')
        #ax.legend(**lgnd_kws)

        ax = axs[2]
        ax.semilogy(t, dat.epsilon, 'k.')
        ax.set_ylabel('$\epsilon\ \mathrm{[W/kg]}$')

        ax = axs[3]
        prod = -(dat.u - dat2.u) / 0.5 * (dat.upwp_ + dat2.upwp_) / 2.
        prod2 = -(dat.u - dat2.u) / 0.5 * (dat.upwp_) / 2.
        inds = prod > 0
        ax.semilogy(t[inds], prod[inds], 'k.')
        ax.semilogy(t[~inds], -prod[~inds], 'r.')
        ax.set_ylim([1e-6, 1e-2])
        #ax.set_ylabel('$\partial u/\partial z\ \mathrm{[s^{-1}]}$')

if flg.get('epsVprod01'):

    with pt.style['onecol']():

        t = (dat.mpltime - (dat.mpltime[0] // 1)) * 24
        fig, axs = pt.newfig(301, 1, 1,
                             figsize=3,
                             hspace=0.17,
                             #left=0.12, right=0.81,
                             #bottom=0.08, top=0.97,
                             sharex=True, sharey=False)

        dnow0, dnow1 = dat3, dat4
        dnow0, dnow1 = dat, dat2
        prod = -(dnow0.u - dnow1.u) / 0.5 * (dnow0.upwp_ + dnow1.upwp_) / 2.
        eps = (dnow0.epsilon + dnow1.epsilon) / 2
        umag = np.abs(dnow0.u)
        iumag = np.abs(dnow0.u) > 1.0
        inds = (prod > 0)

        logbins = np.logspace(-6, -2, 13)

        axs.loglog(eps[inds & iumag], prod[inds & iumag], 'o',
                   label='$P_{uz}>0\ \mathrm{(N = %d)}$' % (inds & iumag).sum(),
                   mfc='k', mec='none', ms=3,
                   zorder=2)
        axs.loglog(eps[~inds & iumag], -prod[~inds & iumag], 'o',
                   label='$P_{uz}<0\ \mathrm{(N = %d)}$' % (~inds & iumag).sum(),
                   mfc='none', mec='k', ms=2.5,
                   zorder=-2)
        axs.legend(loc='upper left', numpoints=1, handlelength=1,
                   handletextpad=0.2,
                   prop=dict(size='medium'))

        # for rng in zip(logbins[:-1], logbins[1:]):
        #     inds = (rng[0] < eps) & (eps < rng[1])
        #     xval = np.mean(rng)
        #     plt.plot(xval, prod[inds].mean(), 'r.')
        #     plt.plot(xval, -prod[inds].mean(), 'ro')

        axs.plot([1e-8, 1], [1e-8, 1], 'k:')
        axs.set_ylim([1e-6, 1e-2])
        axs.set_xlim([1e-6, 1e-2])
        axs.set_xlabel('$\epsilon\ \mathrm{[W/kg]}$')
        axs.set_ylabel('$P_{uz}\ \mathrm{[W/kg]}$')

    if flg.get('save figs'):
        fig.savefig(pt.figdir + 'EpsVProd01.pdf')

if flg.get('epsVprod02'):

    with pt.style['onecol']():

        t = (dat.mpltime - (dat.mpltime[0] // 1)) * 24
        fig, axs = pt.newfig(302, 1, 1,
                             figsize=3,
                             hspace=0.17,
                             sharex=True, sharey=False)

        dX = (gis.diffll2(dat.props['latlon'][::-1],
                          dat3.props['latlon'][::-1])[0][0][0] *
              np.exp(-1j * dat.props['principal_angle']))
        dy = dX.imag  # 32.4 m
        dudy = np.zeros_like(dat.u)
        dudy += (dat4.u - dat2.u) / dy
        dudy += (dat3.u - dat.u) / dy
        dudy *= 0.5

        # dat.u is the top
        dudz = (dat.u - dat2.u) / 0.5

        dwdx = (dat.w[2:] - dat.w[:-2]) / (dat.u[1:-1] * 297.6)
        dudx = (dat.u[2:] - dat.u[:-2]) / (dat.u[1:-1] * 297.6)

        prod_uy = -dudy * (dat.upvp_ + dat2.upvp_ + dat3.upvp_ + dat4.upvp_) / 4.

        # prod_wz = -(dat3.u - dat4.u) / 0.5 * (dat3.upwp_ + dat4.upwp_) / 2.
        # prod_wz += -(dat.u - dat2.u) / 0.5 * (dat.upwp_ + dat2.upwp_) / 2.
        # prod_wz *= 0.5
        prod_wz = -dudz * (dat.upwp_ + dat2.upwp_) / 2.

        # prod = prod_uy + prod_wz
        prod = prod_wz
        iumag = np.abs(dat.u) > 1.0
        #iumag &= (dudz * dat.u) > 0
        inds = (prod > 0)
        axs.loglog(dat.epsilon[~inds & iumag], -prod[~inds & iumag], '.', color='0.6')
        axs.loglog(dat.epsilon[inds & iumag], prod[inds & iumag], 'k.')
        axs.plot([1e-8, 1], [1e-8, 1], 'k-')
        axs.set_ylim([1e-6, 1e-2])
        axs.set_xlim([1e-6, 1e-2])

        print 'good: {}, bad: {} - ({:.1%})'.format(sum(inds & iumag),
                                                    sum(~inds & iumag),
                                                    float(sum(inds & iumag)) / sum(iumag))
