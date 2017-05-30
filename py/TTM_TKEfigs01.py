import ptools as pt
import data.ttm_june2014 as j14
import numpy as np

flg = {}
flg['turb time01'] = True
flg['epsVprod01'] = True
flg['save figs'] = True

lgnd_kws = dict(loc='upper left', bbox_to_anchor=(1.02, 1), )


if 'dat' not in vars():
    dat = j14.load('ttm02b-top', 'pax',
                   bin=True)
    dat2 = j14.load('ttm02b-bot', 'pax',
                    bin=True)

p_angle = (90 - 310) * np.pi / 180
#p_angle = dat.props['principal_angle'] + 16 / 180 * np.pi

#dnow = dat3
#dnow0, dnow1 = dat3, dat4
#dnow2, dnow3 = dat, dat2
dnow = dat
dnow0, dnow1 = dat, dat2

dz = 0.5

tke_eqn = {}
# ###
# Dissipation
tke_eqn['eps'] = (dnow0.epsilon + dnow1.epsilon) / 2

# ###
# Production terms
dudz = (dnow0.u - dnow1.u) / 0.5
tke_eqn['Prod_uz'] = (-dudz * (dnow0.upwp_ + dnow1.upwp_) / 2.)
#tke_eqn['Prod_uz'] = -(dnow.u - dat2.u) / dz * dnow.upwp_
# prod = -(dnow.u - dat2.u) / dz * (dnow.upwp_ + dat2.upwp_) / 2.
# prod2 = -(dnow.u - dat2.u) / dz * (dnow.upwp_) / 2.
# prod = -(dnow.u - dat2.u) / dz * dnow.upwp_

tke_eqn['dt'] = np.pad(np.diff(dnow0.tke + dnow1.tke) / 2,
                       [1, 0], 'edge') / 279.6

# ###
# Advective terms
t1 = (dnow0.tke + dnow1.tke) * 0.5

# ###
# Transport terms
tprd1 = (dnow0.tripprod[:, 1] + dnow1.tripprod[:, 1]).sum(0) / 2

Tw = (dnow0.tripprod[:, 2] - dnow1.tripprod[:, 2]).sum(0) / dz
# Tw += (dnow2.tripprod[:, 2] - dnow3.tripprod[:, 2]).sum(0) / dz
# Tw /= 2
tke_eqn['Trans_z'] = Tw

t_hours = (dnow.mpltime - (dnow.mpltime[0] // 1)) * 24

if flg.get('turb time01'):

    with pt.style['twocol']():

        fig, axs = pt.newfig(201, 4, 1,
                             figsize=6,
                             hspace=0.17,
                             left=0.12, right=0.81,
                             bottom=0.08, top=0.97,
                             sharex=True, sharey=False)

        ax = axs[0]
        ax.plot(t_hours,
                np.sqrt(dnow.u ** 2 + dnow.v ** 2),
                'k-', label=r'$\bar{U}$')
        ax.plot(t_hours, dnow.u, 'b-', label=r'$\bar{u}$')
        ax.plot(t_hours, dnow.v, 'g-', label=r'$\bar{v}$')
        ax.plot(t_hours, dnow.w, 'r-', label=r'$\bar{w}$')
        ax.axhline(0, color='k', linestyle=':')
        ax.set_ylabel('$\mathrm{[m\,s^{-1}]}$')
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
        ax.semilogy(t_hours, dnow.tke, 'k-', label=r"tke")
        ax.semilogy(t_hours, dnow.upup_, 'b-', label=r"$\overline{u^2}$")
        ax.semilogy(t_hours, dnow.vpvp_, 'g-', label=r"$\overline{v^2}$")
        ax.semilogy(t_hours, dnow.wpwp_, 'r-', label=r"$\overline{w^2}$")
        ax.set_ylabel(r'$\mathrm{[m^2\,s^{-2}]}$')
        ax.legend(**lgnd_kws)

        ax = axs[2]
        ax.plot(t_hours, dnow.upvp_,
                '-', c=[0, 0.8, 0.8], label=r"$\overline{uv}$")
        ax.plot(t_hours, dnow.upwp_,
                '-', c=[0.8, 0, 0.8], label=r"$\overline{uw}$")
        ax.plot(t_hours, dnow.vpwp_,
                '-', c=[0.8, 0.8, 0], label=r"$\overline{vw}$")
        ax.set_ylabel(r'$\mathrm{[m^2\,s^{-2}]}$')
        ax.legend(**lgnd_kws)

        ax = axs[3]
        ax.semilogy(t_hours, tke_eqn['eps'], 'ko', ms=3, zorder=4,
                    mfc='k', mec='none',
                    label='$\epsilon$')
        ax.set_ylabel('$\mathrm{[W\,kg^{-1}]}$')

        inds = tke_eqn['Prod_uz'] > 0
        ax.semilogy(t_hours[inds],
                    tke_eqn['Prod_uz'][inds],
                    'ro', ms=3, zorder=5,
                    mfc='r', mec='none', alpha=.6,
                    label='$P_{uz}$')
        ax.semilogy(t_hours[~inds],
                    -tke_eqn['Prod_uz'][~inds],
                    'ro', ms=2.5, zorder=6,
                    mfc='none', mec='r', alpha=0.6,
                    label='$-P_{uz}$')

        ax.set_ylim([1e-6, 1e-2])
        #ax.set_ylabel('$\partial u/\partial z\ \mathrm{[s^{-1}]}$')
        lkw = lgnd_kws.copy()
        lkw['numpoints'] = 1
        lkw['handletextpad'] = 0.2
        ax.legend(**lkw)

        for iax, ax in enumerate(axs):
            ax.text(0.02, 0.93, '({})'.format('ABCDEFG'[iax]),
                    transform=ax.transAxes,
                    va='top', ha='left')
            ax.fill_between(t_hours, 1, 0, where=dnow.u > 1.0,
                            edgecolor='none', facecolor='0.8',
                            zorder=-10,
                            transform=ax.get_xaxis_transform())
            ax.fill_between(t_hours, 1, 0, where=dnow.u < -1.0,
                            edgecolor='none', facecolor='0.9',
                            zorder=-10,
                            transform=ax.get_xaxis_transform())

        axs[2].set_ylim([-0.02, 0.08])
        axs[3].set_ylim([1e-6, 1e-2])
        axs[-1].set_xlim([t_hours[0], t_hours[-1]])
        axs[-1].set_xlabel('Time [Hours since 00:00 June 18, 2014]')

        if flg.get('save figs'):
            fig.savefig(pt.figdir + 'TurbTime_TTM_01.pdf')


if flg.get('epsVprod01'):

    with pt.style['onecol']():

        fig, axs = pt.newfig(301, 1, 1,
                             figsize=3,
                             hspace=0.17,
                             #left=0.12, right=0.81,
                             #bottom=0.08, top=0.97,
                             sharex=True, sharey=False)

        umag = np.abs(dnow0.u)
        iumag = np.abs(dnow0.u) > 1.0
        inds = (tke_eqn['Prod_uz'] > 0)

        alpha = 0.7
        inow = inds & iumag
        axs.loglog(tke_eqn['eps'][inow],
                   tke_eqn['Prod_uz'][inow],
                   'o',
                   label='$P_{uz}>0\ \mathrm{(N = %d)}$' %
                   inow.sum(),
                   mfc='k', mec='none', ms=3,
                   alpha=alpha,
                   zorder=2,)
        inow = ~inds & iumag
        axs.loglog(tke_eqn['eps'][inow],
                   -tke_eqn['Prod_uz'][inow],
                   'o',
                   label='$P_{uz}<0\ \mathrm{(N = %d)}$' %
                   inow.sum(),
                   mfc='none', mec='k', ms=2.5,
                   alpha=alpha,
                   zorder=-2)
        axs.legend(loc='upper left', numpoints=1, handlelength=1,
                   handletextpad=0.2,
                   prop=dict(size='medium'))

        axs.plot([1e-8, 1], [1e-8, 1], 'k:')
        axs.set_ylim([1e-6, 1e-2])
        axs.set_xlim([1e-6, 1e-2])
        axs.set_xlabel('$\epsilon\ \mathrm{[W\,kg^{-1}]}$')
        axs.set_ylabel('$P_{uz}\ \mathrm{[W\,kg^{-1}]}$')

    if flg.get('save figs'):
        fig.savefig(pt.figdir + 'EpsVProd01.pdf')


if flg.get('epsVprod01a'):

    with pt.style['onecol']():

        fig, axs = pt.newfig(301, 1, 1,
                             figsize=3,
                             hspace=0.17,
                             #left=0.12, right=0.81,
                             #bottom=0.08, top=0.97,
                             sharex=True, sharey=False)

        umag = np.abs(dnow0.u)
        iumag = np.abs(dnow0.u) > 1.0
        inds = (tke_eqn['Prod_uz'] > 0)
        iupos = dnow0.u > 0

        alpha = 0.7
        axs.loglog(tke_eqn['eps'][inds & iumag & ~iupos],
                   tke_eqn['Prod_uz'][inds & iumag & ~iupos],
                   'o',
                   label='$P_{uz}>0, u < -1\ \mathrm{(N = %d)}$' %
                   (inds & iumag & ~iupos).sum(),
                   mfc='r', mec='none', ms=3,
                   alpha=alpha,
                   zorder=2,)
        axs.loglog(tke_eqn['eps'][~inds & iumag & ~iupos],
                   -tke_eqn['Prod_uz'][~inds & iumag & ~iupos],
                   'o',
                   label='$P_{uz}<0, u < -1\ \mathrm{(N = %d)}$' %
                   (~inds & iumag & ~iupos).sum(),
                   mfc='none', mec='r', ms=2.5,
                   alpha=alpha,
                   zorder=-2)
        axs.loglog(tke_eqn['eps'][inds & iumag & iupos],
                   tke_eqn['Prod_uz'][inds & iumag & iupos],
                   'o',
                   label='$P_{uz}>0\ \mathrm{(N = %d)}$' %
                   (inds & iumag & iupos).sum(),
                   mfc='b', mec='none', ms=3,
                   alpha=alpha,
                   zorder=2)
        axs.loglog(tke_eqn['eps'][~inds & iumag & iupos],
                   -tke_eqn['Prod_uz'][~inds & iumag & iupos],
                   'o',
                   label='$P_{uz}<0\ \mathrm{(N = %d)}$' %
                   (~inds & iumag & iupos).sum(),
                   mfc='none', mec='b', ms=2.5,
                   alpha=alpha,
                   zorder=-2)

        axs.plot([1e-8, 1], [1e-8, 1], 'k:')
        axs.set_ylim([1e-6, 1e-2])
        axs.set_xlim([1e-6, 1e-2])
        axs.set_xlabel('$\epsilon\ \mathrm{[W\,kg^{-1}]}$')
        axs.set_ylabel('$P_{uz}\ \mathrm{[W\,kg^{-1}]}$')

    if flg.get('save figs'):
        fig.savefig(pt.figdir + 'EpsVProd01a.pdf')

if flg.get('epsVprod02'):

    with pt.style['onecol']():

        fig, axs = pt.newfig(302, 1, 1,
                             figsize=3,
                             hspace=0.17,
                             sharex=True, sharey=False)

        iumag = np.abs(dnow.u) > 1.0
        #iumag &= (dudz * dnow.u) > 0
        inds = (tke_eqn['Prod_uz'] > 0)
        axs.loglog(tke_eqn['eps'][~inds & iumag],
                   -tke_eqn['Prod_uz'][~inds & iumag], '.', color='0.6')
        axs.loglog(tke_eqn['eps'][inds & iumag],
                   tke_eqn['Prod_uz'][inds & iumag], 'k.')
        axs.plot([1e-8, 1], [1e-8, 1], 'k-')
        axs.set_ylim([1e-6, 1e-2])
        axs.set_xlim([1e-6, 1e-2])

        print 'good: {}, bad: {} - ({:.1%})'.format(
            sum(inds & iumag),
            sum(~inds & iumag),
            float(sum(inds & iumag)) / sum(iumag))

if flg.get('epsVprod03'):

    with pt.style['onecol']():

        fig, axs = pt.newfig(303, 2, 1,
                             figsize=5,
                             hspace=0.17,
                             left=0.2, right=0.9,
                             bottom=0.1, top=0.95,
                             sharex=True, sharey=True)

        umag = np.abs(dnow.u)
        iumag = np.abs(dnow.u) > 1.0
        inds = (tke_eqn['Prod_uz'] > 0)
        iupos = dnow.u > 0

        logbins = np.logspace(-6, -2, 13)

        for iax, ax in enumerate(axs):
            if iax == 0:
                inow = iumag & iupos
                ax.set_title('Ebb')
            else:
                ax.set_title('Flood')
                inow = iumag & ~iupos
            ax.loglog(tke_eqn['eps'][inds & inow],
                      tke_eqn['Prod_uz'][inds & inow], 'o',
                      label='$P_{uz}>0,\mathrm{(N = %d)}$' %
                      (inds & inow).sum(),
                      mfc='k', mec='none', ms=3,
                      zorder=2,)
            ax.loglog(tke_eqn['eps'][~inds & inow],
                      -tke_eqn['Prod_uz'][~inds & inow],
                      'o',
                      label='$P_{uz}<0,\mathrm{(N = %d)}$' %
                      (~inds & inow).sum(),
                      mfc='none', mec='k', ms=2.5,
                      zorder=-2)
            ax.plot([1e-8, 1], [1e-8, 1], 'k:')
            ax.set_ylabel('$P_{uz}\ \mathrm{[W\,kg^{-1}]}$')
            ax.legend(loc='upper left', numpoints=1, handlelength=1,
                      handletextpad=0.2,
                      prop=dict(size='medium'))

        # for rng in zip(logbins[:-1], logbins[1:]):
        #     inds = (rng[0] < tke_eqn['eps']) & (tke_eqn['eps'] < rng[1])
        #     xval = np.mean(rng)

        axs[-1].set_xlabel('$\epsilon\ \mathrm{[W\,kg^{-1}]}$')
        axs[0].set_ylim([1e-6, 1e-2])
        axs[0].set_xlim([1e-6, 1e-2])

    if flg.get('save figs'):
        fig.savefig(pt.figdir + 'EpsVProd03.png', dpi=300)
