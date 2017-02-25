import ptools as pt
import ttm.June2014 as j14
import matplotlib.pyplot as plt
import numpy as np
import gis

flg = {}
#flg['turb time01'] = True
#flg['turb time02'] = True
#flg['turb time03'] = True
#flg['epsVprod01'] = True
#flg['epsVprod02'] = True
flg['epsVprod03'] = True
#flg['eqnVeps01'] = True
#flg['eqnVeps02'] = True
flg['save figs'] = True

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

p_angle = (90 - 310) * np.pi / 180
#p_angle = dat.props['principal_angle'] + 16 / 180 * np.pi

#dnow = dat3
#dnow0, dnow1 = dat3, dat4
#dnow2, dnow3 = dat, dat2
dnow = dat
dnow0, dnow1 = dat, dat2
dnow2, dnow3 = dat3, dat4

dX = -(gis.diffll2(dnow0.props['latlon'][::-1],
                   dnow2.props['latlon'][::-1])[0][0][0] *
       np.exp(-1j * p_angle))
dx = dX.real  # -11.0 m
dy = dX.imag  # -40.7 m
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

dudy = np.zeros_like(dnow.u)
dudy += (dnow3.u - dnow1.u) / dy
dudy += (dnow2.u - dnow0.u) / dy
dudy *= 0.5
tke_eqn['Prod_uy'] = (-dudy * (dnow.upvp_ +
                               dnow1.upvp_ +
                               dnow2.upvp_ +
                               dnow3.upvp_) / 4.)

# ###
# Advective terms
t1 = (dnow0.tke + dnow1.tke) * 0.5
t2 = (dnow2.tke + dnow3.tke) * 0.5
v = (dnow.v + dnow1.v + dnow2.v + dnow3.v) * 0.25
tke_eqn['Adv_v'] = v * (t2 - t1) / dy

# ###
# Transport terms
tprd1 = (dnow0.tripprod[:, 1] + dnow1.tripprod[:, 1]).sum(0) / 2
tprd2 = (dnow2.tripprod[:, 1] + dnow3.tripprod[:, 1]).sum(0) / 2
tke_eqn['Trans_y'] = (tprd2 - tprd1) / dy

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
        ax.semilogy(t_hours, dnow.tke, 'k-', label=r"tke")
        ax.semilogy(t_hours, dnow.upup_, 'b-', label=r"$\overline{u^2}$")
        ax.semilogy(t_hours, dnow.vpvp_, 'g-', label=r"$\overline{v^2}$")
        ax.semilogy(t_hours, dnow.wpwp_, 'r-', label=r"$\overline{w^2}$")
        ax.set_ylabel(r'$\mathrm{[m^2/s^2]}$')
        ax.legend(**lgnd_kws)

        ax = axs[2]
        ax.plot(t_hours, dnow.upvp_,
                '-', c=[0, 0.8, 0.8], label=r"$\overline{uv}$")
        ax.plot(t_hours, dnow.upwp_,
                '-', c=[0.8, 0, 0.8], label=r"$\overline{uw}$")
        ax.plot(t_hours, dnow.vpwp_,
                '-', c=[0.8, 0.8, 0], label=r"$\overline{vw}$")
        ax.set_ylabel(r'$\mathrm{[m^2/s^2]}$')
        ax.legend(**lgnd_kws)

        ax = axs[3]
        ax.semilogy(t_hours, tke_eqn['eps'], 'ko', ms=3, zorder=4,
                    mfc='k', mec='none',
                    label='$\epsilon$')
        ax.set_ylabel('$\mathrm{[W/kg]}$')

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


if flg.get('turb time02'):

    with pt.style['twocol']():

        fig, axs = pt.newfig(202, 4, 1,
                             figsize=6,
                             hspace=0.17,
                             left=0.12, right=0.81,
                             bottom=0.08, top=0.97,
                             sharex=True, sharey=False)

        ax = axs[0]
        ax.semilogy(t_hours, dnow.tke, 'k-', label=r"tke")
        ax.semilogy(t_hours, dnow.upup_, 'b-', label=r"$\overline{u^2}$")
        ax.semilogy(t_hours, dnow.vpvp_, 'r-', label=r"$\overline{v^2}$")
        ax.semilogy(t_hours, dnow.wpwp_, 'g-', label=r"$\overline{w^2}$")
        ax.set_ylabel(r'$\mathrm{tke}\ \mathrm{[m^2/s^2]}$')
        ax.legend(**lgnd_kws)

        ax = axs[1]
        ax.plot(t_hours, dnow.upvp_, 'b-', label=r"$\overline{uv}$")
        ax.plot(t_hours, dnow.upwp_, 'r-', label=r"$\overline{uw}$")
        ax.plot(t_hours, dnow.vpwp_, 'g-', label=r"$\overline{vw}$")
        ax.plot(t_hours, dnow1.upvp_, 'b.', label=r"$\overline{uv}$")
        ax.plot(t_hours, dnow1.upwp_, 'r.', label=r"$\overline{uw}$")
        ax.plot(t_hours, dnow1.vpwp_, 'g.', label=r"$\overline{vw}$")
        ax.set_ylabel(r'$\mathrm{stress}\ \mathrm{[m^2/s^2]}$')
        #ax.legend(**lgnd_kws)

        ax = axs[2]
        ax.semilogy(t_hours, tke_eqn['eps'], 'k.')
        ax.set_ylabel('$\epsilon\ \mathrm{[W/kg]}$')

        ax = axs[3]
        inds = tke_eqn['Prod_uz'] > 0
        ax.semilogy(t_hours[inds], tke_eqn['Prod_uz'][inds], 'k.')
        ax.semilogy(t_hours[~inds], -tke_eqn['Prod_uz'][~inds], 'r.')
        ax.set_ylim([1e-6, 1e-2])
        #ax.set_ylabel('$\partial u/\partial z\ \mathrm{[s^{-1}]}$')

if flg.get('turb time03'):

    with pt.style['twocol']():

        fig, axs = pt.newfig(203, 4, 1,
                             figsize=6,
                             hspace=0.17,
                             left=0.12, right=0.81,
                             bottom=0.08, top=0.97,
                             sharex=True, sharey=True)

        ax = axs[0]
        ax.semilogy(t_hours, tke_eqn['eps'], 'ko', ms=3, zorder=4,
                    mfc='k', mec='none',
                    label='$\epsilon$')

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

        lkw = lgnd_kws.copy()
        lkw['numpoints'] = 1
        lkw['handletextpad'] = 0.2
        ax.legend(**lkw)

        ax = axs[1]
        clr = 'r'
        yvar = tke_eqn['Adv_v']
        inds = yvar > 0
        ax.semilogy(t_hours[inds],
                    yvar[inds],
                    'o', ms=3, zorder=5,
                    mfc=clr, mec='none', alpha=.6,
                    label='$A_{y}$')
        ax.semilogy(t_hours[~inds],
                    -yvar[~inds],
                    'o', ms=2.5, zorder=6,
                    mfc='none', mec=clr, alpha=0.6,
                    label='$-A_{y}$')

        clr = 'k'
        yvar = tke_eqn['Trans_y']
        inds = yvar > 0
        ax.semilogy(t_hours[inds],
                    yvar[inds],
                    'o', ms=3, zorder=5,
                    mfc=clr, mec='none', alpha=.6,
                    label=r'$T_{y}$')
        ax.semilogy(t_hours[~inds],
                    -yvar[~inds],
                    'o', ms=2.5, zorder=6,
                    mfc='none', mec=clr, alpha=0.6,
                    label=r'$-T_{y}$')

        # clr = 'b'
        # yvar = tke_eqn['Trans_z']
        # inds = yvar > 0
        # ax.semilogy(t_hours[inds],
        #             yvar[inds],
        #             'o', ms=3, zorder=5,
        #             mfc=clr, mec='none', alpha=.6,
        #             label=r'$T_{z}$')
        # ax.semilogy(t_hours[~inds],
        #             -yvar[~inds],
        #             'ko', ms=2.5, zorder=6,
        #             mfc='none', mec=clr, alpha=0.6,
        #             label=r'$-T_{z}$')

        clr = 'b'
        yvar = tke_eqn['Prod_uy']
        inds = yvar > 0
        ax.semilogy(t_hours[inds],
                    yvar[inds],
                    'o', ms=3, zorder=5,
                    mfc=clr, mec='none', alpha=.6,
                    label=r'$P_{uy}$')
        ax.semilogy(t_hours[~inds],
                    -yvar[~inds],
                    'ko', ms=2.5, zorder=6,
                    mfc='none', mec=clr, alpha=0.6,
                    label=r'$-P_{uy}$')

        clr = 'g'
        yvar = tke_eqn['dt']
        inds = yvar > 0
        ax.semilogy(t_hours[inds],
                    yvar[inds],
                    'o', ms=3, zorder=5,
                    mfc=clr, mec='none', alpha=.6,
                    label=r'$dt$')
        ax.semilogy(t_hours[~inds],
                    -yvar[~inds],
                    'ko', ms=2.5, zorder=6,
                    mfc='none', mec=clr, alpha=0.6,
                    label=r'$-dt$')

        lkw = lgnd_kws.copy()
        lkw['numpoints'] = 1
        lkw['handletextpad'] = 0.2
        ax.legend(**lkw)

        ax = axs[2]
        yvar = (tke_eqn['Trans_y'] +
                tke_eqn['Prod_uz'] +
                tke_eqn['dt'] +
                #tke_eqn['Prod_uy'] +
                tke_eqn['Adv_v']
        )
        inds = yvar > 0
        ax.semilogy(t_hours[inds],
                    yvar[inds],
                    'o', ms=3, zorder=5,
                    mfc='k', mec='none', alpha=.6,
                    label=r'$\sum$')
        ax.semilogy(t_hours[~inds],
                    -yvar[~inds],
                    'o', ms=2.5, zorder=6,
                    mfc='none', mec='k', alpha=0.6,
                    label=r'$-\sum$')

        ax.semilogy(t_hours, tke_eqn['eps'], 'ko', ms=3, zorder=4,
                    mfc='k', mec='none',
                    label='$\epsilon$')

        lkw = lgnd_kws.copy()
        lkw['numpoints'] = 1
        lkw['handletextpad'] = 0.2
        ax.legend(**lkw)

        for iax, ax in enumerate(axs):
            ax.set_ylabel('$\mathrm{[W/kg]}$')
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

        axs[-1].set_ylim([1e-6, 1e-2])
        axs[-1].set_xlim([t_hours[0], t_hours[-1]])
        axs[-1].set_xlabel('Time [Hours since 00:00 June 18, 2014]')

        if flg.get('save figs'):
            fig.savefig(pt.figdir + 'TurbTime_TTM_03.pdf')

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
        # axs.legend(loc='upper left', numpoints=1, handlelength=1,
        #            handletextpad=0.2,
        #            prop=dict(size='medium'))

        axs.plot([1e-8, 1], [1e-8, 1], 'k:')
        axs.set_ylim([1e-6, 1e-2])
        axs.set_xlim([1e-6, 1e-2])
        axs.set_xlabel('$\epsilon\ \mathrm{[W/kg]}$')
        axs.set_ylabel('$P_{uz}\ \mathrm{[W/kg]}$')

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
            ax.set_ylabel('$P_{uz}\ \mathrm{[W/kg]}$')
            ax.legend(loc='upper left', numpoints=1, handlelength=1,
                      handletextpad=0.2,
                      prop=dict(size='medium'))

        # for rng in zip(logbins[:-1], logbins[1:]):
        #     inds = (rng[0] < tke_eqn['eps']) & (tke_eqn['eps'] < rng[1])
        #     xval = np.mean(rng)
        #     plt.plot(xval, tke_eqn['Prod_uz'][inds].mean(), 'r.')
        #     plt.plot(xval, -tke_eqn['Prod_uz'][inds].mean(), 'ro')

        axs[-1].set_xlabel('$\epsilon\ \mathrm{[W/kg]}$')
        axs[0].set_ylim([1e-6, 1e-2])
        axs[0].set_xlim([1e-6, 1e-2])

    if flg.get('save figs'):
        fig.savefig(pt.figdir + 'EpsVProd03.png')

if flg.get('eqnVeps01', False):

    fig, axs = pt.newfig(566, 1, 1,
                         figsize=3,
                         hspace=0.17,
                         # left=0.2, right=0.9,
                         # bottom=0.1, top=0.95,
                         sharex=True, sharey=True)
    ax = axs
    ax.loglog(tke_eqn['eps'], tke_eqn['Adv_v'], 'ko', mfc='k', mec='none')
    ax.loglog(tke_eqn['eps'], -tke_eqn['Adv_v'], 'ko', mfc='none', mec='k')
    ax.loglog(tke_eqn['eps'], tke_eqn['Trans_y'], 'ro', mfc='r', mec='none')
    ax.loglog(tke_eqn['eps'], -tke_eqn['Trans_y'], 'ro', mfc='none', mec='r')
    ax.loglog(tke_eqn['eps'], tke_eqn['Trans_z'], 'bo', mfc='b', mec='none')
    ax.loglog(tke_eqn['eps'], -tke_eqn['Trans_z'], 'bo', mfc='none', mec='b')
    ax.plot([1e-8, 100], [1e-8, 100], 'k--')
    ax.set_xlim([1e-6, 1e-2])
    ax.set_ylim([1e-6, 1e-2])

if flg.get('eqnVeps02'):

    with pt.style['onecol']():

        fig, axs = pt.newfig(306, 2, 1,
                             figsize=5,
                             hspace=0.17,
                             left=0.2, right=0.9,
                             bottom=0.1, top=0.95,
                             sharex=True, sharey=True)

        umag = np.abs(dnow.u)
        iumag = np.abs(dnow.u) > 1.0
        iupos = dnow.u > 0

        yvar = (tke_eqn['Trans_y'] +
                tke_eqn['Prod_uz'] +
                tke_eqn['dt'] +
                tke_eqn['Adv_v'])
        inds = (yvar > 0)
        for iax, ax in enumerate(axs):
            if iax == 0:
                inow = iumag & iupos
                ax.set_title('Ebb')
            else:
                ax.set_title('Flood')
                inow = iumag & ~iupos
            ax.loglog(tke_eqn['eps'][inds & inow],
                      yvar[inds & inow], 'o',
                      label='$\sum>0,\mathrm{(N = %d)}$' %
                      (inds & inow).sum(),
                      mfc='k', mec='none', ms=3,
                      zorder=2,)
            ax.loglog(tke_eqn['eps'][~inds & inow],
                      -yvar[~inds & inow],
                      'o',
                      label='$\sum<0,\mathrm{(N = %d)}$' %
                      (~inds & inow).sum(),
                      mfc='none', mec='k', ms=2.5,
                      zorder=-2)
            ax.plot([1e-8, 1], [1e-8, 1], 'k:')
            ax.set_ylabel('$\sum \ \mathrm{[W/kg]}$')
            ax.legend(loc='upper left', numpoints=1, handlelength=1,
                      handletextpad=0.2,
                      prop=dict(size='medium'))

        axs[-1].set_xlabel('$\epsilon\ \mathrm{[W/kg]}$')
        axs[0].set_ylim([1e-6, 1e-2])
        axs[0].set_xlim([1e-6, 1e-2])

    if flg.get('save figs'):
        fig.savefig(pt.figdir + 'EqnVeps02.pdf')
