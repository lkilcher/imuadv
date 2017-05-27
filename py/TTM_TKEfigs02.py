import ptools as pt
import data.ttm_june2014 as j14
import numpy as np
import gis

flg = {}
#flg['turb time02'] = True
#flg['turb time03'] = True
#flg['epsVprod02'] = True
#flg['epsVprod03'] = True
#flg['eqnVeps01'] = True
#flg['eqnVeps02'] = True

lgnd_kws = dict(loc='upper left', bbox_to_anchor=(1.02, 1), )


if 'dat' not in vars():
    dat = j14.load('ttm02b-top', 'pax',
                   bin=True)
    dat2 = j14.load('ttm02b-bot', 'pax',
                    bin=True)
    # This data isn't downloaded/processed by default. Update the
    # setup_data.py script to get this data.
    dat3 = j14.load('ttm01b-top', 'pax',
                    bin=True).subset(slice(3, -2))
    dat4 = j14.load('ttm01b-bot', 'pax',
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
