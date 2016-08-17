import ptools as pt
import ttm.June2014 as j14
#import matplotlib.pyplot as plt

flg = {}
flg['turb time01'] = True


if 'dat' not in vars():
    dat = j14.load('ttm02b-top', 'pax',
                   bin=True)

if flg.get('turb time01'):

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
