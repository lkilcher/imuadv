import ttm.June2014 as j14
import ptools as pt
import numpy as np
plt = pt.plt

if 'dat' not in vars():
    dat = j14.load('ttm02b-top', 'pax',
                   bin=True)

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
