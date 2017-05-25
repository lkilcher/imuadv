import data.ttm_june2014 as j14
import dolfyn.adv.api as avm
import ptools as pt
import matplotlib.pyplot as plt
import numpy as np

if 'dat' not in vars():
    dat = j14.load('ttm02b-top', 'raw',
                   bin=False)

if 'bindat' not in vars():
    bindat = {}

filters = {
    '5m': 0.003,
    '100s': 0.01,
    '30s': 0.03,
    '10s': 0.1,
    '3s': 0.3,
}

velranges = [(0, 0.5),
             (1, 1.5),
             (2, 2.5)]

noise = [1.5e-4, 1.5e-4, 1.5e-5, ]

binner = avm.TurbBinner(5 * 60 * dat.fs, dat.fs)
for ftag, filt in filters.iteritems():
    if ftag not in bindat:
        dnow = dat.copy()
        mc = avm.motion.CorrectMotion(filt)
        mc(dnow)
        avm.rotate.earth2principal(dnow)
        bd = binner(dnow)
        bd['Spec'] = binner.psd(dnow.vel)
        bd['Spec_umot'] = binner.psd(dnow.velacc + dnow.velrot)
        bindat[ftag] = bd

pt_kws = {
    '100s': dict(color='c', zorder=-2),
    '30s': dict(color='k', lw=3, zorder=-10),
    '10s': dict(color='b', zorder=-4),
    '3s': dict(color='r', zorder=-5)}

with pt.style['twocol']():

    fig, axs = pt.newfig(4101, 3, len(velranges),
                         figsize=5,
                         right=0.86, bottom=0.1,
                         sharex=True, sharey=True)

    freq = binner.calc_omega() / pt.pii

    for ftag in ['100s', '30s', '10s', '3s']:
        filt = filters[ftag]
        kws = pt_kws[ftag]
        bd = bindat[ftag]

        for icol, vr in enumerate(velranges):
            inds = (vr[0] < np.abs(bd.u)) & (np.abs(bd.u) < vr[1])
            for iax, ax in enumerate(axs[:, icol]):
                ax.loglog(freq,
                          bd['Spec'][iax, inds].mean(0) * pt.pii - noise[iax],
                          label=ftag, **kws)
                # ax.loglog(freq,
                #           bd['Spec_umot'][iax, inds].mean(0) * pt.pii,
                #           linestyle=':', **kws)
                ax.axvline(filt, color=kws['color'], linestyle=':',
                           alpha=0.5, zorder=-10)
            if vr[0] == 0:
                axs[0, icol].set_title(r"$ |\bar{u}| < %0.1f$" % vr[1],
                                       fontsize='medium')
            else:
                axs[0, icol].set_title(r"$%0.1f < |\bar{u}| < %0.1f$" % vr,
                                       fontsize='medium')

    for iax, ax in enumerate(axs[:, 0]):
        ax.set_ylabel('$\mathrm{[m^2s^{-2}/Hz]}$')

    ax.set_xlim([1e-3, 1])
    ax.set_ylim([1e-4, 1])
    ax.set_xlabel('$f \mathrm{[Hz]}$')

    #axs[0].legend(loc='lower left', prop=dict(size='medium'))
    axs[0, -1].legend(loc='upper left', bbox_to_anchor=[1.02, 1.0],
                      handlelength=1.4, handletextpad=0.4,
                      prop=dict(size='medium'))

    fig.savefig(pt.figdir + 'SpecFig04_filtering.pdf')


pt_kws = {
    '5m': dict(color='k',
               zorder=-11, label='0.003 Hz', lw=5),
    '100s': dict(color='c', zorder=-2,
                 label='0.01 Hz'),
    '30s': dict(color='b', lw=3, zorder=-10,
                label='0.03 Hz'),
    '10s': dict(color='b', zorder=-4,
                label='0.1 Hz'),
    '3s': dict(color='r', zorder=-5,
               label='0.3 Hz', lw=1)}

filts = [
    '5m',
    #'100s',
    '30s',
    #'10s',
    '3s']

with pt.style['onecol']():

    vr = (2, 2.5)
    fig, axs = pt.newfig(4102, 3, 1,
                         figsize=5,
                         right=0.7, bottom=0.1,
                         hspace=0.2,
                         sharex=True, sharey=True)

    binner = avm.TurbBinner(5 * 60 * dat.fs, dat.fs)
    freq = binner.calc_omega() / pt.pii
    bd = bindat['30s']
    inds = (vr[0] < np.abs(bd.u)) & (np.abs(bd.u) < vr[1])

    for ftag in filts:
        filt = filters[ftag]
        kws = pt_kws[ftag]
        bd = bindat[ftag]
        for iax, ax in enumerate(axs):
            ax.loglog(freq,
                      bd['Spec'][iax, inds].mean(0) * pt.pii - noise[iax],
                      **kws)
            # ax.loglog(freq,
            #           bd['Spec_umot'][iax, inds].mean(0) * pt.pii,
            #           linestyle=':', **kws)
            ax.axvline(filt, color=kws['color'], linestyle=':',
                       alpha=0.5, zorder=-10)
        if vr[0] == 0:
            axs[0].set_title(r"$ |\bar{u}| < %0.1f$" % vr[1],
                             fontsize='medium')
        else:
            axs[0].set_title(r"$%0.1f < |\bar{u}| < %0.1f$" % vr,
                             fontsize='medium')

    for iax, ax in enumerate(axs):
        ax.text(0.07, 0.1, r'$S\{%s\}$' % 'uvw'[iax],
                transform=ax.transAxes)
        ax.set_ylabel('$\mathrm{[m^2s^{-2}/Hz]}$')
        bd = bindat[filts[0]]
        ax.loglog(freq, bd.Spec_umot[iax, inds].mean(0) * pt.pii,
                  lw=3,
                  color='0.7', zorder=-20)
        bd = bindat['3s']
        # ax.loglog(freq, bd.Spec_umot[iax, inds].mean(0) * pt.pii,
        #           lw=1, ls='-',
        #           color='r', zorder=-20, alpha=0.3)
    ax.set_xlim([1e-3, 1])
    ax.set_ylim([1e-3, 1])
    ax.set_xlabel('$f \mathrm{[Hz]}$')

    #axs[0].legend(loc='lower left', prop=dict(size='medium'))
    axs[0].legend(loc='upper left', bbox_to_anchor=[1.02, 1.0],
                  handlelength=1.4, handletextpad=0.4,
                  prop=dict(size='medium'))

    fig.savefig(pt.figdir + 'SpecFig04_filtering2.pdf')
