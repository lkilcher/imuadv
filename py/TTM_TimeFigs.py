import ptools as pt
import dolfyn.adv.api as avm
from dolfyn.tools.misc import delta
import numpy as np
import dolfyn.data.binned as binmod
import data.ttm_june2012 as j12
plt = pt.plt
import scipy.signal as sig
import matplotlib.ticker as tkr

awacr = j12.load('awac', coordsys='pax', data_groups='ALL')
awac = j12.load('awac', coordsys='pax', bin=True)


def rms(arr, ax=-1):
    return np.sqrt(np.mean(arr ** 2))

# This file was copied from ~/work/mhk/ttm/June2012/mets_figs.py

datdir = '/Users/lkilcher/data/ttm/June2012/'

binner = avm.TurbBinner(9600, 32)

flg = {}
# flg['do_time_all'] = True
# flg['do_time_u'] = True
# flg['do_time_u2'] = True
#flg['do_check_awac_amp'] = True
#flg['do_awacVadv_up'] = True
flg['show awac avg fix'] = True
#flg['show awac avg2'] = True
flg['save figs'] = True

if 'dat_mc' not in vars():
    dat_mc = j12.load('adv-nrel', coordsys='pax', bin=True)
    datr = j12.load('adv-nrel', coordsys='pax', bin=False)
    #dat_mc = avm.load(datdir + 'TTM_Vectors/TTM_NRELvector_Jun2012_pax_b5m.h5')
    #datr = avm.load(datdir + 'TTM_Vectors/TTM_NRELvector_Jun2012_pax.h5')
    # binner_nr = avm.turb_binner(
    #     n_bin=10258, fs=datr.fs, n_fft=10240, n_fft_coh=2048)
    # acov = binner_nr.calc_acov(datr.vel)
    # lint = binner_nr.calc_Lint(acov, dat_mc.U_mag)

frmt = 'p'
axsz = .7 * 5
# frmt='h';axsz=.7*5
flsfx = '.png'
dpi = 200

vel = 1.  # This is the default value, it can be customized elsewhere.
dl = 0.1
indu_mc = delta(dat_mc.u, dl)


f_fit = np.logspace(-2, 1, 50)
rad_hz = 2 * np.pi

#bd = np.abs(awac.w[10] - dat_mc.w) > 0.04
#bd = np.abs(awac.w[10] - dat_mc.w) > 0.04

# ucomp=1
idi = 1
lm = {'k': [1e-2, 1e1],
      'f': [1e-3, 2],
      'coh': [0, 1],
      'spec': [1e-5, 1e-1],
      'snr': [1e-2, 10],
      'nsr': [0, 1],
      'spec_hz': [1e-4, 3e-1]}
shp = {'h': [1, 2], 'p': [2, 1]}


ang = 4  # deg
tmp = awac.rotate_var(ang * np.pi / 180)

if flg.get('do_time_all', False):
    fg, axs = pt.newfig(32, 3, 1,
                        sharex=True, )
    ax = axs[0]
    ax.plot(dat_mc.mpltime, dat_mc.u, 'b-', label='ADV')
    #ax.plot(awac.time, np.ma.masked_where(bd, tmp.real[10]), 'r-', label='AWAC')
    ax.set_ylim([-2.5, 2.5])
    c = '\n'
    ax.set_ylabel(r'$\bar{u}\,\mathrm{[ms^{-1}]}$', size='large')
    # ax.set_ylabel(r'$\frac{\bar{u}}{\mathrm{[ms^{-1}]}}$',size='x-large',rotation=0,multialignment='center')
    ax = axs[1]
    ax.plot(dat_mc.mpltime, dat_mc.v, 'b-')
    #ax.plot(awac.time, np.ma.masked_where(bd, tmp.imag[10]), 'r-')
    ax.set_ylim([-1, 1])
    ax.set_ylabel(r'$\bar{v}\,\mathrm{[ms^{-1}]}$', size='large')
    ax = axs[2]
    ax.plot(dat_mc.mpltime, dat_mc.w, 'b-')
    #ax.plot(awac.time, np.ma.masked_where(bd, awac.w[10] + 0.02), 'r-')
    #ax.set_ylim([-0.1, 0.1])
    ax.set_ylabel(r'$\bar{w}\,\mathrm{[ms^{-1}]}$', size='large')
    #ax.set_xlim([12.7, 14.7])
    axs[0].legend(loc='upper left', prop={'size': 'small'},
                  bbox_to_anchor=[0.03, 0.17])
    # ax=axs[3]
    # ax.plot(dat_mc.time,dat_mc.U_mag,'b-',label='ADV')
    # ax.plot(awac.time,np.ma.masked_where(bd,np.abs(tmp)[10]),'r-',label='AWAC')
    # ax.set_ylabel(r'$\bar{U}\,\mathrm{[ms^{-1}]}$',size='large')
    # ax.set_ylim([-2.5,2.5])
    ax.set_xlabel('Day of June, 2012')
    #fg.hide('xticklabels', ax)
    #fg.sax.alphNumAxes('ABCD', prefix='(', suffix=')', pos='ul')
    # fg0.savefig(pt.figdir+'timeFig01.pdf')

point = 35
tr = dat_mc.mpltime[point] + np.array([-0.5, 0.5]) * \
    (dat_mc.mpltime[point + 1] - dat_mc.mpltime[point])
if flg.get('do_time_u', False):
    fg, axs = pt.newfig(33, 2, 1, sharex=False, )
    ax = axs[0]
    ax.plot(dat_mc.time, dat_mc.u, 'b-', label='ADV')
    ax.plot(awac.time, np.ma.masked_where(bd, awac.u[10]), 'r-', label='AWAC')
    ax.set_ylim([-2.5, 2.5])
    ax.set_xlabel('Day of June, 2012', backgroundcolor=[
                  1, 1, 1, 1], zorder=10,)
    ax.hln(0, color='k', linestyle='-', zorder=10)
    c = '\n'
    ax.set_ylabel(r'$\bar{u}\,\mathrm{[ms^{-1}]}$', size='large')
    # ax.plot(dat_mc.time[pt],dat_mc.u[pt],'ko',mfc='none',mew=2,ms=10)
    ax.plot(dat_mc.time[pt], dat_mc.u[pt], 'k.', mec='none', ms=10)
    ax.annotate('', [dat_mc.time[pt], dat_mc.u[pt]],
                [0, -.4],
                textcoords='axes fraction',
                xycoords='data',
                arrowprops=dict(arrowstyle="-",
                                connectionstyle="angle,angleA=53,angleB=-105,rad=6",
                                linewidth=1.6,
                                shrinkA=0, color='0.3'),
                zorder=1)
    # ax.annotate('', [dat_mc.time[pt], dat_mc.u[pt]],
    #             [1, -.4],
    #             textcoords='axes fraction',
    #             xycoords='data',
    #             arrowprops=dict(arrowstyle="-",
    #                             connectionstyle="angle3,angleA=-3,angleB=-50",
    #                             linewidth=1.6,
    #                             shrinkA=0,
    #                             color='0.3'),
    #             zorder=-1)
    ax.annotate('', [dat_mc.time[pt], dat_mc.u[pt]],
                [1, -.4],
                textcoords='axes fraction',
                xycoords='data',
                arrowprops=dict(arrowstyle="-",
                                connectionstyle="angle,angleA=-3,angleB=-100,rad=6",
                                linewidth=1.6,
                                shrinkA=0,
                                color='0.3'),
                zorder=-1)
    # ax.set_ylabel(r'$\frac{\bar{u}}{\mathrm{[ms^{-1}]}}$',size='x-large',rotation=0,multialignment='center')
    ax.set_xlim([12.7, 14.7])
    ax = axs[1]
    iturb = (tr[0] < datr.mpltime) & (datr.mpltime < tr[1])
    ax.plot((datr.time[iturb] - datr.time[iturb][0]) * 24 *
            3600, datr.u[iturb] - datr.u[iturb].mean(), 'b-')
    ax.hln(0, color='k', linestyle='-', zorder=10)
    stdev = datr.u[iturb].std()
    # ax.hln(datr.u[iturb].mean(),color='k',linestyle='-',zorder=10)
    # ax.hln(datr.u[iturb].mean()-stdev,color='k',linestyle='--',zorder=10)
    # ax.hln(datr.u[iturb].mean()+stdev,color='k',linestyle='--',zorder=10)
    ax.hln(-stdev, color='k', linestyle='--', zorder=10)
    ax.hln(stdev, color='k', linestyle='--', zorder=10)
    iturb = (tr[0] < awacr.mpltime) & (awacr.mpltime < tr[1])
    ax.plot((awacr.time[iturb] - awacr.time[iturb][0]) * 24 * 3600,
            awacr.u[10, iturb] - awacr.u[10, iturb].mean(), 'r-', linewidth=.7)
    ax.set_yticks(np.arange(-2, 2, .2))
    # ax.set_ylim([-1.6,-0.8])
    ax.set_ylim([-0.40001, 0.4])
    ax.set_ylabel("$u'\,\mathrm{[ms^{-1}]}$")
    ax.set_xlabel('Time [seconds]')
    ax.set_xlim([0, 300])
    axs[0].legend(loc='upper left', prop={'size': 'small'},
                  bbox_to_anchor=[0.29, .99], handlelength=1.5, handletextpad=.2)
    fg.sax.alphNumAxes('ABC', prefix='(', suffix=')', pos='ul', offset=(5, 5))
    fg.savefig(pt.figdir + 'timeFigU01.pdf')
    fg.savefig(pt.figdir + 'timeFigU01.png', dpi=dpi)


if flg.get('do_time_u2', False):
    fg00 = pt.newfig(34, (2, 1), sharex=True, )
    ax = fg00.ax[0]
    ax.plot(dat_mc.time, dat_mc.u, 'b-', label='ADV')
    #ax.plot(awac.time, np.ma.masked_where(bd, awac.u[10]), 'r-', label='AWAC')
    ax.set_ylim([-2.5, 2.5])
    ax.hln(0, color='k', linestyle='-', zorder=10)
    c = '\n'
    ax.set_ylabel(r'$\bar{u}\,\mathrm{[ms^{-1}]}$', size='large')
    # ax.plot(dat_mc.time[pt],dat_mc.u[pt],'ko',mfc='none',mew=2,ms=10)
    ax = fg00.ax[1]
    ax.plot(dat_mc.time, np.abs(dat_mc.upup_)**0.5 - 0.02, 'b-')
    #ax.plot(awac.time, np.ma.masked_where(bd, np.abs(awac.upup_[10])**0.5), 'r-')
    ax.set_yticks(np.arange(0, 1, .1))
    ax.set_ylim([0, .2])
    ax.set_ylabel("$\mathrm{std}(u')\,\mathrm{[ms^{-1}]}$")
    ax.set_xlim([12.7, 14.7])
    ax.set_xlabel('Day of June, 2012', backgroundcolor=[
                  1, 1, 1, 1], zorder=10,)
    fg00.hide('xticklabels', ax)
    fg00.ax[0].legend(loc='upper left', prop={'size': 'small'}, bbox_to_anchor=[
                      0.29, .99], handlelength=1.5, handletextpad=.2)
    fg00.sax.alphNumAxes('ABC', prefix='(', suffix=')',
                         pos='ul', offset=(5, 5))
    fg00.savefig(pt.figdir + 'timeFigU02.pdf')
    fg00.savefig(pt.figdir + 'timeFigU02.png', dpi=dpi)

if flg.get('do_awacVadv_up', False):

    with pt.style['classic']():
        fg20 = pt.newfig(388, 1, 1, figsize=(8, 8))
        # fg20.plot(np.abs(dat_mc.upup_)**0.5,
        #           np.ma.masked_where(bd, np.abs(awac.upup_[10]))**0.5 + .12 - .1,
        #           'k.')
        fg20.plot([0, 0.2], [0, 0.2], 'k-')

if flg.get('do_check_awac_amp', False):

    def msk(dat):
        return np.ma.masked_where(bd, dat)

    awi = slice(25000, 32000)
    depi = 8
    with pt.style['classic']():
        fg, axs = pt.newfig(40, 4, 1, figsize=[8, 8], sharex=True)
        thresh = 100
        d = awacr.amp[:, depi, awi].mean(0)
        bd = d > thresh
        ax = axs[0]
        ax.plot(msk(awacr.u[depi, awi]))
        ax = axs[1]
        ax.plot(msk(awacr.v[depi, awi]))
        ax = axs[2]
        ax.plot(msk(awacr.w[depi, awi]))
        ax = axs[3]
        ax.plot(d, 'k')
        ax.plot(awacr.amp[0, depi, awi], 'r')
        ax.plot(awacr.amp[1, depi, awi], 'b')
        ax.plot(awacr.amp[2, depi, awi], 'g')
        ax.axhline(thresh, color='r')


def screen(amp, npt=10, thresh=0.05):
    damp = np.diff(np.pad(amp.astype(np.float32),
                          (((0, 1),) + (((0, 0), ) * (amp.ndim - 1))),
                          'wrap'),
                   1, 0)
    damp /= amp
    shape = [1, ] * amp.ndim
    shape[1] = npt
    tmp = np.abs(sig.convolve(damp, np.ones(shape) / npt, 'same'))
    return tmp.max(0) > thresh

if flg.get('do_check_awac_amp2', False):

    def msk(dat):
        return np.ma.masked_where(bd, dat)

    awi = slice(25000, 32000)
    t = pt.dt.num2date(awacr.mpltime[awi])
    avi = slice(np.argmin(np.abs(awacr.mpltime[awi.start] - datr.mpltime)),
                np.argmin(np.abs(awacr.mpltime[awi.stop] - datr.mpltime)), )
    ta = pt.dt.num2date(datr.mpltime[avi])
    depi = 9
    with pt.style['classic']():
        fg, axs = pt.newfig(41, 5, 1, figsize=[8, 8], sharex=True)
        thresh = 100
        d = awacr.amp[:, depi, awi].mean(0)
        #damp = np.diff(awacr.amp[:, depi, awi].astype(np.int16), 1, 0)
        #bd = d > thresh
        damp = np.diff(np.pad(awacr.amp[:, depi, awi].astype(np.int16),
                              ((0, 1), (0, 0)), 'wrap'),
                       1, 0)
        tmp = sig.convolve(damp, np.ones((1, 10)) / 10, 'same')
        bd = screen(awacr.amp[:, depi, awi])
        ax = axs[0]
        ax.plot(t, awacr.u[depi, awi], '0.6')
        ax.plot(ta, datr.u[avi], 'r')
        ax.plot(t, msk(awacr.u[depi, awi]), 'b')
        ax = axs[1]
        ax.plot(t, awacr.v[depi, awi], '0.6')
        ax.plot(ta, datr.v[avi], 'r')
        ax.plot(t, msk(awacr.v[depi, awi]), 'b')
        ax = axs[2]
        ax.plot(t, awacr.w[depi, awi], '0.6')
        ax.plot(ta, datr.w[avi], 'r')
        ax.plot(t, msk(awacr.w[depi, awi]), 'b')
        ax = axs[3]
        ax.plot(t, d, 'k')
        ax.plot(t, awacr.amp[0, depi, awi], 'r')
        ax.plot(t, awacr.amp[1, depi, awi], 'b')
        ax.plot(t, awacr.amp[2, depi, awi], 'g')
        ax.axhline(thresh, color='r')
        ax = axs[4]
        #ax.plot(t, d, 'k')
        ax.plot(t, damp[0], 'r')
        ax.plot(t, damp[1], 'b')
        ax.plot(t, damp[2], 'g')
        #ax.plot(t, np.abs(damp).max(0), 'k')
        ax.plot(t, tmp.max(0), 'k')
        
        #ax.axhline(thresh, color='r')

if flg.get('show awac avg fix', False):

    bnr = binmod.TimeBinner(300, 1)
    t0 = awacr.mpltime[0]
    t = (bnr.mean(awacr.mpltime) - t0) * 24
    dat_mc.time = (dat_mc.mpltime - t0) * 24
    depi = 9
    bd = screen(awacr.amp[:, depi], thresh=0.05)
    u = bnr.mean(np.ma.masked_where(np.tile(bd[None, :], (3, 1)),
                                    awacr.vel[:, depi]),
                 mask_thresh=0.5)
    scale = np.array([[-2.1, 2.1],
                      [-1, 1],
                      [-0.2, 0.2], ])
    height_ratios = np.diff(scale, 1, -1)
    height_ratios[2] *= 2
    height_ratios = np.ones(3)
    offset = [0, 0, 0.02]
    with pt.style['onecol']():
        fg, axs = pt.newfig(42, 3, 1, figsize=4, sharex=True, hspace=0.18,
                            bottom=0.1,
                            gridspec_kw=dict(height_ratios=height_ratios))
        for iax, ax in enumerate(axs):
            ax.text(0.03, 0.92, '({})'.format('ABC'[iax]),
                    transform=ax.transAxes, fontsize='small',
                    ha='left', va='top')
            ax.axhline(0, color='k', linestyle='--', lw=0.5)
            ax.set_ylabel(r'$\bar{%s}\ \mathrm{[m/s]}$' % format('uvw'[iax]))
            ax.plot(dat_mc.time, dat_mc.vel[iax], 'k', lw=2)
            ax.plot(t, u[iax] + offset[iax], 'r', lw=1.0)
        ax = axs[0]
        ax.yaxis.set_ticks(np.arange(-3.0, 3.5, 1.0))
        ax.yaxis.set_ticks(np.arange(-3.5, 3.5, 0.1), minor=True)
        ax.yaxis.set_major_formatter(tkr.FormatStrFormatter('$%0.1f$'))
        ax.set_ylim(scale[0])
        ax = axs[1]
        ax.yaxis.set_ticks(np.arange(-1.5, 1.5, 0.5))
        ax.yaxis.set_ticks(np.arange(-1.5, 1.5, 0.1), minor=True)
        ax.set_ylim(scale[1])
        ax = axs[2]
        ax.yaxis.set_ticks(np.arange(-0.5, 0.5, 0.1))
        #ax.yaxis.set_ticks(np.arange(-0.5, 0.5, 0.05), minor=True)
        ax.set_ylim(scale[2])
        ax.xaxis.set_ticks(np.arange(0, 50, 6))
        # ax.xaxis.set_major_locator(pt.dt.HourLocator(range(0, 24, 6)))
        # ax.xaxis.set_major_formatter(pt.dt.DateFormatter('%H'))
        ax.set_xlim([t[0], t[-1]])
        ax.set_xlabel('Time [Hours]')

        if flg.get('save figs', False):
            fg.savefig(pt.figdir + 'TimeFig02.pdf')

if flg.get('show awac avg2', False):

    def velmag(u):
        return np.sqrt((u ** 2).sum(0))

    def velang(u):
        out = np.angle(u[0] + 1j * u[1]) * 180 / np.pi
        if isinstance(u, np.ma.MaskedArray):
            out = np.ma.masked_where(u.mask[0], out)
        return out

    def dang(u1, u2):
        u1c = u1[0] + 1j * u1[1]
        u2c = u2[0] + 1j * u2[1]
        out = np.angle(u2c / u1c) * 180 / np.pi
        if isinstance(u1, np.ma.MaskedArray):
            out = np.ma.masked_where(u1c.mask, out)
        if isinstance(u2, np.ma.MaskedArray):
            out = np.ma.masked_where(u2c.mask, out)
        return out

    ui = np.array((np.interp(t, dat_mc.time, dat_mc.u),
                   np.interp(t, dat_mc.time, dat_mc.v),
                   np.interp(t, dat_mc.time, dat_mc.w)))
    bnr = binmod.TimeBinner(300, 1)
    t = bnr.mean(awacr.mpltime)
    depi = 9
    bd = screen(awacr.amp[:, depi])
    u = bnr.mean(np.ma.masked_where(np.tile(bd[None, :], (3, 1)), awacr.vel[:, depi]),
                 mask_thresh=0.6)
    umag = velmag(u)
    delta_ang = dang(u, ui)
    delta_ang.mask |= umag < 0.5
    with pt.style['onecol']():
        fg, axs = pt.newfig(43, 3, 1, figsize=5, sharex=True)
        thresh = 100
        ax = axs[0]
        ax.plot(dat_mc.mpltime, velmag(dat_mc.vel), 'b', lw=2)
        ax.plot(t, velmag(u), 'r')
        ax = axs[1]
        ax.plot(dat_mc.mpltime, velang(dat_mc.vel), 'b', lw=2)
        ax.plot(t, velang(u), 'r')
        ax = axs[2]
        ax.plot(t, delta_ang, 'r')
        
