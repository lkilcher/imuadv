#from dolfyn.adv import api as avm
from faretools.ptools import figobj
from dolfyn.tools.misc import delta
import numpy as np
import ttm.June2012 as j12

awacr = j12.load('AWAC', coordsys='pax')
awac = j12.load('AWAC', coordsys='pax', bin=True)

# This file was copied from ~/work/mhk/ttm/June2012/mets_figs.py

datdir = '/Users/lkilcher/data/ttm/June2012/'
figdir = '../fig/'

flg = {}
# flg['do_specfig'] = True
# flg['do_specfig2'] = True
# flg['do_specfig_all'] = True
# flg['do_time_all'] = True
# flg['do_time_u'] = True
# flg['do_time_u2'] = True
flg['do_coh'] = True
# flg['do_specfig_comp'] = True
flg['do_Lcoh'] = True
# flg['do_Lcoh2'] = True

# fctr_lcoh=2*np.pi
fctr_lcoh = 1

if 'dat_mc' not in vars():
    dat_mc = j12.load('NREL', coordsys='pax', bin=True)
    datr = j12.load('NREL', coordsys='pax', bin=False)
    #dat_mc = avm.load(datdir + 'TTM_Vectors/TTM_NRELvector_Jun2012_pax_b5m.h5')
    #datr = avm.load(datdir + 'TTM_Vectors/TTM_NRELvector_Jun2012_pax.h5')
    # binner_nr = avm.turb_binner(
    #     n_bin=10258, fs=datr.fs, n_fft=10240, n_fft_coh=2048)
    # acov = binner_nr.calc_acov(datr._u)
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

bd = np.abs(awac.w[10] - dat_mc.w) > 0.04
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
    fg0 = figobj(32, (3, 1), axsize=[1.6, 4],
                 sharex=True, gap=.24, frame=[.6, .2, 1., .2])
    ax = fg0.ax[0]
    ax.plot(dat_mc.time, dat_mc.u, 'b-', label='ADV')
    ax.plot(awac.time, np.ma.masked_where(
        bd, tmp.real[10]), 'r-', label='AWAC')
    ax.set_ylim([-2.5, 2.5])
    c = '\n'
    ax.set_ylabel(r'$\bar{u}\,\mathrm{[ms^{-1}]}$', size='large')
    # ax.set_ylabel(r'$\frac{\bar{u}}{\mathrm{[ms^{-1}]}}$',size='x-large',rotation=0,multialignment='center')
    ax = fg0.ax[1]
    ax.plot(dat_mc.time, dat_mc.v, 'b-')
    ax.plot(awac.time, np.ma.masked_where(bd, tmp.imag[10]), 'r-')
    ax.set_ylim([-1, 1])
    ax.set_ylabel(r'$\bar{v}\,\mathrm{[ms^{-1}]}$', size='large')
    ax = fg0.ax[2]
    ax.plot(dat_mc.time, dat_mc.w, 'b-')
    ax.plot(awac.time, np.ma.masked_where(bd, awac.w[10] + 0.02), 'r-')
    ax.set_ylim([-0.1, 0.1])
    ax.set_ylabel(r'$\bar{w}\,\mathrm{[ms^{-1}]}$', size='large')
    ax.set_xlim([12.7, 14.7])
    fg0.ax[0].legend(loc='upper left', prop={
                     'size': 'small'}, bbox_to_anchor=[0.03, 0.17])
    # ax=fg0.ax[3]
    # ax.plot(dat_mc.time,dat_mc.U_mag,'b-',label='ADV')
    # ax.plot(awac.time,np.ma.masked_where(bd,np.abs(tmp)[10]),'r-',label='AWAC')
    # ax.set_ylabel(r'$\bar{U}\,\mathrm{[ms^{-1}]}$',size='large')
    # ax.set_ylim([-2.5,2.5])
    ax.set_xlabel('Day of June, 2012')
    fg0.hide('xticklabels', ax)
    fg0.sax.alphNumAxes('ABCD', prefix='(', suffix=')', pos='ul')
    # fg0.savefig(figdir+'timeFig01.pdf')

pt = 35
tr = dat_mc.mpltime[pt] + np.array([-0.5, 0.5]) * \
    (dat_mc.mpltime[pt + 1] - dat_mc.mpltime[pt])
if flg.get('do_time_u', False):
    fg0 = figobj(33, (2, 1), axsize=[1.6, 4],
                 sharex=False, gap=.64, frame=[.6, .2, .8, .2])
    ax = fg0.ax[0]
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
    ax = fg0.ax[1]
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
    fg0.ax[0].legend(loc='upper left', prop={'size': 'small'}, bbox_to_anchor=[
                     0.29, .99], handlelength=1.5, handletextpad=.2)
    fg0.sax.alphNumAxes('ABC', prefix='(', suffix=')', pos='ul', offset=(5, 5))
    fg0.savefig(figdir + 'timeFigU01.pdf')
    fg0.savefig(figdir + 'timeFigU01.png', dpi=dpi)


if flg.get('do_time_u2', False):
    fg00 = figobj(34, (2, 1), axsize=[1.6, 4],
                  sharex=True, gap=.2, frame=[.6, .2, .8, .2])
    ax = fg00.ax[0]
    ax.plot(dat_mc.time, dat_mc.u, 'b-', label='ADV')
    ax.plot(awac.time, np.ma.masked_where(bd, awac.u[10]), 'r-', label='AWAC')
    ax.set_ylim([-2.5, 2.5])
    ax.hln(0, color='k', linestyle='-', zorder=10)
    c = '\n'
    ax.set_ylabel(r'$\bar{u}\,\mathrm{[ms^{-1}]}$', size='large')
    # ax.plot(dat_mc.time[pt],dat_mc.u[pt],'ko',mfc='none',mew=2,ms=10)
    ax = fg00.ax[1]
    ax.plot(dat_mc.time, np.abs(dat_mc.upup_)**0.5 - 0.02, 'b-')
    ax.plot(awac.time, np.ma.masked_where(
        bd, np.abs(awac.upup_[10])**0.5), 'r-')
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
    fg00.savefig(figdir + 'timeFigU02.pdf')
    fg00.savefig(figdir + 'timeFigU02.png', dpi=dpi)

if flg.get('do_awacVadv_up', False):
    fg20 = figobj(388, (1, 1))
    fg20.plot(np.abs(dat_mc.upup_)**0.5,
              np.ma.masked_where(bd, np.abs(awac.upup_[10]))**0.5 + .12 - .1,
              'k.')
    fg20.plot([0, 0.2], [0, 0.2], 'k-')

if flg.get('do_check_awac_amp', False):
    def msk(dat):
        return np.ma.masked_where(bd, dat)

    fga = figobj(40, [4, 1])
    thresh = 100
    d = awacr._amp[:, 10].mean(0)
    bd = d > thresh
    ax = fga.ax[0]
    ax.plot(msk(awacr.u[10]))
    ax = fga.ax[1]
    ax.plot(msk(awacr.v[10]))
    ax = fga.ax[2]
    ax.plot(msk(awacr.w[10]))
    ax = fga.ax[3]
    ax.plot(d)
    ax.hln(thresh, color='r')
