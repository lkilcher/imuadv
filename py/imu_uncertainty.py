from __future__ import division
from dolfyn.adv import api as avm
import nortek_vector.bench_tests as btests
import ptools as pt
import numpy as np


def spec_epsilon(epsilon, tke, k, alpha=0.5, ):
    beta = np.pi * 3 / 5
    a = (np.sin(beta) / beta * tke) ** (5 / 2) / (alpha ** (3 / 2) * epsilon)
    b = a / (alpha * epsilon ** (2 / 3))
    return a / (1 + b * k ** (5 / 3))


flag = {}
flag['plot_spec'] = True
flag['plot_spec2'] = True
flag['plot_spec3'] = True

binner = avm.TurbBinner(9600, 32, )

var_accel_screen_level = 2e-4

mcfilts = {'unfiltered': 0.0,
           '30s': 1. / 30,
           '10s': 1. / 10,
           '5s': 1. / 5, }

if 'rd' not in vars():
    rd = btests.load('C')
    rd.add_data('vel', rd.pop_data('_u'), 'main')
    rd.add_data('corr', rd.pop_data('_corr'), 'main')
    rd.add_data('amp', rd.pop_data('_amp'), 'main')
    rd.props['rotate_vars'].remove('_u')
    rd.props['rotate_vars'].add('vel')
    rd.props['body2head_rotmat'] = np.eye(3)
    rd.props['body2head_vec'] = np.array([1, 1, 1]) * np.sqrt(1. / 3)

    if not hasattr(rd, 'AngRt'):
        rd.add_data('AngRt', rd.DAng * rd.fs, 'orient')
        rd.add_data('Accel', rd.DVel * rd.fs, 'orient')
        rd.del_data(['DVel', 'DAng'])


if 'bindatmc' not in vars():
    bindatmc = {}

    for filt_tag, filt_freq in mcfilts.iteritems():
        mcdat = rd.copy()
        avm.motion.correct_motion(mcdat, accel_filtfreq=filt_freq)

        bdmc = binner(mcdat)
        bdmc.add_data('Spec_urot', binner.psd(mcdat.velrot))
        bdmc.add_data('Spec_uacc', binner.psd(mcdat.velacc))
        bdmc.omega[-1] *= -1
        bindatmc[filt_tag] = bdmc


bd = binner(rd)
bd.add_data('Spec_Accel', binner.psd(rd.Accel))
bd.add_data('Spec_AngRt', binner.psd(rd.AngRt))

bd.add_data('Var_Accel', binner.std(rd.Accel) ** 2)
bd.add_data('Var_AngRt', binner.std(rd.AngRt) ** 2)

gd = bd.Var_Accel.sum(0) < var_accel_screen_level
gd &= bd.mpltime < pt.dt.date2num(pt.datetime(2015, 1, 29, 6, 0, 0))

# tkr = pt.dt.HourLocator(byhour=range(0, 24, 3))
# tkr2 = pt.dt.HourLocator()
# fmt = pt.dt.DateFormatter('%H')

colors = ['r', 'g', 'b']

if flag.get('plot_spec', False):

    fig, axs = pt.newfig(3, 2, 1, sharex=True,
                         gridspec_kw=dict(left=0.2, right=0.94, top=0.95))

    for idx in xrange(3):
        axs[0].loglog(bd.freq, bd.Spec_Accel[idx][gd].mean(0) * 2 * np.pi,
                      color=colors[idx], linestyle='-')
        axs[1].loglog(bd.freq, bd.Spec_AngRt[idx][gd].mean(0) * 2 * np.pi,
                      color=colors[idx], linestyle='-')

    #axs[1].set_xlim([1e-3, 20])
    #axs[0].set_ylim([1e-9, 1e-4])
    #axs[0].set_title('Spectra: Bench Test %s' % tag)
    axs[0].set_ylabel('$a\, \mathrm{[m^2/s^{4}/hz]}$', size='large')
    axs[1].set_ylabel('$\omega\,\mathrm{[rad^2/s^2/hz]}$', size='large')
    axs[1].set_xlabel('$f\, \mathrm{[hz]}$', size='large')

    #fig.saxes.hide('xticklabels', fig.saxes.ax[-1, :])

if flag.get('plot_spec2', False):

    fig, ax = pt.newfig(2, 1, 1, sharex=True,
                        gridspec_kw=dict(left=0.2, right=0.94, top=0.95))

    bdmc = bindatmc['unfiltered']
    for idx in xrange(3):
        ax.loglog(bd.freq, bdmc.Spec_uacc[idx][gd].mean(0) * 2 * np.pi,
                  color=colors[idx], linestyle='-')
        ax.loglog(bd.freq, bdmc.Spec_urot[idx][gd].mean(0) * 2 * np.pi,
                  color=colors[idx], linestyle='-')

    #axs[1].set_xlim([1e-3, 20])
    ax.set_ylim([1e-8, 1e1])
    #axs[0].set_title('Spectra: Bench Test %s' % tag)
    # axs[0].set_ylabel('$u_{acc}\, \mathrm{[m^2/s^{2}/hz]}$', size='large')
    # axs[1].set_ylabel('$u_{rot}\,\mathrm{[m^2/s^2/hz]}$', size='large')
    # axs[1].set_xlabel('$f\, \mathrm{[hz]}$', size='large')

    #fig.saxes.hide('xticklabels', fig.saxes.ax[-1, :])
    fig.savefig(pt.figdir + 'stationary_noise02.pdf')

if flag.get('plot_spec3', False):

    fig, ax = pt.newfig(1, 1, sharex=True, figsize=3,
                        gridspec_kw=dict(left=0.2, right=0.94, top=0.95, bottom=0.15))

    line_f = np.logspace(-3, 2, 50)
    epsilon, tke, U = (5e-4, 0.03, 1.0, )
    line_high = (spec_epsilon(epsilon=epsilon,
                              tke=tke,
                              k=line_f * 2 * np.pi * U) *
                 2 * np.pi * U)
    epsilon, tke, U = (8e-6, 0.0025, 0.5, )
    line_low = (spec_epsilon(epsilon=epsilon,
                             tke=tke,
                             k=line_f * 2 * np.pi * U) *
                2 * np.pi * U)
    # ax.plot(line_f,
    #         line_high,
    #         'k--', zorder=-5)
    # ax.plot(line_f,
    #         line_low,
    #         'k--', zorder=-5)
    ax.fill_between(line_f, line_high, line_low, facecolor='0.8',
                    edgecolor='k',
                    linewidth=0.5, zorder=-8, linestyle='--')

    colors_now = {'unfiltered': 'k',
                  '30s': 'm',
                  '10s': 'b',
                  '5s': 'g', }
    tmp = bdmc.Spec_urot[:, gd].max(0).mean(0) * 2 * np.pi
    ax.loglog(bd.freq, tmp,
              color='y', linestyle='-', label='$n_{rot}$')
    print ''
    print ('u_rot error level: {:0.3f} cm/s'
           .format((np.trapz(tmp, bdmc.freq) ** 0.5) * 100))
    for tag in ['unfiltered', '30s', '10s', '5s']:
        bdmc = bindatmc[tag]
        tmp = bdmc.Spec_uacc[:, gd].max(0).mean(0) * 2 * np.pi
        tmp2 = bdmc.Spec_uacc[:, gd].min(0).mean(0) * 2 * np.pi
        if tag == 'unfiltered':
            label = '$n_{acc}$'
        else:
            label = '$n_{acc}' + '({})$'.format(tag)
            ax.axvline(mcfilts[tag], color=colors_now[tag], linestyle=':')
        ax.loglog(bd.freq, tmp,
                  color=colors_now[tag], linestyle='-',
                  label=label)
        ax.loglog(bd.freq, tmp2,
                  color=colors_now[tag], linestyle='--' )
        print ('u_acc error level ({}): {:0.3f} cm/s'
               .format(tag, (np.trapz(tmp, bdmc.freq) ** 0.5) * 100))
    bdmc = bindatmc['30s']

    #axs[1].set_xlim([1e-3, 20])
    ax.set_ylim([1e-7, 1e1])
    ax.set_xlim([1e-3, 1e1])
    ax.legend(prop=dict(size='small'))
    ax.set_xlabel('$f\ [\mathrm{hz}]$')
    ax.set_ylabel('$\mathrm{[m^2s^{-2}/hz]}$')
    ax.axhline(2e-4, color='k', linestyle=':')
    #axs[0].set_title('Spectra: Bench Test %s' % tag)
    # axs[0].set_ylabel('$u_{acc}\, \mathrm{[m^2/s^{2}/hz]}$', size='large')
    # axs[1].set_ylabel('$u_{rot}\,\mathrm{[m^2/s^2/hz]}$', size='large')
    # axs[1].set_xlabel('$f\, \mathrm{[hz]}$', size='large')

    #fig.saxes.hide('xticklabels', fig.saxes.ax[-1, :])
    fig.savefig(pt.figdir + 'stationary_noise04.pdf')
