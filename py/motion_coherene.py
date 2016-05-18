from scipy.io import loadmat
import scipy.signal as sig
import matplotlib.pyplot as plt
import ttm.sm2015 as data_api
import numpy as np
import dolfyn.adv.api as avm
import dolfyn.tools.psd as psd
import matplotlib.dates as dt
plt.ion()

flag = {}
#flag['bt_basic_time'] = True
#flag['bt_filt_time'] = True
flag['filt spec'] = True
flag['show cohere'] = True
flag['plot time'] = True

filt_freqs = {'unfilt': 0.0,
              '5s': 1. / 5,
              # '10s': 1. / 10,
              # '30s': 1. / 30
}

pii = np.pi * 2

# 4800 points is 5min at 16hz
binner = avm.TurbBinner(4800, 16)

if 'bt' not in vars():
    # This file is from Sam Harding.
    bt = loadmat(data_api.package_root + 'ADCP/BT_IMU_sync_data.mat')
    bt['t_IMU'] = bt['t_IMU'][0] - 366

if 'dat' not in vars():
    # Load the 'raw' (unprocessed) data that corresponds to the file above.
    dat = data_api.load('SM_Nose', coordsys='raw')

    offind = np.abs(bt['t_IMU'][0] - dat.mpltime).argmin()

    # Perform basic motion correction based on IMU data:
    bindat_filt = {}
    dat_filt = {}
    for filt_tag, filt_freq in filt_freqs.iteritems():
        datmc = dat.copy()
        avm.motion.correct_motion(datmc, accel_filtfreq=filt_freq)

        # The ADP data is in the instrument frame, so we make the ADV
        # data consistent.
        avm.rotate.inst2earth(datmc, reverse=True)

        # Crop the data so that it is consistent with the ADP data.
        datmc = datmc.subset(slice(offind, (offind + len(bt['t_IMU']))))

        # Calculate 'rotational' motion of ADP head.
        motcalc = avm.motion.CalcMotion(datmc)
        ur_adp = motcalc.calc_urot(np.array([-1, 0, 0]), to_earth=False)

        ubt = bt['BT_INST_interp']
        ubt2 = bt['BT_INST_interp'] - ur_adp
        if filt_freq > 0:
            filt = sig.butter(2, filt_freq / (dat.fs / 2))
            ubt = sig.filtfilt(filt[0], filt[1], ubt)
            ubt2 = sig.filtfilt(filt[0], filt[1], ubt2)
        datmc.add_data('ubt', ubt, 'orient')
        datmc.add_data('ubt2', ubt2, 'orient')
        if filt_freq > 0:
            datmc._u += datmc.ubt  # Add the bt to the vel
        datmc.rotate_vars.update({'ubt', 'ubt2'})
        avm.rotate.inst2earth(datmc)
        avm.rotate.earth2principal(datmc)
        dat_filt[filt_tag] = datmc

        datnow = datmc

        datbd = binner(datnow)
        datbd.Spec_uraw = binner.psd(datnow.uraw)
        datbd.Spec_uacc = binner.psd(datnow.uacc)
        datbd.Spec_urot = binner.psd(datnow.urot)
        datbd.Spec_ubt = binner.psd(datnow.ubt)
        if filt_freq > 0:
            datbd.Spec_umot = binner.psd(datnow.uacc +
                                         datnow.urot +
                                         datnow.ubt)
        else:
            datbd.Spec_umot = binner.psd(datnow.uacc +
                                         datnow.urot)
        bindat_filt[filt_tag] = datbd

    dat = dat.subset(slice(offind, (offind + len(bt['t_IMU']))))

# The mean velocity is totally contaminated by motion correction, so
# lets fix it here
bindat_filt['unfilt']['_u'] = bindat_filt['5s']['_u']


def offset_scale(dat, offset, scale, ):
    return (dat - offset) * scale


def within(dat, minval, maxval):
    return (minval < dat) & (dat < maxval)


dnow = dat_filt['unfilt']
bnow = bindat_filt['unfilt']
dnow.umot = dnow.urot + dnow.uacc

inds = within(np.abs(bnow.U), 1.0, 1.5)

#n_fft = binner.n_bin / 4
n_fft = binner.n_bin

ubt = dnow.ubt
ubt = dnow.ubt2

if flag.get('show cohere', False):

    cpsd = binner.cpsd(dnow.umot, ubt, n_fft=n_fft)
    sp_bt = np.abs(binner.cpsd(ubt, ubt, n_fft=n_fft))
    sp_mot = np.abs(binner.cpsd(dnow.umot, dnow.umot, n_fft=n_fft))

    cpsdr = binner.cpsd(dnow.urot, ubt, n_fft=n_fft)
    sp_rot = np.abs(binner.cpsd(dnow.urot, dnow.urot, n_fft=n_fft))

    cpsda = binner.cpsd(dnow.uacc, ubt, n_fft=n_fft)
    sp_acc = np.abs(binner.cpsd(dnow.uacc, dnow.uacc, n_fft=n_fft))

    # coh = ((np.abs(cpsd[:, inds]) ** 2).mean(1) /
    #        (sp_mot[:, inds].mean(1) * sp_bt[:, inds].mean(1)))
    coh = ((np.abs(cpsd[:, inds].mean(1)) ** 2) /
           (sp_mot[:, inds].mean(1) * sp_bt[:, inds].mean(1)))
    coha = ((np.abs(cpsda[:, inds].mean(1)) ** 2) /
            (sp_acc[:, inds].mean(1) * sp_bt[:, inds].mean(1)))
    cohr = ((np.abs(cpsdr[:, inds].mean(1)) ** 2) /
            (sp_rot[:, inds].mean(1) * sp_bt[:, inds].mean(1)))

    fig = plt.figure(300, figsize=[6, 9])
    fig.clf()
    fig, axs = plt.subplots(3, 1, num=fig.number,
                            gridspec_kw=dict(right=0.75,
                                             left=0.1,
                                             top=0.96,
                                             bottom=0.06,
                                             hspace=0.08),
                            sharex=True, sharey=True, )

    for iax, ax in enumerate(axs):
        ax.semilogx(bnow.freq, coh[iax],
                    'k', label='$u_{mot}$')
        ax.semilogx(bnow.freq, cohr[iax],
                    'm', label='$u_{rot}$')
        ax.semilogx(bnow.freq, coha[iax],
                    'b', label='$u_{acc}$')
        ax.text(0.02, 0.98, ['u', 'v', 'w'][iax] + '-component',
                transform=ax.transAxes,
                ha='left', va='top')

    axs[0].set_title('Coherence with $u_{bt}$ (%s coord sys)' %
                     dnow.props['coord_sys'], size='large')
    axs[0].legend(loc='upper left', bbox_to_anchor=[1.02, 1])
    ax.set_ylim([0, 1])
    ax.set_xlim([1e-3, 1])
    ax.set_xlabel('freq [hz]')
    fig.savefig('../fig/BT_IMU_Coherence01.pdf')

if flag.get('filt spec', False):

    fig = plt.figure(400, figsize=[6, 9])
    fig.clf()
    fig, axs = plt.subplots(3, 1, num=fig.number,
                            gridspec_kw=dict(right=0.75,
                                             left=0.13,
                                             top=0.96,
                                             bottom=0.06,
                                             hspace=0.08),
                            sharex=True, sharey=True)
    #vars = ['Spec', 'Spec_uraw', 'Spec_uacc', 'Spec_urot', 'Spec_ubt', 'Spec_umot']

    for iax, ax in enumerate(axs):
        ax.loglog(bnow.freq,
                  (bnow.Spec[iax][inds].mean(0) - 2e-5) * pii,
                  'g', label='$u$', linewidth=2, zorder=10)
        ax.loglog(bnow.freq,
                  (bnow.Spec_uraw[iax][inds].mean(0) - 2e-5) * pii,
                  'y', label='$u_{raw}$')
        ax.loglog(bnow.freq,
                  bnow.Spec_umot[iax][inds].mean(0) * pii,
                  'k', label='$u_{mot}$', zorder=8, linewidth=1.5)
        ax.loglog(bnow.freq,
                  bnow.Spec_urot[iax][inds].mean(0) * pii,
                  'm', label='$u_{rot}$')
        ax.loglog(bnow.freq,
                  bnow.Spec_uacc[iax][inds].mean(0) * pii,
                  'b', label='$u_{acc}$')
        ax.loglog(bnow.freq,
                  bnow.Spec_ubt[iax][inds].mean(0) * pii,
                  'r', label='$u_{bt}$')
        ax.text(0.98, 0.98, ['u', 'v', 'w'][iax],
                ha='right', va='top', transform=ax.transAxes)

    ax.set_ylabel('$\mathrm{[m^2/s^2/hz]}$')
    ax.set_ylim([1e-5, 1e1])
    ax.set_xlabel('freq [hz]')
    axs[0].legend(loc='upper left', bbox_to_anchor=[1.02, 1])
    axs[0].set_title('Velocity Spectra (%s coord sys)' %
                     dnow.props['coord_sys'], size='large')

    fig.savefig('../fig/NoMC_Spectra01.pdf')


if flag.get('plot time', False):

    tkr = dt.MinuteLocator(range(0, 60, 10))
    #tkr = dt.HourLocator(range(0, 24, 1))
    tk_fmt = dt.DateFormatter('%H:%M')
    fig = plt.figure(100)
    fig.clf()
    fig, axs = plt.subplots(3, 1, num=fig.number, sharex=True, sharey=True)
    #vars = ['Spec', 'Spec_uraw', 'Spec_uacc', 'Spec_urot', 'Spec_ubt', 'Spec_umot']

    time_inds = slice(10000, 40000)

    for iax, ax in enumerate(axs):
        ax.plot(dnow.mpltime[time_inds], dnow._u[iax][time_inds], 'k-')
        ax.plot(dnow.mpltime[time_inds], dnow.umot[iax][time_inds], 'b-')
        ax.set_ylabel('${}'.format(['u', 'v', 'w'][iax] + '\ \mathrm{[m/s]}$'))

    ax.xaxis.set_major_locator(tkr)
    ax.xaxis.set_major_formatter(tk_fmt)

    axs[0].set_title('Example velocity timeseries No MC (inst coord. sys.)')

    fig.savefig('../fig/NoMC_TimeSeries01.pdf')
    
