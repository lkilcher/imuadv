from dolfyn.tools import psd
import load_binned as bd
import load_rawraw as rd
import load_raw as md
import dolfyn.adv.motion as avmot
import dolfyn.adv.rotate as avrot
import dolfyn.adv.turbulence as avturb
import scipy.signal as sig
import matplotlib.pyplot as plt
from scipy.io import loadmat
import numpy as np
import matplotlib.dates as mpld

plt.show()
plt.ion()
from collections import defaultdict
binner = avturb.TurbBinner(4800, 16)

syncdat = loadmat('/Users/lkilcher/data/ttm/sm2015/ADCP/BT_IMU_sync_data.mat')
syncdat['t_IMU'] = syncdat['t_IMU'][0]
syncdat['t_IMU'] -= 366
time_ind = np.argmin(np.abs((rd.smn.mpltime) - syncdat['t_IMU'][0]))

flag = defaultdict(lambda: False, {})
#flag['SMN'] = True
#flag['SMP'] = True
#flag['noise'] = True
flag['show BT'] = True

if flag['show BT']:

    plot_inds = slice(0, 10000)
    plot_inds2 = slice(plot_inds.start + time_ind, plot_inds.stop + time_ind)
    fignum = 5000
    fig = plt.figure(fignum, figsize=(9, 9))
    fig.clf()
    fig, axs = plt.subplots(3, 1, sharex=True, sharey=True, num=fignum,
                            gridspec_kw=dict(left=0.1, right=0.98, bottom=0.17, top=0.9))

    npole = 2
    fs = 16
    filt_freq = 1. / 5
    flt = sig.butter(npole, filt_freq / (fs / 2))
    btf = sig.filtfilt(flt[0], flt[1], syncdat['BT_INST_interp'])

    rnow = rd.smn.copy()
    avmot.correct_motion(rnow, accel_filtfreq=filt_freq)  # , vel_filtfreq=0.1)
    avrot.inst2earth(rnow, reverse=True)
    tmpd = rnow
    # tmpd = md.smn.copy()
    # avrot.earth2principal(tmpd, reverse=True)
    # avrot.inst2earth(tmpd, reverse=True)
    t0 = syncdat['t_IMU'][0]
    syncdat['time_sec'] = (syncdat['t_IMU'] - t0) * 24 * 3600
    tmpd['time_sec'] = (tmpd.mpltime - t0) * 24 * 3600
    for irow, ax in enumerate(axs):
        ax.plot(syncdat['time_sec'][plot_inds], syncdat['BT_INST_interp'][irow][plot_inds] - btf[irow][plot_inds], 'b')
        ax.plot(syncdat['time_sec'][plot_inds], btf[irow][plot_inds], 'k')
        ax.plot(syncdat['time_sec'][plot_inds], syncdat['IMU_INST'][irow][plot_inds], 'r')
        # ax.plot(tmpd.time_sec[plot_inds2],
        #         tmpd.urot[irow][plot_inds2] + tmpd.uacc[irow][plot_inds2], 'g')
        ax.plot(tmpd.time_sec[plot_inds2], tmpd.uacc[irow][plot_inds2], 'g')
        #ax.plot(tmpd.mpltime[plot_inds2] - t0, tmpd.urot[irow][plot_inds2], 'b')

    # loc = mpld.AutoDateLocator()
    # ax.xaxis.set_major_locator(loc)
    # ax.xaxis.set_major_formatter(mpld.AutoDateFormatter(loc))
    ax.set_xlim([80, 200])

if flag['SMN']:
    ####
    # The Stablemoor Nose figure

    fignum = 1002
    fig = plt.figure(fignum, figsize=(9, 4))
    fig.clf()
    fig, axs = plt.subplots(1, 3, sharex=True, sharey=True, num=fignum,
                            gridspec_kw=dict(left=0.1, right=0.98, bottom=0.17, top=0.9))
    dnow = bd.smn

    for irow, ax in enumerate(axs):
        inds = (1.5 < np.abs(dnow.u)) & (np.abs(dnow.u) < 2.0)
        ax.loglog(dnow.freq, dnow.Spec_uacc[
                  irow, inds].mean(0) / (2 * np.pi), 'b')
        ax.loglog(dnow.freq, dnow.Spec_urot[
                  irow, inds].mean(0) / (2 * np.pi), 'r')
        ax.loglog(dnow.freq, dnow.Spec_umot[
                  irow, inds].mean(0) / (2 * np.pi), 'm')
        ax.loglog(dnow.freq, dnow.Spec_uraw[
                  irow, inds].mean(0) / (2 * np.pi), 'k')
        ax.loglog(dnow.freq, dnow.Spec[irow, inds].mean(0) / (2 * np.pi), 'g')

    ax.set_xlim([1e-3, 1e1])
    ax.set_ylim([1e-6, 1e-1])
    axs[1].set_title('StableMoor Nose')

    rnow = rd.smn.copy()
    avmot.correct_motion(rnow, accel_filtfreq=0.05)  # , vel_filtfreq=0.1)
    dnow = binner(rnow)
    dnow.Spec_uacc = binner.psd(rnow.uacc)
    dnow.Spec_urot = binner.psd(rnow.urot)
    dnow.Spec_umot = binner.psd(rnow.urot + rnow.uacc)
    dnow.Spec_uraw = binner.psd(rnow.uraw)

    fignum = 2002
    fig = plt.figure(fignum, figsize=(9, 4))
    fig.clf()
    fig, axs = plt.subplots(1, 3, sharex=True, sharey=True, num=fignum,
                            gridspec_kw=dict(left=0.1, right=0.98, bottom=0.17, top=0.9))

    for irow, ax in enumerate(axs):
        inds = (1.5 < np.abs(dnow.u)) & (np.abs(dnow.u) < 2.0)
        ax.loglog(dnow.freq, dnow.Spec_uacc[
                  irow, inds].mean(0) / (2 * np.pi), 'b')
        ax.loglog(dnow.freq, dnow.Spec_urot[
                  irow, inds].mean(0) / (2 * np.pi), 'r')
        ax.loglog(dnow.freq, dnow.Spec_umot[
                  irow, inds].mean(0) / (2 * np.pi), 'm')
        ax.loglog(dnow.freq, dnow.Spec_uraw[
                  irow, inds].mean(0) / (2 * np.pi), 'k')
        ax.loglog(dnow.freq, dnow.Spec[irow, inds].mean(0) / (2 * np.pi), 'g')

    ax.set_xlim([1e-3, 1e1])
    ax.set_ylim([1e-6, 1e-1])
    axs[1].set_title('StableMoor Nose2')


if flag['SMP']:
    ####
    # The Stablemoor Nose figure

    fignum = 1003
    fig = plt.figure(fignum, figsize=(9, 4))
    fig.clf()
    fig, axs = plt.subplots(1, 3, sharex=True, sharey=True, num=fignum,
                            gridspec_kw=dict(left=0.1, right=0.98, bottom=0.17, top=0.9))
    dnow = bd.smp

    for irow, ax in enumerate(axs):
        inds = (1.5 < np.abs(dnow.u)) & (np.abs(dnow.u) < 2.0)
        ax.loglog(dnow.freq, dnow.Spec_uacc[
                  irow, inds].mean(0) / (2 * np.pi), 'b')
        ax.loglog(dnow.freq, dnow.Spec_urot[
                  irow, inds].mean(0) / (2 * np.pi), 'r')
        ax.loglog(dnow.freq, dnow.Spec_umot[
                  irow, inds].mean(0) / (2 * np.pi), 'm')
        ax.loglog(dnow.freq, dnow.Spec_uraw[
                  irow, inds].mean(0) / (2 * np.pi), 'k')
        ax.loglog(dnow.freq, dnow.Spec[irow, inds].mean(0) / (2 * np.pi), 'g')

    ax.set_xlim([1e-3, 1e1])
    ax.set_ylim([1e-6, 1e-1])
    axs[1].set_title('StableMoor Port')

if flag['noise']:
    randsig = np.random.randn(int(1e6), )
    nfft = 4800
    fs = 16
    filt_freq = 1. / 30
    npole = 2
    ps = psd.psd(randsig, nfft, fs)
    f = psd.psd_freq(nfft, fs)
    ff = psd.psd_freq(nfft, fs, full=True)

    fignum = 200
    fig = plt.figure(fignum, figsize=(4, 6))
    fig.clf()
    fig, axs = plt.subplots(2, 1, sharex=True, sharey=False, num=fignum,
                            gridspec_kw=dict(left=0.17, right=0.83, bottom=0.05, top=0.97))
    flt = sig.butter(npole, filt_freq / (fs / 2))
    #randsig_f = sig.lfilter(flt[0], flt[1], randsig)
    randsig_f = sig.filtfilt(flt[0], flt[1], randsig)
    ps_f = psd.psd(randsig_f, nfft, fs)

    ax = axs[0]
    ax.loglog(f, ps_f / ps)
    ax.axvline(filt_freq, linestyle='--', color='r')

    ax.set_ylim([1e-10, 1])
    # filt = sig.butter(2, float(filt_freq) / (samp_freq / 2))
    # for idx in range(hp.shape[0]):
    #     bd[idx] = bd[idx] - sig.filtfilt(filt[0], filt[1], bd[idx])

    ang = psd.delta_angle(randsig, randsig_f, nfft)
    ax = axs[1]
    axe = plt.twinx(ax)
    axe.set_ylim([0, 1])
    ax.semilogx(f, 180 / np.pi * np.angle(ang), 'b')
    axe.semilogx(f, np.abs(ang), 'r')
    [tl.set_color('r') for tl in axe.get_yticklabels()]
    [tl.set_color('b') for tl in ax.get_yticklabels()]

    b, a = sig.iirdesign(0.1, 0.3, 5, 50, ftype='butter')
    b, a = sig.butter(npole, filt_freq / (fs / 2))
    #b, a = sig.iirdesign(0.1, 0.3, 5, 50, )
    w, gd = sig.group_delay((b, a))

    fig = plt.figure(202)
    fig.clf()
    plt.title('Digital filter group delay')
    plt.plot(w, gd)
    plt.ylabel('Group delay [samples]')
    plt.xlabel('Frequency [rad/sample]')
    plt.show()
