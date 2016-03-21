import ttm.sm2015 as sm15
import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec
from mpl_toolkits import axes_grid as axg
import numpy as np
plt.show()
plt.ion()

ranges = np.linspace(0, 2., 5)

dat = sm15.load('TTM', coordsys='pax', bin=True)

fig = plt.figure(101, figsize=(10, 8))
fig.clf()
axs = axg.Grid(fig, [.2, .1, .75, .85], (3, len(ranges) - 1), axes_pad=0.3, share_all=True)
#fig, axs = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(4, 8), )

for irow, axes in enumerate(axs.axes_row):
    for icol, ax in enumerate(axes):
        inds = (ranges[icol] < np.abs(dat.u)) & (np.abs(dat.u) < ranges[icol + 1])
        ax.loglog(dat.freq, dat.Spec_umot[irow, inds].mean(0) / (2 * np.pi), 'r')
        ax.loglog(dat.freq, dat.Spec_uraw[irow, inds].mean(0) / (2 * np.pi), 'b')
        ax.loglog(dat.freq, dat.Spec[irow, inds].mean(0) / (2 * np.pi), 'g')

ax.set_ylim([1e-6, 1e-1])
ax.set_xlim([1e-2, 2])
