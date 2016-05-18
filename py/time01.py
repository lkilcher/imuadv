#import load_binned as bd
import load_raw as rd
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
plt.show()
plt.ion()

flag = defaultdict(lambda: False, {})

figs = {}

figs['timeHR01'] = fig = plt.figure(2001, figsize=(9, 9))
fig.clf()
fig, axs = plt.subplots(3, 1, sharex=True, sharey=False, num=fig.number,
                        gridspec_kw=dict(left=0.10, right=0.95, bottom=0.1, top=0.95))

dnow = rd.smn

#rr = 

inds = slice(int(2e5), int(2.1e5))
t0 = dnow.mpltime[inds.start]

# Green = Black + Red + Blue
for idx, ax in enumerate(axs):
    ax.plot((dnow.mpltime[inds] - t0) * 24 * 3600, dnow.uraw[idx][inds], 'k-')
    ax.plot((dnow.mpltime[inds] - t0) * 24 * 3600, dnow._u[idx][inds], 'g-')
    ax.plot((dnow.mpltime[inds] - t0) * 24 * 3600, dnow.urot[idx][inds], 'r-')
    ax.plot((dnow.mpltime[inds] - t0) * 24 * 3600, dnow.uacc[idx][inds], 'b-')
    ax.axhline(0, color='k', linestyle=':')
