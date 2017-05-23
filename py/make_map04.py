#!/usr/local/bin/python
# -*- coding: utf-8 -*-
from mpl_toolkits.basemap import Basemap as bm
from scipy.signal import convolve2d as conv2
import matplotlib.colors as mplc
from matplotlib import cm
import numpy as np
import ptools as pt
import matplotlib.pyplot as plt
from matplotlib.patches import ArrowStyle

#
# From the ./metadata.htm file
proj = 'lcc'  # Lambert Conformal Conic
# proj='omerc'
# proj='cyl'
grid_llcrnr = [-123.043722, 47.978454]
grid_urcrnr = [-122.466557, 48.523023]
grid_center = [-120.833333, 47.]
lat_1 = 47.5
lat_2 = 48.733333
ft_m = (0.3048 ** -1)  # /1.000002

meas_points = {
    'TTM-1': [-122.68581, 48.15285],
    'TTM-2': [-122.68654, 48.15327],
    'SM': [-122.68623, 48.15277],
}

# map_data=bm(projection=proj,lon_0=grid_center[0],lat_0=grid_center[1],llcrnrlon=grid_llcrnr[0],llcrnrlat=grid_llcrnr[1],urcrnrlon=grid_urcrnr[0],urcrnrlat=grid_urcrnr[1],lat_1=lat_1,lat_2=lat_2)
map_data = bm(projection=proj,
              lon_0=grid_center[0], lat_0=grid_center[1],
              llcrnrlon=grid_center[0], llcrnrlat=grid_center[1],
              urcrnrlon=grid_urcrnr[0], urcrnrlat=grid_urcrnr[1],
              lat_1=lat_1, lat_2=lat_2)

# This is for interactive mode
# (I don't want to reload the data every time I modify+run the script)
if 'bdat' not in vars():
    bdat = {}
    dat = np.load('data/bathy/g1230485.npz')
    dparams = {}
    for nm in dat.files:
        if not nm == 'elev':
            dparams[nm] = dat[nm]
    dat = dat['elev'] / ft_m

    grid_spacing = dparams['cellsize'] / ft_m

    xdat = np.arange(0, dat.shape[1]) * grid_spacing + grid_spacing / 2 + (
        dparams['xllcorner'] - 1640416.666667) / ft_m + 400
    ydat = np.arange(dat.shape[0], 0, -1) * grid_spacing - \
        grid_spacing / 2 + (dparams['yllcorner']) / ft_m
    Xdat, Ydat = np.meshgrid(xdat, ydat)

    Lonsdat, Latsdat = map_data(Xdat, Ydat, inverse=True)
    # FOR SOME REASON THE POSITIONS DON'T LINE UP! ARGH!
    #  Using different llcrnr and urcrnr in the map_data projection gives
    #  different results.
    #  Why?
    #  Whatever the case, this must be the source of the misalignment.
    #  Basemap doesn't actually accept a different lon_0 and llcrnrlon
    #  (and lat_0/llcrnrlat).
    # The origin of the coordinate system is always llcrnrlon,llcrnrlat. This
    # is so stupid!

    npt = 5
    smth = np.cos(np.pi * np.arange(-0.5 + 1. / (npt + 1),
                                    0.5 - 1. / (npt + 1) + 0.1 / npt,
                                    1. / (npt + 1))
                  ) ** 2
    smth = smth[:, None] * smth[None, :]
    smth /= smth.sum()
    bdat['z'] = conv2(dat, smth, mode='same')
    bdat['x'] = xdat
    bdat['y'] = ydat
    # grid_llcrnr[1],grid_urcrnr[1]=grid_urcrnr[1],grid_llcrnr[1]

# dlon=(grid_urcrnr[0]-grid_llcrnr[0])/dat.shape[1]
# dlat=(grid_urcrnr[1]-grid_llcrnr[1])/dat.shape[0]
# grid_lons=np.arange(grid_llcrnr[0],grid_urcrnr[0],dlon)
# grid_lats=np.arange(grid_urcrnr[1],grid_llcrnr[1],-dlat)

llcrnr = [-122.80, 48.07]
# urcrnr=[-122.6778,48.1568]
urcrnr = [-122.60, 48.20]
center = [-122.6797, 48.1547]
llcrnr_xy = map_data(*llcrnr)
urcrnr_xy = map_data(*urcrnr)
center_xy = map_data(*center)


# ilon=(llcrnr[0]<grid_lons) & (grid_lons<urcrnr[0])
# ilat=(llcrnr[1]<grid_lats) & (grid_lats<urcrnr[1])
ix = np.nonzero((llcrnr_xy[0] < bdat['x']) & (bdat['x'] < urcrnr_xy[0]))[0]
iy = np.nonzero((llcrnr_xy[1] < bdat['y']) & (bdat['y'] < urcrnr_xy[1]))[0]
ix = slice(ix[0] - 20, ix[-1] + 20)
iy = slice(iy[0] - 20, iy[-1] + 20)
# lons=grid_lons[ilon]
# lats=grid_lats[ilat]
# Lons,Lats=np.meshgrid(lons,lats)

# map=bm(projection=proj,lon_0=center[0],lat_0=center[1],llcrnrlon=center[0],llcrnrlat=center[1],urcrnrlon=urcrnr[0],urcrnrlat=urcrnr[1],lat_1=lat_1,lat_2=lat_2)
# X,Y=map(Lons,Lats)
x0, y0 = map_data(center[0], center[1])

cmap = pt.truncate_colormap(plt.get_cmap('YlGnBu_r'), 0.0, 0.7)
fig, ax = pt.newfig(308, 1, 1, figsize=2.2)
fig.clf()
axrect = [.01, .08, .72, .9]
ax = plt.axes(axrect)

decim = 5

#map_data.drawmapboundary(fill_color='none',zorder=-200, ax=ax)
pcol = ax.pcolor(bdat['x'][ix][::decim] - x0, bdat['y'][iy][::decim] - y0,
                 bdat['z'][iy, ix][::decim, ::decim], cmap=cmap,
                 rasterized=True)
pcol.set_clim([-200, 0])

# cbar_ax = plt.axes([axrect[0] + axrect[2] + .03, .5, .04, .5])
cmap = pt.truncate_colormap(plt.get_cmap('copper_r'), 0.0, 0.66)
cmap.set_over(cmap(1.))
pcon = ax.contourf(bdat['x'][ix][::decim] - x0, bdat['y'][iy][::decim] - y0,
                   bdat['z'][iy, ix][::decim, ::decim], np.arange(0, 100, 10), cmap=cmap,
                   linewidth=1,
                   extend='max',
                   rasterized=True)
# cbar = plt.colorbar(pcol, cax=cbar_ax,
#                     ticks=np.arange(-200, 10, 50))
# pcon.set_clim([0, 90])

erange = [-200, 100]
cmap = mplc.LinearSegmentedColormap.from_list(
    'bathy_land',
    np.vstack((plt.get_cmap('YlGnBu_r')(np.linspace(0, 0.7, 256)),
               plt.get_cmap('copper_r')(np.linspace(0.0, 0.66, 128)), ))
)

pcol = ax.pcolor([0, 0], [0, 0], np.ones((2, 2)) * np.NaN, cmap=cmap,
                 rasterized=True)
pcol.set_clim([-200, 100])
cbar_ax = plt.axes([axrect[0] + axrect[2] + .03, .04, .04, .5])
cbar = plt.colorbar(pcol, cax=cbar_ax,
                    ticks=np.arange(-200, 110, 50))
cbar.mappable.set_rasterized(True)
cbar_ax.set_ylabel('Elevation [m]')


# clf()
dpths = np.arange(-100, 1, 10)
lw = [.6] * len(dpths)
lw[0] = 2

ax.set_xticks([])
ax.set_yticks([])

# ax.set_ylim([-660, 300])
# ax.set_xlim([-1000, 400])
ax.set_aspect('equal', 'box')

ax.set_xlim([-9100, 6100])
ax.set_ylim([-8700, 5000])

pt.drawmapscale(5, (0.1, 0.17), ax=ax)

# bbox_props = dict(boxstyle='darrow,pad=0.3', fc='magenta', ec='none', lw=2)
# t = ax.text(-800, -500, 'Tidal Flow',
#             bbox=bbox_props, zorder=10,
#             ha='center', va='center',
#             rotation=pt.principal_angle,
#             size='large')

point = meas_points['TTM-1']

xd, yd = map_data(*point)
ax.plot(xd - x0, yd - y0, 'o', ms=6, zorder=50, mfc='r', mec='none')

# I think this is not including declination.
# arr = 800 * np.exp(1j * 162.0 * np.pi / 180)
# http://www.ngdc.noaa.gov/geomag-web/
declin = 16.3
arr = 1800 * np.exp(1j * (162.0 - declin) * np.pi / 180)
# ax.arrow(xd - x0, yd - y0, arr.real, arr.imag,
#          head_width=4.0, linewidth=0.5,
#          zorder=100)
ax.annotate("", (xd - x0 + arr.real, yd - y0 + arr.imag), (xd - x0, yd - y0),
            arrowprops=dict(arrowstyle=ArrowStyle('->', head_width=0.3, head_length=0.4),
                            shrinkA=0, shrinkB=0, lw=1.5),
            zorder=200,
            )
# ax.text(xd - x0 + arr.real, yd - y0 + arr.imag, '$u$', ha='right')
arr2 = arr * 2. / 3 * np.exp(1j * np.pi / 2)
# ax.arrow(xd - x0, yd - y0, arr.real, arr.imag,
#          head_width=4.0, linewidth=0.5,
#          zorder=100)
ax.annotate("", (xd - x0 + arr2.real, yd - y0 + arr2.imag), (xd - x0, yd - y0),
            arrowprops=dict(arrowstyle=ArrowStyle('->', head_width=0.2, head_length=0.3),
                            shrinkA=0, shrinkB=0, lw=1.0),
            zorder=200,
            )

lontmp = np.arange(llcrnr[0] - 0.05, urcrnr[0] + .05, 0.01)
lattmp = np.arange(llcrnr[1] - 0.05, urcrnr[1] + .05, 0.01)

for lat in [48.1666666]:
    xd, yd = map_data(lontmp, lontmp * 0 + lat)
    ax.plot(xd - x0, yd - y0, 'k:', lw=0.5)
    ax.text(1.01, yd[-5] - y0, u"48\u00B010'N",
            transform=ax.get_yaxis_transform(),
            ha='left', va='center', clip_on=False)

for lon in [-122.6666666]:
    xd, yd = map_data(lattmp * 0 + lon, lattmp)
    ax.plot(xd - x0, yd - y0, 'k:', lw=0.5)
    ax.text(xd[5] - x0, -0.01, u"122\u00B040'W",
            transform=ax.get_xaxis_transform(),
            ha='center', va='top', clip_on=False)

#fig.savefig(pt.figdir + 'map04.png', dpi=300)
fig.savefig(pt.figdir + 'map04.pdf')
