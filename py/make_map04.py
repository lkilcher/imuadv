from mpl_toolkits.basemap import Basemap as bm
from scipy.signal import convolve2d as conv2
import matplotlib.colors as mplc
from matplotlib import cm
import numpy as np
import ptools as pt
from matplotlib import pyplot as plt


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

# map_data=bm(projection=proj,lon_0=grid_center[0],lat_0=grid_center[1],llcrnrlon=grid_llcrnr[0],llcrnrlat=grid_llcrnr[1],urcrnrlon=grid_urcrnr[0],urcrnrlat=grid_urcrnr[1],lat_1=lat_1,lat_2=lat_2)
map_data = bm(projection=proj,
              lon_0=grid_center[0], lat_0=grid_center[1],
              llcrnrlon=grid_center[0], llcrnrlat=grid_center[1],
              urcrnrlon=grid_urcrnr[0], urcrnrlat=grid_urcrnr[1],
              lat_1=lat_1, lat_2=lat_2)

if 'bdat' not in vars():
    bdat = {}
    dat = np.load('/Users/lkilcher/data/bathy/puget_sound/g1230485/g1230485.npz')
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
decim = 3
ix = slice(ix[0], ix[-1], decim)
iy = slice(iy[0], iy[-1], decim)
# lons=grid_lons[ilon]
# lats=grid_lats[ilat]
# Lons,Lats=np.meshgrid(lons,lats)

# erange = [-200, 100]
# cmap = mplc.LinearSegmentedColormap.from_list(
#     'bathy_land',
#     np.vstack((plt.get_cmap('YlGnBu_r')(np.linspace(0, 0.7, 256)),
#                plt.get_cmap('copper_r')(np.linspace(0.0, 0.66, 128)), ))
# )

# map=bm(projection=proj,lon_0=center[0],lat_0=center[1],llcrnrlon=center[0],llcrnrlat=center[1],urcrnrlon=urcrnr[0],urcrnrlat=urcrnr[1],lat_1=lat_1,lat_2=lat_2)
# X,Y=map(Lons,Lats)
x0, y0 = map_data(center[0], center[1])

cmap = pt.truncate_colormap(plt.get_cmap('YlGnBu_r'), 0.0, 0.7)
fig, ax = pt.newfig(308, 1, 1, figsize=2)
fig.clf()
axrect = [.1, .1, .65, .87]
ax = plt.axes(axrect)
ax.axis('equal')
# map.drawmapboundary(fill_color=[.7,1.,1.],zorder=-200)
pcol = ax.pcolor(bdat['x'][ix] - x0, bdat['y'][iy] - y0, bdat['z'][iy, ix], cmap=cmap,
                 rasterized=True)
cbar_ax = plt.axes([axrect[0] + axrect[2] + .03, axrect[1], .04, .5])
hcb = axrect[1] - 0.5
cbar = plt.colorbar(pcol, cax=cbar_ax,
                    ticks=np.arange(-200, 10, 50))
cbar.mappable.set_rasterized(True)
cbar_ax.set_ylabel('Elevation [m]')
pcol.set_clim([-200, 0])

# cbar_ax = plt.axes([axrect[0] + axrect[2] + .03, .5, .04, .5])
cmap = pt.truncate_colormap(plt.get_cmap('copper_r'), 0.0, 0.66)
# cm.copper_r
# cmap.set_over(cmap(1))
pcon = ax.contourf(bdat['x'][ix] - x0, bdat['y'][iy] - y0,
                   bdat['z'][iy, ix], np.arange(0, 100, 5), cmap=cmap,
                   extend='max',
                   rasterized=True)
# cbar = plt.colorbar(pcol, cax=cbar_ax,
#                     ticks=np.arange(-200, 10, 50))
# pcon.set_clim([0, 90])


# clf()
dpths = np.arange(-100, 1, 10)
lw = [.6] * len(dpths)
lw[0] = 2

pt.drawmapscale(5, (0.2, 0.2), ax=ax)

# fig.savefig(pt.figdir + 'map04.png', dpi=300)
# fig.savefig(pt.figdir + 'map04.pdf')
