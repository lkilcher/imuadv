from mpl_toolkits.basemap import Basemap as bm
from scipy.signal import convolve2d as conv2
from matplotlib import cm
import local_tools as ltbx
import numpy as np
import darrow
from matplotlib import pyplot as pt


if 'veldat' not in vars():
    veldat = {}
    for nm in ['T2b-top', 'T1b-top', 'T1-top']:
        veldat[nm] = ltbx.load(nm, 'mcpax_bin')

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

dat = np.load('/home/lkilcher/data/bathy/puget_sound/g1230485/g1230485.npz')
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
dat = conv2(dat, smth, mode='same')
# grid_llcrnr[1],grid_urcrnr[1]=grid_urcrnr[1],grid_llcrnr[1]

# dlon=(grid_urcrnr[0]-grid_llcrnr[0])/dat.shape[1]
# dlat=(grid_urcrnr[1]-grid_llcrnr[1])/dat.shape[0]
# grid_lons=np.arange(grid_llcrnr[0],grid_urcrnr[0],dlon)
# grid_lats=np.arange(grid_urcrnr[1],grid_llcrnr[1],-dlat)

llcrnr = [-122.6938, 48.146]
# urcrnr=[-122.6778,48.1568]
urcrnr = [-122.674, 48.1590]
center = [-122.6797, 48.1547]
llcrnr_xy = map_data(*llcrnr)
urcrnr_xy = map_data(*urcrnr)
center_xy = map_data(*center)


# ilon=(llcrnr[0]<grid_lons) & (grid_lons<urcrnr[0])
# ilat=(llcrnr[1]<grid_lats) & (grid_lats<urcrnr[1])
ix = (llcrnr_xy[0] < xdat) & (xdat < urcrnr_xy[0])
iy = (llcrnr_xy[1] < ydat) & (ydat < urcrnr_xy[1])
# lons=grid_lons[ilon]
# lats=grid_lats[ilat]
# Lons,Lats=np.meshgrid(lons,lats)


dnow = dat[iy][:, ix]

# map=bm(projection=proj,lon_0=center[0],lat_0=center[1],llcrnrlon=center[0],llcrnrlat=center[1],urcrnrlon=urcrnr[0],urcrnrlat=urcrnr[1],lat_1=lat_1,lat_2=lat_2)
# X,Y=map(Lons,Lats)
x0, y0 = map_data(center[0], center[1])

# veldat=

cmap = ltbx.truncate_colormap(pt.get_cmap('YlGnBu'), 0.3, 1.0)
fig = pt.figure(308, figsize=[8, 5])
pt.clf()
axrect = [.1, .1, .72, .87]
ax = pt.axes(axrect)
ax.axis('equal')
# map.drawmapboundary(fill_color=[.7,1.,1.],zorder=-200)
plt = ax.pcolor(xdat[ix] - x0, ydat[iy] - y0, -dnow, cmap=cmap,
                rasterized=True)
cbar_ax = pt.axes([axrect[0] + axrect[2] + .03, axrect[1], .04, .5])
cbar = pt.colorbar(plt, cax=cbar_ax,
                   ticks=np.arange(0, 100, 10))
cbar_ax.set_ylabel('Depth [m]')
plt.set_clim([0, 70])
#cbar.patch.set_rasterized(True)
plt = ax.contourf(xdat[ix] - x0, ydat[iy] - y0,
                  dnow, np.arange(0, 60, 1), cmap=cm.copper_r,
                  rasterized=True)
plt.set_clim([0, 90])


# clf()
dpths = np.arange(0, -100, -10)
lw = [.6] * len(dpths)
lw[0] = 2
ax.contour(xdat[ix] - x0, ydat[iy] - y0, dnow, dpths, colors=['k']
           * len(dpths), linestyles=['-'] * len(dpths), linewidths=lw)
ax.text(200, 200, 'Admiralty\nHead',
        multialignment='center', ha='center', size='large')
ax.text(200, 80, '(Fort Casey\nState Park)',
        size='medium', multialignment='center', ha='center')

# bbox_props=dict(boxstyle='rarrow,pad=0.3',fc='cyan',ec='none',lw=2)
# t=text(-800,-400,'Tidal
# Flow',bbox=bbox_props,zorder=10,ha='center',va='center',rotation=ltbx.principal_angle,size='large')


bbox_props = dict(boxstyle='darrow,pad=0.3', fc='cyan', ec='none', lw=2)
t = ax.text(-800, -500, 'Tidal Flow',
            bbox=bbox_props, zorder=10,
            ha='center', va='center',
            rotation=ltbx.principal_angle,
            size='large')

ax.set_ylim([-660, 300])
ax.set_xlim([-1000, 400])

kws = {'T1b-top': dict(xytext=(-5, -5),
                       ha='right',
                       va='top'
                       ),
       'T2b-top': dict(xytext=(5, 5),
                       ha='left',
                       va='bottom',
                       )
       }

for nm in ['T1b-top', 'T2b-top']:
    ll = veldat[nm].props['latlon'][::-1]
    xd, yd = map_data(*ll)
    ax.annotate('TTM ' + nm[1], (xd - x0, yd - y0),
                textcoords='offset points', color='k',
                bbox=dict(boxstyle='round,pad=0.2',
                          fc='yellow',
                          alpha=0.6),
                zorder=40,
                #arrowpros=dict(arrowstyle='->'
                **kws[nm]
                )
    ax.plot(xd - x0, yd - y0, 'y.', ms=20, zorder=50)

#fig.savefig(ltbx.figdir + 'map02.png', dpi=300)
#fig.savefig(ltbx.figdir + 'map02.pdf')
