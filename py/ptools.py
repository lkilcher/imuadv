import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mplc
import numpy as np
import matplotlib.dates as dt
import os.path
from matplotlib.path import Path
from matplotlib.patches import BoxStyle
# import matplotlib.gridspec as gridspec

plt.ion()
datetime = dt.datetime.datetime
rcParams = mpl.rcParams


def twocol():
    plt.style.use(['./amspub.mplstyle',
                   './amspub_twocol.mplstyle'])


def onecol():
    plt.style.use('./amspub.mplstyle')

style = dict(
    onecol=lambda: plt.style.context('./amspub.mplstyle'),
    twocol=lambda: plt.style.context(['./amspub.mplstyle',
                                     './amspub_twocol.mplstyle']),
    stylegg=lambda: plt.style.context('ggplot'),
    classic=lambda: plt.style.context('classic'),
)

onecol()


figdir = os.path.abspath('../fig/') + '/'

vel_comps = ['u', 'v', 'w']


def powline(xin=None, factor=1e-3, pow=-5. / 3):
    if xin is None:
        xin = np.logspace(-6, 2, 50)
    return xin, factor * xin ** pow


class label(object):

    def __init__(self, prfx='', sufx='', comp=vel_comps):
        self.prfx = prfx
        self.sufx = sufx
        self.comp = comp

    def __str__(self, ):
        return '$' + self.prfx + r'\vec{' + self.comp[0] + '}' + self.sufx + '$'

    def __getitem__(self, ind):
        return '$' + self.prfx + self.comp[ind] + self.sufx + '$'

    @property
    def spec(self, ):
        return label(prfx=r'S\{', sufx=self.sufx + r'\}', comp=self.comp)


latex = dict(umeas=label(sufx=r'_\mathrm{m}'),
             ue=label(),
             uhead=label(sufx=r'_\mathrm{h}'))


def calcFigSize(n, ax=[1, 0], frm=[.5, .5], norm=False):
    """
    sz, vec = calcFigSize(n, ax, frame) calculates the width (or
    height) of a figure with *n* subplots.  Specify the width (height)
    of each subplot with *ax[0]*, the space between subplots with
    *ax[1]*, and the left / right (bottom / top) spacing with
    *frame[0]* / *frame[1]*.

    calcFigSize returns *sz*, a scalar, which is the width (or height)
    the figure should, and *vec*, which is the three element vector
    for input to saxes.

    See also: saxes, axes, calcAxesSize
    """
    if hasattr(n, '__iter__'):
        n = np.sum(n)
    sz = n * ax[0] + (n - 1) * ax[1] + frm[0] + frm[1]
    # This checks that it is not the default.
    if not (isinstance(norm, bool) and not norm):
        frm = np.array(frm)
        ax = np.array(ax)
        frm = frm / sz * norm
        ax = ax / sz * norm
        sz = norm
    v = np.array([frm[0], (sz - frm[1]), ax[1]]) / sz
    return sz, v


def newfig(fignum, nrows=1, ncols=1,
           figsize=None,  # axsize=None,
           **kwargs):
    gskw = kwargs.setdefault('gridspec_kw', dict())
    # Move these inputs to gskw
    for kw in ['left', 'right', 'bottom', 'top', 'hspace', 'wspace']:
        if kw in kwargs:
            gskw[kw] = kwargs.pop(kw)
        else:
            gskw.setdefault(kw, rcParams['figure.subplot.' + kw])
    # Set the figsize
    if figsize is None:
        figsize = [None, None]
    elif not isinstance(figsize, (tuple, list)):
        figsize = [None, figsize]
    if figsize[0] is None:
        figsize[0] = rcParams['figure.figsize'][0]
    if figsize[1] is None:
        figsize[1] = rcParams['figure.figsize'][1]
    # if axsize is None:
    #     axsize = [None, None]
    # elif not isinstance(axsize, (tuple, list)):
    #     # Assume it is the axheight
    #     axsize = [None, axsize]
    # # axsize overwrites figsize if it is specified.
    # if axsize[0] is not None:
    #     figsize[0] = calcFigSize(ncols, [axsize[0], gskw['wspace']],
    #                              [gskw['left'], gskw['right']])[0]
    # if axsize[1] is not None:
    #     figsize[1] = calcFigSize(nrows, [axsize[1], gskw['wspace']],
    #                              [gskw['left'], gskw['right']])[0]
    fig = plt.figure(fignum, figsize=figsize)
    kwargs['num'] = fig.number
    fig.clf()
    fig, axs = plt.subplots(nrows, ncols, **kwargs)
    return fig, axs


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = mplc.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


# degrees, from several advs, should check this with ADCPs, and JT's data.
principal_angle = - 18

t0 = dt.date2num(datetime(2014, 6, 16, 0, 0, 0))
time_range = {'Xw': [dt.date2num(datetime(2014, 6, 16, 20, 36, 0)),
                     dt.date2num(datetime(2014, 6, 17, 7, 12, 0)), ],
              'T1': [dt.date2num(datetime(2014, 6, 16, 21, 12, 0)),
                     dt.date2num(datetime(2014, 6, 17, 14, 42, 0)), ],
              'T2b': [dt.date2num(datetime(2014, 6, 18, 8, 15, 0)),
                      dt.date2num(datetime(2014, 6, 19, 5, 0, 0)), ],
              'T1b': [dt.date2num(datetime(2014, 6, 18, 8, 0, 0)),
                      dt.date2num(datetime(2014, 6, 19, 5, 12, 0)), ],
              }


class DArrow(BoxStyle.LArrow):

    """
    (Double) Arrow Box
    """

    def __init__(self, pad=0.3):
        self.pad = pad
        super(BoxStyle.DArrow, self).__init__(pad)

    def transmute(self, x0, y0, width, height, mutation_size):

        # padding
        pad = mutation_size * self.pad

        # width and height with padding added.
        width = width + 2. * pad
        width, height = width, height + 2. * pad

        # boundary of the padded box
        x0, y0 = x0 - pad, y0 - pad,
        x1, y1 = x0 + width, y0 + height

        dx = (y1 - y0) / 2.
        dxx = dx * .5
        # adjust x0.  1.4 < - sqrt(2)
        x0 = x0 + pad / 1.4

        cp = [(x0 + dxx, y0), (x1, y0),
              (x1, y0 - dxx), (x1 + dx + dxx, y0 + dx), (
                  x1, y1 + dxx),  # right arrow
              (x1, y1), (x0 + dxx, y1),
              (x0 + dxx, y1 + dxx), (x0 - dx, y0 + dx), (
                  x0 + dxx, y0 - dxx),  # left arrow
              (x0 + dxx, y0), (x0 + dxx, y0)]

        com = [Path.MOVETO, Path.LINETO,
               Path.LINETO, Path.LINETO, Path.LINETO,
               Path.LINETO, Path.LINETO,
               Path.LINETO, Path.LINETO, Path.LINETO,
               Path.LINETO, Path.CLOSEPOLY]

        path = Path(cp, com)

        return path

# NEED TO ADD THIS TO BoxStyle!
BoxStyle.DArrow = DArrow
BoxStyle._style_list['darrow'] = DArrow

latlons = {'T1b': (48.15256666, -122.68678333),
           'T2b': (48.152783333, -122.686316666),
           'T1': (48.1525, -122.6867),
           }

inst = {'Xw-bottom': 'NREL00',
        'Xw-center': 'APL02',
        'Xw-star': 'P01',
        'Xw-port': 'NREL03',
        'T1-bottom': 'NREL01',
        'T1-top': 'NREL02',
        'T1b-top': 'NREL02',
        'T1b-bottom': 'NREL01',
        'T2b-top': 'NREL03',
        'T2b-bottom': 'F01'
        }

filestart = {'Xw': 'Xwing',
             'T1': 'ttm01',
             'T2': 'ttm02',
             'T1b': 'ttm01b',
             'T2b': 'ttm02b'}
