"""

    Plot scripts, Christian Brinch (c) 2015

    This module holds scripts to plot data

"""

from pythonToolkit import standards

try:
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import matplotlib.colors as colors
except:
    print 'Error: plots requires matplotlib'

try:
    import numpy as np
except:
    print 'Error: plots requires numpy'


def createFigure(title='', logaxis=False):
    plt.clf()
    fig = plt.figure(1,figsize=(8, 8), dpi=600)
    ax=fig.add_subplot(111, aspect='equal')
    ax.set_title(title)

    if(logaxis):
        ax.set_xscale("log", nonposx='clip')
        ax.set_yscale("log", nonposx='clip')

    return ax, fig

def createSkyPlot(xaxis,yaxis,xc,yc,header,title,colors='plt.cm.copper'):
    ax, fig = createFigure(title)
    ax.set_xlim(np.max( xaxis[xc] ),np.min( xaxis[xc] ))
    ax.set_ylim(np.min( yaxis[yc] ),np.max( yaxis[yc] ))

    if (0. in xaxis[:] and 0. in yaxis[:]) or 'OBSRA' not in header:
        ax.set_xlabel('R. A. offset (arcsec)')
        ax.set_ylabel('DEC offset (arcsec)')
    else:
        ticks_x = ticker.FuncFormatter(lambda xaxis, pos: '{0:.2f}'.format( ( ((24.*xaxis/360.) - (24.*xaxis/360.).astype(int))*60. - (((24.*xaxis/360.) - (24.*xaxis/360.).astype(int))*60.).astype(int))*60. ))
        ax.xaxis.set_major_formatter(ticks_x)
        ticks_y = ticker.FuncFormatter(lambda yaxis, pos: '{0:.2f}'.format( ((yaxis-yaxis.astype(int))*60. - ((yaxis-yaxis.astype(int))*60.).astype(int))*60. ))
        ax.yaxis.set_major_formatter(ticks_y)
        hours=header['OBSRA']/360.*24
        minutes=np.abs(hours-int(hours))*60.
        degrees=header['OBSDEC']
        arcmin=(np.abs(degrees-int(degrees)))*60.
        ax.set_xlabel('Right Ascension (sec) %dh %dm' % (hours, minutes ))
        ax.set_ylabel('Declination (arcsec) %d$^\circ$ %d\' ' % (degrees, arcmin))

    ax.minorticks_on()
    plt.tick_params(axis='both', which='both', width=0.4)
    plt.subplots_adjust(bottom=0.15,left=0.15)
    cmap = plt.cm.copper
    cmap.set_bad('white',1.)

    return ax,cmap,fig


def beam(ax, corx, cory, posx,posy,beamsize_x,beamsize_y,angle,black):
    from matplotlib.patches import Ellipse
    from matplotlib.collections import PatchCollection
    if(black):
        col='black'
    else:
        col='white'
    beam = Ellipse((posx,posy), beamsize_x,beamsize_y, angle=angle+90)
    patches=[]
    patches.append(beam)
    p = PatchCollection(patches, facecolor=col) #, hatch='//')
    ax.add_collection(p)
    ax.plot([corx,corx-100.],[cory,cory],color=col)
    ax.plot([corx,corx],[cory,cory-100.],color=col)
    return ax
