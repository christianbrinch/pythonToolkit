"""

    Plot scripts, Christian Brinch (c) 2015

    This module holds scripts to plot data

"""

from pythonToolkit import standards

try:
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    from matplotlib.tri import Triangulation, TriAnalyzer, UniformTriRefiner
except:
    print 'Error: plots requires matplotlib'

try:
    import numpy as np
except:
    print 'Error: plots requires numpy'


def createFigure(xaxis,yaxis,xc,yc,header,title):
    plt.clf()
    fig = plt.figure(1,figsize=(8, 8), dpi=600)
    ax=fig.add_subplot(111, aspect='equal')
    if(title):
        ax.set_title(header['OBJECT'])
    ax.set_xlim(np.max( xaxis[xc] )-0.02,np.min( xaxis[xc] )-0.02)
    ax.set_ylim(np.min( yaxis[yc] )-0.3,np.max( yaxis[yc] )-0.3)
    if 0. in xaxis[:] and 0. in yaxis[:]:
        ax.set_xlabel('R. A. offset (arcsec)')
        ax.set_ylabel('DEC offset (arcsec)')
    else:
        hours=header['OBSRA']/360.*24
        minutes=np.abs(hours-int(hours))*60.
        degrees=header['OBSDEC']
        arcmin=(np.abs(degrees-int(degrees)))*(60)
        ax.set_xlabel('Right Ascension (sec) %dh %dm' % (hours, minutes ))
        ax.set_ylabel('Declination (arcsec) %d$^\circ$ %d\' ' % (degrees, arcmin))

    ax.minorticks_on()
    plt.tick_params(axis='both', which='both', width=0.4)
    plt.subplots_adjust(bottom=0.15,left=0.15)
    cmap = plt.cm.copper
    cmap.set_bad('white',1.)

    return ax,cmap,fig

def beam(ax, posx,posy,beamsize_x,beamsize_y,angle):
    from matplotlib.patches import Ellipse
    from matplotlib.collections import PatchCollection
    beam = Ellipse((posx,posy), beamsize_x/15.,beamsize_y, angle=angle+90)
    patches=[]
    patches.append(beam)
    p = PatchCollection(patches, facecolor='white') #, hatch='//') 
    ax.add_collection(p)
    ax.plot([posx+0.04,posx-10]  ,[posy+0.5,posy+0.5],color='white')
    ax.plot([posx+0.04,posx+0.04],[posy+0.5,posy-10.],color='white')
    print "Beam size: ", beamsize_x,"x",beamsize_y, " Angle: ", angle
    return ax




def contour(data, x, y, label=None, log=False):

    tri=Triangulation(x,y)

    plt.close('all')
    plt.figure()
    ax=plt.subplot(111)
    ax.minorticks_on()
    if(log):
        ax.set_xscale("log",nonposx='clip')
        ax.set_yscale("log",nonposy='clip')

    ax.set_xlim([min(x.min(),y.min()),max(x.max(),y.max())])
    ax.set_ylim([min(x.min(),y.min()),max(x.max(),y.max())])
    plt.xlabel('r [AU]')
    plt.ylabel('z [AU]')

    nmax=data.max()
    nmin=data.min()
    levels=np.logspace(np.log10(nmin),np.log10(nmax),num=12)

    plt.tricontourf(tri, data, levels, norm=colors.LogNorm(vmin=nmin, vmax=nmax))
    cbar=plt.colorbar(format='%.2e')
    cbar.set_label(label)

    CS=plt.tricontour(tri, data, levels, colors='black', linewidths=1.5)
    plt.clabel(CS, fontsize=8, inline=1)
    cbar.add_lines(CS)

    plt.triplot(tri, color='black', alpha=0.2)

    plt.show(block=False)
