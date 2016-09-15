"""

fits tools, Christian Brinch (c) 2015

This module contains classes and methods to read and plot FITS files

"""

from pythonToolkit import standards
from pythonToolkit import plots
import sys
try:
    from astropy.io import fits
    from astropy import wcs
except:
    print 'Error: fits requires astropy'
try:
    import numpy as np
except:
    print 'Error: Numpy not found'
try:
    import matplotlib.pyplot as plt
    import matplotlib.colors as cols
    from matplotlib import cm
except:
    print 'Error: Matplotlib not found'
try:
    from scipy import constants
    cc = constants.physical_constants['speed of light in vacuum'][0]
    pi = constants.pi
except:
    print 'Error: Scipy not found'
try:
    import cubehelix
except:
    print 'Error: cubehelix not fount'



class read():
    def __init__(self,fileName, coordSys='rel'):
        try:
            hdulist = fits.open(fileName)
        except:
            sys.exit("Error: FITS file not found")

        self.fileName = fileName
        self.header = hdulist[0].header
        self.image = hdulist[0].data
        hdulist.close()

        # Collapse cube to three dimensions
        if (np.ndim(self.image) == 4):
            self.image=self.image.sum(axis=0)

        # Define xaxis, yaxis, and velocity axis
        # In case of continuum images, collapse frequency axis
        cellsize = self.header['CDELT2'] * pi/180. * 1./constants.arcsecond
        if(coordSys == 'rel'):
          self.xunit = 'arcsec'
          self.xaxis = -(np.arange(self.header['NAXIS1'])
                         -self.header['CRPIX1']) * cellsize
          self.yaxis =  (np.arange(self.header['NAXIS2'])
                         -self.header['CRPIX2']) * cellsize
        elif(coordSys == 'abs'):
            self.xunit = 'seconds'
            w = wcs.WCS(self.header)
            self.xaxis,self.yaxis,dummy,dummy = w.all_pix2world(
                np.linspace(0,self.header['NAXIS1']-1,self.header['NAXIS1']),
                np.linspace(0,self.header['NAXIS2']-1,self.header['NAXIS2']),
                1, 1, 0
            )
        else:
            print "Unknown choice of coordinate system"
            return

        if(self.header['NAXIS3'] > 1):
            if(self.header['CTYPE3'] == 'VELO-LSR'):
                self.vaxis = (np.arange(self.header['NAXIS3'])
                              *self.header['CDELT3'])/1000.
            if(self.header['CTYPE3'] == 'FREQ'):
                self.vaxis = (np.arange(self.header['NAXIS3'])
                              *self.header['CDELT3']
                              *cc/self.header['CRVAL3'])/1000.
                self.vaxis = -self.vaxis - cc*(self.header['CRVAL3']-
                              self.header['RESTFRQ'])/(self.header['RESTFRQ']*
                              1000.)
        else:
            self.vaxis = np.array([])
            if (np.ndim(self.image) == 3):
                self.image = self.image.sum(axis=0)

    def printHeader(self):
        for i in self.header:
            if i != 'HISTORY':
                print i, self.header[i]


        key = [x for x in self.header['HISTORY'] if
               ('Beam' in x or 'restoration' in x ) and 'arcsec' in x]
        try:
            print key[0]
        except:
            print "Header has no beam information in HISTORY"

    def getRMS(self):
        if([x for x in self.header if 'RMS' in x]):
            print "Taking RMS from FITS header"
            rms = self.header['RMS']
        elif([x for x in self.header if 'TELESCOP' in x or 'OBSERVER' in x]):
            print "This seems to be data. Calculating RMS..."
            rms = np.std(self.image[0:self.header['NAXIS2']/4,
                                    0:self.header['NAXIS1']/4])
        else:
            print "This seems to be a model. Contours are 10% of peak value."
            rms = self.image.max()/10.

        print "RMS: ",rms
        return rms

    def _getIndices(self,d,dx,dy,dv=0):
        yc = np.where((self.yaxis > self.yaxis[self.header['NAXIS2']/2]-dy/2.) &
                      (self.yaxis < self.yaxis[self.header['NAXIS2']/2]+dy/2.))
        xc = np.where((self.xaxis > self.xaxis[self.header['NAXIS1']/2]-dx/2.) &
                      (self.xaxis < self.xaxis[self.header['NAXIS1']/2]+dx/2.))
        if(self.vaxis.any()):
            vc = np.where((self.vaxis > -dv) & (self.vaxis < dv))
        else:
            vc = 0
        return xc,yc,vc

    def _plotBeam(self,d,ax,xc,yc,black):
        if ('HISTORY' in self.header) and (not 'BMAJ' in self.header):
            key=[x for x in self.header['HISTORY'] if
                 'Beam' in x and 'arcsec' in x]
            if key:
                beaminfo = key[0].split()
                offset = float(beaminfo[3])/2. + float(beaminfo[3])/4.
                plots.beam(ax, np.min(self.xaxis[xc])+offset,
                           np.min(self.yaxis[yc])+offset, float(beaminfo[3]),
                           float(beaminfo[5]), float(beaminfo[9]))
            else:
                key=[x for x in self.header['HISTORY'] if
                     'restoration' in x and 'arcsec' in x]
                if key:
                    beaminfo = key[0].split()
                    offset = float(beaminfo[2])/2. + float(beaminfo[4])/4.
                    plots.beam(ax, np.min(self.xaxis[xc])+offset,
                               np.min(self.yaxis[yc])+offset,
                               float(beaminfo[2]), float(beaminfo[4]),
                               float(beaminfo[8]))
        elif ('BMAJ' in self.header) and ('BMIN' in self.header):
            bmaj = float(self.header['BMAJ'])*3600.
            bmin = float(self.header['BMIN'])*3600.
            dx = (np.max(self.xaxis[xc])-np.min(self.xaxis[xc]))/10.
            dy = (np.max(self.yaxis[yc])-np.min(self.yaxis[yc]))/10.
            xcor,ycor = np.min(self.xaxis[xc]) + dx, np.min(self.yaxis[yc]) + dy
            xcen,ycen = xcor - dx/2., ycor - dy/2.

            if(self.xunit == 'arcsec'):
                a = 1.
                angle = float(self.header['BPA'])
            else:
                a = 15.
                angle = 90.


            plots.beam(ax, xcor, ycor, xcen, ycen, bmaj/a, bmin,
                       angle, black)

        print "Beam size: ", bmaj,"x", bmin, " Angle: ", float(self.header['BPA'])+90


    def continuum(self, dx=1e30, dy=-1, noBeam=False, rms=-1, immax=-1,
                  solid=True, title=True, wedge=True, lightBackground=True):

        if(dy == -1):
            dy = dx

        if(self.xunit == 'seconds'):
            dy = dy/3600.
            dx = dx/3600.

        if(self.vaxis.any()):
            sys.exit(self.fileName + " is not a continuum image")

        if(rms<0):
            rms = self.getRMS()

        xc,yc,vc = self._getIndices(self,dx,dy)

        ax,cmap,fig = plots.createFigure(self.xaxis, self.yaxis, xc, yc,
                                         self.header, title)

        if(immax == -1):
            immax = np.nanmax(self.image)

        cx = cubehelix.cmap(start=0, rot=-0.5, reverse=lightBackground)
        if(solid):
            im=ax.imshow(self.image, cmap=cx, alpha=1.0,
                        interpolation='nearest', origin='lower',
                        extent=[np.max(self.xaxis), np.min(self.xaxis),
                        np.min(self.yaxis), np.max(self.yaxis)],
                        vmin=np.nanmin(3*rms), vmax=immax, aspect='auto')
            if(wedge):
                cbar = plt.colorbar(im)
                cbar.set_label('Intensity ('+self.header['BUNIT']+')')
        else:
            im=ax.imshow(self.image, cmap=cx, alpha=0.0,
                         interpolation='nearest', origin='lower',
                         extent=[np.max(self.xaxis), np.min(self.xaxis),
                         np.min(self.yaxis), np.max(self.yaxis)],
                         vmin=np.nanmin(3*rms), vmax=immax, aspect='auto')
            im=ax.contour(self.xaxis, self.yaxis, self.image, colors='black',
                          aspect='auto', interpolation='nearest',
                          linewidths=1.5, levels=np.arange(150)*3*rms+3*rms)

        if(not noBeam):
            self._plotBeam(self, ax, xc, yc, lightBackground)

        return ax, fig


    def moment(self, dx=1e30, dv=1e30, sysvel=0, nobeam=False, mom=0, rms=-1,
               overplot=False):

        if not self.vaxis.any():
            sys.exit(self.fileName + " is not a line image")

        if(rms<0):
            rms = self.getRMS()

        self.vaxis = self.vaxis-sysvel
        xc,yc,vc = self._getIndices(self, dx, dv)

        ax,cmap = plots.createFigure(self.xaxis, self.yaxis, xc, yc,
                                     self.header)

        moment0 = (self.image[vc[0]]*(np.abs(self.vaxis[1]
                                             -self.vaxis[0]))).sum(axis=0)
        moment1 = np.zeros((self.header['NAXIS2'], self.header['NAXIS1']),
                           float)
        for j in range(self.header['NAXIS1']):
            for i in range(self.header['NAXIS2']):
                if (moment0[i,j] > 3*rms):
                    moment1[i,j] = (self.image[vc[0],i,j]*
                                    self.vaxis[vc[0]]).sum()/moment0[i,j]
                else:
                    moment1[i,j] = np.nan

        if(mom == 1):
          moment = moment1
        else:
          moment = moment0

        # make custom colormap
        cdict = {'red':  ((0.0, 0.0, 0.0),
                          (1.0, 0.3, 0.3)),
                 'green':((0.0, 0.0, 0.0),
                          (1.0, 0.5, 0.5)),
                 'blue': ((0.0, 0.0, 0.0),
                          (1.0, 1.0, 1.0))
                }
        myBu = cols.LinearSegmentedColormap('myBu', cdict, 256)

        # add ramp alpha channel
        myBuT = cols.LinearSegmentedColormap('myBuT', cdict, 256)
        myBuT._init()
        alphas = np.linspace(0, 1, 256+3)
        myBuT._lut[:,-1] = alphas

        im = ax.imshow(moment, cmap=plt.cm.jet, alpha=1.0,
                       interpolation='nearest', origin='lower',
                       extent=[np.max(self.xaxis), np.min(self.xaxis),
                       np.min(self.yaxis), np.max(self.yaxis)],
                       vmin=np.nanmin(moment), vmax=np.nanmax(moment),
                       aspect='auto')

        plt.subplots_adjust(bottom=0.15, left=0.15)
        cbar = plt.colorbar(im)
        cbar.set_label('Intensity integrated velocity')
        ax.contour(self.xaxis, self.yaxis, moment0, colors='black',
                   levels=np.arange(25)*2*rms+3*rms)



    def veloContour(self, dx=1e30, dv=1e30, sysvel=0, nobeam=False, rms=-1,
                    overplot=False):

        if not self.vaxis.any():
            sys.exit(self.fileName + " is not a line image")

        if(rms<0):
            rms = self.getRMS()

        self.vaxis = self.vaxis-sysvel
        xc,yc,vc = self._getIndices(self, dx, dv)

        i = np.where( elf.vaxis<sysvel)
        moment = (self.image[i[0]]).sum(axis=0)

        rms = 0.2
        ax.contour(self.xaxis, self.yaxis, moment, colors='blue', alpha=0.8,
                   levels=np.arange(25)*2*rms+3*rms, linewidths=3)

        # make custom colormap
        cdict = {'red': ((0.0, 0.0, 0.0),
                         (1.0, 1.0, 1.0)),
                'green':((0.0, 0.0, 0.0),
                         (1.0, 0.5, 0.5)),
                'blue': ((0.0, 0.0, 0.0),
                         (1.0, 0.3, 0.3))
                }
        myBu = cols.LinearSegmentedColormap('myBu', cdict, 256)

        # add ramp alpha channel
        myBuT = cols.LinearSegmentedColormap('myBuT', cdict, 256)
        myBuT._init()
        alphas = np.linspace(0, 1, 256+3)
        myBuT._lut[:,-1] = alphas

        i = np.where(self.vaxis>sysvel)
        moment = (self.image[i[0]]).sum(axis=0)

        rms = 0.5
        ax.contour(self.xaxis, self.yaxis, moment, colors='red', alpha=0.8,
                   levels=np.arange(25)*2*rms+3*rms, linewidths=3)
