"""

fits tools, Christian Brinch (c) 2015

This module contains classes and methods to read and plot FITS files

"""

from pythonToolkit import standards
from pythonToolkit import plots
from pythonToolkit import analysis
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
    print 'Error: cubehelix not found'



class image():
    def __init__(self,fileName=None, header=None, axes=None, image=None,
                 coordSys='relative'):
        if fileName is not None:
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

        else:
            self.header = header
            self.image = image

        self.ax = axes
        self._setAxis(coordSys)









    def _setAxis(self, coordSys):
        # Define xaxis, yaxis, and velocity axis
        # In case of continuum images, collapse frequency axis
        cellsize = self.header['CDELT2'] * pi/180. * 1./constants.arcsecond
        if(coordSys == 'relative'):
          self.xunit = 'arcsec'
          self.xaxis = -(np.arange(self.header['NAXIS1'])
                         -self.header['CRPIX1']) * cellsize
          self.yaxis =  (np.arange(self.header['NAXIS2'])
                         -self.header['CRPIX2']) * cellsize
        elif(coordSys == 'absolute'):
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
            if(int(self.vaxis[0])==0):
                self.vaxis-=self.vaxis[-1]/2.




    def printHeader(self):
        for i in self.header:
            if i != 'HISTORY':
                print i, self.header[i]

        try:
            key = [x for x in self.header['HISTORY'] if
                ('Beam' in x or 'restoration' in x ) and 'arcsec' in x]
        except:
            print "Header has no beam information in HISTORY"



    def getRMS(self):
        if([x for x in self.header if 'RMS' in x]):
            print "Taking RMS from FITS header"
            rms = self.header['RMS']
            print "RMS: ",rms
        elif([x for x in self.header if 'TELESCOP' in x or 'OBSERVER' in x]):
            print "This seems to be data. Calculating RMS..."
            rms = np.std(self.image[0:self.header['NAXIS2']/4,
                                    0:self.header['NAXIS1']/4])
            print "RMS: ",rms
        else:
            print "This seems to be a model. Contours are 10% of peak value."
            rms = self.image.max()/10.
            print "First contour at: ",rms

        return rms



    def continuumSubtract(self):
        if not self.vaxis.any():
            sys.exit(self.fileName + " is not a line image")

        xc,yc,vc = self._getIndices(self, 1e30, 1e30, 1e30)

        tmpimage = np.zeros((len(vc),len(yc), len(xc)))
        for i in range(len(vc)):
            tmpimage[i,:,:]=self.image[i,:,:]-self.image[-1,:,:]

        return image(None, self.header, self.ax, tmpimage)



    def _getIndices(self,d,dx,dy,dv=0,sysvel=0, offset=[0,0] ):
        yc = np.where((self.yaxis > self.yaxis[self.header['NAXIS2']/2]-dy/2.+offset[1]) &
                      (self.yaxis < self.yaxis[self.header['NAXIS2']/2]+dy/2.+offset[1]))
        xc = np.where((self.xaxis > self.xaxis[self.header['NAXIS1']/2]-dx/2.+offset[0]) &
                      (self.xaxis < self.xaxis[self.header['NAXIS1']/2]+dx/2.+offset[0]))
        if(dv>0):
            vc = np.where((self.vaxis-sysvel > -dv) & (self.vaxis-sysvel < dv))
            return xc[0],yc[0],vc[0]
        else:
            return xc[0],yc[0],0


    def _getBeam(self,d,ax,xc,yc,black,silent=False):
        if ('HISTORY' in self.header) and (not 'BMAJ' in self.header):
            key=[x for x in self.header['HISTORY'] if
                 'Beam' in x and 'arcsec' in x]
            if key:
                beaminfo = key[0].split()
                offset = float(beaminfo[3])/2. + float(beaminfo[3])/4.
                #plots.beam(ax, np.min(self.xaxis[xc])+offset,
                #           np.min(self.yaxis[yc])+offset, float(beaminfo[3]),
                #           float(beaminfo[5]), float(beaminfo[9]))
                return ax, np.min(self.xaxis[xc])+offset, np.min(self.yaxis[yc])+offset, float(beaminfo[3]),float(beaminfo[5]), float(beaminfo[9])
            else:
                key=[x for x in self.header['HISTORY'] if
                     'restoration' in x and 'arcsec' in x]
                if key:
                    beaminfo = key[0].split()
                    offset = float(beaminfo[2])/2. + float(beaminfo[4])/4.
                    #plots.beam(ax, np.min(self.xaxis[xc])+offset,
                    #           np.min(self.yaxis[yc])+offset,
                    #           float(beaminfo[2]), float(beaminfo[4]),
                    #           float(beaminfo[8]))
                    return ax, np.min(self.xaxis[xc])+offset,np.min(self.yaxis[yc])+offset,float(beaminfo[2]), float(beaminfo[4]),float(beaminfo[8])
        elif ('BMAJ' in self.header) and ('BMIN' in self.header):
            bmaj = float(self.header['BMAJ'])
            bmin = float(self.header['BMIN'])
            angle = float(self.header['BPA'])
            if not silent:
                print "Beam size: ", bmaj*3600.,"x", bmin*3600., " Angle: ", angle+90
            dx = (np.max(self.xaxis[xc])-np.min(self.xaxis[xc]))/10.
            dy = (np.max(self.yaxis[yc])-np.min(self.yaxis[yc]))/10.
            xcor,ycor = np.min(self.xaxis[xc]) + dx, np.min(self.yaxis[yc]) + dy
            xcen,ycen = xcor - dx/2., ycor - dy/2.
            if(self.xunit == 'arcsec'):
                bmaj=bmaj*3600.
                bmin=bmin*3600.

            return ax, xcor, ycor, xcen, ycen, bmaj, bmin, angle, black
        else:
            print "cannot find any beam info"



    def _initPlot(self,dx,dy,dv,sysvel,rms,overplot,offset):

        if dx <= 0.:
            sys.exit("Delta_x cannot be zero or less")

        if(dy == -1):
            dy = dx

        if(self.xunit == 'seconds'):
            dy = dy/3600.
            dx = dx/3600.

        if(rms<0):
            rms = self.getRMS()

        if(overplot is True):
            alpha=0.4
        else:
            alpha=1.

        xc,yc,vc = self._getIndices(self, dx, dy, dv, sysvel, offset)

        return dy,rms,xc,yc,vc,alpha









    def continuum(self, dx=1e30, dy=-1, offset=[0,0], rms=-1, immax=None,
                  noBeam=False, wedge=True, solid=True, title='',
                  lightBackground=True, overplot=False):
        ''' Plots a map of the continuum emission

            Keyword arguments:
            dx -- delta x in arcseconds (default is entire x axis)
            dy -- delta y in arcseconds (default is same as dx)
            offset -- center offset of the map in arcseconds (default is [0,0])

            rms -- noise level (default is getRMS())
            immax -- Saturation level (default is image maximum)

            noBeam -- if true, do not plot beam (default is false)
            wedge -- if true, plot color bar (default is true)
            solid -- if true, plot filled contours (default is true)
            title -- title string for the plot (default is none)
            lightBackground -- if false, invert colors (default is true)

            overplot -- do not erase previous plot (default is false)
        '''

        dy, rms, xc, yc, vc, alpha = self._initPlot(dx, dy, 0., 0., rms,
                                                    overplot, offset)

        if overplot is not True or self.ax is None:
            self.ax, cmap = plots.createSkyPlot(self.xaxis, self.yaxis, xc, yc,
                                                self.header, title,
                                                cubehelix.cmap(start = 0,
                                                rot = -0.5,
                                                reverse = lightBackground))

        if immax is None:
            immax = np.nanmax(self.image[0,:,:])

        if not solid:
            alpha = 0.0

        im=self.ax.imshow(self.image[0,:,:], cmap=cmap, alpha=alpha,
                          interpolation='nearest', origin='lower',
                          extent=[np.max(self.xaxis-offset[0]),
                                  np.min(self.xaxis-offset[0]),
                                  np.min(self.yaxis-offset[1]),
                                  np.max(self.yaxis-offset[1])],
                          vmin=np.nanmin(3*rms), vmax=immax, aspect='auto')

        if not solid:
            self.ax.contour(self.xaxis, self.yaxis, self.image[0,:,:],
                            colors='black', aspect='auto',
                            interpolation='nearest', lw=1.5,
                            levels=np.arange(150)*3*rms+3*rms)

        if solid and wedge and overplot is not True:
            cbar = plt.colorbar(im)
            cbar.set_label('Intensity ('+self.header['BUNIT']+')')

        if(not noBeam):
            plots.beam(*self._getBeam(self, self.ax, xc, yc, lightBackground))



    def polarization(self, dx=1e30, dy=-1, noBeam=False, rms=-1, immax=-1,
                     title='', wedge=True, lightBackground=True):

        if(dy == -1):
            dy = dx

        if(self.xunit == 'seconds'):
            dy = dy/3600.
            dx = dx/3600.

        if(rms<0):
            rms = self.getRMS()

        xc,yc,vc = self._getIndices(self,dx,dy)

        ax,cmap,fig = plots.createSkyPlot(self.xaxis, self.yaxis, xc, yc,
                                          self.header, title)

        tmpimage = np.sqrt(self.image[1]**2+self.image[2]**2)/self.image[0]

        if(immax == -1):
            immax = np.nanmax(tmpimage)

        cx = cubehelix.cmap(start=0, rot=-0.5, reverse=lightBackground)

        im=ax.imshow(tmpimage, cmap=cx, alpha=1.0,
                     interpolation='nearest', origin='lower',
                     extent=[np.max(self.xaxis), np.min(self.xaxis),
                     np.min(self.yaxis), np.max(self.yaxis)],
                     vmin=np.nanmin(3*rms), vmax=immax, aspect='auto')

        X, Y = np.meshgrid(self.yaxis,self.xaxis, sparse=False)
        skip=(slice(None,None,25),slice(None,None,25))
        tmp2image=np.sqrt(self.image[1]**2+self.image[2]**2)

        plt.quiver(X[skip], Y[skip], self.image[2][skip]/tmp2image[skip],
                   self.image[1][skip]/tmp2image[skip], self.image[0][skip],
                   alpha=.5)
        plt.quiver(X[skip], Y[skip], self.image[2][skip]/tmp2image[skip],
                   self.image[1][skip]/tmp2image[skip], edgecolor='k',
                   facecolor='None', linewidth=.5)
        if(wedge):
            cbar = plt.colorbar(im)
            cbar.set_label('Polarization degree')

        if(not noBeam):
            plots.beam(*self._getBeam(self, self.ax, xc, yc, lightBackground))




    def moment(self, mom=0, dx=1e30, dy=-1, dv=1e30, offset=[0,0], sysvel=0,
               rms=-1, immax=-1, nobeam=False, wedge=True, solid=True,
               title='', overplot=False ):
        ''' Plots a moment map of the line emission

            Keyword arguments:
            mom -- the moment to be displayed (default is 0th moment)

            dx -- delta x in arcseconds (default is entire x axis)
            dy -- delta y in arcseconds (default is same as dx)
            dv -- delta v in km/s (default is entire v axis)
            offset -- center offset of the map in arcseconds (default is [0,0])
            sysvel -- velocity offset in km/s (default is 0)

            rms -- noise level (default is getRMS())
            immax -- Saturation level (default is image maximum)

            noBeam -- if true, do not plot beam (default is false)
            wedge -- if true, plot color bar (default is true)
            solid -- if true, plot filled contours (default is true)
            title -- title string for the plot (default is none)

            overplot -- do not erase previous plot (default is false)
        '''

        if not self.vaxis.any():
            sys.exit(self.fileName + " is not a line image")

        dy, rms, xc, yc, vc, alpha = self._initPlot(dx, dy, dv, sysvel, rms,
                                                    overplot, offset)

        if overplot is not True or self.ax is None:
            self.ax, cmap = plots.createSkyPlot(self.xaxis, self.yaxis, xc, yc,
                                                self.header, title, plt.cm.bwr)

        cmap.set_bad('green',0.03)

        moment0 = (self.image[vc]*(np.abs(self.vaxis[1] -
                   self.vaxis[0]))).sum(axis=0)

        if mom == 1:
            moment1 = np.zeros((self.header['NAXIS2'],self.header['NAXIS1']),
                      float)
            for j in range(self.header['NAXIS1']):
                for i in range(self.header['NAXIS2']):
                    if moment0[i,j] > 3*rms:
                        moment1[i,j] = (self.image[vc,i,j]*(self.vaxis[vc] -
                                        sysvel)).sum()/moment0[i,j]
                    else:
                        moment1[i,j] = np.nan

            moment = moment1
        else:
            moment = moment0

        if(immax == -1):
            immax = np.maximum(np.nanmax(moment), abs(np.nanmin(moment)))

        if mom == 1:
            immin = -immax
        else:
            immin = 0.

        if solid:
            im = self.ax.imshow(moment, cmap=cmap, alpha=alpha,
                                interpolation='nearest', origin='lower',
                                extent=[np.max(self.xaxis-offset[0]),
                                        np.min(self.xaxis-offset[0]),
                                        np.min(self.yaxis-offset[1]),
                                        np.max(self.yaxis-offset[1])],
                                vmin=immin, vmax=immax, aspect='auto')

        self.ax.contour(self.xaxis-offset[0], self.yaxis-offset[1], moment,
                        colors='black', levels=np.arange(500)*2*rms+3*rms)

        if solid and wedge and overplot is not True:
            cbar = plt.colorbar(im)
            cbar.set_label('Intensity integrated velocity')

        if(not nobeam):
            plots.beam(*self._getBeam(self, self.ax, xc, yc, True))


    def pv(self, dx=1e30, dy=-1, dv=1e30, offset=[0,0], sysvel=0., rms=-1,
           immax=-1, wedge=True, solid=True, title='', overplot=False):
        ''' Plots a PV-diagram

            Keyword arguments:
            dx -- delta x in arcseconds (default is 0)
            dy -- delta y in arcseconds (default is same as dx)
            dv -- delta v in km/s (default is entire v axis)
            offset -- center offset of the map in arcseconds (default is [0,0])
            sysvel -- velocity offset in km/s (default is 0)

            rms -- noise level (default is getRMS())

            wedge -- if true, plot color bar (default is true)
            solid -- if true, plot filled contours (default is true)
            title -- title string for the plot (default is none)

            overplot -- do not erase previous plot (default is false)
        '''
        if not self.vaxis.any():
            sys.exit(self.fileName + " is not a line image")

        if dy == -1:
            dy=4*abs(self.yaxis[1]-self.yaxis[0])

        dy, rms, xc, yc, vc, alpha = self._initPlot(dx, dy, dv, sysvel, rms,
                                                    overplot, offset)

        if overplot is not True or self.ax is None:
            self.ax, cmap = plots.createPvPlot(self.xaxis, self.vaxis-sysvel, xc, vc,
                                                self.header, title,
                                                plt.cm.viridis)


        pv = []
        pverr = []
        for i in vc:
            popt = analysis.twoD_gaussFit(self,self.image[i,:,:])

            pv.append(popt[1])
            pverr.append(popt[2])

        image = self.image[:,yc,:].sum(axis=1)

        if(immax == -1):
            immax = np.maximum(np.nanmax(image), abs(np.nanmin(image)))

        immin = 0


        if solid:
            im = self.ax.imshow(image, cmap=cmap, alpha=alpha,
                                interpolation='nearest', origin='lower',
                                extent=[np.max(self.xaxis-offset[0]),
                                        np.min(self.xaxis-offset[0]),
                                        np.min(self.vaxis-sysvel),
                                        np.max(self.vaxis-sysvel)],
                                vmin=immin, vmax=immax, aspect='auto')

        self.ax.contour(self.xaxis-offset[0], self.vaxis-sysvel, image,
                        colors='black', levels=np.arange(500)*2*rms+3*rms)
        self.ax.plot([np.min(self.xaxis[xc]), np.max(self.xaxis[xc])], [0,0],
                     '--', color='black')
        self.ax.plot([0,0],[np.min(self.vaxis[vc]-sysvel),
                     np.max(self.vaxis[vc] - sysvel)], '--', color='black')


        if solid and wedge and overplot is not True:
            cbar = plt.colorbar(im)
            cbar.set_label(self.header['BUNIT'])



    def spectrum(self, dx=-1, dy=-1, dv=1e30, offset=[0,0], sysvel=0., rms=-1,
                 title='', overplot=False, fit=False):
        ''' Plots a spectrum

            Keyword arguments:
            dx -- delta x in arcseconds (default is 0)
            dy -- delta y in arcseconds (default is same as dx)
            dv -- delta v in km/s (default is entire v axis)
            offset -- center offset of the map in arcseconds (default is [0,0])
            sysvel -- velocity offset in km/s (default is 0)

            rms -- noise level (default is getRMS())

            title -- title string for the plot (default is none)

            overplot -- do not erase previous plot (default is false)

            fit -- overplot Gaussian fit to the spectrum (default is no)
        '''

        if not self.vaxis.any():
            sys.exit(self.fileName + " is not a line image")

        if dx == -1:
            dx=abs(self.xaxis[1]-self.xaxis[0])

        dy, rms, xc, yc, vc, alpha = self._initPlot(dx, dy, dv, sysvel, rms,
                                                    overplot, offset)

        tmpspectrum = np.zeros(len(vc))
        count = 0
        for i in xc:
            for j in yc:
                if (self.image[vc, j, i] > 3*rms).any():
                    count += 1
                    tmpspectrum += self.image[vc, j, i]

        if count == 0:
            return "No spectra found."
        else:
            print count, " spectra found."
            spectrum = [k/float(count) for k in tmpspectrum]

        if self.ax is None or overplot is not True:
            self.ax = plots.createSpectrumPlot(self.vaxis-sysvel, vc,
                                               self.header, title,
                                               self.header['BUNIT'])

        self.ax.step(self.vaxis[vc]-sysvel, spectrum, where='mid')
        self.ax.plot([np.min(self.vaxis[vc]-sysvel),
                     np.max(self.vaxis[vc] - sysvel)], [0,0], '--',
                     color='black')
        if(fit):
            popt = analysis.gaussFit(self.vaxis[vc]-sysvel, spectrum,
                                     self.header['BUNIT'] )
            self.ax.plot(self.vaxis[vc]-sysvel,
                         analysis.gauss_function(self.vaxis[vc] - sysvel,
                         *popt), color='red')
