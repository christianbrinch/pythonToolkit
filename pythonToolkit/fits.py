"""

fits tools, Christian Brinch (c) 2015

This module contains classes and methods to read and plot FITS files

"""

from pythonToolkit import standards
from pythonToolkit import plots
import sys
try:
    from astropy.io import fits
except:
    print 'Error: fits requires astropy'

try:
    import numpy as np
except:
    print 'Error: lime requires numpy'

try:
    import matplotlib.pyplot as plt
except:
    print 'Error: lime requires matplotlib'



class read():
    def __init__(self,fileName, coordSys='rel'):
        try:
            hdulist = fits.open(fileName)
        except:
            print 'Error: FITS file not found'

        self.fileName=fileName
        self.header=hdulist[0].header
        self.image =hdulist[0].data
        hdulist.close()

        # Collapse cube to three dimensions
        if (np.ndim(self.image) == 4):
            self.image=self.image.sum(axis=0)

        # Define xaxis, yaxis, and velocity axis
        # In case of continuum images, collapse frequency axis
        cellsize=self.header['CDELT2']*np.pi/180.*206264.806
        if(coordSys == 'rel'):
          self.xaxis=(np.arange(self.header['NAXIS1'])-self.header['CRPIX1']) * cellsize
          self.yaxis=(np.arange(self.header['NAXIS2'])- self.header['CRPIX2']) * cellsize
        elif(coordSys == 'abs'):
          self.xaxis=(((np.arange(self.header['NAXIS1'])-self.header['CRPIX1'])*self.header['CDELT1']/15. + self.header['CRVAL1']  /360*24.-16)*60-27)*60
          self.yaxis=(((np.arange(self.header['NAXIS2'])-self.header['CRPIX2'])*self.header['CDELT2'] + self.header['CRVAL2']       +24)*60+40)*60
        else:
          print "Unknown choice of coordinate system"
          return
        
        if(self.header['NAXIS3'] >1):
            if(self.header['CTYPE3'] == 'VELO-LSR'):
                self.vaxis=(np.arange(self.header['NAXIS3']) * self.header['CDELT3'] ) / 1000.
            if(self.header['CTYPE3'] == 'FREQ'):
                self.vaxis=(np.arange(self.header['NAXIS3']) * self.header['CDELT3']*3e8/self.header['CRVAL3']) / 1000.
        else:
            self.vaxis=np.array([])
            if (np.ndim(self.image) == 3):
                self.image=self.image.sum(axis=0)

    def printHeader(self):
        for i in self.header:
            if i != 'HISTORY':
                print i, self.header[i]

        key=[x for x in self.header['HISTORY'] if ('Beam' in x or 'restoration' in x ) and 'arcsec' in x]
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
            rms = np.std(self.image[0:self.header['NAXIS2']/4.,0:self.header['NAXIS1']/4.])
        else:
            print "This seems to be a model. Contours are 10% of peak value."
            rms = self.image.max()/10./3.

        print "RMS: ",rms
        return rms

    def _plotBeam(self,d,ax,xc,yc):
        if ('HISTORY' in self.header) and (not 'BMAJ' in self.header):
            key=[x for x in self.header['HISTORY'] if 'Beam' in x and 'arcsec' in x]
            if key:
                beaminfo = key[0].split()
                offset=float(beaminfo[3])/2.+float(beaminfo[3])/4.
                plots.beam(ax, np.min(self.xaxis[xc])+offset, np.min(self.yaxis[yc])+offset, float(beaminfo[3]), float(beaminfo[5]), float(beaminfo[9]))
            else:
                key=[x for x in self.header['HISTORY'] if 'restoration' in x  and 'arcsec' in x]
                if key:
                    beaminfo = key[0].split()
                    offset=float(beaminfo[2])/2.+float(beaminfo[4])/4.
                    plots.beam(ax, np.min(self.xaxis[xc])+offset, np.min(self.yaxis[yc])+offset, float(beaminfo[2]), float(beaminfo[4]), float(beaminfo[8]))
        elif ('BMAJ' in self.header) and ('BMIN' in self.header):
            bmaj=float(self.header['BMAJ'])*3600.
            bmin=float(self.header['BMIN'])*3600.
            offset=bmaj/2.+bmin/4.*5.
            plots.beam(ax, np.min(self.xaxis[xc])+offset, np.min(self.yaxis[yc])+offset, bmaj, bmin, float(self.header['BPA']))

    def continuum(self, dx=1e30, nobeam=False, rms=-1):
        if(self.vaxis.any()):
            print self.fileName + " is not a continuum image"
            return

        if(rms<0):
            rms=self.getRMS()

        self.xaxis=(((np.arange(self.header['NAXIS1'])-self.header['CRPIX1'])*self.header['CDELT1']/15. + self.header['CRVAL1']  /360*24.-16)*60-27)*60
        self.yaxis=(((np.arange(self.header['NAXIS2'])-self.header['CRPIX2'])*self.header['CDELT2'] + self.header['CRVAL2']       +24)*60+40)*60

        xc = np.where((self.xaxis > -dx) & (self.xaxis < dx))
        yc = np.where((self.yaxis > -dx) & (self.yaxis < dx))

        ax,cmap=plots.createFigure(self.xaxis,self.yaxis,xc,yc,self.header['OBJECT'])

        im=ax.imshow(self.image, cmap=cmap, alpha=1.0, interpolation='nearest',origin='lower', extent=[np.max( self.xaxis ),np.min( self.xaxis ),np.min( self.yaxis ),np.max( self.yaxis )], vmin=np.nanmin(3*rms),vmax=np.nanmax(self.image), aspect='auto')

        cbar = plt.colorbar(im)
        cbar.set_label('Intensity ('+self.header['BUNIT']+')')

        if(not nobeam):
            self._plotBeam(self,ax,xc,yc)

        return ax


    def moment(self, vmi, vma, dx=1e30, dv=1e30, nobeam=False, mom=0, rms=-1, overplot=False):
        if(self.vaxis == []):
            sys.exit(image + " is not a line image")

        self.xaxis=(((np.arange(self.header['NAXIS1'])-self.header['CRPIX1'])*self.header['CDELT1']/15. + self.header['CRVAL1']  /360*24.-16)*60-27)*60
        self.yaxis=(((np.arange(self.header['NAXIS2'])-self.header['CRPIX2'])*self.header['CDELT2'] + self.header['CRVAL2']       +24)*60+40)*60

        self.vaxis=self.vaxis-self.vaxis[194]
    
#        vc = np.where((self.vaxis > -dv) & (self.vaxis < dv))
        vc = np.where((np.abs(self.vaxis) > vmi) & (np.abs(self.vaxis) < vma))
        xc = np.where((self.xaxis > -dx) & (self.xaxis < dx))
        yc = np.where((self.yaxis > -dx) & (self.yaxis < dx))

        if(rms<0):
            rms=self.getRMS()

        moment0 = (self.image[vc[0]]*(np.abs(self.vaxis[1]-self.vaxis[0]))).sum(axis=0)
        moment1 =np.zeros((self.header['NAXIS2'],self.header['NAXIS1']),float)
        for j in range(self.header['NAXIS1']):
            for i in range(self.header['NAXIS2']):
                if (moment0[i,j] > 3*rms):
                    moment1[i,j]= (self.image[vc[0],i,j]*(self.vaxis[vc[0]])).sum()/moment0[i,j]
                else:
                    moment1[i,j]=np.nan




        plt.clf()
        fig = plt.figure(1)
        ax=fig.add_subplot(111, aspect='equal')
        ax.set_title(self.header['OBJECT'])
        ax.set_xlim(np.max( self.xaxis[xc] ),np.min( self.xaxis[xc] ))
        ax.set_ylim(np.min( self.yaxis[yc] ),np.max( self.yaxis[yc] ))
        ax.set_xlabel('R. A. offset (arcsec)')
        ax.set_ylabel('DEC offset (arcsec)')
        ax.minorticks_on()
        plt.tick_params(axis='both', which='both', width=0.4)
        cmap = plt.cm.jet
        cmap.set_bad('white',1.)

        moment=moment1
        im=ax.imshow(moment, cmap=cmap, alpha=0.9, interpolation='nearest',
                     origin='lower', extent=[np.max( self.xaxis ),np.min( self.xaxis ),np.min( self.yaxis ),np.max(self.yaxis)],
                     vmin=np.nanmin(moment),vmax=np.nanmax(moment), aspect='auto')


        plt.subplots_adjust(bottom=0.15,left=0.15)
        cbar = plt.colorbar(im)
        cbar.set_label('Intensity integrated velocity')
        ax.contour(self.xaxis,self.yaxis,moment0, colors='black', levels=np.arange(25)*2*rms+3*rms)
        plt.plot([26.911676,26.920047],[-50.71791,-51.28717],'o',color='green')
        x=np.linspace(0,2*np.pi,100)
        x0=26.932976536
        y0=-50.44454402
        vx=-0.0005886
        vy=-0.028301
        t=25.25
#        plt.plot(x0+vx*t+0.469/15.*np.cos(x),y0+vy*t+0.469*np.sin(x), color='green', lw=1.3)
#        plt.plot(x0+vx*t+0.132/15.*np.cos(x),y0+vy*t+0.132*np.sin(x), color='green', lw=1.3)


