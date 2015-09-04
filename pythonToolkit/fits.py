"""

fits tools, Christian Brinch (c) 2015

This module contains classes and methods to read and plot FITS files

"""

from pythonToolkit import standards

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
    def __init__(self,file):
        try:
            hdulist = fits.open(file)
        except:
            print 'Error: FITS file not found'

        self.header=hdulist[0].header
        self.image =hdulist[0].data
        hdulist.close()


    def plot(self):
        cellsize=self.header['CDELT2']*np.pi/180.*206264.806
        xaxis=(np.arange(self.header['NAXIS1'])-self.header['CRPIX1']) * cellsize
        yaxis=(np.arange(self.header['NAXIS2'])-self.header['CRPIX2']) * cellsize

        plt.close()
        plt.figure()
        ax=plt.subplot(111)
        ax.minorticks_on()
        plt.xlabel(self.header['CTYPE1'])
        plt.ylabel(self.header['CTYPE2'])

        im=ax.imshow(self.image.sum(axis=0), cmap=plt.cm.jet, alpha=0.9, interpolation='nearest',
             origin='lower', extent=[np.min( xaxis ),np.max( xaxis ),np.min( yaxis ),np.max( yaxis )],
             vmin=np.min(self.image),vmax=np.max(self.image))

        cbar = plt.colorbar(im)
        cbar.set_label(self.header['bunit'])


