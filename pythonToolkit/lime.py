"""

lime tools, Christian Brinch (c) 2015

This module deals with lime models
"""

from pythonToolkit import standards
from pythonToolkit import plots
from pythonToolkit.constants import *


import sys
try:
    import matplotlib.pyplot as plt
    import matplotlib.colors as cols
    import matplotlib.tri as tri
    from matplotlib import cm
except:
    print 'Error: Matplotlib not found'
try:
    import numpy as np
except:
    print 'Error: Numpy not found'
try:
    from astropy.io import ascii
except:
    print 'Error: lime requires astropy'


class read():
    def __init__(self,fileName):
        try:
            self.data=np.loadtxt(fileName, skiprows=1, usecols=(0,1,2,3,4,5), unpack=True)
        except:
            sys.exit("Error: LIME .pop file not found")
        self.r=[np.sqrt(self.data[0][i]**2+self.data[1][i]**2)/AUm for i in range(len(self.data[0]))]
        self.z=[np.abs(self.data[2][i])/AUm for i in range(len(self.data[0]))]
        self.triang = tri.Triangulation(self.r,self.z)

    def plotDensity(self):
        ax,cmap,fig = plots.createModelPlot([0,np.max(self.r),0,np.max(self.z)], 'Gas number density')
        im=plt.tricontourf(self.triang, [np.log10(self.data[3][i]/1e6) for i in range(len(self.data[0]))])
        cbar = plt.colorbar(im)
        cbar.set_label('cm$^{-3}$')

    def plotTemperature(self):
        ax,cmap,fig = plots.createModelPlot([0,np.max(self.r),0,np.max(self.z)], 'Dust temperature')
        im=plt.tricontourf(self.triang, [self.data[4][i] for i in range(len(self.data[0]))])
        cbar = plt.colorbar(im)
        cbar.set_label('K')

    def plotAbundance(self):
        ax,cmap,fig = plots.createModelPlot([0,np.max(self.r),0,np.max(self.z)], 'Relative CO abundance')
        im=plt.tricontourf(self.triang, [np.log(self.data[5][i]) for i in range(len(self.data[0]))])
        cbar = plt.colorbar(im)
        cbar.set_label('')
