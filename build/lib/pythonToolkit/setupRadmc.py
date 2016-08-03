"""

Setup radmc, Christian Brinch (c) 2015

This module reads a radmc3d input model and generates input files 
"""

from pythonToolkit.constants import *

try:
    import numpy as np
except:
    print 'Error: setupRadmc requires numpy' 


def rebin(b, shape):
    a=np.zeros( (b.shape[0] , shape[1]) )

    for i in range(shape[1]):
        a[:,i] = np.asarray(b)

    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)


def setupRadmc(inputFile):
    import sys
    __import__(inputFile)
    par = sys.modules[inputFile]


    lam12   = par.lambda1 * (par.lambda2/par.lambda1)**(np.arange(par.n12*1.)/(1.*par.n12))
    lam23   = par.lambda2 * (par.lambda3/par.lambda2)**(np.arange(par.n23*1.)/(1.*par.n23))
    lam34   = par.lambda3 * (par.lambda4/par.lambda3)**(np.arange(par.n34*1.)/(1.*(par.n34-1.)))
    lambd  = np.concatenate([np.concatenate([lam12,lam23]),lam34])
    nlam    = len(lambd)

    f = open('wavelength_micron.inp', 'w')
    f.write('%s\n' % nlam)
    for ilam in range(nlam):
        f.write('%s\n' % lambd[ilam])
    f.close()

    f = open('stars.inp', 'w') 
    f.write('2\n')
    f.write('1 %s\n' % nlam)
    f.write('\n')
    f.write('%s %s %s %s %s\n' % (par.rstar,par.mstar,par.pstar[0],par.pstar[1],par.pstar[2]))
    f.write('\n')
    for ilam in range(nlam):  
        f.write('%s\n' % lambd[ilam])
    f.write('\n')
    f.write('%s' % -par.tstar)
    f.close()



    f = open('amr_grid.inp', 'w')
    f.write('%s\n' % 1)                      # iformat
    f.write('%s\n' % 0)                      # AMR grid style  (0=regular grid, no AMR)
    f.write('%s\n' % 100)                    # Coordinate system
    f.write('%s\n' % 0)                      # gridinfo
    f.write('%s %s %s\n' % (1,1,0))          # Include x,y,z coordinate
    f.write('%s %s %s\n' % (par.nx,par.ny,1))        # Size of grid
    for i in range(par.nx+1):
        f.write('%s\n' % par.xi[i])              # X coordinates (cell walls)
    for i in range(par.ny+1):
        f.write('%s\n' % par.yi[i])              # Y coordinates (cell walls)
    for i in range(par.nz+1): 
        f.write('%s\n' % par.zi[i])              # Z coordinates (cell walls)
    f.close()

    f = open('dust_density.inp', 'w')
    f.write('%s\n' % 1)                     # Format number
    ncell = par.nx*par.ny*par.nz
    f.write('%s\n' % ncell)               # Nr of cells
    f.write('%s\n' % 1)                      # Nr of dust species
    for iz in range(par.nz):
        for iy in range(par.ny):
            for ix in range(par.nx):
                f.write('%s\n' % par.rhod[ix,iy])
    f.close()

    f = open('dustopac.inp', 'w')
    f.write('2               Format number of this file\n')
    f.write('1               Nr of dust species\n')
    f.write('======================================================================\n')
    f.write('1               Way in which this dust species is read\n')
    f.write('0               0=Thermal grain\n')
    f.write('silicate        Extension of name of dustkappa_***.inp file\n')
    f.write('----------------------------------------------------------------------\n')
    f.close()

    f = open('radmc3d.inp', 'w')
    f.write('nphot = %s\n' % par.nphot)
    f.write('scattering_mode_max = 0')
    f.close()

