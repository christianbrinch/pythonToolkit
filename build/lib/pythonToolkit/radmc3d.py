"""

radmc3d tools, Christian Brinch (c) 2015

This module deals with radmc3d models
"""

from pythonToolkit import standards
from pythonToolkit import plots
from pythonToolkit.constants import *
import sys
import os

try:
    from astropy.io import ascii
except:
    print 'Error: radmc3d requires astropy'

try:
    import numpy as np
except:
    print 'Error: radmc3d requires numpy'

try:
    from scipy.integrate import simps
except:
    print 'Error: radmc3d requires scipy'

try:
    import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation, TriAnalyzer, UniformTriRefiner
except:
    print 'Error: radmc3d requires matplotlib'




class readModel():
    def __init__(self,dir='./'):
        try:
            f=open(dir+'model.py','r')
            f.close()
        except:
            print 'Error: no model is present. Please specify path to model'

        try:
            density=dir+'dust_density.inp'
            gridfile=dir+'amr_grid.inp'
        except:
            print 'Error: density and/or gridfile not found. Please run radmc3d.prepare() first'

        try:
            temperature=dir+'dust_temperature.dat'
        except:
            print 'Warning: temeprature.dat does not exist'

        try:
            dustfile=dir+'dustkappa_silicate.inp'
        except:
            print 'Warning: no dust silicate file found'


        self.T    = ascii.read(temperature, format='basic', data_start=3)
        self.rho  = ascii.read(density,     format='basic', data_start=3)
        self.dust = ascii.read(dustfile,                    data_start=2)

        f = open(gridfile, 'r')
        form  = float(f.readline())
        grid_style  = float(f.readline())
        crd_system  = int(f.readline())
        if crd_system<100:
            self.crd_sys = 'car'
        elif ((crd_system>=100)&(crd_system<200)):
            self.crd_sys = 'sph'
        elif ((crd_system>=200)&(crd_system<300)):
            self.crd_sys = 'cyl'

        grid_info   = float(f.readline())
        dum         = f.readline().split()
        self.act_dim = [int(dum[i]) for i in range(len(dum))]
        dum         = f.readline().split()
        self.nx,self.ny,self.nz    = int(dum[0]), int(dum[1]), int(dum[2])
        self.nxi,self.nyi,self.nzi = self.nx+1, self.ny+1, self.nz+1

        self.xi           = np.zeros(self.nx+1, dtype=np.float64)
        self.yi           = np.zeros(self.ny+1, dtype=np.float64)
        self.zi           = np.zeros(self.nz+1, dtype=np.float64)

        for i in range(self.nxi): self.xi[i] = float(f.readline())
        for i in range(self.nyi): self.yi[i] = float(f.readline())
        for i in range(self.nzi): self.zi[i] = float(f.readline())
        f.close()

        if self.crd_sys=='car':
            self.x = (self.xi[0:self.nx] +  self.xi[1:self.nx+1]) * 0.5
            self.y = (self.yi[0:self.ny] +  self.yi[1:self.ny+1]) * 0.5
            self.z = (self.zi[0:self.nz] +  self.zi[1:self.nz+1]) * 0.5
        else:
            self.x = np.sqrt(self.xi[0:self.nx] * self.xi[1:self.nx+1])
            self.y = (self.yi[0:self.ny] +  self.yi[1:self.ny+1]) * 0.5
            self.z = (self.zi[0:self.nz] +  self.zi[1:self.nz+1]) * 0.5
            r=np.zeros(self.nx*self.ny)
            t=np.zeros(self.nx*self.ny)
            for k in range(self.nx*self.ny):
                r[k] = self.x[k % self.nx]
                t[k] = self.y[k / self.nx]



        self.T = np.reshape(self.T['1'], (self.nx,self.ny), order='C')
        self.rho = np.reshape(self.rho['1'], (self.nx,self.ny), order='C')



    def meanOpacity(self):
        def dB_T(T,freq):
            hnu_kt = hh*freq/(kk*T)
            return 2*hh**2*freq**4/cc**2 * np.exp(hnu_kt)/(T**2*(np.exp(hnu_kt)-1)**2)

        self.kappa  = np.zeros(len(self.x))

        for i in range(len(self.x)):
            f = [dB_T(self.T[i],nu) for nu in cc/self.dust[0][:]]
            nom = simps( self.dust[1][:]**(-1) * f, cc/(1e-4*self.dust[0][:]))
            denom = simps( f, cc/(1e-4*self.dust[0][:]))
            self.kappa[i] = (nom / denom)**(-1.)


    def plot(self, prop='density'):
        r,t = np.meshgrid(self.x/AU, np.pi/2.-self.y)
        ax, fig = plots.createFigure('Radmc3d %s model' % prop, logaxis=True)
        ax.set_xlim(np.min(r*np.cos(t)), np.max(r*np.cos(t)))
        ax.set_ylim(np.min(r*np.sin(t)), np.max(r*np.sin(t)))
        ax.set_xlabel('r (AU)')
        ax.set_ylabel('z (AU)')
        if(prop == 'density'):
            cax = ax.contourf(r*np.cos(t),r*np.sin(t), np.log10(1e-24 +
                              self.rho/(2.37*1.67e-24)*100),
                              levels=np.linspace(1.2,13,100),
                              cmap=plt.cm.spectral)
            cbar = fig.colorbar(cax)
        elif(prop == 'temperature'):
            cax = ax.contourf(r*np.cos(t),r*np.sin(t), np.log10(self.T),
                              levels=np.linspace(0.8,2.8,100),
                              cmap=plt.cm.spectral)
            cbar = fig.colorbar(cax)
            import scipy.ndimage
            T = scipy.ndimage.gaussian_filter(self.T, sigma=3.0, order=0)
            x = scipy.ndimage.zoom(r*np.cos(t), 1)
            y = scipy.ndimage.zoom(r*np.sin(t), 1)
            cax = ax.contour(x,y, T ,
                levels=[10,20,40,100,300], linewidths=3, colors='black')
            ax.clabel(cax, inline=1, fontsize=10)
        elif(prop == 'grid'):
            for r0 in self.x/AU:
                ax.plot(r0*np.cos(np.pi/2.-self.y), r0*np.sin(np.pi/2.-self.y), 'k-')
            for t0 in np.pi/2.-self.y:
                ax.plot(self.x/AU*np.cos(t0), self.x/AU*np.sin(t0), 'k-')
        else:
            print "Unknown property"






def prepare(inputFile='model'):
    try:
        myModel = __import__(inputFile)
        reload(myModel)
    except:
        print 'Error: model file not found or syntax error in model.py'
        return

    par = sys.modules[inputFile]

    if len(par.rhod.shape) == 2:
        par.rhod = np.expand_dims(par.rhod, axis=0)

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
                f.write('%s\n' % par.rhod[iz,iy,ix])
    f.close()

    f = open('dustopac.inp', 'w')
    f.write('2               Format number of this file\n')
    f.write('1               Nr of dust species\n')
    f.write('======================================================================\n')
    f.write('1               Way in which this dust species is read\n')
    f.write('0               0=Thermal grain\n')
    f.write('silbeta         Extension of name of dustkappa_***.inp file\n')
    f.write('----------------------------------------------------------------------\n')
    f.close()

    f = open('radmc3d.inp', 'w')
    f.write('nphot = %s\n' % par.nphot)
    f.write('scattering_mode_max = 0\n')
    f.write('istar_sphere = 1\n')
    f.close()


def runModel(arg='mctherm'):
    os.system('radmc3d.x '+arg)
