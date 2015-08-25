"""

    radmc3d to LIME converter, Christian Brinch (c) 2015

    This module converts radmc3d outputinto LIME readable format.
    Also useful for plotting radmc3d output

"""

from pythonToolkit import standards

try:
  from astropy.io import ascii
except:
  print 'Error: radmc2LIME requires astropy'

try:
  import numpy as np
except:
  print 'Error: radmc2LIME requires numpy' 

try:
  from scipy.integrate import simps
except:
  print 'Error: radmc2LIME requires scipy' 

try:
    import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation, TriAnalyzer, UniformTriRefiner
except:
    print 'Error: radmc2LIME requires matplotlib' 




class readRadmc3d():
  def __init__(self,gridfile,temperature,density):
    self.T    = ascii.read(temperature, format='basic', data_start=3)
    self.rho  = ascii.read(density,     format='basic', data_start=3)
    self.dust = ascii.read('dustkappa_silicate.inp', data_start=2)

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

      self.x= r * np.sin(t)
      self.y= r * np.cos(t)



  def meanOpacity(self):
    def dB_T(T,lam):
      c1 = 2*6.62e-34*(3e8)**2.
      c2 = 6.62e-34*3e8/1.38e-23
      return c1*c2*np.exp(c2/(lam * T)) / (lam**6 * T**2 *(np.exp(c2/(lam * T))-1)**2) 

    self.kappa = np.zeros(len(self.x))
    for i in range(len(self.x)):
      f = [dB_T(self.T[i]['1'],nu) for nu in self.dust[0][:]]
      nom = simps( self.dust[1][:]**(-1) * f, self.dust[0][:])
      denom = simps( f, self.dust[0][:])
      self.kappa[i] = (nom / denom)**(-1.)
    
    self.absorp = self.kappa * self.rho[0:len(self.x)]['1']


  def contour(self, prop=None):
    if(prop=='Tdust'):
      data=self.T[0:len(self.x)]['1']
      label="Dust temperature [K]"
    elif(prop=='rho'):
      data=np.log(self.rho[0:len(self.x)]['1'])
      label="Log some density [?]"
    else:
      try:
        data=prop
	label='Unknown property'
      except:
        print "Please specify property to plot. Valid properties are [Tdust, rho]"
        return

    x_au = self.x/1.49e13
    y_au = self.y/1.49e13

    tri=Triangulation(x_au,y_au)

    plt.close('all')
    plt.figure()
    ax=plt.subplot(111)
    ax.set_xscale("log",nonposx='clip')
    ax.set_yscale("log",nonposy='clip')
    ax.set_xlim([0.1,200])
    ax.set_ylim([0.1,200])
    plt.xlabel('r [AU]')
    plt.ylabel('z [AU]')

    nmax=data.max()
    nmin=data.min()
    levels=np.arange(12) * (nmax-nmin)/12. + nmin
    plt.tricontourf(tri, data, levels)
    cbar=plt.colorbar()
    cbar.set_label(label)

    CS=plt.tricontour(tri, data, levels, colors='black', linewidths=1.5)
    plt.clabel(CS, fontsize=8, inline=1)
    cbar.add_lines(CS)

    plt.triplot(tri, color='black', alpha=0.2)



