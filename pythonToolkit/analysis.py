"""
Analysis functions, Christian Brinch (c) 2016

This module contains various functions for analysis
"""


try:
    import numpy as np
except:
    print 'Error: analysis requires numpy'

try:
    import matplotlib.pyplot as plt
except:
    print 'Error: analysis requires matplotlib'

try:
    from scipy.optimize import curve_fit
except:
    print 'Error: analysis requires scipy'

try:
    from astropy.modeling import models, fitting
except:
    print 'Error: analysis requires astropy'


def gauss_function(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def gaussFit(vaxis, spectrum, unit):
    popt, pcov = curve_fit(gauss_function, vaxis, spectrum, p0 = [np.max(spectrum),
                           vaxis[np.where(spectrum == np.max(spectrum))], 1])
    print "Intensity:     Velocity:    Width:"
    print('{0:3.2f} {1:5s}   {2:3.2f} km/s    {3:3.2f} km/s'.format(popt[0],
           unit, popt[1], popt[2]))

    return popt

def twoD_gauss_function((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta,
                        offset):
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)
                            + c*((y-yo)**2)))
    return g.ravel()

def twoD_gaussFit(fits,image):
    peak = np.max(image)
    pos = np.where(image == peak)
    beam = fits._getBeam(fits, fits.ax, [0], [0], True, True)
    initial_guess = (peak, pos[0][0], pos[1][0], beam[5], beam[6], beam[7], 0.)

    x, y = np.meshgrid(fits.xaxis, fits.yaxis)

    popt, pcov = curve_fit(twoD_gauss_function, (x, y), image.ravel(), p0=initial_guess)

    return popt
