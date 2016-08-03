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
    from astropy.modeling import models, fitting
except:
    print 'Error: analysis requires astropy'


def gaussFit(data, inix=0, iniy=0, amp=1.):
   
    xv = np.linspace(np.max(data.xaxis), np.min(data.xaxis), data.xaxis.shape[0])
    yv = np.linspace(np.min(data.yaxis), np.max(data.yaxis), data.yaxis.shape[0])

    x, y = np.meshgrid(xv, yv)
    
    p_init = models.Gaussian2D(amplitude=amp, x_mean=inix, y_mean=iniy, x_stddev=.01 , y_stddev=.1,theta=0.)
    fit_p = fitting.LevMarLSQFitter()
    
#    for i in range(256):
#        for j in range(256):
#            if (data.yaxis[j]>-51. ): #12*data.xaxis[i]-373.8):
#                data.image[j,i]=0.
    
    
    p = fit_p(p_init, x, y, data.image)
    print p
    
    th=p.theta.value
    a=(np.cos(th)**2)/(2*p.x_stddev.value**2) + (np.sin(th)**2)/(2*p.y_stddev.value**2)
    b=-(np.sin(2*th))/(4*p.x_stddev.value**2) + (np.sin(2*th)) /(4*p.y_stddev.value**2)
    c=(np.sin(th)**2)/(2*p.x_stddev.value**2) + (np.cos(th)**2)/(2*p.y_stddev.value**2)

    z=p.amplitude.value*np.exp(-(a*(x-p.x_mean.value)**2 - 2*b*(x-p.x_mean.value)*(y-p.y_mean.value) + c*(y-p.y_mean.value)**2  ))
    plt.contour(x,y,z, linewidths=2)



