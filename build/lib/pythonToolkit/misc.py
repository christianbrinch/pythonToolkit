"""
Utility functions, Christian Brinch (c) 2015

This module contains various useful functions
"""

try:
    import numpy as np
except:
    print 'Error: misc requires numpy' 


def rebin(b, shape):
    a=np.zeros( (b.shape[0] , shape[1]) )

    for i in range(shape[1]):
        a[:,i] = np.asarray(b)

    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)


