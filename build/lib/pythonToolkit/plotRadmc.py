"""
    
    Plotting tools for radmc3d output, Christian Brinch (c) 2015
    
    This module provides functions for plotting radmc3d model output

"""

try:
    import matplotlib.pyplot as plt
    import matplotlib.delaunay as triang
except:
    print 'Error: plotRadmc requires matplotlib' 



def contour(x,y,data):
    fig=plt.figure()
    ax=fig.add_subplot(111)
 
    cens,edg,tri,neig = triang.delaunay(x,y)
    im=plt.tripcolor(x, y, tri, data, shading='flat')


    




