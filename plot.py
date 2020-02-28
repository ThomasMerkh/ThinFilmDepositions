import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from pylab import *
#import os.path
#import string
import csv
import math
#from scipy import fftpack
#from scipy import signal
#import pyfits
import numpy as np
import pylab as py
#from scipy import weave
#from scipy import interpolate
#from scipy.optimize import curve_fit
#from scipy.interpolate import UnivariateSpline

#from numpy import genfromtxt 

colormap = 'bone' #Options: 'classic', 'hot', 'pink', 'earth', 'heat', 'bone'
surface_datatype = 'zview' #Options: 'xcross', 'ycross', 'zview'
the_dpi = 80

def make_image(filename, data):
    fig = plt.figure() 
    dimensions = float(data_resolution)/float(the_dpi)
    fig.set_size_inches(dimensions, dimensions)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    
    if colormap == 'classic':
        ax.imshow(data, aspect = 'normal', cmap=cm.gray)
    elif colormap == 'hot':
        ax.imshow(data, aspect = 'normal', cmap=cm.hot)
    elif colormap == 'pink':
        ax.imshow(data, aspect = 'normal', cmap=cm.pink)
    elif colormap == 'earth':
        ax.imshow(data, aspect = 'normal', cmap=cm.gist_earth)
    elif colormap == 'heat':
        ax.imshow(data, aspect = 'normal', cmap=cm.gist_heat)
    elif colormap == 'bone':
        ax.imshow(data, aspect = 'normal', cmap=cm.bone)
    else:
        ax.imshow(data, aspect = 'normal', cmap=cm.gray)
        
    plt.savefig(filename + '.png', dpi = the_dpi)

data = np.genfromtxt('KPZ_2d_noise50000.csv', delimiter=',')
fdata = np.flipud(data)
data_resolution = fdata.shape[0]
make_image('50000_', fdata)

print "Program finished"