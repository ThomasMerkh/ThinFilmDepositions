import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from pylab import *
import os.path
import string

#NanoSurfer by Tripp Robert Spivey 2012
#Script for fast automatic analysis of data from full 3-d Monte Carlo Simulation of Thin Film Growth

#Options
colormap = 'hot' #Options: 'classic', 'hot', 'pink', 'earth', 'heat', 'bone'
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

file_front = ''
maxiterationnum = 10000
dir_files = os.listdir('.')    
resolution =  0
components = []
for entry in dir_files:
    if ('Skips' in entry or 'NoDiff' in entry) and (('.cs' in entry) or ('.xcs' in entry) or ('.z' in entry)):
        components = string.split(entry)
        length = len(components)
        resolution = int(components[1])
        maxiterationnum = int(components[5])+1
        break

for i in range(length - 1):
    file_front += components[i]
    file_front += ' '
    
subcomponents = string.split(components[length - 1], '_')
file_front += subcomponents[0]
file_front += '_n'
files = []
file_extension = '.cs'
if surface_datatype == 'ycross':
    file_extension = '.cs'
elif surface_datatype == 'xcross':
    file_extension = '.xcs'
elif surface_datatype == 'zview':
    file_extension = '.z'
else:
    print "invalid surface_datatype"
        
for i in range(maxiterationnum): #checks to large number
    if os.path.isfile(file_front + str(i) + file_extension):
        files.append(i)

   
for i in files:
    data = np.loadtxt(file_front + str(i) + file_extension)
    #pdata = np.loadtxt(file_front + str(i) + '.p')
    ratios = []
    thick = []
    """
    for j in pdata:
        if j[1]!=0:
            ratios.append(float(j[2])/(float(j[1]) + float(j[2])))
            thick.append(j[0])
    figure()
    plot(thick, ratios, '.')
    savefig("p" + str(i) + '.png')
    """
    fdata = np.flipud(data)
    data_resolution = fdata.shape[0]
    
    make_image('n' + str(i), fdata)
    print 'created image for step:' + str(i)
    
files2 = []
for i in range(maxiterationnum): #checks to large number
    if os.path.isfile(file_front + str(i) + '.xcs'):
        files2.append(i)

   
for i in files:
    data = np.loadtxt(file_front + str(i) + '.xcs')

    fdata = np.flipud(data)
    data_resolution = fdata.shape[0]
    
    make_image('nn' + str(i), fdata)
    print 'created image for step:' + str(i)



