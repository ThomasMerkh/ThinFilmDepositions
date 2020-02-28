import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from pylab import *
import os.path
import string
import math
import numpy as np
import pylab as py
from scipy import weave
import time

def generateIntCircle(radius):
    coords = []
    for i in range(radius+1):
        j = sqrt((radius**2) - (i**2))
        sub = array([i, j])
        if floor(j) == j:
            coords.append(sub)
    temp = coords

    for i in range(len(temp)):
        if 0 not in temp[i]:
            coords.append([-1*temp[i][0], -1*temp[i][1]])
            coords.append([temp[i][0], -1*temp[i][1]])
            coords.append([-1*temp[i][0], temp[i][1]])
        elif temp[i][0] == 0:
            coords.append([temp[i][0], -1*temp[i][1]])
        elif temp[i][1] == 0:
            coords.append([-1*temp[i][0], temp[i][1]])
    return array(coords)
    
def autocorrelate(surf, N, w):
    autocor = zeros(N, dtype=float)
    radii = arange(N)
    
    rsum = 0
    rx = 0
    ry = 0
    xysum = 0.0
    xp = 0
    yp = 0
    for r in radii:
        circPoints = array(generateIntCircle(r));
        circPlen = len(circPoints)
        if r == 0:
            circPoints = array([[0.0,0.0]])
            circPlen = len(circPoints)
        rsum = 0   
        #print r
        for points in circPoints:
            rx = int(points[0])
            ry = int(points[1])         
            code = """
            #include <iostream>
            #include <math.h>
            double thexysum = 0;                        
            int xp, yp;

            for (int x = 0; x < N; x++) {
                for (int y = 0; y < N; y++) {
                    xp = x + rx;
                    yp = y + ry;
                    
                    while (xp < 0) { xp = xp + N; }
                    while (xp >= N) { xp = xp - N; }
                    while (yp < 0) { yp = yp + N; }
                    while (yp >= N) { yp = yp - N; }
                    
                    thexysum = thexysum + surf[x*N + y]*surf[xp*N + yp];                                     
                }
            }       
            return_val = thexysum;                      
            """        

            xysum = weave.inline(code, ['rx', 'ry', 'surf', 'N'])
            rsum = rsum + xysum  #this is the total double summation of  h(xi,yj)*h(xi+rx,yj+ry) for a given radius

        acor = float(rsum)/float(circPlen)        
        out = (1.0*acor)/(float(N)**2*float(w)**2)
        autocor[r] = out 
        """
        R(0) ~ 0.99999 as it is supposed to, this is correct.  However usually the corelation at R(r=1) drops to 0.75 and R(r=2) drops tp .54 , .34, .18, .06,....
        This is giving a very low value for alpha.  alpha should be high for smooth, very correlated surfaces, but it is not comming out that way.
        This may be a result of discrete heights.  Example, each '1' is a nm in this 512 nm run, than each particle is 1 nm and the R(r) will only be between 0,1,2,-1 heights
        Rather than heights more realistically such as 0, 0.35, -0.20, etc.  Because of such large particles maintainng a smooth surface, alpha must be low, giving a small slope
        to the H-Ht graph. 
        """
    return (autocor, radii)


#Options
file_front = ''

dir_files = os.listdir('.')    
components = []
statsdata = []
interfaceWidths = []
files = []
XiList = []
AlphaList = []
tList = []

maxiterationnum = 10000
length = 0
resolution =  0
Beta = 0
Alpha = 0
z = 0


#this doesn't loop, this only finds length of simulation and resolution
for entry in dir_files: 
    if ('N ' in entry) and ('.z' in entry):
        components = string.split(entry)
        length = len(components)
        resolution = int(components[1])
        maxiterationnum = int(components[5])+1
        print "Res: " + str(resolution)
        print "MaxIter: " + str(maxiterationnum)
        print "Length: " + str(length)
        break

for entry in dir_files:
    if ('Statistics.txt' in entry and '~' not in entry):
        filen = entry
statsdata = np.loadtxt(filen, delimiter = ' ', skiprows=1, dtype=np.float) ###will throw an error if there are backup tilda files ie: Statistics.txt~ files
interfaceWidths = statsdata[:,2]                                           # slices a column out of array
fig = py.figure()
fig = fig.add_subplot(1,1,1)
fig.plot(statsdata[:,0], statsdata[:,2])
fig.set_yscale('log')
fig.set_xscale('log')
ltransx = log(statsdata[0:20,0])
ltransy = log(statsdata[0:20,2])
fit = polyfit(ltransx, ltransy, 1)
#loglog(exp(ltransx), exp(fit[0]*ltransx +  fit[1]), 'r--', linewidth=3)
Beta = fit[0]
print "Beta: " + str(Beta)
xlabel('Simulation step (time)')    
ylabel('Lattice Units')        
xlim(0,maxiterationnum) 
ylim(0.8,30)
title('')

fig.plot([1,2,6,16,40,63,100,251,398,631,1000] ,[20,20,21,22,22,23,21,19,20,19,19] , '-.o')

show()  #show the interface width over time

for i in range(length-1):
    file_front += components[i]
    file_front += ' '
    
subcomponents = string.split(components[length-1], '_')
file_front += subcomponents[0]
file_front += '_n' #N 512 L 512 T 632 D 0.238 Stick 1.00 NoDiff BA Cos_n is the new form
file_extension = '.z'
      
for i in range(maxiterationnum):   #checks to large number
    if os.path.isfile(file_front + str(i) + file_extension):
        files.append(i)
        #files will only have time step numbers that data was saved at ie: [1,2,4,6,10,16 ... 398 ...]

print "resolution: " + str(resolution)   
MM = 0
filen = open("PeakToValley.txt", "w")
for i in files:
    print 'step:' + str(i)
    print str(file_front + str(i) + file_extension) + " exists:", os.path.isfile(file_front + str(i) + file_extension)
    #print repr(open(str(file_front + str(i) + file_extension), 'rb').read((512*512)))  I used this to find if any \n or \t characters are in data which there may be
    data = np.loadtxt(file_front + str(i) + file_extension, dtype=float)  #if int, then the mean height will not be zero after subtraction
    fdata = np.flipud(data)
    data_resolution = fdata.shape[0]    #size of simulation should match with the resolution found via file name earlier
    themean = np.mean(fdata)
    print themean, "themean"
    print "EXTENT: " + str(amax(fdata)-amin(fdata))  #max() and min() were wrong, they resulted in bools not numbers
    print "Max surface height:", amax(fdata)
    print "Minimum surface height:", amin(fdata)
    filen.write("step"+str(i)+' '+str(amax(fdata)-amin(fdata))+'\n')
    for k in arange(resolution):
        for j in arange(resolution):
            fdata[k][j] = (fdata[k][j] - themean)
    acor = autocorrelate(array(fdata), resolution, interfaceWidths[i-1])

    #HHCOR and ALPHA analysis and Autocorrelation
    """
    loglog(acor[1],acor[0])
    xlabel("Radii")
    ylabel("R(r)")
    title("Autocorrelation function")
    """
    HHcorr = (2*(interfaceWidths[i-1])**2)*(1 - acor[0]) 
    if max(HHcorr) > MM:
        MM = max(HHcorr)
    loglog(acor[1], HHcorr)        
    xlabel('Radii')
    ylabel('Ht-Ht')
    xlim(1,resolution)
    ylim(0.1,MM*1.2)
    title("Ht-Ht function of SiO2") 
    transX1 = log(acor[1][1:4])
    transY1 = log(HHcorr[1:4])
    HHFit = polyfit(transX1, transY1, 1)
    print "HHfitforalpha: " + str(HHFit)
    Alpha = HHFit[0]/2
    print "Alpha: " + str(Alpha) 
filen.close()
show()