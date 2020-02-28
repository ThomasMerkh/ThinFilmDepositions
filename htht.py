#Authors : Thomas Merkh merkht@rpi.edu
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from pylab import *
import os.path
import string
import csv
import math
from scipy import fftpack
from scipy import signal
#import pyfits
import numpy as np
import pylab as py
from scipy import weave
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline


def generateIntCircle(radius):
    coords = []
    for i in range(radius+1):
        j = sqrt(radius**2 - i**2)
        sub = array([i, j])
        if floor(j) == j:
            coords.append(sub)
    temp = coords

    for i in range(len(temp)):
        if 0 not in temp[i]:
            coords.append([-temp[i][0], -temp[i][1]])
            coords.append([temp[i][0], -temp[i][1]])
            coords.append([-temp[i][0], temp[i][1]])
        elif temp[i][0] == 0:
            coords.append([temp[i][0], -temp[i][1]])
        elif temp[i][1] == 0:
            coords.append([-temp[i][0], temp[i][1]])
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

            xylist = weave.inline(code, ['rx', 'ry', 'surf', 'N'])     
            xysum = sum(xylist);
            #print len(xylist);
            rsum = rsum + xysum
        
        acor = float(rsum)/float(circPlen)

        out = 1.0/(float(N)**2*float(w)**2)*acor
        autocor[r] = out 
    
    return (autocor, radii)
    
if __name__ == '__main__':
	
    data = np.loadtxt(open("KPZ_2d_noise100000_v2.csv","rb"),delimiter=",",skiprows=0,dtype=float)
    fdata = np.flipud(data)
    resolution = fdata.shape[0]  

    themean = mean(fdata)
    print "The mean surface height is... " + str(themean)
    
    #get omega
    summm = 0
    for k in arange(resolution):
	for j in arange(resolution):
	    summm = summm + (fdata[k][j] - themean)*(fdata[k][j] - themean)
    summm = summm/(resolution*resolution)
      	    
    omega = summm**(0.5)
    print "Omega = " + str(omega)    

    print "Peak To Valley Range: " + str(amax(fdata)-amin(fdata))

    #Subtract the mean surface height
    for k in arange(resolution):
	for j in arange(resolution):
	    fdata[k][j] = fdata[k][j] - themean
    
    acor = autocorrelate(array(fdata), resolution, omega)   #returns a list [list:acor values, list:radii]
    figure()
    plot(acor[1], acor[0])
    xlim(0,max(acor[1]))
    savefig("AutoCorr_100000.png")
    show()

    HHcorr = 2*(omega)**2*(1 - acor[0])    
    loglog(acor[1], HHcorr)        
    
    transX1 = log(acor[1][1:4])  #fit between what range
    transY1 = log(HHcorr[1:4])   #fit between what range
    
    HHFit = polyfit(transX1, transY1, 1)
    
    Alpha = HHFit[0]/2
    print "Alpha: " + str(Alpha)
    
    loglog(exp(transX1), exp(HHFit[1])*exp(HHFit[0]*transX1)) #shows the alpha fit
    xlim(0,max(acor[1]))
    ylim(0.1,1.2*max(HHcorr))
    savefig("Height_Height100000.png") 


    R = []
    Auto = []
    HtHt = []
    for i in arange(size(acor[1])):
        R.append(acor[1][i])
        Auto.append(acor[0][i])
        HtHt.append(HHcorr[i])


    with open('AutoCorrelated100000.csv', 'wb') as f:
        writer1 = csv.writer(f)
        writer1.writerows([R,Auto,HtHt])

    print "Done"
