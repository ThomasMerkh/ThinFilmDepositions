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
import scipy
from scipy import interpolate
import time

for entry in os.listdir('.'):
    if ('Statistics.txt' in entry and '~' not in entry):
        filen = entry
statsdata = np.loadtxt(filen, delimiter = ' ', skiprows=1, dtype=np.float) ###will throw an error if there are backup tilda files ie: Statistics.txt~ files
interfaceWidths = statsdata[:,2]                                           # slices a column out of array


ltransx1 = log(statsdata[0:20,0])
ltransy1 = log(statsdata[0:20,2])
fit1 = polyfit(ltransx1, ltransy1, 1)

ltransx2 = log(statsdata[100:150,0])
ltransy2 = log(statsdata[100:150,2])
fit2 = polyfit(ltransx2, ltransy2, 1)

ltransx3 = log(statsdata[300:500,0])
ltransy3 = log(statsdata[300:500,2])
fit3 = polyfit(ltransx3, ltransy3, 1)

ltransx4 = log(statsdata[600:800,0])
ltransy4 = log(statsdata[600:800,2])
fit4 = polyfit(ltransx4, ltransy4, 1)

ltransx5 = log(statsdata[900:1000,0])
ltransy5 = log(statsdata[900:1000,2])
fit5 = polyfit(ltransx5, ltransy5, 1)

times = [0,100,300,600,900]
betas = [round(fit1[0],3),round(fit2[0],3),round(fit3[0],3),round(fit4[0],3),round(fit5[0],3)]
print times
print betas
plt.plot(times, betas, 'o')
plt.xlabel("Time step")
plt.ylabel("Beta")
plt.title("Beta measurement vs. time")
plt.ylim(-1,1)
plt.show()

"""
loglog(statsdata[:,0], statsdata[:,2])
loglog(exp(ltransx1), exp(fit1[0]*ltransx1 +  fit1[1]), 'r--', linewidth=3)
xlabel('Simulation step (time)')    
ylabel('Interface Width')        
xlim(0,1000) 
title('RMS roughness over time')
show()  #show the interface width over time
"""


g = scipy.interpolate.UnivariateSpline(times,betas)
xg = times
yg = g(xg)
plt.plot(times,betas, 'o', xg,yg, '-')

plt.xlabel("Time step")
plt.ylabel("Beta")
plt.title("Beta measurement vs. time")
plt.ylim(-1,1)
plt.show()

filen = open("beta_time.txt" , "w")
filen.write( "0 : " + str(round(fit1[0],3)) + "\n" + \
             "100 : " + str(round(fit2[0],3)) + "\n" +\
             "300 : " + str(round(fit3[0],3)) + "\n" +\
             "600 : " + str(round(fit4[0],3)) + "\n" +\
             "900 : " + str(round(fit5[0],3))  )
filen.close()


