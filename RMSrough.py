import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from pylab import *
import os.path
import string

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

import Image

filename = raw_input("Paste .z filename ==> ")
im = loadtxt(filename)
ardat = array(im)
themean = mean(ardat)
rmsrough = std(ardat) #std(ardat*(peaktoval/comppeaktoval))  #standard deviation 
print "RMS Roughness:" + str(rmsrough)
text = open("RMS_roughness.txt",'a')
text.write(str(rmsrough)+',')
text.close()