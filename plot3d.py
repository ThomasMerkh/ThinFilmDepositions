###Written by Thomas Merkh, Class of 2016

import matplotlib
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from pylab import *
import pylab as py
import numpy as np

N = 512

HISTO = True
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

X = np.zeros(N*N, dtype= float)
Y = np.zeros(N*N, dtype = float)
X =X.reshape(N,N)
Y = Y.reshape(N,N)
for i in range(N):
    for j in range(N):
        X[i][j] = j
        Y[j][i] = j
Z = np.loadtxt('N 512 L 512 T 1000 D 0.238 NoDiff SoS RevCos_n1000.z', dtype=float)
print X.shape, Y.shape, Z.shape
ax.plot_surface(X, Y, Z,  rstride=1, cstride=1, cmap=cm.hot,  shade=True,  linewidth=0)   #'classic', 'hot', 'pink', 'earth', 'heat', 'bone', 'coolwarm'
       
plt.show()
if HISTO:
    histo = py.hist(Z)
    title("Surface Profile")
    xlabel("heights")
    ylabel("occurences")
    show()