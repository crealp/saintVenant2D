import math
import time
import datetime
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from matplotlib import cm

#from mpl_toolkits.axes_grid1 import make_axes_locatable
import csv 
import os

# python3 display.py
# latex 
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.xkcd()
# colorbar
cBtype0 = 'Greys_r'
cBtype  = 'Blues'

# import data for post-processing 
D    = np.genfromtxt('./dat/parameters.csv', delimiter=',')
nx   = int(D[1,0])
ny   = int(D[1,1])
dx   =    (D[1,2])
dy   =    (D[1,3])
nsave= int(D[1,6])
D    = np.genfromtxt('./dat/x.csv', delimiter=',')
xc   = D[1:nx*ny+1]

xc   = np.transpose(np.tile(D[1:nx*ny+1],(ny,1)))
D    = np.genfromtxt('./dat/y.csv', delimiter=',')
yc   = np.tile(D[1:nx*ny+1],(nx,1))
n    = "./dat/zhs.csv"
D    = np.genfromtxt(n, delimiter=',')
z    = np.reshape(D[1:nx*ny+1,0],(ny,nx))
hs   = np.reshape(D[1:nx*ny+1,1],(ny,nx))

Z    = np.zeros((ny,nx),dtype=float)
AS   = np.zeros((ny,nx),dtype=float)
HS   = np.zeros((ny,nx),dtype=float)
for i in range(1,ny-1):
	for j in range(1,nx-1):
		A = Z[j-1,i-1]
        B = Z[j-1,i]
        #C = Z[j-1,i+1]
        #D = Z[j,i-1]
        #F = Z[j,i+1]
        #G = Z[j+1,i-1]
        #H = Z[j+1,i]
        #I = Z[j+1,i+1]

        #dzx = ((C+2.0*F+I)-(A+2.0*D+G))/(8.0*dx)
        #dzy = ((G+2.0*H+I)-(A+2.0*B+C))/(8.0*dy)
	#for j in range(0,ny,1):
