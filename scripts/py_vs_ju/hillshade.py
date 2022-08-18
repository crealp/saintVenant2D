import math
from math import cos, sin, atan, atan2, degrees, sqrt, pi
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

Z    = z

tic = time.perf_counter()
AS   = np.zeros((ny,nx),dtype=float)
HS   = np.zeros((ny,nx),dtype=float)
phi  = 45.0*pi/180
theta= 315.0*pi/180
for i in range(1,ny-1):
    for j in range(1,nx-1):
        A = Z[i-1,j-1]
        B = Z[i-0,j-1]
        C = Z[i+1,j-1]

        D = Z[i-1,j]
        F = Z[i+1,j]
        
        G = Z[i-1,j+1]
        H = Z[i,j+1]
        I = Z[i+1,j+1]

        dzx = ((C+2.0*F+I)-(A+2.0*D+G))/(8.0*dx)
        dzy = ((G+2.0*H+I)-(A+2.0*B+C))/(8.0*dy)

        s   = (atan(( sqrt(dzx*dzx+dzy*dzy) ) ))
        a   = (atan2(dzy,-dzx))

        if a<pi/2:
            AS[i,j]=-a+pi/2
        elif a>=pi/2:
            AS[i,j]=2*pi-a+pi/2
        h = 255.0*( (cos(phi)*cos(s))+(sin(phi)*sin(s*cos(theta-AS[i,j]))) )
        if h>=0.0:
            HS[i,j]=abs(h+1.0)
        else:
            HS[i,j]=1.0

toc = time.perf_counter()
print(f"Downloaded the tutorial in {toc - tic:0.4f} seconds")

fig, ax = plt.subplots(figsize=(4,4)) 
im = ax.imshow(hs, cmap=cBtype0, alpha=1.0, interpolation='bicubic')
fig.gca().set_aspect('equal', adjustable='box')
plt.savefig('python_hillshade.png', dpi=300, bbox_inches='tight')