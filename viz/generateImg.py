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
dx   = float(D[1,2])
nsave= int(D[1,6])
D    = np.genfromtxt('./dat/x.csv', delimiter=',')
xc   = D[1:nx*ny+1]

xc   = np.transpose(np.tile(D[1:nx*ny+1],(ny,1)))
D    = np.genfromtxt('./dat/y.csv', delimiter=',')
yc   = np.tile(D[1:nx*ny+1],(nx,1))
print(xc.shape)
print(yc.shape)

print('')
print('o---------------------------------------------o')
print('|               ** Plot data **               |')
print('o---------------------------------------------o')
print('generate fig & export to h_*.png')

#plt.ion()
n  = "./dat/zhs.csv"
D  = np.genfromtxt(n, delimiter=',')
z  = np.reshape(D[1:nx*ny+1,0],(ny,nx))
hs = np.reshape(D[1:nx*ny+1,1],(ny,nx))


if not os.path.exists('./img'):
	os.makedirs('./img') 

lvl = np.linspace(round(np.amin(z)),round(np.amax(z)),20)
lvl = np.arange(round(np.amin(z)),round(np.amax(z)),10)

fig0, ax0 = plt.subplots(figsize=(4,4)) 
im0 = ax0.imshow(z, extent=[0.0, np.amax(xc), 0.0, np.amax(yc)], cmap='gist_earth', alpha=1.0, interpolation='bicubic', vmin=z.min(), vmax=z.max())
fig0.gca().set_aspect('equal', adjustable='box')
plt.xlabel('Easting [m]')
plt.ylabel('Northing [m]')	
cb0=fig0.colorbar(im0, orientation = 'horizontal',extend='max',pad=0.2,label=r'$z(x,y)$ [m]',shrink=0.5)
im0 = ax0.contour(xc,yc, np.transpose(np.flip(z,axis=0)), levels=lvl, colors='black',linewidths=0.5)
ax0.clabel(im0,inline=True, fontsize=3.75)
plt.title('DTM, $\Delta_{x,y}=$ '+str(round(dx))+' [m]')
plt.savefig('./img/DTM.png', dpi=300)
cb0.remove()
plt.draw()
ax0.cla()

fig, ax = plt.subplots(figsize=(4,4)) 
for k in range(0,nsave+1,1):
	# load data
	name = "./dat/tdt_"+str(k)+".csv"
	D    = np.genfromtxt(name, delimiter=',')
	t    = D[1,0]
	name = "./dat/hQxQy_"+str(k)+".csv"
	D    = np.genfromtxt(name, delimiter=',')
	h    = np.reshape(D[1:nx*ny+1,0],(ny,nx))
	# plot data
	im = ax.imshow(hs, extent=[0.0, np.amax(xc), 0.0, np.amax(yc)], cmap=cBtype0, alpha=1.0, interpolation='bicubic')
	im = ax.imshow(h , extent=[0.0, np.amax(xc), 0.0, np.amax(yc)], cmap=cBtype , alpha=0.5, interpolation='bicubic', norm=colors.LogNorm(vmin=1e-4, vmax=1e-2))
	#im = ax.imshow(h , extent=[0.0, np.amax(xc), 0.0, np.amax(yc)], cmap=cBtype , alpha=0.5, interpolation='bicubic', vmin=1e-4, vmax=1e-2)
	fig.gca().set_aspect('equal', adjustable='box')
	plt.xlabel('Easting [m]')
	plt.ylabel('Northing [m]')	
	cb=fig.colorbar(im, orientation = 'horizontal',extend='max',pad=0.2,label=r'$h(x,y)$ [m]',shrink=0.5)
	plt.title("$t_{\mathrm{e}}$ = "+str(time.strftime('%H:%M:%S',time.gmtime(t)))+" [s]")
	# save plot & reinit

	plt.savefig('./img/h_'+str(k).zfill(3)+'.png', dpi=300)
	cb.remove()
	plt.draw()
	ax.cla()
	print(" completion: "+str(round(k/nsave,2))+"\r")
	