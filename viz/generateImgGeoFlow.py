import math
import time
import datetime
import numpy as np 
import matplotlib.pyplot as plt
#from mpl_toolkits.axes_grid1 import make_axes_locatable
import csv 
import os

# https://www.oreilly.com/library/view/python-data-science/9781491912126/ch04.html

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
hs = np.reshape(D[1:nx*ny+1,1],(ny,nx))


if not os.path.exists('./img'):
	os.makedirs('./img') 

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
	im = ax.imshow(h,extent=[0.0, np.amax(xc), 0.0, np.amax(yc)], origin='lower', cmap='viridis' , alpha=1.0, interpolation='bicubic',vmin=0.0,vmax=0.5)
	cb=fig.colorbar(im, orientation = 'horizontal',shrink=0.5,extend='max',pad=0.2,label=r'$h(x,y)$ [m]')
	im = ax.contour(xc,yc, np.transpose(h) , colors='black',linewidths=1.0,levels=[0.031,0.062,0.125,0.25,0.5])
	ax.clabel(im,inline=True, fontsize=5)
	fig.gca().set_aspect('equal', adjustable='box')
	plt.xlabel('Easting [m]')
	plt.ylabel('Northing [m]')	

	plt.title("$t_{\mathrm{e}}$ = "+str(time.strftime('%H:%M:%S',time.gmtime(t)))+" [s]")
	# save plot & reinit

	plt.savefig('./img/h_'+str(k).zfill(3)+'.png', dpi=300)
	cb.remove()
	plt.draw()
	ax.cla()
	print(" completion: "+str(round(k/nsave,2))+"\r")

	