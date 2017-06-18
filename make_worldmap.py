# -*- coding: utf-8 -*-
"""
Created on Wed Jun 07 13:45:12 2017
@author: nederhof
"""

# Make a map!
# Libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from mpl_toolkits.basemap import Basemap, addcyclic
from math import sqrt
from hurricane_functions import *
from scipy import stats
from scipy.stats.stats import pearsonr
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats.distributions import  t
from scipy import interpolate
import cmocean

# Change working directory
import os
os.chdir('/papers/2017/probabilistic_forecasting_tropical_cyclones/figures_r35/')
os.getcwd()

#==============================================================================
# Part 0: settings of this script
#==============================================================================
# always do part 1: loading data
plot_rels = 1               # do part 2
plot_data = 1               # do part 3
# always do part 4: relationships
error_distribution = 1      # do part 5
# always do part 6: showing relationship

#==============================================================================
# Part 1: Loading data 
#==============================================================================
## Retrieve r35
fname = (r'd:\data\hurricanes\analysis\hurricane_vmax_r35.csv')
f = open(fname, 'r')
time=[]
x=[]
y=[]
vt =[]
vmax = []
pc=[]
r35 = []
rmax=[]
dpc=[]
r35_holland1980=[]
AL =[];
for line in f:
    cells = line.split(',')
    time.append(float(cells[8]))
    x.append(float(cells[1]))
    y.append(float(cells[2]))
    vt.append(float(cells[3]))
    vmax.append(float(cells[4]))
    pc.append(float(cells[5]))
    rmax.append(float(cells[6]))
    r35.append(float(cells[7]))
    dpc.append(float(cells[10]))
    r35_holland1980.append(float(cells[9]))
    AL.append(float(cells[12]))
f.close()

# Convert to numpy arrays
time = np.asarray(time)
x = np.asarray(x)
y = np.asarray(y)
vt = np.asarray(vt)
vmax = np.asarray(vmax)
pc = np.asarray(pc)
rmax = np.asarray(rmax)
r35 = np.asarray(r35)
dpc = np.asarray(dpc)
r35_holland1980 = np.asarray(r35_holland1980)
AL = np.asarray(AL)

# R35 is actually delta R35
r35 = r35-rmax
r35_holland1980 = r35_holland1980-rmax
xr35= x
yr35 = y

## Retrieve rmax
fname = (r'd:\data\hurricanes\analysis\hurricane_pc_rmax.csv')
f = open(fname, 'r')
time=[]
x=[]
y=[]
vt =[]
vmax = []
pc=[]
rmax=[]
dpc=[]
AL =[];
for line in f:
    cells = line.split(',')
    x.append(float(cells[1]))
    y.append(float(cells[2]))
    vt.append(float(cells[3]))
    vmax.append(float(cells[4]))
    pc.append(float(cells[5]))
    rmax.append(float(cells[6]))
    time.append(float(cells[7]))
    dpc.append(float(cells[8]))
    AL.append(float(cells[9]))
f.close()

# Convert to numpy arrays
time = np.asarray(time)
x = np.asarray(x)
y = np.asarray(y)
vt = np.asarray(vt)
vmax = np.asarray(vmax)
pc = np.asarray(pc)
rmax = np.asarray(rmax)
dpc = np.asarray(dpc)
AL = np.asarray(AL)
id25dpc = (rmax>0) & (rmax< 100) & (vmax > 0)

# Correlations
xcorrelation        = np.corrcoef(rmax,x)
ycorrelation        = np.corrcoef(rmax,y)
timecorrelation     = np.corrcoef(rmax,time)
vmaxcorrelation     = np.corrcoef(rmax,vmax)
dpccorrelation      = np.corrcoef(rmax,dpc)
vtcorrelation       = np.corrcoef(rmax,vt)



## Make basemap
from matplotlib.patches import Polygon


plt.close('all')
fig, axes = plt.subplots(3, 1, figsize=(10, 15))

axes[0].set_title("Maximum wind speeds and oceanic basins", fontweight='bold')
my_map0 = Basemap(projection='gall',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-360,urcrnrlon=0,ax=axes[0])

#my_map0.drawmapboundary(fill_color='blue')
my_map0.fillcontinents(color='lightgreen',lake_color='lightgreen')
my_map0.drawcoastlines()
my_map0.drawmeridians(np.arange(0, 360, 60),labels=[1,1,0,1])
my_map0.drawparallels(np.arange(-90, 90, 30),labels=[1,0,0,1])

# Wind speeds
xmap,ymap = my_map0(x, y)
fig00 = my_map0.scatter(xmap, ymap, s=5, c=vmax, alpha=0.75,lw = 0,cmap=plt.cm.Reds, vmin=10, vmax=70)
cb = my_map0.colorbar(fig00, extend='both')
cb.set_label('maximum sustainted winds (vmax) [m/s]')

# Colors
colormap = plt.cm.Dark2 
line_colors = colormap(np.linspace(1,1,6))

# Nothern Indian Ocean
x1,y1 = my_map0(-360,0)
x2,y2 = my_map0(-260,0)
x3,y3 = my_map0(-260,45)
x4,y4 = my_map0(-360,45)
polya = Polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4)],facecolor='k',edgecolor='k',linewidth=2, alpha=0.5)
plt.gca().add_patch(polya)

# North
# Western Pacific
x1,y1 = my_map0(-260,0)
x2,y2 = my_map0(-180,0)
x3,y3 = my_map0(-180,60)
x4,y4 = my_map0(-260,60)
polyb = Polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4)],facecolor='k',edgecolor='k',linewidth=2, alpha=0.5)
plt.gca().add_patch(polyb)

# Eastern Pacific
x1,y1 = my_map0(-180,0)
x2,y2 = my_map0(-75,0)
x3,y3 = my_map0(-125,45)
x4,y4 = my_map0(-125,60)
x5,y5 = my_map0(-180,60)
polyc = Polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4), (x5, y5)],facecolor='k',edgecolor='k',linewidth=2, alpha=0.5)
plt.gca().add_patch(polyc)

# Atlantic Ocean
x1,y1 = my_map0(0,0)
x2,y2 = my_map0(-75,0)
x3,y3 = my_map0(-125,45)
x3b,y3b = my_map0(-125,60)
x4,y4 = my_map0(0,60)
polyd = Polygon([(x1,y1),(x2,y2),(x3,y3),(x3b,y3b),(x4,y4)],facecolor='k',edgecolor='k',linewidth=2, alpha=0.5)
plt.gca().add_patch(polyd)

# Southern hemisphere
# Indian Ocean
x1,y1 = my_map0(-360,0)
x2,y2 = my_map0(-225,0)
x3,y3 = my_map0(-225,-45)
x4,y4 = my_map0(-360,-45)
polye = Polygon([(x1,y1),(x2,y2),(x3,y3),(x4,y4)],facecolor='k',edgecolor='k',linewidth=2, alpha=0.5)
plt.gca().add_patch(polye)

# South Pacific
x1c,y1c = my_map0(-225,-0)
x2c,y2c = my_map0(-60,-0)
x3c,y3c = my_map0(-60,-45)
x4c,y4c = my_map0(-225,-45)
polyc = Polygon([(x1c,y1c),(x2c,y2c),(x3c,y3c),(x4c,y4c)],facecolor='k',edgecolor='k',linewidth=2, alpha=0.5)
plt.gca().add_patch(polyc)


# Names
lat = [25, 25, 25, 25, -25, -25]#, 30, 30, -30, -30, -30]
lon = [-320, -230, -160, -90, -355, -200]#, 30 ]
X,Y = my_map0(lon,lat)
labels = ['Northern', 'Western', 'Eastern', 'Atlantic Ocean', 'Southern Indian Ocean', 'Southern Pacific Ocean']
for label, xpt, ypt in zip(labels, X,Y):
    plt.text(xpt, ypt, label, color='0.75')


axes[1].set_title("Inner-core wind structure", fontweight='bold')
my_map1 = Basemap(projection='gall',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-360,urcrnrlon=0,ax=axes[1])
xmap,ymap = my_map1(x, y)
fig1 = my_map1.scatter(xmap, ymap, s=5, c=rmax, alpha=0.75,lw = 0, cmap=plt.cm.OrRd, vmin=0, vmax=100)
#my_map1.drawmapboundary(fill_color='blue')
my_map1.fillcontinents(color='lightgreen',lake_color='lightgreen')
my_map1.drawcoastlines()
my_map1.drawmeridians(np.arange(0, 360, 60),labels=[1,1,0,1])
my_map1.drawparallels(np.arange(-90, 90, 30),labels=[1,0,0,1])
cb = my_map1.colorbar(fig1, extend='max')
cb.set_label('radius of maximum winds (RMW) [km]')

axes[2].set_title("Outer wind structure", fontweight='bold')
my_map2 = Basemap(projection='gall',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-360,urcrnrlon=0,ax=axes[2])
xmap,ymap = my_map2(xr35, yr35)
fig2 = my_map2.scatter(xmap, ymap, s=5, c=r35, alpha=0.75,lw = 0, cmap=plt.cm.YlOrRd, vmin=0, vmax=250)
#my_map2.drawmapboundary(fill_color='blue')
my_map2.fillcontinents(color='lightgreen',lake_color='lightgreen')
my_map2.drawcoastlines()
my_map2.drawmeridians(np.arange(0, 360, 60),labels=[1,1,0,1])
my_map2.drawparallels(np.arange(-90, 90, 30),labels=[1,0,0,1])
cb = my_map2.colorbar(fig2, extend='max')
cb.set_label('average difference radius of 35 knots [km]')

my_map2.drawmapscale(-90, -80, -90, -80, 5000, barstyle='simple', units='km', fontsize=9, yoffset=None, labelstyle='simple', fontcolor='k', fillcolor1='w', fillcolor2='k', ax=None, format='%d', zorder=None)

plt.show()
fnamefig = 'world2.png'
plt.savefig(fnamefig, dpi=300)
plt.close('all')
