# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 16:54:09 2016
Script to analyse the radius of 35 knots (R35)
@author: nederhof
"""

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
plot_data =0               # do part 3
# always do part 4: relationships
error_distribution = 0     # do part 5
# always do part 6: showing relationship

#==============================================================================
# Part 1: Loading data 
#==============================================================================
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

# Correlations
xcorrelation        = np.corrcoef(r35,x)
ycorrelation        = np.corrcoef(r35,y)
timecorrelation     = np.corrcoef(r35,time)
vmaxcorrelation     = np.corrcoef(r35,vmax)
dpccorrelation      = np.corrcoef(r35,dpc)
vtcorrelation       = np.corrcoef(r35,vt)
rmaxcorrelation      = np.corrcoef(r35,rmax)

#==============================================================================
# Part 2: necessity -> other relations result in large error + confidence bounds around coeffcients
#==============================================================================
if plot_rels ==1:
    plt.close('all')
    r35_holland1980a = r35_holland1980[~np.isnan(r35_holland1980)]
    r35a = r35[~np.isnan(r35_holland1980)]
    vmaxa = vmax[~np.isnan(r35_holland1980)]
    [rmse, bias, sci, mae, r2] = modelskill(r35a, r35_holland1980a)
    
    plt.figure(1009,figsize=(8,6))
    sc=plt.scatter(r35a, r35_holland1980a, s=50, c=vmaxa, alpha=0.5,lw = 0, cmap=plt.cm.YlOrRd)
    plt.xlabel('$\Delta AR_{35}$: based on TC database [km]')
    plt.ylabel('$\Delta AR_{35}$: Holland (1980, 2010) wind profile [km]')
    plt.legend(loc = 0)
    text1a =("RMSE: " ); text1b = (str(round(rmse,1)) +" km")
    text2a =("bias: " ); text2b = (str(round(bias,1)) +" km")
    text3a =("SCI: " ); text3b = (str(round(sci,3)))
    plt.text(340,145,text1a); plt.text(405,145,text1b)
    plt.text(340,125,text2a); plt.text(405,125,text2b)
    plt.text(340,105,text3a); plt.text(405,105,text3b)
    cb = plt.colorbar()
    cb.set_label('maximum sustained winds [m/s]')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.clim(20, 80); cb.set_ticks([20, 30, 40, 50, 60, 70, 80])
    plt.ylim(0,500);
    plt.xlim(0,500); 
    plt.grid(True)
    plt.plot([0,500],[0,500],'--k')
    plt.title('Average difference radius of 35 knots: observed \n versus predicted data based on Holland wind profile', fontweight='bold')
    fnamefig = ('r35_holland_vmax.png')
    plt.savefig(fnamefig, dpi=300)
    plt.close(1009)
    
    # Problem per vmax steps of 10 m/s
    nnumber=25;
    vsteps=np.linspace(20,80,nnumber)
    rms2b =[];MAE2=[];numbers=[]
    for i in range(0, nnumber-1):
        # Check value
        id1=vmaxa>=vsteps[i];
        id2=vmaxa<vsteps[i+1];
        id =id1 & id2;
        
        # Calculate MAE and RMSE
        try:
            rmaxtmp =r35a[id];
            rmaxVWtmp = r35_holland1980a[id];
            [rmse, bias, sci, mae,r2] = modelskill(rmaxtmp, rmaxVWtmp)
            rms2b.append(rmse)
            MAE2.append(mae)
        except:
            print 'error/exception'


#==============================================================================
# Part 3: find relations -> simply plotting
#==============================================================================
if plot_data ==1:
    
    # Plotten (vmax,r35,lat)
    plt.figure(1001,figsize=(8,6))
    area = 50
    sc=plt.scatter(vmax, r35, s=area, c=abs(y), alpha=0.5,lw = 0, cmap=plt.cm.YlOrRd)
    plt.xlabel('maximum sustained wind [m/s]')
    plt.ylabel('average difference radius of 35 knots [km]')
    plt.legend(loc = 0)
    cb = plt.colorbar()
    cb.set_label('latitude [$\circ$]')
    plt.clim(0, 60)
    plt.ylim(0,500); 
    plt.grid(True)
    plt.title('Relationships radius of 35 knots',fontweight='bold')
    fnamefig = ('vmax_r35_y.png')
    plt.savefig(fnamefig, dpi=300)
    plt.close(1001)
    
    
    # Plotten (vmax,r35,vt)
    plt.figure(1002,figsize=(8,6))
    area = 50
    sc=plt.scatter(vmax, r35, s=area, c=vt, alpha=0.5,lw = 0, cmap=plt.cm.YlOrRd)
    plt.xlabel('maximum sustained wind [m/s]')
    plt.ylabel('average difference radius of 35 knots [km]')
    plt.legend(loc = 0)
    cb = plt.colorbar()
    cb.set_label('propagation speed storms [m/s]')
    plt.clim(0, 25)
    plt.ylim(0,500); 
    plt.grid(True)
    plt.title('Relationships radius of 35 knots',fontweight='bold')
    fnamefig = ('vmax_r35_vt.png')
    plt.savefig(fnamefig, dpi=300)
    plt.close(1002)
    
    
    # Plotten (x,r35 and lat)
    plt.figure(1003,figsize=(8,6))
    area = 50
    sc=plt.scatter(x, r35, s=area, c=abs(y), alpha=0.5,lw = 0, cmap=plt.cm.YlOrRd)
    plt.xlabel('longitude [$\circ$]')
    plt.ylabel('average difference radius of 35 knots [km]')
    plt.legend(loc = 0)
    cb = plt.colorbar()
    cb.set_label('latitude [$\circ$]')
    plt.clim(0, 60)
    plt.ylim(0,500); 
    plt.grid(True)
    plt.title('Relationships radius of 35 knots',fontweight='bold')
    fnamefig = ('x_r35_y.png')
    plt.savefig(fnamefig, dpi=300)
    plt.close(1003)
    
    
    # Plotten (rmax,r35,y)
    plt.figure(1004,figsize=(8,6))
    area = 50
    sc=plt.scatter(rmax, r35, s=area, c=abs(y), alpha=0.5,lw = 0, cmap=plt.cm.YlOrRd)
    plt.xlabel('radius of maximum winds [km]')
    plt.ylabel('average difference radius of 35 knots [km]')
    plt.legend(loc = 0)
    cb = plt.colorbar()
    cb.set_label('latitude [$\circ$]')
    plt.clim(0, 60)
    plt.ylim(0,500); plt.xlim(0,250); 
    plt.grid(True)
    plt.title('Relationships radius of 35 knots',fontweight='bold')
    fnamefig = ('rmax_r35_y.png')
    plt.savefig(fnamefig, dpi=300)
    plt.close(1004)
    
    
    # Plotten (rmax,r35,y)
    plt.figure(1005,figsize=(8,6))
    area = 50
    sc=plt.scatter(vt, r35, s=area, c=time, alpha=0.5,lw = 0, cmap=plt.cm.YlOrRd)
    plt.xlabel('propagation speed storms [m/s]')
    plt.ylabel('average difference radius of 35 knots [km]')
    plt.legend(loc = 0)
    cb = plt.colorbar()
    cb.set_label('storm duration [hours]')
    plt.clim(0,300)
    plt.ylim(0,500); plt.xlim(0,30); 
    plt.grid(True)
    plt.title('Relationships radius of 35 knots',fontweight='bold')
    fnamefig = ('r35_vt_time.png')
    plt.savefig(fnamefig, dpi=300)
    plt.close(1005)
    
    
    # Plotten (rmax,r35,y)
    plt.figure(1006,figsize=(8,6))
    area = 50
    sc=plt.scatter(time, r35, s=area, c=vt, alpha=0.5,lw = 0, cmap=plt.cm.YlOrRd)
    plt.xlabel('storm duration [hours]')
    plt.ylabel('average difference radius of 35 knots [km]')
    plt.legend(loc = 0)
    cb = plt.colorbar()
    cb.set_label('propagation speed storms [m/s]')
    plt.clim(0, 30)
    plt.ylim(0,500); plt.xlim(0,500) 
    plt.grid(True)
    plt.title('Relationships radius of 35 knots',fontweight='bold')
    fnamefig = ('r35_time_vt.png')
    plt.savefig(fnamefig, dpi=300)
    plt.close(1006)
    
    
    # Plotten
    plt.figure(1007,figsize=(8,6))
    area = 50
    sc=plt.scatter(vt, r35, s=area, c=abs(y), alpha=0.5,lw = 0, cmap=plt.cm.YlOrRd)
    plt.xlabel('propogation speed tropical cyclone [m/s]')
    plt.ylabel('average difference radius of 35 knots [km]')
    plt.legend(loc = 0)
    cb = plt.colorbar()
    cb.set_label('latitude [$\circ$]')
    plt.clim(0, 60)
    plt.ylim(0,500); plt.xlim(0,30); 
    plt.grid(True)
    plt.title('Relationships radius of 35 knots',fontweight='bold')
    fnamefig = ('vt_r35_y.png')
    plt.savefig(fnamefig, dpi=300)
    plt.close(1007)
    plt.close('all')
    
    # Basemap
    plt.close(1008)
    plt.figure(1008,figsize=(10, 5))
    my_map=Basemap(projection='gall',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-360,urcrnrlon=0,resolution='i')
    xmap,ymap = my_map(x, y)
    my_map.scatter(xmap, ymap, s=5, c=r35, alpha=0.75,lw = 0, cmap=plt.cm.YlOrRd)
    my_map.drawcoastlines(linewidth=0.5)
    my_map.drawcountries()
    my_map.fillcontinents(color='lightgreen')
    my_map.drawmapboundary()
    my_map.drawmeridians(np.arange(0, 360, 60),labels=[1,1,0,1])
    my_map.drawparallels(np.arange(-90, 90, 30),labels=[1,0,0,1])
    plt.title('TC database: radius of 35 knots (R35)',fontweight='bold')
    cb = plt.colorbar()
    cb.set_label('average difference radius of 35 knots [km]')
    plt.clim(0, 300)
    my_map.drawmapscale(-90, -80, -90, -80, 5000, barstyle='simple', units='km', fontsize=9, yoffset=None, labelstyle='simple', fontcolor='k', fillcolor1='w', fillcolor2='k', ax=None, format='%d', zorder=None)
    fnamefig = ('world_r35.png')
    plt.savefig(fnamefig, dpi=300)
    plt.close(1008)

#==============================================================================
# Part 4: curve fitting + fit for error
#==============================================================================

# Save main errros
regions = 9;
totalsave=np.zeros(regions)
rmsesave=np.zeros(regions)
biassave=np.zeros(regions)
scisave=np.zeros(regions)

# Save coefficients
asave=np.zeros(regions)
bsave=np.zeros(regions)
csave=np.zeros(regions)
dsave=np.zeros(regions)
esave=np.zeros(regions)

# Largest and smallest error 
rmseAsave=np.zeros(regions)
rmseBsave=np.zeros(regions)

# Holland
rmsesave2=np.zeros(regions)

# Error function
aerror=np.zeros(regions)
berror=np.zeros(regions)

# Regionssss
for xxx in range(0,regions):
    
    # Specify which idregion is used!
    if      xxx == 0:                               # South Indian Ocean
        idregion = (x<-225) & (y < 0)
    elif    xxx == 1:                               # Indian Ocean
        idregion = (x<-260) & (y > 0)
    elif    xxx == 2:                               # Western Pacific
        idregion = (x>-260) & (x<-180) & (y > 0)
    elif    xxx == 3:                               # Central Pacific
        idregion = (x>-180) & (x<-140) & (y > 0)
    elif    xxx == 4:                               # Eastern Pacific Ocean
        idregion = (x>-140) & (AL < 0) & (y > 0)
    elif    xxx == 5:                               # South Pacific Ocean
        idregion = (x>-225) & (y < 0)
    elif    xxx == 6:                               # Atlantic Ocean
        idregion = (x>-140) & (AL>= 0)
    elif    xxx == 7:                               # Central + Eastern Pacific
        idregion = (x>-180) & (AL < 0) & (y > 0)
    elif    xxx == 8:                               # total!
        idregion = r35>0
        
    # Specific variables are being used!
    totalsave[xxx]=sum(idregion)
    X = np.array([vmax[idregion],y[idregion]])
    r35check = r35[idregion]
    vmaxcheck= vmax[idregion]

    # Relationship delta r35 to vmax
    def func(X,a,b,c):
        coriolis = abs(np.sin(X[1] * np.pi / 180))
        yvalue = a*np.power((X[0]-18),0.5)*(1+b*abs(X[1]))
        # some additional skill with progation speed -> not needed
        if np.all(yvalue<0):
            yvalue=0
        return yvalue
    
    # Determine fit + error of the fit
    nnumber=len(r35)
    sigmaval=np.ones(nnumber)
    popt, pcov = curve_fit(func,X,r35check, absolute_sigma=True,method='trf')
    R35_Nederhoff = func(X, *popt)
    [rmse, bias, sci, mae, r2] = modelskill(r35check, R35_Nederhoff)
    rmsesave[xxx]=rmse
    biassave[xxx]=bias
    scisave[xxx]=sci
    asave[xxx] = popt[0]
    bsave[xxx] = popt[1]

    # Check versus Holland
    r35_holland19802 = r35_holland1980[idregion]
    idnan   = ~np.isnan(r35_holland19802)
    [rmse, bias, sci, mae, r2] = modelskill(r35check[idnan], r35_holland19802[idnan])
    rmsesave2[xxx]=rmse
    
    # Problem per dpc steps of 10 m/s
    normalvalue= func(X, *popt)
    nnumber=25;
    vsteps=np.linspace(20,80,nnumber)
    rms2b =[];MAE2=[];numbers=[]
    
    for i in range(0, nnumber-1):
        # Check value
        id1=vmaxcheck>=vsteps[i];
        id2=vmaxcheck<vsteps[i+1];
        id =id1 & id2;
        
        # Calculate MAE and RMSE
        try:
            rmaxtmp =r35check[id];
            rmaxVWtmp = normalvalue[id];
            [rmse, bias, sci, mae,r2] = modelskill(rmaxtmp, rmaxVWtmp)
            rms2b.append(rmse)
            MAE2.append(mae)
        except:
            print 'error/exception'
    
    # Function for the error - option 1
    def func_error(X, a,b,c):
        yvalue = a+b*X+c*np.power(X,2)
        if np.all(yvalue<10):
            yvalue=10
        return yvalue

    # Function for the error - option 2
    xs = (0, 20, 35, 50, 80, 200); 
    ys = (30, 30, 70, 70, 30, 30); 
    func_error2 = interpolate.interp1d(xs, ys)
    errors = func_error2(vmax)
print rmsesave[-1]
print scisave[-1]
print biassave[-1]

# Determine fit + error of the fit
vstepsplot = vsteps[0:nnumber-1]+1.25
param_bounds=([0,-10,-10],[100,10,10])
popt_error, pcov_error = curve_fit(func_error,vstepsplot,rms2b,bounds=param_bounds)
perr_error = np.sqrt(np.diag(pcov_error))
vmaxval = np.linspace(20,80,1000)
errorvalue = func_error(vmaxval, *popt_error)
r35_ERROR = func_error(vmax, *popt_error)
formula="RMSE of $\Delta AR_{35}$ = a+b*$v_{max}$+c*$v_{max}^2$)"
fit1 = "a={:.1f}$\pm${:.1f}$, $b={:.1f}$\pm${:.1f}, c={:.3f}$\pm${:.3f}".format(popt_error[0], perr_error[0], popt_error[1], perr_error[1], popt_error[2], perr_error[2])


# Fit of error 1
fig, ax = plt.subplots(figsize=(15, 7.5))
rects1 = ax.bar(vstepsplot, rms2b, 2.5, color='r', alpha=0.5)
plt.plot(vmaxval, errorvalue)
plt.xlabel('maximum sustained winds [m/s]')
plt.ylabel('root-mean-square error in average difference radius of 35 knots')
plt.grid(True)
plt.text(51,75,formula); 
plt.text(51,72.5,fit1);
plt.xlim(20, 80)
plt.ylim(0, 80)
plt.title('Error distribution relationships average difference radius of 35 knots', fontweight='bold')
fnamefig = ('relationship_error.png')
plt.savefig(fnamefig, dpi=300)
plt.close('all')


# Fit of error 2
fig, ax = plt.subplots(figsize=(15, 7.5))
rects1 = ax.bar(vstepsplot, rms2b, 2.5, color='r', alpha=0.5)
plt.plot(vmaxval, func_error2(vmaxval))
plt.xlabel('maximum sustained winds [m/s]')
plt.ylabel('root-mean-square error in average difference radius of 35 knots')
plt.grid(True)
plt.xlim(20, 80)
plt.ylim(0, 80)
plt.title('Error distribution relationships average difference radius of 35 knots', fontweight='bold')
fnamefig = ('relationship_error2.png')
plt.savefig(fnamefig, dpi=300)
plt.close('all')
r35_ERROR = func_error2(vmaxcheck)

# Save results
import xlwt
book = xlwt.Workbook(encoding="utf-8")
sheet1 = book.add_sheet("Sheet 1")
i = 0
for n in range(0,len(asave)):
    sheet1.write(n, 0,asave[n])
    sheet1.write(n, 1,bsave[n])
    sheet1.write(n, 4,rmsesave[n])
    sheet1.write(n, 5,rmsesave2[n])
book.save("results_r35.xls")


#==============================================================================
# Part 5 Is the total error normal distributed?
#==============================================================================
if error_distribution == 1:

    X = np.array([vmax,y])
    R35_Nederhoff = func(X, *popt)
    idregion = (vmax > 32) & (abs(R35_Nederhoff-r35)<200)
    error = R35_Nederhoff[idregion] - r35[idregion];
    dagostino_results = stats.mstats.normaltest(error)
    ks_results = stats.kstest(error, cdf='norm')
    shapiro_results = stats.shapiro(error)
    
    ## QQplot
    plt.close(1010)
    plt.figure(1010,figsize=(15, 7.5))
    import pylab 
    stats.probplot(error, dist="norm", plot=pylab)
    
    pylab.show()
    pylab.show()
    plt.grid(True)
    fnamefig = ('qqplot.png')
    plt.savefig(fnamefig, dpi=300)
    plt.close('all')

#==============================================================================
# Part 6:  confidence band
#==============================================================================
# Method 3: http://codereview.stackexchange.com/questions/84414/obtaining-prediction-bands-for-regression-model
# https://www.coursera.org/learn/wharton-quantitative-modeling/lecture/Nndhc/4-4-r-squared-and-root-mean-squared-error-rmse
perr = np.sqrt(np.diag(pcov))
pp = (1. + 0.95) / 2.
nstd = stats.norm.ppf(pp)
tval    = nstd
popt_up = popt + tval * perr
popt_dw = popt - tval * perr
X = np.array([vmax, y])
uppervalue = func(X, *popt)+tval*rmsesave[-1]
lowervalue = func(X, *popt)-tval*rmsesave[-1]
normalvalue= func(X, *popt)
r35check    = r35
vmaxcheck   = vmax

#==============================================================================
# Part 7:  show solutions
#==============================================================================
plt.figure(1009,figsize=(8,6))
sc=plt.scatter(r35check, normalvalue, s=50, c=vmaxcheck, alpha=0.5,lw = 0, cmap=plt.cm.YlOrRd)
plt.xlabel('$\Delta AR_{35}$: based on TC database [km]')
plt.ylabel('$\Delta AR_{35}$: proposed [km]')
plt.legend(loc = 0)
text1a =("RMSE: " ); text1b = (str(round(rmsesave[-1],1)) +" km")
text2a =("bias: " ); text2b = (str(round(biassave[-1],1)) +" km")
text3a =("SCI: " ); text3b = (str(round(scisave[-1],3)))
plt.text(50,445,text1a); plt.text(105,445,text1b)
plt.text(50,425,text2a); plt.text(105,425,text2b)
plt.text(50,405,text3a); plt.text(105,405,text3b)
cb = plt.colorbar()
cb.set_label('maximum sustained winds [m/s]')
plt.gca().set_aspect('equal', adjustable='box')
plt.clim(20, 80); cb.set_ticks([20, 30, 40, 50, 60, 70, 80])
plt.ylim(0,500);
plt.xlim(0,500); 
plt.grid(True)
plt.plot([0,500],[0,500],'--k')
plt.title('Average difference radius of 35 knots: observed in data \n versus predicted by proposed relationship', fontweight='bold')
fnamefig = ('r35_nederhoff_vmax.png')
plt.savefig(fnamefig, dpi=300)
plt.close(1009)

# Define confidence interval.Convert to percentile point of the normal distribution. See: https://en.wikipedia.org/wiki/Standard_score and convert to number of standard deviations.
nnumber=20;
vsteps=np.linspace(20, 100,nnumber)
vstepsfound=[];rmax_min=[];rmax_max=[];rmax_mean=[];
for i in range(0, nnumber-1):
    # Check value
    id1=vmaxcheck>=vsteps[i];
    id2=vmaxcheck<vsteps[i+1];
    id =id1 & id2;
    # Calculate MAE and RMSE
    rmaxtmp1a =uppervalue[id];
    rmaxtmp1b =lowervalue[id];
    rmaxtmp1=np.concatenate((rmaxtmp1a, rmaxtmp1b), axis=0)
    try:
        rmax_max.append(np.max(rmaxtmp1))
        rmax_min.append(np.min(rmaxtmp1))
        rmax_mean.append(np.mean(normalvalue[id]))
        vstepsfound.append(vsteps[i])
    except:
        print 'error/exception'

print rmsesave[-1]
print scisave[-1]
print biassave[-1]

# Fit
fit1 = r"fit with: $a={:.3f}\pm{:.3f}$, $b={:.3f}\pm{:.3f}$".format(popt[0], perr[0], popt[1], perr[1])
formula=("$\Delta \overline{R_{35}} = a*(v_{max}-18).^{0.5} * (1+b*lat)$")
fit1 = "a={:.1f}, b={:.3f}".format(popt[0], popt[1])
plt.close(1007)
fig = plt.figure(1007,figsize=(8,6))
h1 = plt.plot(vstepsfound,rmax_mean,label='proposed relationship')
h2 = plt.fill_between(vstepsfound, rmax_max, rmax_min, facecolor='red', alpha=0.5, label='95% prediction interval')
h3 = plt.plot(vmax,r35,'k.',alpha=0.05, label='TC datapoints')
plt.xlim(20,80)
plt.ylim(0,500); 
plt.xlabel('maximum sustained winds [m/s]')
plt.ylabel('average difference radius of 35 knots [km]')
plt.grid(True)
plt.title('Proposed relationship including 95% prediction interval \n compared to data from TC database', fontweight='bold')
plt.legend(loc=0)
plt.text(60.1,75,formula)
plt.text(60.1,55,fit1)
fnamefig = ('relationship1.png')
plt.savefig(fnamefig, dpi=300)
plt.close(1007)

# Smooth version
# Show full solution space
nnumber= 100
vmax2 = np.linspace(20, 100,100)
lat2 = np.linspace(25, 25,100)
rmax2= np.ones((nnumber, nnumber))
rmax2U= np.ones((nnumber, nnumber))
rmax2D= np.ones((nnumber, nnumber))
error = np.ones((nnumber, nnumber))
for i in range(0, nnumber):
    for j in range(0, nnumber):
        rmax2[i,j]=func((vmax2[j],lat2[i]), *popt)
        rmax2U[i,j]=func((vmax2[j],lat2[i]), *popt) + tval*rmsesave[-1]
        rmax2D[i,j]=func((vmax2[j],lat2[i]), *popt) - tval*rmsesave[-1]
        error[i,j]=tval*func_error(vmax2[j], *popt_error)
[vmax2, lat2] = np.meshgrid(vmax2, lat2)

rmax3= np.ones((nnumber))
rmax3U= np.ones((nnumber))
rmax3D= np.ones((nnumber))
for i in range(0, nnumber):
    rvalues=np.concatenate((rmax2[:,i], rmax2U[:,i], rmax2D[:,i]))
    rmax3[i]=np.mean(rvalues)
    rmax3U[i]=np.max(rvalues)
    rmax3D[i]=np.min(rvalues)

idplot1 = (vmax>51) & (vmax < 75) & (r35 > 35) & (r35 < 75)
vmax2   = np.linspace(20, 100,100); lat2 = np.linspace(25, 25,100); X = (vmax2, lat2);
r35mean = func(X, *popt)
plt.close(1007)
fig = plt.figure(1007,figsize=(8,6))
h1 = plt.plot(vmax2,r35mean,label='proposed relationship')
h2 = plt.fill_between(vmax2, rmax3U, rmax3D, facecolor='red', alpha=0.5, label='95% prediction interval')
h3 = plt.plot(vmax[~idplot1],r35[~idplot1],'k.',alpha=0.05, label='TC datapoints')
plt.xlim(20,80)
plt.ylim(0,500); 
plt.xlabel('maximum sustained winds [m/s]')
plt.ylabel('average difference radius of 35 knots [km]')
plt.grid(True)
plt.title('Proposed relationship including 95% prediction interval \n compared to data from TC database', fontweight='bold')
plt.legend(loc=0)
plt.text(51.1,55,formula)
plt.text(51.1,35,fit1)
fnamefig = ('relationship2.png')
plt.savefig(fnamefig, dpi=300)
plt.close(1007)