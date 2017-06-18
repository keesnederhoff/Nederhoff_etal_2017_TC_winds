# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 16:54:09 2016
Script to analyse the radius of maximum winds (rmax)
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
import scipy
import scipy.stats
from scipy.stats.kde import gaussian_kde

# Change working directory
import os
os.chdir('/papers/2017/probabilistic_forecasting_tropical_cyclones/figures_rmax/')
os.getcwd()

#==============================================================================
# Part 0: settings of this script
#==============================================================================
# always do part 1: loading data
plot_rels = 1               # do part 2
plot_data = 0               # do part 3
# always do part 4: relationships
error_distribution = 0      # do part 5
# always do part 6: showing relationship

#==============================================================================
# Part 1: Loading data 
#==============================================================================
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

#==============================================================================
# Part 2: necessity -> other relations result in large error + prediction bounds around coeffcients
#==============================================================================
if plot_rels ==1:

    # Find Vickery and Wadhera value + error
    RMW_VW=vickery_wadhera(dpc,abs(y))
    [rmse, bias, sci, mae,r2] = modelskill(rmax[id25dpc], RMW_VW[id25dpc])
    plt.figure(1001,figsize=(8,6))
    sc=plt.scatter(rmax, RMW_VW, s=50, c=vmax, alpha=0.5,lw = 0, cmap=plt.cm.YlOrRd)
    plt.xlabel('$R_{max}$: based on TC database [km]')
    plt.ylabel('$R_{max}$: based on Vickery & Wadhera (2008) [km]')
    plt.legend(loc = 0)
    text1a =("RMSE: " ); text1b = (str(round(rmse,1)) +" km")
    text2a =("bias: " ); text2b = (str(round(bias,1)) +" km")
    text3a =("SCI: " ); text3b = (str(round(sci,3)) +" [-]")
    plt.text(70,11,text1a); plt.text(82,11,text1b)
    plt.text(70,8,text2a); plt.text(82,8,text2b)
    plt.text(70,5,text3a); plt.text(82,5,text3b)
    cb = plt.colorbar()
    cb.set_label('$v_{max} [m/s]$')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.clim(0, 80)
    plt.ylim(0,100);
    plt.xlim(0,100); 
    plt.grid(True)
    plt.plot([0,100],[0,100],'--k')
    plt.title('Radius of maximum winds: observed \n versus predicted data based on Vickery & Wadhera (2008)', fontweight='bold')
    fnamefig = 'rmax_vickery_vmax.png'
    plt.savefig(fnamefig, dpi=300)
    plt.close(1001)
    
    # Find Knaff and Zehr (2007) + error
    vmaxknots = vmax/0.51444
    RMW_KZ = knaff_zehr(vmaxknots,abs(y))
    [rmse, bias, sci, mae,r2] = modelskill(rmax[id25dpc], RMW_KZ[id25dpc])
    plt.figure(1001,figsize=(8,6))
    sc=plt.scatter(rmax, RMW_KZ, s=50, c=vmax, alpha=0.5,lw = 0, cmap=plt.cm.YlOrRd)
    plt.xlabel('$R_{max}$: based on TC database [km]')
    plt.ylabel('$R_{max}$: based on Knaff and Zehr (2008) [km]')
    plt.legend(loc = 0)
    text1a =("RMSE: " ); text1b = (str(round(rmse,1)) +" km")
    text2a =("bias: " ); text2b = (str(round(bias,1)) +" km")
    text3a =("SCI: " ); text3b = (str(round(sci,3)) +" [-]")
    plt.text(70,11,text1a); plt.text(82,11,text1b)
    plt.text(70,8,text2a); plt.text(82,8,text2b)
    plt.text(70,5,text3a); plt.text(82,5,text3b)
    cb = plt.colorbar()
    cb.set_label('$v_{max} [m/s]$')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.clim(0, 80)
    plt.ylim(0,100);
    plt.xlim(0,100); 
    plt.grid(True)
    plt.plot([0,100],[0,100],'--k')
    plt.title('Radius of maximum winds: observed in data \n versus predicted by Knaff and Zehr (2008)', fontweight='bold')
    fnamefig = ('rmax_knaff_vmax.png')
    plt.savefig(fnamefig, dpi=300)
    plt.close(1001)
    
    # Find Gross (2004) + error
    vmaxknots = vmax/0.51444
    RMW_gross = gross(vmaxknots,abs(y))
    [rmse, bias, sci, mae,r2] = modelskill(rmax[id25dpc], RMW_gross[id25dpc])
    plt.figure(1001,figsize=(8,6))
    sc=plt.scatter(rmax, RMW_gross, s=50, c=vmax, alpha=0.5,lw = 0, cmap=plt.cm.YlOrRd)
    plt.xlabel('$R_{max}$: based on TC database [km]')
    plt.ylabel('$R_{max}$: based on Gross (2004) [km]')
    plt.legend(loc = 0)
    text1a =("RMSE: " ); text1b = (str(round(rmse,1)) +" km")
    text2a =("bias: " ); text2b = (str(round(bias,1)) +" km")
    text3a =("SCI: " ); text3b = (str(round(sci,3)) +" [-]")
    plt.text(70,11,text1a); plt.text(82,11,text1b)
    plt.text(70,8,text2a); plt.text(82,8,text2b)
    plt.text(70,5,text3a); plt.text(82,5,text3b)
    cb = plt.colorbar()
    cb.set_label('$v_{max} [m/s]$')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.clim(0, 80)
    plt.ylim(0,100);
    plt.xlim(0,100); 
    plt.grid(True)
    plt.plot([0,100],[0,100],'--k')
    plt.title('Radius of maximum winds: observed in data \n versus predicted by Gross (2004))', fontweight='bold')
    fnamefig = 'rmax_gross_vmax.png'
    plt.savefig(fnamefig, dpi=300)
    plt.close(1001)
    
    # Problem per vmax steps of 10 m/s
    nnumber=14;
    vsteps=np.linspace(0,130,nnumber)
    rms2a1 =[];biasa1=[];stda1=[]
    rms2a2 =[];biasa2=[];stda2=[]
    rms2a3 =[];biasa3=[];stda3=[]
    
    for i in range(0, nnumber-1):
        # Check value
        id1=dpc>=vsteps[i];
        id2=dpc<vsteps[i+1];
        id =id1 & id2;
        # Calculate MAE and RMSE
        rmaxtmp =rmax[id];
        rmaxVWtmp = RMW_VW[id];
        rmaxGrosstmp = RMW_gross[id];
        rmaxKZtmp = RMW_KZ[id];
        [rmse, bias, sci, mae,r2] = modelskill(rmaxtmp, rmaxVWtmp)
        rms2a1.append(rmse)
        biasa1.append(bias)
        [rmse, bias, sci, mae,r2] = modelskill(rmaxtmp, rmaxGrosstmp)
        rms2a2.append(rmse)
        biasa2.append(bias)
        [rmse, bias, sci, mae,r2] = modelskill(rmaxtmp, rmaxKZtmp)
        rms2a3.append(rmse)
        biasa3.append(bias)

#==============================================================================
# Part 3: find relations -> simply plotting
#==============================================================================
if plot_data ==1:

    # Plotten (vmax,rmax,lat)
    plt.figure(1002,figsize=(8,6))
    sc=plt.scatter(vmax, rmax, s=50, c=abs(y), alpha=0.5,lw = 0, cmap=plt.cm.YlOrRd)
    plt.xlabel('maximum sustained winds [m/s]')
    plt.ylabel('radius of maximum winds [km]')
    plt.legend(loc = 0)
    cb = plt.colorbar()
    cb.set_label('latitude [$\circ$]')
    plt.clim(0, 60)
    plt.ylim(0,100); 
    plt.grid(True)
    plt.title('Relationships radius of maximum winds', fontweight='bold')
    fnamefig = 'vmax_rmax_y.png'
    plt.savefig(fnamefig, dpi=300)
    plt.close(1002)
    
    # Plotten (vmax,rmax,time)
    plt.figure(1003,figsize=(8,6))
    sc=plt.scatter(vmax, rmax, s=50, c=time, alpha=0.5,lw = 0, cmap=plt.cm.YlOrRd)
    plt.xlabel('maximum sustained winds [m/s]')
    plt.ylabel('radius of maximum winds [km]')
    plt.legend(loc = 0)
    cb = plt.colorbar()
    cb.set_label('storm duration / age [hours]')
    plt.clim(0, 300)
    plt.ylim(0,100); 
    plt.grid(True)
    plt.title('Relationships radius of maximum winds', fontweight='bold')
    fnamefig = 'vmax_rmax_time.png'
    plt.savefig(fnamefig, dpi=300)
    plt.close(1003)
    
    # Plotten (dpc,rmax,vt)
    plt.figure(1004,figsize=(8,6))
    area = 50
    sc=plt.scatter(dpc, rmax, s=area, c=vt, alpha=0.5,lw = 0, cmap=plt.cm.YlOrRd)
    plt.xlabel('pressure drop in the eye of the storm [hPa]')
    plt.ylabel('radius of maximum winds [km]')
    plt.legend(loc = 0)
    cb = plt.colorbar()
    cb.set_label('propagation speed storms [m/s]')
    plt.clim(0,15); plt.xlim(0, 120)
    plt.ylim(0,100); 
    plt.grid(True)
    plt.title('Relationships radius of maximum winds', fontweight='bold')
    fnamefig = 'dpc_rmax_vt.png'
    plt.savefig(fnamefig, dpi=300)
    plt.close(1004)
    
    # Plotten (dpc,rmax,vt)
    plt.figure(1004,figsize=(8,6))
    area = 50
    sc=plt.scatter(dpc, rmax, s=area, c=y, alpha=0.5,lw = 0, cmap=plt.cm.YlOrRd)
    plt.xlabel('pressure drop in the eye of the storm [hPa]')
    plt.ylabel('radius of maximum winds [km]')
    plt.legend(loc = 0)
    cb = plt.colorbar()
    cb.set_label('latitude [$\circ$]')
    plt.clim(0,45); plt.xlim(0, 120)
    plt.ylim(0,100); 
    plt.grid(True)
    plt.title('Relationships radius of maximum winds', fontweight='bold')
    fnamefig = 'dpc_rmax_lat.png'
    plt.savefig(fnamefig, dpi=300)
    plt.close(1004)
    
    # Basemap
    plt.close(1010)
    plt.figure(1008,figsize=(10, 5))
    my_map=Basemap(projection='gall',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-360,urcrnrlon=0,resolution='i')
    xmap,ymap = my_map(x, y)
    my_map.scatter(xmap, ymap, s=5, c=rmax, alpha=0.75,lw = 0, cmap=plt.cm.YlOrRd)
    my_map.drawcoastlines(linewidth=0.5)
    my_map.drawcountries()
    my_map.fillcontinents(color='lightgreen')
    my_map.drawmapboundary()
    my_map.drawmeridians(np.arange(0, 360, 60),labels=[1,1,0,1])
    my_map.drawparallels(np.arange(-90, 90, 30),labels=[1,0,0,1])
    plt.title('NHC & JTWC database: radius of maximum winds (RMW)',fontweight='bold')
    cb = plt.colorbar()
    cb.set_label('radius of maximum winds [km]')
    plt.clim(0, 100)
    my_map.drawmapscale(-90, -80, -90, -80, 5000, barstyle='simple', units='km', fontsize=9, yoffset=None, labelstyle='simple', fontcolor='k', fillcolor1='w', fillcolor2='k', ax=None, format='%d', zorder=None)
    fnamefig = 'world_rmax.png'
    plt.savefig(fnamefig, dpi=300)
    plt.close(1010)

#==============================================================================
# Part 4: curve fitting + fit for error -> this is a loop!
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

# Error function
aerror=np.zeros(regions)
berror=np.zeros(regions)

# Vickery save
rmsesave2 =np.zeros(regions)

# Regionss
for xxx in range(0,regions):
    
    # Specify which idregion is used!
    if      xxx == 0:                               # South Indian Ocean
        idregion = (x<-225) & (y < 0) & (rmax>0) & (rmax< 100)
    elif    xxx == 1:                               # Indian Ocean
        idregion = (x<-260) & (y > 0) & (rmax>0) & (rmax< 100)
    elif    xxx == 2:                               # Western Pacific
        idregion = (x>-260) & (x<-180) & (y > 0) & (rmax>0) & (rmax< 100)
    elif    xxx == 3:                               # Central Pacific
        idregion = (x>-180) & (x<-140) & (y > 0) & (rmax>0) & (rmax< 100)
    elif    xxx == 4:                               # Eastern Pacific Ocean
        idregion = (x>-140) & (AL < 0) & (y > 0) & (rmax>0) & (rmax< 100)
    elif    xxx == 5:                               # South Pacific Ocean
        idregion = (x>-225) & (y < 0) & (rmax>0) & (rmax< 100)
    elif    xxx == 6:                               # Atlantic Ocean
        idregion = (x>-140) & (AL> 0) & (rmax>0) & (rmax< 100)
    elif    xxx == 7:                               # Central + Eastern Pacific
        idregion = (x>-180) & (AL < 0) & (y > 0) & (rmax>0) & (rmax< 100)
    elif    xxx == 8:                               # total!
        idregion = (rmax>0) & (rmax< 100)
    
    # Specific variables are being used!
    totalsave[xxx]=sum(idregion)
    
    # New attempt
    value1      = 0;
    value2      = 33;
    def func(X, a, b,c):
        yvalue1 = a + (0*X[0]**1)
        yvalue2 = b + (c*X[0]**1)
        ratio   = (X[0]-value1)/(value2-value1)
        yvalue3 = (1-ratio)*yvalue1 + ratio*yvalue2
        yvalue  = yvalue2
        
        # Small
        id =X[0]<=value1
        yvalue[id]  = yvalue1[id]
    
        # Intermediate
        id = (X[0]>value1) & (X[0]<=value2)
        yvalue[id]  = yvalue3[id]
    
        # Large
        id = X[0]>=value2
        yvalue[id]  = yvalue2[id]
        return yvalue
    
    # Determine fit + error of the fit - ALL
    X = np.array([vmax[idregion],y[idregion]])
    rmaxcheck = rmax[idregion]
    popt, pcov = curve_fit(func,X,rmaxcheck)
    RMW_Nederhoff = func(X, *popt)
    [rmse, bias, sci, mae,r2] = modelskill(rmaxcheck, RMW_Nederhoff)
    asave[xxx]=popt[0]
    bsave[xxx]=popt[1]
    csave[xxx]=popt[2]

    # Determine fit for only TC
    idtc = idregion & (vmax > 33)
    X = np.array([vmax[idtc],y[idtc]])
    rmaxcheck = rmax[idtc]
    popt, pcov = curve_fit(func,X,rmaxcheck)
    RMW_Nederhoff = func(X, *popt)
    [rmse, bias, sci, mae,r2] = modelskill(rmaxcheck, RMW_Nederhoff)
    bsave[xxx]=popt[1]
    csave[xxx]=popt[2]

    # Combine ->  full data set for plotting purposes + part for skill
    X = np.array([vmax[idregion],y[idregion]])
    rmaxcheck = rmax[idregion]
    vmaxcheck = vmax[idregion]
    idreduced = vmaxcheck > 0
    popt[0]=asave[xxx]
    popt[1]=bsave[xxx]
    popt[2]=csave[xxx]
    RMW_Nederhoff = func(X, *popt)
    
    # Save skill
    [rmse, bias, sci, mae, r2] = modelskill(rmaxcheck[idreduced], RMW_Nederhoff[idreduced])
    rmsesave[xxx]=rmse
    biassave[xxx]=bias
    scisave[xxx]=sci

    # Vickery  skill?
    RMW_VW2 = RMW_VW[idregion]
    [rmse, bias, sci, mae,r2] = modelskill(rmaxcheck, RMW_VW2)
    rmsesave2[xxx]=rmse
    
    # Problem per vmax steps of 5 m/s
    normalvalue= func(X, *popt)
    nnumber=14;
    vsteps=np.linspace(10,80,nnumber)
    rms2b =[];bias2b=[];numbers=[];vstepssaved=[];
    for i in range(0, nnumber-1):
        
        # Check value
        try:
            id1=vmaxcheck>=vsteps[i];
            id2=vmaxcheck<vsteps[i+1];
            id =id1 & id2;
            
            # Calculate MAE and RMSE
            rmaxtmp =rmaxcheck[id];
            rmaxVWtmp = normalvalue[id];
            [rmse, bias, sci, mae,r2] = modelskill(rmaxtmp, rmaxVWtmp)
            rms2b.append(rmse)
            bias2b.append(bias)
            vstepssaved.append(vsteps[i]+2.5)
        except:
            print "oops"
    rmseAsave[xxx]=max(rms2b)
    rmseBsave[xxx]=min(rms2b)

##  Plot the error
#fig, ax = plt.subplots(figsize=(15, 7.5))
#rects1 = ax.bar(vstepssaved, rms2b, 5, color='r', alpha=0.5)
#plt.xlabel('maximum sustained winds [m/s]')
#plt.ylabel('root-mean-square error in radius of maximum winds [km]')
#plt.grid(True)
#plt.text(101,45,formula); 
#plt.text(101,42.5,fit1);
#plt.title('Error distribution relationships radius of maximum winds', fontweight='bold')
#fnamefig = 'relationship_error.png'
#plt.savefig(fnamefig, dpi=300)
#plt.close('all')

# Save results
import xlwt
book = xlwt.Workbook(encoding="utf-8")
sheet1 = book.add_sheet("Sheet 1")
i = 0
for n in range(0,len(asave)):
    sheet1.write(n, 0,asave[n])
    sheet1.write(n, 1,bsave[n])
    sheet1.write(n, 2,csave[n])
    sheet1.write(n, 3,rmsesave[n])
    sheet1.write(n,4,totalsave[n])
    sheet1.write(n,5,rmsesave2[n])
book.save("results_rmax.xls")

##==============================================================================
## Part 5 Is the total error normal distributed?
##==============================================================================
#if error_distribution == 1:
#    plt.close(1010)
#    plt.figure(1010,figsize=(8,6))
#    idregion = (vmax > 0)
#    error = RMW_Nederhoff - rmaxcheck
#    print len(error)
#    anderson1_results = stats.anderson(error)
#    dagostino_results = stats.mstats.normaltest(error)
#    ks_results = stats.kstest(error, cdf='norm')
#    shapiro_results = stats.shapiro(error)
#    xvalues =np.linspace(-100,100,50)
#    h = plt.hist(error, xvalues, color='k')
#    dist_names = ['norm', 'cauchy', 't', 'rayleigh'];
#    size = len(error)
#    
#    for dist_name in dist_names:
#        dist = getattr(stats, dist_name)
#        param = dist.fit(error)
#        pdf_fitted = dist.pdf(xvalues, *param[:-2], loc=param[-2], scale=param[-1])*size
#        plt.plot(xvalues,pdf_fitted, label=dist_name)
#        plt.xlim(-100,+100)
#    plt.legend(loc='upper right')
#    plt.xlabel('error in radius of maximum winds (data-relation) [m]')
#    plt.ylabel('probability')
#    plt.grid(True)
#    plt.title('Probability distribution function (fited & data) of the error in RMW',fontweight='bold')
#    plt.show()
#    fnamefig = 'pdfplot1.png'
#    plt.savefig(fnamefig, dpi=300)
#    plt.close('all')
#    
#    # Distribution all
#    idregion = (vmax>0) & (vmax< 10) & (x<-225) & (y < 0)
#    rmaxcheck = rmax[idregion]
#    vmaxcheck= vmax[idregion]
#    kde = gaussian_kde(rmaxcheck)
#    dist_space = np.linspace(0, 150, 500)
#    plt.plot(dist_space, kde(dist_space), label='data')
#    dist = getattr(stats, 'gamma')
#    param = dist.fit(rmaxcheck)
#    pdf_fitted = dist.pdf(dist_space, *param[:-2], loc=param[-2], scale=param[-1])
#    plt.plot(dist_space,pdf_fitted, label=dist_name)
#    
#    # Bootrapsing
#    def bootstrap_resample(X, n=None):
#        """ Bootstrap resample an array_like
#        Parameters
#        ----------
#        X : array_like
#          data to resample
#        n : int, optional
#          length of resampled array, equal to len(X) if n==None
#        Results
#        -------
#        returns X_resamples
#        """
#        if n == None:
#            n = len(X)
#            
#        resample_i = np.floor(np.random.rand(n)*len(X)).astype(int)
#        X_resample = X[resample_i]
#        return X_resample
#    idregion = (vmax>40) & (vmax< 50) & (x<-225) & (y < 0)
#    rmaxcheck = rmax[idregion]
#    rmaxcheck_bootstrap = bootstrap_resample(rmaxcheck, n=10000)
#    rmaxcheck_bootstrap_sort = np.sort(rmaxcheck_bootstrap)
#    rmaxcheck_bootstrap_sort[5000]
#    rmaxcheck_bootstrap_sort[9500]
#    np.mean(rmaxcheck_bootstrap_sort)
#    
#    # Probability density function
#    [tst1, tst2] = stats.ttest_1samp(error, np.mean(error))
#    [tst1, tst2] = stats.kstest(error, 'norm')
#    plt.close(1010)
#    plt.figure(1010,figsize=(15, 7.5))
#    kde = gaussian_kde(error)
#    dist_space = np.linspace( min(error), max(error), 100 )
#    plt.plot(dist_space, kde(dist_space), label='data')
#    for dist_name in dist_names:
#        dist = getattr(stats, dist_name)
#        param = dist.fit(error)
#        pdf_fitted = dist.pdf(xvalues, *param[:-2], loc=param[-2], scale=param[-1])
#        plt.plot(xvalues,pdf_fitted, label=dist_name)
#        plt.xlim(-100,+100)
#    plt.legend(loc='upper right')
#    plt.show()
#    fnamefig = 'pdfplot2.png'
#    plt.savefig(fnamefig, dpi=300)
#    plt.close('all')
#    
#    # Cdfplot 
#    plt.close(1010)
#    plt.figure(1010,figsize=(15, 7.5))
#    dist_space = np.linspace( min(error), max(error), 100 )
#    plt.plot( dist_space, kde(dist_space), label='data')
#    for dist_name in dist_names:
#        dist = getattr(stats, dist_name)
#        param = dist.fit(error)
#        cdf_fitted = dist.cdf(xvalues, *param[:-2], loc=param[-2], scale=param[-1])
#        plt.plot(xvalues,cdf_fitted, label=dist_name)
#        plt.xlim(-100,+100)
#    plt.legend(loc='upper right')
#    plt.show()
#    fnamefig = 'cdfplot.png'
#    plt.savefig(fnamefig, dpi=300)
#    plt.close('all')
#    
#    z95 = 1.96; z50 = 0.674;
#    H,X1 = np.histogram(error, bins = 100, normed = True )
#    dx = X1[1] - X1[0]
#    N = len(error)
#    F1 = np.cumsum(H)*dx
#    X2 = np.sort(error)
#    F2 = np.array(range(N))/float(N)
#    plt.plot(X2, F2, label='data')
#    plt.legend(loc='upper right')
#    plt.ylim([0,1])
#    plt.plot([-z50*rmsesave[-1], -z50*rmsesave[-1]], [0, 1], '--k')
#    plt.plot([z50*rmsesave[-1], z50*rmsesave[-1]], [0, 1], '--k')
#    plt.plot([-z95*rmsesave[-1], -z95*rmsesave[-1]], [0, 1], '--k')
#    plt.plot([z95*rmsesave[-1], z95*rmsesave[-1]], [0, 1], '--k')
#    plt.plot([-100, 100],[0.025, 0.025], '--k')
#    plt.plot([-100, 100],[0.5, 0.5], '--k')
#    plt.plot([-100, 100],[0.975, 0.975], '--k')
#    plt.plot([-100, 100],[0.5, 0.5], '-k')
#    plt.plot([0.0, 0.0],[-100, 100], '-k')
#    
#    ## QQplot
#    plt.close(1010)
#    plt.figure(1010,figsize=(15, 7.5))
#    import pylab 
#    stats.probplot(error, dist="norm", plot=pylab)
#    pylab.show()
#    pylab.show()
#    plt.grid(True)
#    fnamefig = 'qqplot.png'
#    plt.savefig(fnamefig, dpi=300)
#    plt.close('all')

#==============================================================================
# Part 5 confidence and prediction band - based on total
#==============================================================================
# The 95% confidence bands enclose the area that you can be 95% sure contains the true curve. It gives you a visual sense of how well your data define the best-fit curve. SOURCE: http://www.graphpad.com/guides/prism/7/curve-fitting/index.htm?reg_graphing_confidence_and_predic.htm

# Method 1: http://stackoverflow.com/questions/24633664/confidence-interval-for-exponential-curve-fit
# Method 2: http://kitchingroup.cheme.cmu.edu/blog/2013/02/12/Nonlinear-curve-fitting-with-parameter-confidence-intervals/
# Method 3: http://codereview.stackexchange.com/questions/84414/obtaining-prediction-bands-for-regression-model
# https://www.coursera.org/learn/wharton-quantitative-modeling/lecture/Nndhc/4-4-r-squared-and-root-mean-squared-error-rmse
# Limits
idregion    = (rmax>0) & (rmax< 100)
X           = np.array([vmax[idregion],y[idregion]])
rmaxcheck   = rmax[idregion]
vmaxcheck   = vmax[idregion]

# Prediction bands
pp      = (1. + 0.95) / 2.
nstd    = stats.norm.ppf(pp)
tval    = nstd
uppervalue  = func(X, *popt)+tval*rmsesave[-1]
lowervalue  = func(X, *popt)-tval*rmsesave[-1]
normalvalue = func(X, *popt)

#==============================================================================
# Part 6 show solution
#==============================================================================
idreduced = vmaxcheck > 0
[rmse, bias, sci, mae,r2] = modelskill(rmaxcheck[idreduced], normalvalue[idreduced])
bias = 0
print rmse
plt.figure(1006,figsize=(8,6))
sc=plt.scatter(rmaxcheck, normalvalue, s=50, c=vmaxcheck, alpha=0.5,lw = 0, cmap=plt.cm.YlOrRd)
plt.xlabel('$R_{max}$: based on TC database [km]')
plt.ylabel('$R_{max}$: proposed [km]')
plt.legend(loc = 0)
text1a =("RMSE: " ); text1b = (str(round(rmse,1)) +" km")
text2a =("bias: " ); text2b = (str(round(bias,1)) +" km")
text3a =("SCI: " ); text3b = (str(round(sci,3)) +" [-]")
plt.text(70,11,text1a); plt.text(82,11,text1b)
plt.text(70,8,text2a); plt.text(82,8,text2b)
plt.text(70,5,text3a); plt.text(82,5,text3b)
cb = plt.colorbar()
cb.set_label('$v_{max} [m/s]$')
plt.gca().set_aspect('equal', adjustable='box')
plt.clim(0, 80)
plt.ylim(0,100);
plt.xlim(0,100); 
plt.grid(True)
plt.plot([0,100],[0,100],'--k')
plt.title('Radius of maximum winds: observed in data \n versus predicted by proposed relationship',fontweight='bold')
fnamefig = 'rmax_nederhoff_vmax.png'
plt.savefig(fnamefig, dpi=300)
plt.close(1006)

# Difference model and observations
plt.figure(1010,figsize=(8,6))
Z = normalvalue-rmaxcheck
N = len(Z)
H,X1 = np.histogram( Z, bins = 100, normed = True )
dx = X1[1] - X1[0]
F1 = np.cumsum(H)*dx
X2 = np.sort(Z)
F2 = np.array(range(N))/float(N)
plt.plot(X1[1:], F1)
plt.plot(X2, F2)
plt.xlim(-100,100)
plt.ylim(0,1)
plt.grid(True)
plt.xlabel('difference $R_{max}$: TC minus relationship')
plt.ylabel('cumulative probability [-]')
fnamefig = 'cdfplot.png'
plt.savefig(fnamefig, dpi=300)
plt.close(1010)

# Show fit
# Boundary
vmaxs   = np.linspace(0, 70,1000)
ys      = np.linspace(10, 10,1000)
X = np.array([vmaxs, ys])
uppervalue = func(X, *popt)+tval*rmsesave[-1]
lowervalue = func(X, *popt)-tval*rmsesave[-1]
normalvalue= func(X, *popt)
formula1=("$R_{max, storms}$ = a")
formula2=("$R_{max, tropical cyclones}$ = b+c*$v_{max}$")
fit1 = "a={:.1f}, b={:.1f}, c={:.3f}".format(asave[-1], bsave[-1], csave[-1])
plt.close(1007)
idplot1 = (vmax>40) & (vmax < 60) & (rmax > 80) & (rmax < 100)
fig = plt.figure(1007,figsize=(8,6))
h1 = plt.plot(vmaxs,normalvalue,label='proposed relationship')
h2 = plt.fill_between(vmaxs, uppervalue, lowervalue, facecolor='red', alpha=0.5, label='95% prediction interval')
h3 = plt.plot(vmax[~idplot1],rmax[~idplot1],'k.',alpha=0.025, label='TC datapoints')
h4 = plt.plot((0,150), (5.0,5.0),'--k',label='minimum $R_{max} = 5.0 km$', linewidth=4)
plt.xlim(10,70)
plt.ylim(0,150); 
plt.xlabel('maximum sustainted winds [m/s]')
plt.ylabel('radius of maximum winds [km]')
plt.grid(True)
plt.title('Proposed relationship including 95% prediction interval \n compared to data from TC database', fontweight='bold')
plt.legend(loc=0)
plt.text(41,95,formula1)
plt.text(41,89,formula2)
plt.text(41,83,fit1)
fnamefig = 'relationship1.png'
plt.savefig(fnamefig, dpi=300)
plt.close(1007)
rmsesave[-1]