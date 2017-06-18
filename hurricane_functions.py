# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 15:02:48 2016
@author: nederhof
"""
# Functions needed
import numpy as np
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import r2_score


# Standard calulcations
def modelskill(measured,computed):
    
    # Filter: only for real numbers
    id1=np.isreal(measured)
    id2=np.isreal(computed)
    id= id1&id2
    zmt=measured[id];
    zct=computed[id];

    # Compute
    rmse= np.sqrt(mean_squared_error(measured, computed))
    rmsm=np.sqrt(np.mean(np.power(zmt,2)))
    tmp = max(rmsm,(abs(np.mean(zmt))))
    sci =rmse/tmp
    bias=np.mean(zct-zmt);
    mae =mean_absolute_error(zmt,zct)
    r2=r2_score(zmt,zct)
    return rmse,bias,sci,mae,r2


# Radius of Maximum Winds (RMW)
# Calculate Vickery & Wadhera (2008) values
def vickery_wadhera(dpc,y):
    a=3.015;
    b=6.291*pow(10,-5)
    c=0.0337;
    RMW = np.exp(a-b*np.power(dpc,2)+c*y)
    return RMW
# Vickery, P. J., & Wadhera, D. (2008). Winds of Hurricanes from Flight-Level Pressure and H * Wind Data. Journal of Applied Meteorology and Climatology. http://doi.org/10.1175/2008JAMC1837.1

# Calculate Knaff and Zehr (2007) values
def knaff_zehr(vmax,y):
    a=66.785;
    b=0.09102
    c=1.0619;
    RMW = a-b*vmax+c*(y-25)
    return RMW
# Knaff, J. A., & Zehr, R. M. (2007). Reexamination of tropical cyclone wind-pressure relationships. Bulletin of the American Meteorological Society, 88(3), 71–88. http://doi.org/10.1175/WAF965.1


# Calculate Gross (2004) values
def gross(vmax,y):
    a=35.37;
    b=0.111
    c=0.570;
    RMW_NM = a-b*vmax+c*(y-25)
    RMW=RMW_NM*1.852
    return RMW
# Gross, J. M., Demaria, M., Knaff, J. A., & Sampson, C. R. (2004). A new method for determining tropical cyclone wind forecast probabilities, 570(2), 2–4.