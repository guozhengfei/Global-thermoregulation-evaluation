import pandas as pd
import random
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import math
from scipy import io

def smooth(x, window_len=11, window='hanning'):
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len<3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError( "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
    s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    y=np.convolve(w/w.sum(),s,mode='same')
    return y[window_len:-window_len+1]

resultH = pd.read_csv(r'D:\Data\Project-3 homeothermy hypothesis test\csv files\EC_data_halfhour2.csv')

LAI_org = resultH['LAI'].values
LAI_org[LAI_org<2]=np.nan
resultH['LAI']=LAI_org

hour_org = resultH['Hour'].values.astype(float)
hour_org[hour_org<1000]=np.nan
hour_org[hour_org>1400]=np.nan
resultH['Hour'] = hour_org
resultH['date'] = pd.to_datetime(resultH['date'])
resultH = resultH.dropna()
resultD_fromH = resultH.groupby([resultH['site'],resultH['date']],as_index=False).agg('mean')

resultH['dT'] = resultH['Tc']-resultH['Tair']
resultH_site = resultH.groupby(resultH['site'],as_index=False).agg('mean')
dT_site = resultH_site['dT'].values
# dT_site[dT_site>5] = 5
fig, ax = plt.subplots(1, 1, figsize=(4, 4))
n, bins, patches = plt.hist(x=resultH_site['dT'].values, bins=np.linspace(0, 5.5, 20), color=[0.09, 0.79, 0.78],  rwidth=1, density=True,edgecolor = 'k')
figToPath = r'D:\Data\Global Thermoregulation\For RSE\Figures\figure_05_dT_EC'
fig.savefig(figToPath, dpi=600)
