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

siteInfo = pd.read_csv(r'D:\Data\FLUXNET2015\site information1.csv')
igbp_org = siteInfo['IGBP'].values
igbp_new = []
for x in igbp_org:
    if x in (['WET','SNO']):
        x=np.nan
    if x == 'WSA': x = 'SAV'
    if x == 'OSH':  x = 'CSH'
    if x == 'DNF':  x = 'ENF'
    if x == 'DBF':  x = 'EBF'
    igbp_new.append(x)
site_igbp = pd.DataFrame({'site':siteInfo['SITE_ID'].values, 'igbp':igbp_new})
resultH_site = pd.merge(left=resultH_site, right=site_igbp, how='left', on='site', validate='one_to_one').dropna()
resultH_site.to_csv('D:\Data\Global Thermoregulation\For RSE\Data\dT_EC.csv')

order = [2,3,5,6,1,4,0]
dT_EC_PFT_mean = resultH_site.groupby(['igbp'],as_index=False).agg('mean')
dT_EC_PFT_std = resultH_site.groupby(['igbp'],as_index=False).agg('std')

matdata = io.loadmat('D:\Data\Global Thermoregulation\For New Phytologist\Data\dT_satellite_PFT.mat')
dT_satellite_PFT = matdata['dT_satellite_PFT']
dT_satellite_PFT_mean = dT_satellite_PFT[:,0]
dT_satellite_PFT_std = dT_satellite_PFT[:,1]

## slope
slope_site = pd.read_csv('D:\Data\Global Thermoregulation\For New Phytologist\Data\slope_site.csv')
slope_site = pd.merge(left=slope_site, right=site_igbp, how='left', on='site', validate='one_to_one').dropna()
slope_EC_PFT_mean = slope_site.groupby(['igbp'],as_index=False).agg('mean')
slope_EC_PFT_std = slope_site.groupby(['igbp'],as_index=False).agg('std')

matdata = io.loadmat('D:\Data\Global Thermoregulation\For New Phytologist\Data\slope_satellite_PFT.mat')
slope_satellite_PFT = matdata['slope_satellite_PFT']
slope_satellite_PFT_mean = slope_satellite_PFT[:,0]
slope_satellite_PFT_std = slope_satellite_PFT[:,1]

fig, axs = plt.subplots(4, 1, sharex=True, figsize=(7, 6.5))
axs[0].bar(np.arange(start=0, stop=10, step=1.5), dT_EC_PFT_mean['dT'].values[order], width=0.6, yerr=dT_EC_PFT_std['dT'].values[order], fc='#5975A4')
axs[0].set_ylim([0,5])
axs[1].bar(np.arange(start=0, stop=10, step=1.5), dT_satellite_PFT_mean, width=0.6, yerr=dT_satellite_PFT_std, fc='#CC8963')
axs[1].set_ylim([0,5])
axs[2].bar(np.arange(start=0, stop=10, step=1.5), slope_EC_PFT_mean['slope'].values[order], width=0.6, yerr=slope_EC_PFT_std['slope'].values[order], fc='#5975A4',hatch='//')
axs[2].set_ylim([0,1.5])
axs[3].bar(np.arange(start=0, stop=10, step=1.5), slope_satellite_PFT_mean, width=0.6, yerr=slope_satellite_PFT_std, fc='#CC8963',hatch='//')
axs[3].set_ylim([0,1.5])

plt.setp(axs, xticks=np.arange(start=0, stop=10, step=1.5), xticklabels=['BF', 'NF', 'MF', 'SAV', 'SHR', 'GRA', 'CRO'])
fig.tight_layout()
figToPath = r'D:\Data\Global Thermoregulation\For New Phytologist\Figures\PFT_weekly_slope_dT'
fig.savefig(figToPath, dpi=600)
plt.close(fig)

# covariance efficiency
cv_dT_EC = np.std(dT_EC_PFT_mean['dT'])/np.mean(dT_EC_PFT_mean['dT'])
cv_dT_satellite= np.std(dT_satellite_PFT_mean)/np.mean(dT_satellite_PFT_mean)

cv_slope_EC = np.std(slope_EC_PFT_mean['slope'])/np.mean(slope_EC_PFT_mean['slope'])
cv_slope_satellite= np.std(slope_satellite_PFT_mean)/np.mean(slope_satellite_PFT_mean) 

