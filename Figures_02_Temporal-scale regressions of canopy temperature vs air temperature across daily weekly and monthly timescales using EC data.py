# This code for the daily weekly and monthly Tc vs Ta temporal regression analysis
import pandas as pd
import random
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import math

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

#resultH = pd.read_csv(r'D:\Data\Global Thermoregulation\BESS\climate_LAI_PFT_hourly.csv');
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

site_information = pd.read_csv(r'D:\Data\FLUXNET2015\site information1.csv')
site_igbp = site_information[['SITE_ID', 'IGBP']]; site_igbp.columns = ['site', 'igbp']
igbp_org = site_igbp['igbp'].values
igbp_new = []
for x in igbp_org:
    if x in (['WET','SNO']):
        x=np.nan
    igbp_new.append(x)

resultD_fromH = pd.merge(left=resultD_fromH, right=site_igbp, on='site')
resultD_fromH = resultD_fromH.dropna()

siteNames = resultD_fromH['site'].unique()
Coef_D=[]; Coef_W=[]; Coef_M=[]
Trange_D=[]; Trange_W = []; Trange_M=[]
for x in siteNames:
    dataD = resultD_fromH.loc[resultD_fromH['site']==x,:]
    tair = dataD['Tair'].values
    tc = dataD['Tc'].values
    if tair.size == 0 or tc.size == 0:
        daily_coe = [np.nan]*5
    else:
        fig, axs = plt.subplots(1,3,figsize=(9,3),sharey=True,sharex=True) # plot  density point seasonal
        daily_coe = np.squeeze(st.linregress(tair, tc))
        axs[0].plot(tair, tc, 'o', fillstyle='none');
        axs[0].set_title('Daily')
        min_v1 = np.percentile(np.append(tair, tc),2)
        max_v1 = np.percentile(np.append(tair, tc),98)
        axs[0].plot([min_v1, max_v1], [min_v1 * daily_coe[0] + daily_coe[1], max_v1 * daily_coe[0] + daily_coe[1]])
        axs[0].plot([min_v1, max_v1], [min_v1, max_v1], 'r-')
    Coef_D.append(daily_coe)
    Trange_D.append([min_v1,max_v1])

    dataD.set_index('date', inplace=True)
    dataW = dataD.resample('W-MON').mean().dropna()

    tair = dataW['Tair'].values
    tc = dataW['Tc'].values
    if tair.size <5 or tc.size <5:
        weekly_coe = [np.nan]*5
    else:
        weekly_coe = np.squeeze(st.linregress(tair, tc))  # correlation between  EC and Satellite Tair
        axs[1].plot(tair, tc, 'o', fillstyle='none');
        axs[1].set_title('Weekly')
        min_v1 = min(np.append(tair, tc)) - 1;
        max_v1 = max(np.append(tair, tc)) + 1;
        axs[1].plot([min_v1, max_v1], [min_v1 * weekly_coe[0] + weekly_coe[1], max_v1 * weekly_coe[0] + weekly_coe[1]])
        axs[1].plot([min_v1, max_v1], [min_v1, max_v1], 'r-')
    Coef_W.append(weekly_coe)
    Trange_W.append([min_v1, max_v1])

    dataM = dataD.resample('MS').mean().dropna()
    tair = dataM['Tair'].values
    tc = dataM['Tc'].values
    if tair.size <5 or tc.size <5:
        monthly_coe =  [np.nan]*5
    else:
        monthly_coe = np.squeeze(st.linregress(tair, tc))  # correlation between  EC and Satellite Tair
        axs[2].plot(tair, tc, 'o', fillstyle='none');
        axs[2].set_title('Weekly')
        min_v1 = min(np.append(tair, tc)) - 1;
        max_v1 = max(np.append(tair, tc)) + 1;
        axs[2].plot([min_v1, max_v1], [min_v1 * monthly_coe[0] + monthly_coe[1], max_v1 * monthly_coe[0] + monthly_coe[1]])
        axs[2].plot([min_v1, max_v1], [min_v1, max_v1], 'r-')

    Coef_M.append(monthly_coe)
    Trange_M.append([min_v1, max_v1])

    figToPath =r'D:\Data\Global Thermoregulation\For New Phytologist\Figures\Tc vs Ta plots2\\'+ x + '-' + dataD['igbp'][0]
    fig.tight_layout()
    fig.savefig(figToPath, dpi=600)
    plt.close(fig)
    print(x)

# all EC site together except tropical
data = resultD_fromH[:];
site_org = data['site'].values
site_new = []
for x in site_org:
    if x[0:2] in (['ZZ']):
        x=np.nan
    site_new.append(x)

data['site']=site_new
data = data.dropna()
fig, axs = plt.subplots(2,3,figsize=(9,7.5),sharex=True,gridspec_kw={'height_ratios': [1, 1.5]}) # plot  density point seasonal
x = data['Tair'].values
y = data['Tc'].values
axs[0,0].hist2d(x,y, range=[[0, 40],[0,40]], bins=(200, 200), cmap='RdYlBu_r', cmin=5); #axs[0].colorbar()
axs[0,0].plot([0, 40],[0,40],'black')
coe = st.linregress(x, y)  # correlation between  EC and Satellite Tair
r2 = coe.rvalue**2
slope = coe.slope
intercept = coe.intercept
min_v2 = np.percentile(np.append(x, y),2); max_v2 = np.percentile(np.append(x, y),98);
axs[0,0].plot([min_v2, max_v2], [min_v2*slope+intercept, max_v2*slope+intercept],'r')

data.set_index('date',inplace=True)
data['week']= data.index.isocalendar()['week']
dataW = data.groupby([data['site'], data.index.year, data['week']],as_index=False).agg('mean')
x = dataW['Tair'].values
y = dataW['Tc'].values
axs[0,1].hist2d(x,y, range=[[0, 40],[0,40]], bins=(200, 200), cmap='RdYlBu_r', cmin=2); #axs[0].colorbar()
axs[0,1].plot([0, 40],[0,40],'black')
coe = st.linregress(x, y)  # correlation between  EC and Satellite Tair
r2 = coe.rvalue**2
slope = coe.slope
intercept = coe.intercept
min_v2 = np.percentile(np.append(x, y),1); max_v2 = np.percentile(np.append(x, y),99);
axs[0,1].plot([min_v2, max_v2], [min_v2*slope+intercept, max_v2*slope+intercept],'r')

dataM = data.groupby([data['site'], data.index.year, data.index.month],as_index=False).agg('mean')
x = dataM['Tair'].values
y = dataM['Tc'].values
axs[0,2].hist2d(x,y, range=[[0, 40],[0,40]], bins=(200, 200), cmap='RdYlBu_r', cmin=2); #axs[0].colorbar()
axs[0,2].plot([0, 40],[0,40],'black')
coe = st.linregress(x, y)  # correlation between  EC and Satellite Tair
r2 = coe.rvalue**2
slope = coe.slope
intercept = coe.intercept
min_v2 = np.percentile(np.append(x, y),1); max_v2 = np.percentile(np.append(x, y),99);
axs[0,2].plot([min_v2, max_v2], [min_v2*slope+intercept, max_v2*slope+intercept],'r')

# figToPath = r'D:\Data\Global Thermoregulation\For RSE\Figures\Tc_Ta_densitypoint'
# fig.tight_layout()
# fig.savefig(figToPath, dpi=600)
# plt.close(fig)

# plot each site
randomlist = random.sample(range(0, 159), 66)
Coef_D = np.array(Coef_D)
Coef_D[:,0][Coef_D[:,0]>1.3]=np.nan
Trange_D = np.array(Trange_D)
# fig, axs = plt.subplots(1,3,figsize=(9,4.5)) # plot  density point seasonal
for i in randomlist:
    if ~np.isnan(Coef_D[i,0]) :
        min_v = Trange_D[i, 0]*0.8
        if min_v<18:
            max_v = Trange_D[i, 1]*0.8
            axs[1,0].plot([min_v, max_v], [min_v*Coef_D[i,0]+Coef_D[i,1], max_v*Coef_D[i,0]+Coef_D[i,1]],'gray',lw=0.8)
axs[1,0].set_ylim([-20, 40])
axs[1,0].plot([-0, 40], [0, 40],  'k', linewidth=1.5)

Coef_W = np.array(Coef_W)
Coef_W[:,0][Coef_W[:,0]>1.3]=np.nan
Coef_W[Coef_W[:,2]<0.85,:]=np.nan
Trange_W = np.array(Trange_W)
for i in randomlist:
    if ~np.isnan(Coef_W[i,0]) :
        min_v = Trange_W[i, 0]*0.8
        if min_v<18:
            max_v = Trange_W[i, 1]*0.8
            axs[1,1].plot([min_v, max_v], [min_v*Coef_W[i,0]+Coef_W[i,1], max_v*Coef_W[i,0]+Coef_W[i,1]],'gray',lw=0.8)
axs[1,1].set_ylim([-20, 40])
axs[1,1].plot([-0, 40], [0, 40],  'k', linewidth=1.5)

Coef_M = np.array(Coef_M)
Coef_M[:,0][Coef_M[:,0]>1.3]=np.nan
Coef_M[Coef_M[:,2]<0.85,:]=np.nan
Trange_M = np.array(Trange_M)
for i in randomlist:
    if ~np.isnan(Coef_M[i,0]) :
        min_v = Trange_M[i, 0]*0.8
        if min_v<18:
            max_v = Trange_M[i, 1]*0.8
            axs[1,2].plot([min_v, max_v], [min_v*Coef_M[i,0]+Coef_M[i,1], max_v*Coef_M[i,0]+Coef_M[i,1]],'gray',lw=0.8)
axs[1,2].set_ylim([-20, 40])
axs[1,2].plot([-0, 40], [0, 40],  'k', linewidth=1.5)

fig.tight_layout()
fig.subplots_adjust(hspace = 0.01)
figToPath = r'D:\Data\Global Thermoregulation\For New Phytologist\Figures\Tc_Ta_each_site'
fig.savefig(figToPath, dpi=600)

Coef_D[Coef_D[:,2]<0.85,:]=np.nan
Coef_D_nonan = Coef_D[~np.isnan(Coef_D[:,0]),:]
fig, axs = plt.subplots(1, 2, figsize=(2.8, 1.5))
density = st.gaussian_kde(Coef_D_nonan[:,0])  # slope
n, x, _ = plt.hist(Coef_D_nonan[:,0], bins=np.linspace(0.5, 1.5, 20),histtype=u'step', density=True)
axs[0].plot(x, density(x))
density = st.gaussian_kde(Coef_D_nonan[:,2]) # R2
n, x, _ = plt.hist(Coef_D_nonan[:,2], bins=np.linspace(0.75, 1, 20), histtype=u'step', density=True)
axs[1].clear(); axs[1].plot(x, density(x))
axs[0].yaxis.set_visible(False);axs[1].yaxis.set_visible(False)
fig.tight_layout()
fig.subplots_adjust(wspace = 0.1)
figToPath = r'D:\Data\Global Thermoregulation\For New Phytologist\Figures\hist_D'
fig.savefig(figToPath, dpi=600)

Coef_W[Coef_W[:,2]<0.85,:]=np.nan
Coef_W_nonan = Coef_W[~np.isnan(Coef_W[:,0]),:]
fig, axs = plt.subplots(1, 2, figsize=(2.8, 1.5))
density = st.gaussian_kde(Coef_W_nonan[:,0])  # slope
n, x, _ = plt.hist(Coef_W_nonan[:,0], bins=np.linspace(0.5, 1.5, 20),histtype=u'step', density=True)
axs[0].plot(x, density(x))
density = st.gaussian_kde(Coef_W_nonan[:,2]) # R2
n, x, _ = plt.hist(Coef_W_nonan[:,2], bins=np.linspace(0.75, 1, 20), histtype=u'step', density=True)
axs[1].clear(); axs[1].plot(x, density(x))
axs[0].yaxis.set_visible(False);axs[1].yaxis.set_visible(False)
fig.tight_layout()
fig.subplots_adjust(wspace = 0.1)
figToPath = r'D:\Data\Global Thermoregulation\For New Phytologist\Figures\hist_W'
fig.savefig(figToPath, dpi=600)

Coef_M[Coef_M[:,2]<0.85,:]=np.nan
Coef_M_nonan = Coef_M[~np.isnan(Coef_M[:,0]),:]
fig, axs = plt.subplots(1, 2, figsize=(2.8, 1.5))
density = st.gaussian_kde(Coef_M_nonan[:,0])  # slope
n, x, _ = plt.hist(Coef_M_nonan[:,0], bins=np.linspace(0.5, 1.5, 20),histtype=u'step', density=True)
axs[0].plot(x, density(x))
density = st.gaussian_kde(Coef_M_nonan[:,2]) # R2
n, x, _ = plt.hist(Coef_M_nonan[:,2], bins=np.linspace(0.75, 1, 20), histtype=u'step', density=True)
axs[1].clear(); axs[1].plot(x, density(x))
axs[0].yaxis.set_visible(False);axs[1].yaxis.set_visible(False)
fig.tight_layout()
fig.subplots_adjust(wspace = 0.1)
figToPath = r'D:\Data\Global Thermoregulation\For New Phytologist\Figures\hist_M'
fig.savefig(figToPath, dpi=600)


##
coef_D = np.array(Coef_D)[:,0]; coef_W = np.array(Coef_W)[:,0]
st.linregress(coef_D[~np.isnan(coef_D+coef_W)],coef_W[~np.isnan(coef_D+coef_W)])
plt.figure();plt.plot(coef_D[~np.isnan(coef_D+coef_W)],coef_W[~np.isnan(coef_D+coef_W)],'o')

slope_EC_PFT = pd.DataFrame({'site':siteNames, 'slope':np.array(Coef_D)[:,0]})
slope_EC_PFT.to_csv('D:\Data\Global Thermoregulation\For New Phytologist\Data\slope_site')
