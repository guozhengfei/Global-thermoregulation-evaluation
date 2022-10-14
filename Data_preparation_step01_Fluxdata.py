import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import scipy.stats as st
from scipy.stats import gaussian_kde
import seaborn as sns
import time
from datetime import datetime
import scipy.io as io
import random
from scipy.signal import savgol_filter

def smooth(x,window_len=11,window='hanning'):
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
# Canopy temperature inversion
def T_LW(LWout, LWin, emissivity):
    sig = 5.6704 * 10 ** -8
    Tground = ((LWout - (1-emissivity)*LWin) / (emissivity * sig)) ** (1 / 4) - 273.15
    Tground[Tground<-50] = np.nan
    return Tground
# areodynamic temperature inversion
def T_Aero(u, u_star, H, Tair):
    cp = 1.005 * 10**3   # specific heat of air at constant pressure (J kg-1 K-1)
    pa = 0.023*Tair + 1.1139    # air density (kg m-3), Tair (degree C)
    Ra = u/u_star**2 + 6.2*u_star**(-2/3)
    Ra_nonan = Ra[~np.isnan(Ra)]
    min_v = np.percentile(Ra_nonan, 3)
    max_v = min(np.percentile(Ra_nonan, 97), 100)
    Ra_nonan[Ra_nonan>max_v] = max_v
    Ra_nonan[Ra_nonan < min_v] = min_v
    Ra[~np.isnan(Ra)] = Ra_nonan
    return Ra*H/(cp*pa) + Tair
# Outlier remove
def IQR_filter(array):
    p25 = np.percentile(array[~np.isnan(array)], 25)
    p75 = np.percentile(array[~np.isnan(array)],75)
    IQR = p75-p25
    array[array < p25 - 1.5*IQR] = np.nan
    array[array > p75 + 1.5*IQR] = np.nan
    arraynew = array
    return arraynew

def Arr_reverse(x):
    return x.values[::-1]
    # x.loc[:,] = y
    # return

# read and deal emissivity  and LAI data
MYD21_df = pd.read_csv(r'D:\Data\Global Thermoregulation\MODIS LST validation with FLUXNET\Emi_LST_MYD21.csv')
colNameList = list(MYD21_df.columns)
QC = MYD21_df['MYD21A1D_006_QC_Emissivity_accuracy'].values
QC_temp = []
for x in QC:
    QC_temp.append(float(x.split('b')[-1]))

QC_lab = np.array(QC_temp)
QC_lab[QC_lab<10] = np.nan
QC_lab[~np.isnan(QC_lab)] = 1

MODISdata_time = MYD21_df["Date"].values
Emis_29 = MYD21_df['MYD21A1D_006_Emis_29'].values*QC_lab  # *0.002+0.49
Emis_31 = MYD21_df['MYD21A1D_006_Emis_31'].values*QC_lab  # *0.002+0.49
Emis_32 = MYD21_df['MYD21A1D_006_Emis_32'].values*QC_lab  # *0.002+0.49
Emis = (Emis_29 + Emis_31 + Emis_32)/3

MODIS_time = MYD21_df["Date"].values
times = MODIS_time.tolist()
time_list = []
for x in times:
    time_list.append(datetime.strptime(x, "%Y-%m-%d"))
times_extend = pd.date_range('2001-01-01', '2015-12-31',  freq='D')  # for emissivity data
times_month = pd.date_range('2001-01-01', '2015-12-31',  freq='MS')  # for LAI data
Emis_data_0 = pd.DataFrame({'date': times_extend})
Emis_data_0 = Emis_data_0.set_index('date')
siteID_LST = MYD21_df['ID'].values
siteID_unique = MYD21_df['ID'].unique()
IGBPList = MYD21_df['Category'].values

# calculate the emissivity seasonal pattern of each site
Emis_data = pd.DataFrame({'emis': Emis, 'date': time_list, 'site': siteID_LST, 'igbp':IGBPList})
Emis_data = Emis_data.set_index('date')
Emis_data_extend = []

# read LAI
LAI_df = pd.read_csv(r'D:\Data\Global Thermoregulation\MODIS LST validation with FLUXNET\LAI_gapfill_monthly_csv.csv')
colNameList = list(MYD21_df.columns)

for x in siteID_unique:
    lai_site_i = LAI_df[LAI_df['SITE_ID'].isin([x])]
    lai_site_i_values = lai_site_i.iloc[:, :-3].values.T[:,0]
    LAI = pd.DataFrame({'date': times_month,  'LAI':lai_site_i_values}).set_index('date')
    LAI = LAI.reindex(times_extend, axis=0).interpolate(axis=0)
    emis_site_i = Emis_data[Emis_data['site'].isin([x])]
    igbp = emis_site_i['igbp'].unique()[0]
    # emis_site_i['emis'].interpolate(method='linear', inplace=True)
    emis_site_i_values = emis_site_i['emis'].interpolate(method='linear').values
    emis_site_i_values = IQR_filter(emis_site_i_values)
    emis_site_i_values[np.isnan(emis_site_i_values)] = np.nanmean(emis_site_i_values)
    emis_site_i_values = smooth(emis_site_i_values, 8, 'flat')
    # extend from 2000-01-01
    emis_site_i = emis_site_i.drop('emis', 1)
    emis_site_i['emis'] = emis_site_i_values
    emis_mean = emis_site_i.groupby(emis_site_i.index.dayofyear).agg('mean').values.flatten()

    Emis_data_extend_i = pd.merge(left=Emis_data_0, right=emis_site_i, how='left', on='date', validate='one_to_one')
    Emis_data_extend_i['site'] = [x]*Emis_data_extend_i.shape[0]
    Emis_data_extend_i['igbp'] = [igbp]*Emis_data_extend_i.shape[0]
    Emis_data_extend_i.loc[0:emis_mean.shape[0]*2, 'emis'] = np.append(emis_mean,emis_mean)
    nan_df = Emis_data_extend_i[np.isnan(Emis_data_extend_i['emis'])]
    if nan_df.shape[0]>0:
        Emis_data_extend_i.loc[np.isnan(Emis_data_extend_i['emis']),'site'] = [x]*nan_df.shape[0]
        Emis_data_extend_i.loc[np.isnan(Emis_data_extend_i['emis']),'emis'] = [emis_mean.mean()]*nan_df.shape[0]
    Emis_data_extend_i = pd.merge(left=Emis_data_extend_i, right=LAI, how='left', left_index=True, right_index=True, validate='one_to_one')
    Emis_data_extend.append(Emis_data_extend_i)
Emis_LAI_data = pd.concat(Emis_data_extend)
del Emis,Emis_29,Emis_31,Emis_32,MYD21_df

# read FLUXNET data and calculate canopy temperature
siteInfo = pd.read_csv(r'D:\Data\FLUXNET2015\site information1.csv')
FluxFilepath = r'D:\Data\FLUXNET2015\\'
FluxFileFolders = os.listdir(FluxFilepath)
FluxFileFolders = [x for x in FluxFileFolders if x.split('_')[-1].split('-')[-1] == '4']

R2_halfhour = []; R2_daily = []; R2_monthly = [];  R2_yearly = []; # R2
Slope_halfhour = []; Slope_daily = []; Slope_monthly = []; Slope_yearly = []  # Slope
intercept_halfhour = []; intercept_daily = []; intercept_monthly = []; intercept_yearly = []  # Slope
EC_Tc_Ta_halfhour = []; EC_Tc_Ta_daily = []; EC_Tc_Ta_monthly = [];
Trange = []; EC_Tc_Ta_all = [];

for FluxFileFolder in FluxFileFolders:
    FluxFile = FluxFilepath + FluxFileFolder + r'\\' + FluxFileFolder.split('SUBSET')[0] + 'SUBSET_HH' + \
               FluxFileFolder.split('SUBSET')[-1] + '.csv'
    FluxData = pd.read_csv(FluxFile)
    siteName = FluxFileFolder.split('_')[1]

    #  find the corresponding site
    index = Emis_LAI_data['site'] == siteName
    igbp = Emis_LAI_data['igbp'][index][0]
    Emi_df = Emis_LAI_data[index]

    # merge the fluxdata with emissivity data
    if siteName[0:2] == 'ZZ':
        timestr = []
        for i in range(FluxData['Year_LBAMIP'].size):
            timestr.append(datetime.strptime(str(FluxData['Year_LBAMIP'].values[i]) + "-" + str(FluxData['DoY_LBAMIP'].values[i]) ,
                                             "%Y-%j").strftime("%Y%m%d")+str(FluxData['Hour_LBAMIP'].values[i]).zfill(2)+'00')
        FluxData['TIMESTAMP_START'] = timestr
        gppdata = FluxData['GEP_model'].values # umol CO2 m-2 s-1
        gppdata[gppdata<0]=np.nan
        PAR = FluxData['rgs'].values*2  # unit: u mol m-2 s-1
        PAR[PAR<-10]=np.nan
        PAR_org = PAR*1
        PAR[PAR < 20] = np.nan
        RH = FluxData['rh'].values.astype(float)
        RF = FluxData['prec'].values
        LE = FluxData['LE']
    else:
        FluxData['TIMESTAMP_START'] = FluxData['TIMESTAMP_START'].apply(str)
        gppdata = FluxData['GPP_DT_VUT_REF'].values # umol CO2 m-2 s-1
        gppdata[gppdata<0]=np.nan
        PAR = FluxData['SW_IN_F'].values*2  # unit: u mol m-2 s-1
        PAR[PAR<-10]=np.nan
        PAR_org = PAR*1
        PAR[PAR < 20] = np.nan
        if 'RH' in FluxData.columns:
            RH = FluxData['RH'].values
        else:
            Ta = FluxData['TA_F'].values
            Ta[Ta<-100]=np.nan
            VPD =  FluxData['VPD_F'].values
            es = 2.1718e10 * np.exp(-4157./(Ta+273.15-33.91)) #  [Pa]
            RH = (es-VPD*100)/es*100

        RF = FluxData['P_F'].values
        LE = FluxData['LE_CORR']

    Time_flx = FluxData['TIMESTAMP_START'].values
    # apply date to index
    date_flux = []
    Hour_flux = []
    for x in Time_flx:
        date_flux.append(datetime(year=int(x[0:4]), month=int(x[4:6]), day=int(x[6:8])))
        Hour_flux.append(str(x[8:12]))

    FluxData['date'] = date_flux
    FluxData.set_index('date')
    # calculate daily gpp to help derive growth seasaon
    gpp_hour = pd.DataFrame({'gpp': gppdata,'date': date_flux})
    gpp_hour.set_index('date', inplace=True)
    gpp_daily = gpp_hour.resample('D').sum()
    gpp_daily['year'] = gpp_daily.index.year
    y = gpp_daily['gpp'].values
    y1 = np.maximum(y,savgol_filter(y, 31, 2)); gpp_daily['gpp']  = np.maximum(y1,savgol_filter(y1, 31, 2));
    # plt.figure();plt.plot(gpp_daily.index,gpp_daily['gpp'])
    gpp_threshold = 0.25*(gpp_daily.groupby(gpp_daily.index.year).agg('max')-gpp_daily.groupby(gpp_daily.index.year).agg('min'))+gpp_daily.groupby(gpp_daily.index.year).agg('min') # set growth season value as gpp min +15% amplitude
    gpp_threshold = gpp_threshold.rename(columns={'gpp':'gppthreshold'})
    gpp_threshold['year'] = gpp_threshold.index
    gpp_daily['date'] = gpp_daily.index
    gpp_daily = pd.merge(left=gpp_daily, right=gpp_threshold, how='left', on='year', validate='many_to_one')
    gpp_daily['GS_label'] = np.round(smooth((gpp_daily['gpp'] > gpp_daily['gppthreshold']).values*1,10,'flat'))# 1 is in the growth season, 0 is out of growth season
    gpp_daily.set_index('date')
    # plt.figure();plt.plot(gpp_daily.index,gpp_daily['gpp']); plt.plot(gpp_daily.index,gpp_daily['GS_label']*200)
    # gpp_hour = pd.merge(left=gpp_hour, right=gpp_daily, how='left', on='date', validate='many_to_one')
    gpp_daily = gpp_daily.drop(columns=['gpp','year'])
    #  merge Flux data with emissivity for the Tc calculation
    FluxData = pd.merge(left=FluxData, right=Emi_df, how='left', on='date', validate='many_to_one')
    FluxData = pd.merge(left=FluxData, right=gpp_daily, how='left', on='date', validate='many_to_one')

    del gpp_hour, gpp_daily
    Tair = FluxData['TA_F'] # air temperature unit: oC
    Tair[Tair < -1000] = np.nan
    # RH = FluxData['RH'] # relative humidity range :0-100
    RH[RH < -1000] = np.nan
    # RF = FluxData['P_F'] # precipitation, unit: mm/hour or half-hour
    RF[RF < -1000] = np.nan
    LE[LE < -1000] = np.nan
    LAI = FluxData['LAI'].values/10 # leaf area index unit: 0-10
    u = FluxData['WS_F'] # wind speed unit: m/s
    u[u < -1000] = np.nan
    u_star = FluxData['USTAR']
    u_star[u_star < -1000]=np.nan
    H = FluxData['H_CORR']
    H[H < -1000] = np.nan
    if siteName[0:2] != 'ZZ':
        H2 = FluxData['H_F_MDS']
        H2[H2<-1000] = np.nan
        H = pd.DataFrame({'H':H,'H2':H2}).bfill(axis ='columns')['H']
        LE2 = FluxData['LE_F_MDS']
        LE2[LE2<-1000] = np.nan
        LE = pd.DataFrame({'LE':LE,'LE2':LE2}).bfill(axis ='columns')['LE']
    H_QC = FluxData['H_F_MDS_QC'].astype(float)
    H_QC[H_QC == 2] = np.nan # remove the low quality values
    H_QC[H_QC == 3] = np.nan # remove the low quality values
    H_QC[~np.isnan(H_QC)]=1
    # H.values[np.isnan(H_QC)] = np.nan
    H = H*H_QC

    GS_label = FluxData['GS_label'].values
    if igbp == 'EBF': # how about ENF
        GS_label[:] = 1
    else:
        GS_label[GS_label==0] = np.nan
    # plt.figure(); plt.plot(GS_label)
    site =  FluxData['site'].values
    Tareo = T_Aero(u, u_star, H, Tair)
    if np.isnan(Tareo.values).sum() < Tareo.values.size:
        Tareo = IQR_filter(Tareo)
    Tc = Tareo.values
    Tc[np.where((Tc < 0) | (Tc > 45))] = np.nan
    Tc_lab = Tc*1; Tc_lab[~np.isnan(Tc_lab)] = 1;
    Tair[Tair<0] = np.nan

    # average Tcan with equal weight
    # note: force Tair equal to null when Tc is null
    data = pd.DataFrame({'date': date_flux,'Hour': Hour_flux, 'Tair': Tair.values*GS_label*Tc_lab,'Tc': Tc*GS_label,'GS_label': GS_label, 'GPP': gppdata, 'PAR':PAR*GS_label*Tc_lab, 'PAR_org':PAR_org*GS_label*Tc_lab,'site': site, 'RH': RH*GS_label*Tc_lab, 'rainfall': RF*GS_label*Tc_lab,'windspeed':u, 'LAI': LAI*GS_label*Tc_lab,'H':H,'LE':LE,'igbp':FluxData['igbp']})

    data = data.set_index('date')
    EC_Tc_Ta_halfhour.append(data)
    print(FluxFileFolder)

resultH = pd.concat(EC_Tc_Ta_halfhour)

resultH.to_csv(r'D:\Data\Project-3 homeothermy hypothesis test\csv files\EC_data_halfhour2.csv') 
