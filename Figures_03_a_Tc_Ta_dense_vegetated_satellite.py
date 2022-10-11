from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as st

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

Ta_df = pd.read_csv(r'D:\Data\Global Thermoregulation\For New Phytologist\Data\Ta_gs_csv3.csv')
Tc_df = pd.read_csv(r'D:\Data\Global Thermoregulation\For New Phytologist\Data\Tc_gs_csv3.csv')

colNameList_Ta = list(Ta_df.columns)
colNameList_Ta = [x for x in colNameList_Ta if x[0] == 'b']
Ta_df2 = Ta_df[colNameList_Ta]
colNameList_Ta2 = []
for name in colNameList_Ta:
    colNameList_Ta2.append(name[1:].zfill(4))
Ta_df2.columns = colNameList_Ta2
Ta_df2 = Ta_df2.reindex(sorted(Ta_df2.columns), axis=1)

colNameList_Tc = list(Tc_df.columns)
colNameList_Tc = [x for x in colNameList_Tc if x.split('_')[-1][0:3] == '1km']
Tc_df2 = Tc_df[colNameList_Tc]
colNameList_Tc2 = []
for name in colNameList_Tc:
    colNameList_Tc2.append(name.split('_')[0].zfill(4))
Tc_df2.columns = colNameList_Tc2
Tc_df2 = Tc_df2.reindex(sorted(Tc_df2.columns), axis=1)

Tc_fill = []; Ta_gs = []; Tc_gs = []; DT = []; COEF = []; COEF_fit = []; p_95=[]; p_05=[]
for i in range(Ta_df2.shape[0]):
    tc = Tc_df2.loc[i,:].values*1
    if tc[~np.isnan(tc)].size > 80:
        ta = Ta_df2.loc[i,:].values*1; ratio_coef = 0.975;ta =ta*ratio_coef; #ratio_coef from the relationship between EC and reanalysis data
        dt = tc-ta; dt[dt>np.percentile(dt[~np.isnan(dt)],94)] = np.nan; # dt[dt<-3]=np.nan;
        dt[dt<np.percentile(dt[~np.isnan(dt)],6)]=np.nan;
        p_95.append(np.percentile(dt[~np.isnan(dt)],94))
        p_05.append(np.percentile(dt[~np.isnan(dt)], 4))
        tc = tc+dt-dt;
        if tc[~np.isnan(tc)].size > 80:
            perc = np.array([49,50,60,65])
            # ta_rsp = ta.reshape(10,46)
            # ta_mean = smooth(np.mean(ta_rsp,axis=0),3,'flat')
            thresholds = np.abs(np.array([np.percentile(ta,49),np.percentile(ta,50),np.percentile(ta,60),np.percentile(ta,65)])-10.5)
            perc_thre = perc[np.where(thresholds == thresholds.min())]
            threshold = np.percentile(ta, perc_thre)
            # gs_label = ta_mean*1
            # gs_label[gs_label < threshold] = np.nan
            # gs_label[~np.isnan(gs_label)] = 1
            # index = 0
            # for i in range(46):
            #     x = gs_label[(i+1)*-1]
            #     index = index+1
            #     if x == 1: break
            # gs_label[-(index-1):-(index-5)] = 1
            # gs_label_mul = np.tile(gs_label, 10)*0
            # tc = tc+gs_label_mul

            coef = st.linregress(ta[~np.isnan(tc)], tc[~np.isnan(tc)])
            COEF_fit.append(coef.slope)
            tc_gap_fill = ta*(coef.slope)+coef.intercept
            tc_gap_fill[~np.isnan(tc)] = 0
            tc[np.isnan(tc)] = 0
            tc_final = tc_gap_fill+tc

            # threshold = np.maximum(np.percentile(ta,40), 9)
            ta_gs = ta*1; ta_gs[ta_gs<threshold] = np.nan
            if ta_gs[~np.isnan(ta_gs)].size > 80:
                Ta_gs.append(ta_gs)
                ## orginal
                tc[tc==0]=np.nan ; tc = tc+ta_gs-ta_gs;
                coef_org = st.linregress(ta_gs[~np.isnan(tc)],tc[~np.isnan(tc)])
                dT_org = np.mean(tc[~np.isnan(tc)]-ta_gs[~np.isnan(tc)])
                ## gap-filled
                tc = tc_final * 1+ ta_gs - ta_gs;
                coef_fill = st.linregress(ta_gs[~np.isnan(tc)], tc[~np.isnan(tc)])
                dT_fill = np.mean(tc[~np.isnan(tc)] - ta_gs[~np.isnan(tc)])
                Tc_fill.append(tc_final)
                Tc_gs.append(tc)

                DT.append([dT_org, dT_fill])
                COEF.append([coef_org.slope, coef_fill.slope])

        if i==1:
            fig, axs = plt.subplots(1, sharex=True, figsize=(8, 3))
            plt.plot(tc_final[0:230], 'o-',mfc='None')
            plt.plot(Tc_df2.loc[i,:].values[0:230], 'o-',mfc='None')
            axs.set_ylim(-22,37)
            fig.tight_layout()
            figToPath = r'D:\Data\Global Thermoregulation\For New Phytologist\Figures\gap_fill_example2'
            fig.savefig(figToPath, dpi=600)
            plt.close(fig)

Tc_fill = np.array(Tc_fill)
Tc_gs = np.array(Tc_gs)
Ta_gs = np.array(Ta_gs)
DT = np.array(DT)
COEF = np.array(COEF)
p_95 = np.array(p_95)
p_05 = np.array(p_05)

# np.sum(p_95<=8)/p_95.size
# np.sum(p_05>=-3)/p_05.size
# plt.hist(p_05,50)

fig, axs = plt.subplots(1, sharex=True, figsize=(4, 4))
x = Ta_gs.flatten();
y = Tc_gs.flatten()
plt.hist2d(x[~np.isnan(x)], y[~np.isnan(y)]+0.3, bins=(300, 300),vmax=300,cmap='RdYlBu_r', cmin=40)
coef = st.linregress(x[~np.isnan(x)], y[~np.isnan(y)])
axs.set_ylim(0,40); axs.set_xlim(0,40)
axs.plot([0,40], [0,40], 'k-')
axs.plot([8,32.5], [8*coef.slope+coef.intercept,32.5*coef.slope+coef.intercept], 'r-')
fig.tight_layout()
