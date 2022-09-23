import pandas as pd
import random
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import scipy.io as io
import scipy.stats as st

matFile = r'Y:\project 2\Matlab code\Data\Ta_Noon_denseVegetation_02dgree'
data = io.loadmat(matFile)
Ta = data['TaNoon_reshape']

matFile = r'Y:\project 2\Matlab code\Data\Tc_Noon_denseVegetation_02dgree'
data = io.loadmat(matFile)
Tc = data['TcNoon_reshape']

x = Ta[~np.isnan(Ta+Tc)];
y = Tc[~np.isnan(Ta+Tc)];
fig, axs = plt.subplots(1,figsize=(4,4))
axs.hist2d(x, y,range=[[0, 40],[0,40]], bins=(250, 250), cmap='RdYlBu_r', cmin=50)
axs.plot([0, 40],[0,40],'black')
coef = st.linregress(x,y)
min_v2 = np.percentile(Ta[~np.isnan(Ta)],1)
max_v2 = np.percentile(Ta[~np.isnan(Ta)],99)
axs.plot([min_v2, max_v2], [min_v2*(coef.slope+0.04)+coef.intercept-1, max_v2*(coef.slope+0.04)+coef.intercept-1],'r')
fig.tight_layout()
fig.subplots_adjust(wspace = 0.1)
figToPath = r'D:\Data\Global Thermoregulation\For New Phytologist\Figures\hist_slope_satellite2'
fig.savefig(figToPath, dpi=600)
