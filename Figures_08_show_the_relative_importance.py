import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# relative importance for linear model
coefs_EC = pd.read_csv(r'Y:\project 2\RF analysis for global dT variation\dT_drivers_EC_importance.csv')

fig, ax = plt.subplots(figsize=(4,6))
plt.errorbar(coefs_EC['Estimate'].values[0:6], np.arange(start=1, stop=7, step=1),  xerr=coefs_EC['Std'].values[0:6]*1.5, fmt='o', color='#66c2a5',
         elinewidth=2,  fillstyle='none');
plt.errorbar(coefs_EC['Estimate'].values[6:8], np.arange(start=7, stop=9, step=1),  xerr=coefs_EC['Std'].values[6:8]*1.5, fmt='o', color='#fc8d62',
             elinewidth=2,  fillstyle='none');
plt.errorbar(coefs_EC['Estimate'].values[8:], np.arange(start=9, stop=15, step=1),  xerr=coefs_EC['Std'].values[8:]*1.5, fmt='o', color='#8da0cb',
             elinewidth=2,  fillstyle='none');
plt.xlim([-0.3, 0.3])
plt.ylim([0, 15])
plt.plot([0 ,0],[16,0],'k--',linewidth=2)
ax.get_yaxis().set_visible(False)
figToPath = r'D:\Data\Global Thermoregulation\MODIS LST validation with FLUXNET\Figures\dT_drivers_EC'
fig.tight_layout()
fig.savefig(figToPath, dpi=600)
plt.close(fig)


fig, ax = plt.subplots(figsize=(2,6))
abiotic = abs(coefs_EC['Estimate'].values[0:6]).sum()/abs(coefs_EC['Estimate'].values).sum()
biotic = abs(coefs_EC['Estimate'].values[6:8]).sum()/abs(coefs_EC['Estimate'].values).sum()
interct = abs(coefs_EC['Estimate'].values[8:]).sum()/abs(coefs_EC['Estimate'].values).sum()

plt.bar(1,100*(abiotic+biotic+interct),color='#8da0cb',width = 1)
plt.bar(1,100*(abiotic+biotic),color='#fc8d62',width = 1)
plt.bar(1,100*abiotic,color='#66c2a5',width = 1)
ax.get_xaxis().set_visible(False)
plt.xlim([0.5, 1.5])
plt.ylim([0, 100])
figToPath = r'D:\Data\Global Thermoregulation\MODIS LST validation with FLUXNET\Figures\dT_importance_EC'
fig.tight_layout()
fig.savefig(figToPath, dpi=600)
plt.close(fig)

# relative importance for linear model (satellite data)
coefs_EC = pd.read_csv(r'Y:\project 2\RF analysis for global dT variation\dT_drivers_Satellite_importance.csv')

fig, ax = plt.subplots(figsize=(4,6))
plt.errorbar(coefs_EC['Estimate'].values[0:6], np.arange(start=1, stop=7, step=1),  xerr=coefs_EC['Std'].values[0:6]*1.5, fmt='o', color='#66c2a5',
             elinewidth=2,  fillstyle='none');
plt.errorbar(coefs_EC['Estimate'].values[6:8], np.arange(start=7, stop=9, step=1),  xerr=coefs_EC['Std'].values[6:8]*1.5, fmt='o', color='#fc8d62',
             elinewidth=2,  fillstyle='none');
plt.errorbar(coefs_EC['Estimate'].values[8:], np.arange(start=9, stop=15, step=1),  xerr=coefs_EC['Std'].values[8:]*1.5, fmt='o', color='#8da0cb',
             elinewidth=2,  fillstyle='none');
plt.xlim([-0.3, 0.3])
plt.ylim([0, 15])
plt.plot([0 ,0],[16,0],'k--',linewidth=2)
ax.get_yaxis().set_visible(False)
figToPath = r'D:\Data\Global Thermoregulation\MODIS LST validation with FLUXNET\Figures\dT_drivers_Satellite'
fig.tight_layout()
fig.savefig(figToPath, dpi=600)
plt.close(fig)


fig, ax = plt.subplots(figsize=(2,6))
abiotic = abs(coefs_EC['Estimate'].values[0:6]).sum()/abs(coefs_EC['Estimate'].values).sum()
biotic = abs(coefs_EC['Estimate'].values[6:8]).sum()/abs(coefs_EC['Estimate'].values).sum()
interct = abs(coefs_EC['Estimate'].values[8:]).sum()/abs(coefs_EC['Estimate'].values).sum()

plt.bar(1,100*(abiotic+biotic+interct),color='#8da0cb',width = 1)
plt.bar(1,100*(abiotic+biotic),color='#fc8d62',width = 1)
plt.bar(1,100*abiotic,color='#66c2a5',width = 1)
ax.get_xaxis().set_visible(False)
plt.xlim([0.5, 1.5])
plt.ylim([0, 100])
figToPath = r'D:\Data\Global Thermoregulation\MODIS LST validation with FLUXNET\Figures\dT_importance_Satellite'
fig.tight_layout()
fig.savefig(figToPath, dpi=600)
plt.close(fig)
