import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#Table 5 of Holczer et al. 2016, taken from Vizier
data = pd.read_csv('KeplerTTV_amplitudes.txt',sep='|', comment='#')

#creating histogram of ttv amplitudes
fig = plt.figure(figsize=(8,4))

#histogram shows counts rather than probability density, with log-spaced bins
n, bins, patches = plt.hist(np.array(data['Amp']), 
    bins=np.logspace(np.log10(np.min(data['Amp'])), np.log10(np.max(data['Amp'])),25),
    histtype='step',color='#9d98ae',label='Holczer et al., 2016')

#vertical line where KOI-134b's ttv amplitude is located
plt.vlines(20*60, 0,15,linestyle='--',color='#02315d')
plt.semilogx()
plt.xticks(ticks=([1,10,100,1000]),labels=([1,10,100,1000]))
plt.xlabel('TTV Amplitude (Min)',fontsize=14)
plt.ylabel('Counts',fontsize=14)
plt.annotate('KOI-134 b', (550,16),fontsize=14, color='#02315d')
plt.tick_params(axis='both', which='major', labelsize=14)
plt.legend(fontsize=12)
fig.savefig('koi134_ttvdist.png',dpi=300,bbox_inches='tight')