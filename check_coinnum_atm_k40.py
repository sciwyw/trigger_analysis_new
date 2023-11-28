#%%
import sys
sys.path.append("/lustre/neutrino/weizhenyu/trigger/analysis/analyze-trident-mc-master") # 改为自己的git的pull的目录

import os
from pathlib import Path
from typing import Union
import pandas as pd
import numpy as np
import math
from math import pi, cos
import glob
import argparse

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
font = {'family': 'serif',
         'weight': 'normal', 'size': 12}
plt.rc('font', **font)
import domhits_analysis as da
import importlib
importlib.reload(da)
import config as con
con.mcevents_suffix = "mc_events_with_prim.json"
con.raw_power_index = -2
con.sample_area = pi * (2000**2 + 285**2)

plt.clf()
#/lustre/neutrino/zhangfuyudi/hailing/cascade/sy_data/nue_100tev_params/10k_p13/
from matplotlib import font_manager
for i in ['', 'bd', 'bi', 'i']:
    font_manager.fontManager.addfont('/home/zhiwei/fonts/times{}.ttf'.format(i))

plt.rcParams.update({
    'figure.figsize': (8, 5.5),
    'font.size': 20,
    'font.family': 'Times New Roman'
})
#%%

# if __name__ == '__main__':
#     # get samples
pmthits = pd.read_csv('./pmthits.csv').set_index('showerId')
pmthits = pmthits[pmthits.z0==285] # 只有最上面一层的


domhits, coincidence = da.get_coincidence_level_with_showerId(pmthits)
domhits.to_csv('tmp.csv', index_label=['showerId', 'DomId'])
print(domhits)

freq_mean = coincidence.groupby('coincidence_level').mean()
freq_std = coincidence.groupby('coincidence_level').std()
fig, ax = plt.subplots(figsize=(10, 8), dpi=400)
np.save('/lustre/neutrino/weizhenyu/trigger/analysis/result/x_muon_coin.npy',freq_mean.index.unique() )
np.save('/lustre/neutrino/weizhenyu/trigger/analysis/result/y_muon_coin.npy',freq_mean['weight'] )
np.save('/lustre/neutrino/weizhenyu/trigger/analysis/result/yerr_muon_coin.npy',freq_std['weight'] )

#%%
range_ = [0.5,30.5]
bins_ = 31
x_muon = np.load('/lustre/neutrino/weizhenyu/trigger/analysis/result/x_muon_coin.npy')
y_muon = np.load('/lustre/neutrino/weizhenyu/trigger/analysis/result/y_muon_coin.npy')
yerr_muon = np.load('/lustre/neutrino/weizhenyu/trigger/analysis/result/yerr_muon_coin.npy')
coin_pmt_k40 = [1,2,3,4,5,6]
count_coin_k40 = np.array([234600,5851, 306, 22, 4, 1])
# count_coin_k40[:] = count_coin_k40[:]/3.91
data_fake_k40 = []
for i in range(len([1,2,3,4,5,6])):
    data_fake_k40.extend(np.ones(count_coin_k40[i]) * coin_pmt_k40[i])
err_coin_k40 = np.sqrt(count_coin_k40)
counts_k40, ticks = np.histogram(data_fake_k40, density=False,range = range_,bins = bins_)
error_K40 = np.power(count_coin_k40,-0.5) * counts_k40[0:6]
x_axis = (ticks[1:-1] +ticks[0:-2])/2
x_axis = x_axis[0:6]
#%%
plt.figure(dpi = 800)
plt.errorbar(x = x_axis , y = counts_k40[0:6]/3.91, yerr=error_K40[0:6]/3.91, fmt = '.',color = 'orange', markersize=3,label = 'K40')

plt.hist(data_fake_k40,histtype='step', density=False,range = range_,bins = bins_,color = 'orange',weights=np.ones(len(data_fake_k40))/3.91)

plt.errorbar(x_muon, y_muon, yerr=yerr_muon, fmt='.', markersize=3,color = 'green',label = 'Atm Muon')
plt.bar(x = x_muon, height=y_muon, width=np.ones(31),color = 'green')
plt.bar(x = np.linspace(x_muon[0]+0.1,x_muon[-1]-0.1,bins_), height=y_muon*0.9, width=np.ones(31),color = 'white')
plt.semilogy()
plt.legend()
plt.ylabel('Rate [Hz]')
plt.xlabel('Coincidence PMT Number')
plt.axvline(np.quantile(data_fake_k40,0.98)-0.5,color='red',ls='--')
plt.text(np.quantile(data_fake_k40,0.98),3000,'98% K40',color='red')
# plt.set_title('Coincidence Level at Upper Surface (~2615m deep)')
plt.show()
# plt.tight_layout()
# fig.savefig('./save/coincidence_level_uppersurface.pdf')

# %%
