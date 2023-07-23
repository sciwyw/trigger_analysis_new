#%%
import numpy as np
import pandas as pd
import uproot
import os
import sys
import scipy.stats as st
import matplotlib.pyplot as plt
from scipy import optimize as op
from cProfile import label
import matplotlib
from mpl_toolkits import mplot3d
import scipy.interpolate
# import seaborn as sns
plt.clf()
from matplotlib import font_manager
for i in ['', 'bd', 'bi', 'i']:
    font_manager.fontManager.addfont('/home/zhiwei/fonts/times{}.ttf'.format(i))

plt.rcParams.update({
    'figure.figsize': (8, 5.5),
    'font.size': 20,
    'font.family': 'Times New Roman'
})
rate_dark_count_pmt = 300 * 31
rate_dark_count_sipm = 30000 * 24
time_total = 5 

def obtain_trigger(data:pd.DataFrame):
    diff_data_t = np.array(data["t(ns)"][coin_num-1:len(data)]) - np.array(data["t(ns)"][0:len(data)-(coin_num-1)]) # The last term subtracts the first term
    trigger_array = diff_data_t < time_window
    for i0 in range((coin_num-1)):
        trigger_array = np.append(trigger_array, False) # treat the last value as the same with the one before
    True_list = []
    for i in range(len(trigger_array)):
        if trigger_array[i] == True:
            for j0 in range((coin_num-1)):
                True_list.append(i+j0+1)
    for j in True_list:
        trigger_array[j] = True

    data_trigger = data.loc[trigger_array]
    return data_trigger


#%%
time_window_list = np.arange(2,30,2)
coin_list = [2]#range(2,6,1)
rate_list = []

for time_window in time_window_list:
    for coin_num in coin_list:
        x = np.linspace(0,time_window,int(1000 * time_window))
        rate = st.gamma.pdf(x, coin_num - 1, scale = (10**9 / rate_dark_count_pmt))
        p = np.sum(rate) * (x[1] - x[0]) * rate_dark_count_pmt
        rate_list = np.append(rate_list,p) 

#%%
plt.figure(dpi = 800)
plt.scatter(time_window_list,rate_list,label = "coincidence number = 2", color = 'blue')
plt.plot(time_window_list,rate_list, linestyle = '--',color = 'blue')
plt.legend()
plt.xlabel("time window [ns]")
plt.ylabel("Trigger rate [Hz]")
#%%
time_window_list = [10]#np.arange(2,30,2)
coin_list = range(2,9,1)
rate_list = []

for time_window in time_window_list:
    for coin_num in coin_list:
        x = np.linspace(0,time_window,int(1000 * time_window))
        rate = st.gamma.pdf(x, coin_num - 1, scale = (10**9 / rate_dark_count_pmt))
        p = np.sum(rate) * (x[1] - x[0]) * rate_dark_count_pmt
        rate_list = np.append(rate_list,p) 

#%%
plt.figure(dpi = 800)
plt.scatter(coin_list,rate_list,label = "time window = 10 ns", color = 'blue')
plt.semilogy()
plt.plot(coin_list,rate_list,linestyle = '--', color = 'blue')
plt.legend()
plt.xlabel("coincidence photon number ")
# plt.ylabel("Trigger rate (Hz)")
#%%
rate_list_ = np.reshape(rate_list,[14,4])
plt.figure(figsize = (14,6.5),dpi=800)
ax = sns.heatmap(np.array(rate_list_), cmap = "YlGnBu", yticklabels=time_window_list, xticklabels=coin_num)
ax.set_title('Trigger Rate')
ax.set_ylabel('time window(ns)')
ax.set_xlabel('coincidence number')
#ax.set_xticks(np.arange(1.5,58.5,6))
#plt.xticks(np.arange(1.5,58.5,6))
# filename = './heatmap'
# plt.savefig(filename+'.jpg', bbox_inches='tight', dpi=800)
plt.savefig('/lustre/neutrino/weizhenyu/trigger/analysis/dark_rate_trigger.png',bbox_inches='tight', dpi=800)
# %%
