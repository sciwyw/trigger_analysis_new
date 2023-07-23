#%%
import numpy as np
import pandas as pd
import uproot
import os
from matplotlib import rcParams
import sys
import matplotlib.pyplot as plt
from scipy import optimize as op
from cProfile import label
import matplotlib
import seaborn as sns
from mpl_toolkits import mplot3d
import scipy.interpolate
from scipy import optimize as op
import json
import fnmatch
from tqdm import tqdm
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

l2_e_list = np.load('/lustre/neutrino/weizhenyu/analysis_trigger_L2/data_npy/l2_energy_list_cqc_100_2_10_60_330.npy')
l2_e_list_large = np.load('/lustre/neutrino/weizhenyu/analysis_trigger_L2/data_npy/l2_energy_list_cqc_100_2_10_120_660.npy')
total_l2_list = np.load('/lustre/neutrino/weizhenyu/analysis_trigger_L2/data_npy/l2_energy_list_cqc_100_total.npy')
# %%
plt.figure(dpi = 800)
range_ = [1000,100000]
bins_ = np.logspace(np.log10(range_[0]),np.log10(range_[1]),num=20)
counts_total,tick_ = np.histogram(total_l2_list,range=range_, bins=bins_)
counts_L2, tick_ = np.histogram(l2_e_list,range=range_, bins=bins_)
counts_L2_large, tick_ = np.histogram(l2_e_list_large,range=range_, bins=bins_)
# counts_not_L2, tick_ = np.histogram(not_l2_e_list,range=range_, bins=bins_)
energy_ = (tick_[:-1] + tick_[1:])/2
error_L2 = (np.power(counts_L2, 0.5) / counts_total)**2
error_L2 += (np.divide(counts_L2,np.power(counts_total, 2)) * np.power(counts_total, 0.5))**2
error_L2 = np.power(error_L2,0.5)
error_L2_large = (np.power(counts_L2_large, 0.5) / counts_total)**2
error_L2_large += (np.divide(counts_L2_large,np.power(counts_total, 2)) * np.power(counts_total, 0.5))**2
error_L2_large = np.power(error_L2_large,0.5)
plt.errorbar(x = energy_, y = np.divide(counts_L2_large,counts_total), yerr=error_L2,label = '120m interval', fmt='.')
plt.errorbar(x = energy_, y = np.divide(counts_L2,counts_total), yerr=error_L2,label = '60m interval', fmt='.')
plt.semilogx()
plt.xlabel('neutrino energy [GeV]')
# plt.xlabel('CosTheta')
plt.ylim([0,1])
plt.legend()
plt.ylabel('Efficiency')
# %%
