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

def FD_func(x,*param):
    return 1 / (1 + np.exp((x-param[0])/param[1]))
def gauss(x, *param):
    return np.exp(-(x - param[0])**2 / param[1])
# popt2,pcov2 = op.curve_fit(FD_func,dis[0:25], np.divide(counts_L1,counts_total)[0:25], p0=[50,2])
# error_list = np.sqrt(np.diag(pcov2))
# plt.plot(np.linspace(0,100,100), FD_func(np.linspace(0,100,100), *popt2))
def calculate_eff_radius(func,dis, ratio_list):
    popt2,pcov2 = op.curve_fit(func,dis, ratio_list, p0=[50,2])
    # error_list = np.sqrt(np.diag(pcov2))
    if func == FD_func:
        return popt2[0]
    if func == gauss:
        return popt2[0] + np.sqrt(popt2[1]*np.log(2))
#%%
range_ = [0,250]
bins_ = 50
energy_muon_ = []
effective_radius_ = []
effective_radius_upper = []
effective_radius_lower = []
value_ratio = []
for e_id in range(5):
    l1_dis_list = np.load('/lustre/neutrino/weizhenyu/analysis_shower/data_npy/l1_dis_list_{}.npy'.format(e_id))
    total_dis_list = np.load('/lustre/neutrino/weizhenyu/analysis_shower/data_npy/total_dis_list_{}.npy'.format(e_id))
    # plt.hist(l1_dis_list, histtype = 'step', range= range_, bins = bins_)
    # plt.hist(total_dis_list, histtype = 'step', range = range_,bins= bins_)
    counts_total,tick_ = np.histogram(total_dis_list,range=range_, bins=bins_)
    counts_L1, tick_ = np.histogram(l1_dis_list,range=range_, bins=bins_)
    # counts_not_L1, tick_ = np.histogram(not_l1_dis_list,range=range_, bins=bins_)
    dis = (tick_[:-1] + tick_[1:])/2
    popt2,pcov2 = op.curve_fit(FD_func,dis, np.divide(counts_L1,counts_total), p0=[80,2])
    effective_radius_.append(calculate_eff_radius(FD_func,dis, np.divide(counts_L1,counts_total)))
    value_ratio.append(FD_func(np.linspace(0,250,250), *(popt2)))
# %%
plt.figure(dpi = 800)
x_ticks = np.linspace(0,80,10)
y_ticks = np.log(energy_muon_)
ax = sns.heatmap(value_ratio,cmap = "YlGnBu")#, xticklabels=x_ticks, yticklabels=y_ticks)
# plt.xticks(ticks = np.arange(0,120,24))#, labels=['A', 'B', 'C'])
# plt.yticks([])#, labels=['X', 'Y', 'Z'])
line_y = [0.5,1.5,2.5,3.5,4.5]#11.5 / np.log(energy_muon_)[-1] * np.log(energy_muon_)
line_x = effective_radius_
for x_id in range(len(effective_radius_)):
    x_ = line_x[x_id]  # X-coordinates of the line
    y_ = line_y[x_id]
    x__ = np.ones(10) * x_
    y__ = np.linspace(y_-0.5, y_+0.5,10)
    ax.plot(x__, y__,  color = 'red', linestyle = '--', linewidth =  0.8)
    ax.scatter(x_, y_,  color = 'red', s = 2)
ax.set_xlabel('Distance [m]')
ax.set_ylabel(r"$\nu_{E}$ Energy [TeV]")
# %%
