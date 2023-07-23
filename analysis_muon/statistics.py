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
    popt2,pcov2 = op.curve_fit(func,dis[0:25], ratio_list[0:25], p0=[50,2])
    # error_list = np.sqrt(np.diag(pcov2))
    if func == FD_func:
        return popt2[0]
    if func == gauss:
        return popt2[0] + np.sqrt(popt2[1]*np.log(2))

l1_dis_list_tmp = np.load('/lustre/neutrino/weizhenyu/analysis_muon/data_npy/l1_dis_list_tmp.npy')
l1_e_list_tmp = np.load('/lustre/neutrino/weizhenyu/analysis_muon/data_npy/l1_e_list_tmp.npy')

l1_dis_list = np.load('/lustre/neutrino/weizhenyu/analysis_muon/data_npy/l1_dis_list_unique.npy')
l1_e_list = np.load('/lustre/neutrino/weizhenyu/analysis_muon/data_npy/l1_e_list_unique.npy')

dis_list_all = np.load('/lustre/neutrino/weizhenyu/analysis_muon/data_npy/total_dis_list_all.npy')
e_list_all = np.load('/lustre/neutrino/weizhenyu/analysis_muon/data_npy/total_e_list_all.npy')

not_l1_dis_list_tmp = np.load('/lustre/neutrino/weizhenyu/analysis_muon/data_npy/not_l1_dis_list_tmp.npy')
not_l1_e_list_tmp = np.load('/lustre/neutrino/weizhenyu/analysis_muon/data_npy/not_l1_e_list_tmp.npy')
not_l1_dis_list = np.load('/lustre/neutrino/weizhenyu/analysis_muon/data_npy/not_l1_dis_list.npy')
not_l1_e_list = np.load('/lustre/neutrino/weizhenyu/analysis_muon/data_npy/not_l1_e_list.npy')

# l1_total = [*l1_dis_list, *not_l1_dis_list]
l1_total = [*dis_list_all,*l1_dis_list]
e_total = [*e_list_all, *l1_e_list]
# %%
plt.figure(dpi = 800)
range_ = [0,200]
bins_ = 50
# bins_ = np.logspace(np.log10(range_[0]),np.log10(range_[1]),num=20)
plt.hist(l1_total,range=range_, bins=bins_,histtype='step',label='all hDOMs')
plt.hist(l1_dis_list,range=range_, bins=bins_,histtype='step',label='triggered hDOM')
# plt.hist(not_l1_dis_list,range=range_, bins=bins_,histtype='step',label='not-triggered hDOM')
plt.legend()
plt.xlabel('Distance [m]')
# plt.xlabel('CosTheta')
plt.ylabel('number of hDOM')   
plt.semilogy()
#%%
counts_total,tick_ = np.histogram(l1_total,range=range_, bins=bins_)
counts_L1, tick_ = np.histogram(l1_dis_list,range=range_, bins=bins_)
# counts_not_L1, tick_ = np.histogram(not_l1_dis_list,range=range_, bins=bins_)
dis = (tick_[:-1] + tick_[1:])/2
error_L1 = (np.power(counts_L1, 0.5) / counts_total)**2
error_L1 += (np.divide(counts_L1,np.power(counts_total, 2)) * np.power(counts_total, 0.5))**2
error_L1 = np.power(error_L1,0.5)
plt.errorbar(x = dis, y = np.divide(counts_L1,counts_total), yerr=error_L1,label = 'L1 Ratio', fmt='.')
# error_notL1 = (np.power(counts_not_L1, 0.5) / counts_total)**2
# error_notL1 += (np.divide(counts_not_L1,np.power(counts_total, 2)) * np.power(counts_total, 0.5))**2
# error_notL1 = np.power(error_notL1,0.5)
# plt.errorbar(x = dis, y = np.divide(counts_not_L1,counts_total), yerr= error_notL1, label = 'not L1 Ratio', fmt='.')
plt.xlabel('Distance [m]')
# plt.xlabel('CosTheta')
plt.ylim([0,1])
plt.legend()
plt.ylabel('Ratio')
# %%
hdom_trigger_info = pd.DataFrame({'e_l1':l1_e_list,'dis_l1':l1_dis_list})
# hdom_total_trigger_info = pd.DataFrame({'e_total': e_list_all, 'dis_total': dis_list_all})

# hdom_trigger_info = pd.DataFrame({'e_l1':l1_e_list,'dis_l1':l1_dis_list})
hdom_not_trigger_info = pd.DataFrame({'e_l1_not': not_l1_e_list, 'dis_l1_not': not_l1_dis_list})

#%%
energy_muon_ = []
effective_radius_ = []
effective_radius_upper = []
effective_radius_lower = []
value_ratio = []
for e in range(12):
    
    e_upper = np.logspace(np.log10(1),np.log10(100),num=13)[e+1]
    e_upper *= 1000

    e_lower = np.logspace(np.log10(1),np.log10(100),num=13)[e]
    e_lower *= 1000

    hdom_e = hdom_trigger_info.loc[hdom_trigger_info['e_l1'] < e_upper]
    hdom_e = hdom_e.loc[hdom_e['e_l1'] > e_lower]
    l1_dis_list = hdom_e['dis_l1'].to_list()
    energy_muon_.append(np.mean(hdom_e['e_l1']))

    # hdom_e_total = hdom_total_trigger_info.loc[hdom_total_trigger_info['e_total'] < e_upper]
    # hdom_e_total = hdom_e_total.loc[hdom_e_total['e_total'] > e_lower]
    # l1_total = hdom_e_total['dis_total'].to_list()

    hdom_e = hdom_not_trigger_info.loc[hdom_not_trigger_info['e_l1_not'] < e_upper]
    hdom_e = hdom_e.loc[hdom_e['e_l1_not'] > e_lower]
    not_l1_dis_list = hdom_e['dis_l1_not'].to_list()
    l1_total = [*l1_dis_list,*not_l1_dis_list]

    counts_total,tick_ = np.histogram(l1_total,range=range_, bins=bins_)
    counts_L1, tick_ = np.histogram(l1_dis_list,range=range_, bins=bins_)
    # counts_not_L1, tick_ = np.histogram(not_l1_dis_list,range=range_, bins=bins_)
    dis = (tick_[:-1] + tick_[1:])/2
    error_L1 = (np.power(counts_L1, 0.5) / counts_total)**2
    error_L1 += (np.divide(counts_L1,np.power(counts_total, 2)) * np.power(counts_total, 0.5))**2
    error_L1 = np.power(error_L1,0.5)
    plt.errorbar(x = dis, y = np.divide(counts_L1,counts_total), yerr=error_L1,label = 'L1 Ratio', fmt='.')
    effective_radius_.append(calculate_eff_radius(FD_func,dis, np.divide(counts_L1,counts_total)))
    effective_radius_upper.append(calculate_eff_radius(FD_func,dis, (np.divide(counts_L1,counts_total) + error_L1)))
    effective_radius_lower.append(calculate_eff_radius(FD_func, dis, (np.divide(counts_L1,counts_total) - error_L1)))
    popt2,pcov2 = op.curve_fit(FD_func,dis[0:40], np.divide(counts_L1,counts_total)[0:40], p0=[50,2])
    plt.plot(np.linspace(0,80,80), FD_func(np.linspace(0,80,80), *(popt2)),label = 'fitting')
    value_ratio.append(FD_func(np.linspace(0,80,80), *(popt2)))
    # error_notL1 = (np.power(counts_not_L1, 0.5) / counts_total)**2
    # error_notL1 += (np.divide(counts_not_L1,np.power(counts_total, 2)) * np.power(counts_total, 0.5))**2
    # error_notL1 = np.power(error_notL1,0.5)
    # plt.errorbar(x = dis, y = np.divide(counts_not_L1,counts_total), yerr= error_notL1, label = 'not L1 Ratio', fmt='.')
plt.xlabel('Distance [m]')
# plt.xlabel('CosTheta')
plt.ylim([0,1])
# plt.legend()
plt.ylabel('Ratio')
energy_muon_ = np.divide(energy_muon_,1000)

#%%
plt.figure(dpi = 800)
x_ticks = np.linspace(0,80,10)
y_ticks = np.log(energy_muon_)
ax = sns.heatmap(value_ratio,cmap = "YlGnBu")#, xticklabels=x_ticks, yticklabels=y_ticks)
# plt.xticks(ticks = np.arange(0,120,24))#, labels=['A', 'B', 'C'])
# plt.yticks([])#, labels=['X', 'Y', 'Z'])
line_y = 11.5 / np.log(energy_muon_)[-1] * np.log(energy_muon_)
line_x = effective_radius_
for x_id in range(len(effective_radius_)):
    x_ = line_x[x_id]  # X-coordinates of the line
    y_ = line_y[x_id]
    x__ = np.ones(10) * x_
    y__ = np.linspace(y_-0.5, y_+0.5,10)
    ax.plot(x__, y__,  color = 'red', linestyle = '--', linewidth =  0.8)
    ax.scatter(x_, y_,  color = 'red', s = 2)
ax.set_xlabel('Distance [m]')
ax.set_ylabel(r"$\nu_{\mu}$ Energy [TeV]")

# %%

plt.errorbar(x = energy_muon_, y = effective_radius_,yerr=[np.subtract(effective_radius_upper,effective_radius_), np.subtract(effective_radius_, effective_radius_lower)],label = 'L1 Ratio', fmt='.')
plt.xlabel('Energy [TeV]')
# plt.xlabel('CosTheta')
plt.ylim([50,120])
# plt.legend()
plt.ylabel('Effective Radius [m]')
# plt.semilogy()
# plt.plot(energy_muon_, effective_radius_upper)
# plt.plot(energy_muon_, effective_radius_lower)
# plt.semilogx()

# %%
efficiency_list = []
energy_muon_ = []
for e in range(12):
    
    e_upper = np.logspace(np.log10(1),np.log10(100),num=13)[e+1]
    e_upper *= 1000

    e_lower = np.logspace(np.log10(1),np.log10(100),num=13)[e]
    e_lower *= 1000

    hdom_e = hdom_trigger_info.loc[hdom_trigger_info['e_l1'] < e_upper]
    hdom_e = hdom_e.loc[hdom_e['e_l1'] > e_lower]
    l1_dis_list = hdom_e['dis_l1'].to_list()
    energy_muon_.append(np.mean(hdom_e['e_l1']))

    hdom_e = hdom_not_trigger_info.loc[hdom_not_trigger_info['e_l1_not'] < e_upper]
    hdom_e = hdom_e.loc[hdom_e['e_l1_not'] > e_lower]
    not_l1_dis_list = hdom_e['dis_l1_not'].to_list()
    l1_total = [*l1_dis_list,*not_l1_dis_list]

    l1_30 = len(np.array(l1_dis_list)[np.array(l1_dis_list) < 30]) 
    err_l1_30 = l1_30 ** (-0.5)
    l1_tot_30 = len(np.array(l1_total)[np.array(l1_total)<30])
    err_l1_tot_30 = l1_tot_30 ** (-0.5)
    err_eff = np.power(err_l1_30**2 + err_l1_tot_30**2,0.5)
    
    efficiency_list.append(l1_30 / l1_tot_30)

# %%
plt.figure(dpi = 800)
plt.plot(energy_muon_, efficiency_list, linestyle = '--', color = 'blue')
# plt.scatter(energy_muon_, efficiency_list)
plt.errorbar(x = energy_muon_, y = efficiency_list, yerr = err_eff*np.array(efficiency_list), fmt='.', color = 'blue')
plt.xlabel(r"$\nu_{\mu}$ Neutrino Energy [GeV]")
plt.ylabel(r"Triggered hDOM / hit hDOM")
plt.semilogx()
# %%

