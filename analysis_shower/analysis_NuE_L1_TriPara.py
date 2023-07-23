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
import json

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

def qe_module(hit: pd.DataFrame):
    #quantum efficiency
    # 1eV ~ 1240nm
    qe_list = pd.read_csv('/lustre/neutrino/weizhenyu/hdom_signal_simulation/analysis/qe.csv')
    qe_func = scipy.interpolate.interp1d(qe_list['wl'],qe_list['qe'],kind='linear')
    wl_hit = 1240 / np.array(hit["e0"])
    qe_photon = qe_func(wl_hit) / 100
    sample = np.random.sample(len(hit))
    survival_bool = qe_photon > sample
    hit_qe = hit[survival_bool]
    #print(survival_bool)
    hit_qe = hit_qe.reset_index(drop = True)

    return hit_qe

def obtain_trigger(data:pd.DataFrame):
    diff_data_t = np.array(data["t0"][coin_num-1:len(data)]) - np.array(data["t0"][0:len(data)-(coin_num-1)]) # The last term subtracts the first term
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

    data_trigger = data.iloc[trigger_array]
    
    if len(data_trigger) != 0:
        idx_drop = []
        for hits_id in range(len(data_trigger)-3):
            if data_trigger.iloc[hits_id+1]["t0"] - data_trigger.iloc[hits_id]["t0"] < time_window:
                if data_trigger.iloc[hits_id+1]["PmtId"] == data_trigger.iloc[hits_id]["PmtId"]:
                    idx_drop.append(data_trigger.index[hits_id])
                    idx_drop.append(data_trigger.index[hits_id+1])
        data_trigger.drop(idx_drop,inplace = True)
            
        t_trigger_start = data_trigger.iloc[0]["t0"] - 300
        t_trigger_end = data_trigger.iloc[-1]["t0"] + 700
        data_new = data.loc[data["t0"] > t_trigger_start]
        data_new = data_new.loc[data_new["t0"] < t_trigger_end]
        if len(data_trigger) != len(data_new):
            data_trigger = data_trigger.append(data_new)
        data_trigger.drop_duplicates(subset=None, keep= 'first', inplace= True)
    

    return data_trigger
def Is_L1_trigger(data:pd.DataFrame, coin_num):
    if len(data) < coin_num:
        return False
    diff_data_t = np.array(data["t0"][coin_num-1:len(data)]) - np.array(data["t0"][0:len(data)-(coin_num-1)]) # The last term subtracts the first term
    
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
    
    data_trigger = data.iloc[trigger_array]
    if len(np.unique(data_trigger['PmtId'].to_numpy()))<2:
        return False
    else:
        return True

def Is_L1_trigger_loose(data:pd.DataFrame, coin_num):
    diff_data_t = np.array(data["t0"][coin_num-1:len(data)]) - np.array(data["t0"][0:len(data)-(coin_num-1)]) # The last term subtracts the first term
    
    trigger_array = diff_data_t < time_window
    if True in trigger_array:
        return True
    else:
        return False

geo = pd.read_csv("/lustre/neutrino/weizhenyu/trigger/L1_L2/trident/data/penrose.csv",names=['x','y','z'])

file_name_100TeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/sy_data/nue_100tev_params/10k_p13/job_0/'
file_name_10TeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/jobs_10tev_nue_10k_p1/job_0/'
file_json_100TeV = file_name_100TeV + 'mc_events.json'
file_json_10TeV = file_name_10TeV + 'mc_events.json'

#%%
with open(file_json_100TeV) as f:
    mc_events = json.load(f) 

file = uproot.open(file_name_100TeV + 'data/data.root')
tree = file[file.keys()[2]]
hits =tree.arrays(library = 'pd')

#%%
time_window = 10
range_ = [5,300]
bins_ = 50 

coin_num_list = [2,3,4,5,6,7 ]
r_eff_list = []
errorbar_lower_list = []
errorbar_upper_list = []
for coin_num in coin_num_list:
    not_L1_dis_5m = []
    L1_dis_5m = []
    dis_list_5m = []
    for event_id in range(1000):#range(1000):
        hit_event = hits.loc[event_id]
        hit_event = qe_module(hit_event)
        hit_event.sort_values(by=["t0"],inplace = True)
        hit_event.reset_index(drop = True, inplace = True)

        pos_vertex = [mc_events[event_id]['particles_in'][0]['x'],mc_events[event_id]['particles_in'][0]['y'],mc_events[event_id]['particles_in'][0]['z']]
        dir_vertex = [mc_events[event_id]['particles_in'][0]['px'],mc_events[event_id]['particles_in'][0]['py'],mc_events[event_id]['particles_in'][0]['pz']]
        norm_dir = (dir_vertex[0]**2 + dir_vertex[1]**2 + dir_vertex[2]**2)**0.5


        dx_5m = geo.iloc[hit_event['DomId']]['x'] - (pos_vertex[0] + dir_vertex[0] / norm_dir * 5)
        dy_5m = geo.iloc[hit_event['DomId']]['y'] - (pos_vertex[1] + dir_vertex[1] / norm_dir * 5)
        dz_5m = geo.iloc[hit_event['DomId']]['z'] - (pos_vertex[2] + dir_vertex[2] / norm_dir * 5)


        dis_5m = np.power(np.power(dx_5m,2) + np.power(dy_5m,2) + np.power(dz_5m,2), 0.5)
        for domid in np.unique(hit_event['DomId']):
            dom_hit = hit_event.loc[ hit_event['DomId'] == domid]
            dom_photon_number = len(dom_hit)
            dom_pmt_number = len(np.unique(dom_hit['PmtId']))
            dom_dis_5m = dis_5m.iloc[dom_hit.index[0]]
            dis_list_5m.append(dom_dis_5m)
            if len(dom_hit)<2:
                
                not_L1_dis_5m.append(dom_dis_5m)
                        
            else:
                if Is_L1_trigger(dom_hit, coin_num):
                    
                    L1_dis_5m.append(dom_dis_5m)
                    
                    
                else:
                    
                    not_L1_dis_5m.append(dom_dis_5m)
    
    counts_total,tick_ = np.histogram(dis_list_5m,range=range_, bins=bins_)
    counts_L1, tick_ = np.histogram(L1_dis_5m,range=range_, bins=bins_)
    counts_not_L1, tick_ = np.histogram(not_L1_dis_5m,range=range_, bins=bins_)
    dis = (tick_[:-1] + tick_[1:])/2

    error_L1 = (np.power(counts_L1, 0.5) / counts_total)**2
    error_L1 += (np.divide(counts_L1,np.power(counts_total, 2)) * np.power(counts_total, 0.5))**2
    error_L1 = np.power(error_L1,0.5)

    ratio_func = scipy.interpolate.interp1d(x=np.divide(counts_L1,counts_total), y = dis,kind='linear')
    ratio_func_lower = scipy.interpolate.interp1d(x=np.divide(counts_L1,counts_total) - error_L1, y = dis,kind='linear')
    ratio_func_upper = scipy.interpolate.interp1d(x=np.divide(counts_L1,counts_total) + error_L1, y = dis,kind='linear')

    errorbar_lower = ratio_func(0.5) - ratio_func_lower(0.5)
    errorbar_lower_list.append(errorbar_lower)
    errorbar_upper = ratio_func(0.5) - ratio_func_upper(0.5)
    errorbar_upper_list.append(errorbar_upper)

    r_eff = ratio_func(0.5)
    r_eff_list.append(r_eff)
                
                
plt.figure(dpi = 800)
plt.errorbar(x = coin_num_list, y = r_eff_list, yerr = np.array([errorbar_lower_list,np.multiply(-1,errorbar_upper_list)]), fmt='.', capsize=3, color = 'blue', label = 'MC results')
plt.plot(coin_num_list,  r_eff_list, color = 'blue', linestyle = '--')
plt.legend()
plt.xlabel('Coincidence Photon Number')
plt.ylabel('Effective Radius [m]')
#%%
time_window_list = [10,20,30,40,50,60]
r_eff_list = []
errorbar_lower_list = []
errorbar_upper_list = []
coin_num = 2
for time_window in time_window_list:
    not_L1_dis_5m = []
    L1_dis_5m = []
    dis_list_5m = []
    for event_id in range(400):#range(1000):
        hit_event = hits.loc[event_id]
        hit_event = qe_module(hit_event)
        hit_event.sort_values(by=["t0"],inplace = True)
        hit_event.reset_index(drop = True, inplace = True)

        pos_vertex = [mc_events[event_id]['particles_in'][0]['x'],mc_events[event_id]['particles_in'][0]['y'],mc_events[event_id]['particles_in'][0]['z']]
        dir_vertex = [mc_events[event_id]['particles_in'][0]['px'],mc_events[event_id]['particles_in'][0]['py'],mc_events[event_id]['particles_in'][0]['pz']]
        norm_dir = (dir_vertex[0]**2 + dir_vertex[1]**2 + dir_vertex[2]**2)**0.5


        dx_5m = geo.iloc[hit_event['DomId']]['x'] - (pos_vertex[0] + dir_vertex[0] / norm_dir * 5)
        dy_5m = geo.iloc[hit_event['DomId']]['y'] - (pos_vertex[1] + dir_vertex[1] / norm_dir * 5)
        dz_5m = geo.iloc[hit_event['DomId']]['z'] - (pos_vertex[2] + dir_vertex[2] / norm_dir * 5)


        dis_5m = np.power(np.power(dx_5m,2) + np.power(dy_5m,2) + np.power(dz_5m,2), 0.5)
        for domid in np.unique(hit_event['DomId']):
            dom_hit = hit_event.loc[ hit_event['DomId'] == domid]
            dom_photon_number = len(dom_hit)
            dom_pmt_number = len(np.unique(dom_hit['PmtId']))
            dom_dis_5m = dis_5m.iloc[dom_hit.index[0]]
            dis_list_5m.append(dom_dis_5m)
            if len(dom_hit)<2:
                
                not_L1_dis_5m.append(dom_dis_5m)
                        
            else:
                if Is_L1_trigger(dom_hit, coin_num):
                    
                    L1_dis_5m.append(dom_dis_5m)
                    
                    
                else:
                    
                    not_L1_dis_5m.append(dom_dis_5m)
    
    counts_total,tick_ = np.histogram(dis_list_5m,range=range_, bins=bins_)
    counts_L1, tick_ = np.histogram(L1_dis_5m,range=range_, bins=bins_)
    counts_not_L1, tick_ = np.histogram(not_L1_dis_5m,range=range_, bins=bins_)
    dis = (tick_[:-1] + tick_[1:])/2
    
    error_L1 = (np.power(counts_L1, 0.5) / counts_total)**2
    error_L1 += (np.divide(counts_L1,np.power(counts_total, 2)) * np.power(counts_total, 0.5))**2
    error_L1 = np.power(error_L1,0.5)

    ratio_func = scipy.interpolate.interp1d(x=np.divide(counts_L1,counts_total), y = dis,kind='linear')
    ratio_func_lower = scipy.interpolate.interp1d(x=np.divide(counts_L1,counts_total) - error_L1, y = dis,kind='linear')
    ratio_func_upper = scipy.interpolate.interp1d(x=np.divide(counts_L1,counts_total) + error_L1, y = dis,kind='linear')

    errorbar_lower = ratio_func(0.5) - ratio_func_lower(0.5)
    errorbar_lower_list.append(errorbar_lower)
    errorbar_upper = ratio_func(0.5) - ratio_func_upper(0.5)
    errorbar_upper_list.append(errorbar_upper)

    r_eff = ratio_func(0.5)
    r_eff_list.append(r_eff)
                
                
plt.figure(dpi = 800)
plt.errorbar(x = time_window_list, y = r_eff_list, yerr = np.array([errorbar_lower_list,np.multiply(-1,errorbar_upper_list)]), fmt='.', capsize=3, color = 'blue', label = 'MC results')
plt.plot(time_window_list,  r_eff_list, color = 'blue', linestyle = '--')
plt.legend()
plt.ylim([80,160])
plt.xlabel('Time Window [ns]')
plt.ylabel('Effective Radius [m]')
# %%
