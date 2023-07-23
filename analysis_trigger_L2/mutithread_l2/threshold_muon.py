#%%
# mc_events particles_in unit GeV
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
import fnmatch
from tqdm import tqdm
import multiprocessing
import concurrent.futures
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
time_window = 10
coin_num = 2
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
def Is_L1_trigger(data:pd.DataFrame):
    if len(data) < coin_num:
        return False
    else:
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

def Is_L1_trigger_loose(data:pd.DataFrame):
    diff_data_t = np.array(data["t0"][coin_num-1:len(data)]) - np.array(data["t0"][0:len(data)-(coin_num-1)]) # The last term subtracts the first term
    
    trigger_array = diff_data_t < time_window
    if True in trigger_array:
        return True
    else:
        return False

def check_csv_in_list(string_list):
    for string in string_list:
        if fnmatch.fnmatch(string, '*.csv'):
            return True
    return False
def check_json_in_list(string_list):
    for string in string_list:
        if fnmatch.fnmatch(string, '*.json'):
            return True
    return False
def calculate_distance(hit_event: pd.DataFrame, eventinfo: dict):
    dir_vertex = [eventinfo['particles_in'][0]['px'],eventinfo['particles_in'][0]['py'],eventinfo['particles_in'][0]['pz']]

    pos_vertex = [eventinfo['particles_in'][0]['x'],eventinfo['particles_in'][0]['y'],eventinfo['particles_in'][0]['z']]

    dx = geo.iloc[hit_event['DomId']]['x'] - pos_vertex[0]
    dy = geo.iloc[hit_event['DomId']]['y'] - pos_vertex[1]
    dz = geo.iloc[hit_event['DomId']]['z'] - pos_vertex[2]

    norm_dir = (dir_vertex[0]**2 + dir_vertex[1]**2 + dir_vertex[2]**2)**0.5
    para_t = (dir_vertex[0] * dx + dir_vertex[1] * dy + dir_vertex[2] * dz) / norm_dir ** 2
    dis_line = np.power(np.power(dx - dir_vertex[0] * para_t,2) + np.power(dy - dir_vertex[1] * para_t,2) + np.power(dz - dir_vertex[2] * para_t,2), 0.5)

    return list(dis_line)

def Is_L2_trigger(data: pd.DataFrame, tw, distance, geo):
    trigger_dom = pd.DataFrame()
    for domid in np.unique(data['DomId']):
        dom_hit = data.loc[data['DomId']==domid]
        if Is_L1_trigger(dom_hit):
            dom_hit_first = dom_hit.sort_values(by='t0')
            trigger_dom = trigger_dom.append(dom_hit_first.iloc[0])
    if len(trigger_dom)<2:
        return False
    else:
        trigger_dom_time = trigger_dom.sort_values(by = ['t0'])
        
        diff_data_time = np.array(trigger_dom_time["t0"][coin_num-1:len(data)]) - np.array(trigger_dom_time["t0"][0:len(trigger_dom_time)-(coin_num-1)]) # The last term subtracts the first term
        
        diff_data_space = calculate_distance(geo, np.array(trigger_dom_time["DomId"][coin_num-1:len(data)]) , np.array(trigger_dom_time["DomId"][0:len(trigger_dom_time)-(coin_num-1)]))

        trigger_array_t = diff_data_time < tw
        trigger_array_space = np.abs(diff_data_space) < distance
        trigger_array = list(trigger_array_t) and list(trigger_array_space)
        if True in trigger_array:
            return True
        else:
            return False

#fake_data = pd.DataFrame({'t0':[],'PmtId':[],'DomId':[]})

def calculate_distance(geo, id_1, id_2):
    pos1 = geo.iloc[id_1]
    pos2 = geo.iloc[id_2]
    distance = np.power(pos1['x'].to_numpy() - pos2['x'].to_numpy(),2)
    distance += np.power(pos1['y'].to_numpy() - pos2['y'].to_numpy(),2)
    distance += np.power(pos1['z'].to_numpy() - pos2['z'].to_numpy(),2)
    distance = np.power(distance,0.5)
    return distance

def plot_trigger_eff(e_list, not_e_list, label_):
    counts_total,tick_ = np.histogram([*e_list,*not_e_list],range=range_, bins=bins_)
    counts_L2, tick_ = np.histogram(e_list,range=range_, bins=bins_)
    counts_not_L2, tick_ = np.histogram(not_e_list,range=range_, bins=bins_)
    energy_ = (tick_[:-1] + tick_[1:])/2
    error_L2 = (np.power(counts_L2, 0.5) / counts_total)**2
    error_L2 += (np.divide(counts_L2,np.power(counts_total, 2)) * np.power(counts_total, 0.5))**2
    error_L2 = np.power(error_L2,0.5)
    plt.figure(dpi = 800)
    plt.errorbar(x = energy_, y = np.divide(counts_L2,counts_total), yerr=error_L2,label = label_, fmt='.')
    plt.semilogx()
    plt.xlabel('neutrino energy [GeV]')
    # plt.xlabel('CosTheta')
    plt.ylim([0,1])
    plt.legend()
    plt.ylabel('Efficiency')
    return True


file_cqc_root_1 = np.load('/lustre/neutrino/weizhenyu/analysis_muon/data_npy/cqc_root_penrose.npy')
file_cqc_csv_1 = np.load('/lustre/neutrino/weizhenyu/analysis_muon/data_npy/cqc_csv_penrose.npy')
file_cqc_json_1 = np.load('/lustre/neutrino/weizhenyu/analysis_muon/data_npy/cqc_json_penrose.npy')
def obtain_all_energy(file_id):
    file_json = file_cqc_json_1[file_id]
    file_root = file_cqc_root_1[file_id]
    file = uproot.open(file_root)
    tree = file[file.keys()[2]]

    hits =tree.arrays(library = 'pd')
    with open(file_json) as f:
        mc_events = json.load(f) 
    energy_list = []
    for eventid in np.array(hits.index.levels[0]):
        eventinfo = mc_events[eventid]
        e = (eventinfo['particles_in'][0]['px']**2 + eventinfo['particles_in'][0]['py']**2 + eventinfo['particles_in'][0]['pz']**2)**0.5
        energy_list.append(e)
    return energy_list
def obtain_l2_trigger_result(file_id):
    file_root = file_cqc_root_1[file_id]
    file_json = file_cqc_json_1[file_id]
    file_csv = file_cqc_csv_1[file_id]

    with open(file_json) as f:
        mc_events = json.load(f) 
    geo_ = pd.read_csv(file_csv,names=['x','y','z'])
    file = uproot.open(file_root)
    tree = file[file.keys()[2]]

    hits =tree.arrays(library = 'pd')
    l2_e_list = []
    # not_l2_e_list = []
    for eventid in np.array(hits.index.levels[0]):
        event = hits.loc[eventid]
        event = qe_module(event)
        eventinfo = mc_events[eventid]
        e = (eventinfo['particles_in'][0]['px']**2 + eventinfo['particles_in'][0]['py']**2 + eventinfo['particles_in'][0]['pz']**2)**0.5
        
        if len(event) > 3:
            # not_l2_e_list.append(e)
        # else:
            if Is_L2_trigger(event, tw = 660, distance= 120, geo=geo_):
                l2_e_list.append(e)
            # else:
                # not_l2_e_list.append(e)
    return l2_e_list
#%%
if __name__ == '__main__':
    with concurrent.futures.ThreadPoolExecutor() as executor:
        l2_energy_list = executor.map(obtain_l2_trigger_result, range(len(file_cqc_root_1)))
        # total_energy_list = executor.map(obtain_all_energy, range(len(file_cqc_root_1)))

    l2_energy_list = [item for sublist in l2_energy_list for item in sublist]
    # total_energy_list = [item for sublist in total_energy_list for item in sublist]
    range_ = [1000,100000]
    bins_ = np.logspace(np.log10(range_[0]),np.log10(range_[1]),num=20)
    np.save('/lustre/neutrino/weizhenyu/analysis_trigger_L2/data_npy/l2_energy_list_cqc_penrose_2_10_120_660.npy', l2_energy_list)
    # np.save('/lustre/neutrino/weizhenyu/analysis_trigger_L2/data_npy/l2_energy_list_cqc_100_total.npy', total_energy_list)





#%%
# l2_e_list = []
# not_l2_e_list = []
# l2_e_list_tdc = []
# not_l2_e_list_tdc = []
# for id in range(len(file_cqc_root_1)):#tqdm(range(10)):
#     file_root = file_cqc_root_1[id]
#     file_json = file_cqc_json_1[id]
#     file_csv = file_cqc_csv_1[id]

#     with open(file_json) as f:
#         mc_events = json.load(f) 
#     geo = pd.read_csv(file_csv,names=['x','y','z'])
#     file = uproot.open(file_root)
#     tree = file[file.keys()[2]]

#     hits =tree.arrays(library = 'pd')
#     for eventid in np.array(hits.index.levels[0]):
#         event = hits.loc[eventid]
#         event = qe_module(event)
#         eventinfo = mc_events[eventid]
#         e = (eventinfo['particles_in'][0]['px']**2 + eventinfo['particles_in'][0]['py']**2 + eventinfo['particles_in'][0]['pz']**2)**0.5
#         dom_dis = geo.iloc[1]['z'] - geo.iloc[0]['z']  
#         if len(event) < 3:
#             not_l2_e_list.append(e)
#         else:
#             # if Is_L2_trigger(event, tw = dom_dis* 5 + 100):
#                 l2_e_list.append(e)
#             # else:
#                 # not_l2_e_list.append(e)

#         if len(event) < 2:
#             not_l2_e_list_tdc.append(e)
#         else:
#             l2_e_list_tdc.append(e)
        # else:
        #     l2_e_list.append(e)
# range_ = [1000,100000]
# bins_ = np.logspace(np.log10(range_[0]),np.log10(range_[1]),num=20)
# plot_trigger_eff(l2_e_list,not_l2_e_list, '4-photons cut')
# plot_trigger_eff(l2_e_list_tdc, not_l2_e_list_tdc, '3-photons cut')
# # %%
# plt.figure(dpi = 800)
# # range_ = [1000,100000]
# # bins_ = 100
# plt.hist(l2_e_list,range=range_, bins=bins_,histtype='step',label='triggered event')
# plt.hist(not_l2_e_list,range=range_, bins=bins_,histtype='step',label='not triggered event')
# plt.semilogx()
# plt.semilogy()
# plt.legend()
# plt.xlabel('neutrino energy [GeV]')
# # plt.xlabel('CosTheta')
# plt.ylabel('count')
# # %%
# counts_total,tick_ = np.histogram([*l2_e_list,*not_l2_e_list],range=range_, bins=bins_)
# counts_L2, tick_ = np.histogram(l2_e_list,range=range_, bins=bins_)
# counts_not_L2, tick_ = np.histogram(not_l2_e_list,range=range_, bins=bins_)
# energy_ = (tick_[:-1] + tick_[1:])/2
# error_L2 = (np.power(counts_L2, 0.5) / counts_total)**2
# error_L2 += (np.divide(counts_L2,np.power(counts_total, 2)) * np.power(counts_total, 0.5))**2
# error_L2 = np.power(error_L2,0.5)
# plt.errorbar(x = energy_, y = np.divide(counts_L2,counts_total), yerr=error_L2,label = '4-photons cut', fmt='.')
# plt.semilogx()
# plt.xlabel('neutrino energy [GeV]')
# # plt.xlabel('CosTheta')
# plt.ylim([0,1])
# plt.legend()
# plt.ylabel('Efficiency')
# %%
