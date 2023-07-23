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

# %%
file_name_muon = '/lustre/collider/mocen/project/hailing/data/muons/dataStore/penrose/detectorResponse/data.root'
file_json_muon = '/lustre/collider/mocen/project/hailing/data/muons/dataStore/penrose/detectorResponse/mc_events.json'
file_geo_muon = '/lustre/collider/mocen/project/hailing/data/muons/dataStore/penrose/detectorResponse/penrose.csv'
file_muon_cqc_1 = '/lustre/neutrino/changqichao/data/100/'
file_muon_cqc_2 = '/lustre/neutrino/changqichao/data/10/'
file_list_muon_cqc_1 = os.listdir(file_muon_cqc_1)
file_list_muon_cqc_2 = os.listdir(file_muon_cqc_2)
file_muon_cqc_penrose = '/lustre/neutrino/changqichao/icrc2023/telescope/penrose/1.0/'
# os.chdir(file_muon_cqc_1)
# file_dir = os.getcwd()
# file_name_list = []
# for root, dirs,files in os.walk(os.getcwd()):
#     for file in files:
#         if os.path.splitext(file)[-1] == '.root':
#             file_name_list.append(file)
path_cqc_root_list_2 = []
path_cqc_json_list_2 = []
path_cqc_geo_list_2 = []
os.chdir(file_muon_cqc_2)
file_dir = os.getcwd()
path_cqc_root_list_2 = []
for root, dirs,files in os.walk(os.getcwd()):
    for file in files:
        if os.path.splitext(file)[-1] == '.root':
            path_cqc_root_list_2.append(root + '/'+file)


for root_file in path_cqc_root_list_2:
    while not check_csv_in_list(os.listdir(os.path.split(root_file)[0])):
        root_file = os.path.split(root_file)[0]
    for string in os.listdir(os.path.split(root_file)[0]):
        if fnmatch.fnmatch(string, '*.csv'):
            path_cqc_geo_list_2.append(os.path.split(root_file)[0] + '/'+string)
for root_file in path_cqc_root_list_2:
    while not check_json_in_list(os.listdir(os.path.split(root_file)[0])):
        root_file = os.path.split(root_file)[0]
    for string in os.listdir(os.path.split(root_file)[0]):
        if fnmatch.fnmatch(string, '*.json'):
            path_cqc_json_list_2.append(os.path.split(root_file)[0] + '/'+string)
    

path_cqc_root_list_1 = []
path_cqc_json_list_1 = []
path_cqc_geo_list_1 = []
for file_cqc in file_list_muon_cqc_1:
    path_data = file_muon_cqc_1 + file_cqc + '/'
    # file_job = os.listdir(path_data)
    for jobid in range(100):
        path_root = path_data + 'job_{}/'.format(jobid)+'data/data.root'
        path_cqc_root_list_1.append(path_root)

        path_json = path_data + 'job_{}/'.format(jobid)+'mc_events.json'
        path_cqc_json_list_1.append(path_json)

        path_geo = path_data + 'job_{}/'.format(jobid)+'geometry.csv'
        path_cqc_geo_list_1.append(path_geo)

np.save('cqc_root_10',np.array(path_cqc_root_list_2))
np.save('cqc_json_10',np.array(path_cqc_json_list_2))
np.save('cqc_csv_10',np.array(path_cqc_geo_list_2))
np.save('/lustre/neutrino/weizhenyu/analysis_muon/cqc_root_100.npy',np.array(path_cqc_root_list_1))
np.save('/lustre/neutrino/weizhenyu/analysis_muon/cqc_json_100.npy',np.array(path_cqc_json_list_1))
np.save('/lustre/neutrino/weizhenyu/analysis_muon/cqc_csv_100.npy',np.array(path_cqc_geo_list_1))