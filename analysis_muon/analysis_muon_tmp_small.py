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

def FD_func(x,*param):
    return 1 / (1 + np.exp((x-param[0])/param[1]))
# popt2,pcov2 = op.curve_fit(FD_func,dis[0:25], np.divide(counts_L1,counts_total)[0:25], p0=[50,2])
# error_list = np.sqrt(np.diag(pcov2))
# plt.plot(np.linspace(0,100,100), FD_func(np.linspace(0,100,100), *popt2))
# def calculate_eff_radius(array_l1_dis, array_not_l1_dis):


    return eff_radius, error_bar_upper, error_bar_lower
# %%
# file_name_muon = '/lustre/collider/mocen/project/hailing/data/muons/dataStore/penrose/detectorResponse/data.root'
# file_json_muon = '/lustre/collider/mocen/project/hailing/data/muons/dataStore/penrose/detectorResponse/mc_events.json'
# file_geo_muon = '/lustre/collider/mocen/project/hailing/data/muons/dataStore/penrose/detectorResponse/penrose.csv'
# file_muon_cqc_1 = '/lustre/neutrino/changqichao/data/100/'
# file_muon_cqc_2 = '/lustre/neutrino/changqichao/data/10/'
# file_list_muon_cqc_1 = os.listdir(file_muon_cqc_1)
# file_list_muon_cqc_2 = os.listdir(file_muon_cqc_2)

# # os.chdir(file_muon_cqc_1)
# # file_dir = os.getcwd()
# # file_name_list = []
# # for root, dirs,files in os.walk(os.getcwd()):
# #     for file in files:
# #         if os.path.splitext(file)[-1] == '.root':
# #             file_name_list.append(file)
# path_cqc_root_list_2 = []
# path_cqc_json_list_2 = []
# path_cqc_geo_list_2 = []
# os.chdir(file_muon_cqc_2)
# file_dir = os.getcwd()
# path_cqc_root_list_2 = []
# for root, dirs,files in os.walk(os.getcwd()):
#     for file in files:
#         if os.path.splitext(file)[-1] == '.root':
#             path_cqc_root_list_2.append(root + '/'+file)


# for root_file in path_cqc_root_list_2:
#     while not check_csv_in_list(os.listdir(os.path.split(root_file)[0])):
#         root_file = os.path.split(root_file)[0]
#     for string in os.listdir(os.path.split(root_file)[0]):
#         if fnmatch.fnmatch(string, '*.csv'):
#             path_cqc_geo_list_2.append(os.path.split(root_file)[0] + '/'+string)
# for root_file in path_cqc_root_list_2:
#     while not check_json_in_list(os.listdir(os.path.split(root_file)[0])):
#         root_file = os.path.split(root_file)[0]
#     for string in os.listdir(os.path.split(root_file)[0]):
#         if fnmatch.fnmatch(string, '*.json'):
#             path_cqc_json_list_2.append(os.path.split(root_file)[0] + '/'+string)
    

# path_cqc_root_list_1 = []
# path_cqc_json_list_1 = []
# path_cqc_geo_list_1 = []
# for file_cqc in file_list_muon_cqc_1:
#     path_data = file_muon_cqc_1 + file_cqc + '/'
#     # file_job = os.listdir(path_data)
#     for jobid in range(100):
#         path_root = path_data + 'job_{}/'.format(jobid)+'data/data.root'
#         path_cqc_root_list_1.append(path_root)

#         path_json = path_data + 'job_{}/'.format(jobid)+'mc_events.json'
#         path_cqc_json_list_1.append(path_json)

#         path_geo = path_data + 'job_{}/'.format(jobid)+'geometry.csv'
#         path_cqc_geo_list_1.append(path_geo)
file_cqc_root_1 = np.load('/lustre/neutrino/weizhenyu/analysis_muon/cqc_root_100.npy')
file_cqc_csv_1 = np.load('/lustre/neutrino/weizhenyu/analysis_muon/cqc_csv_100.npy')
file_cqc_json_1 = np.load('/lustre/neutrino/weizhenyu/analysis_muon/cqc_json_100.npy')
# %%

l1_e_list = []
not_l1_e_list = []
l1_dis_list = []
not_l1_dis_list = []
for id in tqdm(range(50)):
    file_root = file_cqc_root_1[id]
    file_json = file_cqc_json_1[id]
    file_csv = file_cqc_csv_1[id]

    with open(file_json) as f:
        mc_events = json.load(f) 
    geo = pd.read_csv(file_csv,names=['x','y','z'])
# energy unit is GeV
    file = uproot.open(file_root)
    tree = file[file.keys()[2]]

    hits =tree.arrays(library = 'pd')
    for eventid in np.array(hits.index.levels[0]):
        event = hits.loc[eventid]
        event = qe_module(event)
        eventinfo = mc_events[eventid]
        e = (eventinfo['particles_in'][0]['px']**2 + eventinfo['particles_in'][0]['py']**2 + eventinfo['particles_in'][0]['pz']**2)**0.5
        if len(event) > 0:
            for domid in np.array(event['DomId']):
                event_dom = event.loc[event['DomId']==domid]
                dis = calculate_distance(event_dom, eventinfo)
                if len(event_dom) < 2:
                    not_l1_dis_list.append(dis[0])
                    not_l1_e_list.append(e)
                else:
                    L1_bool = Is_L1_trigger(event_dom)
                    if L1_bool:
                        l1_dis_list.append(dis[0])
                        l1_e_list.append(e)
                    else:
                        not_l1_dis_list.append(dis[0])
                        not_l1_e_list.append(e)
e_total = l1_e_list + not_l1_e_list
l1_total = l1_dis_list + not_l1_dis_list
#%%
np.save('/lustre/neutrino/weizhenyu/analysis_muon/l1_e_list_tmp.npy', l1_e_list)
np.save('/lustre/neutrino/weizhenyu/analysis_muon/not_l1_e_list_tmp.npy', not_l1_e_list)
np.save('/lustre/neutrino/weizhenyu/analysis_muon/l1_dis_list_tmp.npy', l1_dis_list)
np.save('/lustre/neutrino/weizhenyu/analysis_muon/not_l1_dis_list_tmp.npy', not_l1_dis_list)
#%%         
# plt.figure(dpi = 800)
# range_ = [0,200]
# bins_ = 50
# plt.hist(l1_total,range=range_, bins=bins_,histtype='step',label='all hDOMs')
# plt.hist(l1_dis_list,range=range_, bins=bins_,histtype='step',label='L1 hDOM')
# plt.hist(not_l1_dis_list,range=range_, bins=bins_,histtype='step',label='not-L1 hDOM')
# plt.legend()
# plt.xlabel('Distance [m]')
# # plt.xlabel('CosTheta')
# plt.ylabel('number of hDOM')        
# #%%
# counts_total,tick_ = np.histogram(l1_total,range=range_, bins=bins_)
# counts_L1, tick_ = np.histogram(l1_dis_list,range=range_, bins=bins_)
# counts_not_L1, tick_ = np.histogram(not_l1_dis_list,range=range_, bins=bins_)
# dis = (tick_[:-1] + tick_[1:])/2
# error_L1 = (np.power(counts_L1, 0.5) / counts_total)**2
# error_L1 += (np.divide(counts_L1,np.power(counts_total, 2)) * np.power(counts_total, 0.5))**2
# error_L1 = np.power(error_L1,0.5)
# plt.errorbar(x = dis, y = np.divide(counts_L1,counts_total), yerr=error_L1,label = 'L1 Ratio', fmt='.')
# error_notL1 = (np.power(counts_not_L1, 0.5) / counts_total)**2
# error_notL1 += (np.divide(counts_not_L1,np.power(counts_total, 2)) * np.power(counts_total, 0.5))**2
# error_notL1 = np.power(error_notL1,0.5)
# plt.errorbar(x = dis, y = np.divide(counts_not_L1,counts_total), yerr= error_notL1, label = 'not L1 Ratio', fmt='.')
# plt.xlabel('Distance [m]')
# # plt.xlabel('CosTheta')
# plt.ylim([0,1])
# plt.legend()
# plt.ylabel('Ratio')

#%%
hdom_trigger_info = pd.DataFrame({'e_l1':l1_e_list,'dis_l1':l1_dis_list})
hdom_not_trigger_info = pd.DataFrame({'e_not_l1': not_l1_e_list, 'dis_not_l1': not_l1_dis_list})

# for e in np.arange():


