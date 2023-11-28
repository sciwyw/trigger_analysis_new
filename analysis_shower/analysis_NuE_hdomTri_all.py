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

def calculate_distance_5m(hits_event:pd.DataFrame, eventinfo: dict):
    pos_vertex = [eventinfo['particles_in'][0]['x'],eventinfo['particles_in'][0]['y'],eventinfo['particles_in'][0]['z']]
    dir_vertex = [eventinfo['particles_in'][0]['px'],eventinfo['particles_in'][0]['py'],eventinfo['particles_in'][0]['pz']]
    norm_dir = (dir_vertex[0]**2 + dir_vertex[1]**2 + dir_vertex[2]**2)**0.5

    dx_5m = geo.iloc[hits_event['DomId']]['x'] - (pos_vertex[0] + dir_vertex[0] / norm_dir * 5)
    dy_5m = geo.iloc[hits_event['DomId']]['y'] - (pos_vertex[1] + dir_vertex[1] / norm_dir * 5)
    dz_5m = geo.iloc[hits_event['DomId']]['z'] - (pos_vertex[2] + dir_vertex[2] / norm_dir * 5)

    dis_list = np.power(np.power(dx_5m,2) + np.power(dy_5m,2) + np.power(dz_5m,2), 0.5)
    dis_list = dis_list.to_numpy()
    dom_id = hits_event['DomId'].to_numpy()
    dis_result = pd.DataFrame({'DomId':dom_id,'dis': dis_list})
    return dis_result

def calculate_distance_5m_geo_(geo_ : pd.DataFrame, eventinfo: dict):
    pos_vertex = [eventinfo['particles_in'][0]['x'],eventinfo['particles_in'][0]['y'],eventinfo['particles_in'][0]['z']]
    dir_vertex = [eventinfo['particles_in'][0]['px'],eventinfo['particles_in'][0]['py'],eventinfo['particles_in'][0]['pz']]
    norm_dir = (dir_vertex[0]**2 + dir_vertex[1]**2 + dir_vertex[2]**2)**0.5

    dx_5m = geo_['x'] - (pos_vertex[0] + dir_vertex[0] / norm_dir * 5)
    dy_5m = geo_['y'] - (pos_vertex[1] + dir_vertex[1] / norm_dir * 5)
    dz_5m = geo_['z'] - (pos_vertex[2] + dir_vertex[2] / norm_dir * 5)

    dis_list_geo = np.power(np.power(dx_5m,2) + np.power(dy_5m,2) + np.power(dz_5m,2), 0.5)
    return dis_list_geo
#%%
file_name_1TeV = '/lustre/neutrino/weizhenyu/analysis_shower/data_NuE_1TeV/data.root'
file_name_10TeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/jobs_10tev_nue_10k_p1/job_0/data.root'
file_name_100TeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/sy_data/nue_100tev_params/10k_p13/job_0/data/data.root'
file_name_1PeV = '/lustre/neutrino/weizhenyu/analysis_shower/data_NuE_1PeV/data.root'
file_name_10PeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/sy_data/nue_10pev_gv/10k_p0/job_0/data/data.root'

file_json_1TeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/jobs_1tev_nue/mc_events.json'
file_json_10TeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/jobs_10tev_nue_10k_p1/job_0/mc_events.json'
file_json_100TeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/sy_data/nue_100tev_params/10k_p13/job_0/mc_events.json'
file_json_1PeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/jobs_1pev_nue/mc_events.json'
file_json_10PeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/sy_data/nue_10pev_gv/10k_p0/job_0/mc_events.json'

file_geo_1TeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/jobs_1tev_nue/penrose.csv'
file_geo_10TeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/jobs_10tev_nue_10k_p1/penrose.csv'
file_geo_100TeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/sy_data/nue_100tev_params/10k_p13/job_0/penrose.csv' 
file_geo_1PeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/jobs_1pev_nue/penrose.csv'
file_geo_10PeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/sy_data/nue_10pev_gv/10k_p0/job_0/penrose.csv'

file_name_root_list = [file_name_1TeV,file_name_10TeV,file_name_100TeV,file_name_1PeV,file_name_10PeV]
file_name_json_list = [file_json_1TeV,file_json_10TeV,file_json_100TeV,file_json_1PeV,file_json_10PeV]
file_name_geo_list = [file_geo_1TeV,file_geo_10TeV,file_geo_100TeV,file_geo_1PeV,file_geo_10PeV]
#%%

# for e_id in range(3):
#     # l1_e_list = []
#     # not_l1_e_list = []
#     l1_dis_list = []
#     total_dis_list = []
#     root_file = file_name_root_list[e_id]
#     json_file = file_name_json_list[e_id] 
#     geo_file = file_name_geo_list[e_id]
#     with open(json_file) as f:
#             mc_events = json.load(f) 

#     file = uproot.open(root_file)
#     tree = file[file.keys()[2]]
#     hits =tree.arrays(library = 'pd')

#     geo = pd.read_csv(geo_file,names=['x','y','z'])
#     for event_id in range(len(np.unique(hits.index.get_level_values(0)))):#range(1000):
#         hit_event = hits.loc[event_id]
#         hit_event = qe_module(hit_event)
#         hit_event.sort_values(by=["t0"],inplace = True)
#         hit_event.reset_index(drop = True, inplace = True)

#         dis_5m = calculate_distance_5m(hits_event=hit_event, eventinfo=mc_events[event_id])
#         dis_5m_geo = calculate_distance_5m_geo_(geo_=geo, eventinfo= mc_events[event_id])
#         total_dis_list.extend(dis_5m_geo[dis_5m_geo < 250])
#         for domid in np.unique(hit_event['DomId']):
#             dom_hit = hit_event.loc[ hit_event['DomId'] == domid]
#             if len(dom_hit)>1:
#                 L1_bool = Is_L1_trigger(dom_hit)
#                 if L1_bool:
#                 l1_dis_list.append(dis_5m.loc[ dis_5m['DomId'] == dom_hit.iloc[0]['DomId'] ].iloc[0]['dis'])
#     np.save('/lustre/neutrino/weizhenyu/analysis_shower/data_npy/l1_dis_list_{}.npy'.format(e_id), l1_dis_list)
#     np.save('/lustre/neutrino/weizhenyu/analysis_shower/data_npy/total_dis_list_{}.npy'.format(e_id), total_dis_list)

# %%
# 1PeV
# l1_dis_list = []
# total_dis_list = []
# for i in range(100):
#     file_path_root_1PeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/sy_data/bl_1pev_gv/10k_p0/job_{}/data.root'.format(i)
#     file_path_json_1PeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/sy_data/bl_1pev_gv/10k_p0/job_{}/mc_events.json'.format(i)
#     file_path_geo_1PeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/sy_data/bl_1pev_gv/10k_p0/penrose.csv'
#     # root_file = file_name_root_list[e_id]
#     # json_file = file_name_json_list[e_id] 
#     # geo_file = file_name_geo_list[e_id]
#     try:
#         with open(file_path_json_1PeV) as f:
#             mc_events = json.load(f) 
    
#         file = uproot.open(file_path_root_1PeV)
#         tree = file[file.keys()[2]]
#         hits =tree.arrays(library = 'pd')

#         geo = pd.read_csv(file_path_geo_1PeV,names=['x','y','z'])
#         for event_id in range(len(np.unique(hits.index.get_level_values(0)))):#range(1000):
#             hit_event = hits.loc[event_id]
#             hit_event = qe_module(hit_event)
#             hit_event.sort_values(by=["t0"],inplace = True)
#             hit_event.reset_index(drop = True, inplace = True)

#             dis_5m = calculate_distance_5m(hits_event=hit_event, eventinfo=mc_events[event_id])
#             dis_5m_geo = calculate_distance_5m_geo_(geo_=geo, eventinfo= mc_events[event_id])
#             total_dis_list.extend(dis_5m_geo[dis_5m_geo < 250])
#             for domid in np.unique(hit_event['DomId']):
#                 dom_hit = hit_event.loc[ hit_event['DomId'] == domid]
#                 if len(dom_hit)>1:
#                     L1_bool = Is_L1_trigger(dom_hit)
#                     if L1_bool:
#                         l1_dis_list.append(dis_5m.loc[ dis_5m['DomId'] == dom_hit.iloc[0]['DomId'] ].iloc[0]['dis'])
#     except:
#         continue
# np.save('/lustre/neutrino/weizhenyu/analysis_shower/data_npy/l1_dis_list_{}.npy'.format(3), l1_dis_list)
# np.save('/lustre/neutrino/weizhenyu/analysis_shower/data_npy/total_dis_list_{}.npy'.format(3), total_dis_list)

#%%
# 10PeV
l1_dis_list = []
total_dis_list = []
job_10PeV_id_list = [0,2,4,5,6]#,8,9,11,12,13,14,15,16,17,18,20]
for i in job_10PeV_id_list:
    file_path_root_10PeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/sy_data/nue_10pev_gv/10k_p0/job_{}/data/data.root'.format(i)
    file_path_json_10PeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/sy_data/nue_10pev_gv/10k_p0/job_{}/mc_events.json'.format(i)
    file_path_geo_10PeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/sy_data/nue_10pev_gv/10k_p0/job_{}/penrose.csv'.format(i)
    # root_file = file_name_root_list[e_id]
    # json_file = file_name_json_list[e_id] 
    # geo_file = file_name_geo_list[e_id]
    try:
        with open(file_path_json_10PeV) as f:
                mc_events = json.load(f) 

        file = uproot.open(file_path_root_10PeV)
        tree = file[file.keys()[2]]
        hits =tree.arrays(library = 'pd')

        geo = pd.read_csv(file_path_geo_10PeV,names=['x','y','z'])
        for event_id in range(len(np.unique(hits.index.get_level_values(0)))):#range(1000):
            hit_event = hits.loc[event_id]
            hit_event = qe_module(hit_event)
            hit_event.sort_values(by=["t0"],inplace = True)
            hit_event.reset_index(drop = True, inplace = True)

            dis_5m = calculate_distance_5m(hits_event=hit_event, eventinfo=mc_events[event_id])
            dis_5m_geo = calculate_distance_5m_geo_(geo_=geo, eventinfo= mc_events[event_id])
            total_dis_list.extend(dis_5m_geo[dis_5m_geo < 250])
            for domid in np.unique(hit_event['DomId']):
                dom_hit = hit_event.loc[ hit_event['DomId'] == domid]
                if len(dom_hit)>1:
                    L1_bool = Is_L1_trigger(dom_hit)
                    if L1_bool:
                        l1_dis_list.append(dis_5m.loc[ dis_5m['DomId'] == dom_hit.iloc[0]['DomId'] ].iloc[0]['dis'])
    except:
        continue
np.save('/lustre/neutrino/weizhenyu/analysis_shower/data_npy/l1_dis_list_{}.npy'.format(4), l1_dis_list)
np.save('/lustre/neutrino/weizhenyu/analysis_shower/data_npy/total_dis_list_{}.npy'.format(4), total_dis_list)