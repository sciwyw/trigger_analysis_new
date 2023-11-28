import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import uproot
import yaml
import os
import glob
import json
from tqdm import tqdm
import scipy.interpolate
from matplotlib import rcParams
plt.clf()
from matplotlib import font_manager
for i in ['', 'bd', 'bi', 'i']:
    font_manager.fontManager.addfont('/home/zhiwei/fonts/times{}.ttf'.format(i))

plt.rcParams.update({
    'figure.figsize': (8, 5.5),
    'font.size': 20,
    'font.family': 'Times New Roman'
})

def getRealTime_Orb(N_event: int, radius: int):
    '''
    get getRealTime from jingping's result
    '''
    N_event = np.random.poisson(lam = N_event)
    activtyK40 = 10.87
    effective = 1
    volume = 4 * np.pi * radius**3 / 3
    mass = 1.04 * volume * 10**3
    frequncy = mass * activtyK40 * effective
    time = N_event / frequncy * 10 ** 9
    return time
def qe_module(hit: pd.DataFrame):
    #quantum efficiency
    # 1eV ~ 1240nm
    qe_list = pd.read_csv('qe.csv')
    qe_func = scipy.interpolate.interp1d(qe_list['wl'],qe_list['qe'],kind='linear')
    wl_hit = 1240 / np.array(hit["e0"])
    qe_photon = qe_func(wl_hit) / 100
    sample = np.random.sample(len(hit))
    survival_bool = qe_photon > sample
    hit_qe = hit[survival_bool]
    #print(survival_bool)
    hit_qe = hit_qe.reset_index(drop = True)

    return hit_qe
def read_hits(file_name: str) -> pd.DataFrame:
    file = uproot.open(file_name)
    tree = file[file.keys()[0]]
    hits =tree.arrays(library = 'pd')
    for i in range(len(hits["eventID"])):
        # if event id is zero, then the event id is the same with the last one
        if hits["eventID"][i] == 0:
            data_copy_1 = hits.iloc[i - 1].copy(deep=True)
            hits.loc[i, "eventID"] = data_copy_1["eventID"]
            hits.loc[i, "x(mm)"] = data_copy_1["x(mm)"]
            hits.loc[i, "y(mm)"] = data_copy_1["y(mm)"]
            hits.loc[i, "z(mm)"] = data_copy_1["z(mm)"]
    
    hits = qe_module(hits)
    hits["r"] = np.sqrt(
        np.power(hits["x(mm)"], 2) +
        np.power(hits["y(mm)"], 2) +
        np.power(hits["z(mm)"], 2)) / 1000
    return hits
# %%
# onefold
len_root = 0
r_total_list = []
for i in tqdm(range(0,500,1)):#range(500):
    #path = '/lustre/neutrino/weizhenyu/trigger/hdom_k40_new_optical/build_lib/run/{}/'.format(i)
    path = '/lustre/neutrino/weizhenyu/trigger/hdom_k40_new_optical/build_200/run/{}/'.format(i)
    hits = read_hits(get_root_files(path)[0])
    r_total_list = np.append(r_total_list,hits['r'])
    len_root += len(hits)
sim_time = getRealTime_Orb(20000000,200) * 500 / 10 ** 9
rate = len_root / sim_time / 31
print(rate)
# %%
# twofold
result_path = '/lustre/neutrino/weizhenyu/trigger/analysis/result/trigger_hdom_10_2_tts.csv'
result = pd.read_csv(result_path, index_col=0,skiprows=0,names=['event_id', 't', 'copynum', 'energy', 'x', 'y', 'z','batch','r'])
number_list = []
num = 1
for idx in range(len(result)-1):
    
    if result.iloc[idx]["batch"] == result.iloc[idx+1]["batch"]:
        if result.iloc[idx]["copynum"] != result.iloc[idx+1]["copynum"]:
            # if result.iloc[idx]["event_id"] != result.iloc[idx+1]["event_id"]:
                num+=1
    else:
        number_list.append(num)
        num=1
number_list = np.array(number_list)[np.array(number_list)!=1]
print(len(number_list) / 3.91)
# %%
# three fold
result_path = '/lustre/neutrino/weizhenyu/trigger/analysis/result/trigger_hdom_10_2_tts.csv'
result = pd.read_csv(result_path, index_col=0,skiprows=0,names=['event_id', 't', 'copynum', 'energy', 'x', 'y', 'z','batch','r'])
number_list = []
num = 1
for idx in range(len(result)-2):
    
    if result.iloc[idx]["batch"] == result.iloc[idx+2]["batch"]:
        if len(np.unique(result.iloc[idx:idx+3]["copynum"])) > 2:
            # if result.iloc[idx]["event_id"] != result.iloc[idx+1]["event_id"]:
                num+=1
    else:
        number_list.append(num)
        num=1
number_list = np.array(number_list)[np.array(number_list)!=1]
print(len(number_list)/ 3.91)
# %%
# four fold
result_path = '/lustre/neutrino/weizhenyu/trigger/analysis/result/trigger_hdom_10_2_tts.csv'
result = pd.read_csv(result_path, index_col=0,skiprows=0,names=['event_id', 't', 'copynum', 'energy', 'x', 'y', 'z','batch','r'])
number_list = []
num = 1
for idx in range(len(result)-3):
    
    if result.iloc[idx]["batch"] == result.iloc[idx+3]["batch"]:
        if len(np.unique(result.iloc[idx:idx+4]["copynum"])) > 3:
            # if result.iloc[idx]["event_id"] != result.iloc[idx+1]["event_id"]:
                num+=1
    else:
        number_list.append(num)
        num=1
number_list = np.array(number_list)[np.array(number_list)!=1]
print(len(number_list)/ 3.91)
# %%
# five fold
number_list = []
num = 1
for idx in range(len(result)-4):
    
    if result.iloc[idx]["batch"] == result.iloc[idx+4]["batch"]:
        if len(np.unique(result.iloc[idx:idx+5]["copynum"])) > 4:
            # if result.iloc[idx]["event_id"] != result.iloc[idx+1]["event_id"]:
                num+=1
    else:
        number_list.append(num)
        num=1
number_list = np.array(number_list)[np.array(number_list)!=1]
print(len(number_list)/ 3.91)
# %%
# six fold
number_list = []
num = 1
for idx in range(len(result)-5):
    
    if result.iloc[idx]["batch"] == result.iloc[idx+5]["batch"]:
        if len(np.unique(result.iloc[idx:idx+6]["copynum"])) > 5:
            # if result.iloc[idx]["event_id"] != result.iloc[idx+1]["event_id"]:
                num+=1
    else:
        number_list.append(num)
        num=1
number_list = np.array(number_list)[np.array(number_list)!=1]
print(len(number_list)/ 3.91)
# %%
# seven fold
number_list = []
num = 1
for idx in range(len(result)-6):
    
    if result.iloc[idx]["batch"] == result.iloc[idx+6]["batch"]:
        if len(np.unique(result.iloc[idx:idx+7]["copynum"])) > 6:
            # if result.iloc[idx]["event_id"] != result.iloc[idx+1]["event_id"]:
                num+=1
    else:
        number_list.append(num)
        num=1
number_list = np.array(number_list)[np.array(number_list)!=1]
print(len(number_list)/ 3.91)
# %%
# eight fold
number_list = []
num = 1
for idx in range(len(result)-7):
    
    if result.iloc[idx]["batch"] == result.iloc[idx+7]["batch"]:
        if len(np.unique(result.iloc[idx:idx+8]["copynum"])) > 7:
            # if result.iloc[idx]["event_id"] != result.iloc[idx+1]["event_id"]:
                num+=1
    else:
        number_list.append(num)
        num=1
number_list = np.array(number_list)[np.array(number_list)!=1]
print(len(number_list)/ 3.91)
# %%
# night fold
number_list = []
num = 1
for idx in range(len(result)-8):
    
    if result.iloc[idx]["batch"] == result.iloc[idx+8]["batch"]:
        if len(np.unique(result.iloc[idx:idx+9]["copynum"])) > 8:
            # if result.iloc[idx]["event_id"] != result.iloc[idx+1]["event_id"]:
                num+=1
    else:
        number_list.append(num)
        num=1
number_list = np.array(number_list)[np.array(number_list)!=1]
print(len(number_list)/ 3.91)
