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
plt.figure(dpi=800)
time_window_list = [10]
coin_num = 2
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
        for hits_id in range(len(trigger_array)-3):
            if data_trigger.iloc[hits_id+1]["t0"] - data_trigger.iloc[hits_id]["t0"] < time_window:
                if data_trigger.iloc[hits_id+1]["PmtId"] == data_trigger.iloc[hits_id]["PmtId"]:
                    idx_drop.append([data_trigger.index[hits_id],data_trigger.index[hits_id+1]])
        data_trigger.drop(idx_drop,inplace = True)
            
        t_trigger_start = data_trigger.iloc[0]["t0"] - 300
        t_trigger_end = data_trigger.iloc[-1]["t0"] + 700
        data_new = data.loc[data["t0"] > t_trigger_start]
        data_new = data_new.loc[data_new["t0"] < t_trigger_end]
        if len(data_trigger) != len(data_new):
            data_trigger = data_trigger.append(data_new)
        data_trigger.drop_duplicates(subset=None, keep= 'first', inplace= True)
    

    return data_trigger
# file_json = "/lustre/neutrino/zhangfuyudi/hailing/cascade/jobs_10tev_nue_10k_p1/job_0/mc_events.json"
# with open(file_json) as file:
#     mc_events = json.load(file) 
file_name = "/lustre/neutrino/zhangfuyudi/hailing/cascade/jobs_10tev_nue_10k_p1/job_0/data.root"
file = uproot.open(file_name)
tree = file[file.keys()[2]]
hits =tree.arrays(library = 'pd')
for time_window in time_window_list:
    dom_num = []
    trigger_dom_num = []
    for i in range(1000):
        events = hits.loc[i]
        events = qe_module(events)
        dom_num.append(len(np.unique(events["DomId"])))
        count_trigger_dom = 0
        
        for dom_id in np.unique(events["DomId"]):
            try:
                dom_hit = events.loc[events["DomId"]==int(dom_id)]
                dom_hit.sort_values(by=["t0"],inplace = True)
                dom_hit_trigger = obtain_trigger(dom_hit)
                if len(dom_hit_trigger) > 0:
                    count_trigger_dom += 1
            except:
                count_trigger_dom  += 0
                print(dom_id)
        trigger_dom_num.append(count_trigger_dom)
            
    plt.hist(np.divide(trigger_dom_num,dom_num),range=[0,1],bins=20,histtype="step",label="{}ns {}".format(time_window,np.mean(np.divide(trigger_dom_num,dom_num))))
plt.legend()
plt.xlabel("Trigger Ratio")
plt.ylabel("count")
plt.savefig("shower_trigger_10.png")

#%%
plt.figure(dpi=800)
plt.hist(dom_num,histtype="step",label='Hit hDOM',density=True,bins=10)
plt.hist(trigger_dom_num,histtype="step",label='L1-triggered hDOM',density=True,bins=10)
plt.legend()
plt.semilogy()
plt.xlabel("Number of DOM")
plt.ylabel("density")
# %%
