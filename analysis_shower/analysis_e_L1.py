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
coin_num = 3
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


# %%

file_name_1TeV = '/lustre/neutrino/weizhenyu/analysis_shower/data_NuE_1TeV/data.root'
file_name_10TeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/jobs_10tev_nue_10k_p1/job_0/data.root'
file_name_100TeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/sy_data/nue_100tev_params/10k_p13/job_0/data/data.root'
file_name_1PeV = '/lustre/neutrino/weizhenyu/analysis_shower/data_NuE_1PeV/data.root'
file_name_10PeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/jobs_10pev_nue/data/data.root'

file_json_1TeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/jobs_1tev_nue/mc_events.json'
file_json_10TeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/jobs_10tev_nue_10k_p1/job_0/mc_events.json'
file_json_100TeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/sy_data/nue_100tev_params/10k_p13/job_0/mc_events.json'
file_json_1PeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/jobs_1pev_nue/mc_events.json'
file_json_10PeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/jobs_10pev_nue/mc_events.json'

file_geo_1TeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/jobs_1tev_nue/penrose.csv'
file_geo_10TeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/jobs_10tev_nue_10k_p1/penrose.csv'
file_geo_100TeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/sy_data/nue_100tev_params/10k_p13/job_0/penrose.csv' 
file_geo_1PeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/jobs_1pev_nue/penrose.csv'
file_geo_10PeV = ''


# %%
with open(file_json_1TeV) as f:
    mc_events = json.load(f) 

file = uproot.open(file_name_1TeV)
tree = file[file.keys()[2]]
hits =tree.arrays(library = 'pd')

geo = pd.read_csv(file_geo_1TeV,names=['x','y','z'])
# %%
f = open("output.txt", "a")
not_L1_dis = []
L1_dis = []
dis_list = []

not_L1_dis_5m = []
L1_dis_5m = []
dis_list_5m = []

not_L1_dis_line = []
L1_dis_line = []
dis_list_line = []

not_L1_dis_cos = []
L1_dis_cos = []
dis_list_cos = []

not_L1_id = []
photon_number_list = []
pmt_number_list = []
#%%
for event_id in range(100):#range(1000):
    hit_event = hits.loc[event_id]
    hit_event = qe_module(hit_event)
    hit_event.sort_values(by=["t0"],inplace = True)
    hit_event.reset_index(drop = True, inplace = True)

    pos_vertex = [mc_events[event_id]['particles_in'][0]['x'],mc_events[event_id]['particles_in'][0]['y'],mc_events[event_id]['particles_in'][0]['z']]
    dir_vertex = [mc_events[event_id]['particles_in'][0]['px'],mc_events[event_id]['particles_in'][0]['py'],mc_events[event_id]['particles_in'][0]['pz']]
    norm_dir = (dir_vertex[0]**2 + dir_vertex[1]**2 + dir_vertex[2]**2)**0.5

    dx = geo.iloc[hit_event['DomId']]['x'] - pos_vertex[0]
    dy = geo.iloc[hit_event['DomId']]['y'] - pos_vertex[1]
    dz = geo.iloc[hit_event['DomId']]['z'] - pos_vertex[2]
    dx_5m = geo.iloc[hit_event['DomId']]['x'] - (pos_vertex[0] + dir_vertex[0] / norm_dir * 5)
    dy_5m = geo.iloc[hit_event['DomId']]['y'] - (pos_vertex[1] + dir_vertex[1] / norm_dir * 5)
    dz_5m = geo.iloc[hit_event['DomId']]['z'] - (pos_vertex[2] + dir_vertex[2] / norm_dir * 5)

    dis = np.power(np.power(dx,2) + np.power(dy,2) + np.power(dz,2), 0.5)
    dis_5m = np.power(np.power(dx_5m,2) + np.power(dy_5m,2) + np.power(dz_5m,2), 0.5)
    
    para_t = (dir_vertex[0] * dx + dir_vertex[1] * dy + dir_vertex[2] * dz) / norm_dir ** 2
    dis_line = np.power(np.power(dx - dir_vertex[0] * para_t,2) + np.power(dy - dir_vertex[1] * para_t,2) + np.power(dz - dir_vertex[2] * para_t,2), 0.5)

    dis_cos = np.multiply(dir_vertex[0] / norm_dir, np.divide(dx , dis)) + np.multiply(dir_vertex[1] / norm_dir, np.divide(dy , dis)) + np.multiply(dir_vertex[2] / norm_dir, np.divide(dz , dis))

    for domid in np.unique(hit_event['DomId']):
        dom_hit = hit_event.loc[ hit_event['DomId'] == domid]
        dom_photon_number = len(dom_hit)
        dom_pmt_number = len(np.unique(dom_hit['PmtId']))
        photon_number_list.append(dom_photon_number)
        pmt_number_list.append(dom_pmt_number)

        dom_dis = dis.iloc[dom_hit.index[0]]
        dom_dis_5m = dis_5m.iloc[dom_hit.index[0]]
        dom_dis_line = dis_line.iloc[dom_hit.index[0]]
        dom_dis_cos = dis_cos.iloc[dom_hit.index[0]]

        dis_list.append(dom_dis)
        dis_list_5m.append(dom_dis_5m)
        dis_list_line.append(dom_dis_line)
        dis_list_cos.append(dom_dis_cos)
        if len(dom_hit)<2:
            not_L1_dis.append(dom_dis)
            not_L1_dis_5m.append(dom_dis_5m)
            not_L1_dis_line.append(dom_dis_line)
            not_L1_dis_cos.append(dom_dis_cos)

            not_l1 = [event_id, domid]
            not_L1_id.append(not_l1)
            if dom_dis_5m<10:
                if dom_dis_cos > 0.7:
                    dom_hit.to_csv('strange_example.csv', header=True, index=False, sep='\t', mode='a')
                    print(mc_events[event_id]['particles_in'] , file=f)
                    print('\n', file=f)
                    print(geo.iloc[dom_hit['DomId']]['x'] , file=f)
                    print('\n', file=f)
                    print(geo.iloc[dom_hit['DomId']]['y'] , file=f)
                    print('\n', file=f)
                    print(geo.iloc[dom_hit['DomId']]['z'] , file=f)
                    print('\n', file=f)
                    print(dom_dis_5m,file = f)
                    print('\n', file=f)
                    print(event_id,file = f)
                    print('\n', file=f)
                    

        else:
            if Is_L1_trigger(dom_hit):
                L1_dis.append(dom_dis)
                L1_dis_5m.append(dom_dis_5m)
                L1_dis_line.append(dom_dis_line)
                L1_dis_cos.append(dom_dis_cos)
                
            else:
                not_L1_dis.append(dom_dis)
                not_L1_dis_5m.append(dom_dis_5m)
                not_L1_dis_line.append(dom_dis_line)
                not_L1_dis_cos.append(dom_dis_cos)

                not_l1 = [event_id, domid]
                not_L1_id.append(not_l1)
                if dom_dis_5m<10:
                    if dom_dis_cos > 0.7:
                        dom_hit.to_csv('strange_example.csv', header=True, index=False, sep='\t', mode='a')
                        print(mc_events[event_id]['particles_in'] , file=f)
                        print('\n', file=f)
                        print(geo.iloc[dom_hit['DomId']]['x'] , file=f)
                        print('\n', file=f)
                        print(geo.iloc[dom_hit['DomId']]['y'] , file=f)
                        print('\n', file=f)
                        print(geo.iloc[dom_hit['DomId']]['z'] , file=f)
                        print('\n', file=f)
                        print(dom_dis_5m,file = f)
                        print('\n', file=f)
                        print(event_id,file = f)
                        print('\n', file=f)

f.close()
#%%
plt.figure(dpi = 800)
range_ = [5,300]
bins_ = 50
plt.hist(dis_list_5m,range=range_, bins=bins_,histtype='step',label='all hDOMs')
plt.hist(L1_dis_5m,range=range_, bins=bins_,histtype='step',label='L1 hDOM')
plt.hist(not_L1_dis_5m,range=range_, bins=bins_,histtype='step',label='not-L1 hDOM')
plt.legend()
plt.xlabel('Distance [m]')
# plt.xlabel('CosTheta')
plt.ylabel('number of hDOM')

    
# %%
counts_total,tick_ = np.histogram(dis_list_5m,range=range_, bins=bins_)
counts_L1, tick_ = np.histogram(L1_dis_5m,range=range_, bins=bins_)
counts_not_L1, tick_ = np.histogram(not_L1_dis_5m,range=range_, bins=bins_)
dis = (tick_[:-1] + tick_[1:])/2
error_L1 = (np.power(counts_L1, 0.5) / counts_total)**2
error_L1 += (np.divide(counts_L1,np.power(counts_total, 2)) * np.power(counts_total, 0.5))**2
error_L1 = np.power(error_L1,0.5)
plt.errorbar(x = dis, y = np.divide(counts_L1,counts_total), yerr=error_L1,label = 'L1 Ratio', fmt='.')
error_notL1 = (np.power(counts_not_L1, 0.5) / counts_total)**2
error_notL1 += (np.divide(counts_not_L1,np.power(counts_total, 2)) * np.power(counts_total, 0.5))**2
error_notL1 = np.power(error_notL1,0.5)
plt.errorbar(x = dis, y = np.divide(counts_not_L1,counts_total), yerr= error_notL1, label = 'not L1 Ratio', fmt='.')
plt.xlabel('Distance [m]')
# plt.xlabel('CosTheta')
plt.ylim([0,1])
plt.legend()
plt.ylabel('Ratio')
# %%
plt.figure(dpi=800)
plt.scatter(dis_list, photon_number_list,s = 0.1)
plt.scatter(dis_list_cos, pmt_number_list,s = 0.1)
# %%

plt.figure(dpi = 800)
sns.jointplot(x = dis_list_5m, y = photon_number_list, kind = 'kde',ylim=[0,800],xlim=[5,300],palette='Oranges',fill = True)
# plt.semilogy()
plt.xlabel('Distance [m]')
# plt.xlabel('CosTheta')
plt.ylabel('photon number')


# %%
plt.figure(dpi = 800)
sns.jointplot(x = dis_list_5m, y = pmt_number_list, kind = 'kde', xlim=[5,300],fill = True)
plt.xlabel('Distance [m]')
# plt.xlabel('CosTheta')
plt.ylabel('PMT number')

# %%
plt.figure(dpi = 800)
data_plot = pd.DataFrame({'label','Distance [m]','CosTheta'})
L1 = pd.DataFrame({'label':['L1'] * len(L1_dis_5m),'Distance [m]': L1_dis_5m,'CosTheta': L1_dis_cos})

not_L1 = pd.DataFrame({'label':['not L1'] * len(not_L1_dis_5m),'Distance [m]': not_L1_dis_5m,'CosTheta': not_L1_dis_cos})
data_plot = L1.append(not_L1)
sns.jointplot(data=data_plot, x = 'Distance [m]', y = 'CosTheta', kind = 'kde',xlim=[0,300],ylim=[-1,1], hue='label')
# sns.jointplot(x = not_L1_dis_5m, y = not_L1_dis_cos, kind = 'kde',xlim=[0,300],ylim=[-1,1],palette='Oranges',fill = True)
# plt.xlabel('Distance [m]')
# plt.ylabel('CosTheta')
# plt.xlabel('Distance [m]')

#%%
plt.figure(dpi = 800)
plt.hist(photon_number_list, histtype='step', density = True, range=[0,1000], bins=1000, label='photon number')
plt.axvline(x = np.quantile(photon_number_list, 0.5), linestyle = '--', color = 'r', label = 'quantile 50%: 3')
plt.axvline(x = np.quantile(photon_number_list, 0.95), linestyle = '--', color = 'g', label = 'quantile 95%: 202')
plt.legend()
plt.semilogy()
plt.semilogx()
plt.xlabel('Photon Number')
plt.ylabel('Density')
# %%
