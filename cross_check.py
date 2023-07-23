#%%
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
def find_position(copynum:int):
    #input the copy number of a pmt
    #output the direction of the pmt
    theta_array = np.array([55.,73.,107.,125.,150.,180.]) * np.pi / 180
    phi_array_1 = np.array([60., 120., 180., 240., 300., 360.]) * np.pi / 180
    phi_array_2 = np.array([30., 90., 150., 210., 270., 330.]) * np.pi / 180
    phi_array_3 = np.array([60., 120., 180., 240., 300., 360.]) * np.pi / 180
    phi_array_4 = np.array([30., 90., 150., 210., 270., 330.]) * np.pi / 180
    phi_array_5 = np.array([60., 120., 180., 240., 300., 360.]) * np.pi / 180
    phi_array_6 = np.array([0.]) * np.pi / 180
    
    theta = theta_array[int(np.floor(copynum / 6))]
    if copynum < 6:
        phi = phi_array_1[copynum]
    if (copynum > 5) & (copynum < 12):
        phi = phi_array_2[copynum-6]
    if (copynum > 11) & (copynum < 18):
        phi = phi_array_3[copynum-12]
    if (copynum > 17) & (copynum < 24):
        phi = phi_array_4[copynum-18]
    if (copynum > 23) & (copynum < 30):
        phi = phi_array_5[copynum-24]
    if copynum == 30:
        phi = phi_array_6[copynum-30]

    position_vector = [np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)]
    return position_vector
def isNeighbour(copynum1:int, copynum2:int):
    # decide whether two pmts are neighbour
    costheta = np.sum(np.multiply(find_position(copynum1), find_position(copynum2)))
    bool_isNeighbour = (costheta > thres)
    if copynum1 == copynum2:
        bool_isNeighbour = True 
    return bool_isNeighbour
def get_root_files(directory):
    root_files = glob.glob(os.path.join(directory, '*.root'))
    return root_files
rate_coin_xian = 2.6 
rate_hit_xian = 2000
pmt_size_xian = 54.9**2 * (1 - np.cos(58/180*np.pi))
# G4Sphere *solidpmt = new G4Sphere("PMT", 53. * mm, 54.9 * mm, 0, 2 * M_PI, 0, 58. / 180 * M_PI);
pmt_size_wei = 53 ** 2 * (1 - np.cos(45/180*np.pi))
#solidPmt = new G4Sphere("PMT", 52 * mm, 53 * mm, 0, 2 * M_PI, 0, 45. / 180 * M_PI);
rate_coin_xian_to_wei = rate_coin_xian * pmt_size_wei / pmt_size_xian * 31 # per pmt in wei
rate_hit_xian_to_wei = rate_hit_xian * pmt_size_wei / pmt_size_xian 
#%%
sim_time_wei = getRealTime_Orb(10000000,30) * 500 / 10 ** 9
rate_coin_wei = len(pd.read_csv("/lustre/neutrino/weizhenyu/trigger/analysis/result/trigger_hdom.csv")) / sim_time_wei / 31
rate_coin_wei_12 = len(pd.read_csv("/lustre/neutrino/weizhenyu/trigger/analysis/result/trigger_hdom_12.csv")) / sim_time_wei 
#%%
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
#%%
rate_hit_wei = len_root / sim_time_wei / 31
print('rate of coincident hit per pmt(2.6Hz xian to wei)')
print(rate_coin_xian_to_wei)
print('rate of coincident hit per pmt(wei)')
print(rate_coin_wei)
print('rate of coincident hit on pmt1&2 (wei)')
print(rate_coin_wei_12/2)
print('\n')
print('rate of k40 hit (xian_to_wei)')
print(rate_hit_xian_to_wei)
print('rate of k40 hit (wei)')
print(rate_hit_wei)
#%%
result_path = "/lustre/neutrino/weizhenyu/trigger/analysis/result/trigger_hdom_10_2_tts_5.csv"
result = pd.read_csv(result_path, index_col=0,skiprows=0,names=['event_id', 't', 'copynum', 'energy', 'x', 'y', 'z','batch','r'])
result12 = pd.read_csv("/lustre/neutrino/weizhenyu/trigger/analysis/result/trigger_hdom_12.csv", index_col=0,skiprows=0,names=['event_id', 't', 'copynum', 'energy', 'x', 'y', 'z','batch','r'])
result_0_24 = pd.read_csv("/lustre/neutrino/weizhenyu/trigger/analysis/result/trigger_hdom_0_24.csv", index_col=0,skiprows=0,names=['event_id', 't', 'copynum', 'energy', 'x', 'y', 'z','batch','r'])
result_0_12 = pd.read_csv("/lustre/neutrino/weizhenyu/trigger/analysis/result/trigger_hdom_0_12.csv", index_col=0,skiprows=0,names=['event_id', 't', 'copynum', 'energy', 'x', 'y', 'z','batch','r'])
result12.reset_index(drop=True,inplace=True)
#%%
r_list  =  []
for idx in range(len(result)-1):
    
    if result.iloc[idx]["batch"] == result.iloc[idx+1]["batch"]:
        if result.iloc[idx]["copynum"] != result.iloc[idx+1]["copynum"]:
            # if isNeighbour(int(result_0_19.iloc[idx]["copynum"]), int(result_0_19.iloc[idx+1]["copynum"])):
            r_list.append(result.iloc[idx]["r"])
            r_list.append(result.iloc[idx+1]["r"])
r_list_diff  =  []
for idx in range(len(result)-1):
    
    if result.iloc[idx]["batch"] == result.iloc[idx+1]["batch"]:
        if result.iloc[idx]["copynum"] == result.iloc[idx+1]["copynum"]:
            # if isNeighbour(int(result_0_19.iloc[idx]["copynum"]), int(result_0_19.iloc[idx+1]["copynum"])):
            r_list_diff.append(result.iloc[idx]["r"])
            r_list_diff.append(result.iloc[idx+1]["r"])
plt.figure(dpi=800)
plt.hist(r_list,range=[1,30],bins=45,label='different PMT',histtype='step',density=True)
plt.hist(r_total_list,range=[1,30],bins=45,label='raw K-40 decay hit',histtype='step',density=True)
plt.hist(r_list_diff,range=[1,30],bins=45,label='same PMT',histtype='step',density=True)
plt.semilogy()
plt.legend()
plt.xlabel('Radius [m]')
plt.ylabel('density')
#%%
diff = []
cos_dis = []
for idx in range(len(result)-1):
    
    if result.iloc[idx]["batch"] == result.iloc[idx+1]["batch"]:
        diff.append(result.iloc[idx+1]["t"] - result.iloc[idx]["t"])
        cos = np.sum(np.multiply(find_position(np.int(result.iloc[idx]['copynum'])), find_position(np.int(result.iloc[idx+1]['copynum']))))
        cos_dis.append(cos)
plt.figure(dpi=800)
plt.hist(np.arccos(cos_dis),histtype='step',density=True,bins=20)
# plt.semilogy()
plt.xlabel('delta theta [rad]')
plt.ylabel('density')
#%%
plt.figure(dpi=800)
plt.hist(diff,histtype='step',density=True,bins=20)
plt.semilogy()
plt.xlabel('delta time [ns]')
plt.ylabel('density')

#%%
len_event = []
count_len = 1
multiple_hit = []
thres = 0
count_coin = 0
idx_list = []
r_list  =  []
for idx in range(len(result)-1):
    
    if result.iloc[idx]["batch"] == result.iloc[idx+1]["batch"]:
        # if result.iloc[idx]["copynum"] != result.iloc[idx+1]["copynum"]:
            # if isNeighbour(int(result_0_19.iloc[idx]["copynum"]), int(result_0_19.iloc[idx+1]["copynum"])):
            r_list.append(result.iloc[idx]["r"])
            r_list.append(result.iloc[idx+1]["r"])
            # idx_list.append(idx)
            # count_coin += 1
            count_len += 1
    else:
        # count_len += 1
        # len_event.append(count_len)
    #     if count_len > 2:
    #         multiple_hit.append(idx)
        count_len = 1
        # diff.append(result.iloc[idx+1]["t"] - result.iloc[idx]["t"])
        # cos = np.sum(np.multiply(find_position(np.int(result.iloc[idx]['copynum'])), find_position(np.int(result.iloc[idx+1]['copynum']))))
        # if cos < 0:
        #     print(result.iloc[idx+1])
        #     print(result.iloc[idx])
        #     cos_dis.append(cos)
# pmt_0 = result.loc[result['copynum']==0]
# coin_list = []
# neighbour_list_0 = []
# for pmt_id in range(1,31,1):
#     pmt_coin = result.loc[result['copynum']==pmt_id]
#     count = 0 
#     for i in pmt_0['batch']:
#         if i in pmt_coin['batch']:
#             count += 1
#     #print(count)
#     print('coincidence frequency(pmt 0 & {} only)/Hz: '.format(str(pmt_id)) + str(count/sim_time_wei))
#     #bool_neighbour = isNeighbour(0,pmt_id)
#     # neighbour_list_0.append(int(bool_neighbour))
#     coin_list.append(count)

# pmt_costheta_dis_xian = 1 - 3/2 * np.sin(33.96 / 180 *np.pi) ** 2#33.96 #0.5319
#0&12 0.6157
#0&19 0.2521
#0&24 -0.0871
costheta_list = []
for i in range(31):    
    costheta_list.append(np.sum(np.multiply(find_position(i), find_position(0))))
#%%
with open('/lustre/neutrino/weizhenyu/trigger/hdom_k40_new_optical/config/optical_properties.yaml', encoding='utf-8') as f:
    result = yaml.load(f.read(), Loader = yaml.FullLoader)
#'absorption_length', 'energy', 'mie_forward_angle', 'refractive_index', 'scatter_length_mie', 'scatter_length_rayeigh'

# plt.plot(result['energy'], result['absorption_length'])
wl = 1240 / np.array(result['energy'])
eff_att = (np.array(result['absorption_length']) ** (-1) + np.array(result['scatter_length_rayeigh']) ** (-1) + np.array(result['scatter_length_mie']) ** (-1))**(-1)
plt.plot(wl, result['absorption_length'])
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
# %%
coin_pmt_k40 = [1,2,3,4,5,6]
count_coin_k40 = [234600,5851, 306, 22, 4, 1]
err_coin = np.sqrt(count_coin_k40)

# %%
file_cqc_root_1 = np.load('/lustre/neutrino/weizhenyu/analysis_muon/data_npy/cqc_root_100.npy')
file_cqc_csv_1 = np.load('/lustre/neutrino/weizhenyu/analysis_muon/data_npy/cqc_csv_100.npy')
file_cqc_json_1 = np.load('/lustre/neutrino/weizhenyu/analysis_muon/data_npy/cqc_json_100.npy')
pmt_num_list = []
for id in range(100):#range(len(file_cqc_root_1)):#tqdm(range(10)):
    file_root = file_cqc_root_1[id]
    file_json = file_cqc_json_1[id]
    file_csv = file_cqc_csv_1[id]

    with open(file_json) as f:
        mc_events = json.load(f) 
    geo = pd.read_csv(file_csv,names=['x','y','z'])
    file = uproot.open(file_root)# energy unit is GeV
    tree = file[file.keys()[2]]

    hits =tree.arrays(library = 'pd')
    for eventid in np.array(hits.index.levels[0]):
        event = hits.loc[eventid]
        event = qe_module(event)
        if len(event) > 0:
            for domid in np.unique(np.array(event['DomId'])):
                event_dom = event.loc[event['DomId']==domid]
                num_pmt = len(np.unique(event_dom['PmtId']))
                pmt_num_list.append(num_pmt)
# pmt_num_list = np.array(pmt_num_list)[np.array(pmt_num_list)!=1]

# %%
file_name_10TeV = '/lustre/neutrino/zhangfuyudi/hailing/cascade/jobs_10tev_nue_10k_p1/job_0/data.root'
file = uproot.open(file_name_10TeV )# energy unit is GeV
tree = file[file.keys()[2]]
pmt_num_e = []
hits =tree.arrays(library = 'pd')
for eventid in np.array(hits.index.levels[0]):
    event = hits.loc[eventid]
    event = qe_module(event)
    if len(event) > 0:
        for domid in np.unique(np.array(event['DomId'])):
            event_dom = event.loc[event['DomId']==domid]
            num_pmt = len(np.unique(event_dom['PmtId']))
            pmt_num_e.append(num_pmt)
# pmt_num_e = np.array(pmt_num_e)[np.array(pmt_num_e)!=1]
#%%
range_ = [0.5,13.5]
bins_ = 13
plt.figure(dpi = 800, figsize=[8,6])
plt.hist(pmt_num_list,histtype='step', density=True,range = range_,bins = bins_,color = 'blue')
data_fake_k40 = []
for i in range(len([1,2,3,4,5,6])):
    data_fake_k40.extend(np.ones(count_coin_k40[i]) * coin_pmt_k40[i])
# plt.errorbar(x=coin_pmt_k40, y = err_coin, fmt='.', label = 'K40')
plt.hist(pmt_num_e,histtype='step', density=True,range = range_,bins = bins_,color = 'red')
plt.hist(data_fake_k40,histtype='step', density=True,range = range_,bins = bins_,color = 'green')
counts_k40, ticks = np.histogram(data_fake_k40, density=True,range = range_,bins = bins_)
counts_e, ticks_e = np.histogram(pmt_num_e, density=True,range = range_,bins = bins_)
counts_e_total, ticks_e = np.histogram(pmt_num_e, density=False,range = range_,bins = bins_)
error_e = np.power(counts_e_total,-0.5) * counts_e

x_axis = (ticks[1:-1] +ticks[0:-2])/2
x_axis = x_axis[0:6]

counts_mu, ticks_mu = np.histogram(pmt_num_list, density=True,range = range_,bins = bins_)
counts_mu_total, ticks_mu = np.histogram(pmt_num_list, density=False,range = range_,bins = bins_)
error_mu = np.power(counts_mu_total,-0.5) * counts_mu

error_K40 = np.power(count_coin_k40,-0.5) * counts_k40[0:6]
x_axis = (ticks[1:-1] +ticks[0:-2])/2
x_axis = x_axis[0:6]
plt.errorbar(x = x_axis , y = counts_k40[0:6], yerr=error_K40[0:6], fmt = '.',color = 'green',label = 'K40')
plt.errorbar(x = (ticks_mu[:-1] +ticks_mu[1:])/2 , y = counts_mu, yerr=error_mu, fmt = '.',color = 'blue',label = r"$\nu_{mu}$")
plt.errorbar(x = (ticks_e[:-1] +ticks_e[1:])/2 , y = counts_e, yerr=error_e, fmt = '.',color = 'red',label = r"$\nu_{E}$")
plt.semilogy()
plt.legend()
plt.xlabel('Coincidence PMT Number')
plt.ylabel('Density')
plt.xlim([0.5,13])
plt.axvline(np.quantile(data_fake_k40,0.98)-0.5,color='red',ls='--')
plt.text(np.quantile(data_fake_k40,0.98),0.3,'98% K40',color='red')
# %%
