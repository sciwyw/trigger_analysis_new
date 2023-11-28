#31 same dark count rate
#only adjacent pmt hits are considered as triggered
#only two neighbour hits are recorded 

import numpy as np
import pandas as pd
import uproot
import os
import sys
import glob
import matplotlib.pyplot as plt
from scipy import optimize as op
from cProfile import label
import matplotlib
from mpl_toolkits import mplot3d
import scipy.interpolate
import argparse
# check the number in the primarygenerator action
parser = argparse.ArgumentParser()
parser.add_argument(
    'time_window', type=int,
    help='time_window')
parser.add_argument(
    'coin_num',  type=int,
    help='coincidence number')
args = parser.parse_args()
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

def addGlobalTime(hits: pd.DataFrame, time_per_event: float):
    event_id_pre = 0
    event_time_pre = 0
    for idx in range(len(hits)):
        event_id_cur = hits.loc[idx, "eventID"]
        if event_id_cur == event_id_pre:
            hits.loc[idx, "t(ns)"] += event_time_pre
        else:
            delta_e = event_id_cur - event_id_pre
            delta_t = np.random.exponential(delta_e * time_per_event)
            event_time_pre += delta_t
            event_id_pre = event_id_cur
            hits.loc[idx, "t(ns)"] += event_time_pre
    return hits

def read_hits(file_name: str,bool_tts,tts) -> pd.DataFrame:
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
    
    if bool_tts:
        tts_ = np.random.normal(loc = 0,scale=tts / 2.3548,size=len(hits))
        hits["t(ns)"] = hits["t(ns)"] + tts_
    hits["r"] = np.sqrt(
        np.power(hits["x(mm)"], 2) +
        np.power(hits["y(mm)"], 2) +
        np.power(hits["z(mm)"], 2)) / 1000
    return hits

def qe_module(hit: pd.DataFrame):
    #quantum efficiency
    # 1eV ~ 1240nm
    qe_list = pd.read_csv('qe.csv')
    qe_func = scipy.interpolate.interp1d(qe_list['wl'],qe_list['qe'],kind='linear')
    wl_hit = 1240 / np.array(hit["energy(eV)"])
    qe_photon = qe_func(wl_hit) / 100
    sample = np.random.sample(len(hit))
    survival_bool = qe_photon > sample
    hit_qe = hit[survival_bool]
    #print(survival_bool)
    hit_qe = hit_qe.reset_index(drop = True)

    return hit_qe

def test_radius_distribution(hits: pd.DataFrame):
    radius_list = []
    event_id_pre = 0
    time_pre = -1000
    for i, hit in hits.iterrows():
        if hit["eventID"] == event_id_pre:
            if hit["t(ns)"] - time_pre < 10:
                radius_list.append(hit["r"])
            else:
                time_pre = hit["t(ns)"]
        else:
            event_id_pre = hit["eventID"]
            time_pre = -1000
    plt.hist(radius_list, bins=np.arange(0, 30, 1))
    plt.show()

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
if __name__=="__main__":
    # N_event = 10000000
    N_event = 10000000
    dcr = 1000 #dark count rate for pmts is 1000Hz
    radius_sim = 30
    TTS = 1.3
    for i in range(0,500,1):
        path_data = "/lustre/neutrino/weizhenyu/trigger/hdom_k40/build/lib/{}/data.root".format(i)
        #path_data_gamma = "build_iso/data_nt_1.csv"
        # path_data = "/lustre/neutrino/weizhenyu/trigger/hdom_k40_new_optical/build_lib/run/{}/".format(i)
        # path_data = '/lustre/neutrino/weizhenyu/trigger/hdom_k40_new_optical/build_200/run/{}/'.format(i)
        path_save_data = "/lustre/neutrino/weizhenyu/trigger/analysis/result/trigger_hdom_{}_{}_tts_0711.csv".format(args.time_window, args.coin_num)

        time_total_sample = getRealTime_Orb(N_event, radius_sim) # check the radius in the .mac file
        # time_series_e = pd.DataFrame(
        #     {'shift': np.random.uniform(0,time_total_sample,react_e)}).sort_values(by = 'shift', ascending=True).reset_index(drop=True)
        # time_series_gamma = pd.DataFrame({'shift': np.random.uniform(0,time_total_sample,react_gamma)}).sort_values(by = 'shift', ascending=True).reset_index(drop=True)
        #sample reaction time

        data = read_hits(path_data,True,TTS)
        #data_2 = read_hits(path_data_gamma,True)


        #shift time according to poisson and unifrom distribution
        #sample the time of k-40 decay reaction
        addGlobalTime(data, time_total_sample/N_event)
        # addGlobalTime(data_2, time_total_sample/react_gamma)

        # data = data_1.append(data_2)
        # data = data_1.append(data_2)

        #data.insert(7,'r',(data["x(mm)"]**2 + data["y(mm)"]**2 +data["z(mm)"]**2)**0.5)
        # data.insert(7,'r',(data["x(mm)"]**2 + data["y(mm)"]**2 +data["z(mm)"]**2)**0.5)
        # ratio = len(data.loc[data['r']<30000])/len(data)
        # plt.axvline(x=30,label=str(np.round(ratio*100,2))+'%'+' (0,30m)',color='r')
        # plt.hist(data['r']/1000,range=[0.215,60],histtype='step',label='reaction distance',linewidth=2)
        # plt.legend()
        # plt.xlabel('r [m]')
        # plt.ylabel('Count')

        # sample dark noise
        # for i in range(31): #31pmt in total
        #     dc = np.random.poisson(dcr * time_total_sample / 10 ** 9) #sample the number of total dark photon
        #     dp = np.sort(np.random.uniform(0,time_total_sample, dc)) #sample the time of dark photon
        #     copy_num = np.ones(dc).astype(int) * i 
        #     dark_count = pd.DataFrame({"t(ns)":dp,"copynumber":copy_num})
        #     data = data.append(dark_count)
        # data_0 = data.loc[data["copynumber"]==0]
        # data_12 = data.loc[data["copynumber"]==12]
        # data = data_0.append(data_12)
        # data.reset_index(drop = True, inplace = True)
        # data.sort_values(by = ["t(ns)"], ascending=True, inplace=True)
        
        thres = 0 # costheta threshold to decide if it's neighbour
        time_window = args.time_window # time window to decide if it's coincidence
        coin_num = args.coin_num
        diff_data_t = np.array(data["t(ns)"][coin_num-1:len(data)]) - np.array(data["t(ns)"][0:len(data)-(coin_num-1)]) # The last term subtracts the first term
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

        data_trigger = data.loc[trigger_array]
        data_trigger.insert(7, 'batch', np.ones(len(data_trigger)))
        data_trigger = data_trigger.reset_index(drop = True)
    
        for ii in range(1, len(data_trigger)-(coin_num-1), 1):
            if (data_trigger["t(ns)"][ii+(coin_num-1)] - data_trigger["t(ns)"][ii]) < args.time_window:
                temp = data_trigger['batch'][ii] 
                for k0 in range((coin_num-1)):
                    data_trigger.loc[ii+k0+1,'batch'] = temp
            else:
                temp = data_trigger['batch'][ii] + 1
                data_trigger.loc[ii+1,'batch'] = temp

        # for i in range(1, np.int(np.max(data_trigger['batch'])),1):
            
        #     data_batch = data_trigger.loc[data_trigger['batch']==i]
        #     copynum_list = np.array(data_batch["copynumber"])
        #     if not isNeighbour(int(copynum_list[0]),int(copynum_list[1])):
        #         data_trigger.drop(data_batch.index)

        data_trigger.to_csv(path_save_data,header=False, index=True,mode = 'a')