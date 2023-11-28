
#%%
import os
import sys
sys.path.append("/lustre/neutrino/weizhenyu/trigger/L1_L2")
#sys.path.remove('/lustre/neutrino/hufan/trident-analysis-py')
sys.path.append(os.getcwd())
from trident import Hits, Track, plot_track_3d, geomtry
import pkg_resources
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
from mpl_toolkits import mplot3d
from tqdm import tqdm
import json
import uproot
plt.clf()
import matplotlib as mpl
#/lustre/neutrino/zhangfuyudi/hailing/cascade/sy_data/nue_100tev_params/10k_p13/
from matplotlib import font_manager
for i in ['', 'bd', 'bi', 'i']:
    font_manager.fontManager.addfont('/home/zhiwei/fonts/times{}.ttf'.format(i))

plt.rcParams.update({
    'figure.figsize': (8, 5.5),
    'font.size': 20,
    'font.family': 'Times New Roman'
})
geo = pd.read_csv("/lustre/neutrino/weizhenyu/trigger/L1_L2/trident/data/penrose.csv",names=['x','y','z'])
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
size_x = 2200
size_y = 2200
size_z = 400
space = 100
time_window = 3000
K40_trigger_rate_series = [1500, 500, 50, 5]
raw_hit_fq = 5000*31
delta_t_mean = 0.6

light_dom_num = np.random.poisson(lam=time_window * raw_hit_fq * 10**(-9) * 24220)
k40_dom_id = np.random.randint(low=1, high=24221,size=light_dom_num)
fake_data = pd.DataFrame({'t0': np.random.uniform(low=0,high=time_window, size=len(k40_dom_id)),'DomId': k40_dom_id})
fake_data_ = fake_data
# print(len(fake_data))
for i in range(len(K40_trigger_rate_series)):
    
    N_dom = np.random.poisson(lam=time_window * K40_trigger_rate_series[i] * 10**(-9) * 24220)
    dom_id = np.random.randint(low=1, high=len(fake_data_),size=N_dom)
    fake_data_ = pd.DataFrame({'t0': np.random.exponential(scale=delta_t_mean, size=len(dom_id)) + np.array(fake_data_.iloc[dom_id]['t0']),'DomId': fake_data_.iloc[dom_id]['DomId']})
    fake_data = fake_data.append(fake_data_)
    # print(N_dom)
fake_data.sort_values(by = 'DomId', inplace=True)
fake_data.reset_index(drop=True,inplace=True)

#%%
size_x = 1400
size_y = 1400
size_z = 400
space = 100
fake_data_merge = getattr(fake_data.groupby("DomId"), "min")()
fake_data_merge['n'] = fake_data['DomId'].value_counts()
fake_data_merge.reset_index(drop = False,inplace = True)
# fake_data_merge_L1 = fake_data_merge.loc[fake_data_merge['n']>2]
fig = plt.figure(figsize=(18, 15), dpi=600)
ax = fig.add_axes([0.07, 0.3, 0.8, 0.6], projection='3d')
ax.set_box_aspect(aspect = (5, 5, 3))
ax.scatter(geo.iloc[fake_data_merge['DomId']]["x"], geo.iloc[fake_data_merge['DomId']]["y"], geo.iloc[fake_data_merge['DomId']]["z"],
            s=np.power(fake_data_merge["n"], .5) * 9.,
            c=fake_data_merge["t0"],
            marker='.',
            alpha=1.,
            cmap='rainbow_r')
ax.scatter(geo["x"], geo["y"], geo["z"],
            s=0.1,
            color = 'grey',
            marker='.',
            alpha=0.5,)
            # cmap='rainbow_r')
cmap_ = plt.get_cmap('rainbow_r')
ax.set_xlim(-size_x-space, size_x+space)
ax.set_ylim(-size_y-space, size_y+space)
ax.set_zlim(-size_z-space, size_z+space)
ax.set_xticks(np.linspace(-size_x, size_x, 4))
ax.set_yticks(np.linspace(-size_y, size_y, 4))
ax.set_zticks(np.linspace(-size_z, size_z, 2))
ax.set_xlabel("X [m]")
ax.set_ylabel("Y [m]")
ax.set_zlabel("Z [m]")
ax.grid(False)

    # fig.add_axes([left, bottom, width, height])
cax = fig.add_axes([0.2, 0.1, 0.6, 0.05])
length = int(len(fake_data_merge) / 20)
fake_data_merge.sort_values("t0", inplace=True)
minT = fake_data_merge["t0"].iloc[length]
maxT = fake_data_merge["t0"].iloc[-length-1]
norm = mpl.colors.Normalize(vmin=minT, vmax=maxT)
cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap_,
                                norm=norm,
                                orientation='horizontal')
cax.set_xlabel('Time [ns]')
plt.show()
#%%
fake_data_merge_L1 = fake_data_merge.loc[fake_data_merge['n']>1]
fig = plt.figure(figsize=(18, 15), dpi=600)
ax = fig.add_axes([0.07, 0.3, 0.8, 0.6], projection='3d')
ax.scatter(geo.iloc[fake_data_merge_L1['DomId']]["x"], geo.iloc[fake_data_merge_L1['DomId']]["y"], geo.iloc[fake_data_merge_L1['DomId']]["z"],
            s=np.power(fake_data_merge_L1["n"], .5) * 9.,
            c=fake_data_merge_L1["t0"],
            marker='.',
            alpha=1.,
            cmap='rainbow_r')
# ax.scatter(geo["x"], geo["y"], geo["z"],
#             s=0.1,
#             color = 'grey',
#             marker='.',
#             alpha=0.5,)
# ax.scatter(geo["x"], geo["y"], geo["z"],
#             s=0.1,
#             color = 'grey',
#             marker='.',
#             alpha=0.5,)
t = np.linspace(-size_x-space, size_x + space, 100)

# ax.plot3D(xs = eventinfo['x'] + eventinfo['px']/10000 * t, ys = eventinfo['y'] + eventinfo['py']/10000 * t,zs = eventinfo['z'] + eventinfo['pz']/10000 * t, linewidth = 0.7)
ax.set_box_aspect(aspect = (5, 5, 3))
cmap_ = plt.get_cmap('rainbow_r')
ax.set_xlim(-size_x-space, size_x+space)
ax.set_ylim(-size_y-space, size_y+space)
ax.set_zlim(-size_z-space, size_z+space)
ax.set_xticks(np.linspace(-size_x, size_x, 4))
ax.set_yticks(np.linspace(-size_y, size_y, 4))
ax.set_zticks(np.linspace(-size_z, size_z, 2))
ax.set_xlabel("X [m]")
ax.set_ylabel("Y [m]")
ax.set_zlabel("Z [m]")
ax.grid(False)

    # fig.add_axes([left, bottom, width, height])
cax = fig.add_axes([0.2, 0.1, 0.6, 0.05])
length = int(len(fake_data_merge_L1) / 20)
fake_data_merge_L1.sort_values("t0", inplace=True)
minT = fake_data_merge_L1["t0"].iloc[length]
maxT = fake_data_merge_L1["t0"].iloc[-length-1]
norm = mpl.colors.Normalize(vmin=minT, vmax=maxT)
cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap_,
                                norm=norm,
                                orientation='horizontal')
cax.set_xlabel('Time [ns]')
plt.show()
# %%
muon_data_path = '/lustre/neutrino/changqichao/icrc2023/telescope/penrose/1.0/job_0/data/data.root'
with open('/lustre/neutrino/changqichao/icrc2023/telescope/penrose/1.0/job_0/mc_events.json') as f:
            mc_events = json.load(f) 
file = uproot.open(muon_data_path)
tree = file[file.keys()[2]]
hits =tree.arrays(library = 'pd')
num_list = []
for i in np.unique(hits.index.get_level_values(0)):
    hits_event = hits.loc[i]
    num = len(np.unique(hits_event['DomId']))
    num_list.append(num)
event_idx = np.argmax(num_list)
event_id =  np.unique(hits.index.get_level_values(0))[event_idx]
eventinfo = mc_events[event_id]['particles_out'][0]
fake_data = hits.loc[event_id]
fake_data['t0'] = fake_data['t0'] - np.min(fake_data['t0'])
#%%

fake_data_merge = getattr(fake_data.groupby("DomId"), "min")()
fake_data_merge['n'] = fake_data['DomId'].value_counts()
fake_data_merge.reset_index(drop = False,inplace = True)
size_x = np.max(geo.iloc[fake_data_merge['DomId']]["x"])
size_y = np.max(geo.iloc[fake_data_merge['DomId']]["y"])
size_z = np.max(geo.iloc[fake_data_merge['DomId']]["z"])
space = 100
# fake_data_merge_L1 = fake_data_merge.loc[fake_data_merge['n']>2]
fig = plt.figure(figsize=(18, 15), dpi=600)
ax = fig.add_axes([0.07, 0.3, 0.8, 0.6], projection='3d')
ax.set_box_aspect(aspect = (5, 5, 3))
ax.scatter(geo.iloc[fake_data_merge['DomId']]["x"], geo.iloc[fake_data_merge['DomId']]["y"], geo.iloc[fake_data_merge['DomId']]["z"],
            s=np.power(fake_data_merge["n"], .5) * 9.,
            c=fake_data_merge["t0"],
            marker='.',
            alpha=1.,
            cmap='rainbow_r')
ax.scatter(geo["x"], geo["y"], geo["z"],
            s=0.1,
            color = 'grey',
            marker='.',
            alpha=0.5,)
t = np.linspace(-size_x-space, size_x + space, 100)
ax.plot3D(xs = eventinfo['x'] + eventinfo['px']/10000 * t, ys = eventinfo['y'] + eventinfo['py']/10000 * t,zs = eventinfo['z'] + eventinfo['pz']/10000 * t, linewidth = 0.7)
            # cmap='rainbow_r')
cmap_ = plt.get_cmap('rainbow_r')
ax.set_xlim(-size_x-space, size_x+space)
ax.set_ylim(-size_y-space, size_y+space)
ax.set_zlim(-size_z-space, size_z+space)
ax.set_xticks(np.linspace(-size_x, size_x, 4))
ax.set_yticks(np.linspace(-size_y, size_y, 4))
ax.set_zticks(np.linspace(-size_z, size_z, 2))
ax.set_xlabel("X [m]")
ax.set_ylabel("Y [m]")
ax.set_zlabel("Z [m]")
ax.grid(False)

    # fig.add_axes([left, bottom, width, height])
cax = fig.add_axes([0.2, 0.1, 0.6, 0.05])
length = int(len(fake_data_merge) / 20)
fake_data_merge.sort_values("t0", inplace=True)
minT = fake_data_merge["t0"].iloc[length]
maxT = fake_data_merge["t0"].iloc[-length-1]
norm = mpl.colors.Normalize(vmin=minT, vmax=maxT)
cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap_,
                                norm=norm,
                                orientation='horizontal')
cax.set_xlabel('Time [ns]')
plt.show()
# %%
fake_data_merge_L1 = fake_data_merge.loc[fake_data_merge['n']>1]
fig = plt.figure(figsize=(18, 15), dpi=600)
ax = fig.add_axes([0.07, 0.3, 0.8, 0.6], projection='3d')
ax.scatter(geo.iloc[fake_data_merge_L1['DomId']]["x"], geo.iloc[fake_data_merge_L1['DomId']]["y"], geo.iloc[fake_data_merge_L1['DomId']]["z"],
            s=np.power(fake_data_merge_L1["n"], .5) * 9.,
            c=fake_data_merge_L1["t0"],
            marker='.',
            alpha=1.,
            cmap='rainbow_r')
ax.scatter(geo["x"], geo["y"], geo["z"],
            s=0.1,
            color = 'grey',
            marker='.',
            alpha=0.5,)
ax.scatter(geo["x"], geo["y"], geo["z"],
            s=0.1,
            color = 'grey',
            marker='.',
            alpha=0.5,)
t = np.linspace(-size_x-space, size_x + space, 100)

ax.plot3D(xs = eventinfo['x'] + eventinfo['px']/10000 * t, ys = eventinfo['y'] + eventinfo['py']/10000 * t,zs = eventinfo['z'] + eventinfo['pz']/10000 * t, linewidth = 0.7)
ax.set_box_aspect(aspect = (5, 5, 3))
cmap_ = plt.get_cmap('rainbow_r')
ax.set_xlim(-size_x-space, size_x+space)
ax.set_ylim(-size_y-space, size_y+space)
ax.set_zlim(-size_z-space, size_z+space)
ax.set_xticks(np.linspace(-size_x, size_x, 4))
ax.set_yticks(np.linspace(-size_y, size_y, 4))
ax.set_zticks(np.linspace(-size_z, size_z, 2))
ax.set_xlabel("X [m]")
ax.set_ylabel("Y [m]")
ax.set_zlabel("Z [m]")
ax.grid(False)

    # fig.add_axes([left, bottom, width, height])
cax = fig.add_axes([0.2, 0.1, 0.6, 0.05])
length = int(len(fake_data_merge_L1) / 20)
fake_data_merge_L1.sort_values("t0", inplace=True)
minT = fake_data_merge_L1["t0"].iloc[length]
maxT = fake_data_merge_L1["t0"].iloc[-length-1]
norm = mpl.colors.Normalize(vmin=minT, vmax=maxT)
cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap_,
                                norm=norm,
                                orientation='horizontal')
cax.set_xlabel('Time [ns]')
plt.show()
# %%
