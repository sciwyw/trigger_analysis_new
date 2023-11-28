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
import scipy.interpolate
from mpl_toolkits import mplot3d
plt.clf()
from matplotlib import font_manager
for i in ['', 'bd', 'bi', 'i']:
    font_manager.fontManager.addfont('/home/zhiwei/fonts/times{}.ttf'.format(i))

plt.rcParams.update({
    'figure.figsize': (8, 5.5),
    'font.size': 20,
    'font.family': 'Times New Roman'
})

data_10 = '/lustre/neutrino/weizhenyu/trigger/check_k40_photon_new_optical/build/10/data.root'
data_20 = '/lustre/neutrino/weizhenyu/trigger/check_k40_photon_new_optical/build/20/data.root'
#%%
file_10 = uproot.open(data_10)
tree_10 = file_10[file_10.keys()[0]]
hits_10 = tree_10.arrays(library = 'pd')

file_20 = uproot.open(data_20)
tree_20 = file_20[file_20.keys()[0]]
hits_20 = tree_20.arrays(library = 'pd')
# %%
range_ = [40,400]
bins_ = 50
plt.hist(hits_10['t(ns)'], histtype='step',range=range_, bins= bins_)
plt.hist(hits_20['t(ns)'], histtype='step',range=range_, bins= bins_)
plt.semilogy()
plt.xlabel('Time [ns]')
plt.ylabel('Count')
# %%
range_ = [200,800]
bins_ = 50
plt.hist(np.divide(1240, hits_10['energy(eV)']), histtype='step',range=range_, bins= bins_)
plt.hist(np.divide(1240, hits_20['energy(eV)']), histtype='step',range=range_, bins= bins_)
plt.xlabel('Wavelength [nm]')
plt.ylabel('Count')
# %%
N_event = 10000
N_10 = len(hits_10)
N_20 = len(hits_20)
lambda_eff = 10 / np.log(N_10 / N_20)
N_abs = N_10 * N_10 / N_20
print('photon_num at 10m is')
print(N_10)
print('photon_num at 20m is')
print(N_20)
print('lambda is' )
print(lambda_eff)
print('N_abs is' )
print(N_abs / N_event)
# %%
