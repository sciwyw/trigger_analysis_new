#%%
import numpy as np
import matplotlib.pyplot as plt
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
wf_hz2MB = 31 * 500 * 16 /1024/1024/8
tdc_hz2MB = 24 * 100/8/1024/1024
# Bandwidth = (50k + 9kHz * 31pmt) * 100bit/Hz / 1024 / 1024/8 = 4 MB/s 
#%%
time_sim = 3.91
coin_num = 2
time_window_tts = [1,2,3,4,6,7,8,9,10,14,15,18,22,25,30]
hit_num_tts = np.array([16072,21778,22691,22620,22786,22756,22641,23284,22813,23339,23538,23379,23747,23505,23858]) 
coin_hit_rate_list_tts = np.divide(hit_num_tts,2*3.91)
errorbar_tts_rate = np.sqrt(hit_num_tts  / (3.91**2) + 400000 * 0.2**2 /3.91 )* (wf_hz2MB + tdc_hz2MB)
bd_k40_per_hdom = coin_hit_rate_list_tts * (wf_hz2MB + tdc_hz2MB)
plt.figure(dpi=800)
plt.errorbar(time_window_tts,bd_k40_per_hdom,yerr=errorbar_tts_rate,fmt='.',color = 'blue' )
plt.plot([1,2,3,30],np.array([2.05524297, 2.78491049,2.9016624,3.05089514])*1000*(wf_hz2MB + tdc_hz2MB),"--",color = 'blue')
plt.xlabel("time window [ns]")
plt.ylabel("MB/s per hDOM")
#%%
coin_num_list = [2,3,4,5,6]
coin_hit_rate_list = np.array([2.93,0.63,0.08,0.03,0.02])

bd_k40_per_hdom_coin_num = coin_hit_rate_list * 1000 * (wf_hz2MB + tdc_hz2MB)
errorbar_tts_rate_bd = np.sqrt((coin_hit_rate_list*1000)  / (3.91) + 400000 * 0.2**2 /3.91 ) * (wf_hz2MB + tdc_hz2MB)
rate_muon_total_array = 1400
N_photon_per_muon = 500

bd_muon_per_hDOM = N_photon_per_muon * rate_muon_total_array * (wf_hz2MB + tdc_hz2MB) / 22420
plt.figure(dpi=800)

plt.errorbar(coin_num_list,bd_k40_per_hdom_coin_num,yerr=errorbar_tts_rate_bd,color = 'blue',fmt='.' )
plt.plot(coin_num_list,bd_k40_per_hdom_coin_num,'--',color = 'blue' )
plt.xlabel("coincidence photon number in 10 ns")
plt.ylabel("MB/s per hDOM")
# %%
