import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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
time_window = [1,2,6,10,14,18,22]#,300]
rate = [0.103,0.152,0.240,0.277,0.304,0.317,0.334]#,0.410]
error_bar = np.sqrt(np.array(rate)*1000) / 1000 
plt.figure(dpi = 800)
plt.plot(time_window,rate,'--',color='blue')
plt.errorbar(time_window,rate,error_bar,fmt='.',color='blue')
# plt.semilogx()
# plt.ylim([0,0.6])
plt.xlabel("time window [ns]")
plt.ylabel("mean trigger ratio")
