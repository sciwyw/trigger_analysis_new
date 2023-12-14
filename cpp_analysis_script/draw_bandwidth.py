#! /usr/bin/env python3
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def tim(radius):
    N_event = 20000000 * 20 #BeamOn * Run number 
    N_event = np.random.poisson(lam = N_event)
    activtyK40 = 10.87
    effective = 1
    volume = 4 * np.pi * radius**3 / 3
    mass = 1.04* volume * 10**3
    frequncy = mass * activtyK40 * effective
    time = N_event / frequncy #* 10 ** 9
    print(time)
    return time

if __name__ == "__main__":
    wf_hz2MB = 31 * 500 * 16 / 1024 / 1024 / 8
    tdc_hz2MB = 24 * 100 / 8 / 1024 / 1024
    radius = 30
    sim_time = tim(radius)  # s  
    coin_num = [1, 2, 3, 4, 5, 6, 7]
    time_window = [1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]
    csvfile = "trigger_num.csv"
    trigger_list = pd.read_csv(csvfile, header=None)
    # run_num = len(trigger_list[0])
    # print("run %d times" % run_num)
    coin2 = []
    for i in range(16, 32):
        coin2.append(sum(trigger_list[i][:]))

    window20 = []
    for j in range(0, 7):
        window20.append(sum(trigger_list[16 * j + 10][:]))

    coin2_trigger_rate = np.divide(coin2, sim_time)
    window20_trigger_rate = np.divide(window20, sim_time)
    print(window20)

    coin22MB = coin2_trigger_rate * (wf_hz2MB + tdc_hz2MB)
    window202MB = window20_trigger_rate * (wf_hz2MB + tdc_hz2MB)

    plt.figure(dpi=800)
    plt.plot(time_window, coin22MB, "o--", color='blue')
    plt.xlabel("time window [ns]")
    plt.ylabel("MB/s per hDOM")
    plt.savefig("coin2.png")

    plt.figure(dpi=800)
    plt.plot(coin_num, window20_trigger_rate, "o--", color='blue')
    plt.xlabel("coincidence photon number")
    plt.ylabel("rate")
    plt.yscale('log')
    plt.savefig("window20.png")
