# trigger_analysis
Three parts are included. 
1. K40 and dark noise simulation (Geant4) & analysis in analysis_trigger_root_tts.py, bandwidth_estimation.py and analysis_trigger_dark_noise.py.

   Analysis target: Count coincidence photon events in hDOM. 
   
3. L1 analysis(muon and shower events are included) in analysis_muon and analysis_shower

   Analysis target: Count hDOM trigger efficiency for certain energy events at certain distance.
   
5. L2 analysis(only muon threshold analysis and fake K40 noise generator are included) in analysis_trigger_L2

   Analysis target: Obtain muon trigger efficiency for certain energy events.

Trigger is a game about counting some number in certain time window.
Code here is a little sloppy, but the basic workflow is simple and similar. 

step 1: Read information from data.root

example:
```
file = uproot.open(file_name)
tree = file[file.keys()[0]]
hits =tree.arrays(library = 'pd')
```

step 2: count number in certain time window (Just define your own ways to count photons. This is trigger.)

example: this funtion can decide if a hDOM is triggered.
```
def Is_L1_trigger(data:pd.DataFrame):
    if len(data) < coin_num:
        return False
    else:
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
```
step 3: make some plots from step 2

example:
```
plt.hist(data, range, bins,histtype = 'step')
```
All analysis is written in python.
