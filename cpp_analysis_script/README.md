### Simulation 
The simulation code: `k40_oneDOM_execu_version`  
Run the simulation by: 
```shell 
mkdir build 
cd build 
cmake Path/To/k40_oneDOM_execu_version 
make             
./main [config.yaml] [optical_properties.yaml] 
```
the `.yaml` files could be found in `k40_oneDOM_execu_version/config`
#### Output
the output is a ROOT file `data.root`

#### Default 
- BeamOn = 20000000 (modify it in config.yaml) 
- world length = 30m (modify it in `k40_oneDOM_execu_version/src/DetetctorConstruction.cc&PrimaryGenerator.cc`
- stopandkill 

### Analysis 
The analysis code: `cpp_analysis_script` 
Run analysis by:
```shell
mkdir build
cd build
cmake Path/To/cpp_analysis_script
make
hadd copy.root data.root
./analysis copy.root
```
#### Output              
The output is a csv file (**remember to change the path of the .csv file in main() of analysis_trigger_root_tts.cpp**) 
#### Default
- only analyze time window = 20ns (add or edit in the .cpp file)
#### Draw
`python draw_bandwidth.py`

