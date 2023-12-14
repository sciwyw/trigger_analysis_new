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

#### Default 
- BeamOn = 20000000 (modify it in config.yaml) 
- world length = 30m (modify it in `k40_oneDOM_execu_version/src/DetetctorConstruction.cc&PrimaryGenerator.cc`
- stopandkill 

### Analysis 
The
