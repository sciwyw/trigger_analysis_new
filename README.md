# trigger_analysis
Three parts are included. 
1. K40 simulation (Geant4) & analysis
   Analysis target: Count coincidence photon events in hDOM. 
   
3. L1 analysis(muon and shower events are included)
   Analysis target: Count hDOM trigger efficiency for certain energy events at certain distance.
   
5. L2 analysis(only muon threshold analysis is included)
   Analysis target: Obtain muon trigger efficiency for certain energy events.

Trigger is a game about counting some number in certain time window.

Code here is mazy and sloppy, but the basic workflow is simple and similar. 

step 1: Read information from data.root

step 2: count number in certain time window (Just define your own ways to count photons. This is trigger.)

step 3: make some plots from step 2

All analysis is written in python.
