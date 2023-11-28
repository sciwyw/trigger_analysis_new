##### Default:
BeamOn = 20000000 (see config/config.yaml)

absorb photon (see `aStep->GetTrack()->SetTrackStatus(fStopAndKill);` in TrackerSD.cc)

##### ROOT file naming convention:
Tree: `Hit0`
Branch: `eventid` `hittime` `PMTid` `energy` `x` `y` `z` `decaytime`

##### Notice:
`hittime` includes decaytime. If you want the time interval between primary vertex and hit, `hittime` - `decaytime`. 

Doesn't place Gel in DetectorConstruction
