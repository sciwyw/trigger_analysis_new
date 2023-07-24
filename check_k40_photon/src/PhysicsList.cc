/// \file Physicslist
/// \brief allowed physical process

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4StoppingPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4IonPhysics.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmParameters.hh"

#include "PhysicsList.hh"

PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
    SetVerboseLevel(1);

    //defult physics
    RegisterPhysics(new G4DecayPhysics());
    RegisterPhysics(new G4EmStandardPhysics_option1());
    RegisterPhysics(new G4HadronElasticPhysics());
    RegisterPhysics(new G4HadronPhysicsFTFP_BERT());
    RegisterPhysics(new G4StoppingPhysics());
    RegisterPhysics(new G4IonPhysics());

    //neutron tracking cut
    G4NeutronTrackingCut *neutronTrackingCut = new G4NeutronTrackingCut();
    neutronTrackingCut->SetKineticEnergyLimit(1 * GeV);
    RegisterPhysics(neutronTrackingCut);
    RegisterPhysics(new G4RadioactiveDecayPhysics());
    RegisterPhysics(new G4OpticalPhysics());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts()
{
    G4double cutForGamma = 0.3 * mm;
    G4double cutForElectron = 0.3 * mm;  // 1cm for 2 MeV cut raughly
    G4double cutForPositron = 0.3 * mm;
    SetCutValue(cutForGamma,"gamma");
    SetCutValue(cutForElectron,"e-");
    SetCutValue(cutForPositron,"e+");
    DumpCutValuesTable();
}