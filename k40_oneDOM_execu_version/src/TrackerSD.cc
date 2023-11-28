//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file TrackerSD.cc
/// \brief Implementation of the TrackerSD class

#include "g4root.hh"
#include "TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
//#include "G4Event.hh"
//#include "G4EventManager.hh"
//#include "G4UserEventAction.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TrackerSD::TrackerSD(const G4String &name,
                         const G4String &hitsCollectionName)
    : G4VSensitiveDetector(name),
      fHitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerSD::~TrackerSD()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::Initialize(G4HCofThisEvent *hce)
{
  // Create hits collection

  fHitsCollection = new TrackerHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce

  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection(hcID, fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TrackerSD::ProcessHits(G4Step *aStep,
                                G4TouchableHistory *)
{
  // paticle find
  //const G4Event *event;
  const G4Event* event = G4RunManager::GetRunManager()->GetCurrentEvent();
  G4Track *pTrack;
  pTrack = aStep->GetTrack();
  G4String ParticleName;
  ParticleName = pTrack->GetParticleDefinition()->GetParticleName();
  G4int eventid = event->GetEventID();
  if (ParticleName == "opticalphoton")
  {
    //G4cout << "Oops!Hit!"<< G4endl;
    G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
    G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
    G4int copyNo = theTouchable->GetReplicaNumber();
    TrackerHit *newHit = new TrackerHit();

    newHit->SetTime(aStep->GetPreStepPoint()->GetGlobalTime());
    newHit->SetcpN(copyNo);//newHit->SetcpN_phi(copyNo);
    newHit->SetEnergy(aStep->GetPreStepPoint()->GetKineticEnergy());
    newHit->SetEventid(eventid);
    newHit->SetDecaytime(event->GetPrimaryVertex()->GetT0());
    newHit->SetPosx(event->GetPrimaryVertex()->GetPosition()[0]);
    newHit->SetPosy(event->GetPrimaryVertex()->GetPosition()[1]);
    newHit->SetPosz(event->GetPrimaryVertex()->GetPosition()[2]);


    fHitsCollection->insert(newHit);
    newHit->Print();
    aStep->GetTrack()->SetTrackStatus(fStopAndKill);
  }
  
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::EndOfEvent(G4HCofThisEvent *)
{
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
