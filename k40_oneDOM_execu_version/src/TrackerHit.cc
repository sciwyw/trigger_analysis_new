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
/// \file TrackerHit.cc
/// \brief Implementation of the TrackerHit class

#include "TrackerHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
// #include "MyAnalysis.hh"
//#include "G4AnalysisManager.hh"

#include "g4root.hh"
#include <iomanip>
#include "G4SystemOfUnits.hh"
//#include "DataCollector.hh"
G4ThreadLocal G4Allocator<TrackerHit> *TrackerHitAllocator = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerHit::TrackerHit()
    : G4VHit(),
      fEnergy(0), fTime(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerHit::~TrackerHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerHit::TrackerHit(const TrackerHit &right)
    : G4VHit()
{
  fTime = right.fTime;
  fEnergy = right.fEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const TrackerHit &TrackerHit::operator=(const TrackerHit &right)
{
  fTime = right.fTime;
  fEnergy = right.fEnergy;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TrackerHit::operator==(const TrackerHit &right) const
{
  return (this == &right) ? true : false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerHit::Print()
{
  auto *analysisManager = G4AnalysisManager::Instance();
  analysisManager ->FillNtupleDColumn(1,fTime / ns);
  analysisManager ->FillNtupleDColumn(2,cp_number);
  analysisManager ->FillNtupleDColumn(3,fEnergy / eV);
  analysisManager->FillNtupleIColumn(0, feventID);
  analysisManager->FillNtupleDColumn(4,fpos_x / mm);
  analysisManager->FillNtupleDColumn(5,fpos_y / mm);
  analysisManager->FillNtupleDColumn(6,fpos_z / mm);
  analysisManager->FillNtupleDColumn(7,fdecaytime / ns);
  // record information 
  analysisManager->AddNtupleRow();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
