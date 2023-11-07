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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeneralParticleSource.hh"
#include "Randomize.hh"
#include "G4SPSEneDistribution.hh"
#include "g4root.hh"
#include <iostream>
#include <fstream>
// #include "G4INCLRandom.hh"
//#include "G4AnalysisManager.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction()
{
  //set electron source
  // G4cout << "PrimaryGeneratorAction::PrimaryGeneratorAction():G4VUserPrimaryGeneratorAction()" << G4endl;
  fParticleSource = new G4GeneralParticleSource;
  G4ParticleDefinition* particleDefinition=G4ParticleTable::GetParticleTable()->FindParticle("e-");
  fParticleSource->SetParticleDefinition(particleDefinition);
  InitFunction();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  //G4cout << "delete primary generator" << G4endl;
  delete fParticleSource;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  //G4cout << "GeneratePrimaries(G4Event* anEvent)" << G4endl;
  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
  G4int event_num = 1;
  G4double pos_x, pos_y, pos_z;
  G4double e_branch_ratio = 0.893;
  G4double energy_gamma = 1.459 * MeV;
  //G4SPSEneDistribution energy_e_spec;// = new G4SPSEneDistribution();
  G4double probability;
  G4double sample_r;
  G4ThreeVector dir;
  G4double px_e, py_e, pz_e;
  G4double cos_theta_r, sin_theta_r, phi_r;
  G4double energy_e; 
  G4double cos_theta_p, sin_theta_p, phi_p;
  G4double px_gamma, py_gamma,pz_gamma;
  G4PrimaryParticle *particle;
  
  G4PrimaryVertex **fVecPrimaryVertex;
  fVecPrimaryVertex = new G4PrimaryVertex *[event_num];
  //G4cout << "G4PrimaryVertex **fVecPrimaryVertex;" << G4endl;
  for(G4int event_id = 0; event_id < event_num; event_id++){
    fVecPrimaryVertex[event_id] = new G4PrimaryVertex;
    //G4cout << "event_id " << event_id << G4endl;
    //G4cout << "fVecPrimaryVertex[event_id] = new G4PrimaryVertex;" << G4endl;
    sample_r = cbrt(G4RandFlat::shoot(0.215 * 0.215 * 0.215, 30. * 30. * 30)) * m;
    G4cout << sample_r << G4endl;
    //theta_r = G4RandFlat::shoot(0., M_PI);
    phi_r = G4RandFlat::shoot(0.,2 * M_PI);
    cos_theta_r = G4RandFlat::shoot(-1, 1.);
    sin_theta_r = sqrt(1 - cos_theta_r * cos_theta_r);
    pos_x = sample_r * sin_theta_r * cos(phi_r);
    pos_y = sample_r * sin_theta_r * sin(phi_r);
    pos_z = sample_r * cos_theta_r; 
    G4cout << pos_x << G4endl;
    G4cout << pos_y << G4endl;
    G4cout << pos_z << G4endl;
    analysisManager->FillNtupleDColumn(4,pos_x);
    analysisManager->FillNtupleDColumn(5,pos_y);
    analysisManager->FillNtupleDColumn(6,pos_z);
    probability = G4RandFlat::shoot(0., 1.);
    if (probability < e_branch_ratio) {
      // energy_e = energy_e_spec.GenerateOne(G4ParticleTable::GetParticleTable()->FindParticle("e-"));
      energy_e = InverseCumul();
      
      outfile.open("energy_spec_e.csv", std::ios_base::app);
      outfile << energy_e << G4endl;
      
      fVecPrimaryVertex[event_id]->SetPosition(pos_x, pos_y, pos_z);

      cos_theta_p = G4RandFlat::shoot(-1, 1.);
      sin_theta_p = sqrt(1 - cos_theta_p * cos_theta_p);
      phi_p = G4RandFlat::shoot(0.,2 * M_PI);
      dir = G4ThreeVector(sin_theta_p * cos(phi_p),sin_theta_p * sin(phi_p),cos_theta_p);
      px_e = (0.511 + energy_e) * sqrt(1-(0.511 / (energy_e + 0.511)) * (0.511 / (energy_e + 0.511))) * sin_theta_p * cos(phi_p);
      py_e = (0.511 + energy_e) * sqrt(1-(0.511 / (energy_e + 0.511)) * (0.511 / (energy_e + 0.511))) * sin_theta_p * sin(phi_p);
      pz_e = (0.511 + energy_e) * sqrt(1-(0.511 / (energy_e + 0.511)) * (0.511 / (energy_e + 0.511))) * cos_theta_p;
      G4cout << px_e << G4endl;
      G4cout << py_e << G4endl;
      G4cout << pz_e << G4endl;
      particle =
          new G4PrimaryParticle(G4ParticleTable::GetParticleTable()->FindParticle("e-"), px_e, py_e, pz_e);
      // G4cout << energy_e << G4endl;
      //G4cout << "dir " << dir << G4endl;
      fVecPrimaryVertex[event_id]->SetT0(0. * ns);
      fVecPrimaryVertex[event_id]->SetPrimary(particle);
      // fParticleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("e-"));
      // fParticleGun->SetParticleMomentumDirection(dir);
      // fParticleGun->SetParticleEnergy(energy_e);    
      // fParticleGun->GeneratePrimaryVertex(anEvent);
      //G4cout << " fParticleGun->GeneratePrimaryVertex(anEvent); " << G4endl;
    }
    else {
      cos_theta_p = G4RandFlat::shoot(-1, 1.);
      sin_theta_p = sqrt(1 - cos_theta_p * cos_theta_p);
      phi_p = G4RandFlat::shoot(0.,2 * M_PI);
      dir = G4ThreeVector(sin_theta_p * cos(phi_p),sin_theta_p * sin(phi_p),cos_theta_p);
      px_gamma = energy_gamma * sin_theta_p * cos(phi_p);
      py_gamma = energy_gamma * sin_theta_p * sin(phi_p);
      pz_gamma = energy_gamma * cos_theta_p;
      G4cout << px_gamma << G4endl;
      G4cout << py_gamma << G4endl;
      G4cout << pz_gamma << G4endl;
      particle =
          new G4PrimaryParticle(G4ParticleTable::GetParticleTable()->FindParticle("gamma"), px_gamma, py_gamma, pz_gamma);

      fVecPrimaryVertex[event_id]->SetPosition(pos_x, pos_y, pos_z);

      fVecPrimaryVertex[event_id]->SetT0(0. * ns);
      fVecPrimaryVertex[event_id]->SetPrimary(particle);
      // fParticleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("gamma"));
      // fParticleGun->SetParticleMomentumDirection(dir);
      // fParticleGun->SetParticleEnergy(energy_gamma);    
      // fParticleGun->GeneratePrimaryVertex(anEvent);
      //G4cout << " fParticleGun->GeneratePrimaryVertex(anEvent); " << G4endl;
      //G4cout << "fVecPrimaryVertex[event_id]->SetPrimary(particle);" << G4endl;
    }
  }
  for (int i = 0; i < event_num; i++){
    //G4cout << "anEvent->AddPrimaryVertex(fVecPrimaryVertex[i]);" << G4endl;
    anEvent->AddPrimaryVertex(fVecPrimaryVertex[i]);
  }


  // fParticleSource->GeneratePrimaryVertex(anEvent);
  // G4ThreeVector position = fParticleSource->GetParticlePosition();
  // G4double radius = sqrt(position.getX()*position.getX() + 
  //   position.getY()*position.getY() + 
  //   position.getZ()*position.getZ());
  // G4cout << position.getX() << " " << position.getY() << " " << position.getZ()  << G4endl;
  
  // // if (radius > 0.215 * m) {
  // G4cout << " radius: " << radius / m << G4endl;
  
  // // record source partcle information
  // // G4ThreeVector position = fParticleSource->GetParticlePosition();
  
  
  


  // }
  // else{
  // G4cout << " GG " << G4endl;
  // // anEvent->SetEventAborted();
  // delete anEvent;
  // // delete fParticleSource;
  //   return;

  // }


  delete[] fVecPrimaryVertex;
  // G4cout << "**************write the position***************" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::InitFunction()
{
  // tabulated function 
  // Y is assumed positive, linear per segment, continuous
  //
  fNPoints = 21;
  const G4double xx[] = 
    { 0.10 * MeV,
      0.15 * MeV,
      0.20 * MeV,
      0.25 * MeV,
      0.30 * MeV,
      0.35 * MeV,
      0.40 * MeV,
      0.45 * MeV,
      0.50 * MeV,
      0.55 * MeV,
      0.60 * MeV,
      0.65 * MeV,
      0.70 * MeV,
      0.75 * MeV,
      0.80 * MeV,
      0.85 * MeV,
      0.90 * MeV,
      0.95 * MeV,
      1.00 * MeV,
      1.05 * MeV,
      1.10 * MeV }; 
      
  const G4double yy[] =
    { 0.0363,
      0.0402,
      0.0440,
      0.0478,
      0.0507,
      0.0532,
      0.0550,
      0.0562,
      0.0568,
      0.0571,
      0.0569,
      0.0562,
      0.0550,
      0.0535,
      0.0514,
      0.0488,
      0.0455,
      0.0416,
      0.0367,
      0.0314,
      0.0257 };
  
  //copy arrays in std::vector and compute fMax
  //
  fX.resize(fNPoints); fY.resize(fNPoints);
  fYmax = 0.;
  for (G4int j=0; j<fNPoints; j++) {
    fX[j] = xx[j]; fY[j] = yy[j];
    if (fYmax < fY[j]) fYmax = fY[j];
  };
     
  //compute slopes
  //
  fSlp.resize(fNPoints);
  for (G4int j=0; j<fNPoints-1; j++) { 
    fSlp[j] = (fY[j+1] - fY[j])/(fX[j+1] - fX[j]);
  };
  
  //compute cumulative function
  //
  fYC.resize(fNPoints);  
  fYC[0] = 0.;
  for (G4int j=1; j<fNPoints; j++) {
    fYC[j] = fYC[j-1] + 0.5*(fY[j] + fY[j-1])*(fX[j] - fX[j-1]);
  };     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double PrimaryGeneratorAction::InverseCumul()
{
  // tabulated function
  // Y is assumed positive, linear per segment, continuous 
  // --> cumulative function is second order polynomial
  
  //choose y randomly
  G4double Yrndm = G4UniformRand()*fYC[fNPoints-1];
  //find bin
  G4int j = fNPoints-2;
  while ((fYC[j] > Yrndm) && (j > 0)) j--;
  //y_rndm --> x_rndm :  fYC(x) is second order polynomial
  G4double Xrndm = fX[j];
  G4double a = fSlp[j];
  if (a != 0.) {
    G4double b = fY[j]/a, c = 2*(Yrndm - fYC[j])/a;
    G4double delta = b*b + c;
    G4int sign = 1; if (a < 0.) sign = -1;
    Xrndm += sign*std::sqrt(delta) - b;    
  } else if (fY[j] > 0.) {
    Xrndm += (Yrndm - fYC[j])/fY[j];
  };
  return Xrndm;
}
