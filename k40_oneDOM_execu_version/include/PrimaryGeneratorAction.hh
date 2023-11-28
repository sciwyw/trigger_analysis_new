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
/// \file PrimaryGeneratorAction.hh
/// \brief Definition of the PrimaryGeneratorAction class

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include <vector>
#include "G4ParticleGun.hh"
#include <iostream>
#include <fstream>
#include "yaml-cpp/yaml.h"
class G4ParticleGun;
class G4Event;
class G4GeneralParticleSource;
/// The primary generator action class with particle gum.
///
/// It defines a single particle which hits the Tracker 
/// perpendicular to the input face. The type of the particle
/// can be changed via the G4 build-in commands of G4ParticleGun class 
/// (see the macros provided with this example).

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();    
    virtual ~PrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event* );

    G4GeneralParticleSource* GetParticleSource() {return fParticleSource;}

  public:        
    //G4double RejectAccept();
    G4double InverseCumul();
    G4double decaytime;

  private:
    std::ofstream outfile;
    G4GeneralParticleSource*          fParticleSource; // G4 ParticleSource
    G4ParticleGun*         fParticleGun;
    G4int                  fNPoints;    //nb of points
    std::vector<G4double>  fX;          //abscisses X
    std::vector<G4double>  fY;          //values of Y(X)
    std::vector<G4double>  fSlp;        //slopes
    std::vector<G4double>  fYC;         //cumulative function of Y
    G4double               fYmax;       //max(Y)
  private:
    void InitFunction();      
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
