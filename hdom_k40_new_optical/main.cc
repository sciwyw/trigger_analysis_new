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
/// \file main.cc
/// \brief Main program of the coincidence event

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

// #include "G4RunManagerFactory.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "G4StepLimiterPhysics.hh"

#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
// #include "MyAnalysis.hh"
#include "PhysicsList.hh"
#include "time.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Detect interactive mode (if no arguments) and define UI session
  //
  //char* fileG4 = argv[1];
  //char* fileGeometry= argv[2];
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
  G4long seed = time(NULL);
  CLHEP::HepRandom::setTheSeed(seed);
  // G4UIExecutive* ui = 0;
  if (argc != 3) {
    std::cerr << "Wrong input arguments! The right way is: ./main [config/config.yaml] [config/optical_property.yaml]";
    return -1;
  }
  // } else {
  //   //fileG4 = argv[1];
  //   //fileGeometry = argv2;
  //   ui = new G4UIExecutive(2, &fileG4);
  // }
  
 

  // Optionally: choose a different Random engine...
  // G4Random::setTheEngine(new CLHEP::MTwistEngine);

  // Construct the default run manager
  //
  // auto* runManager =
  //   G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);
  //runManager->SetNumberOfThreads((G4Threading::G4GetNumberOfCores()));

  G4RunManager* runManager = new G4RunManager;
  // Set mandatory initialization classes
  //
  runManager->SetUserInitialization(new DetectorConstruction(argv[1],argv[2]));

  runManager->SetUserInitialization(new PhysicsList());
    
  // Set user action classes
  runManager->SetUserInitialization(new ActionInitialization());


 // Initialize visualization
  //
  //G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  //visManager->Initialize();
  // Process macro or start UI session
  //
  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand("/tracking/verbose 0");
  UImanager->ApplyCommand("/run/verbose      0");
  UImanager->ApplyCommand("/control/verbose  0");
  UImanager->ApplyCommand("/analysis/verbose 0");
  UImanager->ApplyCommand("/event/verbose    0");
  UImanager->ApplyCommand("/process/verbose  0");
  runManager->Initialize(); 
  runManager->BeamOn(20000000);//(250000000)
  //runManager->/gps/particle gamma
  ///gps/energy 1.459 MeV

  ///run/beamOn 29950000

  // Process macro or start UI session
  //
  // std::string command = "/control/execute ";
  // std::string fileName = argv[1];
  // UImanager->ApplyCommand(command+fileName);
  // if ( ! ui ) {
  //   // batch mode
  //   std::string command = "/control/execute ";
  //   std::string fileName = argv[1];
  //   UImanager->ApplyCommand(command+fileName);
  // }
  // else {
  //   // interactive mode
    
  //   ui->SessionStart();
  //   delete ui;
  // }
  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !
  //
  //delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
