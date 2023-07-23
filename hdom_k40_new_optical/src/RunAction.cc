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
/// \file RunAction.cc
/// \brief Implementation of the RunAction class
#include "g4root.hh"
#include "RunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
// #include "MyAnalysis.hh"
#include "Randomize.hh"
// s#include "G4AnalysisManager.hh"
// #include "G4VAnalysisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
 : G4UserRunAction()
{ 
  // set printing event number per each 100 events
  G4RunManager::GetRunManager()->SetPrintProgress(50000);
  
  //auto analysisManager = G4AnalysisManager::Instance();
  
  // //make a new folder
  // folderPath = std::to_string(G4UniformRand());
  // G4String command = "mkdir -p " + folderPath;  
  
  // int systemRet = system(command.c_str());
  // if(systemRet == -1){
  //   G4cout << "make a new folder : failed" << G4endl;
  // }
  // }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  //inform the runManager to save random number seed
  // G4cout << "++++ Initialize Analysis " << G4endl;
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  
  
  auto *analysisManager = G4AnalysisManager::Instance();
  G4int runid =  aRun->GetRunID();
  G4String run_num = std::to_string(runid);
  // G4String OpenPath = folderPath + "/data.csv";
  G4String OpenPath =  std::to_string(G4UniformRand()) + "_data.root";
  analysisManager->OpenFile(OpenPath);
  analysisManager->SetVerboseLevel( 1 );
  analysisManager->CreateNtuple(run_num ,"time");
  
  analysisManager->CreateNtupleIColumn("eventID");
  analysisManager->CreateNtupleDColumn("t(ns)");
  analysisManager->CreateNtupleDColumn("copynumber");
  //analysisManager->CreateNtupleDColumn("copynumber_phi");
  analysisManager->CreateNtupleDColumn("energy(eV)");
  analysisManager->CreateNtupleDColumn("x(mm)");
  analysisManager->CreateNtupleDColumn("y(mm)");
  analysisManager->CreateNtupleDColumn("z(mm)");
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* )
{ 
  // G4cout << "*****************end of run action ******************" << G4endl;
  auto *analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  //analysisManager->CloseFile();

  delete analysisManager;
   //analysisManager->Clear();
  //  delete analysisManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
