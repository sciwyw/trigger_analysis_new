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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "TrackerSD.hh"
#include "G4Material.hh"
#include "G4SDManager.hh"
#include "G4Sphere.hh"
#include "G4Box.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "yaml-cpp/yaml.h"
#include <vector>
#include "G4GeometryManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4UserLimits.hh"
#include "G4Orb.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction(std::string fileName)
    : G4VUserDetectorConstruction(),
      solidPmt(NULL), logicPmt(NULL),
      Seawater_Material(NULL),
      fStepLimit(NULL),
      fCheckOverlaps(true), fileGeometry(fileName)
{
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  //delete fLogichDom;
  delete solidPmt;
  delete fStepLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // Material definition

  //define the element of seawater
  G4double a; // Zeff
  a = 1.01 * g / mole;
  G4Element *elH = new G4Element("Hydrogen", "H", 1., a);
  a = 12.01 * g / mole;
  G4Element *elC = new G4Element("Carbon", "C", 6., a);
  a = 16.00 * g / mole;
  G4Element *elO = new G4Element("Oxygen", "O", 8., a);
  a = 28.00 * g / mole;
  G4Element *elSi = new G4Element("Silicon", "Si", 14., a);
  a = 22.99 * g / mole;
  G4Element *elNa = new G4Element("Sodium", "Na", 11, a);
  a = 35.453 * g / mole;
  G4Element *elCl = new G4Element("Chlorine", "Cl", 17., a);
  G4Material *NaCl = new G4Material("Sodium Chlorure", 2.16 * g / cm3, 2);

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
  NaCl->AddElement(elNa, 1);
  NaCl->AddElement(elCl, 1);

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  G4Material *H2O = new G4Material("Water", 1.000 * g / cm3, 2);
  H2O->AddElement(elH, 2);
  H2O->AddElement(elO, 1);

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
  //define the seawater property: density,state,press,temperture

  Seawater_Material = new G4Material("SeaWater", 1.04 * g / cm3, 2, kStateLiquid, 300. * atmosphere, 275. * kelvin);
  Seawater_Material->AddMaterial(NaCl, 3.5 * perCent);
  Seawater_Material->AddMaterial(H2O, 96.5 * perCent);

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  G4MaterialPropertiesTable *myMPT1 = new G4MaterialPropertiesTable();
  
  // Three wavelength LEDs in TRIDENT Pathfinder
  G4double wl[3] = { 525 * nm, 450 * nm, 405 * nm };
  G4double photonEnergy[3];

  for (int i = 0; i < 3; i++)
      photonEnergy[i] = 1240 * nm * eV / wl[i];
  
  // Set the refractiveIndex, From Bjorn Herold's PhD thesis Equation A 1
  G4double refractiveIndex1[3];
  G4double preasure = 3400; // in unit of bar
  for (int i = 0; i < 3; i++)
  {
      refractiveIndex1[i] = 1.3201 + 1.4e-5 * preasure + 16.2566 * std::pow(wl[i] / nm, -1) - 4383.0 * std::pow(wl[i] / nm, -2) + 1.1455e6 * std::pow(wl[i] / nm, -3);
  }
  myMPT1->AddProperty("RINDEX", photonEnergy, refractiveIndex1, 3);

  // Set ABSLENGTH from TRIDENT Pathfinder paper: https://arxiv.org/abs/2207.04519 
  G4double absorb[3] = {18.5 * m, 26.4 * m, 19.2 * m};
  myMPT1->AddProperty("ABSLENGTH", photonEnergy, absorb, 3);

  // Set RAYLEIGH scatter length 
  G4double Rey_scatter[3] = {300. * m, 203 * m, 114 * m};
  myMPT1->AddProperty("RAYLEIGH", photonEnergy, Rey_scatter, 3);

  // Set MIE scatter length 
  G4double Mie_scatter[3] = {127 * m, 64 * m, 46 * m};
  myMPT1->AddProperty("MIEHG", photonEnergy, Mie_scatter, 3);

  // Set MIE angle
  G4double mie_water_const[3]={0.98,0.98,0.8};

  myMPT1->AddConstProperty("MIEHG_FORWARD",mie_water_const[0]);
  myMPT1->AddConstProperty("MIEHG_BACKWARD",mie_water_const[1]);
  myMPT1->AddConstProperty("MIEHG_FORWARD_RATIO",mie_water_const[2]);

  Seawater_Material->SetMaterialPropertiesTable(myMPT1);
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  plastic_Material = new G4Material("plastic", 1.19 * g / cm3, 3);
  plastic_Material->AddElement(elC, 5);
  plastic_Material->AddElement(elH, 8);
  plastic_Material->AddElement(elO, 2);
  plastic_Material->GetIonisation()->SetBirksConstant(0.01 * m);
  // Plastic_Material->GetIonisation()->s
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
  Gel_Material = new G4Material("Gel", 1.20 * g / cm3, 3);
  Gel_Material->AddElement(elC, 4);
  Gel_Material->AddElement(elH, 8);
  Gel_Material->AddElement(elO, 2);
  // Gel_Material->GetIonisation()->
  // build optical property
  G4MaterialPropertiesTable *mpt = new G4MaterialPropertiesTable();
  
  const G4int num = 2;
  G4double photon_Energy1[num] = {2.06667 * eV, 4.13333 * eV};
  G4double refractive_Index1[num] = {1.41, 1.41}; // SilGel 601 A/B by Wacker
  G4double absorption_Length1[num] = {10. * m, 10. * m};
  mpt->AddProperty("RINDEX", photon_Energy1, refractive_Index1, num);
  mpt->AddProperty("ABSLENGTH", photon_Energy1, absorption_Length1, num);
  Gel_Material->SetMaterialPropertiesTable(mpt);
  
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
  // build Glass material component
  Glass_Material = new G4Material("Glass", 1.19 * g / cm3, 2);
  Glass_Material->AddElement(elSi, 1);
  Glass_Material->AddElement(elO, 2);

  // build optical property
  G4MaterialPropertiesTable *mpt2 = new G4MaterialPropertiesTable();
  
  G4double photon_Energy2[num] = {2.06667 * eV, 4.13333 * eV};
  G4double refractive_Index2[num] = {1.50, 1.50};
  G4double absorption_Length2[num] = {10. * m, 10. * m};
  mpt2->AddProperty("RINDEX", photon_Energy2, refractive_Index2, num);
  mpt2->AddProperty("ABSLENGTH", photon_Energy2, absorption_Length2, num);
  Glass_Material->SetMaterialPropertiesTable(mpt2);

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  PMTandSiPM_Material = new G4Material("PMTandSiPM", 1.20 * g / cm3, 1);
  PMTandSiPM_Material->AddElement(elSi, 1);

  // build optical property
  G4MaterialPropertiesTable *mpt_PS = new G4MaterialPropertiesTable();
  const G4int num3 = 2;
  G4double photon_Energy_PS[num3] = {2.06667 * eV, 4.13333 * eV};
  G4double refractive_Index_PS[num3] = {1.34, 1.34};
  mpt_PS->AddProperty("RINDEX", photon_Energy_PS, refractive_Index_PS, num3);
  PMTandSiPM_Material->SetMaterialPropertiesTable(mpt2);

  //build pmt material component
  Pmt_Material = new G4Material("Pmt_Material", 1.19 * g / cm3, 2);
  Pmt_Material ->AddElement(elSi, 1);
  Pmt_Material ->AddElement(elO, 2);

  G4MaterialPropertiesTable *mpt3 = new G4MaterialPropertiesTable();
  
  G4double photon_Energy3[num] = {2.06667 * eV, 4.13333 * eV};
  G4double refractive_Index3[num] = {1.50, 1.50};
  G4double absorption_Length3[num] = {0.1 * mm, 0.1 * mm};
  mpt3->AddProperty("RINDEX", photon_Energy3, refractive_Index3, num);
  mpt3->AddProperty("ABSLENGTH", photon_Energy3, absorption_Length3, num);
  Pmt_Material->SetMaterialPropertiesTable(mpt3);
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
  // Print materials
  //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::DefineVolumes()
{
  
  rootNode = YAML::LoadFile(fileGeometry);
  Shell_R = rootNode["shell_radius"].as<double>();
  //world size
  G4double worldlength = (Shell_R + 1)* m;
  G4Orb* worldS = new G4Orb("world", worldlength);//, worldlength, worldlength);
  G4LogicalVolume *worldLV = new G4LogicalVolume(
                                                worldS,            //its solid
                                                Seawater_Material, //its material
                                                "World");          //its name
  G4VPhysicalVolume *worldPV = new G4PVPlacement(
                                                  0,               // no rotation
                                                  G4ThreeVector(), // at (0,0,0)
                                                  worldLV,         // its logical volume
                                                  "World",         // its name
                                                  0,               // its mother  volume
                                                  false,           // no boolean operations
                                                  0);              // checking overlaps
  solidPmt =
      new G4Sphere("PMT", Shell_R * m, (Shell_R + 0.1) * m, 0, 2 * M_PI, 0, M_PI);
  logicPmt = new G4LogicalVolume(solidPmt, PMTandSiPM_Material, "PMT");
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicPmt, "PMT",
                    worldLV, false, 0, true);
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors
  // G4cout << "Set sensitive det" << G4endl;
  G4String trackerChamberSDname = "PMT";
  TrackerSD *aTrackerSD = new TrackerSD(trackerChamberSDname,
                                            "TrackerHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  // Setting aTrackerSD to all logical volumes with the same name
  SetSensitiveDetector("PMT", aTrackerSD, true);
  //SetSensitiveDetector("SiPMs", aTrackerSD, true);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetMaxStep(G4double maxStep) {
  if ((fStepLimit) && (maxStep > 0.))
    fStepLimit->SetMaxAllowedStep(maxStep);
}
void DetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
  fCheckOverlaps = checkOverlaps;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

