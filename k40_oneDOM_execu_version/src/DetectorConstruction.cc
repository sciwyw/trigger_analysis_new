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
using std::vector;
DetectorConstruction::DetectorConstruction(std::string fileName,std::string fileName_1)
    : G4VUserDetectorConstruction(),
      solidPmt(NULL), logicPmt(NULL),
      solidSipm(NULL), logicSipm(NULL), fLogichDom(NULL),
      Seawater_Material(NULL), Gel_Material(NULL), Glass_Material(NULL),
      plastic_Material(NULL), PMTandSiPM_Material(NULL), fStepLimit(NULL),
      fCheckOverlaps(true), fileGeometry(fileName),fileOptical(fileName_1)
{
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  //delete fLogichDom;
  delete solidPmt;
  delete solidSipm;
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
  
  vector<double> energy_;
  vector<double> refracIdx_;
  vector<double> absLen_;
  vector<double> scaLenRay_;
  vector<double> scaLenMie_;
  rootNode_water_property = YAML::LoadFile(fileOptical);
  energy_ = rootNode_water_property["energy"].as<vector<double>>();
  geoOptical.num = energy_.size();
  geoOptical.energy = new double[geoOptical.num];
  // Set the refractiveIndex, From Bjorn Herold's PhD thesis Equation A 1
  geoOptical.energy = new double[geoOptical.num];
  for (int i = 0; i < geoOptical.num; i++)
    geoOptical.energy[i] = energy_[i] * eV;

  refracIdx_ = rootNode_water_property["refractive_index"].as<vector<double>>();
  geoOptical.refracIdx = new double[geoOptical.num];
  for (int i = 0; i < geoOptical.num; i++)
    geoOptical.refracIdx[i] = refracIdx_[i];

  absLen_ = rootNode_water_property["absorption_length"].as<vector<double>>();
  geoOptical.absLen = new double[geoOptical.num];
  for (int i = 0; i < geoOptical.num; i++)
    geoOptical.absLen[i] = absLen_[i] * m;

  scaLenRay_ = rootNode_water_property["scatter_length_rayeigh"].as<vector<double>>();
  geoOptical.scaLenRay = new double[geoOptical.num];
  for (int i = 0; i < geoOptical.num; i++)
    geoOptical.scaLenRay[i] = scaLenRay_[i] * m;

  scaLenMie_ = rootNode_water_property["scatter_length_mie"].as<vector<double>>();
  geoOptical.scaLenMie = new double[geoOptical.num];
  for (int i = 0; i < geoOptical.num; i++)
    geoOptical.scaLenMie[i] = scaLenMie_[i] * m;

  geoOptical.mieForward = rootNode_water_property["mie_forward_angle"].as<double>();

  myMPT1->AddProperty("RINDEX", geoOptical.energy, geoOptical.refracIdx, geoOptical.num);

  myMPT1->AddProperty("ABSLENGTH", geoOptical.energy, geoOptical.absLen, geoOptical.num);

  myMPT1->AddProperty("RAYLEIGH", geoOptical.energy, geoOptical.scaLenRay, geoOptical.num);

  myMPT1->AddProperty("MIEHG", geoOptical.energy, geoOptical.scaLenMie, geoOptical.num);

  myMPT1->AddConstProperty("MIEHG_FORWARD",geoOptical.mieForward);
  myMPT1->AddConstProperty("MIEHG_BACKWARD",geoOptical.mieBackward);
  myMPT1->AddConstProperty("MIEHG_FORWARD_RATIO",geoOptical.mieRatio);

  Seawater_Material->SetMaterialPropertiesTable(myMPT1);
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  plastic_Material = new G4Material("plastic", 1.19 * g / cm3, 3);
  plastic_Material->AddElement(elC, 5);
  plastic_Material->AddElement(elH, 8);
  plastic_Material->AddElement(elO, 2);
  plastic_Material->GetIonisation()->SetBirksConstant(0.01 * m);
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
  Gel_Material = new G4Material("Gel", 1.20 * g / cm3, 3);
  Gel_Material->AddElement(elC, 4);
  Gel_Material->AddElement(elH, 8);
  Gel_Material->AddElement(elO, 2);
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::DefineVolumes()
{
  SetPSVolumes();
  //world size
  G4double worldlength = 30 * m;
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
  G4Sphere *solidDomGlass =
      new G4Sphere("DOM", 0, 21.5 * cm, 0.f, 2 * M_PI, 0.f, M_PI);
  // G4Sphere *solidGel = new G4Sphere("Gel",
  //                                   0,
  //                                   20.2 * cm,
  //                                   0,
  //                                   2 * M_PI,
  //                                   0,
  //                                   M_PI);
  G4Sphere *solidSupportor =
      new G4Sphere("Supportor", 0, 16.4 * cm, 0, 2 * M_PI, 0, M_PI);

  G4LogicalVolume *logicDom =
      new G4LogicalVolume(solidDomGlass, Glass_Material, "DOM");
  // G4LogicalVolume *logicGel = new G4LogicalVolume(solidGel, Gel_Material,
  // "Gel");
  G4LogicalVolume *logicSupportor =
      new G4LogicalVolume(solidSupportor, plastic_Material, "Supportor");
  // G4OpticalSurface *opGlassSurface = new G4OpticalSurface("DOM");
  // opGlassSurface->SetType(dielectric_dielectric);
  // opGlassSurface->SetFinish(polished);
  // opGlassSurface->SetModel(glisur);
  // new G4LogicalSkinSurface("DOM",
  //                          logicDom,
  //                          opGlassSurface);

  // G4OpticalSurface *opGelSurface = new G4OpticalSurface("Gel");
  // opGelSurface->SetType(dielectric_dielectric);
  // opGelSurface->SetFinish(polished);
  // opGelSurface->SetModel(glisur);
  // new G4LogicalSkinSurface("Gel",
  //                          logicGel,
  //                          opGelSurface);
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicDom, "DOM", worldLV, false,
                    0, false);
  // new G4PVPlacement(0,
  //                   G4ThreeVector(0, 0, 0),
  //                   logicGel,
  //                   "Gel",
  //                   logicDom,
  //                   false,
  //                   0,
  //                   true);
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicSupportor, "Supportor",
                    logicDom, false, 0, true);
  rootNode = YAML::LoadFile(fileGeometry);
  // read the yaml file
  // build PMT
  numOftheta0 = rootNode["numOftheta0"].as<int>();
  for (G4int i = 0; i < numOftheta0; i++) { 
    theta_array0.push_back(rootNode["theta_array0"][i].as<double>());
    phi_array0.push_back(rootNode["phi_array0"][i].as<std::vector<double>>());
  }
  // get the theta and phi parameter
  std::vector<G4Transform3D> transform0;
  // set the transform matrix
  G4Transform3D rot1, rot2;
  G4Transform3D tran_PMT = G4Translate3D(G4ThreeVector(0, 0, 149. * mm));
  for (int i = 0; i < numOftheta0; i++) {
    for (int j = 0; j < int(phi_array0[i].size()); j++) {

      rot1 = G4RotateX3D(theta_array0[i] / 180. * M_PI);
      rot2 = G4RotateZ3D(phi_array0[i][j] / 180. * M_PI);
      transform0.push_back(rot2 * rot1 * tran_PMT);

    }
  }

  // transform0.push_back(G4RotateX3D(M_PI) * tran_PMT);
  for (uint i = 0; i < transform0.size(); i++) {
    new G4PVPlacement(transform0[i], logicPmt, "PMT", logicDom, false, i, true);
  }



  // build Sipm same as PMTs
  numOftheta1 = rootNode["numOftheta1"].as<int>();
  for (G4int i = 0; i < numOftheta1; i++) {
    theta_array1.push_back(rootNode["theta_array1"][i].as<double>());
    phi_array1.push_back(rootNode["phi_array1"][i].as<std::vector<double>>());
  }
  std::vector<G4Transform3D> transform1;
  G4Transform3D rot3, rot4;
  G4Transform3D tranSipm = G4Translate3D(G4ThreeVector(0, 0, 180 * mm));
  for (int i = 0; i < numOftheta1; i++) {
    for (int j = 0; j < int(phi_array1[i].size()); j++) {
      rot3 = G4RotateX3D(theta_array1[i] / 180. * M_PI);
      rot4 = G4RotateZ3D(phi_array1[i][j] / 180. * M_PI);
      transform1.push_back(rot4 * rot3 * tranSipm);
    }
  }
  // transform1.push_back(G4RotateX3D(M_PI) * tranSipm);
  for (uint i = 0; i < transform1.size(); i++) {
    new G4PVPlacement(transform1[i], logicSipm, "SiPM", logicDom, false, i,
                      true);
  }
  // Always return the physical world
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetPSVolumes() {
  solidPmt =
      new G4Sphere("PMT", 52 * mm, 53 * mm, 0, 2 * M_PI, 0, 45. / 180 * M_PI);
  logicPmt = new G4LogicalVolume(solidPmt, PMTandSiPM_Material, "PMT");
  G4double chip_size_sipm = 2.7 * cm;
  solidSipm = new G4Box("SiPM", chip_size_sipm / 2, chip_size_sipm / 2,
                        chip_size_sipm / 2);
  logicSipm = new G4LogicalVolume(solidSipm, PMTandSiPM_Material, "SiPMs");
}

void DetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors
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

