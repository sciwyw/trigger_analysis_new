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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1
#include "G4AutoDelete.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "tls.hh"
#include "yaml-cpp/yaml.h"
#include <vector>
#include "G4Sphere.hh"
#include "G4Box.hh"
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4UserLimits;

/// Detector construction class to define materials, geometry


class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction(std::string fileGeometry);
    virtual ~DetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();
    void SetCheckOverlaps(G4bool);
    
    YAML::Node rootNode;
    void SetMaxStep(G4double);
  private:
    // methods
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();
    // data members
    G4Material*        PMTandSiPM_Material;  
    G4LogicalVolume**   fLogicPMT;           // pointer to the logical PMT
    
    //G4bool  fCheckOverlaps;                  // option to activate checking of volumes overlaps 
    //G4VPhysicalVolume *DefineVolumes();
  // define the solid and logic vloume of 
    G4Sphere *solidPmt;
    G4LogicalVolume *logicPmt;
    G4Material*        Seawater_Material;    // pointer to the seawater material
    G4Material*        Gel_Material;         // optical glue
    G4Material*        Glass_Material;       // glass ball
    G4Material*        Pmt_Material;         // pmt
    //G4Material*        Plastic_Material;     // plastic suppoter
    G4Material*        plastic_Material;
    // G4Material*        PMTandSiPM_Material;            
    
    G4UserLimits *fStepLimit;  // pointer to user step limits
    G4bool fCheckOverlaps;     // option to activate checking of volumes overlaps
    G4double Shell_R;
    std::string fileGeometry;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
