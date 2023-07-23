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
/// \file TrackerHit.hh
/// \brief Definition of the TrackerHit class

#ifndef TrackerHit_h
#define TrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

/// Tracker hit class
///
/// It defines data members to store the trackID, chamberNb, energy deposit,
/// and position of charged particles in a selected volume:
/// - fTrackID, fChamberNB, fEdep, fPos

class TrackerHit : public G4VHit
{
  public:
    TrackerHit();
    TrackerHit(const TrackerHit&);
    virtual ~TrackerHit();

    // operators
    const TrackerHit& operator=(const TrackerHit&);
    G4bool operator==(const TrackerHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    //virtual void Draw();
    virtual void Print();

    // Set methods
    //void SetTrackID  (G4int track)      { fTrackID = track; };
    //void SetChamberNb(G4int chamb)      { fChamberNb = chamb; };
    //void SetEdep     (G4double de)      { fEdep = de; };
    //void SetPos      (G4ThreeVector xyz){ fPos = xyz; };

    void SetTime (G4double time){fTime =time;};
    void SetcpN  (G4int cpN){cp_number=cpN;};
    //void SetcpN_theta  (G4int cpN){cp_number_theta=cpN;};
    void SetEnergy (G4double energy){fEnergy = energy;};

    // Get methods
    //G4int GetTrackID() const     { return fTrackID; };
    //G4int GetChamberNb() const   { return fChamberNb; };
    //G4double GetEdep() const     { return fEdep; };
    //G4ThreeVector GetPos() const { return fPos; };
    
    G4double GetTime()const {return fTime;};
    G4int GetcpN()const {return cp_number;};
    //G4int GetcpN_theta()const {return cp_number_theta;};
    G4double GetEnergy ()const {return fEnergy;};
  private:
      /*
      G4int         fTrackID;
      G4int         fChamberNb;
      G4double      fEdep;
      G4ThreeVector fPos;
      */
     G4double      fEnergy;
     G4double fTime;
     G4int cp_number;
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<TrackerHit> TrackerHitsCollection;

extern G4ThreadLocal G4Allocator<TrackerHit>* TrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* TrackerHit::operator new(size_t)
{
  if(!TrackerHitAllocator)
      TrackerHitAllocator = new G4Allocator<TrackerHit>;
  return (void *) TrackerHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void TrackerHit::operator delete(void *hit)
{
  TrackerHitAllocator->FreeSingle((TrackerHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
