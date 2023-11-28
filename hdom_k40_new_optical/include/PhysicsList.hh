/// \file Physicslist
/// \brief allowed physical process

#ifndef PhysicsList_h
#define PhysicsList_h

#include "globals.hh"
#include "G4VModularPhysicsList.hh"

class PhysicsList : public G4VModularPhysicsList
{
   public:
       PhysicsList();
       ~PhysicsList();
       virtual void SetCuts();
};

#endif

