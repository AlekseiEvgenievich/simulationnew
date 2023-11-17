// **************************************************************
//
//      ------ CrystalEyePhysicsList.hh  ------
//           Modified by Libo Wu (23 Sep 2021)
//           Email: libo.wu@gssi.it
//
// **************************************************************

#ifndef CrystalEyePhysicsList_h
#define CrystalEyePhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "G4EmConfigurator.hh"
#include "globals.hh"

class G4VPhysicsConstructor;
class CrystalEyePhysicsListMessenger;

class CrystalEyePhysicsList: public G4VModularPhysicsList
{
public:

  CrystalEyePhysicsList();
  virtual ~CrystalEyePhysicsList();

  void ConstructParticle();

  void SetCuts();
  void SetCutForGamma(G4double);
  void SetCutForElectron(G4double);
  void SetCutForPositron(G4double);

  void AddPhysicsList(const G4String& name);
  void ConstructProcess();

  void AddPackage(const G4String& name);

private:

  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;

  G4EmConfigurator em_config;
  G4VPhysicsConstructor*               fDecaPhys;
	std::vector<G4VPhysicsConstructor*>  fElecPhys;
  std::vector<G4VPhysicsConstructor*>  fHadrPhys;

  CrystalEyePhysicsListMessenger* pMessenger;
};

#endif
