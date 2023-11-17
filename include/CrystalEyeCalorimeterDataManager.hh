
// ------------------------------------------------------------
//
//      ------------ CrystalEyeCalorimeterDataManager.hh  ------
//           by Libo Wu (15 Feb 2022)
//           Email: libo.wu@gssi.it
//
// ------------------------------------------------------------

#ifndef CrystalEyeCalorimeterDataManager_h
#define CrystalEyeCalorimeterDataManager_h 1

#include "globals.hh"

class G4Run;
class G4Event;
class CrystalEyeDetectorConstruction;

class CrystalEyeCalorimeterDataManager
{
	public:
		CrystalEyeCalorimeterDataManager();
		~CrystalEyeCalorimeterDataManager();

	public:
		void BeginRun(const G4Run* aRun);
		void BeginEvent(const G4Event* aEvt);
		void FillEvent(const G4Event* aEvt);//EndEvent

	private:
		CrystalEyeDetectorConstruction* fDetector;

		G4int fNofCells;
		G4int CollectionID;

		G4int    fnHits;
		G4int    fNsteps;
		G4double fTotalEdep;
		G4double fTotalEdepDelay;
		G4double *fEdep;
		G4double *fEdepDelay;
		G4int    *fOrder;
	
};


#endif
