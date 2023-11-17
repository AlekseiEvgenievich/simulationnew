// **************************************************************
//
//      ------------ CrystalEyeCalorimeterSD  ------
//           by Libo Wu (20 Sep 2021)
//           Email: libo.wu@gssi.it
//
// **************************************************************

#ifndef CrystalEyeCalorimeterSD_h
#define CrystalEyeCalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"
#include "CrystalEyeCalorimeterHit.hh"

class G4Step;
class G4HCofThisEvent;
class CrystalEyeDetectorConstruction;

class CrystalEyeCalorimeterSD : public G4VSensitiveDetector
{
		public:
				CrystalEyeCalorimeterSD(const G4String name);
				~CrystalEyeCalorimeterSD();

				void Initialize(G4HCofThisEvent*);
				G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* history);
				void EndOfEvent(G4HCofThisEvent* HCE);
				void clear();
				void DrawAll();
				void PrintAll();

		private:
				CrystalEyeCalorimeterHitsCollection* CalorimeterCollection;
				CrystalEyeDetectorConstruction* Detector;

				G4int* fCHitID;
				G4int  fNofCells;

				G4int fhit_order;
				G4int tmp_pixelID;
				G4bool *flag_order;

				G4double fBkgTimeU;
				G4double fBkgTimeD;
				G4double fThrTime;
};

#endif



