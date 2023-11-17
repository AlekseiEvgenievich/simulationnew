// ------------------------------------------------------------
//
//      ------------ CrystalEyeCalorimeterDataManager.cc  ------
//           by Libo Wu (15 Feb 2022)
//           Email: libo.wu@gssi.it
//
// ------------------------------------------------------------

#include "CrystalEyeCalorimeterDataManager.hh"
#include "CrystalEyeDetectorConstruction.hh"
#include "CrystalEyeCalorimeterHit.hh"
#include "CrystalEyeAnalysis.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"


// ------------------------------------------------------------
CrystalEyeCalorimeterDataManager::CrystalEyeCalorimeterDataManager()
	: CollectionID(-1)
{
	G4RunManager* runManager = G4RunManager::GetRunManager();
	fDetector = (CrystalEyeDetectorConstruction*)(runManager->GetUserDetectorConstruction());
	fNofCells = 2 * ( fDetector -> GetNbOfPixel() );

	fEdep       = new G4double[fNofCells];
	fEdepDelay  = new G4double[fNofCells];
	fOrder      = new G4int[fNofCells];
}

// ------------------------------------------------------------
CrystalEyeCalorimeterDataManager::~CrystalEyeCalorimeterDataManager()
{
	delete [] fEdep;      fEdep  = NULL;
	delete [] fEdepDelay; fEdepDelay = NULL;
	delete [] fOrder;     fOrder = NULL;
}

// ------------------------------------------------------------
void CrystalEyeCalorimeterDataManager::BeginRun(const G4Run* aRun)
{
	;
}

// ------------------------------------------------------------
void CrystalEyeCalorimeterDataManager::BeginEvent(const G4Event* aEvt)
{

	fnHits          = 0; 
	fNsteps         = 0;
	fTotalEdep      = 0;
	fTotalEdepDelay = 0;
	

	for(int i=0; i<fNofCells; i++){
		fEdep[i]      = 0;
		fEdepDelay[i] = 0;
		fOrder[i]     = 0;
	}

}

// ------------------------------------------------------------
void CrystalEyeCalorimeterDataManager::FillEvent(const G4Event* aEvt)
{

	G4HCofThisEvent* HCE = aEvt->GetHCofThisEvent();
	CrystalEyeCalorimeterHitsCollection* CHC = 0;

	if( HCE ){

		G4SDManager* SDman = G4SDManager::GetSDMpointer();
		if( CollectionID==-1 )
			CollectionID = SDman->GetCollectionID("CalorimeterCollection");

		CHC = (CrystalEyeCalorimeterHitsCollection*)(HCE->GetHC(CollectionID));
		if( CHC ){

			fnHits = CHC->entries();
			for(int i=0; i<fnHits; i++){
				G4int fID        = (*CHC)[i]->GetCrystalID();
				fEdep[fID]       = (*CHC)[i]->GetEdep();
				fEdepDelay[fID]  = (*CHC)[i]->GetEdepDelay();
				fOrder[fID]      = (*CHC)[i]->GetCrystalOrder();
				fNsteps          = (*CHC)[i]->GetnSteps();

				fTotalEdep      = fTotalEdep + fEdep[fID];
				fTotalEdepDelay = fTotalEdepDelay + fEdepDelay[fID];
			}
		}
	}

	CrystalEyeAnalysis* fanalysis = CrystalEyeAnalysis::getInstance();
	fanalysis -> SetDataCAL(
			fTotalEdep,
			fnHits,
			fEdep,
			fOrder);
	fanalysis -> SetDataCAL_Delay(
			fTotalEdepDelay,
			fEdep
			);

}

// ------------------------------------------------------------
// ------------------------------------------------------------
