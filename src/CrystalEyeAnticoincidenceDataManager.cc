// ------------------------------------------------------------
//
//      ------------ CrystalEyeAnticoincidenceDataManager.cc  ------
//           by Libo Wu (15 Feb 2022)
//           Email: libo.wu@gssi.it
//
// ------------------------------------------------------------

#include "CrystalEyeAnticoincidenceDataManager.hh"
#include "CrystalEyeDetectorConstruction.hh"
#include "CrystalEyeAnticoincidenceHit.hh"
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
CrystalEyeAnticoincidenceDataManager::CrystalEyeAnticoincidenceDataManager()
	: CollectionID(-1)
{
	G4RunManager* runManager = G4RunManager::GetRunManager();
	fDetector = (CrystalEyeDetectorConstruction*)(runManager->GetUserDetectorConstruction());
	fNofCells = fDetector -> GetNbOfPixel()+1;//UpACDs+DownACD

	fEdep       = new G4double[fNofCells];
	fEdepDelay  = new G4double[fNofCells];
	fOrder      = new G4int[fNofCells];
}

// ------------------------------------------------------------
CrystalEyeAnticoincidenceDataManager::~CrystalEyeAnticoincidenceDataManager()
{
	delete [] fEdep;      fEdep  = NULL;
	delete [] fEdepDelay; fEdepDelay = NULL;
	delete [] fOrder;     fOrder = NULL;
}

// ------------------------------------------------------------
void CrystalEyeAnticoincidenceDataManager::BeginRun(const G4Run* aRun)
{
	;
}

// ------------------------------------------------------------
void CrystalEyeAnticoincidenceDataManager::BeginEvent(const G4Event* aEvt)
{

	fnHits          = 0; 
	fTotalEdep      = 0;
	fTotalEdepDelay = 0;
	for(int i=0; i<fNofCells; i++){
		fEdep[i]      = 0;
		fEdepDelay[i] = 0;
		fOrder[i]     = 0;
	}

}

// ------------------------------------------------------------
void CrystalEyeAnticoincidenceDataManager::FillEvent(const G4Event* aEvt)
{

	G4HCofThisEvent* HCE = aEvt->GetHCofThisEvent();
	CrystalEyeAnticoincidenceHitsCollection* AHC = 0;

	if( HCE ){

		G4SDManager* SDman = G4SDManager::GetSDMpointer();
		if( CollectionID==-1 )
			CollectionID = SDman->GetCollectionID("AnticoincidenceCollection");

		AHC = (CrystalEyeAnticoincidenceHitsCollection*)(HCE->GetHC(CollectionID));
		if( AHC ){

			fnHits = AHC->entries();
			for(int i=0; i<fnHits; i++){
				G4int fID       = (*AHC)[i]->GetACDID();
				fEdep[fID]      = (*AHC)[i]->GetEdep();
				fEdepDelay[fID] = (*AHC)[i]->GetEdepDelay();
				fOrder[fID]     = (*AHC)[i]->GetACDOrder();

				fTotalEdep      = fTotalEdep + fEdep[fID];
				fTotalEdepDelay = fTotalEdepDelay + fEdepDelay[fID];
			}
		}
	}

	CrystalEyeAnalysis* fanalysis = CrystalEyeAnalysis::getInstance();
	fanalysis -> SetDataACD(
			fTotalEdep,
			fnHits,
			fEdep,
			fOrder);
	fanalysis -> SetDataACD_Delay(
			fTotalEdepDelay,
			fEdepDelay
			);

}

// ------------------------------------------------------------
// ------------------------------------------------------------
