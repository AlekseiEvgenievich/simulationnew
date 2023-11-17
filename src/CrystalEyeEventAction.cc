// **************************************************************
//
//      ------------ CrystalEyeEventAction  -----------
//           by Libo Wu (23 Sep 2021)
//           Email: libo.wu@gssi.it
//
// **************************************************************

#include "CrystalEyeEventAction.hh"
#include "CrystalEyeAnalysis.hh"
#include "CrystalEyeAnticoincidenceDataManager.hh"
#include "CrystalEyeCalorimeterDataManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"


//------------------------------------------------------------
CrystalEyeEventAction::CrystalEyeEventAction(
		CrystalEyeRunAction*                  aRunAction,
		CrystalEyeAnticoincidenceDataManager* aDataManagerAcd,
		CrystalEyeCalorimeterDataManager*     aDataManagerCal) 
:
	fEventID(-1),
	theRunAction(aRunAction),
	fDataManagerAcd(aDataManagerAcd),
	fDataManagerCal(aDataManagerCal)
{
	;
}

//------------------------------------------------------------
CrystalEyeEventAction::~CrystalEyeEventAction()
{
	;
}

//------------------------------------------------------------
void CrystalEyeEventAction::BeginOfEventAction(const G4Event* evt)
{
	fStruEdep = 0;
	fGammaNum = 0;
	firHitsVolume = 0;
	fInteraction  = 0;
	fEvtInFOV     = 0;
	fFirstPos     = G4ThreeVector(-999,-999,-999);

	fEventID = evt->GetEventID();
	if(fEventID%1000==0)
		G4cout<<"EventID: "<< fEventID << G4endl;

	fDataManagerAcd -> BeginEvent( evt );
	fDataManagerCal -> BeginEvent( evt );

	CrystalEyeAnalysis* fanalysis = CrystalEyeAnalysis::getInstance();
	fanalysis->BeginOfEvent();

}

//------------------------------------------------------------
void CrystalEyeEventAction::EndOfEventAction(const G4Event* evt)
{

	if(!theRunAction){
		G4Exception("CrystalEyeEventAction::BeginOfEventAction()",
				"CRYE0001",FatalException,"Null pointer to Run Action!");
	}

	fDataManagerCal -> FillEvent( evt );
	fDataManagerAcd -> FillEvent( evt );

	CrystalEyeAnalysis* fanalysis = CrystalEyeAnalysis::getInstance();
	fanalysis -> SetEventData(
			fGammaNum,
			firHitsVolume,
			fInteraction,
			fEvtInFOV,
			fStruEdep,
			fFirstPos);

	fanalysis -> EndOfEvent(fEventID);


}

//------------------------------------------------------------
//------------------------------------------------------------
//------------------------------------------------------------
//------------------------------------------------------------
//------------------------------------------------------------
//------------------------------------------------------------
//------------------------------------------------------------
