// **************************************************************
//
//      ------------ CrystalEyeEventAction  -----------
//           by Libo Wu (23 Sep 2021)
//           Email: libo.wu@gssi.it
//
// **************************************************************

#include "CrystalEyeRunAction.hh"
#include "CrystalEyeAnalysis.hh"
#include "CrystalEyeAnticoincidenceDataManager.hh"
#include "CrystalEyeCalorimeterDataManager.hh"

#include "G4Run.hh"
#include "G4UImanager.hh" 
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4Threading.hh"


//------------------------------------------------------------
CrystalEyeRunAction::CrystalEyeRunAction(
		CrystalEyeAnticoincidenceDataManager* aDataManagerAcd,
		CrystalEyeCalorimeterDataManager*     aDataManagerCal) 
	: 
		fRunID(-1),
		fDataManagerAcd(aDataManagerAcd),
		fDataManagerCal(aDataManagerCal)
{
		;
}

//------------------------------------------------------------
CrystalEyeRunAction::~CrystalEyeRunAction()
{
		;
}

//------------------------------------------------------------
void CrystalEyeRunAction::BeginOfRunAction(const G4Run* aRun)
{
		fRunID = aRun->GetRunID();
		fDataManagerAcd->BeginRun( aRun );
		fDataManagerCal->BeginRun( aRun );

		if (IsMaster())
				G4cout << "### Run " << aRun->GetRunID() << " starts (master)." << G4endl;
		else 
				G4cout << "### Run " << aRun->GetRunID() << " starts (worker)." << G4endl;

		if (G4VVisManager::GetConcreteInstance()){
				G4UImanager* UI = G4UImanager::GetUIpointer();
				UI->ApplyCommand("/vis/scene/notifyHandlers");
		}

		CrystalEyeAnalysis* analysis = CrystalEyeAnalysis::getInstance();
		analysis->BeginOfRun();

}

//------------------------------------------------------------
void CrystalEyeRunAction::EndOfRunAction(const G4Run* aRun)
{
		G4cout << "End of Run "<< aRun->GetRunID() << G4endl;

		CrystalEyeAnalysis* analysis = CrystalEyeAnalysis::getInstance();
		analysis->EndOfRun();

}

//------------------------------------------------------------
//------------------------------------------------------------
//------------------------------------------------------------




