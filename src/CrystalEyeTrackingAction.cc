// **************************************************************
//
//      ------ CrystalEyeTrackingAction  ------
//           by Libo Wu (23 Oct 2021)
//           Email: libo.wu@gssi.it
//
// **************************************************************

#include "CrystalEyeTrackingAction.hh"
#include "CrystalEyeEventAction.hh"
#include "G4Track.hh"


//------------------------------------------------------------
CrystalEyeTrackingAction::CrystalEyeTrackingAction(CrystalEyeEventAction* evt)
	:G4UserTrackingAction()
{
	eventAct = evt;
}

//------------------------------------------------------------
CrystalEyeTrackingAction::~CrystalEyeTrackingAction()
{
}

//------------------------------------------------------------
void CrystalEyeTrackingAction::PreUserTrackingAction(const G4Track* track)
{

	G4int ftrackid = track->GetTrackID();
	G4int fpdg = track->GetDefinition()->GetPDGEncoding();
	//G4String fparticleName = track->GetDefinition()->GetParticleName();

	if(fpdg==22)		eventAct->AddGammaNumber();

	/*
		 G4cout
		 <<"trackID: "<<ftrackid<<"  "
		 <<"PDGCode: "<<fpdg<<"  "
	//<<"Particle: "<<fparticleName.c_str()<<"  "
	<<G4endl;
	 */
}

//------------------------------------------------------------
void CrystalEyeTrackingAction::PostUserTrackingAction(const G4Track* track)
{
}

//------------------------------------------------------------
//------------------------------------------------------------
//------------------------------------------------------------
