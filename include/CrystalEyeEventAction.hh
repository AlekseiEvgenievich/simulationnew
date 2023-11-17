// **************************************************************
//
//      ------------ CrystalEyeEventAction  -----------
//           by Libo Wu (23 Sep 2021)
//           Email: libo.wu@gssi.it
//
// **************************************************************

#ifndef CrystalEyeEventAction_h
#define CrystalEyeEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class CrystalEyeRunAction;
class CrystalEyeAnticoincidenceDataManager;
class CrystalEyeCalorimeterDataManager;

class CrystalEyeEventAction : public G4UserEventAction
{
	public:
		CrystalEyeEventAction(
				CrystalEyeRunAction* aRunAction,
				CrystalEyeAnticoincidenceDataManager* aDataManagerAcd,
				CrystalEyeCalorimeterDataManager*     aDataManagerCal);

		virtual ~CrystalEyeEventAction();

		virtual void BeginOfEventAction(const G4Event*);
		virtual void EndOfEventAction(const G4Event*);

		void  AddStructureEdep(G4double edep){ fStruEdep+=edep; }
		void  AddGammaNumber(){	fGammaNum++; }
		void  SetFirstHisVolume(G4int i){		firHitsVolume=i; }
		void  SetInteraction(G4int i){ fInteraction=i; }
		G4int GetInteraction(){ return fInteraction; }
		void  SetEvtInFOV(G4int i){ fEvtInFOV=i; }
		int   GetEvtInFOV(){ return fEvtInFOV;}
		void  SetFirstPos(G4ThreeVector xyz){ fFirstPos = xyz;};

		G4int GetEventID(){ return fEventID; }

	private:
		G4int			fEventID;
		G4double  fStruEdep;
		G4int			fGammaNum;
		G4int		  firHitsVolume;//0/1: Strcuture/Senstive Detector
		G4int     fInteraction;//0/1: No/Yes Interact with Detector
		G4int     fEvtInFOV;
		G4ThreeVector fFirstPos;

		CrystalEyeRunAction* theRunAction;
		CrystalEyeAnticoincidenceDataManager* fDataManagerAcd;
		CrystalEyeCalorimeterDataManager* fDataManagerCal;

};

#endif

