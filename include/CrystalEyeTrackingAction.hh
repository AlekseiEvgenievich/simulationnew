// **************************************************************
//
//      ------ CrystalEyeTrackingAction  ------
//           by Libo Wu (23 Oct 2021)
//           Email: libo.wu@gssi.it
//
// **************************************************************

#ifndef CrystalEyeTrackingAction_h
#define CrystalEyeTrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "globals.hh"


class CrystalEyeEventAction;

class CrystalEyeTrackingAction : public G4UserTrackingAction{

		public:
				CrystalEyeTrackingAction(CrystalEyeEventAction* evt);
				~CrystalEyeTrackingAction();

				virtual void PreUserTrackingAction(const G4Track*);
				virtual void PostUserTrackingAction(const G4Track*);

		private:
				CrystalEyeEventAction* eventAct;

};

#endif
