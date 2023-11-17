// **************************************************************
//
//      ------------ CrystalEyeEventAction  -----------
//           by Libo Wu (23 Sep 2021)
//           Email: libo.wu@gssi.it
//
// **************************************************************

#ifndef CrystalEyeRunAction_h
#define CrystalEyeRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;
class CrystalEyeAnticoincidenceDataManager;
class CrystalEyeCalorimeterDataManager;

class CrystalEyeRunAction : public G4UserRunAction
{
		public:
				CrystalEyeRunAction(
						CrystalEyeAnticoincidenceDataManager* aDataManagerAcd,
						CrystalEyeCalorimeterDataManager* aDataManagerCal);
				~CrystalEyeRunAction();
				
				void	BeginOfRunAction(const G4Run*);
				void  EndOfRunAction(const G4Run*);

		private:
				G4int fRunID; 

				CrystalEyeAnticoincidenceDataManager* fDataManagerAcd;
				CrystalEyeCalorimeterDataManager* fDataManagerCal;

};

#endif

