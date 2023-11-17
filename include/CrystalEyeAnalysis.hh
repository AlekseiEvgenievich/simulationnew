// **************************************************************
//
//      ------ CrystalEyeAnalysis  ------
//           by Libo Wu (23 Sep 2021)
//           Email: libo.wu@gssi.it
//
// **************************************************************


#ifndef CrystalEyeAnalysis_h
#define CrystalEyeAnalysis_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <vector>

#include "g4root.hh"
//#include "g4xml.hh" 

//class CrystalEyeAnalysisMessenger;
class CrystalEyeDetectorConstruction;

class CrystalEyeAnalysis
{
	private:
		CrystalEyeAnalysis(G4String fname);
	public:
		virtual ~CrystalEyeAnalysis();

		void BeginOfRun();
		void EndOfRun();
		void BeginOfEvent();
		void EndOfEvent(G4int evtid);

		void SetEventData(
				G4int fgammaN,
			  G4int tag_firstHit,
			  G4int tag_interaction,
				G4int tag_evtInFov,
			  G4double fstructureE,	
				G4ThreeVector fpos
				);

		void SetPrimaryData(
				G4int  fpdg,
				G4double fenergy,
				G4ThreeVector fdirection,
				G4ThreeVector fposition
				);

		void SetDataCAL(
				G4double ftotedep,
				G4int    fnhits,
				G4double* fedep,
		    G4int*    forder		
				);

		void SetDataACD(
				G4double  ftotedep,
				G4int     fnhits,
				G4double* fedep,
		    G4int*    forder		
				);

		void SetDataCAL_Delay(
				G4double  ftotedep,
				G4double* fedep
				);

		void SetDataACD_Delay(
				G4double  ftotedep,
				G4double* fedep
				);

		static CrystalEyeAnalysis* getInstance(G4String fname);
		static CrystalEyeAnalysis* getInstance();

		void SetRootFilePath(G4String fname){
			rootFilePath = fname;
		}

	private:
		static CrystalEyeAnalysis* AnalysisInstance;
		const CrystalEyeDetectorConstruction* fDetector;

		G4String histFileName;
		G4String rootFilePath;
		
		G4int    fNofCells;
		G4int    fNofCellsACD;

		G4ThreeVector fPriPosition;
		G4ThreeVector fPriDirection;
		G4int    fPDG;
		G4float  fPrimaryE;
		G4int    fEvtID_ACD;
		G4int    fEvtID_CAL;

		G4float fTotEdep;
		G4float fTotEdepAcd;
		G4float fTotEdepAcd_Bottom;
		G4float fTotEdepCalUp;
		G4float fTotEdepCalDown;

		G4float fEdepStructure;
		G4int    fnHitsAcd;
		G4int    fnHitsCalUp;
		G4int    fnHitsCalDown;
		G4int    fGammaNum; // Amount of Gamma in TrackAction
		G4int    fTagFirstHit;
		G4int    fTagInteraction;
		G4int    fTagEvtInFOV;
		G4ThreeVector fFirstPos;

		G4int    fMaxIdAcd;
		G4float fMaxEdepAcd;
		G4int    fMaxIdUp;
		G4float fMaxEdepUp;
		G4int    fMaxIdDown;
		G4float fMaxEdepDown;

		G4float fTotEdepCal;
		G4int    fnHitsCal;
		G4float *fEdepCal;
		G4float *fEdepCalDelay;
		G4int    *fOrderCal;

		G4float *fEdepAcd;
		G4int    *fOrderAcd;
		G4float *fEdepAcdDelay;

		G4float fTotEdepDelay;
		G4float fTotEdepAcdDelay;
		G4float fTotEdepAcdDelay_Bottom;
		G4float fTotEdepCalDelay;
		G4float fTotEdepCalUpDelay;
		G4float fTotEdepCalDownDelay;



};

#endif

