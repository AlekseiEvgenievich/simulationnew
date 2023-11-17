// **************************************************************
//
//      ------ CrystalEyeAnalysis.cc  ------
//           by Libo Wu (23 Sep 2021)
//           Email: libo.wu@gssi.it
//
// **************************************************************


#include <fstream>
#include <iomanip>

#include "G4RunManager.hh"

#include "CrystalEyeAnalysis.hh"
#include "CrystalEyeDetectorConstruction.hh"
//#include "CrystalEyeAnalysisMessenger.hh"

CrystalEyeAnalysis* CrystalEyeAnalysis::AnalysisInstance = 0;

//------------------------------------------------------------
CrystalEyeAnalysis::CrystalEyeAnalysis(G4String fname)
	: fDetector(0), rootFilePath("")
{
	fDetector = static_cast<const CrystalEyeDetectorConstruction*>
		(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

	fNofCells = fDetector->GetNbOfPixel();
	//fNofCellsACD = fNofCells + 1;//UpACDs+DownACD
	fNofCellsACD = fNofCells ;//UpACDs+DownACD


	fEdepCal       = new G4float[2*fNofCells];
	fEdepCalDelay  = new G4float[2*fNofCells];
	fOrderCal      = new G4int[2*fNofCells];
	  
	fEdepAcd       = new G4float[fNofCellsACD];
	fOrderAcd      = new G4int[fNofCellsACD]; 
	fEdepAcdDelay  = new G4float[fNofCellsACD];
	
	histFileName = "CrystalEyeMC_" + fname + ".root";
}

//------------------------------------------------------------
CrystalEyeAnalysis::~CrystalEyeAnalysis()
{
	delete [] fEdepCal;  fEdepCal  = NULL;
	delete [] fOrderCal; fOrderCal = NULL;
	delete [] fEdepAcd;  fEdepAcd  = NULL;
	delete [] fOrderAcd; fOrderAcd = NULL;

	delete [] fEdepCalDelay; fEdepCalDelay  = NULL;
	delete [] fEdepAcdDelay; fEdepAcdDelay = NULL;

	delete G4AnalysisManager::Instance();
}

//------------------------------------------------------------
CrystalEyeAnalysis* CrystalEyeAnalysis::getInstance(G4String fname)
{
	if(AnalysisInstance == 0) 
		AnalysisInstance = new CrystalEyeAnalysis(fname);
	return AnalysisInstance;
}

//------------------------------------------------------------
CrystalEyeAnalysis* CrystalEyeAnalysis::getInstance()
{
	if(AnalysisInstance == 0){
		G4cout << "CrystalEyeAnalysis::getInstance(): Please Pass The ROOT File Name!!!" << G4endl;
		exit(0);
	}
	else
		return AnalysisInstance;
}

//------------------------------------------------------------
void CrystalEyeAnalysis::BeginOfRun()
{
	G4AnalysisManager* man = G4AnalysisManager::Instance();
	man->OpenFile( rootFilePath + histFileName );

	man->SetFirstNtupleId(1);

	// ---  1  ---
	man->CreateNtuple("Total","Total");                   // 
	man->CreateNtupleFColumn("TotalEdep");			          // 0
	man->CreateNtupleFColumn("TotalEdep_ACD");            // 1
//	man->CreateNtupleFColumn("TotalEdep_BottomACD");      // 2
	man->CreateNtupleFColumn("TotalEdep_UP");	            // 3
	man->CreateNtupleFColumn("TotalEdep_DOWN");           // 4

	man->CreateNtupleFColumn("StructureEdep");            // 5 
	man->CreateNtupleIColumn("Interaction");              // 6  
	man->CreateNtupleIColumn("EventInFOV");               // 7  
	man->CreateNtupleFColumn("TotalEdepDelay");			      // 8  
	man->CreateNtupleFColumn("TotalEdepDelay_ACD");       // 9 
	man->CreateNtupleFColumn("TotalEdepDelay_CAL");       // 10
//	man->CreateNtupleFColumn("TotalEdepDelay_BottomACD"); // 11
	man->FinishNtuple();                                  
                                                        
	// ---  2  ---                                        
	man->CreateNtuple("Anticoincidence","Anticoincidence");
	man->CreateNtupleIColumn("PixelID_ACD");              
	man->CreateNtupleFColumn("Edep_ACD");
	man->CreateNtupleIColumn("Order_ACD");
	man->CreateNtupleFColumn("EdepDelay_ACD");
	man->CreateNtupleIColumn("EventID");
	man->FinishNtuple();

	// ---  3  ---
	man->CreateNtuple("Calorimeter","Calorimeter");
	man->CreateNtupleIColumn("PixelID_CAL");
	man->CreateNtupleFColumn("Edep_CAL");
	man->CreateNtupleIColumn("Order_CAL");
	man->CreateNtupleFColumn("EdepDelay_CAL");
	man->CreateNtupleIColumn("EventID");
	man->FinishNtuple();

	// ---  4  ---
	man->CreateNtuple("PrimaryParameter","PrimaryParameter");
	man->CreateNtupleIColumn("PDG");    //0
	man->CreateNtupleFColumn("Kinetic");//1
	man->CreateNtupleFColumn("Dir_x");	//2
	man->CreateNtupleFColumn("Dir_y");	//3
	man->CreateNtupleFColumn("Dir_z");	//4
	man->CreateNtupleFColumn("Pos_x");	//5
	man->CreateNtupleFColumn("Pos_y");	//6
	man->CreateNtupleFColumn("Pos_z");  //7	
	man->CreateNtupleIColumn("EventID_ACD"); //8
	man->CreateNtupleIColumn("Hits_ACD");//9
	man->CreateNtupleIColumn("EventID_CAL"); //10
	man->CreateNtupleIColumn("Hits_CAL");//11
	man->CreateNtupleIColumn("Int_x");//12
	man->CreateNtupleIColumn("Int_y");//13
	man->CreateNtupleIColumn("Int_z");//14
	man->FinishNtuple();

	fEvtID_ACD=0;
	fEvtID_CAL=0;
}

//------------------------------------------------------------
void CrystalEyeAnalysis::EndOfRun()
{
	G4AnalysisManager* man = G4AnalysisManager::Instance();
	man->Write();
	man->CloseFile();

}

//------------------------------------------------------------
void CrystalEyeAnalysis::BeginOfEvent(){

	fTotEdep           = 0;
	fTotEdepAcd        = 0;
//	fTotEdepAcd_Bottom = 0;
	fTotEdepCalUp      = 0;
	fTotEdepCalDown    = 0;

	fTotEdepDelay           = 0;
	fTotEdepAcdDelay        = 0;
//	fTotEdepAcdDelay_Bottom = 0;
	fTotEdepCalUpDelay      = 0;
	fTotEdepCalDownDelay    = 0;

	fEdepStructure  = 0;
	fnHitsAcd       = 0;
	fnHitsCalUp     = 0;
	fnHitsCalDown   = 0;
	fGammaNum       = 0;
	fTagFirstHit    = 0;
	fTagInteraction = 0;
	fTagEvtInFOV    = 0;
	fFirstPos       = G4ThreeVector(-999,-999,-999);

	fMaxIdAcd   = -1;
	fMaxEdepAcd =  0;
  fMaxIdUp    = -1;
  fMaxEdepUp  =  0;
  fMaxIdDown  = -1;	
	fMaxEdepDown=  0;

	fTotEdepCal =  0;
	fnHitsCal   =  0;
	for(int i=0; i<2*fNofCells; i++){
		fEdepCal[i]  = 0;
		fOrderCal[i] = -1;

		fEdepCalDelay[i] = 0;
	}

	fTotEdepAcd        = 0;
	fTotEdepAcdDelay   = 0;
//	fTotEdepAcd_Bottom = 0;
	fnHitsAcd   = 0;
	for(int i=0; i<fNofCellsACD; i++){
		fEdepAcd[i]  = 0;
		fOrderAcd[i] = -1;

		fEdepAcdDelay[i] = 0;
	}


}

//------------------------------------------------------------
void CrystalEyeAnalysis::EndOfEvent(G4int evtid)
{
//	fTotEdep = fTotEdepAcd + fTotEdepCal + fTotEdepAcd_Bottom;
	fTotEdep = fTotEdepAcd + fTotEdepCal;
//	fTotEdepDelay = fTotEdepAcdDelay + fTotEdepCalDelay + fTotEdepAcdDelay_Bottom;
	fTotEdepDelay = fTotEdepAcdDelay + fTotEdepCalDelay ;

	G4AnalysisManager* man = G4AnalysisManager::Instance();

	// ---  1: TotalInfo  ---
	man->FillNtupleFColumn(1, 0, fTotEdep);
	man->FillNtupleFColumn(1, 1, fTotEdepAcd);
//	man->FillNtupleFColumn(1, 2, fTotEdepAcd_Bottom);
	man->FillNtupleFColumn(1, 2, fTotEdepCalUp);
	man->FillNtupleFColumn(1, 3, fTotEdepCalDown);
	man->FillNtupleFColumn(1, 4, fEdepStructure);
	man->FillNtupleIColumn(1, 5, fTagInteraction);
	man->FillNtupleIColumn(1, 6, fTagEvtInFOV);
	man->FillNtupleFColumn(1, 7, fTotEdepDelay);
	man->FillNtupleFColumn(1, 8, fTotEdepAcdDelay);
	man->FillNtupleFColumn(1, 9, fTotEdepCalDelay);
//	man->FillNtupleFColumn(1, 11,fTotEdepAcdDelay_Bottom);
	man->AddNtupleRow(1);

	// ---  2: AnticoincidenceInfo  ---
	int FirstID_ACD = -1;
	int fHits_ACD    = -1;
	for(int i=0; i<fNofCellsACD; i++){
		if(fEdepAcd[i]>0.0000001){
			man->FillNtupleIColumn(2, 0, i);
			man->FillNtupleFColumn(2, 1, fEdepAcd[i] );
			man->FillNtupleIColumn(2, 2, fOrderAcd[i]);
			man->FillNtupleFColumn(2, 3, fEdepAcdDelay[i] );
			man->FillNtupleIColumn(2, 4, evtid);
			man->AddNtupleRow(2);

			if(FirstID_ACD==-1){
				FirstID_ACD = fEvtID_ACD;
				fHits_ACD = 0;
			}
			fEvtID_ACD++;
			fHits_ACD++;
		}
	}

	// ---  3: CalorimeterInfo  ---
	int FirstID_CAL = -1;
	int fHits_CAL    = -1;
	for(int i=0; i<2*fNofCells; i++){
		if(fEdepCal[i]>0.000001){
			man->FillNtupleIColumn(3, 0, i);
			man->FillNtupleFColumn(3, 1, fEdepCal[i] );
			man->FillNtupleIColumn(3, 2, fOrderCal[i]);
			man->FillNtupleFColumn(3, 3, fEdepCalDelay[i] );
			man->FillNtupleIColumn(3, 4, evtid);
			man->AddNtupleRow(3);

			if(FirstID_CAL==-1){
				FirstID_CAL = fEvtID_CAL;
				fHits_CAL    = 0;
			}
			fEvtID_CAL++;
			fHits_CAL++;
		}
	}

	// ---  4: PrimaryInfo  ---
	man->FillNtupleIColumn(4, 0,  fPDG);
	man->FillNtupleFColumn(4, 1,  fPrimaryE );
	man->FillNtupleFColumn(4, 2,  fPriDirection.x() );
	man->FillNtupleFColumn(4, 3,  fPriDirection.y() );
	man->FillNtupleFColumn(4, 4,  fPriDirection.z() );
	man->FillNtupleFColumn(4, 5,  fPriPosition.x() );
	man->FillNtupleFColumn(4, 6,  fPriPosition.y() );
	man->FillNtupleFColumn(4, 7,  fPriPosition.z() );
	man->FillNtupleIColumn(4, 8,  FirstID_ACD);
	man->FillNtupleIColumn(4, 9,  fHits_ACD);
	man->FillNtupleIColumn(4, 10, FirstID_CAL);
	man->FillNtupleIColumn(4, 11, fHits_CAL);
	man->FillNtupleIColumn(4, 12, fFirstPos.x());
	man->FillNtupleIColumn(4, 13, fFirstPos.y());
	man->FillNtupleIColumn(4, 14, fFirstPos.z()); 
	man->AddNtupleRow(4);

}

//------------------------------------------------------------
void CrystalEyeAnalysis::SetEventData(
		G4int fgammaN, G4int tag_firstHit, G4int tag_interaction, G4int tag_evtInFov, G4double fstructureE, G4ThreeVector fpos ){
	fGammaNum       = fgammaN;
	fTagFirstHit    = tag_firstHit;
	fTagInteraction = tag_interaction;
	fTagEvtInFOV    = tag_evtInFov;

	fEdepStructure  = fstructureE;
	fFirstPos       = fpos;
}

//------------------------------------------------------------
void CrystalEyeAnalysis::SetPrimaryData(
		G4int fpdg, G4double fenergy, G4ThreeVector fdirection, G4ThreeVector fposition){
	fPDG      = 0;
	fPrimaryE = 0;
	fPriDirection = G4ThreeVector(0,0,0);
	fPriPosition  = G4ThreeVector(0,0,0);

	fPDG      = fpdg;
	fPrimaryE = fenergy;
	fPriDirection = fdirection;
	fPriPosition  = fposition;

}

//------------------------------------------------------------
void CrystalEyeAnalysis::SetDataCAL(
		G4double ftotedep, G4int fnhits, G4double* fedep, G4int* forder){

	fTotEdepCal = ftotedep;
	fnHitsCal   = fnhits;
	for(int i=0; i<2*fNofCells; i++){
		fEdepCal[i]  = fedep[i];
		fOrderCal[i] = forder[i];

		if( i<fNofCells ){
			if( fEdepCal[i]>fMaxEdepUp ){
				fMaxEdepUp = fEdepCal[i];
				fMaxIdUp   = i;
			}
			if( fEdepCal[i]>0 ){
				fTotEdepCalUp += fEdepCal[i];
				fnHitsCalUp++;
			}
		}
		else{
			if( fEdepCal[i]>fMaxEdepDown ){
				fMaxEdepDown = fEdepCal[i];
				fMaxIdDown   = i;
			}
			if( fEdepCal[i]>0 ){
				fTotEdepCalDown += fEdepCal[i];
				fnHitsCalDown++;
			}
		}

	}

}

//------------------------------------------------------------
void CrystalEyeAnalysis::SetDataCAL_Delay(G4double ftotedep, G4double* fedep){

	fTotEdepCalDelay = ftotedep;
	for(int i=0; i<2*fNofCells; i++){
		fEdepCal[i]  = fedep[i];
	}

}

//------------------------------------------------------------
void CrystalEyeAnalysis::SetDataACD(
		G4double ftotedep, G4int fnhits, G4double* fedep, G4int* forder){

	fTotEdepAcd        = ftotedep - fedep[fNofCellsACD-1];
//	fTotEdepAcd_Bottom = fedep[fNofCellsACD-1]; 
	fnHitsAcd          = fnhits;
	for(int i=0; i<fNofCellsACD; i++){
		fEdepAcd[i]  = fedep[i];
		fOrderAcd[i] = forder[i];

		if(fEdepAcd[i]>fMaxEdepAcd){
			fMaxEdepAcd = fEdepAcd[i];
			fMaxIdAcd   = i;
		}
	}

}

//------------------------------------------------------------
void CrystalEyeAnalysis::SetDataACD_Delay( G4double ftotedep, G4double* fedep ){

	fTotEdepAcdDelay        = ftotedep - fedep[fNofCellsACD-1];
//	fTotEdepAcdDelay_Bottom = fedep[fNofCellsACD-1];

	for(int i=0; i<fNofCellsACD; i++){
		fEdepAcdDelay[i] = fedep[i];
	}

}

//------------------------------------------------------------
