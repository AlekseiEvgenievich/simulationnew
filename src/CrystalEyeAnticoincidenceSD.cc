// **************************************************************
//
//      ------------ CrystalEyeAnticoincidenceSD  ------
//           by Libo Wu (20 Sep 2021)
//           Email: libo.wu@gssi.it
//
// **************************************************************

#include "G4RunManager.hh"
#include "CrystalEyeAnticoincidenceHit.hh"
#include "CrystalEyeAnticoincidenceSD.hh"
#include "CrystalEyeDetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"



//---------------------------------------------
CrystalEyeAnticoincidenceSD::CrystalEyeAnticoincidenceSD(G4String name):G4VSensitiveDetector(name)
{
	G4RunManager* runManager = G4RunManager::GetRunManager();
	Detector = (CrystalEyeDetectorConstruction*)(runManager->GetUserDetectorConstruction());

	fNofCells  = Detector->GetNbOfPixel()+1;//UpACDs + DownACD
	fBkgTimeU  = Detector->GetBackgroundTimeU();//Second
	fBkgTimeD  = Detector->GetBackgroundTimeD();//Second
	fThrTime   = Detector->GetThresholdTime();//Microsecond

	fCHitID = new G4int[fNofCells];

	collectionName.insert("AnticoincidenceCollection");

	fhit_order = 0;
	tmp_pixelID= -1;
	flag_order = new G4bool[fNofCells];
}

//---------------------------------------------
CrystalEyeAnticoincidenceSD::~CrystalEyeAnticoincidenceSD()
{
	delete [] fCHitID;
	delete [] flag_order;
}

//---------------------------------------------
void CrystalEyeAnticoincidenceSD::Initialize(G4HCofThisEvent*)
{
	AnticoincidenceCollection = new CrystalEyeAnticoincidenceHitsCollection(
			SensitiveDetectorName, collectionName[0]);

	for(G4int i=0; i<fNofCells; i++){
		fCHitID[i] = -1;
		flag_order[i] = false;
	}

}

//---------------------------------------------
G4bool CrystalEyeAnticoincidenceSD::ProcessHits(G4Step *aStep, G4TouchableHistory* )
{
	G4double edep = 0.;
	edep  = aStep->GetTotalEnergyDeposit();

	G4TouchableHistory* theTouchable = 
		(G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());

	G4VPhysicalVolume* crystal_bar = theTouchable->GetVolume();
	G4String pixeName = crystal_bar->GetName();
	G4int pixelId = crystal_bar->GetCopyNo();

	G4StepPoint* fpreStepPoint = aStep->GetPreStepPoint();
	G4StepPoint* fpostStepPoint= aStep->GetPostStepPoint();
	G4double fedep    = aStep->GetTotalEnergyDeposit();
	G4double ftime    = fpreStepPoint->GetGlobalTime();
	G4int    ftrackid = aStep->GetTrack()->GetTrackID();
	G4int    fstepid  = aStep->GetTrack()->GetCurrentStepNumber();
	G4LogicalVolume* fpreVolum = fpreStepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
	G4String fvolumName = fpreVolum->GetName();
	G4String fparticleName = aStep->GetTrack()->GetDefinition()->GetParticleName();

	//G4cout
	//	<<"ftrackid: "<<ftrackid<<"  "
	//	<<"fstepid: "<<fstepid<<"  "
	//	<<"Volume Name: "<<fvolumName.c_str()<<"  "
	//	<<"Partcle Name: "<<fparticleName.c_str()<<"  "
	//	//<<"Parent ID: "<<fparentid<<"  "
	//	<<"Step Energy: "<<fedep*1000.<<"  "
	//	<<G4endl;

	if (edep == 0.) return false;


	if( fCHitID[pixelId]==-1 ){
		CrystalEyeAnticoincidenceHit* AnticoincidenceHit = new CrystalEyeAnticoincidenceHit;
		//AnticoincidenceHit->SetPos(aStep->GetPreStepPoint()->GetPosition());
		fCHitID[pixelId] = AnticoincidenceCollection->insert(AnticoincidenceHit)-1;

		//G4cout
		//	<<"------  ACDCollectionID = "
		//	<<fCHitID[pixelId]
		//	<<"  pixelId = "<<pixelId
		//	<<"  ------"
		//	<<G4endl;

	}
	if( fCHitID[pixelId]!=-1 ){
		if( ftime < fThrTime ){
			(*AnticoincidenceCollection)[fCHitID[pixelId]]->AddEnergy(edep);
			(*AnticoincidenceCollection)[fCHitID[pixelId]]->SetACDID(pixelId);
		}

		if( ftime>=fBkgTimeD && ftime<=fBkgTimeU ){
			(*AnticoincidenceCollection)[fCHitID[pixelId]]->AddEnergyDelay( edep );
		}

		if( tmp_pixelID!=pixelId && flag_order[pixelId]==false ){
			tmp_pixelID = pixelId;
			(*AnticoincidenceCollection)[fCHitID[pixelId]]->SetACDOrder(fhit_order);
			fhit_order++;
			flag_order[pixelId] = true;
		}

	}

	return true;
}

//---------------------------------------------
void CrystalEyeAnticoincidenceSD::EndOfEvent(G4HCofThisEvent* HCE){

	static G4int HCID = -1;
	if(HCID<0){
		HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
	}
	HCE -> AddHitsCollection(HCID, AnticoincidenceCollection);

	for(G4int i=0;i<fNofCells;i++){
		fCHitID[i]  = -1;
		flag_order[i] = false;
	}
	fhit_order  = 0;
	tmp_pixelID = -1;
}

//---------------------------------------------
void CrystalEyeAnticoincidenceSD::clear(){
}

//---------------------------------------------
void CrystalEyeAnticoincidenceSD::DrawAll(){
}

//---------------------------------------------
void CrystalEyeAnticoincidenceSD::PrintAll(){
}

//---------------------------------------------
