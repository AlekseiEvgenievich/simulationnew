// **************************************************************
//
//      ------ CrystalEyePhysicsList.cc  ------
//           Modified by Libo Wu (23 Sep 2021)
//           Email: libo.wu@gssi.it
//
// **************************************************************

#include "CrystalEyePhysicsList.hh"
#include "CrystalEyePhysicsListMessenger.hh"
#include "G4PhysListFactory.hh"
#include "G4VPhysicsConstructor.hh"

#include "G4Decay.hh"
#include "G4StepLimiter.hh"
#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ProcessManager.hh"
#include "G4IonFluctuations.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4EmProcessOptions.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmLivermorePolarizedPhysics.hh" // TBC

#include "G4HadronElasticPhysics.hh"
#include "G4HadronDElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

// ------------------------------------------------------------
CrystalEyePhysicsList::CrystalEyePhysicsList() : G4VModularPhysicsList()
{
	G4LossTableManager::Instance();

	verboseLevel    = 1;
	defaultCutValue = 100*micrometer;
	SetVerboseLevel( verboseLevel );

	cutForGamma     = defaultCutValue;
	cutForElectron  = defaultCutValue;
	cutForPositron  = defaultCutValue;

	pMessenger = new CrystalEyePhysicsListMessenger(this);

	//---  Decay physics  ---
	fDecaPhys = new G4DecayPhysics("decays");
	G4cout << ">>> PhysicsList:: G4DecayPhysics is activated! <<<" << G4endl;
}

// ------------------------------------------------------------
CrystalEyePhysicsList::~CrystalEyePhysicsList()
{
	delete pMessenger;

	delete fDecaPhys;
	for(size_t i=0; i<fElecPhys.size(); i++){
		delete fElecPhys[i];
	}
	for(size_t i=0; i<fHadrPhys.size(); i++){
		delete fHadrPhys[i];
	}

}

// ------------------------------------------------------------
void CrystalEyePhysicsList::AddPackage(const G4String& name)
{
	//---   For PhysicsListMessenger  ----
	G4PhysListFactory factory;
	G4VModularPhysicsList* phys =factory.GetReferencePhysList(name);
	G4int i=0;
	const G4VPhysicsConstructor* elem= phys->GetPhysics(i);
	G4VPhysicsConstructor* tmp = const_cast<G4VPhysicsConstructor*> (elem);
	while (elem !=0)
	{
		RegisterPhysics(tmp);
		elem= phys->GetPhysics(++i) ;
		tmp = const_cast<G4VPhysicsConstructor*> (elem);
	}
}

// ------------------------------------------------------------
void CrystalEyePhysicsList::ConstructParticle()
{
	G4BosonConstructor  pBosonConstructor;
	pBosonConstructor.ConstructParticle();

	G4LeptonConstructor pLeptonConstructor;
	pLeptonConstructor.ConstructParticle();

	G4MesonConstructor pMesonConstructor;
	pMesonConstructor.ConstructParticle();

	G4BaryonConstructor pBaryonConstructor;
	pBaryonConstructor.ConstructParticle();

	G4IonConstructor pIonConstructor;
	pIonConstructor.ConstructParticle();

	G4ShortLivedConstructor pShortLivedConstructor;
	pShortLivedConstructor.ConstructParticle();
}

// ------------------------------------------------------------
void CrystalEyePhysicsList::ConstructProcess()
{
	// transportation
	AddTransportation();

	// decay physics list
	fDecaPhys->ConstructProcess();

	// electromagnetic physics list
	for(size_t i=0; i<fElecPhys.size(); i++){
		fElecPhys[i]->ConstructProcess();
		em_config.AddModels();
	}

	// hadronic physics lists
	for(size_t i=0; i<fHadrPhys.size(); i++) {
		fHadrPhys[i]->ConstructProcess();
	}

}

// ------------------------------------------------------------
void CrystalEyePhysicsList::AddPhysicsList(const G4String& name)
{

	if (verboseLevel>1) {
		G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
	}

	//---  ELECTROMAGNETIC MODELS  ---
	if (name == "Em_Standard_opt3") {
		fElecPhys.push_back( new G4EmStandardPhysics_option3( verboseLevel ) );
		G4cout << ">>> PhysicsList:: G4EmStandardPhysics_option3 is activated! <<<" << G4endl;
	} 
	else if (name == "Em_Livermore") {
		fElecPhys.push_back( new G4EmLivermorePhysics( verboseLevel ) );
		G4cout << ">>> PhysicsList:: G4EmLivermorePhysics is activated! <<<" << G4endl;
	} 
	else if (name == "Em_Penelope") {
		fElecPhys.push_back( new G4EmPenelopePhysics( verboseLevel ) );
		G4cout << ">>> PhysicsList:: G4EmPenelopePhysics is activated! <<<" << G4endl;
	} 
	else if (name == "Em_Polarized") {
		fElecPhys.push_back( new G4EmLivermorePolarizedPhysics( verboseLevel ) );
		G4cout << ">>> PhysicsList:: G4EmPenelopePhysics is activated! <<<" << G4endl;
	} 

	//---  HADRONIC MODELS  ---
	//Hadron Elastic Physics List
	else if (name == "Hd_Elastic") {
		fHadrPhys.push_back( new G4HadronElasticPhysics( verboseLevel ));
		G4cout << ">>> PhysicsList:: G4HadronElasticPhysics is activated! <<<" << G4endl;
	} 
	else if (name == "Hd_DElastic" ) {
		fHadrPhys.push_back( new G4HadronDElasticPhysics( verboseLevel ));
		G4cout << ">>> PhysicsList:: G4HadronDElasticPhysics is activated! <<<" << G4endl;
	} 
	else if (name == "Hd_ElasticHP" ) {
		fHadrPhys.push_back( new G4HadronElasticPhysicsHP( verboseLevel ));
		G4cout << ">>> PhysicsList:: G4HadronElasticPhysicsHP is activated! <<<" << G4endl;
	} 
	else if (name == "QGSP_BIC_HP" ) {
		fHadrPhys.push_back(new G4HadronPhysicsQGSP_BIC_HP( verboseLevel ));
		G4cout << ">>> PhysicsList:: G4HadronPhysicsQGSP_BIC_HP is activated! <<<" << G4endl;
	} 
	else if (name == "Hd_IonBinary" ) {
		fHadrPhys.push_back(new G4IonBinaryCascadePhysics( verboseLevel ));
		G4cout << ">>> PhysicsList:: G4IonBinaryCascadePhysics is activated! <<<" << G4endl;
	} 
	else if (name == "Hd_RadioactiveDecay" ) {
		fHadrPhys.push_back(new G4RadioactiveDecayPhysics( verboseLevel ));
		G4cout << ">>> PhysicsList:: G4RadioactiveDecayPhysics is activated! <<<" << G4endl;
	} 
	else {
		G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
			<< " is not defined!  The choice list:"
			<< " Em_Standard_opt3, Em_Livermore, Em_Penelope, Em_Polarized,"
			<< " Hd_Elastic, Hd_DElastic, Hd_ElasticHP, Hd_IonBinary, Hd_RadioactiveDecay, QGSP_BIC_HP"
			<<"  \n"
			<< G4endl;

		exit(0);
	}
}

// ------------------------------------------------------------
void CrystalEyePhysicsList::SetCuts()
{

	if (verboseLevel >0){
		G4cout << "PhysicsList::SetCuts:";
		G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
	}

	G4double lowLimit  = 250. * eV;
	G4double highLimit = 100. * GeV;
	G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowLimit, highLimit);

	// set cut values for gamma at first and for e- second and next for e+,
	// because some processes for e+/e- need cut values for gamma
	SetCutValue(cutForGamma, "gamma");
	SetCutValue(cutForElectron, "e-");
	SetCutValue(cutForPositron, "e+");

	if (verboseLevel>0) DumpCutValuesTable();
}

void CrystalEyePhysicsList::SetCutForGamma(G4double cut)
{
	cutForGamma = cut;
	SetParticleCuts(cutForGamma, G4Gamma::Gamma());
}

void CrystalEyePhysicsList::SetCutForElectron(G4double cut)
{
	cutForElectron = cut;
	SetParticleCuts(cutForElectron, G4Electron::Electron());
}

void CrystalEyePhysicsList::SetCutForPositron(G4double cut)
{
	cutForPositron = cut;
	SetParticleCuts(cutForPositron, G4Positron::Positron());
}

