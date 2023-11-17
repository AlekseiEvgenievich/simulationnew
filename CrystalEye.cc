// **************************************************************
//
//      ------ CrystalEye.cc  ------
//           by Libo Wu (23 Sep 2021)
//           Email: libo.wu@gssi.it
//
// **************************************************************


//#include "G4Types.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif
  
#include "G4UImanager.hh"
#include "Randomize.hh"

#include "CrystalEyeDetectorConstruction.hh"
#include "CrystalEyePhysicsList.hh"
#include "CrystalEyeActionInitialization.hh"
#include "CrystalEyeAnalysis.hh"

#include "QGSP_BIC.hh"
#include "FTFP_BERT.hh"

using namespace std;

// ------------------------------------------------------------
int main(int argc, char** argv)
{
	G4cout << "=====  argc = " << argc << "  =====" << G4endl;
	for(int i=0; i<argc; i++){
		G4cout << " argc["<<i<<"] = "<<argv[i]<<G4endl;
	}
	G4cout<<"==========================\n"<<G4endl;

	//---  Choose the random engine  ----
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	G4long fRandomSeed = time(0)+getpid();
	if(argc==4){
		G4long fseed = atoi( argv[3] );
		fRandomSeed = fseed;
	}
	G4cout<<"===>  Seed = "<<fRandomSeed<<"  <==="<<G4endl;
	CLHEP::HepRandom::setTheSeed( fRandomSeed );
	CLHEP::HepRandom::showEngineStatus();

#ifdef G4MULTITHREADED
	G4MTRunManager* runManager = new G4MTRunManager;
	runManager->SetNumberOfThreads(G4Threading::G4GetNumberOfCores());
#else
	G4RunManager* runManager = new G4RunManager;
#endif  

#ifdef G4VIS_USE
	G4VisManager* visManager = new G4VisExecutive;
	visManager->Initialize();
#endif

	runManager->SetUserInitialization(new CrystalEyeDetectorConstruction);
	runManager->SetUserInitialization(new CrystalEyePhysicsList);
	//runManager->SetUserInitialization(new FTFP_BERT);//QGSP_BIC
	runManager->SetUserInitialization(new CrystalEyeActionInitialization());

	CrystalEyeAnalysis* analysis;
	if(argc!=1)	analysis = CrystalEyeAnalysis::getInstance(argv[2]);
	else        analysis = CrystalEyeAnalysis::getInstance("run");
	//analysis -> SetRootFilePath("../result/");
	analysis -> SetRootFilePath("../result/expo/");

	//---  Get the pointer to the UI manager  ---
	G4UImanager* UImanager = G4UImanager::GetUIpointer();
	if (argc!=1)  
	{
		G4String command = "/control/execute ";
		G4String fileName = argv[1];
		UImanager->ApplyCommand("/control/execute ../macros/physics.mac");
		UImanager->ApplyCommand(command+fileName);

	}
	else
	{
		G4UIExecutive* ui = new G4UIExecutive(argc, argv);
		if (ui->IsGUI())
		{
			UImanager->ApplyCommand("/control/execute ../macros/physics.mac");
			UImanager->ApplyCommand("/control/execute ../macros/vis.mac");
//			UImanager->ApplyCommand("/control/execute ../macros/onepixel/10MeV.mac");
//			UImanager->ApplyCommand("/control/execute ../macros/circle_theta_phi_pow/gamma_circle_Angle0_0_E10_50000.mac");
			//UImanager->ApplyCommand("/control/execute ../macros/run.mac ");
//			UImanager->ApplyCommand("/control/execute ../macros/test.mac ");

			ui->SessionStart();
		}
		delete ui;
	}

#ifdef G4VIS_USE
	delete visManager;
#endif
	delete analysis;
	delete runManager;

	return 0;
}
