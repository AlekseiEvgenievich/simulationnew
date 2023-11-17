// **************************************************************
//
//      ------ CrystalEyeAnalysis  ------
//           Modified by Libo Wu (06 May 2022)
//           Email: libo.wu@gssi.it
//
// **************************************************************

#include "CrystalEyeDetectorConstruction.hh"
#include "CrystalEyeDetectorMessenger.hh"
#include "CrystalEyeAnticoincidenceSD.hh"
#include "CrystalEyeCalorimeterSD.hh"
#include "CADMesh.hh"
 
#include "G4NistManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UImanager.hh"

#include "G4Transform3D.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"

#include <string>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <fstream> //header necessario per ofstream!
#include <cstring>


#include "G4RegionStore.hh"

G4ThreadLocal G4GlobalMagFieldMessenger* 
CrystalEyeDetectorConstruction::fMagFieldMessenger = 0;

//------------------------------------------------------------

CrystalEyeDetectorConstruction::CrystalEyeDetectorConstruction()
	:solidWorld(0),logicWorld(0),physiWorld(0),
	solidUpStructure(0),logicUpStructure(0),physiUpStructure(0),
	solidUpStructure2(0),
	uppix(0),logicUpPixel(0),physiUpPixel(0),
	solidDownStructure(0),logicDownStructure(0),physiDownStructure(0),
	downpix(0),logicDownPixel(0),physiDownPixel(0),
	solidACD(0),logicACD(0),physiACD(0)
{

	//StructurePath = "/Users/lb/work/CrystalEye/CrystalEyeCode/"
	StructurePath = "../../";
	checkOverlaps = false;
	fCadUnit = 0.1*mm;

	// default parameter values of the payload
	NbOfPixel = 3;
	NbOfLayers = 2;
	ACDThickness = 0.5*cm;

	fBackgroundTimeD = 30*CLHEP::s;
	fBackgroundTimeU = 100*24*60*60*CLHEP::s;
	fThresholdTime   = 50*CLHEP::us;

	G4cout
		<<"  \n@@@@@"
		<<"  ACDThickness = "<<ACDThickness<<"mm"
		<<"  fBkgTimeUp = "<<fBackgroundTimeU<<"ns"
		<<"  fBkgTimeDown = "<<fBackgroundTimeD<<"ns"
		<<"  fThresholdTime = "<<fThresholdTime<<"ns"
		<<"  @@@@@@@\n"<<G4endl;

	calorimeterSD.Put(0);
	anticoincidenceSD.Put(0);

	detectorMessenger = new CrystalEyeDetectorMessenger(this);
}

//------------------------------------------------------------

CrystalEyeDetectorConstruction::~CrystalEyeDetectorConstruction()
{ delete detectorMessenger;}

//------------------------------------------------------------

G4VPhysicalVolume* CrystalEyeDetectorConstruction::Construct()
{
	DefineMaterials();
	return ConstructPayload();
}

//------------------------------------------------------------

void CrystalEyeDetectorConstruction::DefineMaterials()
{ 
	G4NistManager* man = G4NistManager::Instance();

	G4String name, symbol;    
	G4double a, z, density;            

	G4int ncomponents;
	G4int natoms;
	G4double fractionmass;
	G4double temperature, pressure;

	//
	// define Elements
	//

	/*
		 a= 55.845*g/mole;
		 G4Element* Iron = new G4Element(name="Iron", symbol="Fe", z=26.,a);

		 a = 126.904*g/mole;
		 G4Element* I  = new G4Element(name="Iodine"  ,symbol="I" , z= 53., a);

		 a = 132.905*g/mole;
		 G4Element* Cs  = new G4Element(name="Cesium"  ,symbol="Cs" , z= 55., a);

		 a = 207.19*g/mole;
		 G4Element* Lead = new G4Element(name="Lead", symbol="Pb", z=82., a);
	 */

	//a = 174.97*g/mole;
	//G4Element* Lu = new G4Element(name="Luthetium", symbol="Lu", z=71. , a);
	//a = 88.9059*g/mole;
	//G4Element* Y = new G4Element(name="Yttrium", symbol="Y", z=39. , a);
	//a = 28.09*g/mole;
	//G4Element* Silicon = new G4Element(name="Silicon", symbol="Si", z=14., a);
	//a = 15.99*g/mole;
	//G4Element* O  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);
	//a = 14.006*g/mole;
	//G4Element* N  = new G4Element(name="Nitrogen"  ,symbol="N" , z= 7., a);
	//a = 1.01*g/mole;
	//G4Element* H  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);
	//a = 12.01*g/mole;
	//G4Element* C  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);
	//a = 26.98*g/mole;
	//G4Element* Alumin = new G4Element(name="Aluminum"  ,symbol="A G4Element* Alumin = man->FindOrBuildElement("Al"); G4Element* Alumin = man->FindOrBuildElement("Al"); G4Element* Alumin = man->FindOrBuildElement("Al"); G4Element* Alumin = man->FindOrBuildElement("Al");l" , z= 13., a);

	G4Element* Lu = man->FindOrBuildElement("Lu");
	G4Element* Y  = man->FindOrBuildElement("Y");
	G4Element* Si = man->FindOrBuildElement("Si");
	G4Element* O  = man->FindOrBuildElement("O");
	G4Element* N  = man->FindOrBuildElement("N");
	G4Element* H  = man->FindOrBuildElement("H");
	G4Element* C  = man->FindOrBuildElement("C");
	G4Element* Alumin = man->FindOrBuildElement("Al");

	//
	// define simple materials
	//

	G4Material* LSO = new G4Material("LSO",7.4*g/cm3,3);
	LSO->AddElement(Lu,2);
	LSO->AddElement(Si,1);
	LSO->AddElement(O, 5);
	G4Material* YSO = new G4Material("YSO",4.5*g/cm3,3);
	YSO->AddElement(Y, 2);
	YSO->AddElement(Si,1);
	YSO->AddElement(O, 5);
	G4Material* lyso = new G4Material("LYSO",7.1*g/cm3,2);
	lyso->AddMaterial(LSO,90*perCent);
	lyso->AddMaterial(YSO,10*perCent);

	//density = 7.3*g/cm3;
	//G4Material* lyso = new G4Material(name="Lu2SiO5", density, ncomponents=4);
	//lyso->AddElement(Lu, natoms=1.8);//????
	//lyso->AddElement(Y, fractionmass=0.2);
	//lyso->AddElement(Silicon, natoms=1.);
	//lyso->AddElement(O , natoms=5.);

	//
	// define a material from elements.   case 1: chemical molecule
	//

	density = 1.032*g/cm3;
	G4Material* Sci = new G4Material(name="Scintillator", density, ncomponents=2);
	Sci->AddElement(C, natoms=9);
	Sci->AddElement(H, natoms=10);

	//
	// define a material from elements.   case 2: mixture by fractional mass
	//

	density = 1.290*mg/cm3;
	G4Material* Air = new G4Material(name="Air"  , density, ncomponents=2);
	Air->AddElement(N, fractionmass=0.7);
	Air->AddElement(O, fractionmass=0.3);

	density = 2.700*g/cm3;
	G4Material* Al = new G4Material(name="Aluminum", density, ncomponents=1);
	Al->AddElement(Alumin, fractionmass=1.);

	//
	// examples of vacuum
	//

	density     = universe_mean_density;    //from PhysicalConstants.h
	pressure    = 3.e-18*pascal;
	temperature = 2.73*kelvin;
	G4Material* vacuum = new G4Material(name="Galactic", z=1., a=1.01*g/mole, density,kStateGas,temperature,pressure);

	density     = 1.e-5*g/cm3;
	pressure    = 2.e-2*bar;
	temperature = STP_Temperature;         //from PhysicalConstants.h
	G4Material* beam = new G4Material(name="Beam", density, ncomponents=1,
			kStateGas,temperature,pressure);
	beam->AddMaterial(Air, fractionmass=1.);
	//G4cout << *(G4Material::GetMaterialTable()) << G4endl;

	//default materials of the payload

	defaultMaterial  = vacuum;
	ACDMaterial = Sci;
	PixelMaterial = lyso;
	structureMaterial = Al;

}

//------------------------------------------------------------

G4VPhysicalVolume* CrystalEyeDetectorConstruction::ConstructPayload()
{
	// Some visualization styles
	// Old Version
	//G4Colour(0.3,0.8,0.1), G4Colour(0.2,0.3,0.8), G4Colour(0.8,0.2,0.3)

	G4VisAttributes* VisAtt_ACD= new G4VisAttributes();
	VisAtt_ACD->SetVisibility(true);
	VisAtt_ACD->SetForceSolid(TRUE);
	VisAtt_ACD->SetColour(G4Colour::Cyan());

	G4VisAttributes* VisAtt_BottomACD= new G4VisAttributes();
	VisAtt_BottomACD->SetVisibility(true);
	VisAtt_BottomACD->SetForceSolid(TRUE);
	VisAtt_BottomACD->SetColour(G4Colour::Green());

	G4VisAttributes* VisAtt_Up= new G4VisAttributes();
	VisAtt_Up->SetVisibility(true);
	VisAtt_Up->SetForceSolid(TRUE);
	VisAtt_Up->SetColour(G4Colour::White());

	G4VisAttributes* VisAtt_Down= new G4VisAttributes();
	VisAtt_Down->SetVisibility(true);
	VisAtt_Down->SetForceSolid(TRUE);
	VisAtt_Down->SetColour(G4Colour::Yellow());

	G4VisAttributes* VisAtt_Structure= new G4VisAttributes();
	VisAtt_Structure->SetVisibility(true);
	VisAtt_Structure->SetForceWireframe(TRUE);
	VisAtt_Structure->SetColour(G4Colour::Red());

	// complete the Payload parameters definition
	G4ThreeVector offset;

	WorldSizeXY=80*cm;
	WorldSizeZ=80*cm;

	//     
	// World
	//

	solidWorld = new G4Box("World",WorldSizeXY/2,WorldSizeXY/2,WorldSizeZ/2);  
	logicWorld = new G4LogicalVolume(solidWorld,defaultMaterial,"World");		
	physiWorld = new G4PVPlacement(0,G4ThreeVector(0,0,0),"World",logicWorld,0,false,0, checkOverlaps);			

	//logicWorld->SetVisAttributes (G4VisAttributes::GetInvisible());

	//                                 
	// UP STRUCTURE
	//

	// fascia modello nuovo con anti a piastrelle
	std::string file_name;
	char* tmp_char;

	file_name.clear();
	//file_name = StructurePath + "DetectorStructure/Structure_V1/Structure/Mechanics.stl";
	//file_name = "/Users/aleksejsmirnov/Downloads/double/detector3/CRE_STRUCTURE.stl";
	file_name = "../detector3/CRE_STRUCTURE.stl";
	tmp_char  = const_cast<char*>(file_name.c_str());

	CADMesh * UpStructureMesh = new CADMesh(tmp_char, mm, offset, false);
	solidUpStructure = UpStructureMesh->TessellatedMesh( fCadUnit );
	logicUpStructure = new G4LogicalVolume(
			solidUpStructure,
			structureMaterial, 
			"UpStructure",0, 0, 0);
	logicUpStructure->SetVisAttributes(VisAtt_Structure);

	physiUpStructure = new G4PVPlacement(
			0,
			G4ThreeVector(),
			"UpStructure",
			logicUpStructure,
			physiWorld,
			false,	
			0,
			checkOverlaps);

	//
	// ------------UP PIXEL------------
	//

	for (int i=0; i<NbOfPixel; i++) {
    // for (int i=91; i<92; i++) {/DetectorStructure/Structure_V1/Cal
		file_name.clear();
		//file_name = StructurePath + "DetectorStructure/Structure_V1/Cal/Up"+std::to_string(i+1)+".stl";
	//	file_name = "/Users/aleksejsmirnov/Downloads/double/detector3/CRISTALLO_INF"+std::to_string(i+1)+".stl";
	//	file_name = "/Users/aleksejsmirnov/Downloads/double/detector3/CRISTALLO_SUP"+std::to_string(i+1)+".stl";
		file_name = "../detector3/CRISTALLO_SUP"+std::to_string(i+1)+".stl";
		tmp_char  = const_cast<char*>(file_name.c_str());

		CADMesh * UpPixelmesh = new CADMesh(tmp_char ,mm, offset, false);
		uppix = UpPixelmesh->TessellatedMesh( fCadUnit );
		logicUpPixel = new G4LogicalVolume(
				uppix,
				PixelMaterial,
				"UpPixel"+std::to_string(i),0, 0, 0);
		logicUpPixel->SetVisAttributes(VisAtt_Up);

		physiUpPixel = new G4PVPlacement(
				0,
				G4ThreeVector(),
				"UpPixel"+std::to_string(i),
				logicUpPixel,
				physiWorld,
				false,	
				i,
				checkOverlaps);

		/*
			 G4ThreeVector xyz = physiUpPixel->GetTranslation();
			 G4cout<<" ===== PixelID =  "<<i
			 <<"  X = "<<xyz.x()
			 <<"  Y = "<<xyz.y()
			 <<"  Z = "<<xyz.z()
			 <<"    ======"
			 <<G4endl;
		 */
	}

	//
	// ------------DOWN STRUCTURE ------------
	//

/*	file_name.clear();
	file_name = StructurePath + "DetectorStructure/Pixel_Down/07-FASCIA_INT.stl";
	tmp_char  = const_cast<char*>(file_name.c_str());

	CADMesh * DownStructureMesh = new CADMesh(tmp_char, mm, offset, false);
	solidDownStructure = DownStructureMesh->TessellatedMesh( fCadUnit );

	logicDownStructure = new G4LogicalVolume(
			solidDownStructure,
			structureMaterial, 
			"DownStructure",0, 0, 0);
	logicDownStructure->SetVisAttributes(VisAtt_Structure);
	physiDownStructure = new G4PVPlacement(0,
			G4ThreeVector(),
			"DownStructure",
			logicDownStructure,
			physiWorld,
			false,
			0,
			checkOverlaps);*/

	//
	// ------------DOWN PIXEL------------
	//

	for (int i=0; i<NbOfPixel; i++) {
	//	for (int i=109; i<110; i++) {

		file_name.clear();
		//file_name =  StructurePath + "DetectorStructure/Structure_V1/Cal/Dw"+std::to_string(i+1)+".stl";
		//file_name = "/Users/aleksejsmirnov/Downloads/double/detector3/CRISTALLO_INF"+std::to_string(i+1)+".stl";
		file_name = "../detector3/CRISTALLO_INF"+std::to_string(i+1)+".stl";

		tmp_char  = const_cast<char*>(file_name.c_str());

		CADMesh * DownPixelmesh = new CADMesh(tmp_char, mm, offset, false);
		downpix = DownPixelmesh->TessellatedMesh( fCadUnit );

		logicDownPixel = new G4LogicalVolume(
				downpix,
				PixelMaterial,
				"DownPixel"+std::to_string(i),0, 0, 0);
		logicDownPixel->SetVisAttributes(VisAtt_Down);

		physiDownPixel = new G4PVPlacement(
				0,
				G4ThreeVector(),
				"DownPixel"+std::to_string(i),
				logicDownPixel,
				physiWorld,
				false,
				NbOfPixel+i,
				checkOverlaps);// Please NOTE the CopyNo!!!
	}

	//                               
	// ---------  Anticoincidence (ACD)  -------
	//

	for (int i=0; i<1; i++) {
    //  for (int i=91; i<92; i++) {
		file_name.clear();
		//file_name = "/Users/aleksejsmirnov/Downloads/double/detector3/ACD1.stl";
		file_name = "../detector3/ACD1.stl";
		//file_name =  StructurePath + "DetectorStructure/Structure_V1/Veto/V"+std::to_string(i+1)+".stl";;
		tmp_char  = const_cast<char*>(file_name.c_str());

		CADMesh * ACDmesh = new CADMesh(tmp_char, mm, offset, false);
		solidACD = ACDmesh->TessellatedMesh( fCadUnit );

		logicACD = new G4LogicalVolume(
				solidACD,
				ACDMaterial,
				"ACD"+std::to_string(i),
				0, 0, 0);
		logicACD->SetVisAttributes(VisAtt_ACD);
		physiACD = new G4PVPlacement(
				0,
				G4ThreeVector(),
				"ACD"+std::to_string(i),
				logicACD,
				physiWorld,
				false,
				i,
				checkOverlaps);
	}
	
	/*G4Tubs* solidACD_Bottom = new G4Tubs(
			"solidBottomACD", // name
			0.*cm,            // innerRadius
			14.5*cm,          // outerRadius
			0.5*cm,           // half length
			0.*deg,           // starting Phi
			360.*deg          // segment angle
			);
	G4LogicalVolume* logicACD_Bottom = new G4LogicalVolume(
			solidACD_Bottom,
			ACDMaterial,
			"ACD"+std::to_string(NbOfPixel),
			0, 0, 0);
	logicACD_Bottom->SetVisAttributes(VisAtt_BottomACD);
	physiACD_Bottom = new G4PVPlacement(
			0,
			G4ThreeVector(0., 0., -0.5*cm),
			"ACD"+std::to_string(NbOfPixel),
			logicACD_Bottom,
			physiWorld,
			false,
			NbOfPixel,
			checkOverlaps);*/

	//
	//always return the physical World
	//
	return physiWorld;
}
 
//------------------------------------------------------------  

void CrystalEyeDetectorConstruction::ConstructSDandField()
{ 
	//  -------  Set Sensitive Detectors: Anticoincidence --------
	if( anticoincidenceSD.Get()==0 ){
		CrystalEyeAnticoincidenceSD* SD = new CrystalEyeAnticoincidenceSD("AnticoincidenceSD");
		anticoincidenceSD.Put(SD);
	}
	G4SDManager::GetSDMpointer()->AddNewDetector(anticoincidenceSD.Get());

	//for(int i=0; i<NbOfPixel; i++){
	for (int i=0; i<1; i++) {
		SetSensitiveDetector("ACD"+std::to_string(i),anticoincidenceSD.Get());
	}
	//SetSensitiveDetector("ACD"+std::to_string(NbOfPixel), anticoincidenceSD.Get());

	//  -------  Set Sensitive Detectors: Calorimeter --------
	if( calorimeterSD.Get()==0 ){
		CrystalEyeCalorimeterSD* SD = new CrystalEyeCalorimeterSD("CalorimeterSD");
		calorimeterSD.Put(SD);
	}
	G4SDManager::GetSDMpointer()->AddNewDetector( calorimeterSD.Get() );

	for(int i=0; i<NbOfPixel; i++){
	//for (int i=91; i<92; i++) {
		SetSensitiveDetector("DownPixel"+std::to_string(i), calorimeterSD.Get());
		SetSensitiveDetector("UpPixel"+std::to_string(i), calorimeterSD.Get());
	}
	//for (int i=109; i<110; i++) {
	//	SetSensitiveDetector("DownPixel"+std::to_string(i), calorimeterSD.Get());
	//}

	// Create global magnetic field messenger.
	// Uniform magnetic field is then created automatically if
	// the field value is not zero.
	G4ThreeVector fieldValue = G4ThreeVector();
	fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
	fMagFieldMessenger->SetVerboseLevel(1);

	// Register the field messenger for deleting
	G4AutoDelete::Register(fMagFieldMessenger);

	return;

}


void CrystalEyeDetectorConstruction::SetNbOfLayers(G4int val)
{
	NbOfLayers = val;
}

void CrystalEyeDetectorConstruction::SetNbOfPixel(G4int val)
{
	NbOfPixel = val;
}


void CrystalEyeDetectorConstruction::SetACDThickness(G4double val)
{
	ACDThickness = val;
}


void CrystalEyeDetectorConstruction::SetMagField(G4double fieldValue)
{
	// Just invoke manually the MT-safe command 
	// /globalField/setValue 
	// instantiated by the GlobalFieldMessenger
	std::stringstream sss;
	sss << "/globalField/setValue 0 0 " << fieldValue/tesla << " tesla";

	G4String command = sss.str();
	G4cout << "Going to execute: " << command << G4endl;

	G4UImanager* UImanager = G4UImanager::GetUIpointer();  
	UImanager->ApplyCommand(command);

}

//------------------------------------------------------------

void CrystalEyeDetectorConstruction::UpdateGeometry()
{
	//  delete payloadSD;
	G4RunManager::GetRunManager()->DefineWorldVolume(ConstructPayload());
	G4RunManager::GetRunManager()->PhysicsHasBeenModified();
	G4RegionStore::GetInstance()->UpdateMaterialList(physiWorld);

	G4RunManager::GetRunManager()->ReinitializeGeometry();

}

//------------------------------------------------------------

















