// **************************************************************
//
//      ------ CrystalEyeAnalysis  ------
//           Modified by Libo Wu (06 May 2022)
//           Email: libo.wu@gssi.it
//
// **************************************************************

#ifndef CrystalEyeDetectorConstruction_h
#define CrystalEyeDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4Cache.hh"
#include "globals.hh"

class G4Box;
class G4Tubs;
class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class G4Region;
class CrystalEyeDetectorMessenger;
class CrystalEyeAnticoincidenceSD;
class CrystalEyeCalorimeterSD;
class G4GlobalMagFieldMessenger;

//class CrystalEyeTrackerROGeometry;

//------------------------------------------------------------

class CrystalEyeDetectorConstruction : public G4VUserDetectorConstruction
{
		public:

				CrystalEyeDetectorConstruction();
				~CrystalEyeDetectorConstruction();

		public:

				void SetNbOfLayers (G4int); // CAL material, lenght, thickness
				void SetNbOfPixel (G4int);

				void SetACDThickness (G4double); //ACD Thickness

				void SetMagField(G4double); // Magnetic Field


				G4VPhysicalVolume* Construct();
				void UpdateGeometry();
				void ConstructSDandField();

		public:

				void PrintPayloadParameters();

				G4double GetWorldSizeZ()  const  {return WorldSizeZ;}; 
				G4double GetWorldSizeXY() const  {return WorldSizeXY;};

				G4int GetNbOfLayers() const      {return NbOfLayers;};
				G4int GetNbOfPixel()  const      {return NbOfPixel;};

				G4double GetBackgroundTimeU() const {return fBackgroundTimeU;};
				G4double GetBackgroundTimeD() const {return fBackgroundTimeD;};
				G4double GetThresholdTime() const {return fThresholdTime;};
				G4double GetACDThickness() const {return ACDThickness;};

		private:
				std::string StructurePath;
				bool checkOverlaps;
				G4double fCadUnit;

				G4int NbOfLayers;
				G4int NbOfPixel;
				G4double ACDThickness;
				G4double fBackgroundTimeD;
				G4double fBackgroundTimeU;
				G4double fThresholdTime;
				G4bool  fCheckOverlaps;

				G4Material*        defaultMaterial;
				G4Material*        PixelMaterial;
				G4Material*        ACDMaterial;
				G4Material*        structureMaterial;

				G4double           WorldSizeXY;
				G4double           WorldSizeZ;

				G4Box*             solidWorld;// World 
				G4LogicalVolume*   logicWorld;    
				G4VPhysicalVolume* physiWorld;    

				G4VSolid*          solidUpStructure;// Tracker
				G4LogicalVolume*   logicUpStructure;
				G4VPhysicalVolume* physiUpStructure;

				G4VSolid*          solidUpStructure2;// Tracker
				G4LogicalVolume*   logicUpStructure2;

				G4VSolid*          uppix;// Calorimeter
				G4LogicalVolume*   logicUpPixel;
				G4VPhysicalVolume* physiUpPixel;

				G4VSolid*          solidDownStructure;// Top Anticoincidence
				G4LogicalVolume*   logicDownStructure;
				G4VPhysicalVolume* physiDownStructure;

				G4VSolid*          downpix;// Lateral Anticoincidence
				G4LogicalVolume*   logicDownPixel;
				G4VPhysicalVolume* physiDownPixel;

				G4VSolid*          solidACD;
				G4LogicalVolume*   logicACD;
				G4VPhysicalVolume* physiACD;

				G4Tubs*            solidACD_Bottom;
				G4LogicalVolume*   logicACD_Bottom;
				G4VPhysicalVolume* physiACD_Bottom;


				// magnetic field messenger
				static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger; 

				CrystalEyeDetectorMessenger* detectorMessenger;  //pointer to the Messenger

				G4Cache<CrystalEyeAnticoincidenceSD*> anticoincidenceSD;
				G4Cache<CrystalEyeCalorimeterSD*>     calorimeterSD;


		private:

				void DefineMaterials();
				G4VPhysicalVolume* ConstructPayload();
};

//------------------------------------------------------------


#endif









