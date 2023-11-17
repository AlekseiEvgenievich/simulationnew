// **************************************************************
//
//      ------ CrystalEyePrimaryGeneratorAction.cc  ------
//           by Libo Wu (23 Sep 2021)
//           Email: libo.wu@gssi.it
//
// **************************************************************

#include "G4RunManager.hh"
#include "CrystalEyePrimaryGeneratorAction.hh"

#include "CrystalEyeDetectorConstruction.hh"
#include "CrystalEyePrimaryGeneratorMessenger.hh"
#include "CrystalEyeAnalysis.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

//------------------------------------------------------------
CrystalEyePrimaryGeneratorAction::CrystalEyePrimaryGeneratorAction()
{
	fDetector = static_cast<const CrystalEyeDetectorConstruction*>
		(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

	//---  Create a messenger for this class  ---
	gunMessenger = new CrystalEyePrimaryGeneratorMessenger(this);

	rndmFlag = "off";
	nSourceType = 0;
	nSpectrumType = 0;
	sourceGun = false;
	dVertexRadius = 15.*cm; 

	G4int n_particle = 1;

	fParticleGun  = new G4ParticleGun(n_particle);     
	// default particle kinematic

	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;
	G4ParticleDefinition* particle = particleTable->FindParticle(particleName="e-");
	//fParticleGun->SetParticleDefinition(particle);
	//fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
	//fParticleGun->SetParticleEnergy(1.*MeV);
	//G4double position = 0.5*(fDetector->GetWorldSizeZ());
	//fParticleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,position));

	fParticleSource = new G4GeneralParticleSource();

	//fPosDist = new G4SPSPosDistribution();
	//fPosDist = fParticleSource->GetCurrentSource()->GetPosDist();
	//fAngDist = new G4SPSAngDistribution();
	//fAngDist = fParticleSource->GetCurrentSource()->GetAngDist();
	fEngDist = new G4SPSEneDistribution();
	fEngDist = fParticleSource->GetCurrentSource()->GetEneDist();

	//default setting
	//fPosDist->SetPosDisType("Point");
	//fPosDist->SetCentreCoords( G4ThreeVector(0, 0, 20) );
	//fAngDist->SetParticleMomentumDirection( G4ThreeVector(0., 0., -1));
	//fEngDist->SetEnergyDisType("Mono");
	//fEngDist->SetMonoEnergy(1.*MeV);
}

//------------------------------------------------------------

CrystalEyePrimaryGeneratorAction::~CrystalEyePrimaryGeneratorAction()
{

	delete fParticleGun;
	delete fParticleSource;
	delete gunMessenger;
}

/*double CrystalEyePrimaryGeneratorAction::Func_DiffuseGamma(double eng) {
    double f;
    if(eng <= 60.0) {
        f = 7.877 * std::pow(eng, -1.29) * std::exp(-eng / 41.13);
    } else {
        f = 4.32e-4 * std::pow(eng / 60.0, -6.5)
            + 8.4e-3 * std::pow(eng / 60.0, -2.58)
            + 4.8e-4 * std::pow(eng / 60.0, -2.05);
    }
    return f;
}

// Function to generate a particle with energy according to Func_DiffuseGamma
double CrystalEyePrimaryGeneratorAction::GenerateParticleEnergy() {
    double eng_min = 1.0; // Minimum energy in keV
    double eng_max = 100.0; // Maximum energy in keV
    double max_pdf = 8; // Calculate the maximum value of Func_DiffuseGamma in the range

    while (true) {
        double eng = G4UniformRand() * (eng_max - eng_min) + eng_min; // Random energy
        double pdf_val = Func_DiffuseGamma(eng);
        if (G4UniformRand() * max_pdf < pdf_val) {
            return eng * keV; // Convert to Geant4 units and return
        }
    }
}

void CrystalEyePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

    // Generate energy using the Func_DiffuseGamma distribution
    double energy = GenerateParticleEnergy();
    fParticleGun->SetParticleEnergy(energy);

    // Generate random direction within the upper hemisphere
    double theta = std::acos(G4UniformRand()); // Polar angle between 0 and π/2
    double phi = 2.0 * M_PI * G4UniformRand(); // Azimuthal angle between 0 and 2π

    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(std::sin(theta) * std::cos(phi), 
                                                            std::sin(theta) * std::sin(phi), 
                                                            std::cos(theta)));

    fParticleGun->GeneratePrimaryVertex(anEvent);
}*/


//------------------------------------------------------------

void CrystalEyePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	//---  For Particle Gun  ---
	if (sourceGun)
	{
		G4double z0 = 0.5*(fDetector->GetWorldSizeZ());
		G4double x0 = 0.*cm, y0 = 0.*cm;

		G4ThreeVector pos0;
		G4ThreeVector dir0;
		G4ThreeVector vertex0 = G4ThreeVector(x0,y0,z0);

		dir0 = G4ThreeVector(0.,0.,-1.);

		G4double theta, phi, y, f;
		G4double theta0=0.;
		G4double phi0=0.;

		switch(nSourceType) {
			case 0:
				fParticleGun->SetParticlePosition(vertex0);
				fParticleGun->SetParticleMomentumDirection(dir0);
				break;

			case 1:
				// GS: Generate random position on the 4PIsphere to create a unif. distrib.
				// GS: on the sphere
				phi = G4UniformRand() * twopi;
				do {
					y = G4UniformRand()*1.0;
					theta = G4UniformRand() * pi;
					f = std::sin(theta);
				} while (y > f);
				vertex0 = G4ThreeVector(1.,0.,0.);
				vertex0.setMag(dVertexRadius);
				vertex0.setTheta(theta);
				vertex0.setPhi(phi);
				fParticleGun->SetParticlePosition(vertex0);

				dir0 = G4ThreeVector(1.,0.,0.);
				do {
					phi = G4UniformRand() * twopi;
					do {
						y = G4UniformRand()*1.0;
						theta = G4UniformRand() * pi;
						f = std::sin(theta);
					} while (y > f);
					dir0.setPhi(phi);
					dir0.setTheta(theta);
				} while (vertex0.dot(dir0) >= -0.7 * vertex0.mag());
				fParticleGun->SetParticleMomentumDirection((G4ParticleMomentum)dir0);
				break;

			case 2:
				// GS: Generate random position on the upper semi-sphere z>0 to create a unif. distrib.
				// GS: on a plane
				phi = G4UniformRand() * twopi;
				do {
					y = G4UniformRand()*1.0;
					theta = G4UniformRand() * halfpi;
					f = std::sin(theta) * std::cos(theta);
				} while (y > f);
				vertex0 = G4ThreeVector(1.,0.,0.);

				G4double xy = fDetector->GetWorldSizeXY();
				G4double z  = fDetector->GetWorldSizeZ();

				if (dVertexRadius > xy*0.5)
				{ 
					G4cout << "vertexRadius too big " << G4endl;
					G4cout << "vertexRadius setted to " << xy*0.45 << G4endl;
					dVertexRadius = xy*0.45;
				}

				if (dVertexRadius > z*0.5)
				{ 
					G4cout << "vertexRadius too high " << G4endl;
					G4cout << "vertexRadius setted to " << z*0.45 << G4endl;
					dVertexRadius = z*0.45;
				}

				vertex0.setMag(dVertexRadius);
				vertex0.setTheta(theta);
				vertex0.setPhi(phi);

				// GS: Get the user defined direction for the primaries and
				// GS: Rotate the random position according to the user defined direction for the particle

				dir0 = fParticleGun->GetParticleMomentumDirection();
				if (dir0.mag() > 0.001) 
				{
					theta0 = dir0.theta();
					phi0   = dir0.phi();   
				}

				if (theta0!=0.) 
				{
					G4ThreeVector rotationAxis(1.,0.,0.);
					rotationAxis.setPhi(phi0+halfpi);
					vertex0.rotate(theta0+pi,rotationAxis);
				}
				fParticleGun->SetParticlePosition(vertex0);
				break;
		}//Source Type

		G4double pEnergy;
		switch(nSpectrumType) {
			case 0:
				fEngDist->SetEnergyDisType("Pow");
				fEngDist->SetAlpha(-1);
				fEngDist->SetEmin(10*keV);
				fEngDist->SetEmax(100*MeV);
				break;

			case 1:
				fEngDist->SetEnergyDisType("Cdg");
				fEngDist->SetEmin(10*keV);
				fEngDist->SetEmax(100*MeV);
				break;

			case 2:
				do {
					y = G4UniformRand()*100000.0;
					pEnergy = G4UniformRand() * 10. * GeV;
					f = std::pow(pEnergy * (1/GeV), -4.);
				} while (y > f);

				fParticleGun->SetParticleEnergy(pEnergy);
				break;

			case 3:
				break;
		}//Spectrum Type

		fParticleGun->GeneratePrimaryVertex(anEvent);
	}
	//---  For Particle Source  ---
	else
	{
		fParticleSource->GeneratePrimaryVertex(anEvent);

		G4int		      pdg			  = fParticleSource->GetParticleDefinition()->GetPDGEncoding();
		G4double      priEnergy = fParticleSource->GetParticleEnergy();
		G4ThreeVector direction = fParticleSource->GetParticleMomentumDirection();
		G4ThreeVector position  = fParticleSource->GetParticlePosition();

		CrystalEyeAnalysis* fanalysis = CrystalEyeAnalysis::getInstance();
		fanalysis -> SetPrimaryData(
				pdg, 
				priEnergy, 
				direction,
			 	position);

	}

}

//------------------------------------------------------------
//------------------------------------------------------------
//------------------------------------------------------------








