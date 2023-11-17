//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: CrystalEyeActionInitialization.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file CrystalEyeActionInitialization.cc
/// \brief Implementation of the CrystalEyeActionInitialization class

#include "CrystalEyeActionInitialization.hh"
#include "CrystalEyePrimaryGeneratorAction.hh"
#include "CrystalEyeRunAction.hh"
#include "CrystalEyeEventAction.hh"
#include "CrystalEyeTrackingAction.hh"
#include "CrystalEyeSteppingAction.hh"
#include "CrystalEyeAnticoincidenceDataManager.hh"
#include "CrystalEyeCalorimeterDataManager.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CrystalEyeActionInitialization::CrystalEyeActionInitialization()
 : G4VUserActionInitialization()
{
	fDataManagerAcd = new CrystalEyeAnticoincidenceDataManager();
	fDataManagerCal = new CrystalEyeCalorimeterDataManager();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CrystalEyeActionInitialization::~CrystalEyeActionInitialization()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CrystalEyeActionInitialization::Build() const
{
		CrystalEyeRunAction*      theRunAction = new CrystalEyeRunAction(
				fDataManagerAcd, fDataManagerCal);
		CrystalEyeEventAction*    theEvtAction = new CrystalEyeEventAction(
				theRunAction, fDataManagerAcd, fDataManagerCal);
		CrystalEyeTrackingAction*	theTrkAction = new CrystalEyeTrackingAction(theEvtAction);
		CrystalEyeSteppingAction* theStpAction = new CrystalEyeSteppingAction(theEvtAction);

  SetUserAction(new CrystalEyePrimaryGeneratorAction);
  SetUserAction(theRunAction);
	SetUserAction(theEvtAction);
	SetUserAction(theTrkAction);
	SetUserAction(theStpAction);
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CrystalEyeActionInitialization::BuildForMaster() const
{
  SetUserAction(new CrystalEyeRunAction(
				fDataManagerAcd,
				fDataManagerCal)
			);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

