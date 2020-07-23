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
// $Id$
//
/// \file exampleB2b.cc
/// \brief Main program of the B2b example

#include "OFOS_DetectorConstruction.h"
#include "OFOS_ActionInitialization.h"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
//#include "FTFP_BERT.hh"
//#include "G4StepLimiterPhysics.hh"
#include "OFOS_PhysicsList.h"

#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4PhysListFactory.hh"

#include "OFOS_Verbosity.h"
#include "OFOS_OutputNtuples.h"
#include "OFOS_GlobalNtuplesPtr.h"


OFOS_OutputNtuples *global_ntuples_ptr;

int main(int argc,char** argv)
{
  // Detect interactive mode (if no arguments) and define UI session
  //


  OFOS_Verbosity::level = 0;
  global_ntuples_ptr = nullptr;


  G4UIExecutive* ui = nullptr;
  if ( argc == 1  || argc == 3 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // Construct the default run manager
  //

  auto* runManager = new G4RunManager;

  auto * det = new OFOS_DetectorConstruction();

  // Set mandatory initialization classes
  //
  runManager->SetUserInitialization(det);

//G4VModularPhysicsList* physicsList = new FTFP_BERT;
//physicsList->RegisterPhysics(new G4StepLimiterPhysics());

  /// customized physics list, inherited by Stefano
  auto *physicsList = new OFOS_PhysicsList;
  runManager->SetUserInitialization(physicsList);

  /// Now using standard physics list
  /// LBE seems the be the only suitable to simulate optical photons and radioactive decays
  /// cf geant4.in2p3.fr/IMG/pdf_PhysicsLists.pdf
//G4PhysListFactory *physListFactory = new G4PhysListFactory();
//G4VUserPhysicsList *physicsList = physListFactory->GetReferencePhysList("LBE");
//runManager->SetUserInitialization(physicsList);


  // Set user action classes
  runManager->SetUserInitialization(new OFOS_ActionInitialization(det));

  // >>> MY DEBUG <<<<
//G4cout << "runManager->InitializeGeometry()" << G4endl;
//runManager->InitializeGeometry();

//G4cout << "runManager->InitializePhysics()" << G4endl;
//runManager->InitializePhysics();

//G4cout << "runManager->Initialize()" << G4endl;
//runManager->Initialize();

//G4cout << G4endl << G4endl;





  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive ();
  if(OFOS_Verbosity::level==2)
      visManager->SetVerboseLevel("Quiet");

  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  G4cout << "Initialize visual manager" << G4endl;
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( argc == 1 ) {
    // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    if (ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }
    ui->SessionStart();
    delete ui;
  }
  if ( argc == 2 ) {
    // batch mode
    G4cout << "Exectute " << argv[1] << G4endl;
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
    G4cout << argv[1] << " Execution : Done" << G4endl;
  }
  if (argc == 3 ) {
    // show interactive mode for the configuration specified by the macro
    G4cout << "Exectute " << argv[1] << G4endl;
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
    // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    if (ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !

  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
