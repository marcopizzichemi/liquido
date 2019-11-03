
#include "OFOS_ActionInitialization.h"
#include "OFOS_PrimaryGeneratorAction.h"
#include "OFOS_SteppingAction.h"
#include "OFOS_RunAction.h"
#include "OFOS_EventAction.h"
#include "OFOS_TrackingAction.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OFOS_ActionInitialization::OFOS_ActionInitialization(OFOS_DetectorConstruction* det)
 : G4VUserActionInitialization(),
    fDetector(det)
{
G4cout << "OFOS_ActionInitialization::" <<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OFOS_ActionInitialization::~OFOS_ActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OFOS_ActionInitialization::BuildForMaster() const
{
  SetUserAction(new OFOS_RunAction(fDetector));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OFOS_ActionInitialization::Build() const
{
    G4cout << "OFOS_ActionInitialization::Build()" << G4endl;

/// OFOS_PrimaryGeneratorAction is needed but parameters are overridden by mac file
/// Possibly the way to retrieve them is through a messenger class 
    SetUserAction(new OFOS_PrimaryGeneratorAction);
                                                   
/// Needed to init root file
    SetUserAction(new OFOS_RunAction(fDetector) ); 
    SetUserAction(new OFOS_EventAction);
    SetUserAction(new OFOS_SteppingAction(fDetector));
    SetUserAction(new OFOS_TrackingAction);
}  

