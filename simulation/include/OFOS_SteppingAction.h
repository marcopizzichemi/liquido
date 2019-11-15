
#ifndef OFOS_SteppingAction_h
#define OFOS_SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "OFOS_DetectorConstruction.h"

class OFOS_DetectorConstruction;
//class OFOS_SteppingActionMessenger;

class G4Track;
class G4StepPoint;
class G4OpBoundaryProcess;

//enum g4_subprocesses{ OpAbsorption=31, Compton=13, Photoelectric=12, Bremsstrahlung=3, Cherenkov=21, Scintillation=22 };

class OFOS_SteppingAction : public G4UserSteppingAction
{
  public:

    explicit OFOS_SteppingAction(OFOS_DetectorConstruction*);
    virtual ~OFOS_SteppingAction();

    virtual void UserSteppingAction(const G4Step*);
 
  private:

    // maximum number of save states
    static G4int fMaxRndmSave;
 
    OFOS_DetectorConstruction* fDetector;
    

};

#endif
