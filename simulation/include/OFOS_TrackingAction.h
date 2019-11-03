#ifndef OFOS_TrackingInformation
#define OFOS_TrackingInformation 1

#include "G4UserTrackingAction.hh"
#include "globals.hh"

class OFOS_TrackingAction : public G4UserTrackingAction {

  public:  
    OFOS_TrackingAction( );
   ~OFOS_TrackingAction() {};
   
    virtual void PostUserTrackingAction(const G4Track*);
    
//private:
//  PrimaryGeneratorAction* fPrimary;
};

#endif
