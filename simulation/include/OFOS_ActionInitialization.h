
#ifndef OFOS_ActionInitialization_h
#define OFOS_ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "OFOS_DetectorConstruction.h"

//class B4DetectorConstruction;

/// Action initialization class.
///

class OFOS_ActionInitialization : public G4VUserActionInitialization
{
  public:
    explicit OFOS_ActionInitialization( OFOS_DetectorConstruction* det );
    virtual ~OFOS_ActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;

  private :
    OFOS_DetectorConstruction* fDetector;
};

#endif

    
