#ifndef OFOS_EventAction_h
#define OFOS_EventAction_h 1

#include <vector>

#include "G4UserEventAction.hh"

#include "globals.hh"

//#include "OFOS_SiPmtHit.h"
//#include "OFOS_FiberHit.h"



class OFOS_EventAction : public G4UserEventAction
{

public:
  OFOS_EventAction();
  virtual ~OFOS_EventAction();
  

  // data members
   G4int  vert_fiber_hit_collection_id_;
   G4int  hori_fiber_hit_collection_id_;
   G4int  vert_sipmt_hit_collection_id_;
   G4int  hori_sipmt_hit_collection_id_;
   G4int  vessel_hit_collection_id_;

  virtual void  BeginOfEventAction(const G4Event* );
  virtual void    EndOfEventAction(const G4Event* );
  

private:
 
    G4int process_hit_collections ( const G4Event* evt );
    static void fetch_event_info( const G4Event* evt );

};



#endif
