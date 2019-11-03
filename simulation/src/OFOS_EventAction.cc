#include <vector>

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

#include "g4root.hh" 

#include "OFOS_EventAction.h"
//#include "OFOS_SiPmtHit.h"
//#include "OFOS_FiberHit.h"
//#include "OFOS_VesselHit.h"
#include "OFOS_Verbosity.h"
//#include "OFOS_FiberSD.h"
#include "OFOS_OpticalPhotonSD.h"
#include "OFOS_OutputNtuples.h"
#include "OFOS_GlobalNtuplesPtr.h" 



OFOS_EventAction::OFOS_EventAction() : G4UserEventAction(),
                                       vert_fiber_hit_collection_id_ (-1),
                                       hori_fiber_hit_collection_id_ (-1),
                                       vert_sipmt_hit_collection_id_ (-1),
                                       hori_sipmt_hit_collection_id_ (-1),
                                       vessel_hit_collection_id_ (-1)
{
}



OFOS_EventAction::~OFOS_EventAction()
{
    if (OFOS_Verbosity::level>1) 
        G4cout << "~OFOS_EventAction()" << G4endl;
}



/***********************************************************************************/
/***********************************************************************************/
 



void OFOS_EventAction::BeginOfEventAction(const G4Event*)
{
    if ( OFOS_Verbosity::level>1) 
        G4cout << "Begin of Event Action" << G4endl;
}



/***********************************************************************************/
/***********************************************************************************/
 



void 
OFOS_EventAction::EndOfEventAction(const G4Event* event)
{
    
    G4int eventID = event->GetEventID();

  //if( eventID < 100 || eventID%100 == 0) 
    if ( OFOS_Verbosity::level>1) 
        G4cout << "EndOfEventAction  ev:" << eventID << G4endl;
    

    /// fetch collection IDs
    if ( vert_fiber_hit_collection_id_ == -1 ) vert_fiber_hit_collection_id_ = G4SDManager::GetSDMpointer()->GetCollectionID("VertFiberHitsCollection");
    if ( hori_fiber_hit_collection_id_ == -1 ) hori_fiber_hit_collection_id_ = G4SDManager::GetSDMpointer()->GetCollectionID("HoriFiberHitsCollection");
    if ( vessel_hit_collection_id_ == -1 )     vessel_hit_collection_id_     = G4SDManager::GetSDMpointer()->GetCollectionID("VesselHitsCollection");
    if ( vert_sipmt_hit_collection_id_ == -1 ) vert_sipmt_hit_collection_id_ = G4SDManager::GetSDMpointer()->GetCollectionID("VertSiPmtHitsCollection");
    if ( hori_sipmt_hit_collection_id_ == -1 ) hori_sipmt_hit_collection_id_ = G4SDManager::GetSDMpointer()->GetCollectionID("HoriSiPmtHitsCollection");

    process_hit_collections( event );
    fetch_event_info(event);
    global_ntuples_ptr->store_event();

}





/***********************************************************************************/
/***********************************************************************************/
 



void
OFOS_EventAction::fetch_event_info( const G4Event* evt )
{
    /// TODO
    //  To be implemented
}




/***********************************************************************************/
/***********************************************************************************/
 




G4int 
OFOS_EventAction::process_hit_collections ( const G4Event* evt )
{
    if ( OFOS_Verbosity::level>1) 
        G4cout << "EventAction::process_hit_collections()" << G4endl;

    /// setup container of hit collections
    G4HCofThisEvent* HCofEvent = evt->GetHCofThisEvent();

    const int n_collections = 5;

    OFOS_OPHitCollection* collection_container[n_collections];


    /// fetch hit collections 
    collection_container[0] = (OFOS_OPHitCollection*) (HCofEvent->GetHC(vert_fiber_hit_collection_id_));
    collection_container[1] = (OFOS_OPHitCollection*) (HCofEvent->GetHC(hori_fiber_hit_collection_id_));
    collection_container[2] = (OFOS_OPHitCollection*) (HCofEvent->GetHC(vessel_hit_collection_id_    ));
    collection_container[3] = (OFOS_OPHitCollection*) (HCofEvent->GetHC(vert_sipmt_hit_collection_id_));
    collection_container[4] = (OFOS_OPHitCollection*) (HCofEvent->GetHC(hori_sipmt_hit_collection_id_));

    /// determine collection size
    G4int n_vert_fiber_hits      = collection_container[0]->entries();
    G4int n_hori_fiber_hits      = collection_container[1]->entries();
    G4int n_vessel_hits          = collection_container[2]->entries();
    G4int n_vert_sipmt_hits      = collection_container[3]->entries();
    G4int n_hori_sipmt_hits      = collection_container[4]->entries();
    G4int n_tot_hits             = n_vert_fiber_hits + n_hori_fiber_hits + n_vessel_hits + n_vert_sipmt_hits + n_hori_sipmt_hits ; 


    if ( OFOS_Verbosity::level>0) 
    {
        G4cout << "EndOfEventAction :: n_vert_fiber_hits = " << n_vert_fiber_hits << G4endl; 
        G4cout << "EndOfEventAction :: n_hori_fiber_hits = " << n_hori_fiber_hits << G4endl; 
        G4cout << "EndOfEventAction :: n_vess_fiber_hits = " << n_vessel_hits << G4endl; 
    }

    OFOS_OPHit *a_hit;

    for (auto & iC : collection_container)
    {
        int n_hits = iC->entries();

        for(int iH=0; iH<n_hits; ++ iH)
        {
            a_hit = (*iC)[iH];
            global_ntuples_ptr->fill_hit( a_hit );

        }
    }
    return n_tot_hits;
}  


