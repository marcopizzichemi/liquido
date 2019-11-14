

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4TouchableHandle.hh"
#include "G4StepPoint.hh"
#include "G4VProcess.hh"
#include "G4OpticalProcessIndex.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmProcessSubType.hh"

#include "OFOS_OpticalPhotonSD.h"
#include "OFOS_Verbosity.h"
#include "OFOS_OutputLog.h"



OFOS_OpticalPhotonSD::OFOS_OpticalPhotonSD ( const G4String& name, const G4String& hitCollectionName,  SDType type) :
        G4VSensitiveDetector(name),
        hit_collection_(nullptr),
        type_(type),
        sd_name_(name),
        eV_to_nm_(1239.84193) // nm
{
    collectionName.insert(hitCollectionName);


    G4cout << "OFOS_OpticalPhotonSD :: name(" << name << ")  collection(" << hitCollectionName << ")"
           << "   type (" << (G4int)type_ << ")" << G4endl;

    OFOS_OutputLog::log_cache << name << " id: " << (G4int)type_ << "  collection: " << hitCollectionName << G4endl;
}


/***********************************************************************************/
/***********************************************************************************/


OFOS_OpticalPhotonSD::~OFOS_OpticalPhotonSD()
{}



/***********************************************************************************/
/***********************************************************************************/


void OFOS_OpticalPhotonSD::Initialize(G4HCofThisEvent* event_container)
{
    // Create hits collection
    hit_collection_ = new OFOS_OPHitCollection(SensitiveDetectorName, collectionName[0]);

    /// Get collection ID
    G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);

    // Add collection to event_container
    event_container->AddHitsCollection( hcID, hit_collection_ );

    if( OFOS_Verbosity::level > 1 )
        G4cout << "OFOS_OpticalPhotonSD::Initialize()  "
               << "SD Name: "   << SensitiveDetectorName << "   "
               << "Coll Name: " << collectionName[0] << "   "
               << "Coll ID: "   <<  hcID << G4endl;
}



/***********************************************************************************/
/***********************************************************************************/


G4bool OFOS_OpticalPhotonSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    // if( OFOS_Verbosity::level > 2 )
    //     G4cout << "OFOS_OpticalPhotonSD::ProcessHits >>> DOING NOTHING" << G4endl;

    return false;
}


/***********************************************************************************/
/***********************************************************************************/


G4bool OFOS_OpticalPhotonSD::process_hit( const G4Step* aStep, const G4TouchableHistory* )
{
    if( OFOS_Verbosity::level > 3 )
    {
        G4cout << "OFOS_OpticalPhotonSD::process_hit >>> CREATING HIT" << G4endl;
        //  G4cout << "wavelength = " << eV_to_nm_ / (thePostPoint->GetTotalEnergy()/eV) << " nm" << G4endl;
    }

    /// create new hit
    auto aHit = new OFOS_OPHit();


    /// fetch hit properties
    G4Track*      theTrack         = aStep->GetTrack();
    G4StepPoint*  thePrePoint      = aStep->GetPreStepPoint();
    G4StepPoint*  thePostPoint     = aStep->GetPostStepPoint();
    G4ThreeVector polarization     = thePostPoint->GetPolarization();
    G4ThreeVector position         = thePostPoint->GetPosition();
    G4int         gen_proc_type    = (theTrack->GetCreatorProcess() ? theTrack->GetCreatorProcess()->GetProcessType()    : my_undefined_proc );
    G4int         gen_proc_subtype = (theTrack->GetCreatorProcess() ? theTrack->GetCreatorProcess()->GetProcessSubType() : my_undefined_proc );



    /// needed to know what is the volume being hit
    const G4TouchableHandle& th = thePostPoint->GetTouchableHandle();

    /// IDs to be stored in the Hit
    G4int pid = 0; // primary
    G4int sid = 0; // secondary
    G4int my_proc_id = my_undefined_proc;

    /// determine the ID of the volume being hit
    if (this->type_ & SDType::Fiber)
    {
        pid = th->GetCopyNumber(4); // readout_unit_id
        sid = th->GetCopyNumber(3); // fiber_id
    }
    else if( this->type_ & SDType::Pmt)
    {
        pid = th->GetCopyNumber(3); // readout_unit_id
        sid = 0; /// needed to implement top/bottom separation
    }

    /// determine what process created the optical photon being processed
    if(gen_proc_type == G4ProcessType::fElectromagnetic)
    {
        switch (gen_proc_subtype)
        {
            // case g4_scintillation : my_proc_id = my_scintillation; break;
            // case g4_cherenkov     : my_proc_id = my_cherenkov;     break;
            case G4EmProcessSubType::fScintillation : my_proc_id = my_scintillation; break;
            case G4EmProcessSubType::fCerenkov      : my_proc_id = my_cherenkov;     break;
            default: G4cout << "Unknown optical photon creation process" << G4endl;  break;
        }
    }

    /// fill hit
    aHit->set_track_id     (theTrack->GetTrackID());
//    aHit->set_parent_id    (theTrack->GetParentID());  // Josh addition
    aHit->set_sd_type      ( static_cast<G4int>(type_) );
    aHit->set_primary_id   ( pid );
    aHit->set_secondary_id ( sid );
    aHit->set_gen_proc     ( my_proc_id );
    aHit->set_time         ( thePostPoint->GetGlobalTime() );
    aHit->set_wavelength   ( eV_to_nm_ / (thePostPoint->GetTotalEnergy()/eV) );
    aHit->set_position     ( thePostPoint->GetPosition() );
    aHit->set_polarization ( thePostPoint->GetPolarization() );

    /// store hit
    hit_collection_->insert(aHit);

    return true;

}


/***********************************************************************************/
/***********************************************************************************/



void OFOS_OpticalPhotonSD::EndOfEvent(G4HCofThisEvent*)
{
    if (OFOS_Verbosity::level>1)
        G4cout << "OFOS_OpticalPhotonSD::EndOfEvent() :: number of hits = " << hit_collection_->entries() << G4endl;
}

