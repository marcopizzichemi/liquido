#include <sstream>
#include <cmath>

#include "G4Run.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ProcessManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4Scintillation.hh"
#include "G4Cerenkov.hh"

#include "G4OpProcessSubType.hh"

#include "OFOS_GlobalNtuplesPtr.h"
#include "OFOS_SteppingAction.h"
#include "OFOS_DetectorConstruction.h"
#include "OFOS_Verbosity.h"
#include "OFOS_OpticalPhotonSD.h"

#include "G4EmProcessSubType.hh"
#include "G4SteppingManager.hh"
#include "G4EventManager.hh"

#include "TString.h"

//G4int OFOS_SteppingAction::fMaxRndmSave = 10000;


OFOS_SteppingAction::OFOS_SteppingAction(OFOS_DetectorConstruction *detector) :
        fDetector(detector) {

}


/***********************************************************************************/
/***********************************************************************************/



OFOS_SteppingAction::~OFOS_SteppingAction() = default;


/***********************************************************************************/
/***********************************************************************************/



void OFOS_SteppingAction::UserSteppingAction(const G4Step *theStep) {
    G4Track *theTrack = theStep->GetTrack();
    G4ParticleDefinition *particleDefinition = theTrack->GetDefinition();
    auto PDGEEncoding = particleDefinition->GetPDGEncoding();
    auto trackID = theTrack->GetTrackID();

    /// post_step is where absorption takes place
    G4StepPoint *thePrePoint = theStep->GetPreStepPoint();
    G4StepPoint *thePostPoint = theStep->GetPostStepPoint();

    G4VPhysicalVolume *thePrePV = thePrePoint->GetPhysicalVolume();
    G4VPhysicalVolume *thePostPV = thePostPoint->GetPhysicalVolume();

    G4String thePrePVname = "";
    G4String thePostPVname = "";
    G4int n_scin_phot = 0;
    G4int n_cher_phot = 0;
    G4String post_step_proc_list = "";
    G4String along_step_proc_list = "";


    thePrePVname = thePrePV->GetName();

    thePostPVname = thePostPV->GetName();


    G4SteppingManager *fpSteppingManager = G4EventManager::GetEventManager()->GetTrackingManager()->GetSteppingManager();
    G4StepStatus stepStatus = fpSteppingManager->GetfStepStatus();


    if (stepStatus == fWorldBoundary) {
        theTrack->SetTrackStatus(fStopAndKill);
        G4cout << G4endl << "*** WARNING :: Particle touching World Volume >>> Killing Track " << G4endl;
        return;
    }
    if (thePostPVname.data()[0] == 'B' // shortcut for Buffer
        && thePrePVname == "LsVessel"  // only if particle is coming from inside liquido vessel (i.e. pointing outwards)
        )
    {
        theTrack->SetTrackStatus(fStopAndKill);
        global_ntuples_ptr->add_lost_phot(1);
        global_ntuples_ptr->set_lost_containment();
        G4cout << G4endl << "*** WARNING :: Particle in Buffer >>> Killing Track " << G4endl;
        return;
    }


    /// Steps with status fAtRestDoItProc have zero length
    /// and zero deposited energy. Not clear how to use them.
    /// In this  an example is Scintillation with
    /// zero emitted photons;
//  if (stepStatus == fAtRestDoItProc) return;


    /// Scintillation and Cherenkov are processes which do not play a role
    /// in defining the step length (no max number of photons specified)
    /// in the physics list. As a consequence they cannot be fetched through
    /// the post point and must be retrieved from the vector of postStep
    /// processes
    size_t MAXofPostStepLoops = fpSteppingManager->GetMAXofPostStepLoops();
    G4ProcessVector *procPost = fpSteppingManager->GetfPostStepDoItVector();


    /// Loop over all the postStep processes available for this particle
    /// Avoid picking up the scintilation process defined at rest
    ///       since it does not contain the number of generated
    ///       scintillation photons. The true number of scintillation
    ///       photons is defined at the postStep level.
    for (size_t iProc = 0; iProc < MAXofPostStepLoops && stepStatus != fAtRestDoItProc; iProc++) {
        G4ProcessType proc_type = (*procPost)[iProc]->GetProcessType();
        G4int proc_subtype = (*procPost)[iProc]->GetProcessSubType();

        /// debug output
        if (OFOS_Verbosity::level > 2) {
            post_step_proc_list += G4String(TString::Format(">p>  [%d] %x ", iProc, (*procPost)[iProc]).Data());
            post_step_proc_list += (*procPost)[iProc]->GetProcessName();
            post_step_proc_list += G4String(TString::Format(" (%d / %d)\n",
                                                            static_cast<G4int>(proc_type),
                                                            proc_subtype).Data());
        }

        /// not able to store photons in output storage >>> do nothing
        if (!global_ntuples_ptr) continue;

        /// Scintillation and Cherenkov light are emitted only in the case of a EM process
        if (proc_type != G4ProcessType::fElectromagnetic) continue;

        switch (proc_subtype) {
            case G4EmProcessSubType::fScintillation : {
                auto *proc1 = dynamic_cast<G4Scintillation *>((*procPost)[iProc] );
                n_scin_phot += proc1->GetNumPhotons();
                global_ntuples_ptr->add_scint_phot(proc1->GetNumPhotons());
                break;
            }
            case G4EmProcessSubType::fCerenkov : {
                auto *proc1 = dynamic_cast<G4Cerenkov *>((*procPost)[iProc]);
                n_cher_phot += proc1->GetNumPhotons();
                global_ntuples_ptr->add_chere_phot(proc1->GetNumPhotons());
                break;
            }
            default:
                break;
        }
    }


    /// get process which defined the step length
    G4String pp_proc_name = thePostPoint->GetProcessDefinedStep()->GetProcessName();
    G4ProcessType pp_proc_type = thePostPoint->GetProcessDefinedStep()->GetProcessType();
    G4int pp_proc_subtype = thePostPoint->GetProcessDefinedStep()->GetProcessSubType();


    /// Check if a sensitive detector is being crossed
    /// Needed in order to store hits
    G4VSensitiveDetector *thisSD = thePostPoint->GetSensitiveDetector();
    G4int sd_type = -9;
    if (thisSD)
        sd_type = static_cast<G4int>( dynamic_cast<OFOS_OpticalPhotonSD *>(thisSD)->get_sd_type());


    /// check material being crossed
    if (OFOS_Verbosity::level > 3) {
        G4cout << G4endl << "pre: " << thePrePVname << "  >>  post: " << thePostPVname << G4endl;
        // G4cout << "time: " << theTrack->GetGlobalTime()/ns << G4endl ;

        G4Material *current_material = thePrePoint->GetMaterial();

        /// check if optical properties of the material being crossed are correct
        if (current_material && thePrePVname == "LS") {
            G4cout << current_material->GetName();
            G4cout << "  ( a: "
                   << current_material->GetMaterialPropertiesTable()->GetProperty("ABSLENGTH")->GetMaxValue() / mm
                   << " mm ";
            G4cout << "    s: "
                   << current_material->GetMaterialPropertiesTable()->GetProperty("RAYLEIGH")->GetMaxValue() / mm
                   << " mm )";
            G4cout << G4endl;

            //G4cout << G4endl  << "----------------" << G4endl;
            //current_material->GetMaterialPropertiesTable()->DumpTable();
            //G4cout << "----------------" << G4endl;
        }
    }


    if (OFOS_Verbosity::level > 2 && pp_proc_type != G4ProcessType::fTransportation &&
        particleDefinition->GetPDGEncoding() != 0)
        // if(pp_proc_type != G4ProcessType::fTransportation)
        // if (OFOS_Verbosity::level>2 && pp_proc_type == G4ProcessType::fElectromagnetic) //&& pp_proc_type != G4ProcessType::fTransportation  && pp_proc_type != G4ProcessType::fOptical)
        //if ( (OFOS_Verbosity::level>2 && pp_proc_type != G4ProcessType::fTransportation && pp_proc_type != G4ProcessType::fOptical)
        //     || (thePostPVname.data()[0] == 'B') ) /// Buffer
    {

        const std::vector<const G4Track *> *sec_tracks = theStep->GetSecondaryInCurrentStep();
        G4double sec_energy(0.0);
        G4String list_of_secondaries = "";
        for (auto sec_track : *sec_tracks) {
            if (sec_track->GetDefinition()->GetPDGEncoding() != 0) {
                sec_energy += sec_track->GetKineticEnergy();
                list_of_secondaries += "sec: trk_id:";
                list_of_secondaries += G4String(TString::Format("%d", sec_track->GetTrackID()).Data());
                list_of_secondaries += " (";
                list_of_secondaries += G4String(
                        TString::Format("%d", sec_track->GetDefinition()->GetPDGEncoding()).Data());
                list_of_secondaries += ")  ";
                list_of_secondaries += "mom_id: ";
                list_of_secondaries += G4String(TString::Format("%d\n", sec_track->GetParentID()).Data());
            }
        }


        G4cout << G4endl;

        G4cout << "particle name: " << particleDefinition->GetParticleName() << G4endl;
        G4cout << "pdgid: " << particleDefinition->GetPDGEncoding() << G4endl;
        G4cout << "Track ID: " << theTrack->GetTrackID() << G4endl;
        G4cout << "Parent ID: " << theTrack->GetParentID() << G4endl;
        G4cout << "Prepoint: (" << thePrePoint->GetPosition().x() << ", "
               << thePrePoint->GetPosition().y() << ", "
               << thePrePoint->GetPosition().z() << ") " << G4endl;
        G4cout << "Postpoint: (" << thePostPoint->GetPosition().x() << ", "
               << thePostPoint->GetPosition().y() << ", "
               << thePostPoint->GetPosition().z() << ") " << G4endl;
        G4cout << "Energy (pre)  : " << thePrePoint->GetTotalEnergy() / MeV << " MeV   "
               << "k = " << thePrePoint->GetKineticEnergy() / MeV << G4endl;
        G4cout << "Energy (post) : " << thePostPoint->GetTotalEnergy() / MeV << " MeV   "
               << "k = " << thePostPoint->GetKineticEnergy() / MeV << G4endl;
        G4cout << "Energy (track): " << theTrack->GetTotalEnergy() / MeV << " MeV   "
               << "k = " << theTrack->GetKineticEnergy() / MeV << G4endl;
        G4cout << "deltaE  : " << (thePrePoint->GetTotalEnergy() - thePostPoint->GetTotalEnergy()) / MeV << G4endl;
        G4cout << "deltaK  : " << (thePrePoint->GetKineticEnergy() - thePostPoint->GetKineticEnergy()) / MeV << G4endl;
        G4cout << "Deposited Energy " << theStep->GetTotalEnergyDeposit() / MeV << G4endl;
        G4cout << "Missing Energy "
               << (thePrePoint->GetTotalEnergy() - thePostPoint->GetTotalEnergy() - theStep->GetTotalEnergyDeposit()) /
                  MeV
               << "  (i.e. non deposited)" << G4endl;
        G4cout << "Non Ionizing En  " << theStep->GetNonIonizingEnergyDeposit() / MeV << G4endl;
        G4cout << "process name : " << pp_proc_name << "  (" << (G4int) pp_proc_type << " / " << pp_proc_subtype << ")"
               << G4endl;
        G4cout << "sensitive detector address : " << thisSD << G4endl;
        G4cout << "step status: " << stepStatus << G4endl;
        G4cout << post_step_proc_list;
        G4cout << "n scint photons: " << n_scin_phot << G4endl;
        G4cout << "n chere photons: " << n_cher_phot << G4endl;
        //  G4Scintillation *proc = (G4Scintillation*) theStep->GetPostStepPoint()->GetProcessDefinedStep();
        //  G4int n_phot =proc->GetNumPhotons();

        if (thisSD) {
            G4cout << "sensitive detector type: " << sd_type << G4endl;
            G4cout << "sensitive detector name: "
                   << dynamic_cast<OFOS_OpticalPhotonSD *>(thisSD)->get_name() << G4endl;
        }


        G4cout << "Number of secondaries: " << theStep->GetNumberOfSecondariesInCurrentStep() << G4endl;
        G4cout << list_of_secondaries;
        //const std::vector<const G4Track*> *tracks = theStep->GetSecondaryInCurrentStep();
        //for(size_t i=0; i<tracks->size(); ++i)
        //{
        //    if(tracks->at(i)->GetDefinition()->GetPDGEncoding() != 0)
        //    {
        //        G4cout << "sec: trk_id:" << tracks->at(i)->GetTrackID();
        //        G4cout << "(" << tracks->at(i)->GetDefinition()->GetPDGEncoding()<< ")  ";
        //        G4cout << " mom_id:" << tracks->at(i)->GetParentID() << G4endl;
        //    }
        //}
        G4cout << "Energy of Secondaries: " << sec_energy / MeV << "  MeV " << G4endl;
        G4cout << G4endl << G4endl;

    }





    /// if output storage is unavaiable, it is pointless to create a hit
    if (!global_ntuples_ptr) return;



    /// track energy loss by primary particle
    if (theTrack->GetTrackID() == 1) {
        const std::vector<const G4Track *> *sec_tracks = theStep->GetSecondaryInCurrentStep();
//        G4double sec_energy(0.0);
        double en_loss_via_deltas = 0;
        double radiative_en_loss = 0;

        for (auto sec_track : *sec_tracks) {
            if (sec_track->GetDefinition()->GetPDGEncoding() == 22 &&
                /// if the secondary particle is a gamma
                !(pp_proc_type == G4ProcessType::fElectromagnetic &&           /// that has not been emitted in a
                  pp_proc_subtype == G4EmProcessSubType::fAnnihilation))       /// annihilation process
                radiative_en_loss += sec_track->GetKineticEnergy();            /// consider the gamma energy as radiative energy loss

            else if (sec_track->GetDefinition()->GetPDGEncoding() == 11) /// if the secondary particle is an electron
                en_loss_via_deltas += sec_track->GetKineticEnergy();      /// consider its energy as delta-ray energy loss
        }
        global_ntuples_ptr->add_en_loss_via_deltas(en_loss_via_deltas);
        global_ntuples_ptr->add_radiative_en_loss(radiative_en_loss);

        /// save steps along the track (proxy for ionization energy loss)
        /// only in the case of charged particles
        if (fabs(particleDefinition->GetPDGCharge()) > 0) {
            global_ntuples_ptr->fill_interaction(theTrack->GetTrackID(),
                                                 theTrack->GetParentID(),
                                                 100,
                                                 particleDefinition->GetPDGEncoding(),
                                                 thePostPoint->GetPosition(),
                                                 thePrePoint->GetKineticEnergy() - thePostPoint->GetKineticEnergy(),
                                                 thePrePoint->GetTotalEnergy(),
                                                 theTrack->GetGlobalTime());

        }
    }




    /// store properties common to all steps, regardless of the process
    if (theStep->GetTotalEnergyDeposit() > 0) {
        global_ntuples_ptr->add_edep(theStep->GetTotalEnergyDeposit());

        if (thePrePVname == "LS" || thePostPVname == "LS")
            global_ntuples_ptr->add_ls_ion(theStep->GetTotalEnergyDeposit());
        // global_ntuples_ptr->add_ls_ion( theStep->GetTotalEnergyDeposit() - theStep->GetNonIonizingEnergyDeposit ());

        if (theTrack->GetTrackID() == 1) {
            global_ntuples_ptr->add_prim_ion(theStep->GetTotalEnergyDeposit() - theStep->GetNonIonizingEnergyDeposit());
//          global_ntuples_ptr->add_prim_ion_incl_deltas(theStep->GetTotalEnergyDeposit() - theStep->GetNonIonizingEnergyDeposit ());
        }

//      if( theTrack->GetParentID() == 1 && abs(particleDefinition->GetPDGEncoding()) == 11 )
//          global_ntuples_ptr->add_prim_ion_incl_deltas(theStep->GetTotalEnergyDeposit() - theStep->GetNonIonizingEnergyDeposit ());
    }

    switch (pp_proc_type) {
        case G4ProcessType::fOptical : {
            /// optical photon absorbed (i.e. detected) by a fiber
            if (pp_proc_subtype == G4OpProcessSubType::fOpAbsorption && thisSD) {
                if (OFOS_Verbosity::level > 3) {
                    G4cout << "Filling Hit :: sensitive_det_type = " << sd_type << " (add: " << thisSD << ")" << G4endl;
                }

                /// >>> Process hit in sensitive detector <<<
                /// However a hit could be generated here (after fetching the SD type,
                /// i.e. thisSD->get_sd_type()), and be sent directly to OutputNtuples
                /// through global_ntuples_ptr->fill_hit(). In this case a loop
                /// over the hit collection could be avoided at each end of event
                dynamic_cast<OFOS_OpticalPhotonSD *>(thisSD)->process_hit(theStep, nullptr);
                theTrack->SetTrackStatus(fStopAndKill);

                if (OFOS_Verbosity::level > 3)
                    G4cout << "Hit Filled" << G4endl;
            }
            break; // case G4ProcessType::fOptical
        }
        case G4ProcessType::fElectromagnetic : {
            if (pp_proc_subtype == G4EmProcessSubType::fBremsstrahlung ||       // 3
                pp_proc_subtype == G4EmProcessSubType::fPhotoElectricEffect ||  // 12
                pp_proc_subtype == G4EmProcessSubType::fComptonScattering ||    // 13
                pp_proc_subtype == G4EmProcessSubType::fAnnihilation ||         // 5
                pp_proc_subtype == G4EmProcessSubType::fIonisation ||           // 2
                pp_proc_subtype == G4EmProcessSubType::fCerenkov)               // 21
            {
                // TODO : why annihilation is never printed out in the case of a positron???
                /// At the moment annihilations are stored thanks to the TrackingAction invoked
                ///     on gammas, but ideally they should be visible also here.
                if (pp_proc_subtype == G4EmProcessSubType::fAnnihilation) {
                    G4cout << "ANNIHILATION!!! (originated by pdgid " << theTrack->GetDefinition()->GetPDGEncoding()
                           << ")"
                           << G4endl;
                }
                /// energy inherited by secondary particles
                G4double dE = (thePrePoint->GetTotalEnergy() - thePostPoint->GetTotalEnergy());
//                dE -= theStep->GetTotalEnergyDeposit();
                G4int interaction_id = G4ProcessType::fElectromagnetic * 1000 + pp_proc_subtype;

                if (OFOS_Verbosity::level > 2)
                    G4cout << "Filling Interaction: trk_id: " << theTrack->GetTrackID()
                           << "  int_id: " << interaction_id
                           << "  dE: " << dE / MeV << G4endl;

                global_ntuples_ptr->fill_interaction(theTrack->GetTrackID(),
                                                     theTrack->GetParentID(),
                                                     interaction_id,
                                                     particleDefinition->GetPDGEncoding(),
                                                     thePostPoint->GetPosition(),
                                                     dE,
                                                     thePrePoint->GetTotalEnergy(),
                                                     theTrack->GetGlobalTime());

            }
            break; // case G4ProcessType::fElectromagnetic
        }

        case G4ProcessType::fTransportation : {
            if (PDGEEncoding != 0) {
                G4double dE = (thePrePoint->GetTotalEnergy() - thePostPoint->GetTotalEnergy());
//                dE -= theStep->GetTotalEnergyDeposit();
                G4int interaction_id = G4ProcessType::fHadronic * 1000 + pp_proc_subtype;

                if (OFOS_Verbosity::level > 2)
                    G4cout << "Filling Interaction: trk_id: " << trackID
                           << "  int_id: " << interaction_id
                           << "  dE: " << dE / MeV << G4endl;

                global_ntuples_ptr->fill_interaction(trackID,
                                                     theTrack->GetParentID(),
                                                     interaction_id,
                                                     PDGEEncoding,
                                                     thePostPoint->GetPosition(),
                                                     static_cast<float>(dE),
                                                     thePrePoint->GetTotalEnergy(),
                                                     static_cast<float>(theTrack->GetGlobalTime()));
            }
            break;
        }
        default:
            break;
    }

    /// formating debug output
    // if (OFOS_Verbosity::level>2 && pp_proc_type != G4ProcessType::fTransportation)
    // if (OFOS_Verbosity::level>2 && pp_proc_type != G4ProcessType::fTransportation && pp_proc_type != G4ProcessType::fOptical)
    //     G4cout << G4endl;


//  G4double x = theTrack->GetVertexPosition().x();
//  G4double y = theTrack->GetVertexPosition().y();
//  G4double z = theTrack->GetVertexPosition().z();


}
