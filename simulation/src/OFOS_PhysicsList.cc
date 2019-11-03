// Physics List
//this
#include "OFOS_PhysicsList.h"
#include "OFOS_Verbosity.h"

//std
#include <iomanip>

//Geant 4
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4Material.hh"
#include "G4ios.hh"
#include "G4UImanager.hh"
#include "globals.hh"

//MBA
#include "G4MuonDecayChannel.hh"
#include "G4SystemOfUnits.hh"
#include "LBE.hh"

/**************************
 *     EM PROCESSES      *
 *************************/


OFOS_PhysicsList::OFOS_PhysicsList() : G4VUserPhysicsList()
                                    /* physicsMessenger(0)*/ 
{
    SetVerboseLevel(0);

    if(OFOS_Verbosity::level>1)
        G4cout << "OFOS_PhysicsList()" << G4endl;

 // RANGE CUTS are defined in the SetCuts() method
 //
 // defaultCutValue = 0.1*cm;
 // electronCutValue = 5.0*mm;
 // gammaCutValue = 1*mm;

//  physicsMessenger = new PhysicsMessenger(this);
}


OFOS_PhysicsList::~OFOS_PhysicsList() 
{
    //delete physicsMessenger;
    //physicsMessenger = 0;
}


void OFOS_PhysicsList::ConstructParticle() 
{
    if(OFOS_Verbosity::level>1) G4cout << "OFOS_PhysicsList::ConstructParticle()" << G4endl;

    /// basically all the families except for the short-lived one
    /// http://geant4.cern.ch/G4UsersDocuments/UsersGuides/ForApplicationDeveloper/html/GettingStarted/particleDef.html
    G4LeptonConstructor leptonConstructor;
                        leptonConstructor.ConstructParticle();
    G4MesonConstructor  mesonConstructor;
                        mesonConstructor.ConstructParticle();
    G4BaryonConstructor baryonConstructor;
                        baryonConstructor.ConstructParticle();
    G4BosonConstructor  bosonConstructor;
                        bosonConstructor.ConstructParticle();
    G4IonConstructor    ionConstructor;
                        ionConstructor.ConstructParticle();


    /// not present in OpNovicePhysicsList
    G4OpticalPhoton::OpticalPhotonDefinition();
    G4Geantino::GeantinoDefinition();

    if(OFOS_Verbosity::level>1) G4cout << "OFOS_PhysicsList::ConstructParticle() :: Done" << G4endl;
}


void OFOS_PhysicsList::ConstructProcess() 
{
    if(OFOS_Verbosity::level>1) G4cout << "OFOS_PhysicsList::ConstructProcess()" << G4endl;


    AddTransportation();
    ConstructDecay();
    ConstructEM();
    ConstructOp();

 // ConstructGeneral();
 // ConstructHad();

 // ConstructHad_QSGP_BIC_HP();//MBA 4/12/15 for neutrons
 // ConstructHadEl();//MBA 4/12/15 for neutrons

    if(OFOS_Verbosity::level>1)  G4cout << "OFOS_PhysicsList::ConstructProcess() :: Done" << G4endl;
}


#include "G4EmStandardPhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4LossTableManager.hh"
#include "G4EmProcessOptions.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

//#include "G4MultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4eMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"

void OFOS_PhysicsList::ConstructEM() 
{
    if(OFOS_Verbosity::level>1) G4cout << "OFOS_PhysicsList::ConstructEM()" << G4endl;
    
    auto theParticleIterator=GetParticleIterator();
    theParticleIterator->reset();

    while( (*theParticleIterator)() )
    {
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pmanager     = particle->GetProcessManager();
        G4String particleName          = particle->GetParticleName();
        G4double charge                = particle->GetPDGCharge();
        
        if (particleName == "gamma") 
        {
            pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
            pmanager->AddDiscreteProcess(new G4ComptonScattering);
            pmanager->AddDiscreteProcess(new G4GammaConversion);
        } 
        else if (particleName == "e-") 
        {
           // G4int G4ProcessManager::AddProcess ( G4VProcess *    aProcess,
           //                                      G4int   ordAtRestDoIt = ordInActive,
           //                                      G4int   ordAlongSteptDoIt = ordInActive,
           //                                      G4int   ordPostStepDoIt = ordInActive    
           //                                     )   

            pmanager->AddProcess(new G4eMultipleScattering,-1, 1, 1);
            pmanager->AddProcess(new G4eIonisation,        -1, 2, 2);
         // DoIt is inactive if ordering parameter is negative (cf G4ProcessManager @ 474)
         // pmanager->AddProcess(new G4eBremsstrahlung(),  -1,-3, 3);  // legacy from Stefano
            pmanager->AddProcess(new G4eBremsstrahlung(),  -1, 3, 3);  // from OpNovice
            
        } 
        else if (particleName == "e+") 
        {
            
            /* G4VProcess* theeplusMultipleScattering = new G4eMultipleScattering();
             G4VProcess* theeplusIonisation         = new G4eIonisation();
             G4VProcess* theeplusBremsstrahlung     = new G4eBremsstrahlung();
             G4VProcess* theeplusAnnihilation       = new G4eplusAnnihilation();
             
             pmanager->AddProcess(theeplusMultipleScattering);
             pmanager->AddProcess(theeplusIonisation);
             pmanager->AddProcess(theeplusBremsstrahlung);
             pmanager->AddProcess(theeplusAnnihilation);
             //
             // set ordering for AtRestDoIt
             pmanager->SetProcessOrderingToFirst(theeplusAnnihilation, idxAtRest);
             //
             // set ordering for AlongStepDoIt
             pmanager->SetProcessOrdering(theeplusMultipleScattering, idxAlongStep,1);
             pmanager->SetProcessOrdering(theeplusIonisation,         idxAlongStep,2);
             //
             // set ordering for PostStepDoIt
             pmanager->SetProcessOrdering(theeplusMultipleScattering, idxPostStep,1);
             pmanager->SetProcessOrdering(theeplusIonisation,         idxPostStep,2);
             pmanager->SetProcessOrdering(theeplusBremsstrahlung,     idxPostStep,3);
             pmanager->SetProcessOrdering(theeplusAnnihilation,       idxPostStep,4);*/
            
            pmanager->AddProcess(new G4eMultipleScattering,-1, 1, 1);
            pmanager->AddProcess(new G4eIonisation,        -1, 2, 2);
         // pmanager->AddProcess(new G4eBremsstrahlung,    -1,-3, 3); // legacy from Stefano
            pmanager->AddProcess(new G4eBremsstrahlung,    -1, 3, 3); // from OpNovice
            pmanager->AddProcess(new G4eplusAnnihilation,   0,-1, 4);
            
        } else if (particleName == "mu+" ||
                   particleName == "mu-"    ) {
            
            pmanager->AddProcess(new G4hMultipleScattering,-1, 1, 1);
            pmanager->AddProcess(new G4MuIonisation,       -1, 2, 2);
         // pmanager->AddProcess(new G4MuBremsstrahlung,   -1,-3, 3); // legacy from Stefano
         // pmanager->AddProcess(new G4MuPairProduction,   -1,-4, 4); // legacy from Stefano
            pmanager->AddProcess(new G4MuBremsstrahlung,   -1, 3, 3); // from OpNovice
            pmanager->AddProcess(new G4MuPairProduction,   -1, 4, 4); // from OpNovice
            
            
        } else if (particleName == "alpha" ||
                   particleName == "He3" ||
                   particleName == "GenericIon") {
            
            pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
            pmanager->AddProcess(new G4ionIonisation,       -1, 2, 2);
            //     pmanager->AddProcess(new G4hLowEnergyIonisation  -1, 2, 2);
            
        } else if (particleName == "pi+" ||
                   particleName == "pi-" ||
                   particleName == "proton" ) {
            
            pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
            pmanager->AddProcess(new G4hIonisation,         -1, 2, 2);
            pmanager->AddProcess(new G4hBremsstrahlung,     -1,-3, 3);
            pmanager->AddProcess(new G4hPairProduction,     -1,-4, 4);
            
        } else if (particleName == "B+" ||
                   particleName == "B-" ||
                   particleName == "D+" ||
                   particleName == "D-" ||
                   particleName == "Ds+" ||
                   particleName == "Ds-" ||
                   particleName == "anti_lambda_c+" ||
                   particleName == "anti_omega-" ||
                   particleName == "anti_proton" ||
                   particleName == "anti_sigma_c+" ||
                   particleName == "anti_sigma_c++" ||
                   particleName == "anti_sigma+" ||
                   particleName == "anti_sigma-" ||
                   particleName == "anti_xi_c+" ||
                   particleName == "anti_xi-" ||
                   particleName == "deuteron" ||
                   particleName == "kaon+" ||
                   particleName == "kaon-" ||
                   particleName == "lambda_c+" ||
                   particleName == "omega-" ||
                   particleName == "sigma_c+" ||
                   particleName == "sigma_c++" ||
                   particleName == "sigma+" ||
                   particleName == "sigma-" ||
                   particleName == "tau+" ||
                   particleName == "tau-" ||
                   particleName == "triton" ||
                   particleName == "xi_c+" ||
                   particleName == "xi-" ||
                   (particleName == "nucleus" && charge != 0)) {
            
            pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
            pmanager->AddProcess(new G4hIonisation,         -1, 2, 2);
            
        } else if  ((!particle->IsShortLived()) &&
                    (particle->GetPDGCharge() != 0.0) &&
                    (particle->GetParticleName() != "chargedgeantino")) {
            
            pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
            pmanager->AddProcess(new G4hIonisation,         -1, 2, 2);
            
        }
    }
    
    
   
    if(OFOS_Verbosity::level>1) G4cout << "OFOS_PhysicsList::ConstructEM() :: Done" << G4endl;
}






//////////////////////////////////////////////////////////////////////////////
// Optics             ////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include "G4Cerenkov.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4OpBoundaryProcess.hh"

//#include "BxOpAbsorptionReemission.hh"
#include "G4OpAbsorption.hh"

//#include "BxScintillation.hh"
#include "G4Scintillation.hh"
#include "G4OpWLS.hh"//MBA 03/10/15

void OFOS_PhysicsList::ConstructOp()
{
    
    if(OFOS_Verbosity::level>1) G4cout << "OFOS_PhysicsList::ConstructOp()" << G4endl;
    
    auto* theScintProcess = new G4Scintillation("Scintillation");
                     theScintProcess->SetVerboseLevel(0);
                     theScintProcess -> SetTrackSecondariesFirst(true);
                     /// meant to have a particle-dependent light yield 
                     theScintProcess -> SetScintillationYieldFactor(1.);


    auto*     theCerenkovProcess = new G4Cerenkov("Cerenkov");
                    theCerenkovProcess->SetTrackSecondariesFirst(true);
                    /// affect the step definition, not the cherenkov yield
                    theCerenkovProcess->SetMaxNumPhotonsPerStep(30);
                    theCerenkovProcess->SetMaxBetaChangePerStep(10.0);


    auto*        theRayleighScatteringProcess = new G4OpRayleigh();
    auto*      theAbsorptionProcess         = new G4OpAbsorption();
    auto* theBoundaryProcess           = new G4OpBoundaryProcess();
    auto*             theReemissionProcess         = new G4OpWLS();
    auto*           theMieHGScatteringProcess    = new G4OpMieHG();
    
    
    //theBoundaryProcess->DumpPhysicsTable();
    //theAbsorptionProcess->DumpPhysicsTable();
    
    G4cout << "PhysicsList :: Scintillation Physics Table" << G4endl;
    theScintProcess->DumpPhysicsTable();
    G4cout << "PhysicsList :: End of Scintillation Physics Table" << G4endl;
    
    
    theAbsorptionProcess        ->SetVerboseLevel(0);
    theRayleighScatteringProcess->SetVerboseLevel(0);
    theBoundaryProcess          ->SetVerboseLevel(0);
    theReemissionProcess        ->SetVerboseLevel(0);
    theMieHGScatteringProcess   ->SetVerboseLevel(0);

    
 // Birks constant is applied directly to LS in detector construction 
 // G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
 //                 theScintProcess->AddSaturation(emSaturation);
    
    
    auto theParticleIterator=GetParticleIterator();
    theParticleIterator->reset();
    while( (*theParticleIterator)() )
    {
        G4ParticleDefinition* particle     = theParticleIterator->value();
        G4ProcessManager*     pmanager     = particle->GetProcessManager();
        G4String              particleName = particle->GetParticleName();
        
        if (theCerenkovProcess->IsApplicable(*particle)) 
        {
            pmanager->AddProcess(theCerenkovProcess);
            pmanager->SetProcessOrdering(theCerenkovProcess,idxPostStep);
        }
        if (theScintProcess->IsApplicable(*particle)) 
        {
            pmanager->AddProcess(theScintProcess);
            pmanager->SetProcessOrderingToLast(theScintProcess, idxAtRest);
            pmanager->SetProcessOrderingToLast(theScintProcess, idxPostStep);
        }
        /*  if (particleName == "mu+" || particleName == "mu-")   {
         pmanager->AddProcess(theMuScintillationProcess);
         pmanager->SetProcessOrdering(theMuScintillationProcess, idxAtRest);
         pmanager->SetProcessOrdering(theMuScintillationProcess, idxPostStep);
         }
         if (particleName == "proton")   {
         pmanager->AddProcess(theProtonScintillationProcess);
         pmanager->SetProcessOrdering(theProtonScintillationProcess, idxAtRest);
         pmanager->SetProcessOrdering(theProtonScintillationProcess, idxPostStep);
         }
         if (particleName == "alpha")   {
         pmanager->AddProcess(theAlphaScintillationProcess);
         pmanager->SetProcessOrdering(theAlphaScintillationProcess, idxAtRest);
         pmanager->SetProcessOrdering(theAlphaScintillationProcess, idxPostStep);
         }*/
        //}
        
        
        if (particleName == "opticalphoton") 
        {
            G4cout << " AddDiscreteProcess to OpticalPhoton " << G4endl;
            pmanager->AddDiscreteProcess(theAbsorptionProcess);
            pmanager->AddDiscreteProcess(theRayleighScatteringProcess);
            pmanager->AddDiscreteProcess(theMieHGScatteringProcess);
            pmanager->AddDiscreteProcess(theReemissionProcess);//MBA 3/10/15
            pmanager->AddDiscreteProcess(theBoundaryProcess);
        }
    }
    G4cout << "OFOS_PhysicsList::ConstructOp() :: Done" << G4endl;
}


//MBA 4/12/15   elastic processes for hadrons
#include "G4HadronElasticPhysics.hh"
void OFOS_PhysicsList::ConstructHadEl() {
    G4cout << "OFOS_PhysicsList:: Hadronic Elastic Scattering Active" << G4endl;
    auto* hElPhysicsList = new G4HadronElasticPhysics(1);//verbose
    hElPhysicsList->SetVerboseLevel(0);
    hElPhysicsList->ConstructProcess();
    G4cout << "OFOS_PhysicsList:: Hadronic Elastic Scattering :: Done" << G4endl;
}



//QGSP BINARY cascade NEUTRONHP model
//MBA nel 2013 era la meglio per i neutroni
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
void OFOS_PhysicsList::ConstructHad_QSGP_BIC_HP() {
    G4cout << "Hadronic Physics Active: QSGP_BIC_HP model" << G4endl;
    auto *hadPhysicsList = new G4HadronPhysicsQGSP_BIC_HP;
    hadPhysicsList->SetVerboseLevel(0);
    hadPhysicsList->ConstructProcess();
    G4cout << "Hadronic Physics Active: QSGP_BIC_HP model :: Done" << G4endl;
    
}


#include "G4Decay.hh"

void OFOS_PhysicsList::ConstructDecay() {
    G4cout << "PhysicsList: decays "<<G4endl;
    // Add Decay Process
    auto* theDecayProcess = new G4Decay();

    auto theParticleIterator=GetParticleIterator();
    theParticleIterator->reset();

    while( (*theParticleIterator)() )
    {
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4String              particleName      = particle->GetParticleName();
 //     if(particleName != "mu+")
 //     {  
            G4ProcessManager* pmanager = particle->GetProcessManager();
            if (theDecayProcess->IsApplicable(*particle)) 
            { 
                pmanager ->AddProcess(theDecayProcess);
                // set ordering for PostStepDoIt and AtRestDoIt
                pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
                pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
            }
 //     }
    }
}//ConstructGeneral



void OFOS_PhysicsList::SetCuts()
{
    
    ///special for low energy physics
    G4double lowlimit=250*eV;  
    G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowlimit,100.*GeV);
    
    /// set cut values for gamma at first and for e- second and next for e+,
    /// because some processes for e+/e- need cut values for gamma
     
    double cut_em = defaultCutValue;
    G4double cut_had= defaultCutValue;
    G4double cut_ion= 1000*mm;
    
   
    /// NO REASON AT THE MOMENT TO OVERRIDE DEFULT CUT VALUES
 // SetCutValue(gammaCutValue, "gamma");
 // SetCutValue(electronCutValue, "e-");
 // SetCutValue(electronCutValue, "e+");
 // 
 // SetCutValue(cut_had, "proton");
 // SetCutValue(cut_had, "anti_proton");
 // SetCutValue(cut_had, "neutron");
 // 
 // SetCutValue(cut_ion, "alpha");
 // SetCutValue(cut_ion, "GenericIon");

    SetCutsWithDefault(); 
    
    
}
