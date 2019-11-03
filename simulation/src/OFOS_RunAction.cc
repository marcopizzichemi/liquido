#include <iostream>
#include <string>
#include <cstdio>
#include <time.h>
#include <fstream>

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "g4root.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmCalculator.hh"
#include "G4Electron.hh"
#include "G4UnitsTable.hh"

#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"

#include "TFile.h"

#include "OFOS_RunAction.h"
#include "OFOS_EventAction.h"
#include "OFOS_Verbosity.h"
#include "OFOS_OutputLog.h"
#include "OFOS_OutputNtuples.h"
#include "OFOS_GlobalNtuplesPtr.h"

#include "OFOS_DetectorConstruction.h"
#include "OFOS_PrimaryGeneratorAction.h"



std::string 
OFOS_RunAction::get_current_time()
{
    time_t raw_time;
    struct tm * time_info;

    time(&raw_time);  /* get current time */
    time_info = localtime (&raw_time);

    int year  = time_info->tm_year - 100;
    int month = time_info->tm_mon + 1;
    int day   = time_info->tm_mday;
    int hour  = time_info->tm_hour;
    int min   = time_info->tm_min;
    int sec   = time_info->tm_sec;


    char date[100];
    sprintf(date, "%d%s%d%s%d%s%d%s%d%s%d", year,  (month<10 ? "0" : "" ) , 
                                            month, (day<10 ? "0": ""), 
                                            day,   (hour<10 ? "0" : ""), 
                                            hour,  (min<10 ? "0" : ""),
                                            min ,  (sec<10 ? "0": ""),
                                            sec);

    std::string out(date);
    return out;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OFOS_RunAction::OFOS_RunAction( OFOS_DetectorConstruction *det ):
    G4UserRunAction(), // ,messenger(0)
    output_hit_file_(0),
    detector_( det )
{ 
    G4cout << "OFOS_RunAction::OFOS_RunAction()" << G4endl;
    // set printing event number per each 100 events
    G4RunManager::GetRunManager()->SetPrintProgress(1);     

    // messenger = new LIMbuSRunActionMessenger(this); 
    
    if(OFOS_Verbosity::level>0)
        G4cout << G4endl << "OFOS_RunAction::OFOS_RunAction::CreateNtuple" << G4endl << G4endl;


}



OFOS_RunAction::~OFOS_RunAction()
{

    delete global_ntuples_ptr; 
    global_ntuples_ptr = 0;

    if(OFOS_Verbosity::level>0)
	G4cout << "OFOS_RunAction deleted" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OFOS_RunAction::BeginOfRunAction(const G4Run* a_run)
{ 

    std::string current_time = get_current_time();

    char run_name[100];
    sprintf (run_name, "_%d", a_run->GetRunID());

    out_filename = "OFOS_";
    out_filename += current_time;
    out_filename += run_name;
    out_filename += "_data.root";


    output_hit_file_ = new TFile(out_filename.data(), "recreate");
    G4cout << "RunAction :: filename = " << out_filename << G4endl;

    G4cout << "Init Ntuples" << G4endl;
    global_ntuples_ptr = new OFOS_OutputNtuples( "op_hits"  , "Optical Photon Hits"          , int(1e6),
                                                 "mc_truth" , "Electromagnetic Interactions" , int(2e3)); 

    
    out_filename = "OFOS_output_";
    out_filename += current_time;
    out_filename += run_name;
    out_filename += ".root";



    log_filename = "OFOS_";
    log_filename += current_time;
    log_filename += run_name;
    log_filename += "_log.txt";

    OFOS_OutputLog::logfile = new std::ofstream(log_filename); 
 // OutputLog::logfile = 0;


    geom_filename = "OFOS_";
    geom_filename += current_time;
    geom_filename += run_name;
    geom_filename += "_geom.txt";

    OFOS_OutputLog::geom_logfile = new std::ofstream(geom_filename); 


    OFOS_OutputLog::log_cache << G4endl;
    OFOS_OutputLog::log_cache << "***************************************" << G4endl;
    OFOS_OutputLog::log_cache << "*        BEGIN OF RUN ACTION          *" << G4endl;
    OFOS_OutputLog::log_cache << "***************************************" << G4endl;


    G4double density     = detector_->get_ls_material()->GetDensity();  
    G4double radl        = detector_->get_ls_material()->GetRadlen();

    OFOS_OutputLog::log_cache << "LS Material: " <<  detector_->get_ls_material()->GetName() << G4endl;
    OFOS_OutputLog::log_cache << "Density:                             " << G4BestUnit(density,"Volumic Mass") << G4endl;
    OFOS_OutputLog::log_cache << "Rad length:                          " << G4BestUnit(radl,   "Length") << G4endl;

    ComputeElectronCriticalEnergy();



    
    G4Element* elementXe = new G4Element("Xenon","Xe",54.,131.29*g/mole);
    G4Material* LXe = new G4Material ("LXeâ",3.02*g/cm3, 1, kStateLiquid, 173.15*kelvin, 1.5*atmosphere);
    LXe -> AddElement(elementXe, 1);



    // PRINT GAMMA CROSS SECTIONS AT 2 MEV
    double gamma_en = 2.0 * MeV;
    G4ParticleDefinition *gamma     = G4Gamma::Gamma();
    G4ProcessManager     *proc_mng  = gamma->GetProcessManager();
    G4ProcessVector      *proc_list = proc_mng->GetProcessList () ;

    G4EmCalculator emCalculator;

    /***********
     *   LAB   *
     ***********/

    double tot_sigma = 0;
    double phot_sigma = 0;

    OFOS_OutputLog::log_cache << G4endl << "Gamma cross section per process (at 2MeV) " << G4endl;
    OFOS_OutputLog::log_cache << G4endl << "***** LAB *****"  << G4endl;
    for(int i=0; i<proc_list->size(); ++i)
    {
        G4VProcess *a_proc = (*proc_list)[i];

        double sigma = emCalculator.ComputeCrossSectionPerVolume( gamma_en, gamma , a_proc->GetProcessName() , detector_->get_ls_material() )/density;
        tot_sigma += sigma;
        OFOS_OutputLog::log_cache << a_proc->GetProcessName(); 

        if( a_proc->GetProcessName() == "phot" ) phot_sigma = sigma;

        /// formatting output
        int s_size = (a_proc->GetProcessName()).size();
        for( int j=0; j<(37-s_size); ++j) OFOS_OutputLog::log_cache << " ";

        OFOS_OutputLog::log_cache << G4BestUnit(sigma, "Surface/Mass") << G4endl; 
    }

    OFOS_OutputLog::log_cache << "Total                                " 
                              << G4BestUnit(tot_sigma, "Surface/Mass") << G4endl;
    OFOS_OutputLog::log_cache << "Photofraction                        " 
                              << phot_sigma / tot_sigma << G4endl;



    /****************
     *     LXe      *
     ***************/
    tot_sigma = 0;
    phot_sigma = 0;

    OFOS_OutputLog::log_cache << G4endl << G4endl << "***** LXe ***** (for comparison)"  << G4endl;
    for(int i=0; i<proc_list->size(); ++i)
    {
        G4VProcess *a_proc = (*proc_list)[i];

        double sigma = emCalculator.ComputeCrossSectionPerVolume( gamma_en, gamma , a_proc->GetProcessName() , LXe )/LXe->GetDensity();;
        tot_sigma += sigma;
        OFOS_OutputLog::log_cache << a_proc->GetProcessName(); 

        if( a_proc->GetProcessName() == "phot" ) phot_sigma = sigma;

        /// formatting output
        int s_size = (a_proc->GetProcessName()).size();
        for( int j=0; j<(37-s_size); ++j) OFOS_OutputLog::log_cache << " ";

        OFOS_OutputLog::log_cache << G4BestUnit(sigma, "Surface/Mass") << G4endl; 
    }
    OFOS_OutputLog::log_cache << "Total                                " 
                              << G4BestUnit(tot_sigma, "Surface/Mass") << G4endl;
    OFOS_OutputLog::log_cache << "Photofraction                        " 
                              << phot_sigma / tot_sigma << G4endl;




    /// dump all the info that was stored in the buffer up to now
    *(OFOS_OutputLog::logfile) << OFOS_OutputLog::log_cache.rdbuf();
    OFOS_OutputLog::logfile->flush();

    *(OFOS_OutputLog::geom_logfile) << OFOS_OutputLog::geom_cache.rdbuf();
    OFOS_OutputLog::geom_logfile->flush();


}



/// from Geant4-10.4.0/examples/extended/electromagnetic/TestEm0
void OFOS_RunAction::ComputeElectronCriticalEnergy()
{

    // compute e- critical energy (Rossi definition) and Moliere radius.
    // Review of Particle Physics - Eur. Phys. J. C3 (1998) page 147
    //
    G4EmCalculator emCal;
      
    const G4Material* material = detector_->get_ls_material();
    const G4double radl = material->GetRadlen();
    G4double ekin = 5*MeV;
    G4double deioni;
    G4double err  = 1., errmax = 0.001;
    G4int    iter = 0 , itermax = 10;  
    while (err > errmax && iter < itermax) 
    {
        iter++;          
        deioni  = radl*
                  emCal.ComputeDEDX(ekin,G4Electron::Electron(),"eIoni",material);
        err = std::abs(deioni - ekin)/ekin;
        ekin = deioni;
    }

    OFOS_OutputLog::log_cache << "Critical Energy (Rossi):            "  << std::setw(8) << G4BestUnit(ekin,"Energy") << G4endl;
           
    //Pdg formula (only for single material)
    G4double pdga[2] = { 610*MeV, 710*MeV };
    G4double pdgb[2] = { 1.24, 0.92 };
    G4double EcPdg = 0.;
    
    if (material->GetNumberOfElements() == 1) {
      G4int istat = 0;
      if (material->GetState() == kStateGas) istat = 1;  
      G4double Zeff = material->GetZ() + pdgb[istat];
      EcPdg = pdga[istat]/Zeff;
      G4cout << "\t\t\t (from Pdg formula : " 
             << std::setw(8) << G4BestUnit(EcPdg,"Energy") << ")";    
    }
       

    const G4double Es = 21.2052*MeV;
    G4double rMolier1 = Es/ekin, rMolier2 = rMolier1*radl;

 // OFOS_OutputLog::log_cache << "Moliere Radius:                      "  << std::setw(8) << rMolier1 << G4endl
 //                           << "X0 = "  << std::setw(8) << G4BestUnit(rMolier2,"Length") << G4endl;
           
    if (material->GetNumberOfElements() == 1) 
    {
       G4double rMPdg = radl*Es/EcPdg;
       G4cout << "\t (from Pdg formula : " 
              << std::setw(8) << G4BestUnit(rMPdg,"Length") << ")";    
     }         


}




void OFOS_RunAction::EndOfRunAction(const G4Run* )
{
 // global_ntuples_ptr->write(); <-- It causes the tree headers to be written twice since
 //                                  tree->write() method is implicitly called by file->write()
    output_hit_file_->Write();
    output_hit_file_->Close();


    OFOS_OutputLog::logfile->close();

    *(OFOS_OutputLog::geom_logfile) << OFOS_OutputLog::geom_cache.rdbuf();
    OFOS_OutputLog::geom_logfile->flush();
    OFOS_OutputLog::geom_logfile->close();

    G4cout << "EndOfRunAction :: all files closed " << G4endl;

 // delete global_ntuples_ptr; <-- NOT ALLOWED
    global_ntuples_ptr = 0;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
