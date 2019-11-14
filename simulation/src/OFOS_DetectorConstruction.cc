#include <algorithm>

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4RunManager.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4SDManager.hh"
#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "OFOS_DetectorConstruction.h"
#include "OFOS_DetectorMessenger.h"
#include "OFOS_OpticalPhotonSD.h"
#include "OFOS_Verbosity.h"
#include "OFOS_LsMatProperties.h"
#include "OFOS_OutputLog.h"


/******************
 *   
 *   Readout Unit Geometry
 *
 *   Horizontal Readout Units [new]:
 *   The Readout takes place on the XZ plane
 *   They are needed to reconstruct the Z coordinate
 *
 *   Vertical Readout Units:
 *   The readout takes place on the XY plane
 *   They are used for the main topology on the horizontal plane
 *
 ********************/

G4ThreadLocal G4GlobalMagFieldMessenger *OFOS_DetectorConstruction::fMagFieldMessenger = nullptr;

OFOS_DetectorConstruction::OFOS_DetectorConstruction() : G4VUserDetectorConstruction(),
                                                         logic_vert_ru_sipmt_(nullptr),
                                                         logic_hori_ru_sipmt_(nullptr),
                                                         logic_vert_fiber_(nullptr),
                                                         logic_hori_fiber_(nullptr),
                                                         logic_vessel_(nullptr),
                                                         fCheckOverlaps(false),
                                                         logic_ls_(nullptr),
                                                         ls_mpt(nullptr),
                                                         ls(nullptr),
                                                         air(nullptr),
                                                         is_geom_built(false),
                                                         vert_sipmt_SD_(nullptr),
                                                         vert_fiber_SD_(nullptr),
                                                         hori_sipmt_SD_(nullptr),
                                                         hori_fiber_SD_(nullptr),
                                                         vessel_SD_(nullptr)
//  optical_photon_sd_(0)
{
    OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::OFOS_DetectorConstruction()" << G4endl;

    /// detector parameter initialization
    vert_ru_geometry_ = ReadoutUnitGeometry::Square;
    hori_ru_geometry_ = ReadoutUnitGeometry::Square;
    number_ru_x_ = 10;
    number_ru_y_ = 15;
    number_ru_z_ = 12;
    vert_ru_distance_ = 1.0 * cm;
    vert_ru_size_ = 0.2 * cm;
    hori_ru_distance_ = 10.0 * cm;
    hori_ru_size_ = 0.2 * cm;
    fiber_radius_ = 0.25 * mm;
    outer_cladding_fractional_radius_ = 0.03;
    inner_cladding_fractional_radius_ = 0.03;
    beta_emitter_wall_thickness_ = 1e-8 * cm;
    air_buffer_thickness_ = 20 * cm;
    outer_vessel_thickness_ = 1. * cm;
    ls_vessel_thickness_ = 1. * cm;


    messenger_ = new OFOS_DetectorMessenger(this);


    //define_materials();


    /// LS properties are instantiated together with the detector 
    /// so that they can be modified at any time
    my_ls_properties = new OFOS_LsMatProperties();

    /// LS material property table is instantiated on demand in
    /// deploy_ls_properties(), so that several tables can be used
    /// at runtime
    //ls_mpt = new G4MaterialPropertiesTable();
}



/***********************************************************************************/
/***********************************************************************************/


OFOS_DetectorConstruction::~OFOS_DetectorConstruction() {
    delete my_ls_properties;
    delete ls_mpt;
    delete messenger_;
}


/***********************************************************************************/
/***********************************************************************************/


G4VPhysicalVolume *OFOS_DetectorConstruction::Construct() {
    OFOS_OutputLog::log_cache << "***************************************" << G4endl;
    OFOS_OutputLog::log_cache << "*   DEFAULT GEOMETRY IMPLEMENTATION   *" << G4endl;
    OFOS_OutputLog::log_cache << "***************************************" << G4endl;

    OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::Construct()" << G4endl;

    /// define all materials but LS
    define_materials();

    /// define LS properties
    deploy_ls_properties();

    set_field();

    /// build geometry and return pointer to World volume
    return build_geom();
}



/***********************************************************************************/
/***********************************************************************************/


void OFOS_DetectorConstruction::define_materials() {
    if (OFOS_Verbosity::level > 1)
        OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::DefineMaterials()" << G4endl;


    auto *C = new G4Element("Carbon", "C", 6., 12.01 * g / mole);
    auto *H = new G4Element("Hydrogen", "H", 1., 1.01 * g / mole);
    auto *O = new G4Element("Oxygen", "O", 8., 16.00 * g / mole);
    auto *N = new G4Element("Nitrogen", "N", 7., 14.01 * g / mole);
    auto *Si = new G4Element("Silicon", "Si", 14., 28.09 * g / mole);
    auto *Al = new G4Element("Aluminium", "Al", 13., 26.98 * g / mole);
    auto *S = new G4Element("Sulfur", "S", 16., 32.066 * g / mole);
//  G4Element *Pb = new G4Element("Lead", "Pb", 82., 207.2*g/mole);



    /* DUMMY FIBER OPTICAL PROPERTIES */
    G4double energy_10[10] =
            {1.0 * eV,
             2.0 * eV,
             3.0 * eV,
             4.0 * eV,
             5.0 * eV,
             6.0 * eV,
             7.0 * eV,
             8.0 * eV,
             9.0 * eV,
             10.0 * eV};

    G4double out_clad_n[10] =
            {1.42,
             1.42,
             1.42,
             1.42,
             1.42,
             1.42,
             1.42,
             1.42,
             1.42,
             1.42};

    G4double inn_clad_n[10] =
            {1.49,
             1.49,
             1.49,
             1.49,
             1.49,
             1.49,
             1.49,
             1.49,
             1.49,
             1.49};

    G4double core_n[10] =
            {1.59,
             1.59,
             1.59,
             1.59,
             1.59,
             1.59,
             1.59,
             1.59,
             1.59,
             1.59};

    G4double fib_att_len[10] =
            {16 * m,
             16 * m,
             16 * m,
             16 * m,
             16 * m,
             16 * m,
             16 * m,
             16 * m,
             16 * m,
             16 * m};

    G4double core_att_len[10] =
            {7.0E-1 * mm,
             7.0E-1 * mm,
             7.0E-1 * mm,
             7.0E-1 * mm,
             7.0E-1 * mm,
             7.0E-1 * mm,
             7.0E-1 * mm,
             7.0E-1 * mm,
             7.0E-1 * mm,
             7.0E-1 * mm};


    G4double black_10[10] =
            {1.0E-9 * mm,
             1.0E-9 * mm,
             1.0E-9 * mm,
             1.0E-9 * mm,
             1.0E-9 * mm,
             1.0E-9 * mm,
             1.0E-9 * mm,
             1.0E-9 * mm,
             1.0E-9 * mm,
             1.0E-9 * mm};
    /* end of DUMMYFIBER OPTICAL PROPERTIES */



//  /*** DEFINE LAB-BASED LS ***/
    lab = new G4Material("LAB", 0.859 * g / cm3, 5);
    lab->AddElement(C, 0.87924);
    lab->AddElement(H, 0.1201);
    lab->AddElement(O, 0.00034);
    lab->AddElement(N, 0.00027);
    lab->AddElement(S, 0.00005);


    /// Start with undoped Liquid Scintillator
    ls = lab;


    /// LS properties are initialized at the moment of geoemtry building
    /// ( method Construct() and update_geom() )
    // deploy_ls_properties();



    G4NistManager *man = G4NistManager::Instance();
    air = man->FindOrBuildMaterial("G4_AIR");


    black_acrylic = new G4Material("BlackAcrylic", 1.18 * g / cm3, 3);
    black_acrylic->AddElement(C, 0.59984);
    black_acrylic->AddElement(H, 0.08055);
    black_acrylic->AddElement(O, 0.31961);

    auto *black_mpt = new G4MaterialPropertiesTable();


    polystyrene = new G4Material("Polystyrene", 1.05 * g / cm3, 2);
    polystyrene->AddElement(C, 0.5);
    polystyrene->AddElement(H, 0.5);

    auto *polystyrene_mpt = new G4MaterialPropertiesTable();
    polystyrene_mpt->AddProperty("RINDEX", energy_10, core_n, 10);
    polystyrene_mpt->AddProperty("ABSLENGTH", energy_10, core_att_len, 10);
    polystyrene->SetMaterialPropertiesTable(polystyrene_mpt);


    polymethylmethacrylate = new G4Material("PMMA", 1.19 * g / cm3, 3);
    polymethylmethacrylate->AddElement(C, 0.336449);
    polymethylmethacrylate->AddElement(H, 0.532710);
    polymethylmethacrylate->AddElement(O, 0.130841);

    auto *polymethylmethacrylate_mpt = new G4MaterialPropertiesTable();
    polymethylmethacrylate_mpt->AddProperty("RINDEX", energy_10, inn_clad_n, 10);
    polymethylmethacrylate_mpt->AddProperty("ABSLENGTH", energy_10, fib_att_len, 10);
    polymethylmethacrylate->SetMaterialPropertiesTable(polymethylmethacrylate_mpt);

    out_cladding_mat = new G4Material("OuterCladding", 1.43 * g / cm3, 3);
    out_cladding_mat->AddElement(C, 0.336449);  /// WRONG COMPOSITION
    out_cladding_mat->AddElement(H, 0.532710);  /// JUST A PLACEHOLDER
    out_cladding_mat->AddElement(O, 0.130841);

    auto *out_cladding_mpt = new G4MaterialPropertiesTable();
    out_cladding_mpt->AddProperty("RINDEX", energy_10, out_clad_n, 10);
    out_cladding_mpt->AddProperty("ABSLENGTH", energy_10, fib_att_len, 10);
    out_cladding_mat->SetMaterialPropertiesTable(out_cladding_mpt);


    black_mpt->AddProperty("RINDEX", energy_10, core_n, 10);  /// PLACEHOLDER
    black_mpt->AddProperty("ABSLENGTH", energy_10, black_10, 10);
    black_acrylic->SetMaterialPropertiesTable(black_mpt);


    // Print materials
    if (OFOS_Verbosity::level > 2)
        OFOS_OutputLog::log_cache << G4Material::GetMaterialTable() << G4endl;
}



/***********************************************************************************/
/***********************************************************************************/








void
OFOS_DetectorConstruction::ConstructSDandField() {
    // Sensitive detectors

    G4String vert_sipmt_SDname = "OFOS/VertSiPmtSD";
    G4String vert_fiber_SDname = "OFOS/VertFiberSD";
    G4String hori_sipmt_SDname = "OFOS/HoriSiPmtSD";
    G4String hori_fiber_SDname = "OFOS/HoriFiberSD";
    G4String vessel_SDname = "OFOS/VesselSD";

    OFOS_OutputLog::log_cache << "---------------------------------------" << G4endl;
    OFOS_OutputLog::log_cache << "***       sensitive detectors       ***" << G4endl;

    vert_sipmt_SD_ = new OFOS_OpticalPhotonSD(vert_sipmt_SDname, "VertSiPmtHitsCollection", SDType::Pmt);
    vert_fiber_SD_ = new OFOS_OpticalPhotonSD(vert_fiber_SDname, "VertFiberHitsCollection", SDType::VerticalFiber);
    hori_sipmt_SD_ = new OFOS_OpticalPhotonSD(hori_sipmt_SDname, "HoriSiPmtHitsCollection", SDType::Pmt);
    hori_fiber_SD_ = new OFOS_OpticalPhotonSD(hori_fiber_SDname, "HoriFiberHitsCollection", SDType::HorizontalFiber);
    vessel_SD_ = new OFOS_OpticalPhotonSD(vessel_SDname, "VesselHitsCollection", SDType::Vessel);

    /// new since Geant4 v10.3
    G4SDManager::GetSDMpointer()->AddNewDetector(vert_sipmt_SD_);
    G4SDManager::GetSDMpointer()->AddNewDetector(vert_fiber_SD_);
    G4SDManager::GetSDMpointer()->AddNewDetector(hori_sipmt_SD_);
    G4SDManager::GetSDMpointer()->AddNewDetector(hori_fiber_SD_);
    G4SDManager::GetSDMpointer()->AddNewDetector(vessel_SD_);


    if (is_geom_built) {
        if (OFOS_Verbosity::level > 1)
            OFOS_OutputLog::log_cache
                    << "OFOS_DetectorConstruction::ConstructSDandField() :: Setting Sensitive Detectors" << G4endl;

        SetSensitiveDetector(logic_vert_ru_sipmt_, vert_sipmt_SD_);
        SetSensitiveDetector(logic_hori_ru_sipmt_, hori_sipmt_SD_);
        SetSensitiveDetector(logic_vert_fiber_, vert_fiber_SD_);
        SetSensitiveDetector(logic_hori_fiber_, hori_fiber_SD_);
        SetSensitiveDetector(logic_vessel_, vessel_SD_);
    }
}



/***********************************************************************************/
/***********************************************************************************/



G4Material *
OFOS_DetectorConstruction::make_cocktail(double density, double loading_fraction, G4String &material) {
    /// this method is meant to load the liquid scintillator with a non-native isotope
    if (loading_fraction < 1e-6) return lab;


    /// define here the list of all material suitable to load LSi at runtime
    G4NistManager *man = G4NistManager::Instance();


    G4Material *Pb = man->FindOrBuildMaterial("G4_Pb");
    G4Material *Te = man->FindOrBuildMaterial("G4_Te");
    G4Material *In = man->FindOrBuildMaterial("G4_In");

    G4double lab_density = lab->GetDensity() / (g / cm3);
    G4double pb_density = 11.34; //  g/cm3 --> units are applied when defining the material
    G4double te_density = 6.24; //  g/cm3 --> units are applied when defining the material
    G4double in_density = 7.31; //  g/cm3 --> units are applied when defining the material

    double tot_mass = 1.0 + loading_fraction; // adding loading_fraction[g] to 1g of LAB

    /// start the cocktail with simple LAB
    G4Material *a_isotope = lab;

    /// list here all the available options in terms of loading isotope
    if (material.contains("Pb") && loading_fraction > 1e-6) {
        a_isotope = Pb;

        if (density < 1e-6) {
            /// volume of 1g of LAB + loading_fraction[g] of loading_material
            double tot_volume = 1.0 / lab_density + loading_fraction / pb_density;
            density = tot_mass / tot_volume;
        }
    } else if (material.contains("Te") && loading_fraction > 1e-6) {
        a_isotope = Te;

        if (density < 1e-6) {
            /// volume of 1g of LAB + loading_fraction[g] of loading_material
            double tot_volume = 1.0 / lab_density + loading_fraction / te_density;
            density = tot_mass / tot_volume;
        }
    } else if (material.contains("In") && loading_fraction > 1e-6) {
        a_isotope = In;

        if (density < 1e-6) {
            /// volume of 1g of LAB + loading_fraction[g] of loading_material
            double tot_volume = 1.0 / lab_density + loading_fraction / in_density;
            density = tot_mass / tot_volume;
        }
    } else {
        OFOS_OutputLog::log_cache << G4endl << G4endl;
        OFOS_OutputLog::log_cache << "******************************************" << G4endl;
        OFOS_OutputLog::log_cache << "*                 WARNING                *" << G4endl;
        OFOS_OutputLog::log_cache << "******************************************" << G4endl;

        OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::make_cocktail()" << G4endl;
        OFOS_OutputLog::log_cache << "Loading material (" << material << ") not found" << G4endl;
        OFOS_OutputLog::log_cache << "Using undoped LAB" << G4endl;
        OFOS_OutputLog::log_cache << "******************************************" << G4endl << G4endl << G4endl;
    }

    if (OFOS_Verbosity::level > 1) {
        OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::make_cocktail()" << G4endl;
        OFOS_OutputLog::log_cache << "                           density = " << density << G4endl;
        OFOS_OutputLog::log_cache << "                           loading fraction = " << loading_fraction << G4endl;
        OFOS_OutputLog::log_cache << "                           loading material = " << a_isotope->GetName() << G4endl;
    }


    /// define a new name for the doped LS to avoid conflict with previously defined materials
    G4String ls_name("LS");
    ls_name += "_doped_w_";
    ls_name += material;
    ls_name += "_at_";
    ls_name += G4UIcommand::ConvertToString(loading_fraction);

    auto *a_ls = new G4Material(ls_name, density * g / cm3, 2);
    a_ls->AddMaterial(lab, 1. / (1. + loading_fraction));
    a_ls->AddMaterial(a_isotope, loading_fraction / (1. + loading_fraction));

    return a_ls;
}



/***********************************************************************************/
/***********************************************************************************/


void
OFOS_DetectorConstruction::deploy_ls_properties() {
    if (OFOS_Verbosity::level > 1)
        OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::deploy_ls_properties()" << G4endl;


    /****** DEFINE NEW LS ******/
    /// if(ls) delete ls; /// ownership of materials is not mine !!!
    ls = make_cocktail(my_ls_properties->get_density(),
                       my_ls_properties->get_loading_fraction(),
                       my_ls_properties->get_loading_material());

    OFOS_OutputLog::log_cache << G4endl;
    OFOS_OutputLog::log_cache << "LS PROPERTIES" << G4endl;
    OFOS_OutputLog::log_cache << "Cocktail Name:                       " << ls->GetName() << G4endl;
    OFOS_OutputLog::log_cache << "LS Loading Fraction:                 " << my_ls_properties->get_loading_fraction();
    OFOS_OutputLog::log_cache << (ls == lab ? " (below threshold)" : "") << G4endl;

    OFOS_OutputLog::log_cache << "LS Density:                          " << ls->GetDensity() / (g / cm3);
    OFOS_OutputLog::log_cache << (ls == lab ? " (LAB default)" : "") << G4endl;
    OFOS_OutputLog::log_cache << "                      ext input was: " << my_ls_properties->get_density();
    OFOS_OutputLog::log_cache << (ls == lab || my_ls_properties->get_density() < 1e-6 ? " (ignored)" : "") << G4endl;

    OFOS_OutputLog::log_cache << "LS Loading Material:                 " << my_ls_properties->get_loading_material();
    OFOS_OutputLog::log_cache << (ls == lab ? " (ignored)" : "") << G4endl;

    OFOS_OutputLog::log_cache << "LS Light Yield:                      " << my_ls_properties->get_light_yield()
                              << " / MeV" << G4endl;
    OFOS_OutputLog::log_cache << OFOS_OutputLog::ls_cache.str();
    OFOS_OutputLog::log_cache << "---------------------------------------" << G4endl;

    /// resetting LS log cache
    OFOS_OutputLog::ls_cache.str("");
    OFOS_OutputLog::ls_cache.clear();


    if (OFOS_Verbosity::level > 1) {
        OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::deploy_ls_properties() :: cocktail prepared" << G4endl;
        OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::deploy_ls_properties() :: cocktail name: "
                                  << ls->GetName() << G4endl;
    }

    // if(ls_mpt) delete ls_mpt; /// if deleted causes seg fault -> G4 ownership?

    auto *a_ls_mpt = new G4MaterialPropertiesTable();

    a_ls_mpt->AddProperty("RINDEX", my_ls_properties->get_energy(), my_ls_properties->get_ref_index(),
                          my_ls_properties->get_n_data());
    a_ls_mpt->AddProperty("ABSLENGTH", my_ls_properties->get_energy(), my_ls_properties->get_absorption(),
                          my_ls_properties->get_n_data());
    a_ls_mpt->AddProperty("FASTCOMPONENT", my_ls_properties->get_energy(), my_ls_properties->get_fast_spectrum(),
                          my_ls_properties->get_n_data());
    a_ls_mpt->AddProperty("SLOWCOMPONENT", my_ls_properties->get_energy(), my_ls_properties->get_slow_spectrum(),
                          my_ls_properties->get_n_data());
    a_ls_mpt->AddProperty("RAYLEIGH", my_ls_properties->get_energy(), my_ls_properties->get_scattering(),
                          my_ls_properties->get_n_data());
    a_ls_mpt->AddConstProperty("SCINTILLATIONYIELD", my_ls_properties->get_light_yield());
    a_ls_mpt->AddConstProperty("RESOLUTIONSCALE", my_ls_properties->get_resolution_scale());
    a_ls_mpt->AddConstProperty("FASTTIMECONSTANT", my_ls_properties->get_fast_time_const());
    a_ls_mpt->AddConstProperty("SLOWTIMECONSTANT", my_ls_properties->get_slow_time_const());
    a_ls_mpt->AddConstProperty("YIELDRATIO", my_ls_properties->get_yield_ratio());

    ls->SetMaterialPropertiesTable(a_ls_mpt);

    ls->GetIonisation()->SetBirksConstant(my_ls_properties->get_birks_const() * CLHEP::mm / CLHEP::MeV);

    if (OFOS_Verbosity::level > 1) {
        OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::deploy_ls_properties() :: LS properties updated"
                                  << G4endl;
    }


    /// composition of logical volumes relying on new LS definition
    /// is updated at the moment of building new geometry
}





/***********************************************************************************/
/***********************************************************************************/


void
OFOS_DetectorConstruction::update_geom() {
    OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::update_geom()" << G4endl << G4endl;
    OFOS_OutputLog::log_cache << "***************************************" << G4endl;
    OFOS_OutputLog::log_cache << "*   UPDATED GEOMETRY IMPLEMENTATION   *" << G4endl;
    OFOS_OutputLog::log_cache << "***************************************" << G4endl;

    if (OFOS_Verbosity::level > 1)
        OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::update_geom() :: opening geometry" << G4endl;

    /// reset geometry
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

    /// NB geometry needs to be opened not just to modify positioning of detector's components
    /// bul also in order to modify materials
    /// http://geant4.cern.ch/G4UsersDocuments/UsersGuides/ForApplicationDeveloper/html/Detector/geomDynamic.html


    /// reset geometry log file to dump only final configuration
    OFOS_OutputLog::geom_cache.str("");
    OFOS_OutputLog::geom_cache.clear();


    if (OFOS_Verbosity::level > 1)
        OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::update_geom() :: updating LS properties" << G4endl;

    /// update LS properites
    deploy_ls_properties();

    if (OFOS_Verbosity::level > 1)
        OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::update_geom() :: updating geometry" << G4endl;

    /// update detector geometry
    G4RunManager::GetRunManager()->DefineWorldVolume(build_geom());

    if (OFOS_Verbosity::level > 1)
        OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::update_geom() :: updating sensitive volumes" << G4endl;

    /// update sensitive volumes
    /// http://hypernews.slac.stanford.edu/HyperNews/geant4/get/geometry/1090/2.html
    logic_vert_fiber_->SetSensitiveDetector(vert_fiber_SD_);
    logic_hori_fiber_->SetSensitiveDetector(hori_fiber_SD_);
    logic_vessel_->SetSensitiveDetector(vessel_SD_);
    logic_vert_ru_sipmt_->SetSensitiveDetector(vert_sipmt_SD_);
    logic_hori_ru_sipmt_->SetSensitiveDetector(hori_sipmt_SD_);



    /// issue flag to re-optimize geometry and physics before to run
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();

    if (OFOS_Verbosity::level > 1)
        OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::update_geom() :: geometry now closed" << G4endl;

    /// debug: check the material associated to LS logical volume
    G4LogicalVolumeStore *pLVStore = G4LogicalVolumeStore::GetInstance();
    G4int nLV = pLVStore->size();
    G4int iLV = 0;
    G4LogicalVolume *pLV = nullptr;

    for (iLV = 0; iLV < nLV && OFOS_Verbosity::level > 1; ++iLV) {
        pLV = (*pLVStore)[iLV];
        if (pLV->GetName() == "LS") {
            OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::update_geom() :: LS Material : "
                                      << pLV->GetMaterial()->GetName() << G4endl;
            // pLV->GetMaterial()->GetMaterialPropertiesTable()->DumpTable();
        }
    }
}




/***********************************************************************************/
/***********************************************************************************/




G4VPhysicalVolume *
OFOS_DetectorConstruction::build_geom() {

    OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::build_geom()" << G4endl;

    check_geom_params();

    print_geom_params();


    /// TODO check if this radius is actually containing the RU when Traingle geometry is used
    G4double vert_ru_radius = sqrt(2.) * 0.5 * vert_ru_size_ + fiber_radius_; /// needed for G4 implementation of the RU
    G4double hori_ru_radius = sqrt(2.) * 0.5 * hori_ru_size_ + fiber_radius_;

    /// pmt geometry
    G4double sipmt_radius = fiber_radius_;
    G4double sipmt_sens_radius = fiber_radius_ * (1.0 - outer_cladding_fractional_radius_);
    G4double sipm_thickness = 1. * mm;


    /// lattice geometry
    G4double instrumented_x_size = vert_ru_distance_ * (number_ru_x_ - 1) + 2.0 * vert_ru_radius;
    G4double instrumented_y_size = vert_ru_distance_ * (number_ru_y_ - 1) + 2.0 * vert_ru_radius;
    G4double instrumented_z_size = hori_ru_distance_ * (number_ru_z_ - 1) + 2.0 * hori_ru_radius;


    /// ls volume - adding some extra space between the last fiber and the LS vessel
    G4double ls_x_size = instrumented_x_size + 10 * cm;
    G4double ls_y_size = instrumented_y_size + 10 * cm;
    G4double ls_z_size = instrumented_z_size + 10 * cm;


    /// ls vessel
    G4double ls_vessel_x_size = ls_x_size + 2. * ls_vessel_thickness_;
    G4double ls_vessel_y_size = ls_y_size + 2. * ls_vessel_thickness_;
    G4double ls_vessel_z_size = ls_z_size + 2. * ls_vessel_thickness_;


    /// vertical fibers (running along z)
    G4double out_of_vessel_fiber_len = 0. * cm;
    G4double vert_fiber_length = ls_z_size + 2 * out_of_vessel_fiber_len - 2. * sipm_thickness;
    G4double vert_fiber_env_length = ls_z_size + 2 * out_of_vessel_fiber_len;

    /// horizontal fibers (running along y)
    G4double hori_fiber_length = ls_y_size + 2 * out_of_vessel_fiber_len - 2. * sipm_thickness;
    G4double hori_fiber_env_length = ls_y_size + 2 * out_of_vessel_fiber_len;


    /// outer vessel
    G4double outer_vessel_x_size = ls_vessel_x_size + 2. * outer_vessel_thickness_ + 2. * air_buffer_thickness_;
    G4double outer_vessel_y_size = ls_vessel_y_size + 2. * outer_vessel_thickness_ + 2. * air_buffer_thickness_;
    G4double outer_vessel_z_size = ls_vessel_z_size + 2. * outer_vessel_thickness_ + 2. * air_buffer_thickness_;


    /// external world
    G4double world_x_size = outer_vessel_x_size + 1.0 * m; //size_x_LiquidO_MC_New;
    G4double world_y_size = outer_vessel_y_size + 1.0 * m; //size_y_LiquidO_MC_New;
    G4double world_z_size = outer_vessel_z_size + 1.0 * m; //size_z_LiquidO_MC_New;


    /// info dump
    OFOS_OutputLog::log_cache << "---------------------------------------" << G4endl;
    OFOS_OutputLog::log_cache << "***   geometry derived parameters   ***" << G4endl;
    OFOS_OutputLog::log_cache << "LS size (x,y,z): (" << ls_vessel_x_size / m << ", " << ls_vessel_y_size / m << ", "
                              << ls_vessel_z_size / m << ") m" << G4endl;


    auto *world_s = new G4Box("World", 0.5 * world_x_size, 0.5 * world_y_size, 0.5 * world_z_size);
    auto *world_l = new G4LogicalVolume(world_s, air, "World");
    G4VPhysicalVolume *world_p = new G4PVPlacement(nullptr,                     //no rotation
                                                   G4ThreeVector(),       //at (0,0,0)
                                                   world_l,          //its logical volume
                                                   "World",               //its name
                                                   nullptr,                     //its mother  volume
                                                   false,                 //no boolean operation
                                                   0,                     //copy number
                                                   fCheckOverlaps);       //overlaps checking


    // Outer Vessel
    auto *outer_vessel_s = new G4Box("OuterVessel", 0.5 * outer_vessel_x_size, 0.5 * outer_vessel_y_size,
                                     0.5 * outer_vessel_z_size);
    auto *outer_vessel_l = new G4LogicalVolume(outer_vessel_s, black_acrylic, "OuterVessel");
    G4VPhysicalVolume *outer_vessel_p = new G4PVPlacement(nullptr, G4ThreeVector(), outer_vessel_l, "OuterVessel",
                                                          world_l, false, 0, fCheckOverlaps);


    /// Buffer
    auto *buffer_s = new G4Box("Buffer", 0.5 * (outer_vessel_x_size - 2. * outer_vessel_thickness_),
                               0.5 * (outer_vessel_y_size - 2. * outer_vessel_thickness_),
                               0.5 * (outer_vessel_z_size - 2. * outer_vessel_thickness_));
    auto *buffer_l = new G4LogicalVolume(buffer_s, air, "Buffer");
    G4VPhysicalVolume *buffer_p = new G4PVPlacement(nullptr, G4ThreeVector(), buffer_l, "Buffer", outer_vessel_l, false,
                                                    0, fCheckOverlaps);


    // LS Vessel
    auto *ls_vessel_s = new G4Box("LsVessel", 0.5 * ls_vessel_x_size, 0.5 * ls_vessel_y_size, 0.5 * ls_vessel_z_size);
    auto *ls_vessel_l = new G4LogicalVolume(ls_vessel_s, black_acrylic, "LsVessel");
    G4VPhysicalVolume *ls_vessel_p = new G4PVPlacement(nullptr, G4ThreeVector(), ls_vessel_l, "LsVessel", buffer_l,
                                                       false, 0, fCheckOverlaps);


    // LS
    auto *ls_s = new G4Box("LS", 0.5 * ls_x_size, 0.5 * ls_y_size, 0.5 * ls_z_size);
    auto *ls_l = new G4LogicalVolume(ls_s, ls, "LS");
    logic_ls_ = ls_l;
    G4VPhysicalVolume *ls_p = new G4PVPlacement(nullptr, G4ThreeVector(), ls_l, "LS", ls_vessel_l, false, 0,
                                                fCheckOverlaps);


    /// 2beta source
    auto *wall_s = new G4Box("Wall", 0.5 * beta_emitter_wall_thickness_, 0.5 * ls_y_size, 0.5 * ls_z_size);
    auto *wall_l = new G4LogicalVolume(wall_s, black_acrylic, "Wall");
    if (beta_emitter_wall_thickness_ > 1e-6) {
        G4VPhysicalVolume *wall_p = new G4PVPlacement(nullptr, G4ThreeVector(), wall_l, "Wall", ls_l, false, 0, true);
        OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::build_geom() :: 2beta source foil deployed" << G4endl;
    }


    /// Horizontal Readout Unit
    auto *hori_ru_s = new G4Tubs("hori_ru_s", 0, hori_ru_radius, 0.5 * hori_fiber_env_length, 0, twopi);
    auto *hori_ru_l = new G4LogicalVolume(hori_ru_s, ls, "z_readout_unit");

    // Horizontal Fiber Envelope
    auto *hori_fiber_env_s = new G4Tubs("hori_fiber_env_s", 0., fiber_radius_, 0.5 * hori_fiber_env_length, 0 * deg,
                                        360 * deg);
    auto *hori_fiber_env_l = new G4LogicalVolume(hori_fiber_env_s, out_cladding_mat, "Horizontal Fiber Envelope");

    //  Horizontal Fiber outer cladding
    auto *hori_fiber_outer_cladding_s = new G4Tubs("hori_fiber_outer_cladding_s", 0., fiber_radius_,
                                                   0.5 * hori_fiber_length, 0 * deg, 360 * deg);
    auto *hori_fiber_outer_cladding_l = new G4LogicalVolume(hori_fiber_outer_cladding_s, out_cladding_mat,
                                                            "Horizontal Fiber Outer Cladding");

    //  Horizontal Fiber inner cladding
    auto *hori_fiber_inner_cladding_s = new G4Tubs("hori_fiber_inner_cladding_s", 0.,
                                                   (1. - outer_cladding_fractional_radius_) * fiber_radius_,
                                                   0.5 * hori_fiber_length, 0 * deg, 360 * deg);
    auto *hori_fiber_inner_cladding_l = new G4LogicalVolume(hori_fiber_inner_cladding_s, polymethylmethacrylate,
                                                            "Horizontal Fiber Inner Cladding");

    // Horizontal Fiber Core
    auto *hori_fiber_core_s = new G4Tubs("hori_fiber_core_s", 0.,
                                         (1. - outer_cladding_fractional_radius_ - inner_cladding_fractional_radius_) *
                                         fiber_radius_,
                                         0.5 * hori_fiber_length, 0 * deg, 360 * deg);
    auto *hori_fiber_core_l = new G4LogicalVolume(hori_fiber_core_s, polystyrene, "Horizontal Fiber Core");

    // SiPMT within Horizontal RU
    auto *hori_ru_sipmt_s = new G4Tubs("hori_ru_sipmt_s", 0., sipmt_radius, 0.5 * sipm_thickness, 0 * deg, 360 * deg);
    auto *hori_ru_sipmt_l = new G4LogicalVolume(hori_ru_sipmt_s, black_acrylic, "Horizontal RU SiPMT");

    // SiPMT Sensitive Volume within Horizontal RU
    auto *hori_ru_sipmt_sens_s = new G4Tubs("hori_ru_sipmt_sens_s", 0., sipmt_sens_radius, (0.5 * sipm_thickness) / 2.,
                                            0 * deg, 360 * deg);
    auto *hori_ru_sipmt_sens_l = new G4LogicalVolume(hori_ru_sipmt_sens_s, polystyrene,
                                                     "Horizontal RU Sensitive SiPMT");
    G4VPhysicalVolume *hori_ru_sipmt_sens_p = new G4PVPlacement(nullptr, G4ThreeVector(0., 0., -0.5 * sipm_thickness +
                                                                                               (0.5 * sipm_thickness) /
                                                                                               2.),
                                                                hori_ru_sipmt_sens_l, "SiPMTsens", hori_ru_sipmt_l,
                                                                false, 0, fCheckOverlaps);


    /// Fill fiber envelops with actual fibers and SiPm
    G4VPhysicalVolume *hori_fiber_outer_cladding_p = new G4PVPlacement(nullptr, G4ThreeVector(),
                                                                       hori_fiber_outer_cladding_l, "OuterCladding",
                                                                       hori_fiber_env_l, false, 0, fCheckOverlaps);
    G4VPhysicalVolume *hori_fiber_inner_cladding_p = new G4PVPlacement(nullptr, G4ThreeVector(),
                                                                       hori_fiber_inner_cladding_l, "InnerCladding",
                                                                       hori_fiber_outer_cladding_l, false, 0,
                                                                       fCheckOverlaps);
    G4VPhysicalVolume *hori_fiber_core_p = new G4PVPlacement(nullptr, G4ThreeVector(), hori_fiber_core_l, "FiberCore",
                                                             hori_fiber_inner_cladding_l, false, 0, fCheckOverlaps);
    G4VPhysicalVolume *hori_ru_sipmt_top_p = new G4PVPlacement(nullptr, G4ThreeVector(0., 0.,
                                                                                      0.5 * hori_fiber_env_length -
                                                                                      0.5 * sipm_thickness),
                                                               hori_ru_sipmt_l, "SiPMT", hori_fiber_env_l, false, 0,
                                                               fCheckOverlaps);

    G4RotationMatrix rot_sipmt;
    rot_sipmt.rotateX(180. * deg);

    G4VPhysicalVolume *hori_ru_sipmt_bot_p = new G4PVPlacement(G4Transform3D(rot_sipmt,
                                                                             G4ThreeVector(0., 0., -0.5 *
                                                                                                   hori_fiber_env_length +
                                                                                                   0.5 *
                                                                                                   sipm_thickness)),
                                                               hori_ru_sipmt_l, "SiPMT",
                                                               hori_fiber_env_l, false, 1, fCheckOverlaps);



    /// Deploy Horizontal Readout Units Along Z
    G4double hori_ru_z0 = -0.5 * hori_ru_distance_ * (number_ru_z_ - 1);
    G4double vert_ru_x0 = -0.5 * vert_ru_distance_ * (number_ru_x_ - 1) + 0.5 * vert_ru_distance_;

    G4RotationMatrix rot_ru_x;
    rot_ru_x.rotateX(90. * deg);

    OFOS_OutputLog::geom_cache << "Horizontal Readout Units (RU position is defined by X,Z couple)" << G4endl;
    OFOS_OutputLog::geom_cache << "RU id\tX[mm]\tZ[mm]" << G4endl;

    for (int ix = 0; ix < number_ru_x_ - 1; ix++) {
        for (int iz = 0; iz < number_ru_z_; iz++) {

            G4double xpos = vert_ru_x0 + static_cast<G4double>(ix) * vert_ru_distance_;
            G4double ypos = 0.;
            G4double zpos = hori_ru_z0 + static_cast<G4double>(iz) * hori_ru_distance_;

            new G4PVPlacement(G4Transform3D(rot_ru_x, G4ThreeVector(xpos, ypos, zpos)),
                              hori_ru_l,
                              "readout_unit",
                              ls_l,
                              false,
                    // iz + ix*vert_ru_x0,
                              number_ru_x_ * number_ru_y_ + iz + ix * number_ru_x_,
                              fCheckOverlaps);

            OFOS_OutputLog::geom_cache << number_ru_x_ * number_ru_y_ + iz + ix * number_ru_x_ << "\t" << xpos << "\t"
                                       << zpos << G4endl;
        }
    }


    /// Fill a readout unit with fiber envelopes
    G4VPhysicalVolume *physFiberEn1 = nullptr;
    G4VPhysicalVolume *physFiberEn2 = nullptr;
    G4VPhysicalVolume *physFiberEn3 = nullptr;
    G4VPhysicalVolume *physFiberEn4 = nullptr;
    G4double r_circ = hori_ru_size_ / sqrt(3.0);


    switch (hori_ru_geometry_) {
        case ReadoutUnitGeometry::SingleFiber:
            physFiberEn1 = new G4PVPlacement(nullptr, G4ThreeVector(), hori_fiber_env_l, "FiberEnv", hori_ru_l, false,
                                             0, fCheckOverlaps);
            break;

        case ReadoutUnitGeometry::Triangle:
            physFiberEn1 = new G4PVPlacement(nullptr, G4ThreeVector(-0.5 * hori_ru_size_, -0.5 * r_circ, 0),
                                             hori_fiber_env_l, "FiberEnv", hori_ru_l, false, 0, fCheckOverlaps);
            physFiberEn2 = new G4PVPlacement(nullptr, G4ThreeVector(0.5 * hori_ru_size_, -0.5 * r_circ, 0),
                                             hori_fiber_env_l, "FiberEnv", hori_ru_l, false, 1, fCheckOverlaps);
            physFiberEn3 = new G4PVPlacement(nullptr, G4ThreeVector(0, r_circ, 0), hori_fiber_env_l, "FiberEnv",
                                             hori_ru_l, false, 2, fCheckOverlaps);
            break;

        case ReadoutUnitGeometry::Square:
            physFiberEn1 = new G4PVPlacement(nullptr, G4ThreeVector(-0.5 * hori_ru_size_, 0.5 * hori_ru_size_, 0),
                                             hori_fiber_env_l, "FiberEnv", hori_ru_l, false, 0, fCheckOverlaps);
            physFiberEn2 = new G4PVPlacement(nullptr, G4ThreeVector(0.5 * hori_ru_size_, 0.5 * hori_ru_size_, 0),
                                             hori_fiber_env_l, "FiberEnv", hori_ru_l, false, 1, fCheckOverlaps);
            physFiberEn3 = new G4PVPlacement(nullptr, G4ThreeVector(-0.5 * hori_ru_size_, -0.5 * hori_ru_size_, 0),
                                             hori_fiber_env_l, "FiberEnv", hori_ru_l, false, 2, fCheckOverlaps);
            physFiberEn4 = new G4PVPlacement(nullptr, G4ThreeVector(0.5 * hori_ru_size_, -0.5 * hori_ru_size_, 0),
                                             hori_fiber_env_l, "FiberEnv", hori_ru_l, false, 3, fCheckOverlaps);
            break;

        default:
            physFiberEn1 = new G4PVPlacement(nullptr, G4ThreeVector(-0.5 * hori_ru_size_, 0.5 * hori_ru_size_, 0),
                                             hori_fiber_env_l, "FiberEnv", hori_ru_l, false, 0, fCheckOverlaps);
            physFiberEn2 = new G4PVPlacement(nullptr, G4ThreeVector(0.5 * hori_ru_size_, 0.5 * hori_ru_size_, 0),
                                             hori_fiber_env_l, "FiberEnv", hori_ru_l, false, 1, fCheckOverlaps);
            physFiberEn3 = new G4PVPlacement(nullptr, G4ThreeVector(-0.5 * hori_ru_size_, -0.5 * hori_ru_size_, 0),
                                             hori_fiber_env_l, "FiberEnv", hori_ru_l, false, 2, fCheckOverlaps);
            physFiberEn4 = new G4PVPlacement(nullptr, G4ThreeVector(0.5 * hori_ru_size_, -0.5 * hori_ru_size_, 0),
                                             hori_fiber_env_l, "FiberEnv", hori_ru_l, false, 3, fCheckOverlaps);
            break;
    }



    /********************  VERTICAL RU  **************************/


    /// Vertical Readout Unit
    auto *vert_ru_s = new G4Tubs("vert_ru_s", 0, vert_ru_radius, 0.5 * vert_fiber_env_length, 0, twopi);
    auto *vert_ru_l = new G4LogicalVolume(vert_ru_s, ls, "readout_unit");

    // Vertical Fiber Envelope
    auto *vert_fiber_env_s = new G4Tubs("vert_fiber_env_s", 0., fiber_radius_, 0.5 * vert_fiber_env_length, 0 * deg,
                                        360 * deg);
    auto *vert_fiber_env_l = new G4LogicalVolume(vert_fiber_env_s, out_cladding_mat, "vert_fiber_env_l");

    //  Fiber outer cladding
    auto *vert_fiber_outer_cladding_s = new G4Tubs("vert_fiber_outer_cladding_s", 0., fiber_radius_,
                                                   0.5 * vert_fiber_length, 0 * deg, 360 * deg);
    auto *vert_fiber_outer_cladding_l = new G4LogicalVolume(vert_fiber_outer_cladding_s, out_cladding_mat,
                                                            "OuterCladding");

    if (OFOS_Verbosity::level > 1)
        OFOS_OutputLog::log_cache << "vert_fiber_outer_cladding_l :: " << vert_fiber_outer_cladding_s->GetOuterRadius()
                                  << G4endl;

    //  Fiber inner cladding
    auto *vert_fiber_inner_cladding_s = new G4Tubs("vert_fiber_inner_cladding_s", 0.,
                                                   (1. - outer_cladding_fractional_radius_) * fiber_radius_,
                                                   0.5 * vert_fiber_length, 0 * deg, 360 * deg);
    auto *vert_fiber_inner_cladding_l = new G4LogicalVolume(vert_fiber_inner_cladding_s, polymethylmethacrylate,
                                                            "vert_fiber_inner_cladding_l");
    if (OFOS_Verbosity::level > 1)
        OFOS_OutputLog::log_cache << "vert_fiber_inner_cladding_l :: " << vert_fiber_inner_cladding_s->GetOuterRadius()
                                  << G4endl;


    //  Fiber Core
    auto *vert_fiber_core_s = new G4Tubs("vert_fiber_core_s", 0.,
                                         (1. - outer_cladding_fractional_radius_ - inner_cladding_fractional_radius_) *
                                         fiber_radius_,
                                         0.5 * vert_fiber_length, 0 * deg, 360 * deg);
    auto *vert_fiber_core_l = new G4LogicalVolume(vert_fiber_core_s, polystyrene, "vert_fiber_core_l");
    if (OFOS_Verbosity::level > 1)
        OFOS_OutputLog::log_cache << "vert_fiber_core_l :: " << vert_fiber_core_s->GetOuterRadius() << G4endl;

    // SiPMT
    auto *vert_ru_sipmt_s = new G4Tubs("vert_ru_sipmt_s", 0., sipmt_radius, 0.5 * sipm_thickness, 0 * deg, 360 * deg);
    auto *vert_ru_sipmt_l = new G4LogicalVolume(vert_ru_sipmt_s, black_acrylic, "SiPMT");

    // SiPMT Sensitive Volume
    auto *vert_ru_sipmt_sens_s = new G4Tubs("vert_ru_sipmt_sens_s", 0., sipmt_sens_radius, (0.5 * sipm_thickness) / 2.,
                                            0 * deg, 360 * deg);
    auto *vert_ru_sipmt_sens_l = new G4LogicalVolume(vert_ru_sipmt_sens_s, polystyrene, "SiPMTsens");
    G4VPhysicalVolume *vert_ru_sipmt_sens_p = new G4PVPlacement(nullptr, G4ThreeVector(0., 0., -0.5 * sipm_thickness +
                                                                                               (0.5 * sipm_thickness) /
                                                                                               2.),
                                                                vert_ru_sipmt_sens_l, "SiPMTsens", vert_ru_sipmt_l,
                                                                false, 0, fCheckOverlaps);

    /// Fill fiber envelops with actual fibers and SiPm
    G4VPhysicalVolume *vert_fiber_outer_cladding_p = new G4PVPlacement(nullptr, G4ThreeVector(),
                                                                       vert_fiber_outer_cladding_l, "OuterCladding",
                                                                       vert_fiber_env_l, false, 0, fCheckOverlaps);
    G4VPhysicalVolume *vert_fiber_inner_cladding_p = new G4PVPlacement(nullptr, G4ThreeVector(),
                                                                       vert_fiber_inner_cladding_l, "InnerCladding",
                                                                       vert_fiber_outer_cladding_l, false, 0,
                                                                       fCheckOverlaps);
    G4VPhysicalVolume *vert_fiber_core_p = new G4PVPlacement(nullptr, G4ThreeVector(), vert_fiber_core_l, "FiberCore",
                                                             vert_fiber_inner_cladding_l, false, 0, fCheckOverlaps);
    G4VPhysicalVolume *vert_ru_sipmt_top_p = new G4PVPlacement(nullptr, G4ThreeVector(0., 0.,
                                                                                      0.5 * vert_fiber_env_length -
                                                                                      0.5 * sipm_thickness),
                                                               vert_ru_sipmt_l, "SiPMT", vert_fiber_env_l, false, 0,
                                                               fCheckOverlaps);

    G4VPhysicalVolume *vert_ru_sipmt_bot_p = new G4PVPlacement(G4Transform3D(rot_sipmt,
                                                                             G4ThreeVector(0., 0., -0.5 *
                                                                                                   vert_fiber_env_length +
                                                                                                   0.5 *
                                                                                                   sipm_thickness)),
                                                               vert_ru_sipmt_l, "SiPMT", vert_fiber_env_l, false, 1,
                                                               fCheckOverlaps);



    /// Deploy Vertical Readout Units Along X and Y
    G4double first_ru_x = -0.5 * vert_ru_distance_ * (number_ru_x_ - 1);
    G4double first_ru_y = -0.5 * vert_ru_distance_ * (number_ru_y_ - 1);

    G4RotationMatrix rot_ru_xy;
    rot_ru_xy.rotateZ(180. * deg);

    OFOS_OutputLog::geom_cache << G4endl << G4endl;
    OFOS_OutputLog::geom_cache << "Vertical Readout Units (RU position is defined by X,Y couple)" << G4endl;
    OFOS_OutputLog::geom_cache << "RU id\tX[mm]\tY[mm]" << G4endl;

    for (int ix = 0; ix < number_ru_x_; ix++) {
        for (int iy = 0; iy < number_ru_y_; iy++) {
            G4double xpos = first_ru_x + static_cast<G4double>(ix) * vert_ru_distance_;
            G4double ypos = first_ru_y + static_cast<G4double>(iy) * vert_ru_distance_;
            G4double zpos = 0.;

            OFOS_OutputLog::geom_cache << iy + ix * number_ru_x_ << "\t" << xpos << "\t" << ypos << G4endl;

            /// triangle arrangement forsees the readout unit to be rotated every other unit
            if (ix % 2 && vert_ru_geometry_ == ReadoutUnitGeometry::Triangle) {
                new G4PVPlacement(G4Transform3D(rot_ru_xy, G4ThreeVector(xpos, ypos, zpos)),
                                  vert_ru_l,
                                  "readout_unit_xy",
                                  ls_l,
                                  false,
                                  iy + ix * number_ru_x_,
                                  fCheckOverlaps);
            } else {
                new G4PVPlacement(nullptr,
                                  G4ThreeVector(xpos, ypos, zpos),
                                  vert_ru_l,
                                  "readout_unit",
                                  ls_l,
                                  false,
                                  iy + ix * number_ru_x_,
                                  fCheckOverlaps);
            }
        }
    }

    /// Fill a readout unit with fiber envelopes
    physFiberEn1 = nullptr;
    physFiberEn2 = nullptr;
    physFiberEn3 = nullptr;
    physFiberEn4 = nullptr;

    r_circ = vert_ru_size_ / sqrt(3.0);

    switch (vert_ru_geometry_) {
        case ReadoutUnitGeometry::SingleFiber:
            physFiberEn1 = new G4PVPlacement(nullptr, G4ThreeVector(), vert_fiber_env_l, "FiberEnv", vert_ru_l, false,
                                             0, fCheckOverlaps);
            break;

        case ReadoutUnitGeometry::Triangle:
            physFiberEn1 = new G4PVPlacement(nullptr, G4ThreeVector(-0.5 * vert_ru_size_, -0.5 * r_circ, 0),
                                             vert_fiber_env_l, "FiberEnv", vert_ru_l, false, 0, fCheckOverlaps);
            physFiberEn2 = new G4PVPlacement(nullptr, G4ThreeVector(0.5 * vert_ru_size_, -0.5 * r_circ, 0),
                                             vert_fiber_env_l, "FiberEnv", vert_ru_l, false, 1, fCheckOverlaps);
            physFiberEn3 = new G4PVPlacement(nullptr, G4ThreeVector(0, r_circ, 0), vert_fiber_env_l, "FiberEnv",
                                             vert_ru_l, false, 2, fCheckOverlaps);
            break;

        case ReadoutUnitGeometry::Square:
            physFiberEn1 = new G4PVPlacement(nullptr, G4ThreeVector(-0.5 * vert_ru_size_, 0.5 * vert_ru_size_, 0),
                                             vert_fiber_env_l, "FiberEnv", vert_ru_l, false, 0, fCheckOverlaps);
            physFiberEn2 = new G4PVPlacement(nullptr, G4ThreeVector(0.5 * vert_ru_size_, 0.5 * vert_ru_size_, 0),
                                             vert_fiber_env_l, "FiberEnv", vert_ru_l, false, 1, fCheckOverlaps);
            physFiberEn3 = new G4PVPlacement(nullptr, G4ThreeVector(-0.5 * vert_ru_size_, -0.5 * vert_ru_size_, 0),
                                             vert_fiber_env_l, "FiberEnv", vert_ru_l, false, 2, fCheckOverlaps);
            physFiberEn4 = new G4PVPlacement(nullptr, G4ThreeVector(0.5 * vert_ru_size_, -0.5 * vert_ru_size_, 0),
                                             vert_fiber_env_l, "FiberEnv", vert_ru_l, false, 3, fCheckOverlaps);
            break;

        default:
            physFiberEn1 = new G4PVPlacement(nullptr, G4ThreeVector(-0.5 * vert_ru_size_, 0.5 * vert_ru_size_, 0),
                                             vert_fiber_env_l, "FiberEnv", vert_ru_l, false, 0, fCheckOverlaps);
            physFiberEn2 = new G4PVPlacement(nullptr, G4ThreeVector(0.5 * vert_ru_size_, 0.5 * vert_ru_size_, 0),
                                             vert_fiber_env_l, "FiberEnv", vert_ru_l, false, 1, fCheckOverlaps);
            physFiberEn3 = new G4PVPlacement(nullptr, G4ThreeVector(-0.5 * vert_ru_size_, -0.5 * vert_ru_size_, 0),
                                             vert_fiber_env_l, "FiberEnv", vert_ru_l, false, 2, fCheckOverlaps);
            physFiberEn4 = new G4PVPlacement(nullptr, G4ThreeVector(0.5 * vert_ru_size_, -0.5 * vert_ru_size_, 0),
                                             vert_fiber_env_l, "FiberEnv", vert_ru_l, false, 3, fCheckOverlaps);
            break;
    }



    /// needed for scoring
    logic_vert_fiber_ = vert_fiber_core_l;
    logic_hori_fiber_ = hori_fiber_core_l;
    logic_vert_ru_sipmt_ = vert_ru_sipmt_sens_l;
    logic_hori_ru_sipmt_ = hori_ru_sipmt_sens_l;
    logic_vessel_ = ls_vessel_l;


    /// visualization attributes
    G4Colour white(G4Colour::White());
    G4Colour grey(G4Colour::Gray());
    G4Colour black(G4Colour::Black());
    G4Colour red(G4Colour::Red());
    G4Colour green(G4Colour::Green());
    G4Colour blue(G4Colour::Blue());
    G4Colour cyan(G4Colour::Cyan());
    G4Colour magenta(G4Colour::Magenta());
    G4Colour yellow(G4Colour::Yellow());

    // Visualization attributes
    auto *vis_white = new G4VisAttributes(white);
    auto *vis_gray = new G4VisAttributes(gray);
    auto *vis_black = new G4VisAttributes(black);
    auto *vis_red = new G4VisAttributes(red);
    auto *vis_green = new G4VisAttributes(green);
    auto *vis_blue = new G4VisAttributes(blue);
    auto *vis_cyan = new G4VisAttributes(cyan);
    auto *vis_magenta = new G4VisAttributes(magenta);
    auto *vis_yellow = new G4VisAttributes(yellow);

    world_l->SetVisAttributes(vis_gray);
    outer_vessel_l->SetVisAttributes(vis_green);
    buffer_l->SetVisAttributes(vis_blue);
    ls_vessel_l->SetVisAttributes(vis_red);
    ls_l->SetVisAttributes(vis_yellow);
    vert_ru_l->SetVisAttributes(vis_cyan);
    hori_ru_l->SetVisAttributes(vis_magenta);
    wall_l->SetVisAttributes(vis_gray);


    auto *fiber_inn_vis_att = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
    auto *fiber_out_vis_att = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
    auto *fiber_cor_vis_att = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));

    vert_fiber_outer_cladding_l->SetVisAttributes(fiber_out_vis_att);
    vert_fiber_inner_cladding_l->SetVisAttributes(fiber_inn_vis_att);
    vert_fiber_core_l->SetVisAttributes(fiber_cor_vis_att);
    vert_fiber_env_l->SetVisAttributes(vis_black);
    vert_ru_sipmt_l->SetVisAttributes(vis_black);

    is_geom_built = true;

    if (OFOS_Verbosity::level > 1)
        OFOS_OutputLog::log_cache << G4endl;

    return world_p;
}




/***********************************************************************************/
/***********************************************************************************/



void
OFOS_DetectorConstruction::print_geom_params() {
    OFOS_OutputLog::log_cache << G4endl << "***    geometry input parameters    ***" << G4endl;

    OFOS_OutputLog::log_cache << "VERTICAL RU (along z) Geometry:      ";
    switch (vert_ru_geometry_) {
        case ReadoutUnitGeometry::SingleFiber:
            OFOS_OutputLog::log_cache << "Single Fiber" << G4endl;
            break;
        case ReadoutUnitGeometry::Triangle:
            OFOS_OutputLog::log_cache << "Triangle (3 Fibers)" << G4endl;
            break;
        case ReadoutUnitGeometry::Square:
            OFOS_OutputLog::log_cache << "Square (4 Fibers)" << G4endl;
            break;
        default:
            OFOS_OutputLog::log_cache << "Undefined" << G4endl;
    }


    OFOS_OutputLog::log_cache << "HORIZONTAL RU (along y) Geometry:    ";
    switch (hori_ru_geometry_) {
        case ReadoutUnitGeometry::SingleFiber:
            OFOS_OutputLog::log_cache << "Single Fiber" << G4endl;
            break;
        case ReadoutUnitGeometry::Triangle:
            OFOS_OutputLog::log_cache << "Triangle (3 Fibers)" << G4endl;
            break;
        case ReadoutUnitGeometry::Square:
            OFOS_OutputLog::log_cache << "Square (4 Fibers)" << G4endl;
            break;
        default:
            OFOS_OutputLog::log_cache << "Undefined" << G4endl;
    }

    OFOS_OutputLog::log_cache << "Number of RU along X:                " << number_ru_x_ << G4endl;
    OFOS_OutputLog::log_cache << "Number of RU along Y:                " << number_ru_y_ << G4endl;
    OFOS_OutputLog::log_cache << "Number of RU along Z:                " << number_ru_z_ << G4endl;
    OFOS_OutputLog::log_cache << "Distance among RU on XY plane:       " << vert_ru_distance_ / mm << " mm" << G4endl;
    OFOS_OutputLog::log_cache << "Size of a RU on the XY plane:        " << vert_ru_size_ / mm << " mm" << G4endl;
    OFOS_OutputLog::log_cache << "Distance among RU along Z direction: " << hori_ru_distance_ / mm << " mm" << G4endl;
    OFOS_OutputLog::log_cache << "Size of a RU along Z direction:      " << hori_ru_size_ / mm << " mm" << G4endl;

    OFOS_OutputLog::log_cache << "Fiber Radius:                        " << fiber_radius_ / mm << " mm" << G4endl;
    OFOS_OutputLog::log_cache << "Fraction of Outer Cladding:          " << outer_cladding_fractional_radius_ << G4endl;
    OFOS_OutputLog::log_cache << "Fraction of Inner Cladding:          " << inner_cladding_fractional_radius_ << G4endl;


}



/***********************************************************************************/
/***********************************************************************************/



void
OFOS_DetectorConstruction::check_geom_params() {
    char line[500];

    if (vert_ru_geometry_ == ReadoutUnitGeometry::Undefined) {
        sprintf(line, "Invalid vert_ru_geometry_(%d)", vert_ru_geometry_);
        G4Exception("OFOS_DetectorConstruction::check_geom_params()",
                    "001",
                    G4ExceptionSeverity::FatalException,
                    line);
    }


    if (hori_ru_geometry_ == ReadoutUnitGeometry::Undefined) {
        sprintf(line, "Invalid ru_z_geometry_(%d)", hori_ru_geometry_);
        G4Exception("OFOS_DetectorConstruction::check_geom_params()",
                    "001",
                    G4ExceptionSeverity::FatalException,
                    line);
    }


    if (not(number_ru_x_ > 0)) {
        sprintf(line, "Invalid number_ru_x_(%d)", number_ru_x_);
        G4Exception("OFOS_DetectorConstruction::check_geom_params()",
                    "001",
                    G4ExceptionSeverity::FatalException,
                    line);
    }

    if (not(number_ru_y_ > 0)) {
        sprintf(line, "Invalid number_ru_y_(%d)", number_ru_y_);
        G4Exception("OFOS_DetectorConstruction::check_geom_params()",
                    "001",
                    G4ExceptionSeverity::FatalException,
                    line);
    }

    if (not(number_ru_z_ > 0)) {
        sprintf(line, "Invalid number_ru_z_(%d)", number_ru_z_);
        G4Exception("OFOS_DetectorConstruction::check_geom_params()",
                    "001",
                    G4ExceptionSeverity::FatalException,
                    line);
    }


    if (not(vert_ru_distance_ > 0.)) {
        sprintf(line, "Invalid vert_ru_distance_(%f mm)", vert_ru_distance_ / mm);
        G4Exception("OFOS_DetectorConstruction::check_geom_params()",
                    "001",
                    G4ExceptionSeverity::FatalException,
                    line);
    }

    if (vert_ru_size_ < 0. || (vert_ru_size_ + fiber_radius_) > vert_ru_distance_) {
        sprintf(line, "Invalid vert_ru_size_(%f mm) - likely exceeding readout unit spacing", vert_ru_size_ / mm);
        G4Exception("OFOS_DetectorConstruction::check_geom_params()",
                    "001",
                    G4ExceptionSeverity::FatalException,
                    line);
    }

    if (not(hori_ru_distance_ > 0.)) {
        sprintf(line, "Invalid hori_ru_distance_(%f mm)", hori_ru_distance_ / mm);
        G4Exception("OFOS_DetectorConstruction::check_geom_params()",
                    "001",
                    G4ExceptionSeverity::FatalException,
                    line);
    }

    if (hori_ru_size_ < 0. || (hori_ru_size_ + fiber_radius_) > hori_ru_distance_) {
        sprintf(line, "Invalid ru_z_size_(%f mm) - likely exceeding readout unit spacing", hori_ru_size_ / mm);
        G4Exception("OFOS_DetectorConstruction::check_geom_params()",
                    "001",
                    G4ExceptionSeverity::FatalException,
                    line);
    }

    if (fiber_radius_ < 0. || fiber_radius_ > hori_ru_size_ || fiber_radius_ > vert_ru_size_) {
        sprintf(line, "Invalid fiber_radius_(%f mm) - likely exceeding readout unit size", fiber_radius_ / mm);
        G4Exception("OFOS_DetectorConstruction::check_geom_params()",
                    "001",
                    G4ExceptionSeverity::FatalException,
                    line);
    }

    if (outer_cladding_fractional_radius_ < 0. || outer_cladding_fractional_radius_ > 1.0) {
        sprintf(line, "Invalid outer_cladding_fractional_radius_(%f mm)", outer_cladding_fractional_radius_ / mm);
        G4Exception("OFOS_DetectorConstruction::check_geom_params()",
                    "001",
                    G4ExceptionSeverity::FatalException,
                    line);
    }

    if (inner_cladding_fractional_radius_ < 0. ||
        (outer_cladding_fractional_radius_ + inner_cladding_fractional_radius_) > 1.0) {
        sprintf(line, "Invalid inner_cladding_fractional_radius_(%f mm) - outer cladding is %f mm",
                inner_cladding_fractional_radius_ / mm,
                outer_cladding_fractional_radius_ / mm);
        G4Exception("OFOS_DetectorConstruction::check_geom_params()",
                    "001",
                    G4ExceptionSeverity::FatalException,
                    line);
    }
}





/***********************************************************************************/
/***********************************************************************************/



void
OFOS_DetectorConstruction::set_ls_dummy_absorption(double value) {
    if (OFOS_Verbosity::level > 1)
        OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::set_ls_dummy_absorption" << G4endl;
    my_ls_properties->set_dummy_absorption(value);

    OFOS_OutputLog::ls_cache << "LS Absorption:                       " << value / m << " m" << G4endl;
}


void
OFOS_DetectorConstruction::set_ls_dummy_scattering(double value) {
    if (OFOS_Verbosity::level > 1)
        OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::set_ls_dummy_scattering" << G4endl;
    my_ls_properties->set_dummy_scattering(value);

    OFOS_OutputLog::ls_cache << "LS Scattering:                       " << value / cm << " cm" << G4endl;
}


void
OFOS_DetectorConstruction::set_ls_density(double value) {
    if (OFOS_Verbosity::level > 1)
        OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::set_ls_density" << G4endl;

    my_ls_properties->set_density(value);
}


void
OFOS_DetectorConstruction::set_ls_birks(double value) {
    if (OFOS_Verbosity::level > 1)
        OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::set_ls_birks" << G4endl;

    my_ls_properties->set_birks_const(value);
}


void
OFOS_DetectorConstruction::set_ls_light_yield(double value) {
    if (OFOS_Verbosity::level > 1)
        OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::set_ls_light_yield" << G4endl;

    my_ls_properties->set_light_yield(value / MeV);
}


void
OFOS_DetectorConstruction::set_ls_loading_fraction(double value) {
    if (OFOS_Verbosity::level > 1)
        OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::set_ls_loading_fraction" << G4endl;

    my_ls_properties->set_loading_fraction(value);
}


void
OFOS_DetectorConstruction::set_ls_loading_material(G4String &value) {
    if (OFOS_Verbosity::level > 1)
        OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::set_ls_loading_material" << G4endl;

    my_ls_properties->set_loading_material(value);
}




/***********************************************************************************/
/***********************************************************************************/





void
OFOS_DetectorConstruction::set_hori_ru_geometry(const G4String &value) {
    if (value == "SingleFiber")
        hori_ru_geometry_ = ReadoutUnitGeometry::SingleFiber;
    else if (value == "Square")
        hori_ru_geometry_ = ReadoutUnitGeometry::Square;
    else if (value == "Triangle")
        hori_ru_geometry_ = ReadoutUnitGeometry::Triangle;
    else {
        G4cerr << "ERROR:: Horizontal Readout Unit Type (along y) not recognized" << G4endl;
        G4cerr << "ERROR:: hori_ru_geometry_ = " << value << G4endl;
        G4Exception("OFOS_DetectorConstruction::set_hori_ru_geometry()",
                    "001",
                    G4ExceptionSeverity::FatalException,
                    "Horizontal Readout Unit Type not recognized");
    }
}


void
OFOS_DetectorConstruction::set_vert_ru_geometry(const G4String &value) {
    if (value == "SingleFiber")
        vert_ru_geometry_ = ReadoutUnitGeometry::SingleFiber;
    else if (value == "Square")
        vert_ru_geometry_ = ReadoutUnitGeometry::Square;
    else if (value == "Triangle")
        vert_ru_geometry_ = ReadoutUnitGeometry::Triangle;
    else {
        G4cerr << "ERROR:: Vertical Readout Unit Type (along z) not recognized" << G4endl;
        G4cerr << "ERROR:: vert_ru_geometry_ = " << value << G4endl;
        G4Exception("OFOS_DetectorConstruction::set_vert_ru_geometry()",
                    "001",
                    G4ExceptionSeverity::FatalException,
                    "Vertical Readout Unit Type not recognized");
    }
}


void
OFOS_DetectorConstruction::set_field() {

    if (OFOS_Verbosity::level > 1)
        OFOS_OutputLog::log_cache << "OFOS_DetectorConstruction::set_field" << G4endl;

    G4ThreeVector fieldValue = G4ThreeVector();
    fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
    fMagFieldMessenger->SetVerboseLevel(1);
}