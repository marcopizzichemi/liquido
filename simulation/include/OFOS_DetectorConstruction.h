#ifndef OFOS_DetectorConstruction_h
#define OFOS_DetectorConstruction_h 1

#include <Geant4/G4TransportationManager.hh>
#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4MaterialPropertiesTable.hh"
#include "tls.hh"

#include "OFOS_LsMatProperties.h"
#include "OFOS_DetectorMessenger.h"
#include "OFOS_ReadoutUnitGeometry.h"

//#include "OFOS_SiPmtSD.h"
//#include "OFOS_FiberSD.h"
//#include "OFOS_VesselSD.h"
#include "OFOS_OpticalPhotonSD.h"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4UserLimits;
class G4GlobalMagFieldMessenger;

//class B2bDetectorMessenger;


/// Detector construction class to define materials, geometry
/// and global uniform magnetic field.


class OFOS_DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    OFOS_DetectorConstruction();
    virtual ~OFOS_DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    // Set methods

    void SetCheckOverlaps(G4bool );

    /// define ls composition
    G4Material* make_cocktail( double density, double loading_fraction, G4String& material );

    /// act on OFOS_LsMatProperties
    void set_ls_dummy_absorption(double value); 
    void set_ls_dummy_scattering(double value); 
    void set_ls_birks           ( double value );
    void set_ls_light_yield     ( double value );
    void set_ls_density         ( double value );
    void set_ls_loading_fraction( double value );
    void set_ls_loading_material( G4String & value );

    /// act on G4MaterialPropertiesTable
    void deploy_ls_properties();

    /// setters
    void set_vert_ru_geometry ( const G4String& value );
    void set_hori_ru_geometry ( const G4String& value );
    void set_number_ru_x      ( G4int    value ) { number_ru_x_      = value; }
    void set_number_ru_y      ( G4int    value ) { number_ru_y_      = value; }
    void set_number_ru_z      ( G4int    value ) { number_ru_z_      = value; }
    void set_vert_ru_dist     ( G4double value ) { vert_ru_distance_ = value; }  
    void set_vert_ru_size     ( G4double value ) { vert_ru_size_     = value; }  
    void set_hori_ru_dist     ( G4double value ) { hori_ru_distance_ = value; }  
    void set_hori_ru_size     ( G4double value ) { hori_ru_size_     = value; }  
    void set_fiber_radius     ( G4double value ) { fiber_radius_     = value; }

    void set_outer_cladding_fractional_radius  (G4double value ) { outer_cladding_fractional_radius_ = value; }
    void set_inner_cladding_fractional_radius  (G4double value ) { inner_cladding_fractional_radius_ = value; }

    void update_geom();

    inline G4Material* get_ls_material() { return (logic_ls_ ? logic_ls_->GetMaterial() : nullptr);}


  private:
    // methods
    void define_materials();
 // G4VPhysicalVolume* DefineVolumes();
    G4VPhysicalVolume* build_geom();
    void check_geom_params();
    void print_geom_params();
  
    // data members
    G4LogicalVolume*   logic_vert_ru_sipmt_;    // pointer to the logical SiPmt for scoring
    G4LogicalVolume*   logic_hori_ru_sipmt_;    // pointer to the logical SiPmt for scoring
    G4LogicalVolume*   logic_vert_fiber_;    // pointer to the logical Fiber for scoring
    G4LogicalVolume*   logic_hori_fiber_;    // pointer to the logical Fiber for scoring
    G4LogicalVolume*   logic_vessel_;   // pointer to the logical Vessel for scoring
    G4LogicalVolume*   logic_ls_;   // pointer to the logical LS to allow live material modification

   /***************
    *  MATERIALS  *
    ***************/

    G4Material                *ls; /// need to be data member in order to change its properties at run time
    G4Material                *lab;
    OFOS_LsMatProperties      *my_ls_properties;
    G4MaterialPropertiesTable *ls_mpt;

    G4Material *out_cladding_mat;
    G4Material *polymethylmethacrylate; // fiber inner cladding
    G4Material *polystyrene; // fiber core
    G4Material *black_acrylic;
    G4Material *air;

   /***************
    *  GEOMETRY   *
    ***************/
    G4bool is_geom_built;
    ReadoutUnitGeometry vert_ru_geometry_;
    ReadoutUnitGeometry hori_ru_geometry_;
    G4int               number_ru_x_;
    G4int               number_ru_y_;
    G4int               number_ru_z_;
    G4double            vert_ru_distance_;  // csi
    G4double            vert_ru_size_;   // psi (fiber dist within a readout unit)
    G4double            hori_ru_distance_;  // csi
    G4double            hori_ru_size_;   // psi (fiber dist within a readout unit)
    G4double            fiber_radius_;
    G4double            outer_cladding_fractional_radius_;
    G4double            inner_cladding_fractional_radius_;
    G4double            beta_emitter_wall_thickness_;
    G4double            air_buffer_thickness_;
    G4double            outer_vessel_thickness_;
    G4double            ls_vessel_thickness_;


   /**************************
    *  SENSITIVE DETECTORS   *
    **************************/

    OFOS_OpticalPhotonSD *vert_sipmt_SD_;
    OFOS_OpticalPhotonSD *vert_fiber_SD_;
    OFOS_OpticalPhotonSD *hori_sipmt_SD_;
    OFOS_OpticalPhotonSD *hori_fiber_SD_;
    OFOS_OpticalPhotonSD *vessel_SD_; 

 // OFOS_SiPmtSD  *vert_sipmt_SD_;
 // OFOS_FiberSD  *vert_fiber_SD_;
 // OFOS_SiPmtSD  *hori_sipmt_SD_;
 // OFOS_FiberSD  *hori_fiber_SD_;
 // OFOS_VesselSD *vessel_SD_; 

 // OFOS_OpticalPhotonSD *optical_photon_sd_;


    OFOS_DetectorMessenger *messenger_;

    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger;
                                         // magnetic field messenger


    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps 

    void set_field();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
