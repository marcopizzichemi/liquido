
#ifndef OFOS_DetectorMessenger_h
#define OFOS_DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

//class OFOS_LsMatProperties;
class OFOS_DetectorConstruction; 
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;

class OFOS_DetectorMessenger: public G4UImessenger
{
  public:
    explicit OFOS_DetectorMessenger( OFOS_DetectorConstruction *det );
    virtual ~OFOS_DetectorMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    
  private:
    OFOS_DetectorConstruction *det_;

    G4UIdirectory*           ofos_directory_;
    G4UIdirectory*           ls_directory_;
    G4UIdirectory*           geom_directory_;
    G4UIdirectory*           detector_directory_;

    /// liquid scintillator commands
    G4UIcmdWithADoubleAndUnit *absorption_length_cmd_;
    G4UIcmdWithADoubleAndUnit *scattering_length_cmd_;
    G4UIcmdWithADouble        *light_yield_cmd_;
    G4UIcmdWithADouble        *birks_cmd_;
    G4UIcmdWithADouble        *density_cmd_;
    G4UIcmdWithADouble        *loading_fraction_cmd_;
    G4UIcmdWithAString        *loading_material_cmd_;
 // G4UIcmdWithoutParameter   *deploy_properties_cmd_;

    /// geometry commands
    G4UIcmdWithAnInteger      *num_ru_x_cmd_      ;
    G4UIcmdWithAnInteger      *num_ru_y_cmd_      ;
    G4UIcmdWithAnInteger      *num_ru_z_cmd_      ;
    G4UIcmdWithAString        *vert_ru_geom_cmd_    ;
    G4UIcmdWithAString        *hori_ru_geom_cmd_     ;
    G4UIcmdWithADoubleAndUnit *vert_ru_distance_cmd_   ;
    G4UIcmdWithADoubleAndUnit *hori_ru_distance_cmd_    ;
    G4UIcmdWithADoubleAndUnit *vert_ru_size_cmd_    ;
    G4UIcmdWithADoubleAndUnit *hori_ru_size_cmd_     ;
    G4UIcmdWithADoubleAndUnit *fiber_rad_cmd_     ;
    G4UIcmdWithADouble        *outer_cladding_cmd_;
    G4UIcmdWithADouble        *inner_cladding_cmd_;
    G4UIcmdWithoutParameter   *build_geometry_cmd_;

};



#endif
