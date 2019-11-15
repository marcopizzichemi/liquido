#include "OFOS_LsMatProperties.h" 
#include "OFOS_DetectorMessenger.h"
#include "OFOS_Verbosity.h"


OFOS_LsMatProperties::OFOS_LsMatProperties() : n_data_(10),
                                               energy_(nullptr),
                                               ref_index_(nullptr),
                                               fast_spectrum_(nullptr),
                                               slow_spectrum_(nullptr),
                                               scattering_(nullptr),
                                               absorption_(nullptr),

                                               light_yield_(0.),
                                               resolution_scale_(0.),
                                               fast_time_const_(0.),
                                               slow_time_const_(0.),
                                               yield_ratio_(0.),
                                               birks_const_(1e-9),

                                               density_(0.),
                                               loading_fraction_(0.),
                                               loading_material_("")
{
    if(OFOS_Verbosity::level>1)
        G4cout << "OFOS_LsMatProperties :: constructor " << G4endl;


    energy_        = new double[n_data_];
    ref_index_     = new double[n_data_];
    fast_spectrum_ = new double[n_data_];
    slow_spectrum_ = new double[n_data_];
    scattering_    = new double[n_data_];
    absorption_    = new double[n_data_];

 // messenger_     = new OFOS_DetectorMessenger(this);


    set_dummy_energy();
    set_dummy_ref_index    ( 1.50 );
    set_dummy_fast_spectrum( 2 );
    set_dummy_slow_spectrum( 2 );
    set_dummy_scattering   ( 1.0 * m );
    set_dummy_absorption   ( 1.0 * m );
    set_light_yield        ( 10000.0 / MeV );
    set_resolution_scale   ( 1.0 );
    set_fast_time_const    ( 5.0 * ns);
    set_slow_time_const    ( 20.0 * ns );
    set_yield_ratio        ( 0.8 );

}



OFOS_LsMatProperties::~OFOS_LsMatProperties()
{
    if(OFOS_Verbosity::level>1)
        G4cout << "OFOS_LsMatProperties :: destructor " << G4endl;

                     delete energy_;
                  delete ref_index_;
              delete fast_spectrum_;
              delete slow_spectrum_;
                 delete scattering_;
                 delete absorption_;

    energy_        = nullptr;
    ref_index_     = nullptr;
    fast_spectrum_ = nullptr;
    slow_spectrum_ = nullptr;
    scattering_    = nullptr;
    absorption_    = nullptr;
}


void
OFOS_LsMatProperties::set_dummy_energy()
{
    for(int i=0; i<n_data_; ++i)
    {
        energy_[i] = (i+1.0) * eV;
    }
}



void
OFOS_LsMatProperties::set_dummy_ref_index( double value )
{
    for(int i=0; i<n_data_; ++i)
    {
        ref_index_[i] = value;
    }
}



void
OFOS_LsMatProperties::set_dummy_fast_spectrum( int non_empty_bin )
{
    for(int i=0; i<n_data_; ++i)
    {
        fast_spectrum_[i] = (i==non_empty_bin ? 1.0 : 0.0 );
    }
}



void
OFOS_LsMatProperties::set_dummy_slow_spectrum( int non_empty_bin )
{
    for(int i=0; i<n_data_; ++i)
    {
        slow_spectrum_[i] = (i==non_empty_bin ? 1.0 : 0.0 );
    }
}



void
OFOS_LsMatProperties::set_dummy_scattering( double value )
{
    if(OFOS_Verbosity::level>1)
        G4cout << "OFOS_LsMatProperties::set_dummy_scattering( " << value << " )" << G4endl;

    for(int i=0; i<n_data_; ++i)
    {
        scattering_[i] = value;
    }
}



void
OFOS_LsMatProperties::set_dummy_absorption( double value )
{
    if(OFOS_Verbosity::level>1)
        G4cout << "OFOS_LsMatProperties::set_dummy_absorption( " << value << " )" << G4endl;
    for(int i=0; i<n_data_; ++i)
    {
        absorption_[i] = value;
    }
}

