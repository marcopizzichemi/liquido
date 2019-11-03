#include "OFOS_OutputNtuples.h"


OFOS_OutputNtuples::OFOS_OutputNtuples(const char *h_name, const char *h_title, int n_max_hits, const char *i_name, const char *i_title, int n_max_interactions) :
        n_max_hits_(n_max_hits),
        n_hits_(0),
        is_nh_overflow_(false),
        collection_eff_(0.0),
        n_max_interactions_(n_max_interactions),
        n_interactions_(0),
        is_ni_overflow_(false),
        is_evt_contained_(true),

        tot_en_dep_(static_cast<float>(0.)),
        ion_en_in_ls_(static_cast<float>(0.)),
        primary_ionization_(static_cast<float>(0.)),
//  primary_ionization_incl_deltas_(0),
        en_loss_via_deltas_(static_cast<float>(0.)),
        radiative_en_loss_ (static_cast<float>(0.)),

        n_scint_phot_(0),
        n_chere_phot_(0),
        n_lost_phot_(0)
{
    hit_tree_   = new TTree(h_name, h_title);
    truth_tree_ = new TTree(i_name, i_title);

    init_arrays();
    set_branches();
}

void
OFOS_OutputNtuples::write()
{
    hit_tree_->Write();
    truth_tree_->Write();
}


void
OFOS_OutputNtuples::init_arrays()
{
    h_track_id_    = new int    [n_max_hits_];
    h_parent_id_   = new int    [n_max_hits_];
    h_sd_type_     = new int    [n_max_hits_];
    h_primary_id_  = new int    [n_max_hits_];
    h_secondary_id_= new int    [n_max_hits_];
    h_gen_proc_    = new int    [n_max_hits_];
    h_time_        = new double [n_max_hits_];
    h_wavelength_  = new float  [n_max_hits_];
    h_pos_x_       = new float  [n_max_hits_];
    h_pos_y_       = new float  [n_max_hits_];
    h_pos_z_       = new float  [n_max_hits_];
    h_pol_x_       = new float  [n_max_hits_];
    h_pol_y_       = new float  [n_max_hits_];
    h_pol_z_       = new float  [n_max_hits_];


    t_id_          = new int   [n_max_interactions_];
    p_id_          = new int   [n_max_interactions_];
    i_id_          = new int   [n_max_interactions_];
    i_pos_x_       = new float [n_max_interactions_];
    i_pos_y_       = new float [n_max_interactions_];
    i_pos_z_       = new float [n_max_interactions_];
    i_particle_    = new int   [n_max_interactions_];
    i_dE_          = new float [n_max_interactions_];
    i_time_        = new float [n_max_interactions_];

}


void
OFOS_OutputNtuples::set_branches()
{
    hit_tree_->Branch( "n_hits"         , &n_hits_         , "n_hits/I");
    hit_tree_->Branch( "is_nh_overflow" , &is_nh_overflow_ , "is_nh_overflow/O");
    hit_tree_->Branch( "h_track_id"     , h_track_id_      , "h_track_id[n_hits]/I");
    hit_tree_->Branch( "h_parent_id"    , h_parent_id_     , "h_parent_id[n_hits]/I");
    hit_tree_->Branch( "h_sd_type"      , h_sd_type_       , "h_sd_type[n_hits]/I");
    hit_tree_->Branch( "h_primary_id"   , h_primary_id_    , "h_primary_id[n_hits]/I");
    hit_tree_->Branch( "h_secondary_id" , h_secondary_id_  , "h_secondary_id[n_hits]/I");
    hit_tree_->Branch( "h_gen_proc"     , h_gen_proc_      , "h_gen_proc[n_hits]/I");
    hit_tree_->Branch( "h_time"         , h_time_          , "h_time[n_hits]/D");
    hit_tree_->Branch( "h_wavelength"   , h_wavelength_    , "h_wavelength[n_hits]/F");
    hit_tree_->Branch( "h_pos_x"        , h_pos_x_         , "h_pos_x[n_hits]/F");
    hit_tree_->Branch( "h_pos_y"        , h_pos_y_         , "h_pos_y[n_hits]/F");
    hit_tree_->Branch( "h_pos_z"        , h_pos_z_         , "h_pos_z[n_hits]/F");
    hit_tree_->Branch( "h_pol_x"        , h_pol_x_         , "h_pol_x[n_hits]/F");
    hit_tree_->Branch( "h_pol_y"        , h_pol_y_         , "h_pol_y[n_hits]/F");
    hit_tree_->Branch( "h_pol_z"        , h_pol_z_         , "h_pol_z[n_hits]/F");



    truth_tree_->Branch( "n_interactions"    , &n_interactions_     , "n_interactions/I");
    truth_tree_->Branch( "is_ni_overflow"    , &is_ni_overflow_     , "is_ni_overflow/O");
    truth_tree_->Branch( "tot_en_dep"        , &tot_en_dep_         , "tot_en_dep_/F"      );
    truth_tree_->Branch( "primary_ionization", &primary_ionization_ , "primary_ionization_/F"      );
    truth_tree_->Branch( "en_loss_via_deltas", &en_loss_via_deltas_ , "en_loss_via_deltas_/F");
    truth_tree_->Branch( "radiative_en_loss",  &radiative_en_loss_  , "radiative_en_loss_/F");
//  truth_tree_->Branch( "ion_en_in_ls"      , &ion_en_in_ls_       , "ion_en_in_ls_/F"      );
    truth_tree_->Branch( "n_scint_phot"      , &n_scint_phot_       , "n_scint_phot/I");
    truth_tree_->Branch( "n_chere_phot"      , &n_chere_phot_       , "n_chere_phot/I");
    truth_tree_->Branch( "track_id"          , t_id_                , "track_id[n_interactions]/I"      );
    truth_tree_->Branch( "parent_id"         , p_id_                , "parent_id[n_interactions]/I"      );
    truth_tree_->Branch( "interaction_id"    , i_id_                , "interaction_id[n_interactions]/I"      );
    truth_tree_->Branch( "i_particle"        , i_particle_          , "i_particle[n_interactions]/I");
    truth_tree_->Branch( "i_pos_x"           , i_pos_x_             , "i_pos_x[n_interactions]/F"   );
    truth_tree_->Branch( "i_pos_y"           , i_pos_y_             , "i_pos_y[n_interactions]/F"   );
    truth_tree_->Branch( "i_pos_z"           , i_pos_z_             , "i_pos_z[n_interactions]/F"   );
    truth_tree_->Branch( "i_dE"              , i_dE_                , "i_dE[n_interactions]/F"      );
    truth_tree_->Branch( "i_time"            , i_time_              , "i_time[n_interactions]/F"    );
    truth_tree_->Branch( "is_evt_contained"  , &is_evt_contained_   , "is_evt_contained_/O"         );
}



OFOS_OutputNtuples::~OFOS_OutputNtuples()
{
    G4cout << "Deleting OFOS_OutputNtuples" << G4endl;

    delete h_track_id_    ;
    delete h_parent_id_   ;
    delete h_sd_type_     ;
    delete h_primary_id_  ;
    delete h_secondary_id_;
    delete h_gen_proc_    ;
    delete h_time_        ;
    delete h_wavelength_  ;
    delete h_pos_x_       ;
    delete h_pos_y_       ;
    delete h_pos_z_       ;
    delete h_pol_x_       ;
    delete h_pol_y_       ;
    delete h_pol_z_       ;

    delete t_id_;
    delete p_id_;
    delete i_id_;
    delete i_pos_x_;
    delete i_pos_y_;
    delete i_pos_z_;
    delete i_particle_;
    delete i_dE_;
    delete i_time_;

    delete hit_tree_ ;
    delete truth_tree_ ;

    G4cout << "Deleting OFOS_OutputNtuples >>> DONE" << G4endl;
}



bool
OFOS_OutputNtuples::fill_interaction( int track_id, int parent_id, int interaction_id, int particle,  const G4ThreeVector &pos, float dE, float time )
{
    if(n_interactions_ == n_max_interactions_)
    {
        is_ni_overflow_ = true;
        return false;
    }

    t_id_      [n_interactions_] = track_id;
    i_id_      [n_interactions_] = interaction_id;
    p_id_      [n_interactions_] = parent_id;
    i_particle_[n_interactions_] = particle;
    i_pos_x_   [n_interactions_] = static_cast<float>(pos.x() / mm);
    i_pos_y_   [n_interactions_] = static_cast<float>(pos.y() / mm);
    i_pos_z_   [n_interactions_] = static_cast<float>(pos.z() / mm);
    i_dE_      [n_interactions_] = static_cast<float>(dE / MeV);
    i_time_    [n_interactions_] = static_cast<float>(time / ns);


    n_interactions_++;

    return true;
}



bool
OFOS_OutputNtuples::fill_hit( const OFOS_OPHit *a_hit )
{
    if(n_hits_ == n_max_hits_)
    {
        is_nh_overflow_ = true;
        return false;
    }

    h_track_id_    [n_hits_] = static_cast<int>   (a_hit->get_track_id());
    h_parent_id_   [n_hits_] = static_cast<int>   (a_hit->get_parent_id());
    h_sd_type_     [n_hits_] = static_cast<int>   (a_hit->get_sd_type());
    h_primary_id_  [n_hits_] = static_cast<int>   (a_hit->get_primary_id());
    h_secondary_id_[n_hits_] = static_cast<int>   (a_hit->get_secondary_id());
    h_gen_proc_    [n_hits_] = static_cast<int>   (a_hit->get_gen_proc());
    h_time_        [n_hits_] = static_cast<double>(a_hit->get_time() / ns);
    h_wavelength_  [n_hits_] = static_cast<float> (a_hit->get_wavelength());
    h_pos_x_       [n_hits_] = static_cast<float> (a_hit->get_position().x() / mm);
    h_pos_y_       [n_hits_] = static_cast<float> (a_hit->get_position().y() / mm);
    h_pos_z_       [n_hits_] = static_cast<float> (a_hit->get_position().z() / mm);
    h_pol_x_       [n_hits_] = static_cast<float> (a_hit->get_polarization().x());
    h_pol_y_       [n_hits_] = static_cast<float> (a_hit->get_polarization().y());
    h_pol_z_       [n_hits_] = static_cast<float> (a_hit->get_polarization().z());

    n_hits_++;

    return true;
}

bool
OFOS_OutputNtuples::fill_all_branches( double value )
{
    if(n_hits_ == n_max_hits_)
    {
        is_nh_overflow_ = true;
        return false;
    }

    h_track_id_    [n_hits_] = static_cast<int>   (value);
    h_parent_id_   [n_hits_] = static_cast<int>   (value);
    h_sd_type_     [n_hits_] = static_cast<int>   (value);
    h_primary_id_  [n_hits_] = static_cast<int>   (value);
    h_secondary_id_[n_hits_] = static_cast<int>   (value);
    h_gen_proc_    [n_hits_] = static_cast<int>   (value);
    h_time_        [n_hits_] = value;
    h_wavelength_  [n_hits_] = static_cast<float> (value);
    h_pos_x_       [n_hits_] = static_cast<float> (value);
    h_pos_y_       [n_hits_] = static_cast<float> (value);
    h_pos_z_       [n_hits_] = static_cast<float> (value);
    h_pol_x_       [n_hits_] = static_cast<float> (value);
    h_pol_y_       [n_hits_] = static_cast<float> (value);
    h_pol_z_       [n_hits_] = static_cast<float> (value);

    n_hits_++;

    return true;
}


void
OFOS_OutputNtuples::store_event()
{

    hit_tree_->Fill();
//  G4cout << "store_event() :: n_scin: " << n_scint_phot_ << "   n_chere: " << n_chere_phot_ 
//         << "  n_lost: " << n_lost_phot_ << G4endl;
    n_hits_ = 0;
    // n_lost_phot_ = 0;

    truth_tree_->Fill();
//  G4cout << "tot_edep = " << tot_en_dep_ << G4endl
//         << "prim_ion = " << primary_ionization_ << G4endl
//         << "p_ion_wd = " << primary_ionization_incl_deltas_ << G4endl 
//         << "ion_ls   = " << ion_en_in_ls_ << G4endl
//         << "delta_en = " << en_loss_via_deltas_ << G4endl
//         << "rad_en   = " << radiative_en_loss_ << G4endl
//         << "sum(idr) = " << primary_ionization_ + en_loss_via_deltas_ + radiative_en_loss_ << G4endl;

    n_interactions_ = 0;
    n_scint_phot_ = 0;
    n_chere_phot_ = 0;
    tot_en_dep_ = 0;
    ion_en_in_ls_ = 0;
    primary_ionization_ = 0;
//  primary_ionization_incl_deltas_ = 0;
    en_loss_via_deltas_ = 0;
    radiative_en_loss_ = 0;

    n_lost_phot_ = 0;
    is_ni_overflow_ = false;
    is_evt_contained_ = true;


}
