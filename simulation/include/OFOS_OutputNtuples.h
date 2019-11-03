#ifndef OFOS_OutputNtuples_h
#define OFOS_OutputNtuples_h 1

#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "TTree.h"
#include "OFOS_OPHit.h"

class OFOS_OutputNtuples
{
    public:
        OFOS_OutputNtuples(const char *h_name, const char *h_title, int n_max_hits, const char *i_name, const char *i_title, int n_max_interactions); 
        ~OFOS_OutputNtuples();
        void init_arrays();
        void set_branches();
        void write() ;
        bool fill_all_branches( double value);
        void store_event() ;

        bool fill_hit( const OFOS_OPHit *a_hit );
        bool fill_interaction( int track_id, int parent_id, int interaction_id, int particle,  const G4ThreeVector &pos, float dE, float time ); 

        inline void add_edep( double e )    { tot_en_dep_   += (e/MeV); };
        inline void add_ls_ion (double e)   { ion_en_in_ls_       += (e/MeV); };
        inline void add_prim_ion( double e ){ primary_ionization_ += (e/MeV); }; 
//      inline void add_prim_ion_incl_deltas( double e ){ primary_ionization_incl_deltas_ += (e/MeV); }; 
        inline void add_en_loss_via_deltas(double e)    { en_loss_via_deltas_  += (e/MeV); };
        inline void add_radiative_en_loss (double e)    { radiative_en_loss_   += (e/MeV); };

        inline void add_scint_phot( int n ) { n_scint_phot_ += n;};
        inline void add_chere_phot( int n ) { n_chere_phot_ += n;};
        inline void add_lost_phot ( int n ) { n_lost_phot_ += n;};
        inline void set_lost_containment()  { is_evt_contained_ = false; };

        TTree *hit_tree_;
        TTree *truth_tree_;


    private:
        float       tot_en_dep_;
        float       ion_en_in_ls_;
        float       primary_ionization_;
//      float       primary_ionization_incl_deltas_;
        float       en_loss_via_deltas_;
        float       radiative_en_loss_;
        int         n_scint_phot_;
        int         n_chere_phot_;
        int         n_lost_phot_;
        int         n_max_interactions_;
        int         n_interactions_;
        bool        is_ni_overflow_;
        int        *t_id_;
        int        *p_id_;
        int        *i_id_;
        int        *i_particle_;
        float      *i_pos_x_;
        float      *i_pos_y_;
        float      *i_pos_z_;
        float      *i_dE_;
        float      *i_time_;
        bool        is_evt_contained_;


        float       collection_eff_;
        int         n_max_hits_;
        bool        is_nh_overflow_;
        int         n_hits_;
        int        *h_track_id_    ;
        int        *h_parent_id_;  // Josh addition
        int        *h_sd_type_     ;
        int        *h_primary_id_  ;
        int        *h_secondary_id_;
        int        *h_gen_proc_;
        double     *h_time_        ;
        float      *h_wavelength_;
        float      *h_pos_x_;
        float      *h_pos_y_;
        float      *h_pos_z_;
        float      *h_pol_x_;
        float      *h_pol_y_;
        float      *h_pol_z_;
};


#endif
