#include "OFOS_TrackInformation.h"
#include "G4ios.hh"

G4Allocator<OFOS_TrackInformation> aTrackInfoAllocator;


OFOS_TrackInformation::OFOS_TrackInformation(const G4Track*)
{
  saveit = false;
  drawit = false;
  status = active;
}//Ctor


//------------------------------------------------------------------------
void OFOS_TrackInformation::AddTrackStatusFlag(int s)
{
  if(s&active) //track is now active
    status&=~inactive; //remove any flags indicating it is inactive 
  else if(s&inactive) //track is now inactive
    status&=~active; //remove any flags indicating it is active
  status|=s; //add new flags
}//AddTrackStatusFlag

//------------------------------------------------------------------------
void OFOS_TrackInformation::Print() const
{
  G4cout << "TrackInformation : [" << saveit 
	 << "]-[" << drawit 
	 << "]-[" << status << "]"
	 << G4endl;
}//Print


