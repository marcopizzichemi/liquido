#include "OFOS_LsVesselParametrization.h"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4SystemOfUnits.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OFOS_LsVesselParametrization::OFOS_LsVesselParametrization(G4int    nModules,
                                             std::vector<G4double> x,
                                             std::vector<G4double> y,
                                             std::vector<G4double> z,
                                             G4ThreeVector dimensions
                                            )
 : G4VPVParameterisation()
{
  fnModules = nModules;
  fx = x;
  fy = y;
  fz = z;
  fDimensions = dimensions;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OFOS_LsVesselParametrization::~OFOS_LsVesselParametrization()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OFOS_LsVesselParametrization::ComputeTransformation
(const G4int copyNo,
 G4VPhysicalVolume* physVol) const
{
  // // Note: copyNo will start with zero!
  G4double x = fx[copyNo];
  G4double y = fy[copyNo];
  G4double z = fz[copyNo];
  G4ThreeVector origin(x,y,z);
  physVol->SetTranslation(origin);
  physVol->SetRotation(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OFOS_LsVesselParametrization::ComputeDimensions
(G4Box &module,
 const G4int copyNo,
 const G4VPhysicalVolume* physVol) const
{
  module.SetXHalfLength (fDimensions.getX()/2.0);
  module.SetYHalfLength (fDimensions.getY()/2.0);
  module.SetZHalfLength (fDimensions.getZ()/2.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
