#include "llh/public/DecBkgProb.h"

/* Things to take care of:

OLD:  psData.SetDecMapParameters(180,1) for IC40, (180,4) for IC9
NOW:  scripts have to set up BkgSpaceProb first

OLD:  psData.SetBaseEvents would also execute
         baseDecMap_ = DecMapFromCatalogAndEventVector(mapCatalog_, baseEvents_)
NOW:  this has to be done before giving bkgSpaceProb to psData
      (psData should no longer be in charge of setting this)
*/


// recall, ic9 default values: 180 declination maps, smoothing 4 deg
// ic40 default values: 180 declination maps, smoothing 1 deg

void DecBkgProb::Initialize(int nBins, double sigmaSmooth) {
  if (mapCatalog_.GetNMaps() != 0) {
    log_fatal("Error: Cannot re-initialize DecBkgProb.\n");
  }
  mapCatalog_.MakeCatalog(nBins, sigmaSmooth);
  // initialized to zero everywhere
  baseDecMap_.SetNBinsDec(nBins);
  outDecMap_.SetNBinsDec(nBins);
}


void DecBkgProb::FixToBasePlusEvents(const EventPtrList& evList) {
  outDecMap_ = baseDecMap_ + 
    DecMapFromCatalogAndEventPtrList(mapCatalog_, evList);
}
