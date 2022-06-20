#ifndef LLH_DECBKGPROB_H_
#define LLH_DECBKGPROB_H_

#include "llh/public/BkgSpaceProb.h"
#include "llh/public/DeclinationDensityMap.h"
#include "llh/public/I3Event.h"


class DecBkgProb : public BkgSpaceProb {
 private:
  DecMapCatalog mapCatalog_;
  DeclinationDensityMap baseDecMap_; 
  DeclinationDensityMap outDecMap_; // density for baseEvents + any srcEvents

 public:

  DecBkgProb() { }
  DecBkgProb(int nBins, double sigmaSmooth) {
    Initialize(nBins, sigmaSmooth);
  }
  ~DecBkgProb() { }

  // BASE-CLASS FUNCTIONS

  BkgSpaceProb* Clone() const { return new DecBkgProb(*this); }
  // This allows us to copy a derived class object w/o knowing what it is.
  // THIS MUST BE DEFINED SEPARATELY for each derived class.
  // For simple class (not requiring a deep copy) this should suffice:
  // { return new DerivedClass(*this); }

  void FixToBase() { outDecMap_ = baseDecMap_; }

  void FixToBasePlusEvents(const EventPtrList& evList);
    
  double GetBkgProbDensity(const Event& event) const {
    return GetBkgProbDensity(event.GetCoord());
  }

  double GetBkgProbDensity(const Coord& coord) const {
    return outDecMap_.GetProbDensity(coord);
  }


  // SPECIFIC FUNCTIONS FOR DecBkgProb CLASS

  void Initialize(int nMaps, double sigmaSmooth);

  void SetBaseDecMap(const vector<I3Event>& eVect) {
    baseDecMap_ = DecMapFromCatalogAndEventVector(mapCatalog_, eVect);
  }

  // To Do: make a EventPtrList version of this, if/when ready
};


#endif // LLH_DECBKGPROB_H_
