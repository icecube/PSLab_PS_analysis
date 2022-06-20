#ifndef BKGSPACEPROB_H_
#define BKGSPACEPROB_H_

#include "llh/public/classes.h"
#include "llh/public/DeclinationDensityMap.h"
#include "TH2D.h"

//#include "llh/public/I3Event.h"

// Forward Declarations (when feasible, more efficient than including headers)
//class I3Event;



class BkgSpaceProb {
 public:
  virtual ~BkgSpaceProb() { }

  virtual BkgSpaceProb* Clone() const = 0;
  // This allows us to copy a derived class object w/o knowing what it is.
  // THIS MUST BE DEFINED SEPARATELY for each derived class.
  // For simple class (not requiring a deep copy) this should suffice:
  // { return new DerivedClass(*this); }

  virtual void FixToBase() = 0;
  virtual void FixToBasePlusEvents(const EventPtrList& evList) = 0;

  virtual double GetBkgProbDensity(const Event& event) const = 0;
  virtual double GetBkgProbDensity(const Coord& coord) const = 0;
};

#endif // BKGSPACEPROB_H_
