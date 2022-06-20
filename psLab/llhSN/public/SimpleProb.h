#ifndef LLH_SIMPLEPROB_H_
#define LLH_SIMPLEPROB_H_

#include <vector>

#include "TH1D.h"
//#include "TH2D.h"
#include "TString.h"
#include "TTree.h"

#include "llhSN/public/SNEvent.h"


// What is this?
// Made a basic container class to hold onto histogram which
// Supernova triggers can have. It is designed to take in the
// strength of the alert and deliver the fraction of alerts
// with strength above that seen. For SN alerts, we presently
// have no signal hypothesis.


class SimpleProb {
 protected:

  //bool optConstrainSignal_;
  bool rangeIsSet_;

//  TH2D hProbGamma_;
  TH1D hProbBkg_;
  int nBinStartBackFill_;  // lower limit on back-filling


  int GetBin(double energy) const;

 public:

  SimpleProb() { };
  virtual ~SimpleProb() { }
  
  void SetRangeAndBackFill(int nBins, double sMin, double sMax,
	 int nBinStartBackFill);

  // This determines whether E> ELastBin gets moved back to last bin or not
  //virtual void SetConstrainSignal(bool opt) { optConstrainSignal_ = opt; }

  virtual void SetTableBkg(const vector<SNEvent> eventVect);
  virtual void SetTableBkg(const EventPtrList& evList);

  // const Functions

  int GetBins() const;
  double GetMin() const;
  double GetMax() const;

  int GetNBinStartBackFill() const { return nBinStartBackFill_; }

  double GetProbBkg(const SNEvent * event) const ;
  double GetProbBkg(const Event * event) const {
    return GetProbBkg( dynamic_cast<const SNEvent*> (event));
  }
  
  TH1D GetHistProbBkg() const { return hProbBkg_; }
  
};

// Forward Declarations (when feasible, more efficient than including headers)
class TH1;

void HistFillIn(TH1* h, int nBinStart=1);


#endif // LLH_SIMPLEPROB_H_
