#ifndef PLOTSCRIPTS_EFFECTIVEAREA_H_
#define PLOTSCRIPTS_EFFECTIVEAREA_H_

#include "TH1D.h"
#include "TString.h"
#include "TTree.h"


class EffectiveArea {
 public:
  EffectiveArea() : nBins_(40) , logEMin_(1.) , logEMax_(9.) { }

  ~EffectiveArea() { }

  void SetEnergyRange(int nBins, double logEMin, double logEMax);

  TH1D* Calculate(TTree *tree, double zenMinDeg, double zenMaxDeg, 
		  TCut cut, TString hname="hEffArea");

 private:
  int nBins_;
  double logEMin_;
  double logEMax_;
};



#endif // PLOTSCRIPTS_EFFECTIVEAREA_H_
