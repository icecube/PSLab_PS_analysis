#ifndef LLH_SIMPLEENERGYPROB_H_
#define LLH_SIMPLEENERGYPROB_H_

#include <vector>

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TTree.h"

#include "llh/public/EnergyProb.h"
#include "llh/public/I3Event.h"


// TO DO:
// Probably just back-filling data but not simulation will not be enough,
// once we start classifying events in different declination bands.
// Need to smooth *both* distributions....


class SimpleEnergyProb : public EnergyProb {
 protected:
  bool optConstrainSignal_;
  bool optLoadModeNew_;

  TH2D hProbGamma_;
  TH1D hProbBkg_;
  TH1D hProbGammaMax_;
  bool rangeIsSet_;
  int energyNBinStartBackFill_;  // lower limit on back-filling


  int GetEnergyBin(double energy) const;

 public:

  SimpleEnergyProb();
  virtual ~SimpleEnergyProb() { }

  // This determines whether E> ELastBin gets moved back to last bin or not
  virtual void SetConstrainSignal(bool opt) { optConstrainSignal_ = opt; }
  virtual void SetLoadModeNew(bool opt) { optLoadModeNew_ = opt; }

  virtual void SetEnergyGammaRangeAndBackFill(
	    int nBinsEnergy, double energyMin, double energyMax,
	    int nBinsGamma, double gammaMin, double gammaMax,
	    int energyNBinStartBackFill);

  TH1D* CalculateEnergyPDF(TH2D* effectiveArea, double gamma, double decMinDeg, double decMaxDeg);
  void SetTableGamma(TH2D *effectiveArea, vector<TH1D*> logEproxy, double decMinDeg, double decMaxDeg);
  virtual void SetTableGamma(TTree *srcTree, TCut cut, TString enString);
  virtual void SetTableBkg(const vector<I3Event>& eventVect);
  virtual void SetTableBkg(const EventPtrList& evList);

  // const Functions

  int GetEnergyBins() const;
  double GetEnergyMin() const;
  double GetEnergyMax() const;

  int GetGammaBins() const;
  double GetGammaMin() const;
  double GetGammaMax() const;

  int GetEnergyNBinStartBackFill() const { return energyNBinStartBackFill_; }

  virtual double GetEnergyProbGamma(const Event& event, double gamma) const ;
  virtual double GetEnergyProbBkg(const Event& event) const ;
  virtual double GetEnergyMaxRatio(const Event& event) const ;

  virtual TH2D GetHistProbGamma() const { return hProbGamma_; }
  virtual TH1D* GetHistProbGamma(double gamma) const;
  virtual TH1D GetHistProbBkg() const { return hProbBkg_; }
  virtual TH1D GetHistProbGammaMax() const { return hProbGammaMax_; }
  
};

// Forward Declarations (when feasible, more efficient than including headers)
class TH1;

void HistFillIn(TH1* h, int nBinStart=1);


#endif // LLH_SIMPLEENERGYPROB_H_
