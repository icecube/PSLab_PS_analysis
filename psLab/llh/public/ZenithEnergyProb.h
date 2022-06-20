#ifndef LLH_ZENITHENERGYPROB_H_
#define LLH_ZENITHENERGYPROB_H_

#include <vector>

#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TString.h"

#include "llh/public/I3Event.h"
#include "llh/public/SimpleEnergyProb.h"



class ZenithEnergyProb : public SimpleEnergyProb {
 protected:
  // THESE INHERITED DATA MEMBERS ARE USED FOR SETTING INDIVIDUAL EPROBS

  //  bool optConstrainSignal_;
  //  bool optLoadModeNew_;

  //  TH2D hProbGamma_;
  //  TH1D hProbBkg_;
  //  TH1D hProbGammaMax_;
  //  bool rangeIsSet_;
  //  int energyNBinStartBackFill_;  // lower limit on back-filling


  //  int GetEnergyBin(double energy) const;

  TH1D hZenDegBands_;
  double sourceZenWidthDeg_;
  vector<SimpleEnergyProb> eProbVect_;
  TString recoZenRadName_;
  int selectZenithBand_;

 public:

  ZenithEnergyProb();
  virtual ~ZenithEnergyProb() { }

  // this extends the boundaries of each zenith range when gamma table
  // is being filled; since the EventLoader may select signal events up to this
  // far away, and relocate them inside the nominal range as injected signal
  void SetSourceZenWidthDeg(double zenWidthDeg) {
    sourceZenWidthDeg_ = zenWidthDeg;
  }

  void SetZenithBandsDeg(const vector<double>& zenMinDegVect);
  void SetZenithBandsRad(const vector<double>& zenMinRadVect);

  // Create vector of SimpleEnergyProbs, and pass along the settings
  void CreateBands();

  // This determines whether E> ELastBin gets moved back to last bin or not
  virtual void SetConstrainSignal(bool opt);
  virtual void SetLoadModeNew(bool opt);

  virtual void SetEnergyGammaRangeAndBackFill(
	    int nBinsEnergy, double energyMin, double energyMax,
	    int nBinsGamma, double gammaMin, double gammaMax,
	    int energyNBinStartBackFill);

  void SetName_recoZenith_rad(TString varexp) { recoZenRadName_ = varexp; }
  void SetTableGamma(TH2D* Aeff, vector< vector<TH1D*> > logEproxy);
  virtual void SetTableGamma(TTree *srcTree, TCut cut, TString enString);
  virtual void SetTableBkg(const vector<I3Event>& eventVect);
  virtual void SetTableBkg(const EventPtrList& evList);

  // base class functions to retrieve histograms will operate on this band
  void SelectZenithBand(int band);  

  // const Functions

  int GetZenDegBand(double zenDeg) const;
  int GetZenDegBand(const Event& event) const;
  int GetNZenithBands() const { return hZenDegBands_.GetNbinsX(); }
  int ValidateZenithBand(int band) const ;

  //  int GetEnergyBins() const;
  //  double GetEnergyMin() const;
  //  double GetEnergyMax() const;
  //
  //  int GetGammaBins() const;
  //  double GetGammaMin() const;
  //  double GetGammaMax() const;
  //
  //  int GetEnergyNBinStartBackFill() const { return energyNBinStartBackFill_; }

  virtual double GetEnergyProbGamma(const Event& event, double gamma) const ;
  virtual double GetEnergyProbBkg(const Event& event) const ;
  virtual double GetEnergyMaxRatio(const Event& event) const ;

  virtual TH2D GetHistProbGamma() const { 
    //Printf("selected zenith band = %d", selectZenithBand_);
    return eProbVect_[selectZenithBand_].GetHistProbGamma(); 
  }
  virtual TH2D GetHistProbGamma(int zenBand) const {
    assert(zenBand<=eProbVect_.size());
    //Printf("selected zenith band = %d", zenBand);
    return eProbVect_[zenBand].GetHistProbGamma();
  }
  virtual TH1D* GetHistProbGamma(double gamma) const {
    return eProbVect_[selectZenithBand_].GetHistProbGamma(gamma);
  }
  virtual TH1D GetHistProbBkg() const {
    return eProbVect_[selectZenithBand_].GetHistProbBkg();
  }
  virtual TH1D GetHistProbGammaMax() const {
    return eProbVect_[selectZenithBand_].GetHistProbGammaMax();
  }
  

  TH1D GetHistZenDegBands() const { return hZenDegBands_; }

  const SimpleEnergyProb& GetSimpleEnergyProb(int band) const {
    return eProbVect_[ValidateZenithBand(band)];
  }

};


#endif // LLH_ZENITHENERGYPROB_H_
