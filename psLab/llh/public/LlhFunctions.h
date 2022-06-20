#ifndef LLH_LLHFUNCTIONS_H_
#define LLH_LLHFUNCTIONS_H_

#include "TProfile.h"

#include "llh/public/classes.h"
#include "llh/public/LlhFunctionsBase.h"



class BinnedAnalysis : public AnalysisLlh {
 private:
  double binSize_;
  double psfFractionInBin_;
  TProfile* profileLogProb_;

  vector<double> distanceVector_;
  int counts_;
  double bkgDensity_;
  double nBkgMean_;
  double logLambdaBest_;
  double nSrcMin_;
  double nSrcBest_;

 public:

  BinnedAnalysis();
  virtual ~BinnedAnalysis() { }


  void SetBinSize(double binSize) { 
    binSize_ = binSize; 
    psfFractionInBin_ = 1.;
  }
  void SetBinSizeAndPsfFractionInBin(double binSize, double fraction) { 
    binSize_ = binSize;
    psfFractionInBin_ = fraction;
  }

  void SetNSrcMin(double nSrcMin) { nSrcMin_ = nSrcMin; }

  // Optional: if you want stats for a range of bin sizes
  void SetTProfile(TProfile* profileLogProb) 
    { profileLogProb_ = profileLogProb; }


  void PrepareAnalysis();

  void DoAnalysis();

  // for bins, it's not really maximizing!  it's just doing the analysis
  void MaximizeLlh() { DoAnalysis(); }

  double EvaluateLlh(double nSrc);
  virtual double EvaluateLlh(double *parValueArray) {
    return EvaluateLlh(parValueArray[0]);
  }
    

  double GetTestStatistic() const { return counts_; }
  //  double GetTestStatistic() const { return Get_counts(); }
  double GetEstProb() const { return GetPoissonProb(); }
  double GetPar(int i) const { 
    if (i==0) { return nSrcBest_; }
    else { return 0.; }
  }

  double Get_logLambdaBest() const {return logLambdaBest_;}
  double Get_binSize() const {return binSize_;}
  double Get_counts() const {return counts_;}
  double Get_bkgDensity() const {return bkgDensity_;}
  double Get_nBkgMean() const {return nBkgMean_;}
  double GetPoissonProb() const;

};




class SimpleLlh : public AnalysisLlh {

 private:
  enum Method { SIMPLE, MINUIT_MIGRAD };
  Method method_;

  vector<double> ProbRatio_;

  double nSrcInc_;
  double nSrcMin_;
  double nSrcMax_;
  bool opt_Default_SrcRange_;

  double logLambdaBest_;
  double chiSq_;
  double chiSqProb_;
  double nSrcBest_;
  int nEvents_;

  virtual void EvalMinuitFCN(int &npar, double *gin, double &f, 
			     double *par, int iflag);

  void MaximizeLlh_Simple();
  void MaximizeLlh_MinuitMigrad();

 public:

  SimpleLlh();
  virtual ~SimpleLlh() { }

  void SetMethodSimple() { method_ = SIMPLE;}
  void SetMethodMinuitMigrad() { method_ = MINUIT_MIGRAD;}

  void Set_nSrcInc(double nSrcInc) {  nSrcInc_ = nSrcInc; }
  void Set_nSrcMin(double nSrcMin) {  
    nSrcMin_ = nSrcMin; 
    opt_Default_SrcRange_ = false;
  }
  void Set_nSrcMax(double nSrcMax) {  
    nSrcMax_ = nSrcMax; 
    opt_Default_SrcRange_ = false;
  }

  void PrepareAnalysis();

  void MaximizeLlh() {
    if (method_ == SIMPLE) {
      MaximizeLlh_Simple();
    } 
    else if (method_ == MINUIT_MIGRAD) {
      MaximizeLlh_MinuitMigrad();
    }
  }

  double EvaluateLlh(double nSrc);

  virtual double EvaluateLlh(double *parValueArray) {
    return EvaluateLlh(parValueArray[0]);
  }

  double GetTestStatistic() const { return Get_logLambdaBest(); }
  double GetEstProb() const { return Get_oneSided_chiSqProb(); }
  double GetPar(int i) const { 
    if (i==0) { return Get_nSrcBest(); } 
    else { return 0.; }
  }

  double Get_logLambdaBest() const {return logLambdaBest_;}
  double Get_chiSq() const { return chiSq_;}
  double Get_sigma() const { return sqrt(chiSq_);}
  double Get_twoSided_chiSqProb() const {return chiSqProb_;}
  double Get_oneSided_chiSqProb() const {return chiSqProb_/2.;}
  double Get_nEvents() const {return nEvents_;}
  double Get_nSrcBest() const {return nSrcBest_;}
  double Get_weightBest() const {return nSrcBest_/nEvents_;}

};


double llhRatio(const vector<double>& ProbRatio, double WeightS);

#endif //  LLH_LLHFUNCTIONS_H_
