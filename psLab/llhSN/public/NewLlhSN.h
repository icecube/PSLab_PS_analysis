#ifndef NEWLLH_LLHSN_H_
#define NEWLLH_LLHSN_H_

#include "TFitterMinuit.h"
#include "TStopwatch.h"

#include "llh/public/classes.h"
#include "llh/public/MinuitAnalysisFn.h"

#include "llhSN/public/SNEvent.h"
#include "llhSN/public/SNAnalysis.h"

// Forward Declarations (when feasible, more efficient than including headers)
//class SimpleProb;
class NewLlhSNFCN;


class NewLlhSN : public AnalysisFn {

 private:
  TFitterMinuit* minuit_;
  NewLlhSNFCN* fcn_;
  //SNAnalysis* sna_;

  EventPtrList selectedList_;
  
  vector<const SimpleProb*> sProbVect_;
  vector<double> powRatioVect_;

  int nEventsTot_;
  int nEvents_;

  vector<SNEvent> eVect_;
  string lcFile_;

  bool optParAuto_[3];

  int icstatWarnLevel_;
  // print warning if, after minimization, icstat<=icstatWarnLevel_.
  //   -1 never prints anything
  //    0 = default (no covar matrix, something wrong)
  //    3 prints something always) 
  
  double nSrcMin_;
  double srcFracMax_;

  double logLambdaBest_;
  double estProb_;
  double livetime_;

  double nSrcBest_;
  double lagBest_;
  double threshBest_;

  double chiSq_;
  double chiSqProb_;

  int monitorLevel_;

  void StoreLogLambdaBest();

 public:
  // 'friend' gives this function access to the data
  friend class NewLlhSNFCN;

  double ndof_;

  NewLlhSN();
  virtual ~NewLlhSN() { 
    if (minuit_) delete minuit_;
  }

  // not sure if this is a good idea to allow direct access... we'll try it
  TFitterMinuit* Minuit() { return minuit_; }
  NewLlhSNFCN* GetFCN() { return fcn_; }

  // level 0: no logging info
  // level 1: event optimization results for each trial
  // level 2: maximum likelihood stuff for each iteration of minimizer
  void SetMonitorLevel(int monitor) { monitorLevel_ = monitor; }

  // under optimization, the absolute level of error in logLambda tolerated.
  // if zero, no optimization (i.e. zero error tolerated)
//  void SetOptimizeTolerance(double tolerance) { optimizeTolerance_ =tolerance;}
//  void SetOptimizeSrcFracMax(double srcFracMax) { srcFracMax_ = srcFracMax; }
//  void SetOptimizeAngleDeg(double angleDeg) { optimizeAngleDeg_ = angleDeg; }

  void PrepareAnalysis();

  void MaximizeLlh();

  // pa1=0, pa2=1, means x=nsrc, y=gamma
  //  TGraph* GetContour(double sigma, int npoints, int pa1=0, int pa2=1);

  double EvalFCN(const vector<double>& parVect) const;

  double EvaluateLlh(double nSrc, double lag, double threshold);

  double EvaluateLlh(double* parArray) {
    return EvaluateLlh(parArray[0], parArray[1], parArray[2]);
  }

  //  vector<int> GetEventNumbers() { return eventNumber_;}
  vector<double> GetEventRatios() { return powRatioVect_;}
  vector<SNEvent> GetEvents() { return eVect_;}
  
  void SetBlocks(string infile) { lcFile_ = infile; }
  void SetOptParAuto(int i, int b=0) { optParAuto_[i] = b; }

  virtual double GetPar(int i) const { return minuit_->GetParameter(i); }
  virtual double Get_logLambdaBest() const { return logLambdaBest_; }
  double GetTestStatistic() const { return Get_logLambdaBest(); }
  double GetEstProb() const { return Get_oneSided_chiSqProb(); }

//  bool GetUseEnergy() const { return optUseEnergy_; }
//  double GetOptimizeTolerance() const { return optimizeTolerance_; }
//  double GetOptimizeSrcFracMax() const { return srcFracMax_; }
//  double GetOptimizeAngleDeg() const { return optimizeAngleDeg_; }

  double Get_chiSq() const { return chiSq_;}
  double Get_sigma() const { return sqrt(chiSq_);}
  double Get_twoSided_chiSqProb() const {return chiSqProb_;}
  double Get_oneSided_chiSqProb() const {return chiSqProb_/2.;}
  double Get_nEvents() const {return nEventsTot_;}
  double Get_nSrcBest() const {return nSrcBest_;}
  double Get_weightBest() const {return nSrcBest_/nEventsTot_;}
//  double Get_gammaBest() const {return gammaBest_;}

  void SetLivetime(double t) { livetime_ = t; }
  double GetLivetime() { return livetime_; }
  void ScanParams(double & nsGuess, double & lagGuess);

  //void SetAnalysis(SNAnalysis ss) { sna_ = &ss; }

  // For time-testing.  These continue and stop with each execution,
  // and can be accessed directly by the user after many executuions for
  // a summary of time used.
  TStopwatch stopwatch_MaximizeLlh_;
  TStopwatch stopwatch_optimize_;
  TStopwatch stopwatch_minuitMigrad_;
  TStopwatch stopwatch_minuitFCN_;
  // DONT FORGET TO *STOP* THESE WATCHES AFTER THEY ARE CREATED IN 
  // NewLlhEnergy() CONSTRUCTOR !!!

};


class NewLlhSNFCN : public ROOT::Minuit2::FCNBase {
 private:
  NewLlhSN* ptr;
 public:
  // Pure Virtual Fn inherited from ROOT::Minuit2::FCNBase
  // This is related to how "errors" are defined, see Minuit2 for documentation
  virtual double Up() const { return 0.5; }

  // Pure Virtual Fn inherited from ROOT::Minuit2.:FCNBase
  // This is what gets minimized: you have to define your likelihood
  virtual double operator() (const vector<double>& par) const;

  // Here's where the connection is made to FCN can access the data
  virtual void Point(AnalysisFn* fn) {
    ptr = dynamic_cast<NewLlhSN*>(fn); 
  }
};


#endif // NEWLLH_LLHSN_H_


