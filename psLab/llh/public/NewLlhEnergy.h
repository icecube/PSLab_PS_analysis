#ifndef NEWLLH_LLHENERGY_H_
#define NEWLLH_LLHENERGY_H_

#include "TStopwatch.h" // This may be temporary, for optimization purposes

#include "llh/public/classes.h"
#include "llh/public/MinuitAnalysisFn.h"
#include "llh/public/MultiAnalysisSet.h"
#include "llh/public/LlhEnergy.h"

// Forward Declarations (when feasible, more efficient than including headers)
class EnergyProb;
//class NewLlhEnergyFCN; 

class NewLlhEnergy : public AnalysisFn {

 protected:
  TMinuit* minuit_;
  //NewLlhEnergyFCN* fcn_; 

  EventPtrList selectedList_;

  vector<const EnergyProb*> eProbVect_;
  vector<double> spaceRatioVect_;

  vector<double> eventRatioVect_;
  // mutable, because we want to be able to store the result of the 
  // llhEvalulation, and Minuit2::FCNBase::operator() is const

  bool optUseEnergy_;
  int eMaxRatioWarnStatus_;
  bool optStoreRatios_;
  bool optParAuto_[2];

  int icstatWarnLevel_;
  // print warning if, after minimization, icstat<=icstatWarnLevel_.
  //   -1 never prints anything
  //    0 = default (no covar matrix, something wrong)
  //    3 prints something always) 
  
  double nSrcMin_;
  double srcFracMax_;

  double gammaMin_;
  double gammaMax_;

  double logLambdaBest_;
  double estProb_;

  double nSrcBest_;
  double gammaBest_;

  double chiSq_;
  double chiSqProb_;

  int nEventsTot_;

  int monitorLevel_;

  double optimizeTolerance_;
  double optimizeAngleDeg_;

  void OptimizeEventSelection();

  int  minuitOut_;

  void StoreLogLambdaBest();

 public:
  // 'friend' gives this function access to the data
  //friend class NewLlhEnergyFCN;

  NewLlhEnergy();
  virtual ~NewLlhEnergy();

  // not sure if this is a good idea to allow direct access... we'll try it
  TMinuit* Minuit() { return minuit_; }
  //NewLlhEnergyFCN* GetFCN() { return fcn_; }

  void SetUseEnergy(bool use) { optUseEnergy_ = use; }

  void SetStoreRatios(bool store) { optStoreRatios_ = store; }

  void SetEMaxRatioWarnOnlyOnce(bool opt) {
    if (opt) { eMaxRatioWarnStatus_ = 0; } // goes to 1 with first warning
    else { eMaxRatioWarnStatus_ = -1; } // stays -1 for all warnings
  }

  // level 0: no logging info
  // level 1: event optimization results for each trial
  // level 2: maximum likelihood stuff for each iteration of minimizer
  void SetMonitorLevel(int monitor) { monitorLevel_ = monitor; }

  // under optimization, the absolute level of error in logLambda tolerated.
  // if zero, no optimization (i.e. zero error tolerated)
  void SetOptimizeTolerance(double tolerance) { optimizeTolerance_ =tolerance;}
  void SetOptimizeSrcFracMax(double srcFracMax) { srcFracMax_ = srcFracMax; }
  void SetOptimizeAngleDeg(double angleDeg) { optimizeAngleDeg_ = angleDeg; }
  void SetEventRatioVect(vector<double> evtRatios) { eventRatioVect_ = evtRatios; }

  void PrepareAnalysis();

  void MaximizeLlh();

  // pa1=0, pa2=1, means x=nsrc, y=gamma
  //  TGraph* GetContour(double sigma, int npoints, int pa1=0, int pa2=1);

  virtual double EvalFCN(const vector<double>& parVect) const;

  double EvaluateLlh(double nSrc, double gamma);

  virtual double EvaluateLlh(double* parArray) {
    return EvaluateLlh(parArray[0], parArray[1]);
  }

  //  vector<int> GetEventNumbers() { return eventNumber_;}
  vector<double>* GetEventRatios() { return &eventRatioVect_;}

  virtual double GetPar(int i) const { double par=0., err=0.; return minuit_->GetParameter(i, par, err); }
  virtual double Get_logLambdaBest() const { return logLambdaBest_; }
  double GetTestStatistic() const { return Get_logLambdaBest(); }
  double GetEstProb() const { return Get_oneSided_chiSqProb(); }

  bool GetUseEnergy() const { return optUseEnergy_; }
  double GetOptimizeTolerance() const { return optimizeTolerance_; }
  double GetOptimizeSrcFracMax() const { return srcFracMax_; }
  double GetOptimizeAngleDeg() const { return optimizeAngleDeg_; }
  bool   GetOptStoreRatios() const { return optStoreRatios_; }

  int GetMonitorLevel() const { return monitorLevel_; }

  double Get_chiSq() const { return chiSq_;}
  double Get_sigma() const { return sqrt(chiSq_);}
  double Get_twoSided_chiSqProb() const {return chiSqProb_;}
  double Get_oneSided_chiSqProb() const {return chiSqProb_/2.;}
  double Get_nEvents() const {return nEventsTot_;}
  double Get_nSrcBest() const {return nSrcBest_;}
  double Get_weightBest() const {return nSrcBest_/nEventsTot_;}
  double Get_gammaBest() const {return gammaBest_;}
  EventPtrList GetSelectedList() const { return selectedList_; }
  vector<const EnergyProb*> GetEProbVect() const { return eProbVect_; }
  vector<double> GetSpaceRatioVect() const { return spaceRatioVect_; }

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


/*
class NewLlhEnergyFCN : public ROOT::Minuit2::FCNBase {
 private:
  NewLlhEnergy* ptr;
 public:
  // Pure Virtual Fn inherited from ROOT::Minuit2::FCNBase
  // This is related to how "errors" are defined, see Minuit2 for documentation
  virtual double Up() const { return 0.5; }

  // Pure Virtual Fn inherited from ROOT::Minuit2.:FCNBase
  // This is what gets minimized: you have to define your likelihood
  virtual double operator() (const vector<double>& par) const;

  // Here's where the connection is made to FCN can access the data
  virtual void Point(AnalysisFn* fn) {
    ptr = dynamic_cast<NewLlhEnergy*>(fn); 
  }
};
*/


// For MultiAnalysis: Translate Global Par values to Individual Par Values:


class NewLlhEnergy_ParTranslator : public ParTranslator {
 protected:
  vector< vector<double> > srcWeightVect_;
  double gammaMin_;
  double gammaMax_;
  int nStopsGamma_;

 public:
  virtual ~NewLlhEnergy_ParTranslator() { }

  virtual const vector<double>
    Translate(int llhIndex, const vector<double>& parGlobal) const;
  
  virtual void SetRange(double gammaMin, double gammaMax, int nStopsGamma) {
    gammaMin_ = gammaMin;
    gammaMax_ = gammaMax;
    nStopsGamma_ = nStopsGamma;
    srcWeightVect_.clear();
  }
  
  virtual void SetTranslator(MultiAnalysisSet* mASet) {
    SetTranslator(mASet->GetAnalysisSetVect());
  }
  virtual void SetTranslator(const vector<AnalysisSet*>& aSetVect);
};

/*
double LinearInterpolate(double x, 
		   double xmin, double xmax, int nstops, const double* yarray)
{ 
  if (x<xmin) { x=xmin; } 
  if (x>xmax) { x=xmax; }
  

  double ix = (nstops-1) * (x-xmin) / (xmax-xmin);
  
  int i = int(ix);
  
  if (i==nstops-1) { 
      return yarray[nstops-1]; 
      
    } // last element, no interp.
  double y0 = yarray[i];
  
  double y1 = yarray[i+1];
  return y0 + (y1-y0) * (ix-i);
}

double LinearInterpolate(double x, 
		   double xmin, double xmax, const vector<double> yvect) {
  return LinearInterpolate(x, xmin, xmax, yvect.size(), &yvect[0]);
}
*/
NewLlhEnergy *llh;
#endif // NEWLLH_LLHENERGY_H_

