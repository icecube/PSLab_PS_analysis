#ifndef LLH_MULTIANALYSISFN_H_
#define LLH_MULTIANALYSISFN_H_

#include "TGraph.h"
#include "TStopwatch.h" // This may be temporary, for optimization purposes

#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"

#include "rootExt/public/generalfunctions.h"
#include "llh/public/LlhFunctionsBase.h"
#include "llh/public/NewLlhEnergy.h"
#include "llh/public/MinuitAnalysisFn.h"
#include "llh/public/MultiAnalysisSet.h"
#include "llh/public/I3Analysis.h"
#include "llh/public/I3Event.h"

class FluxBase;

class MultiAnalysisFn : public AnalysisFn {

 protected:

  TMinuit* minuit_;
  vector<AnalysisFn*> analysisFnVect_;

  const NewLlhEnergy_ParTranslator* parTrans_;
  vector<MinuitParDef> parDefVect_;

  int status_;

  vector<const EnergyProb*> eProbVect_;
  vector<double> energyRatioVect_; 
  vector<double> spaceRatioVect_;
  vector<double> eventRatioVect_;
  vector<const I3Event*> events_;

  bool useEnergy_;
  int eMaxRatioWarnStatus_;
  bool storeRatios_;

  int icstatWarnLevel_;
  // print warning if, after minimization, icstat<=icstatWarnLevel_.
  //   -1 never prints anything
  //    0 = default (no covar matrix, something wrong)
  //    3 prints something always) 
  
  double nSrcMin_;
  double srcFracMax_;
  double nSrcGuess_;

  double gammaMin_;
  double gammaMax_;
  double gammaGuess_;

  double logLambdaBest_;
  double chiSq_;
  double chiSqProb_;
  double nSrcBest_;
  double gammaBest_;

  int minuitOut_;

  int nEventsTot_;

  int monitorLevel_;

  double optimizeTolerance_;
  double optimizeAngleDeg_;

  void OptimizeEventSelection();

 public:

  int nPar_;

  MultiAnalysisFn();

  virtual ~MultiAnalysisFn() {
    if(minuit_) delete minuit_;
  }
 
  void AddAnalysisFn(AnalysisFn* llh) {
    analysisFnVect_.push_back(llh);
  }

  void SetAnalysisSet(AnalysisSet*) {
    log_error("use AddAnalysis  instead of  SetAnalysisSet(aSet)\n");
  }

  void SetSearchCoord(const Coord& coord) {
    srcCoord_ = &coord;
    for (int i=0; i<int(analysisFnVect_.size()); ++i) {
      analysisFnVect_[i]->SetSearchCoord(coord);
    }
  }

  void SetParDefs(vector<MinuitParDef>& parDefVect) {
    parDefVect_ = parDefVect;
    nPar_ = parDefVect_.size();
  }

  void SetParTranslator(const NewLlhEnergy_ParTranslator* pt) { parTrans_ = pt; }

  void PrepareAnalysis(); 
  void MaximizeLlh();

  virtual double EvaluateLlh(double ns, double gamma);
  double EvaluateLlh(double *parValueArray) {
    return EvaluateLlh(parValueArray[0], parValueArray[1]);
  }


  virtual double GetPar(int i) const {
    double par=0., err=0.;
    minuit_->GetParameter(i, par, err);
    return par;
  }

  virtual double EvalFCN(const vector<double>&) const {
    Printf("EvalFCN not implemented. Return 1e50 by default");
    return 1e50;
  }

  void SetUseEnergy(bool use) { useEnergy_ = use; }
  void SetStoreRatios(bool use) { storeRatios_ = use; }
  void SetEMaxRatioWarnOnlyOnce(bool opt) {
    if (opt) { eMaxRatioWarnStatus_ = 0; } // goes to 1 with first warning
    else { eMaxRatioWarnStatus_ = -1; } // stays -1 for all warnings
  }
  
  // level 0: no logging info
  // level 1: event optimization results for each trial
  // level 2: maximum likelihood stuff for each iteration of minimizer
  // level 3: details of selected events
  void SetMonitorLevel(int monitor) { monitorLevel_ = monitor; }

  void SetIcstatWarnLevel(int level) { icstatWarnLevel_ = level; }

  // under optimization, the absolute level of error in logLambda tolerated.
  // if zero, no optimization (i.e. zero error tolerated)
  void SetOptimizeTolerance(double tolerance) { optimizeTolerance_ =tolerance;}
  void SetOptimizeSrcFracMax(double srcFracMax) { srcFracMax_ = srcFracMax; }
  void SetOptimizeAngleDeg(double angleDeg) { optimizeAngleDeg_ = angleDeg; }


  int GetMinimizationResult() {return minuitOut_;}

  double GetParsGuess(double & Guess_nsrc, double & Guess_gamma);
  
  vector<double> GetSpaceRatios() { return spaceRatioVect_;}
  vector<double> GetEnergyRatios() { return energyRatioVect_;}
  vector<double> GetEventRatios() { return eventRatioVect_;}
  vector<const I3Event*> GetEventsVector() { return events_;} 

  bool GetUseEnergy() const { return useEnergy_; }
  double GetOptimizeTolerance() const { return optimizeTolerance_; }
  double GetOptimizeSrcFracMax() const { return srcFracMax_; }
  double GetOptimizeAngleDeg() const { return optimizeAngleDeg_; }


  double GetTestStatistic() const { return Get_logLambdaBest(); }
  double GetEstProb() const { return Get_oneSided_chiSqProb(); }

  double Get_logLambdaBest() const {return logLambdaBest_;}
  double Get_chiSq() const { return chiSq_;}
  double Get_sigma() const { return sqrt(chiSq_);}
  double Get_twoSided_chiSqProb() const {return chiSqProb_;}
  double Get_oneSided_chiSqProb() const {return chiSqProb_/2.;}
  double Get_nEvents() const {return nEventsTot_;}
  double Get_nSrcBest() const {return nSrcBest_;}
  double Get_weightBest() const {return nSrcBest_/nEventsTot_;}
  double Get_gammaBest() const {return gammaBest_;}
  vector<AnalysisFn*> GetAnalysisFnVect() const { return analysisFnVect_; }
  const NewLlhEnergy_ParTranslator* GetParTrans() const { return parTrans_; }

  // For time-testing.  These continue and stop with each execution,
  // and can be accessed directly by the user after many executuions for
  // a summary of time used.
  TStopwatch stopwatch_MaximizeLlh_;
  TStopwatch stopwatch_optimize_;
  TStopwatch stopwatch_minuitMigrad_;
  TStopwatch stopwatch_minuitFCN_;
  // DONT FORGET TO *STOP* THESE WATCHES AFTER THEY ARE CREATED IN 
  // LlhEnergy() CONSTRUCTOR !!!
};




#endif // LLH_MULTIANALYSISFN_H_

