#ifndef LLH_LLHENERGY_H_
#define LLH_LLHENERGY_H_

#include "TGraph.h"
#include "TStopwatch.h" // This may be temporary, for optimization purposes

#include "llh/public/Interpolate.h"
#include "llh/public/classes.h"
#include "llh/public/LlhFunctionsBase.h"
#include "llh/public/I3Event.h"
#include "llh/public/MinuitAnalysisFn.h"
#include "llh/public/MultiAnalysisSet.h"

// Forward Declarations (when feasible, more efficient than including headers)
class EnergyProb;


class LlhEnergy : public AnalysisLlh {

 private:
  EventPtrList selectedList_;

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

  virtual void EvalMinuitFCN(int &npar, double *gin, double &f, 
			     double *par, int iflag);

 public:

  LlhEnergy();

  virtual ~LlhEnergy() { }

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

  void SetAnalysis(AnalysisSet *aSet, Coord& sourceCoord) {
    aSet_  = aSet;
    srcCoord_ = &sourceCoord;
  }

  void PrepareAnalysis(); // for now, this is all done in MaximizeLlh();
  void MaximizeLlh();

  // pa1=0, pa2=1, means x=nsrc, y=gamma
  TGraph* GetContour(double sigma, int npoints, int pa1=0, int pa2=1);

  double EvaluateLlh(double ns, double gamma);

  //  vector<int> GetEventNumbers() { return eventNumber_;}
  vector<double> GetSpaceRatios() { return spaceRatioVect_;}
  vector<double> GetEnergyRatios() { return energyRatioVect_;}
  vector<double> GetEventRatios() { return eventRatioVect_;}
  vector<const I3Event*> GetEventsVector() { return events_;} 

  virtual double EvaluateLlh(double *parValueArray) {
    return EvaluateLlh(parValueArray[0], parValueArray[1]);
  }

  bool GetUseEnergy() const { return useEnergy_; }
  double GetOptimizeTolerance() const { return optimizeTolerance_; }
  double GetOptimizeSrcFracMax() const { return srcFracMax_; }
  double GetOptimizeAngleDeg() const { return optimizeAngleDeg_; }


  double GetTestStatistic() const { return Get_logLambdaBest(); }
  double GetEstProb() const { return Get_oneSided_chiSqProb(); }
  double GetParsGuess(double & Guess_nsrc, double & Guess_gamma);

  double Get_logLambdaBest() const {return logLambdaBest_;}
  double Get_chiSq() const { return chiSq_;}
  double Get_sigma() const { return sqrt(chiSq_);}
  double Get_twoSided_chiSqProb() const {return chiSqProb_;}
  double Get_oneSided_chiSqProb() const {return chiSqProb_/2.;}
  double Get_nEvents() const {return nEventsTot_;}
  double Get_nSrcBest() const {return nSrcBest_;}
  double Get_weightBest() const {return nSrcBest_/nEventsTot_;}
  double Get_gammaBest() const {return gammaBest_;}
  int GetMinimizationResult() {return minuitOut_;}

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


// For MultiAnalysis: Translate Global Par values to Individual Par Values:
class LlhEnergy_ParTranslator : public ParTranslator {
 protected:
  vector< vector<double> > srcWeightVect_;
  double gammaMin_;
  double gammaMax_;
  int nStopsGamma_;

 public:
  virtual ~LlhEnergy_ParTranslator() { }

  virtual const vector<double> Translate(int llhIndex, const vector<double>& parGlobal) const;

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
double LinearInterpolate(double x, double xmin, double xmax, int nstops, const double* yarray)
{
  if (x<xmin) { x=xmin; }
  if (x>xmax) { x=xmax; }


  double ix = (nstops-1) * (x-xmin) / (xmax-xmin);
  //cout << "nstops=" << nstops << "\tx=" << x << "\txmin=" << xmin << "\txmax=" << xmax << "\tix = " << ix << endl;
  int i = int(ix);

  if (i==nstops-1) {
      return yarray[nstops-1];

    } // last element, no interp.
  double y0 = yarray[i];

  double y1 = yarray[i+1];
  return y0 + (y1-y0) * (ix-i);
}

double LinearInterpolate(double x, double xmin, double xmax, const vector<double> yvect) {
  return LinearInterpolate(x, xmin, xmax, yvect.size(), &yvect[0]);
}
*/
#endif // LLH_LLHENERGY_H_
