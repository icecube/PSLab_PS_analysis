#ifndef NEWLLH_LLHBLOCKTIME_H_
#define NEWLLH_LLHBLOCKTIME_H_

#include "TStopwatch.h" // This may be temporary, for optimization purposes

#include "llh/public/Interpolate.h"
#include "llh/public/classes.h"
#include "llh/public/MinuitAnalysisFn.h"
#include "llh/public/MultiAnalysisSet.h"

#include "llh/public/I3Analysis.h"
#include "llhTimeDep/public/BlockLevel.h"
//#include "TF1.h"

// Forward Declarations (when feasible, more efficient than including headers)
class EnergyProb;

  // Welcome to NewLlhBlockTime.h!
  
  // Have a look in the .C file for more interesting information!
  // I really hope that things here are somewhat explanatory.


class NewLlhBlockTime : public AnalysisFn {

 private:
  TMinuit* minuit_;

  I3Analysis *i3Set_; //AnalysisSet inherits from AnalysisFn
  EventPtrList selectedList_;

  vector<const EnergyProb*> eProbVect_;
  vector<double> spaceRatioVect_;
  
  TH1 * nullTestStat_;
  TH1 * pvalHisto_;
  //TF1 * fitfn_;
  bool histoForProb_;
  
  vector<double> eRatio_;
  vector<I3Event> eVect_;
  vector<I3Event> srcVect_;
  vector<double> eventRatio_;

  mutable vector<double> eventRatioVect_;

  bool optUseEnergy_;
  int eMaxRatioWarnStatus_;
  bool optStoreRatios_;
  bool optParAuto_[4];

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
  double meanBest_;
  double sigmaBest_;

  double chiSq_;
  double chiSqProb_;
  
  int nEventsSelected_;
  int nEventsTot_;
  
  int monitorLevel_;
  double optimizeTolerance_;
  double optimizeAngleDeg_;

  double margValue_;
  int timePdfType_;
  
  double livetime_;
  double sigmamin_;
  int nevs_; //number of events which contribute S/B > 1 in the minimizer
  
  double nsrcGuess_;
  double gammaGuess_;
  double meanGuess_;
  double sigmaGuess_;

  void OptimizeEventSelection();
  void StoreLogLambdaBest();

  int  minuitOut_;

  //virtual void EvalMinuitFCN(int &npar, double *gin, double &f, 
  //		     double *par, int iflag);

 public: 

  NewLlhBlockTime();
  virtual ~NewLlhBlockTime() { 
    if (minuit_) delete minuit_;
  }

  void SetAnalysis(AnalysisSet *aSet, Coord& sourceCoord) {
    aSet_  = aSet;
    i3Set_ = dynamic_cast<I3Analysis*> (aSet);
    srcCoord_ = &sourceCoord;
  }

  // not sure if this is a good idea to allow direct access... we'll try it
  //NewLlhBlockTimeFCN* GetFCN() { return fcn_; }

  double tmin_;
  double tmax_;
  
  string BlocksTimeFile_;
  double blocksThreshold_;
  double blocksMinThreshold_;
  TH1 * pvalHisto;
  double close_;
  bool JimsTerm_;
  bool SpectralPenalty_;
  bool useBkgPdfFunc_;
  double maxClusterLength_;
  int ndof_;
  double laglimit_;
  double fitfrac_;

  void SetUseEnergy(bool use) { optUseEnergy_ = use; }

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

  void PrepareAnalysis();

  void MaximizeLlh();

  // pa1=0, pa2=1, means x=nsrc, y=gamma
  //  TGraph* GetContour(double sigma, int npoints, int pa1=0, int pa2=1);

  virtual double EvalFCN(const vector<double>& parVect) const;

  double EvaluateLlh(double ns, double gamma, double mean, double sigma);
  
  virtual double EvaluateLlh(double *parValueArray) {
    return EvaluateLlh(parValueArray[0], 
                       parValueArray[1], 
                       parValueArray[2], 
                       parValueArray[3]);
  }

  //  vector<int> GetEventNumbers() { return eventNumber_;}
  vector<double> GetEventRatios() { return eventRatioVect_;}

  virtual double GetPar(int i) const { double par=0., err=0.; return minuit_->GetParameter(i, par, err); }
  virtual double Get_logLambdaBest() const { return logLambdaBest_; }
  double GetTestStatistic() const { return Get_logLambdaBest(); }
  double GetEstProb() const { return Get_oneSided_chiSqProb(); }

  double GetProbFromHisto(double teststat);//, bool useFit=false);
  void SetNullTestStat(TH1D * inputhisto);

  bool GetUseEnergy() const { return optUseEnergy_; }
  double GetOptimizeTolerance() const { return optimizeTolerance_; }
  double GetOptimizeSrcFracMax() const { return srcFracMax_; }
  double GetOptimizeAngleDeg() const { return optimizeAngleDeg_; }
  bool   GetOptStoreRatios() const { return optStoreRatios_; }

  bool GetSpectralPenalty() const { return SpectralPenalty_;}
  EventPtrList GetSelectedList() const { return selectedList_; }
  vector<const EnergyProb*> GetEProbVect() const { return eProbVect_; }
  vector<double> GetSpaceRatioVect() const { return spaceRatioVect_; }

  double Get_chiSq() const { return chiSq_;}
  double Get_sigma() const { return sqrt(chiSq_);}
  double Get_twoSided_chiSqProb() const {return chiSqProb_;}
  double Get_oneSided_chiSqProb() const {return chiSqProb_/2.;}
  double Get_nEvents() const {return nEventsTot_;}
  double Get_nSrcBest() const {return nSrcBest_;}
  double Get_weightBest() const {return nSrcBest_/nEventsTot_;}
  double Get_gammaBest() const {return gammaBest_;}
  double Get_meanBest() const {return meanBest_;}
  double Get_sigmaBest() const {return sigmaBest_;}  
  double Get_margValue() const {return margValue_;}
  int    Get_nEventsSelected() const {return nEventsSelected_;}
  int GetMonitorLevel() const { return monitorLevel_; } 
  int Get_nevs() {return nevs_;}

  double Get_nsrcGuess(){ return nsrcGuess_;}
  double Get_meanGuess() { return meanGuess_; }
  double Get_sigmaGuess(){ return sigmaGuess_;}

  double GetSigmaMin() { return sigmamin_;}
  
  void SetTimeBounds(TimePdf * tPdf){
    tmin_ = tPdf->GetTmin();
    tmax_ = tPdf->GetTmax();
//    cout << "Time covered: " << tmin_ << " to: " << tmax_ << endl;
  }
  
  void SetTimeBounds(double tmin, double tmax) {
    tmin_ = tmin;
    tmax_ = tmax;
  }
      
  vector<double> GetProbPairRatios() {return spaceRatioVect_;}
  vector<I3Event> GetEventVector() {return eVect_;}
  vector<I3Event> GetEventVectorSrc() {return srcVect_;}

  void SetStatWarnLevel(int warnlevel) { icstatWarnLevel_ = warnlevel; }
  void SetLivetime(double t) { livetime_ = t; }
  double GetLivetime() { return livetime_; }
  
  int GetTimePdfType() { return timePdfType_; }
    
  void SetBlocks(string str, double t, double offset=0) {
    BlocksTimeFile_ = str;
    blocksThreshold_ = t;
    if (offset) { }
  }
  void SetMinThreshold(float minth){
      blocksMinThreshold_=minth;
  }
  string GetBlocksFile() { return BlocksTimeFile_; }
  double GetBlocksThreshold() { return blocksThreshold_; }
  double SearchForLag();
  void SearchBlockSpace(double lag, double & initV, double & maxT);

  void Increment_nevs(int n) { nevs_+=n; }

  void SetMargValue(double m) {margValue_ = m;}
  void Set_nevs(int n)        {nevs_ = n;}

  // For time-testing.  These continue and stop with each execution,
  // and can be accessed directly by the user after many executuions for
  // a summary of time used.
  TStopwatch stopwatch_MaximizeLlh_;
  TStopwatch stopwatch_optimize_;
  TStopwatch stopwatch_minuitMigrad_;
  TStopwatch stopwatch_minuitFCN_;
  // DONT FORGET TO *STOP* THESE WATCHES AFTER THEY ARE CREATED IN 
  // NewLlhBlockTime() CONSTRUCTOR !!!

};


// For MultiAnalysis: Translate Global Par values to Individual Par Values:


class NewLlhBlockTime_ParTranslator : public ParTranslator {
 protected:
  vector< vector<double> > srcWeightVect_;
  vector< vector<double> > blockWeightVect_;
  
  vector<double> dataTMinVect_; // a pair of vectors for the
  vector<double> dataTMaxVect_;   // start and end times of each dataset
                                          
  vector<double> localweight_; // time weights mean a new normalization for the
                               // ns weight which needs to be recalculated each time
                               // so this holds onto the information
  double gammaMin_;
  double gammaMax_;
  int nStopsGamma_;

  double threshMin_;
  double threshMax_;
  
  string BlocksTimeFile_;  //Adding in where the blocks can be found.

 public:
  NewLlhBlockTime_ParTranslator() { }
  virtual ~NewLlhBlockTime_ParTranslator() { }
  //void SetUpTranslate(const vector<double>& parGlobal);

  void SetUpBlocks(string s) {
    BlocksTimeFile_ = s;
  }

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
  //cout << ix << " " << flush;
  int i = int(ix);
  //cout << i << " " << flush;
  if (i==nstops-1) { return yarray[nstops-1]; } // last element, no interp.
  //cout << "| " << flush;
  double y0 = yarray[i];
  //cout << y0 << " " << flush;
  double y1 = yarray[i+1];
  //cout << y1 << " " << endl;
  return y0 + (y1-y0) * (ix-i);
}

double LinearInterpolate(double x, 
		   double xmin, double xmax, const vector<double> yvect) {
		   //cout << endl << "Interpolate! " << yvect.size() << " " << flush;
  return LinearInterpolate(x, xmin, xmax, yvect.size(), &yvect[0]);
} */


NewLlhBlockTime *tdLlhBlock;

#endif // NEWLLH_LLHBLOCKTIME_H_
