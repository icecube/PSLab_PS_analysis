#ifndef NEWLLH_LLHBOXTIME_H_
#define NEWLLH_LLHBOXTIME_H_

#include "TStopwatch.h" // This may be temporary, for optimization purposes

#include "llh/public/Interpolate.h"
#include "llh/public/classes.h"
#include "llh/public/MinuitAnalysisFn.h"
#include "llh/public/MultiAnalysisSet.h"

#include "llh/public/I3Analysis.h"
#include "llh/public/TimePdf.h"
#include "llhTimeDep/public/LocalCoordBkgProb.h"
//#include "TF1.h"

// Forward Declarations (when feasible, more efficient than including headers)
class EnergyProb;

  // Welcome to LlhGausTime.h!
  
  // Have a look in the .C file for more interesting information!
  // I really hope that things here are somewhat explanatory.


class NewLlhBoxTime : public AnalysisFn {

 public:
  TMinuit* minuit_;

  I3Analysis *i3Set_; //AnalysisSet inherits from AnalysisFn
  LocalCoordBkgProb *lcBkgProb_;
  bool useLCBkgProb_;

  vector<const EnergyProb*> eProbVect_;
  vector<double> spaceRatioVect_;
  vector<double> sigmaSpaceVect_;
  vector<double> spaceWeightVect_;
  vector<double> enWeightVect_;
  //vector<double> enMaxWeightVect_;
  vector<double> timeWeightVect_;
  vector<double> raVect_;
  vector<double> decVect_;
  vector<double> angErrVect_;
  vector<double> eneVect_;
  vector<double> timeVect_;
  vector<int> eventID_;
  
  TH1 * nullTestStat_;
  TH1 * pvalHisto_;
  //TF1 * fitfn_;
  bool histoForProb_;
  
  vector<double> eRatio_;
  vector<I3Event> eVect_;
  vector<double> eventRatio_;

  mutable vector<double> eventRatioVect_;
  // mutable, because we want to be able to store the result of the 
  // llhEvalulation, and Minuit2::FCNBase::operator() is const

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
  double boxMinBest_;
  double boxMaxBest_;
  
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
  int nevs_; //number of events which contribute S/B > 1 in the minimizer
  
  double nsrcGuess_;
  double gammaGuess_;
  double minGuess_;
  double maxGuess_;

  double gammaFixed_;
  double minFixed_;
  double maxFixed_;

  int minuitOut_;

  void OptimizeEventSelection();
  void StoreLogLambdaBest(TMinuit *minuit);

  //virtual void EvalMinuitFCN(int &npar, double *gin, double &f, 
  //		     double *par, int iflag);

  EventPtrList selectedList_;

  NewLlhBoxTime();
  virtual ~NewLlhBoxTime() { 
    //if (minuit_) delete minuit_;
  }

  void SetAnalysis(AnalysisSet *aSet, Coord& sourceCoord) {
    aSet_  = aSet;
    i3Set_ = dynamic_cast<I3Analysis*> (aSet);
    srcCoord_ = &sourceCoord;
  }
  
  void SetLocalCoordBkgProb( LocalCoordBkgProb lcBP ) {
    lcBkgProb_ = lcBP.Clone();
    useLCBkgProb_ = true;
  }

  // not sure if this is a good idea to allow direct access... we'll try it
  TMinuit* Minuit() { return minuit_; }


  double tmin_;
  double tmax_;

  TH1 * pvalHisto;
  bool JimsTerm_;
  bool SpectralPenalty_;
  bool useBkgPdfFunc_;
  double maxClusterLength;
  int ndof_;
  double fitfrac_;
  double seedWtMin;

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
  void SetOptStoreRatios(bool b) {optStoreRatios_ = b;}

  void PrepareAnalysis();

  void MaximizeLlh();

  // pa1=0, pa2=1, means x=nsrc, y=gamma
  //  TGraph* GetContour(double sigma, int npoints, int pa1=0, int pa2=1);

  virtual double EvalFCN(const vector<double>& parVect) const;

  double EvaluateLlh(double ns, double gamma, double boxmin, double boxmax);
  
  virtual double EvaluateLlh(double *parValueArray) {
    return EvaluateLlh(parValueArray[0], 
                       parValueArray[1], 
                       parValueArray[2], 
                       parValueArray[3]);
  }

  vector<const EnergyProb*> GetEProbVect() { return eProbVect_; }
  vector<double>* GetEventRatios()         { return &eventRatioVect_; }
  vector<double>* GetSpatialWeights()      { return &spaceWeightVect_;}
  vector<double>* GetEnergyWeights()       { return &enWeightVect_;   }
  //vector<double>* GetEnergyMaxWeights()  { return &enMaxWeightVect_;}
  vector<double>* GetTimeWeights()         { return &timeWeightVect_; }
  vector<double>* GetraVect()              { return &raVect_;         }
  vector<double>* GetdecVect()             { return &decVect_;        }
  vector<double>* GetAngErrVect()          { return &angErrVect_;     }
  vector<double>* GeteneVect()             { return &eneVect_;        }
  vector<double>* GettimeVect()            { return &timeVect_;       }
  vector<int>* GetEventID()                { return &eventID_;        }


  virtual double GetPar(int i) const {
    double par=0., err=0.;
    return minuit_->GetParameter(i, par, err);
  } 
  virtual double Get_logLambdaBest() const { return logLambdaBest_; }
  double GetTestStatistic() const { return Get_logLambdaBest(); }
  double GetEstProb() const { return Get_oneSided_chiSqProb(); }

  double GetProbFromHisto(double teststat);//, bool useFit=false);
  void SetNullTestStat(TH1D * inputhisto);

  int GetMonitorLevel() { return monitorLevel_; }
  int GetNEventsSelected() { return nEventsSelected_; }
  bool GetUseEnergy() const { return optUseEnergy_; }
  double GetOptimizeTolerance() const { return optimizeTolerance_; }
  double GetOptimizeSrcFracMax() const { return srcFracMax_; }
  double GetOptimizeAngleDeg() const { return optimizeAngleDeg_; }

  double Get_chiSq() const { return chiSq_;}
  double Get_sigma() const { return sqrt(chiSq_);}
  double Get_twoSided_chiSqProb() const {return chiSqProb_;}
  double Get_oneSided_chiSqProb() const {return chiSqProb_/2.;}
  double Get_nEvents() const {return nEventsTot_;}
  double Get_nSrcBest() const {return nSrcBest_;}
  double Get_weightBest() const {return nSrcBest_/nEventsTot_;}
  double Get_gammaBest() const {return gammaBest_;}
  double Get_BoxMinBest() const {return boxMinBest_;}
  double Get_BoxMaxBest() const {return boxMaxBest_;}
  bool GetOptStoreRatios() { return optStoreRatios_; }
  int Get_nevs() {return nevs_;}

  double GetNsrcGuess(){ return nsrcGuess_;}
  double GetMinGuess() { return minGuess_; }
  double GetMaxGuess() { return maxGuess_; }
  AnalysisSet* GetAnalysisSet() {return aSet_;}
 
  void SetTimeBounds(TimePdf * tPdf){
    tmin_ = tPdf->GetTmin();
    tmax_ = tPdf->GetTmax();
//    cout << "Time covered: " << tmin_ << " to: " << tmax_ << endl;
  }

  void SetTimeBounds(double timemin, double timemax){
    tmin_ = timemin;
    tmax_ = timemax;
//    cout << "Time covered: " << tmin_ << " to: " << tmax_ << endl;
  }
  
  void SetTimeWindow(double bmin, double bmax) {
    // this fixes the min and max of the box if we don't want to fit.
    // it should break or be screwy if you give a min or max outside of
    // the dataset.
    minGuess_ = bmin;
    maxGuess_ = bmax;
    Int_t ierflg=0;
    
    optParAuto_[2] = false;
    optParAuto_[3] = false;
    minuit_->mnparm(2, "boxmin", bmin, 0., bmin, bmin, ierflg);
    minuit_->mnparm(3, "boxmax", bmax, 0., bmax, bmax, ierflg);
    
    cout << "Analyzing data from " << tmin_ << " to " << tmax_ << ". Flare from " << minGuess_ << " to " << maxGuess_ << "." << endl;
  }

  void SetSeedWtMin(double value) {seedWtMin=value;}
      
  vector<double> GetProbPairRatios() {return spaceRatioVect_;}
  vector<I3Event> GetEventVector() {return eVect_;}
     
  void SetStatWarnLevel(int warnlevel) { icstatWarnLevel_ = warnlevel; }
  void SetLivetime(double t) { livetime_ = t; }
  void Set_nevs(int n) { nevs_ = n; }
  void SetMargValue(double margVal) { margValue_ = margVal; }
  
  void SetGammaFixed(double g) { gammaFixed_ = g; optParAuto_[1] = false; }
  void SetBoxMinFixed(double t) { minFixed_ = t; optParAuto_[2] = false; }
  void SetBoxMaxFixed(double t) { maxFixed_ = t; optParAuto_[3] = false; }
  void SetnSrcMin(double n) { nSrcMin_=n;}

  double GetLivetime() { return livetime_; }
  
  void GetFlareGuess(bool UseE, double & Guess_nsrc, double & Guess_gamma, double & Guess_boxmin, double & Guess_boxmax);
  double Get_margWeight() const {return margValue_;}

  // For time-testing.  These continue and stop with each execution,
  // and can be accessed directly by the user after many executuions for
  // a summary of time used.
  TStopwatch stopwatch_MaximizeLlh_;
  TStopwatch stopwatch_optimize_;
  TStopwatch stopwatch_minuitMigrad_;
  TStopwatch stopwatch_minuitFCN_;
  // DONT FORGET TO *STOP* THESE WATCHES AFTER THEY ARE CREATED IN 
  // NewLlhBoxTime() CONSTRUCTOR !!!

};


// For MultiAnalysis: Translate Global Par values to Individual Par Values:


class NewLlhBoxTime_ParTranslator : public ParTranslator {
 protected:
  vector< vector<double> > srcWeightVect_;
  vector< vector<double> > dataTimeVect_; // a pair of doubles for the
                                          // start and end times of each dataset
  vector< vector<double> > sourcesTimeVect_; // a pair of doubles for the
                                             // start and end times of each source
                                             //ONLY FOR STACKING
  vector<double> enhanceFactor_; //source weights, ONLY FOR STACKING
                                          
  vector<double> localweight_; // time weights mean a new normalization for the
                               // ns weight which needs to be recalculated each time
                               // so this holds onto the information

  double gammaMin_;
  double gammaMax_;
  int nStopsGamma_;

 public:
  virtual ~NewLlhBoxTime_ParTranslator() { }
  void SetUpTranslate(const vector<double>& parGlobal);

  virtual const vector<double>
    Translate(int llhIndex, const vector<double>& parGlobal) const;
  const vector<double>
    TranslateStacking(int llhIndex, const vector<double>& parGlobal) const;
  
  virtual void SetRange(double gammaMin, double gammaMax, int nStopsGamma) {
    gammaMin_ = gammaMin;
    gammaMax_ = gammaMax;
    nStopsGamma_ = nStopsGamma;
    srcWeightVect_.clear();
  }

  void SetPeriodsTimeBounds(vector<double> tmin, vector<double> tmax);
  void SetSourceTimeBounds(vector<double> tmin, vector<double> tmax);
  void SetEnhanceFactor(vector<double> weights);

  virtual void SetTranslator(MultiAnalysisSet* mASet) {
    SetTranslator(mASet->GetAnalysisSetVect());
  }
  virtual void SetTranslator(const vector<AnalysisSet*>& aSetVect);

  void SetTranslatorStacking(MultiAnalysisSet* mASet) {
    SetTranslatorStacking(mASet->GetAnalysisSetVect());
  }
  void SetTranslatorStacking(const vector<AnalysisSet*>& aSetVect);
};

NewLlhBoxTime *tdBLlh;
#endif // NEWLLH_LLHBOXTIME_H_

