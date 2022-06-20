#ifndef NEWLLH_LLHGAUSTIME_H_
#define NEWLLH_LLHGAUSTIME_H_

#include "TMinuit.h"
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
class NewLlhGausTimeFCN;

// Welcome to LlhGausTime.h!
  
// Have a look in the .C file for more interesting information!
// I really hope that things here are somewhat explanatory.

bool sortinrev(const pair<double,int> &a, const pair<double,int> &b)
{
  return (a.first > b.first);
}

class NewLlhGausTime : public AnalysisFn {

 private:
  TMinuit* minuit4_;
  TMinuit* minuit6_;
  NewLlhGausTimeFCN* fcn_; 
  
  I3Analysis *i3Set_; //AnalysisSet inherits from AnalysisFn
  LocalCoordBkgProb *lcBkgProb_;
  bool useLCBkgProb_;
  bool UseFolded_;
  bool fitSrcCoord_;
   
  vector<const EnergyProb*> eProbVect_;
  vector<double> spaceRatioVect_;
  vector<double> bkgSpaceProbVect_;
  vector<double> lcBkgProbVect_;
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
  vector<int> runID_;
  
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
  bool optParAuto_[6];

  int icstatWarnLevel_;
  // print warning if, after minimization, icstat<=icstatWarnLevel_.
  //   -1 never prints anything
  //    0 = default (no covar matrix, something wrong)
  //    3 prints something always) 
  
  double nSrcMin_;
  double srcFracMax_;

  double gammaMin_;
  double gammaMax_;

  int nFlareGuess_;
  double sigmamin_;
  double sigmamax_;
  bool isSetLowLimSigma_;
  bool isSetUpLimSigma_;
  double TSthr_;

  double logLambdaBest_;
  double estProb_;

  double nSigmaTrunc_;

  vector<pair<double,double> > nSrcBest_;
  vector<pair<double,double> > gammaBest_;
  vector<pair<double,double> > meanBest_;
  vector<pair<double,double> > sigmaBest_;
  vector<pair<double,double> > raSrcBest_;
  vector<pair<double,double> > decSrcBest_;

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
  
  vector<double> nsrcGuess_;
  vector<double> gammaGuess_;
  vector<double> meanGuess_;
  vector<double> sigmaGuess_;

  bool doInitGauss_; //flag to use the class method to initialize parameters
  bool optUsePrior_; //flag to choose between standard/bayesian likelihood
                     //with a prior term modelled with a 2D (a)symmetric Gaussian PDF
  bool setNormFromGRL_;
  bool isGaussTrunc_;
  const Coord* priorCoord_; //coord  of the event used to give a prior term
  double sigmaX_;   //X extension of the area associated to prior term
  double sigmaY_;   //Y extension of the area associated to prior term
  double theta_;    //rotation angle of area associated to prior term
                    //(in degrees, from X axis counterclockwise) 

  int  minuitOut_;   //store minimization result (0=fit converged, <0=fit not converged)

  vector<double> startMissRunVect_;
  vector<double> stopMissRunVect_;

  void ClearBestFitParams() {
    nSrcBest_.clear();
    gammaBest_.clear();
    meanBest_.clear();
    sigmaBest_.clear();
    raSrcBest_.clear();
    decSrcBest_.clear();
  }

  void OptimizeEventSelection();
  void StoreLogLambdaBest(TMinuit *minuit);

  //virtual void EvalMinuitFCN(int &npar, double *gin, double &f, 
  //		     double *par, int iflag);

 public:
  EventPtrList selectedList_;
  // 'friend' gives this function access to the data
  friend class NewLlhGausTimeFCN;

  NewLlhGausTime();
  virtual ~NewLlhGausTime() { 
  }
  EquatorialDeg* srcCoordGuess_;	
  double srcCoordUncRA_;
  double srcCoordUncDEC_;
  double stepSizeSrcCoord_;

  void Set_nevs(int nevs) { nevs_ = nevs; }
  void SetMargValue(double margVal) { margValue_ = margVal; }
  
  void SetOptStoreRatios(bool flag) {optStoreRatios_ = flag;}
  void SetAnalysis(AnalysisSet *aSet, Coord& sourceCoord) {
    aSet_  = aSet;
    i3Set_ = dynamic_cast<I3Analysis*> (aSet);
    srcCoord_ = &sourceCoord;
  }
  AnalysisSet* GetAnalysisSet() {return aSet_;}
  
  void SetLocalCoordBkgProb( LocalCoordBkgProb lcBP , bool folded=false) {
    lcBkgProb_ = lcBP.Clone();
    useLCBkgProb_ = true;
    UseFolded_ = folded;
  }
  void SetUseFitSrc(bool flag) {//to activate the fit of the src position in the llh
    fitSrcCoord_ = flag;
    optParAuto_[4] = flag;
    optParAuto_[5] = flag;
  }
  void SetFitSrc(EquatorialDeg& srcCoord, double uncRA, double uncDEC, 
		 double stepSize) {
    //to activate the fit of the src position in the llh
    //and set their init values, ranges and minuit step
    fitSrcCoord_   = true;
    optParAuto_[4] = true;
    optParAuto_[5] = true;
    srcCoordGuess_ = &srcCoord;
    srcCoordUncRA_    = uncRA;
    srcCoordUncDEC_   = uncDEC;
    stepSizeSrcCoord_ = stepSize;
  }
  void SetParamGuess(vector<double> nsGuess, vector<double> gammaGuess, 
         vector<double> timeGuess, vector<double> sigmaGuess,
		     bool optParAuto) {
    //in case the initialization of the (4) llh parameters is done manually
    doInitGauss_= false;
    nsrcGuess_  = nsGuess;
    gammaGuess_ = gammaGuess;
    meanGuess_  = timeGuess;
    size_t nel = sigmaGuess.size();
    for(size_t i=0; i<nel; i++) sigmaGuess_[i] = pow(10,sigmaGuess[i]);
    //optParAuto_[i] = true  : param. i free to vary in the fit minimization
    //optParAuto_[i] = false : param. i kept fixed in the fit minimization
    optParAuto_[0] = optParAuto;
    optParAuto_[1] = optParAuto;
    optParAuto_[2] = optParAuto;
    optParAuto_[3] = optParAuto;
  }
  void SetParamGuess(vector<double> nsGuess, vector<double> gammaGuess, 
         vector<double> timeGuess, vector<double> sigmaGuess,
		     double srcRAGuess, double srcDECGuess,
		     bool optParAuto) {
    //in case the initialization of the (6) llh parameters is done manually
    doInitGauss_= false;
    nsrcGuess_  = nsGuess;
    gammaGuess_ = gammaGuess;
    meanGuess_  = timeGuess;
    size_t nel = sigmaGuess.size();
    for(size_t i=0; i<nel; i++) sigmaGuess_[i] = pow(10,sigmaGuess[i]);
    //optParAuto_[i] = true  : param. i free to vary in the fit minimization
    //optParAuto_[i] = false : param. i kept fixed in the fit minimization
    optParAuto_[0] = optParAuto;
    optParAuto_[1] = optParAuto;
    optParAuto_[2] = optParAuto;
    optParAuto_[3] = optParAuto;
    SetUseFitSrc(true);
    srcCoordGuess_->SetCoords(srcRAGuess,srcDECGuess);
  }  
  void SetOptParAuto (bool flag) {
    //optParAuto_[i] = true  : param. i free to vary in the fit minimization
    //optParAuto_[i] = false : param. i kept fixed in the fit minimization
    optParAuto_[0] = flag;
    optParAuto_[1] = flag;
    optParAuto_[2] = flag;
    optParAuto_[3] = flag;
  }
  void SetOptParAuto (int i, bool flag) {
    optParAuto_[i] = flag;
  }
  void SetSpatialPrior(const Coord& priorCoord, double sigmaX, double sigmaY, double theta/*in degree*/) {
    optUsePrior_ = true;
    priorCoord_  = &priorCoord;
    sigmaX_ = sigmaX;
    sigmaY_ = sigmaY;
    theta_  = theta;
  }

  void SetLowLimitSigma(double limit) {sigmamin_ = limit; isSetLowLimSigma_ = true;}
  void SetUpLimitSigma(double limit) {sigmamax_ = limit; isSetUpLimSigma_ = true;}
  void SetNormFromGRL(bool flag) {setNormFromGRL_ = flag;}

  void UnsetParamGuess() {
    doInitGauss_= true;
  }
  bool IsFitSrc() {return fitSrcCoord_;}
  
  NewLlhGausTimeFCN* GetFCN() { return fcn_; }

  double tmin_;
  double tmax_;
 
  TH1 * pvalHisto;
  double close_;
  bool JimsTerm_;
  bool SpectralPenalty_;
  bool useBkgPdfFunc_;
  double maxClusterLength_;
  int ndof_;
  double laglimit_;
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

  vector<const EnergyProb*> GetEProbVect() { return eProbVect_;       }
  vector<double>  GetBkgSpaceProbVect()    { return bkgSpaceProbVect_;}
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
  vector<int>* GetEventID()	           { return &eventID_;        }
  vector<int>* GetRunID()             { return &runID_;          }
  bool GetIsGaussTrunc()              { return isGaussTrunc_;    }
  int GetNSigTrunc()                  { return nSigmaTrunc_;     }
 
  LocalCoordBkgProb* GetLcBkgProb() { return lcBkgProb_; }
  EventPtrList GetSelectedList() { return selectedList_; }
  virtual double GetPar(int i) const { 
    double par=0., err=0.;
    minuit4_->GetParameter(i, par, err);
    return par;
  }

  virtual double Get_logLambdaBest() const { return logLambdaBest_; }
  double GetTestStatistic() const { return Get_logLambdaBest(); }
  double GetEstProb() const { return Get_oneSided_chiSqProb(); }

  double GetProbFromHisto(double teststat);//, bool useFit=false);
  void SetNullTestStat(TH1D * inputhisto);

  int GetMonitorLevel() { return monitorLevel_; }
  int GetNEventsSelected() { return nEventsSelected_; }
  bool GetUseEnergy() const { return optUseEnergy_; }
  bool GetUsePrior() const { return optUsePrior_; }
  bool GetOptStoreRatios() { return optStoreRatios_; }
  bool GetUseLCBkgProb() { return useLCBkgProb_; }
  bool GetUseFolded() { return UseFolded_; }
  bool IsSetNormFromGRL() { return setNormFromGRL_; }
  double GetOptimizeTolerance() const { return optimizeTolerance_; }
  double GetOptimizeSrcFracMax() const { return srcFracMax_; }
  double GetOptimizeAngleDeg() const { return optimizeAngleDeg_; }

  double Get_chiSq() const { return chiSq_;}
  double Get_sigma() const { return sqrt(chiSq_);}
  double Get_twoSided_chiSqProb() const {return chiSqProb_;}
  double Get_oneSided_chiSqProb() const {return chiSqProb_/2.;}
  double Get_nEvents() const {return nEventsTot_;}
  vector<double> Get_weightBest() const {
    vector<double> fracVect; 
    size_t nel = nSrcBest_.size();
    for(size_t i=0; i<nel; i++) fracVect.push_back(nSrcBest_[i].first/nEventsTot_); 
    return fracVect;
  }
  vector<pair<double,double> > Get_nSrcBest() const {return nSrcBest_;}
  vector<pair<double,double> > Get_gammaBest() const {return gammaBest_;}
  vector<pair<double,double> > Get_meanBest() const {return meanBest_;}
  vector<pair<double,double> > Get_sigmaBest() const {return sigmaBest_;}  
  vector<pair<double,double> > Get_raSrcBest() const {return raSrcBest_;}  
  vector<pair<double,double> > Get_decSrcBest() const {return decSrcBest_;}  
  double GetSigmaMax() { return sigmamax_;}
  double GetNSigmaTrunc() { return nSigmaTrunc_; }

  double Get_margWeight() const {return margValue_;}

  int Get_nevs() {return nevs_;}
  int GetNFlareGuess() {return nFlareGuess_;}
  void SetNFlareGuess(int nFl) {nFlareGuess_=nFl;}
  void SetTSthr(double TSthr) {TSthr_=TSthr;}
  void SetSeedWtMin(double value) {seedWtMin=value;}

  vector<double> GetNsrcGuess(){ return nsrcGuess_;}
  vector<double> GetMeanGuess() { return meanGuess_; }
  vector<double> GetSigmaGuess(){ return sigmaGuess_;}

  double GetSigmaMin() { return sigmamin_;}
  double GetGammaMin() {return gammaMin_;}
  double GetGammaMax() {return gammaMax_;}  

  void SetTimeBounds(TimePdf * tPdf){
    tmin_ = tPdf->GetTmin();
    tmax_ = tPdf->GetTmax();
    //    cout << "Time covered: " << tmin_ << " to: " << tmax_ << endl;
  }
  void SetTimeBounds(double tmin, double tmax){
    tmin_ = tmin;
    tmax_ = tmax;
    //    cout << "Time covered: " << tmin_ << " to: " << tmax_ << endl;
  }
      
  void SetMissingRuns(vector<double> &start, vector<double> &stop){
    startMissRunVect_.clear();
    stopMissRunVect_.clear();
    if(start.size() != stop.size()) log_error("ERROR: start and stop of missing runs have different size\n"); 
    for(unsigned int i=0; i<start.size(); ++i){ 
      startMissRunVect_.push_back(start[i]);
      stopMissRunVect_.push_back(stop[i]);
    } 
  }

  vector<double> GetStartMissRunVect() {return startMissRunVect_;}
  vector<double> GetStopMissRunVect() {return stopMissRunVect_;}
  vector<double> GetProbPairRatios() {return spaceRatioVect_;}
  vector<I3Event> GetEventVector() {return eVect_;}
     
  void SetStatWarnLevel(int warnlevel) { icstatWarnLevel_ = warnlevel; }
  void SetLivetime(double t) { livetime_ = t; }
  double GetLivetime() { return livetime_; }
  double calcSpatialPrior(double testRA, double testDEC);

  void GetFlareGuess(bool UseE, int UseSpace, vector<double> & Guess_nsrc, vector<double> & Guess_gamma, vector<double> & Guess_mean, vector<double> & rms);

  int GetMinimizationResult() {return minuitOut_;}

  void SetNSigmaTrunc(double nsig) {
    nSigmaTrunc_ = nsig;
    isGaussTrunc_ = true;
  }

  // For time-testing.  These continue and stop with each execution,
  // and can be accessed directly by the user after many executuions for
  // a summary of time used.
  TStopwatch stopwatch_MaximizeLlh_;
  TStopwatch stopwatch_optimize_;
  TStopwatch stopwatch_minuitMigrad_;
  TStopwatch stopwatch_minuitFCN_;
  // DONT FORGET TO *STOP* THESE WATCHES AFTER THEY ARE CREATED IN 
  // NewLlhGausTime() CONSTRUCTOR !!!

};


// For MultiAnalysis: Translate Global Par values to Individual Par Values:
class NewLlhGausTime_ParTranslator : public ParTranslator {
 protected:
  vector< vector<double> > srcWeightVect_;
  vector< vector<double> > dataTimeVect_; // a pair of doubles for the
  // start and end times of each dataset
                                          
  vector<double> localweight_; // time weights mean a new normalization for the
  // ns weight which needs to be recalculated each time
  // so this holds onto the information

  double gammaMin_;
  double gammaMax_;
  int nStopsGamma_;

 public:
  virtual ~NewLlhGausTime_ParTranslator() { }

  virtual const vector<double> Translate(int llhIndex, const vector<double>& parGlobal) const;
  virtual const vector<double> Translate(int llhIndex, const vector<double>& parGlobal, vector<vector<double> > weights) const;

  virtual void SetRange(double gammaMin, double gammaMax, int nStopsGamma) {
    gammaMin_ = gammaMin;
    gammaMax_ = gammaMax;
    nStopsGamma_ = nStopsGamma;
    srcWeightVect_.clear();
  }

  virtual vector<vector<double> > MakeWeights(vector<NewLlhGausTime*> llhVect, const vector<double>& parGlobal, int nFlares) const;

  virtual void SetTranslator(MultiAnalysisSet* mASet) {    SetTranslator(mASet->GetAnalysisSetVect());  }
  virtual void SetTranslator(const vector<AnalysisSet*>& aSetVect);
};

NewLlhGausTime *tdLlh;
#endif 
