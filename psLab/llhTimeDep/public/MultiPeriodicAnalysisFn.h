#ifndef LLH_MULTIPERIODICANALYSISFN_H_
#define LLH_MULTIPERIODICANALYSISFN_H_


#include "TMinuit.h"

#include "rootExt/public/generalfunctions.h"
#include "llh/public/LlhFunctionsBase.h"
#include "llh/public/MinuitAnalysisFn.h"
#include "llh/public/MultiAnalysisSet.h"
#include "llh/public/I3Event.h"
#include "llhTimeDep/public/BlockLevel.h"
#include "llhTimeDep/public/NewLlhPeriodicTime.h"

class FluxBase;
void llhFncMP(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

class MultiPeriodicAnalysisFn : public AnalysisFn {
 public:
  TMinuit* minuit_;
  Int_t minuitOut_;

  vector<NewLlhPeriodicTime*> analysisFnVect_; // this needs to be specifically
                                              // for periodic analyses, not an <I3Analysis>
  const NewLlhPeriodicTime_ParTranslator* parTrans_;
  vector<MinuitParDef> parDefVect_;
  double tmin_;
  double tmax_;
  
  double period_; // period and ephemeris for
  double t0_;     // the binary system.

  TH1 * nullTestStat_; // a few parameters for using a custom teststat
  TH1 * pvalHisto_;   // distribution to estimate p-values
  //TF1 * fitfn_;            // <- this seems to be broken in newer versions of ROOT
  bool histoForProb_;

  double logLambdaBest_;
  void StoreLogLambdaBest();

  double nSrcBest_;
  double gammaBest_;
  double meanBest_;
  double sigmaBest_;

  double nsrcGuess_;
  double gammaGuess_;
  double meanGuess_;
  double sigmaGuess_;
  double sigmamin_;
  
  double gammaMin_;
  double gammaMax_;
  double nEventsTot_;

  int nPar;
  double seedWtMin;

  MultiPeriodicAnalysisFn();
  virtual ~MultiPeriodicAnalysisFn() { 
    if (minuit_) delete minuit_;
  }

  TMinuit* Minuit() { return minuit_; }

  virtual void AddAnalysisFn(AnalysisFn* llh) {
//  virtual void AddAnalysisFn(I3Analysis* llh) {
    NewLlhPeriodicTime* llh1 = dynamic_cast<NewLlhPeriodicTime*>(llh);
    analysisFnVect_.push_back(llh1);
  }

  virtual void SetAnalysisSet(AnalysisSet*) {
    log_error("use AddAnalysis  instead of  SetAnalysisSet(aSet)\n");
  }

  virtual void SetSearchCoord(const Coord& coord) {
    srcCoord_ = &coord;
    for (int i=0; i<int(analysisFnVect_.size()); ++i) {
      analysisFnVect_[i]->SetSearchCoord(coord);
    }
  }
  
  void SetTimePeriodic(double p, double t0=0) { //setting system period+ephemeris
    period_ = p;
    t0_ = t0;
    for (int i=0; i<int(analysisFnVect_.size()); ++i) {
      analysisFnVect_[i]->SetTimePeriodic(p,t0);
    }
  }

  virtual void SetParDefs(vector<MinuitParDef>& parDefVect) {
    parDefVect_ = parDefVect;
  }

  virtual void SetParTranslator(const NewLlhPeriodicTime_ParTranslator* pt) { parTrans_ = pt; }

  virtual void PrepareAnalysis() {
    nEventsTot_ = 0.;
    for (int i=0; i<int(analysisFnVect_.size()); ++i) {
      analysisFnVect_[i]->PrepareAnalysis();
      nEventsTot_ += analysisFnVect_[i]->Get_nEvents();
    }
  }
    
  virtual void MaximizeLlh() {
  
    // Wipes out parameters (minuit crashes easily when changing existing pars)
    minuit_->Clear();
    Double_t arglist[1];
    arglist[0] = 0.5; //0.5 for likelihood fit
    Int_t ierflg=0;
    minuit_->SetFCN(llhFncMP);
    minuit_->mnexcm("SET ERR",arglist,1,ierflg);
    minuit_->SetPrintLevel(-1);
 
    PrepareAnalysis();
    parDefVect_.clear();
    
    //gets our guess parameters (all of them), plus minimum width
    GetFlareGuessPeriodic(nsrcGuess_, gammaGuess_, meanGuess_, sigmaGuess_, sigmamin_);
    
    //find the bounds on spectral index  
    gammaMin_ = analysisFnVect_[0]->GetGammaMin();
    gammaMax_ = analysisFnVect_[0]->GetGammaMax();
   
    // Loop over samples and spectral indices to ensure ns for MESE sample does not exceed sample size
    double nsrcMax=nEventsTot_;
    for (unsigned int i=0; i<analysisFnVect_.size(); ++i) {
        double maxw=-1;
        for (unsigned int k=0;k<(parTrans_->srcWeightVect_[i]).size();k++) {
            if (maxw<parTrans_->srcWeightVect_[i][k]) maxw=parTrans_->srcWeightVect_[i][k];
        }
        nsrcMax= min(analysisFnVect_[i]->nEventsTot_/maxw,nsrcMax);
        cout<<"Check weights : "<< nsrcMax<<" "<< nsrcGuess_<<" weight max: "<<maxw<<endl; 
    } 
//    cout<<"Check weights : "<< nsrcMax<<" "<< nsrcGuess_<<" weight max: "<<maxw<<endl; 
    parDefVect_.push_back( MinuitParDef("nSrc",nsrcGuess_,0.1, 0.,0.5*nsrcMax) );
    parDefVect_.push_back( MinuitParDef("gamma",gammaGuess_,0.1, gammaMin_, gammaMax_) );   
    parDefVect_.push_back( MinuitParDef("mean", meanGuess_, sigmaGuess_, -0.5, 1.5) );
    
    // The width is fit in log10 space, following the untriggered flare search.
    // It probably doesn't make too much of a difference here, since it only searches over a few
    // orders of magnitude in width.
    parDefVect_.push_back( MinuitParDef("sigma",log10(sigmaGuess_), 1., log10(sigmamin_), log10(1.) ) );
    
    for (int i=0; i<int(parDefVect_.size()); ++i) {
      const MinuitParDef& pd = parDefVect_[i];
      minuit_->mnparm(i, pd.name.c_str(), pd.initValue, pd.initStepSize, pd.lowLimit, pd.upLimit, ierflg);
    }
    double arglist2[2];
    arglist2[0] = 500;
    arglist2[1] = 0.1;
    minuit_->mnexcm("MIGRAD", arglist2, 2, minuitOut_);
    StoreLogLambdaBest();
   
    double par=0., err=0.;
    nSrcBest_  = minuit_->GetParameter(0, par, err);
    gammaBest_ = minuit_->GetParameter(1, par, err);
    meanBest_  = minuit_->GetParameter(2, par, err); // The mean is allowed to float between (-1,2) to prevent
                                           // it from getting stuck at 0 or 1.
    while(meanBest_<0.) {meanBest_ = meanBest_+1.;} // Then we make sure the mean is between 0 and 1.
    while(meanBest_>1.) {meanBest_ = meanBest_-1.;}
    sigmaBest_ = pow( 10., minuit_->GetParameter(3, par, err) );
    
  }

  double Get_nSrcBest()  { return nSrcBest_;}
  double Get_gammaBest() { return gammaBest_;}
  double Get_meanBest()  { return meanBest_;}
  double Get_sigmaBest() { return sigmaBest_;}

  double GetNsrcGuess() { return nsrcGuess_; }
  double GetMeanGuess() { return meanGuess_; }
  double GetSigmaGuess(){ return sigmaGuess_;}

  virtual double EvalFCN(const vector<double>& parVect) const;

  virtual double EvaluateLlh(double *parValueArray);
  
  double EvaluateLlh( double a, double b, double c, double d) {    
    double dd[] = {a, b, c, log10(d) }; //use log10(sigma) in minimizer
    double Llh = EvaluateLlh(dd);
    return Llh;
  }

  double GetProbFromHisto(double teststat) const;
  void SetNullTestStat(TH1D * inputhisto);

  virtual double GetPar(int i) const { double par=0., err=0.; return minuit_->GetParameter(i, par, err); }
  virtual double Get_logLambdaBest() const { return logLambdaBest_; }
  virtual double GetTestStatistic() const { return Get_logLambdaBest(); }
  virtual double GetSigmaMin() { return sigmamin_;}
  virtual double GetEstProb() const {
    //cout << "aa " << histoForProb_ << endl;
    if (histoForProb_) { return GetProbFromHisto( Get_logLambdaBest() ); }
    double chiSq = 2. * Get_logLambdaBest();
    //    cout << "chiSq: " << chiSq << endl;
    if(chiSq<0) { chiSq=0.; }
    double p_temp, p;
    //int nDoF = minuit_->GetNumberFreeParameters();
    int nDoF = analysisFnVect_[0]->ndof_;
    //    cout << " nDof: " << nDoF << endl;
    chisq_prob(chiSq, nDoF, &p_temp, &p);
    return p / 2.;  // one-sided chi-sq prob
  }
  
  vector<I3Event> GetAllEvents();
  void GetFlareGuessPeriodic(double & Guess_nsrc, double & Guess_gamma, double & Guess_mean, double & Guess_rms, double & sigmamin);
 
  void SetTimeBounds(double tmin, double tmax) {
    tmin_ = tmin;
    tmax_ = tmax;
  }
  
};

MultiPeriodicAnalysisFn* ptr;
#endif // LLH_MultiPeriodicANALYSISFN_H_

