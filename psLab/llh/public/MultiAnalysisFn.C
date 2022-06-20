#include "llh/public/LlhEnergy.h"
#include "llh/public/MultiAnalysisFn.h"
#include "iostream"

#include "TMinuit.h"
#include "TGraph.h"

#include "rootExt/public/FunctionsRoot.h"
#include "rootExt/public/generalfunctions.h"
#include "rootExt/public/randomfunctions.h"
#include "rootExt/public/log_report.h"
#include "rootExt/public/ModDistanceFn.h"

#include "llh/public/BkgSpaceProb.h"
#include "llh/public/EnergyProb.h"
#include "llh/public/I3Event.h"

MultiAnalysisFn *Mllh;

void EvalMinuitFCN(int &npar, double *gin, double &f, double *par, int iflag)
{
  vector<AnalysisFn*> analysisFnVect = Mllh->GetAnalysisFnVect();
  f=0;
  for (int i = 0; i < int(analysisFnVect.size()); ++i) {
    for (int j = 0; j < npar; j++)
      {
        if (par[j] != par[j]){
          log_warn("MultiAnalysisFn: Global parameter is NaN. Returning 1e50 in the minimizer\n");

          return;
        }
      }

    //Create a vector<double> it is needed for the Translator
    vector<double> parVect;
    for (int j = 0; j <  npar; j++){
      parVect.push_back(par[j]);
    }
    vector<double> individualParVect = Mllh->GetParTrans()->Translate(i, parVect);

    double pars[npar];
    for (int j = 0; j < npar; j++){
      pars[j] = individualParVect[j];
    }
    
    f -= analysisFnVect[i]->EvaluateLlh(pars);
  }
  return;
}

MultiAnalysisFn::MultiAnalysisFn() :
  //AnalysisFn(2),  // initialize minuit_ with 2 parameters for minimizing
  useEnergy_(true),
  eMaxRatioWarnStatus_(0),  // default -1 means warn every time, 0 only once, >0 never
  storeRatios_(false),
  icstatWarnLevel_(0),
  nSrcMin_(0.),
  srcFracMax_(0.5), // default is half of total nEvents
  gammaMin_(0.),
  gammaMax_(0.),
  logLambdaBest_(0.),
  chiSq_(0.),
  chiSqProb_(0.),
  nSrcBest_(0.),
  gammaBest_(0.),
  nEventsTot_(0),
  monitorLevel_(0),
  optimizeTolerance_(0.),
  optimizeAngleDeg_(0.)
{
  //parDefVect_[0].name = "nSrc";
  //parDefVect_[1].name = "gamma";
  Mllh = this;
  minuit_ = new TMinuit(2);
  minuit_->SetFCN(EvalMinuitFCN);
  minuit_->SetPrintLevel(-1);
  // These start when created, so stop immediately
  stopwatch_MaximizeLlh_.Stop();
  stopwatch_optimize_.Stop();
  stopwatch_minuitMigrad_.Stop();
  stopwatch_minuitFCN_.Stop();
}

void MultiAnalysisFn::PrepareAnalysis(){
  nEventsTot_ = 0.;
  for (int i=0; i<int(analysisFnVect_.size()); ++i) {
    analysisFnVect_[i]->PrepareAnalysis();
    nEventsTot_ += analysisFnVect_[i]->Get_nEvents();
  }
} 

void MultiAnalysisFn::MaximizeLlh()
{
  stopwatch_MaximizeLlh_.Start(false);  //  false = don't reset

  PrepareAnalysis();

  // DEFINE PARAMETERS
  vector<MinuitParDef> pars;
  double nSrcMax = 0;
  GetParsGuess(nSrcGuess_, gammaGuess_);

  //setting ns
  nSrcMax = srcFracMax_*nEventsTot_;
  nSrcMin_ = 0.;
  pars.push_back( MinuitParDef("nSrc",nSrcGuess_,0.2, nSrcMin_, nSrcMax) );

  //setting gamma
  gammaMin_ = 1.;
  gammaMax_ = 4.; 
  pars.push_back( MinuitParDef("gamma",gammaGuess_,0.5, gammaMin_, gammaMax_) );

  SetParDefs(pars);

  Int_t ierflg = 0;
  for (int i=0; i<int(parDefVect_.size()); ++i) {
    const MinuitParDef& pd = parDefVect_[i];
    minuit_->mnparm(i, pd.name.c_str(), pd.initValue, pd.initStepSize,pd.lowLimit, pd.upLimit,ierflg);
  }

  // MINIMIZE

  // Set the Minuit FCN indirectly, by using a wrapper fcn which 
  // internally will point to 'this' object and its own EvalMinuitFCN
  stopwatch_minuitMigrad_.Start(false);  //  false = don't reset
  //minuit_->Migrad();          // do the minimization
  double arglist[2];
  arglist[0] = 500;
  arglist[1] = 0.1;
  minuit_->mnexcm("MIGRAD", arglist ,2,minuitOut_);
  stopwatch_minuitMigrad_.Stop();
  

  // GET RESULTS 
  if(minuitOut_!=0) {
    Printf("I'm trying to make the fit converge...");
    int itrial = 0;
    double tmpnsrcGuess  = 0.;
    double tmpgammaGuess = 0.;
    //try to change ns value
    while(minuitOut_!=0 && itrial<100) {
      tmpnsrcGuess = random_uniform(0.,2*nSrcGuess_);
      Printf("ns guess = %f\n", tmpnsrcGuess);
      minuit_->mnparm(0, "nSrc", tmpnsrcGuess, 0.1, nSrcMin_, nSrcMax, ierflg);
      minuit_->mnexcm("MIGRAD", arglist ,2,minuitOut_);
      if(minuitOut_ == 0 ) Printf("ns was the problem. Now solved. ns=%f, gamma=%f", GetPar(0), GetPar(1));
      itrial++;
    }
    itrial = 0;
    if(minuitOut_!=0){
      minuit_->mnparm(0, "nSrc", nSrcGuess_, 0.1, nSrcMin_, nSrcMax, ierflg); 
      while(minuitOut_!=0 && itrial<100) {
        tmpgammaGuess = random_uniform(gammaGuess_-0.7, gammaGuess_+0.7);
        if(tmpgammaGuess<gammaMin_)tmpgammaGuess=gammaMin_+0.1;
        if(tmpgammaGuess>gammaMax_)tmpgammaGuess=gammaMax_-0.1;
        minuit_->mnparm(1, "gamma", tmpgammaGuess, 0.1, gammaMin_, gammaMax_, ierflg);
        minuit_->mnexcm("MIGRAD", arglist ,2,minuitOut_);
        if(minuitOut_ == 0 ) Printf("gamma was the problem. Now solved. ns=%f, gamma=%f", GetPar(0), GetPar(1));
        itrial++;
      }
    }
  }
  
  if(minuitOut_==0){
    double fmin, fedm, errdef;
    int nvpar, nparx, icstat;
    minuit_->mnstat(fmin, fedm, errdef, nvpar, nparx, icstat);

    if (icstat <= icstatWarnLevel_) {
      log_warn("icstat=%d : ",icstat);
      if (icstat == 0) { log_warn("covar. matrix not calculated at all.\n"); }
      if (icstat == 1) { log_warn("covar. matrix is approx. only\n"); }
      if (icstat == 2) { log_warn("full covar.matrix, but forced pos-def.\n");}
      if (icstat == 3) { log_warn("full, accurate covariance matrix.\n");}
    }

    logLambdaBest_ = -fmin;
    nSrcBest_ = GetPar(0);
    gammaBest_ = GetPar(1);

    if (logLambdaBest_ < 0.) {
      // we can always do better, since LogLambda(ns=0,gamma) = 0.
      //nSrcBest_ = 0.;
      //gammaBest_ = 0.;  // i.e., no fit makes sense, if nSrc=0.
      logLambdaBest_ = 0.;
    }

    chiSq_ = 2.*logLambdaBest_;
    Printf("nsBest=%f, gammaBest=%f, TS=%f", nSrcBest_, gammaBest_, chiSq_);

    double p_temp;
    if (useEnergy_) {
      chisq_prob(chiSq_, 2., &p_temp, &chiSqProb_);
    }
    else {
      chisq_prob(chiSq_, 1., &p_temp, &chiSqProb_);
    }
  }
  else{Printf("fit problem not solved!");} 
  stopwatch_MaximizeLlh_.Stop();
}

/*
TGraph* MultiAnalysisFn::GetContour(double sigma, int npoints, int pa1, int pa2) {
  minuit_->SetFCN(MinuitWrapperFCN);
  SetMinuitWrapperPtr(this);  // point to *our* EvalMinuitFCN

  minuit_->SetErrorDef(sigma*sigma);
  TGraph* g = dynamic_cast<TGraph*> (minuit_->Contour(npoints, pa1, pa2));

  SetMinuitWrapperPtr(NULL);  // reset pointer
  return g;
}
*/

double MultiAnalysisFn::GetParsGuess(double & Guess_nsrc, double & Guess_gamma) {

  double nsStep = 1;
  double llhMax=-100.;
  double llhTmp;

  for(double ns=0.; ns<=20.; ns += nsStep){
    for(double gamma=1.1; gamma<=3.9; gamma+=0.4){
      llhTmp = EvaluateLlh(ns, gamma);
      //Printf("Searching for pars guess: ns=%f, gamma=%f, Fnc=%.2e", ns, gamma, llhTmp);
      if(llhTmp > llhMax){
        llhMax = llhTmp;
        Guess_nsrc = ns;
        Guess_gamma = gamma;
      }
    }
  }
  return llhMax;
}


double MultiAnalysisFn::EvaluateLlh(double ns, double gamma) {
  assert(nEventsTot_);
  int npar = 2;
  double gin[2];
  double f;
  double par[2];
  int iflag = 0;
  par[0] = ns;
  par[1] = gamma;
  EvalMinuitFCN(npar, gin, f, par, iflag);
  return (-f);  // We want Llh, not -Llh
}

