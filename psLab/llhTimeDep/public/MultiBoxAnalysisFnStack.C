#include "llhTimeDep/public/MultiBoxAnalysisFnStack.h"
#include "rootExt/public/FunctionsRoot.h"
#include "rootExt/public/randomfunctions.h"
MultiBoxAnalysisFnStack *tdMBoxStackLlh;

void llhFncMultiBoxStack(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  f = 0.;
  double rr;
  vector<double> parVect;
  for(int i=0; i<npar; i++) parVect.push_back(par[i]);
  //cout << "parVect: len = " << npar <<", pars = " << parVect[0] << " " << parVect[1] << endl;

  for (int i=0; i<int(tdMBoxStackLlh->analysisFnVect_.size()); ++i) {
    vector<double> individualParVect = tdMBoxStackLlh->parTrans_->TranslateStacking(i, parVect);
    //cout << "individualParVect " << individualParVect[0] << " " << individualParVect[1] << " " << individualParVect[2] << " " << individualParVect[3] << endl;
    tdBoxStack = tdMBoxStackLlh->analysisFnVect_[i];
    rr=tdMBoxStackLlh->analysisFnVect_[i]->EvalFCN(individualParVect);
    f += rr;
    //Now we correct for the marginalization. 
    //Since it was applied N times on each individual llh
    //we have to remove it all the times. 
    //Since result = -LogLikelihood it needs to be added (+). 
    //In the end a single marginalization term will be added,
    //containing the integral of the Gaussian time PDF
    //across the full period of the analysis
  
    if (tdMBoxStackLlh->analysisFnVect_[i]->JimsTerm_ && tdMBoxStackLlh->analysisFnVect_[i]->Get_margWeight() < 0.) { 
      f += tdMBoxStackLlh->analysisFnVect_[i]->Get_margWeight();
    }
  }

  //f += par[0];
  //if ( tdMBoxStackLlh->JimsTerm_ && margValue < 0. ) {
  //  f -=  margValue;
  //}
  if(f!=f) f=1e50;
  //cout << "parVect: len = " << npar <<", pars = " << parVect[0] << " " << parVect[1] << " f = " << f << endl;  
  return;
}

MultiBoxAnalysisFnStack::MultiBoxAnalysisFnStack():
  nEventsTot_(0),
  boxMinGuess_(54500.),
  boxMaxGuess_(54501.),
  gammaMin_(0.),
  gammaMax_(0.),
  maxClusterLength(200.),
  srcFracMax_(0.5), // default is half of total nEvents
  logLambdaBest_(0.),
  nSrcBest_(0.),
  gammaBest_(0.),
  monitorLevel_(0),
  JimsTerm_(true),
  optUseEnergy_(true),
  optStoreRatios_(false)
{
  tdMBoxStackLlh = this;
    
  nPar = 2;
  histoForProb_ = false;
  seedWtMin=1000;
  optParAuto_[0] = true;
  optParAuto_[1] = true;
  optParAuto_[2] = false;
  optParAuto_[3] = false;
}

void MultiBoxAnalysisFnStack::StoreLogLambdaBest(TMinuit *minuit) {
    Double_t amin, edm, errdef;
    Int_t nvpar, nparx, istat;
    minuit->mnstat(amin, edm, errdef, nvpar, nparx,istat);
    logLambdaBest_ = -amin;  // that is, max llh = - (minimizer result) 

    if ( analysisFnVect_[0]->monitorLevel_ > 1 ) { 
        printf("LogLambdaBest=%12.6lg\n",logLambdaBest_ );
    } 

    // The *worst* logLambdaBest should be zero (i.e. null hypothesis).
    // But minimizer will miss exact zero, leading logLambdaBest_ to be slightly
    // negative. Since this will cause probability calculation to choke, we 
    // fix it here
  
    if (logLambdaBest_ < 0.) { logLambdaBest_ = 0.;}
}


vector<I3Event> MultiBoxAnalysisFnStack::GetAllEvents() {
    vector<I3Event> allEvents;
    vector<I3Event> tempVect;

    for (int i=0; i<int(analysisFnVect_.size()); ++i) {
        NewLlhBoxTime* i3an = dynamic_cast<NewLlhBoxTime*>(analysisFnVect_[i]);
    
        if (i3an->Get_nEvents()==0) { i3an->PrepareAnalysis(); }
        tempVect.clear();
        tempVect = i3an->GetEventVector();
    
        for (int j=0; j<int(tempVect.size()); j++) {
        allEvents.push_back(tempVect[j]);
        }
    }
    return allEvents;
}


void MultiBoxAnalysisFnStack::MaximizeLlh() { 
  PrepareAnalysis();
  parDefVect_.clear();

  GetFlareGuess(optUseEnergy_, nsrcGuess_, gammaGuess_);

  TMinuit *minuit = new TMinuit(4);
  Double_t arglist[1];
  arglist[0] = 0.5; //0.5 for likelihood fit
  Int_t ierflg=0;
  minuit->SetFCN(llhFncMultiBoxStack);
  minuit->mnexcm("SET ERR",arglist,1,ierflg);
  minuit->SetPrintLevel(-1);

  if(gammaGuess_<=1.) gammaGuess_ = 1.1;
  if(gammaGuess_>=4.) gammaGuess_ = 3.9;
  if (!nsrcGuess_) { nsrcGuess_=0.5; }
  
  double nsrcMax=nEventsTot_;
    
  for (unsigned int i=0; i<analysisFnVect_.size(); ++i) {
    if(analysisFnVect_[i]->nEventsTot_<nsrcMax) nsrcMax = analysisFnVect_[i]->nEventsTot_;
  }

  ///////  ns ///////                 
  double nSrcMax = srcFracMax_*nsrcMax;
  double nSrcMin = 0;
  double initStep = 0.1;
  if(nsrcGuess_>nSrcMax) { nsrcGuess_ = nSrcMax-0.1; initStep = 1;}
  minuit->mnparm(0, Form("nSrc"), nsrcGuess_, initStep, nSrcMin, nSrcMax, ierflg);
  if (!optParAuto_[0])  minuit->FixParameter(0);

  ///////  gamma ///////
  gammaMin_ = 1.;
  gammaMax_ = 4.;
  if(!optParAuto_[1]) {gammaGuess_ = gammaFixed_; gammaMin_ = gammaFixed_-1e-4; gammaMax_=gammaFixed_+1e-4; }
  minuit->mnparm(1, Form("gamma"), gammaGuess_, 0.1, gammaMin_, gammaMax_, ierflg);

  /*
  ///////  boxmin time ///////
  double  initStepSize = 1.;
  double  initValue = boxMinGuess_;
  double lowLimit = tmin_;
  double upLimit = tmax_;
  if (!optParAuto_[2]) { initValue = boxMinFixed_; lowLimit = boxMinFixed_-1e-6; lowLimit=boxMinFixed_-1e-6; upLimit=boxMinFixed_+1e-6;}
  minuit->mnparm(2, Form("boxMin"), initValue, initStepSize, lowLimit, upLimit, ierflg);

  ///////  boxmax time ///////
  double initValueS = boxMaxGuess_;
  double lowLimitS = tmin_;
  double upLimitS = tmax_;
  double initStepSizeS = 1.;
  if(!optParAuto_[3]) { initValueS = boxMaxFixed_; lowLimitS = boxMaxFixed_-1e-6; upLimitS = boxMaxFixed_+1e-6; }
  minuit->mnparm(3, "boxMax", initValueS, initStepSizeS, lowLimitS, upLimitS, ierflg);   
  */

  double arglist2[2];
  arglist2[0] = 500;
  arglist2[1] = 0.1;
  minuit->mnexcm("MIGRAD", arglist2, 2, minuitOut_);

  // GET RESULTS                                                             
  StoreLogLambdaBest(minuit);   // Set logLambdaBest_
  double par=0., err=0.;
  minuit->GetParameter(0, par, err);
  nSrcBest_  = par;
  minuit->GetParameter(1, par, err);
  gammaBest_ = par;
  //minuit->GetParameter(2, par, err);
  //boxMinBest_  = par;
  //minuit->GetParameter(3, par, err);
  //boxMaxBest_ =  par;

  delete minuit;
}

void MultiBoxAnalysisFnStack::GetFlareGuess(bool useE, double & Guess_nsrc, double & Guess_gamma){
  double nsStep = 1;
  double llhMax=-100.;
  double llhTmp;

  if(useE){
    for(double ns=0.; ns<=20.; ns += nsStep){
      if(optParAuto_[1]){
        for(double gamma=1.1; gamma<=3.9; gamma+=0.4){
          llhTmp = EvaluateLlh(ns, gamma);
          if(llhTmp > llhMax){
            llhMax = llhTmp;
            Guess_nsrc = ns;
            Guess_gamma = gamma;
          }
        }
      }
      else{
        double gamma = gammaFixed_;
        llhTmp = EvaluateLlh(ns, gamma);
        if(llhTmp > llhMax){
          llhMax = llhTmp;
          Guess_nsrc = ns;
          Guess_gamma = gamma;
        }
      }
    }
  }
  else{
    for(double ns=0.; ns<=20.; ns += nsStep){
      llhTmp = EvaluateLlh(ns, 0);
      if(llhTmp > llhMax){
        llhMax = llhTmp;
        Guess_nsrc = ns;
        Guess_gamma = 0;
      }
    }
  }
}

    
double MultiBoxAnalysisFnStack::EvaluateLlh(double *parValueArray) {
    vector<double> parVect;
    for (int i=0; i<nPar; ++i) {
        parVect.push_back(parValueArray[i]);
    }
    double minusLlh = EvalFCN(parVect);
    return -minusLlh;   // that is, max llh = - (minimizer result)
}

double MultiBoxAnalysisFnStack::EvalFCN(const vector<double>& parVect) const {
  double f;
  const int npar = parVect.size();
  double gin[npar];
  double par[npar];
  int iflag = 0;
  for(int i=0; i<npar; i++) par[i] = parVect[i];
  int npars=npar;
  llhFncMultiBoxStack(npars, gin, f, par, iflag);

  return f;
}

double MultiBoxAnalysisFnStack::GetProbFromHisto(double teststat) const { 
    int bin = pvalHisto_->FindBin(teststat);
    double ptemp = pvalHisto_->GetBinContent(bin);
    return ptemp;
}

void MultiBoxAnalysisFnStack::SetNullTestStat(TH1D * inputhisto) {
  
  // This reads in a specific TH1D as the null test statistic
  // to use for p-values instead of using a chisquare distribution.
  // It also used to fit an exponential to the upper tail (default set to top 0.01%).
  
  // It takes the raw test statistic distribution (pdf) as input
  // and then makes the cdf to work with later.
  
    histoForProb_ = true;
  
    pvalHisto_ = new TH1D();
    nullTestStat_ = new TH1D();
  
    char sc[] = "Scale";
    pvalHisto_ = DescendingCumulate(inputhisto,sc);
    nullTestStat_ = (TH1*)inputhisto->Clone("hnew");  
}
