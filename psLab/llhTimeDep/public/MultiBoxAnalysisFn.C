#include "llhTimeDep/public/MultiBoxAnalysisFn.h"
#include "rootExt/public/FunctionsRoot.h"
#include "rootExt/public/randomfunctions.h"
MultiBoxAnalysisFn *tdMBoxLlh;

void llhFncMultiBox(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  f = 0.;
  double rr;
  vector<double> parVect;
  for(int i=0; i<npar; i++) parVect.push_back(par[i]);
  //cout << "parVect " << parVect[0] << " " << parVect[1] << " " << parVect[2] << " " << parVect[3] << endl;
  double boxmin = parVect[2];
  double boxmax = parVect[3];

  for (int i=0; i<int(tdMBoxLlh->analysisFnVect_.size()); ++i) {
    vector<double> individualParVect = tdMBoxLlh->parTrans_->Translate(i, parVect);
    //cout << "individualParVect " << individualParVect[0] << " " << individualParVect[1] << " " << individualParVect[2] << " " << individualParVect[3] << endl;
    tdBLlh = tdMBoxLlh->analysisFnVect_[i];
    rr=tdMBoxLlh->analysisFnVect_[i]->EvalFCN(individualParVect);
    f += rr;
    //Now we correct for the marginalization. 
    //Since it was applied N times on each individual llh
    //we have to remove it all the times. 
    //Since result = -LogLikelihood it needs to be added (+). 
    //In the end a single marginalization term will be added,
    //containing the integral of the Gaussian time PDF
    //across the full period of the analysis
  
    if (tdMBoxLlh->analysisFnVect_[i]->JimsTerm_ && tdMBoxLlh->analysisFnVect_[i]->Get_margWeight() < 0.) { 
      f += tdMBoxLlh->analysisFnVect_[i]->Get_margWeight();
    }
  }

  //f += par[0];

  //add marginalization term
  double margValue = log( (fabs(boxmax-boxmin)) / ( tdMBoxLlh->tmax_ - tdMBoxLlh->tmin_ ) );

  if ( tdMBoxLlh->JimsTerm_ && margValue < 0. ) {
    f -=  margValue;
  }
  if(f!=f) f=1e50;
  return;
}

MultiBoxAnalysisFn::MultiBoxAnalysisFn():
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
  tdMBoxLlh = this;
    
  nPar = 4;
  histoForProb_ = false;
  seedWtMin=1000;
  optParAuto_[0] = true;
  optParAuto_[1] = true;
  optParAuto_[2] = true;
  optParAuto_[3] = true;
}

void MultiBoxAnalysisFn::StoreLogLambdaBest(TMinuit *minuit) {
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


vector<I3Event> MultiBoxAnalysisFn::GetAllEvents() {
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


void MultiBoxAnalysisFn::MaximizeLlh() { 
  PrepareAnalysis();
  parDefVect_.clear();

  GetFlareGuess(optUseEnergy_, nsrcGuess_, gammaGuess_, boxMinGuess_, boxMaxGuess_);

  TMinuit *minuit = new TMinuit(4);
  Double_t arglist[1];
  arglist[0] = 0.5; //0.5 for likelihood fit
  Int_t ierflg=0;
  minuit->SetFCN(llhFncMultiBox);
  minuit->mnexcm("SET ERR",arglist,1,ierflg);
  minuit->SetPrintLevel(-1);

  if(gammaGuess_<1.) gammaGuess_ = 1.1;
  if(gammaGuess_>4.) gammaGuess_ = 3.9;
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
  //if (!optParAuto_[1]) minuit->FixParameter(1+i*4);//minuit->mnfixp(1+i*4,ierflg);//minuit->mnexcm("FIX ", arglist ,1, ierflg);//FixParameter(1+i*4); 

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
  minuit->GetParameter(2, par, err);
  boxMinBest_  = par;
  minuit->GetParameter(3, par, err);
  boxMaxBest_ =  par;

  delete minuit;
}

void MultiBoxAnalysisFn::GetFlareGuess(bool useE, double & Guess_nsrc, double & Guess_gamma, double & Guess_boxmin, double & Guess_boxmax){

  vector<double> tVectorclose;
  double llhMax= -100.;
  double llhtemp;
  vector<I3Event> evs;

  for (int i=0; i<int(analysisFnVect_.size()); ++i) {
      vector<I3Event> eVect=analysisFnVect_[i]->GetEventVector();
      evs.insert(evs.end(), eVect.begin(), eVect.end());
  }

  double sProb, bProb, spaceRatio, eMaxRatio;
  for (int j=0; j<int(evs.size()); j++) {
    sProb = evs[j].ProbFrom(*srcCoord_);
    bProb = evs[j].GetBkgSpaceProbFn()->GetBkgProbDensity(evs[j]);
    spaceRatio = sProb/bProb;
    eMaxRatio = evs[j].GetEnergyProbFn()->GetEnergyMaxRatio(evs[j]);
    if ( (eMaxRatio*spaceRatio) > seedWtMin ) {
        tVectorclose.push_back(evs[j].GetMJD());
    }
  }

  sort(tVectorclose.begin(),tVectorclose.end());

  if (monitorLevel_ > 0) { cout << tVectorclose.size() << " events with S/B > " << seedWtMin << "." << endl; }

//  const int ngamma=2;
//  double g[ngamma] = {2.0, 3.5};
  const int ngamma=2;
  double g[ngamma] = {2.0, 3.0};
  double boxmin, boxmax;
  for (unsigned int i=0; i<(tVectorclose.size()); i++) {
    boxmin = tVectorclose[i];
    for (unsigned int j=(i+1); j<(tVectorclose.size()); j++) {
      if ( (tVectorclose[j]-boxmin) < maxClusterLength ) {
        boxmax = tVectorclose[j];
      } else {
        continue;
      }
      if (useE) {
        for (int h=0; h<ngamma; h++) {
          llhtemp = EvaluateLlh( j-i, g[h], boxmin, boxmax);
          if (llhtemp > llhMax) {
            llhMax = llhtemp;
            Guess_boxmax = boxmax;
            Guess_boxmin = boxmin;
            Guess_nsrc = j-i;
          }
        }
      } else { // I have no idea why you'd want to do the analysis without energy, but here you go...
        llhtemp = EvaluateLlh( j-i, 0., boxmin, boxmax );
        if (llhtemp > llhMax) {
            llhMax = llhtemp;
            Guess_boxmax = boxmax;
            Guess_boxmin = boxmin;
            Guess_nsrc = j-i;
        }
      }
    }
  }

  tVectorclose.clear();

  double sllhmax = -100;
  double sllhtemp;
  for (double d=1.; d<4.; d+=0.2) { // loop over gamma with 15 steps for best seed
    sllhtemp = EvaluateLlh( Guess_nsrc, d, Guess_boxmin, Guess_boxmax );
    if (sllhtemp > sllhmax ) {
      Guess_gamma = d;
      sllhmax = sllhtemp;
    }
  }
  //if (useE) { } //nip compiler complaints

  //Printf("finished GetFlareGuess: guess boxmin = %f", Guess_boxmin);
}

    
double MultiBoxAnalysisFn::EvaluateLlh(double *parValueArray) {
    vector<double> parVect;
    for (int i=0; i<nPar; ++i) {
        parVect.push_back(parValueArray[i]);
    }
    double minusLlh = EvalFCN(parVect);
    return -minusLlh;   // that is, max llh = - (minimizer result)
}

double MultiBoxAnalysisFn::EvalFCN(const vector<double>& parVect) const {
  double f;
  const int npar = parVect.size();
  double gin[npar];
  double par[npar];
  int iflag = 0;
  for(int i=0; i<npar; i++) par[i] = parVect[i];
  int npars=npar;
  llhFncMultiBox(npars, gin, f, par, iflag);

  return f;
}

double MultiBoxAnalysisFn::GetProbFromHisto(double teststat) const { 
    int bin = pvalHisto_->FindBin(teststat);
    double ptemp = pvalHisto_->GetBinContent(bin);
    return ptemp;
}

void MultiBoxAnalysisFn::SetNullTestStat(TH1D * inputhisto) {
  
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
