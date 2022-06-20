#include "llhTimeDep/public/MultiGaussAnalysisFn.h"
#include "rootExt/public/FunctionsRoot.h"
#include "rootExt/public/randomfunctions.h"
MultiGaussAnalysisFn *tdMLlh;

void llhFncMS(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  f = 0.;
  double rr;
  int nFlares = tdMLlh->GetNFlareGuess();
  vector<double> parVect;
  for(int i=0; i<npar; i++) parVect.push_back(par[i]);

  vector< vector<double> > weights;
  weights.clear();
  weights = tdMLlh->GetParTranslator()->MakeWeights(tdMLlh->GetAnalysisFn(), parVect, nFlares);

  for (int i=0; i<int(tdMLlh->GetAnalysisFn().size()); ++i) {
    vector<double> individualParVect = tdMLlh->GetParTranslator()->Translate(i, parVect, weights);
    tdLlh = tdMLlh->GetAnalysisFn().at(i);
    tdLlh->SetNFlareGuess(nFlares);
    rr=tdMLlh->GetAnalysisFn().at(i)->EvalFCN(individualParVect);
    f += rr;
    //Now we correct for the marginalization. 
    //Since it was applied N times on each individual llh
    //we have to remove it all the times. 
    //Since result = -LogLikelihood it needs to be added (+). 
    //In the end a single marginalization term will be added,
    //containing the integral of the Gaussian time PDF
    //across the full period of the analysis
  
    if (tdMLlh->GetAnalysisFn().at(i)->JimsTerm_ && tdMLlh->GetAnalysisFn().at(i)->Get_margWeight() < 0.) { 
      f += tdMLlh->GetAnalysisFn().at(i)->Get_margWeight();
    }
  }

  //add marginalization term
  double norm=0.;
  double margValue=log(2);

  for(int iFl=0; iFl<nFlares; iFl++){
    double zlow  = TMath::Erf( (tdMLlh->GetTmin()-par[2+iFl*4])/(TMath::Power(10, par[3+iFl*4])*TMath::Sqrt(2.)) );
    double zhigh = TMath::Erf( (tdMLlh->GetTmax()-par[2+iFl*4])/(TMath::Power(10, par[3+iFl*4])*TMath::Sqrt(2.)) );
    norm  = ( (zhigh-zlow)/2. );
    margValue += log( TMath::Power(10, par[3+iFl*4]) * norm / ( tdMLlh->GetTmax() - tdMLlh->GetTmin() ) );
  }

  if ( tdMLlh->GetJimsTerm() && margValue < 0. ) {
    f -=  margValue;
  }
  if(f!=f) f=1e50;
  return;
}

MultiGaussAnalysisFn::MultiGaussAnalysisFn():
  srcCoordUncRA_(0.),
  srcCoordUncDEC_(0.),
  stepSizeSrcCoord_(0.1),
  TSthr_(0),
  sigmamin_(10000),
  sigmamax_(-999),
  isSetLowLimSigma_(false),
  isSetUpLimSigma_(false),
  srcFracMax_(0.5), // default is half of total nEvents
  nSigmaTrunc_(4),
  fitSrcCoord_(false),
  doInitGauss_(true),
  optUsePrior_(false),
  JimsTerm_(true),
  optStoreRatios_(false),
  isSingleFlare_(true)
{
  tdMLlh = this;
    
  nPar = 4;
  histoForProb_ = false;
  seedWtMin=1000;
  optParAuto_[0] = true;
  optParAuto_[1] = true;
  optParAuto_[2] = true;
  optParAuto_[3] = true;
  if(IsFitSrc()) {
    optParAuto_[4] = true;
    optParAuto_[5] = true;
  }
  else {
    optParAuto_[4] = false;
    optParAuto_[5] = false;
  }
}

void MultiGaussAnalysisFn::StoreLogLambdaBest(TMinuit *minuit) {
    Double_t amin, edm, errdef;
    Int_t nvpar, nparx, istat;
    minuit->mnstat(amin, edm, errdef, nvpar, nparx,istat);
    logLambdaBest_ = -amin;  // that is, max llh = - (minimizer result)

    if(optUsePrior_) {
      double testRa  = 0.;
      double testDec = 0.;
      if(IsFitSrc()) {
	testRa    = GetPar(4);
	testDec   = GetPar(5);
      }
      else {
	const EquatorialDeg *srcEqDeg = ( dynamic_cast<const EquatorialDeg*>(srcCoord_) );
	testRa    = srcEqDeg->GetRa();
	testDec   = srcEqDeg->GetDec();
      }
      logLambdaBest_ += log(calcSpatialPrior(testRa,testDec));
    }

    if ( analysisFnVect_[0]->GetMonitorLevel() > 1 ) { 
        printf("LogLambdaBest=%12.6lg\n",logLambdaBest_ );
    } 

    // The *worst* logLambdaBest should be zero (i.e. null hypothesis).
    // But minimizer will miss exact zero, leading logLambdaBest_ to be slightly
    // negative. Since this will cause probability calculation to choke, we 
    // fix it here
  
    if (logLambdaBest_ < 0.) { logLambdaBest_ = 0.;}
}
double MultiGaussAnalysisFn::calcSpatialPrior(double testRA, double testDEC) {

  const EquatorialDeg *priorEqDeg =  ( dynamic_cast<const EquatorialDeg*>(priorCoord_) );
  if(!priorEqDeg) {
    Printf("ERROR, prior coordinates not set!");
    return -99.;
  }
  theta_ *= TMath::DegToRad();

  //center of the Gaussian                                                                                                                                                                   
  double centerRa  = priorEqDeg->GetRa();
  double centerDec = priorEqDeg->GetDec();

  //rotation of the point                                                                                                                                                                    
  double Ra  = testRA  - centerRa;
  double Dec = testDEC - centerDec;

  double xm = (Ra)*TMath::Cos(theta_) - (Dec)*TMath::Sin(theta_);
  double ym = (Ra)*TMath::Sin(theta_) + (Dec)*TMath::Cos(theta_);
  double u  = (xm/sigmaX_)*(xm/sigmaX_) + (ym/sigmaY_)*(ym/sigmaY_);

  double prior = TMath::Exp(-u/2.);

  return prior;
}

vector<I3Event> MultiGaussAnalysisFn::GetAllEvents() {
    vector<I3Event> allEvents;
    vector<I3Event> tempVect;

    for (int i=0; i<int(analysisFnVect_.size()); ++i) {
        NewLlhGausTime* i3an = dynamic_cast<NewLlhGausTime*>(analysisFnVect_[i]);
    
        if (i3an->Get_nEvents()==0) { i3an->PrepareAnalysis(); }
        tempVect.clear();
        tempVect = i3an->GetEventVector();
    
        for (int j=0; j<int(tempVect.size()); j++) {
        allEvents.push_back(tempVect[j]);
        }
    }
    return allEvents;
}


void MultiGaussAnalysisFn::MaximizeLlh() { 
  PrepareAnalysis();
  parDefVect_.clear();

  nFlareGuess_ = 1;
  if(doInitGauss_) GetFlareGuessGauss(nsrcGuess_, gammaGuess_, meanGuess_, sigmaGuess_, sigmamin_);
  if(isSingleFlare_) nFlareGuess_ = 1;
  else nFlareGuess_ = nsrcGuess_.size();

  TMinuit *minuit = new TMinuit(4*nFlareGuess_);
  Double_t arglist[1];
  arglist[0] = 0.5; //0.5 for likelihood fit
  Int_t ierflg=0;
  minuit->SetFCN(llhFncMS);
  minuit->mnexcm("SET ERR",arglist,1,ierflg);
  minuit->SetPrintLevel(-1);

  for(int i=0; i<nFlareGuess_; i++) {
    if(gammaGuess_[i]<1.) gammaGuess_[i] = 1.1;
    if(gammaGuess_[i]>4.) gammaGuess_[i] = 3.9;
    if (!nsrcGuess_[i]) { nsrcGuess_[i]=0.5; }
    if(!isSetUpLimSigma_) sigmamax_ = (tmax_-tmin_)/2.;
    if(sigmaGuess_[i]>=sigmamax_) sigmaGuess_[i] = 0.9*sigmamax_;
    if(sigmaGuess_[i]<=sigmamin_) sigmaGuess_[i] = 1.1*sigmamin_;
  }
  double nsrcMax=nEventsTot_;
    
  for (unsigned int i=0; i<analysisFnVect_.size(); ++i) {
    if(analysisFnVect_[i]->Get_nEvents()<nsrcMax) nsrcMax = analysisFnVect_[i]->Get_nEvents();
  }     

  ///////  ns ///////                 
  double nSrcMax = srcFracMax_*nsrcMax;
  double nSrcMin = 0;
  double initStep = 0.1;
  for(int i=0; i<nFlareGuess_; i++)  {
    if(nsrcGuess_[i]>nSrcMax) { nsrcGuess_[i] = nSrcMax-0.1; initStep = 1;}
    minuit->mnparm(0+i*4, Form("nSrc_f%d",i), nsrcGuess_[i], initStep, nSrcMin, nSrcMax, ierflg);
    if (!optParAuto_[0])  minuit->FixParameter(0+i*4);
  }

  ///////  gamma ///////
  gammaMin_ = analysisFnVect_[0]->GetGammaMin();
  gammaMax_ = analysisFnVect_[0]->GetGammaMax();
  gammaMin_ = 1.;
  gammaMax_ = 4.;
  for(int i=0; i<nFlareGuess_; i++)  {
    if(!optParAuto_[1]) {gammaGuess_[i] = gammaFixed_; gammaMin_ = gammaFixed_-1e-4; gammaMax_=gammaFixed_+1e-4; }
    minuit->mnparm(1+i*4, Form("gamma_f%d",i), gammaGuess_[i], 0.1, gammaMin_, gammaMax_, ierflg);
    //if (!optParAuto_[1]) minuit->FixParameter(1+i*4);//minuit->mnfixp(1+i*4,ierflg);//minuit->mnexcm("FIX ", arglist ,1, ierflg);//FixParameter(1+i*4); 
  } 

  ///////  mean time ///////
  for(int i=0; i<nFlareGuess_; i++)  {
    double  initStepSize = sigmaGuess_[i];
    double  initValue = meanGuess_[i];
    double lowLimit;
    if( nFlareGuess_>1 && (meanGuess_[i]-nSigmaTrunc_*sigmaGuess_[i] > tmin_) ) lowLimit = meanGuess_[i]-nSigmaTrunc_*sigmaGuess_[i];
    else lowLimit = tmin_;
    double upLimit;
    if( nFlareGuess_>1 && (meanGuess_[i]+nSigmaTrunc_*sigmaGuess_[i] < tmax_) ) upLimit = meanGuess_[i]+nSigmaTrunc_*sigmaGuess_[i];
    else upLimit = tmax_;
    if (!optParAuto_[2]) { initValue = meanFixed_; initStepSize=0.01; lowLimit=meanFixed_-1e-6; upLimit=meanFixed_+1e-6;}
    minuit->mnparm(2+i*4, Form("mean_f%d",i), initValue, initStepSize, lowLimit, upLimit, ierflg);
    //if (!optParAuto_[2]) minuit->FixParameter(2+i*4);//minuit->mnfixp(2+i*4,ierflg); //minuit->mnexcm("FIX ", arglist ,2+i*4, ierflg);
  }

  ///////  sigma time ///////
  // Using the log10 of sigma as the fit parameter to make things more smooth.
  if(!isSetUpLimSigma_) sigmamax_ = (tmax_-tmin_)/2.; 
  double lowLimitS = log10(sigmamin_);
  double upLimitS = log10(sigmamax_); 
  double initStepSizeS = 1.;
  for(int i=0; i<nFlareGuess_; i++)  {
    if(sigmaGuess_[i]>=sigmamax_) sigmaGuess_[i] = sigmamax_-0.1;
    if(sigmaGuess_[i]<=sigmamin_) sigmaGuess_[i] = sigmamin_+0.1;
    double initValueS = log10(sigmaGuess_[i]);
    if (!optParAuto_[3]) { initValueS = log10(sigmaFixed_); initStepSizeS=1e-6; lowLimitS = log10(sigmaFixed_-1e-6); upLimitS = log10(sigmaFixed_+1e-6);}
    minuit->mnparm(3+i*4, Form("sigma_f%d",i), initValueS, initStepSizeS, lowLimitS, upLimitS, ierflg);
    //cout << "flare "<< i << ": par set " << 3<< " " << "sigma"<< " " << initValueS<< " " << initStepSizeS<< " " << lowLimitS<< " " << upLimitS << endl;
    //if (!optParAuto_[3]) minuit->FixParameter(3+i*4); //minuit->mnexcm("FIX ", arglist ,3+i*4, ierflg);;
  }

  /*for(int i=nFlareGuess_-1; i>=0; i--){
    if (!optParAuto_[3]) minuit->mnfixp(3+i*4,ierflg);
    if (!optParAuto_[2]) minuit->mnfixp(2+i*4,ierflg);
    if (!optParAuto_[1]) minuit->mnfixp(1+i*4,ierflg);
    if (!optParAuto_[0]) minuit->mnfixp(0+i*4,ierflg);
  }*/
        
  if (optParAuto_[4] && optParAuto_[5]) {
    // here we fit src position coordinates                                
    double initValueRA  = srcCoordGuess_->GetRa();
    double initValueDEC = srcCoordGuess_->GetDec();
    double lowLimitRA   = initValueRA  - srcCoordUncRA_;
    double lowLimitDEC  = initValueDEC - srcCoordUncDEC_;
    double upLimitRA    = initValueRA  + srcCoordUncRA_;
    double upLimitDEC   = initValueDEC + srcCoordUncDEC_;
    parDefVect_.push_back( MinuitParDef("Ra",  initValueRA,  stepSizeSrcCoord_, lowLimitRA,  upLimitRA) );
    parDefVect_.push_back( MinuitParDef("Dec", initValueDEC, stepSizeSrcCoord_, lowLimitDEC, upLimitDEC) );
    //cout << "par set " << 4<< " " << "Ra"<< " " << initValueRA << " " << stepSizeSrcCoord_ << " " << lowLimitRA << " " << upLimitRA << endl;                  
    //cout << "par set " << 5<< " " << "Dec"<< " " << initValueDEC << " " << stepSizeSrcCoord_<< " " << lowLimitDEC << " " << upLimitDEC << endl;                
  }
  
  for (int i=0; i<int(parDefVect_.size()); ++i) {
    const MinuitParDef& pd = parDefVect_[i];
    minuit->mnparm(i, pd.name.c_str(), pd.initValue, pd.initStepSize,pd.lowLimit, pd.upLimit,ierflg);
    if(!optParAuto_[i]) minuit->FixParameter(i);
    //cout << "par set " << i<< " " << pd.name.c_str()<< " " << pd.initValue<< " " << pd.initStepSize<< " " <<pd.lowLimit<< " " << pd.upLimit << endl;
  }
  
  double arglist2[2];
  arglist2[0] = 500;
  arglist2[1] = 0.1;
  minuit->mnexcm("MIGRAD", arglist2, 2, minuitOut_);
  
  // GET RESULTS                                                             
  if(minuitOut_!=0) {
    Printf("Fit did not converge, trying to change initial parameters");
    int itrial = 0;
    double tmpnsrcGuess  = 0.;
    double tmpgammaGuess = 0.;
    double tmpsigmaGuess = 0.;
    //try to change init ns value                          
    while(minuitOut_!=0 && itrial<100) {
      for(int i=0; i<nFlareGuess_; i++)  {
        tmpnsrcGuess = random_uniform(0.,2*nsrcGuess_[i]);
        minuit->mnparm(0+i*4, Form("nSrc_f%d",i), tmpnsrcGuess, 0.1, nSrcMin, nSrcMax, ierflg);
      }
      minuit->mnexcm("MIGRAD", arglist2 ,2,minuitOut_);
      itrial++;
    }
    itrial = 0;
    if(minuitOut_!=0) {
      if(optParAuto_[1]){
        //reset init ns value to nsrcGuess_            
        for(int i=0; i<nFlareGuess_; i++)  minuit->mnparm(0+i*4, Form("nSrc_f%d",i), nsrcGuess_[i], 0.1, nSrcMin, nSrcMax, ierflg);
        //try to change init value of gamma                             
        while(minuitOut_!=0 && itrial<100) {
          for(int i=0; i<nFlareGuess_; i++)  {
            tmpgammaGuess = random_uniform(gammaGuess_[i]-0.7,gammaGuess_[i]+0.7);
            if(tmpgammaGuess<gammaMin_)tmpgammaGuess=gammaMin_+0.1;
            if(tmpgammaGuess>gammaMax_)tmpgammaGuess=gammaMax_-0.1;
            minuit->mnparm(1+i*4, Form("gamma_f%d",i), tmpgammaGuess, 0.1, gammaMin_, gammaMax_, ierflg); 
          }
          minuit->mnexcm("MIGRAD", arglist2 ,2,minuitOut_); 
          itrial++;
        }
      }
      itrial = 0;
      if(minuitOut_!=0) {
        if(optParAuto_[3]){
          //reset init gamma value to gammaGuess_                          
          for(int i=0; i<nFlareGuess_; i++)  minuit->mnparm(1+i*4, Form("gamma_f%d",i), gammaGuess_[i], 0.1, gammaMin_, gammaMax_, ierflg); 
          //try to change init value of sigma time       
          while(minuitOut_!=0 && itrial<100) {
            for(int i=0; i<nFlareGuess_; i++)  { 
              tmpsigmaGuess = log10(random_uniform(sigmaGuess_[i]/5.,sigmaGuess_[i]*5.));
              if(tmpsigmaGuess<lowLimitS) tmpgammaGuess=lowLimitS+0.1;
              if(tmpsigmaGuess>upLimitS)  tmpgammaGuess=upLimitS-0.1;
              minuit->mnparm(3+i*4, Form("sigma_f%d",i), tmpsigmaGuess, initStepSizeS, lowLimitS, upLimitS, ierflg);
            }
            minuit->mnexcm("MIGRAD", arglist2 ,2,minuitOut_);
            itrial++;
          }
        }
      }
    }
  }
  if(minuitOut_==0) {
    double par=0., err=0.;
    StoreLogLambdaBest(minuit);   // Set logLambdaBest_   
    ClearBestFitParams();
    vector<double> parVect;
    for(int iflare=0; iflare<nFlareGuess_; iflare++) {
      minuit->GetParameter(0+iflare*4,par,err);
      nSrcBest_.push_back(make_pair(par,err));
      parVect.push_back(par);
      minuit->GetParameter(1+iflare*4,par,err);
      gammaBest_.push_back(make_pair(par,err));
      parVect.push_back(par);
      minuit->GetParameter(2+iflare*4,par,err);
      meanBest_.push_back(make_pair(par,err));
      parVect.push_back(par);
      minuit->GetParameter(3+iflare*4,par,err);
      sigmaBest_.push_back(make_pair(pow(10.,par),pow(10.,err)));
      parVect.push_back(par);
    }

    spaceWeightVect_.clear();
    enWeightVect_.clear();
    timeWeightVect_.clear();
    raVect_.clear();
    decVect_.clear();
    angErrVect_.clear();
    timeVect_.clear();
    eneVect_.clear();
    eventID_.clear();
    runID_.clear();

    if(optStoreRatios_){
      for(int i=0; i<int(analysisFnVect_.size()); ++i) analysisFnVect_[i]->SetOptStoreRatios(true);
    }

    double minusLlh = EvalFCN(parVect);

    for(int i=0; i<int(analysisFnVect_.size()); ++i){
      if(analysisFnVect_[i]->GetOptStoreRatios()){
        spaceWeightVect_.insert(spaceWeightVect_.end(), analysisFnVect_[i]->GetSpatialWeights()->begin(), analysisFnVect_[i]->GetSpatialWeights()->end());
        enWeightVect_.insert(enWeightVect_.end(), analysisFnVect_[i]->GetEnergyWeights()->begin(), analysisFnVect_[i]->GetEnergyWeights()->end());
        timeWeightVect_.insert(timeWeightVect_.end(), analysisFnVect_[i]->GetTimeWeights()->begin(), analysisFnVect_[i]->GetTimeWeights()->end());
        raVect_.insert(raVect_.end(), analysisFnVect_[i]->GetraVect()->begin(), analysisFnVect_[i]->GetraVect()->end());
        decVect_.insert(decVect_.end(), analysisFnVect_[i]->GetdecVect()->begin(), analysisFnVect_[i]->GetdecVect()->end());
        timeVect_.insert(timeVect_.end(), analysisFnVect_[i]->GettimeVect()->begin(), analysisFnVect_[i]->GettimeVect()->end());
        eneVect_.insert(eneVect_.end(), analysisFnVect_[i]->GeteneVect()->begin(), analysisFnVect_[i]->GeteneVect()->end());
        eventID_.insert(eventID_.end(), analysisFnVect_[i]->GetEventID()->begin(), analysisFnVect_[i]->GetEventID()->end());
        runID_.insert(runID_.end(), analysisFnVect_[i]->GetRunID()->begin(), analysisFnVect_[i]->GetRunID()->end());
      }
    }

    for(int i=0; i<int(analysisFnVect_.size()); ++i) analysisFnVect_[i]->SetOptStoreRatios(false);
  }

  delete minuit;
}

void MultiGaussAnalysisFn::GetFlareGuessGauss(vector<double> & Guess_nsrc, vector<double> & Guess_gamma, vector<double> & Guess_mean, vector<double> & Guess_rms, double & sigmamin){

    nFlareGuess_ = 1;
    Guess_mean.clear();
    Guess_rms.clear();
    Guess_nsrc.clear();
    Guess_gamma.clear();

    vector<double> tVectortemp;
    vector<double> tVectorclose;
    vector<I3Event> evs;
    vector<double> LCBkgProb; //first it will be set to 1. if no LCBkgProb to be used and else -1 for normal LCBkgProb, -2 for folded
    for (int i=0; i<int(analysisFnVect_.size()); ++i) {
        vector<I3Event> eVect=analysisFnVect_[i]->GetEventVector();
        evs.insert(evs.end(), eVect.begin(), eVect.end());
        for (unsigned int ii=0;ii<eVect.size();ii++){
            if (analysisFnVect_[i]->GetUseLCBkgProb()) {
                if ( analysisFnVect_[i]->GetUseFolded() ) LCBkgProb.push_back(analysisFnVect_[i]->GetLcBkgProb()->BackgroundLCProb_folded(eVect[ii]));
                else LCBkgProb.push_back(analysisFnVect_[i]->GetLcBkgProb()->BackgroundLCProb(eVect[ii]));
            }                                                  
            else LCBkgProb.push_back(1.);
        }
    }
    
    if(!isSetLowLimSigma_){
      vector<double> ttt;
      for (unsigned int k=0;k<evs.size();k++) ttt.push_back( evs[k].GetMJD());
      sort(ttt.begin(),ttt.end());
      double timed;
      sigmamin=1000.;
      for (unsigned int i=1; i<ttt.size(); i++) {//This assumes events are time-ordered.
          timed = ttt[i] - ttt[i-1];
          if (analysisFnVect_[0]->GetMonitorLevel() > 2)  cout << ttt[i] << " "; 
          if ( (fabs(timed) < sigmamin) && (fabs(timed)>0) ) sigmamin = fabs(timed);
      }
      if (analysisFnVect_[0]->GetMonitorLevel() > 2) cout << endl;
      sigmamin /= sqrt(2.);
      if (analysisFnVect_[0]->GetMonitorLevel() > 2) { cout << " sigmamin= " << sigmamin << endl; }
      ttt.clear();
    }
    
    double llhtemp, rms, avgtime;
    double sProb,bProb,eMaxRatio;
   
    for (int j=0; j<int(evs.size()); j++) { 
        
        sProb = evs[j].ProbFrom(*srcCoord_);
        bProb = evs[j].GetBkgSpaceProbFn()->GetBkgProbDensity(evs[j]);
        eMaxRatio = evs[j].GetEnergyProbFn()->GetEnergyMaxRatio(evs[j]);
        bProb = bProb*LCBkgProb[j];
        if ( (eMaxRatio * sProb/bProb) > seedWtMin && evs[j].GetMJD() > tmin_ && evs[j].GetMJD()<tmax_) tVectorclose.push_back(evs[j].GetMJD());           
    }
    /*if (tVectorclose.size() <10) {
        tVectorclose.clear();
        for (int j=0; j<int(evs.size()); j++) tVectorclose.push_back(evs[j].GetMJD());
    }*/

    sort(tVectorclose.begin(),tVectorclose.end());

    if (analysisFnVect_[0]->GetMonitorLevel() > 0) { cout << tVectorclose.size() << " events with S/B > " << seedWtMin << ". Sigmamin is " << sigmamin << endl; }
  
    const int ngamma=2;
    double g[ngamma] = {2.0, 3.0};

    int kmax = 11; // MUST BE <=11

    int n[] =  {1,2,3,4,5,10,15,20,25,30,35,40};
    int ns[] = {1,2,3,4,5, 9, 5, 5, 5, 5, 5, 5};

    int kmin=0;
  
    while ( int(tVectorclose.size()) < n[kmax] ) { kmax--; }
 
    vector<double> tmpNsVec, nsVec; //tmp vector storing mean, rms, ns from tested flare
    vector<double> tmpGammaVec, gammaVec; //tmp vector storing mean, rms, ns from tested flare
    vector<double> tmpT0Vec, t0Vec; //tmp vector storing mean, rms, ns from tested flare
    vector<double> tmpSigVec, sigVec; //tmp vector storing mean, rms, ns from tested flare
    vector<pair<double,int> > tmpTSVec, TSVec; //tmp vector storing mean, rms, ns from tested flare

    vector<double> tminVec;
    vector<double> tmaxVec;
  
    vector<double> Guess_TS; // vector storing TS from tested flare
  
    int ind = 0;
    double llhMax= -100.;
    double meanBckup, rmsBckup, nsrcBckup, gammaBckup;

    //cout << "kmax " << kmax << endl;
    
    for (int k=kmin; k<kmax; k++){ //k==1 for pairs, k==2 for triples, etc.
        for (unsigned int i=0; i<(tVectorclose.size()-n[k]); i++) {
            for (int j=0;j<=n[k];j++){
                tVectortemp.push_back(tVectorclose[i+j]);
            }
            //get mean of the times to seed this gaussian
            double sum=0, sumsq=0;
            for (unsigned int ii=0;ii<tVectortemp.size();ii++) sum += tVectortemp[ii];
            avgtime = sum/tVectortemp.size();

            for (unsigned int ii=0;ii<tVectortemp.size();ii++) sumsq += (tVectortemp[ii] - avgtime)*(tVectortemp[ii] - avgtime);
            //calculate the rms to use as the sigma
            rms = sqrt( sumsq/tVectortemp.size() );

            //cout << " GUESS k " << k << " i " << i << " avgtime " << avgtime << " rms " << rms << endl;

            llhtemp = 0.;

            for (int h=0; h<ngamma; h++) {
                llhtemp = EvaluateLlh( ns[k]+1., g[h], avgtime, rms );
                if (llhtemp > llhMax) {
                  llhMax = llhtemp;
                  meanBckup  = avgtime;
                  rmsBckup   = rms;
                  nsrcBckup  = ns[k]+1.;
                  gammaBckup = g[h];
                }
                if(llhtemp>TSthr_) {
                  tmpNsVec.push_back(ns[k]+1.);  
                  tmpGammaVec.push_back(g[h]);  
                  tmpT0Vec.push_back(avgtime);  
                  tmpSigVec.push_back(rms);  
                  tmpTSVec.push_back(make_pair(llhtemp,ind));
                  ind++;  
                }
            }
            tVectortemp.clear();
        }
    }

    if(tmpTSVec.size()==0) {
      tmpNsVec.push_back(nsrcBckup);  
      tmpGammaVec.push_back(gammaBckup);  
      tmpT0Vec.push_back(meanBckup);  
      tmpSigVec.push_back(rmsBckup);  
      tmpTSVec.push_back(make_pair(llhMax,0));
    }

    tVectorclose.clear();
    ind=0;
    //let's do a finer scan through gamma values and save results in vectors
    for(unsigned int i=0; i<tmpTSVec.size(); i++) {
      for (double d=1.; d<4.; d+=0.2) { // loop over gamma with 15 steps for best seed
        llhtemp = EvaluateLlh( tmpNsVec[i], d, tmpT0Vec[i], tmpSigVec[i] );
        if(llhtemp>0.) {
    nsVec.push_back(tmpNsVec[i]);
          gammaVec.push_back(d);
          t0Vec.push_back(tmpT0Vec[i]);
          sigVec.push_back(tmpSigVec[i]);
          TSVec.push_back(make_pair(llhtemp,ind));
      ind++;
        }           
        if(llhtemp>llhMax) {
          llhMax = llhtemp;
          gammaBckup = d;
        }
      }
    }
  
    if(TSVec.size()==0) {
      nsVec.push_back(nsrcBckup);  
      gammaVec.push_back(gammaBckup);  
      t0Vec.push_back(meanBckup);  
      sigVec.push_back(rmsBckup);  
      TSVec.push_back(make_pair(llhMax,0));
    }

    sort(TSVec.begin(), TSVec.end(), sortinrev);
    tmpTSVec.clear();
    tmpNsVec.clear();
    tmpGammaVec.clear();
    tmpSigVec.clear();
    tmpT0Vec.clear();

    //now select only non-overlapping flares
    Guess_mean.push_back(t0Vec[TSVec[0].second]);
    Guess_rms.push_back(sigVec[TSVec[0].second]);
    Guess_nsrc.push_back(nsVec[TSVec[0].second]);
    Guess_gamma.push_back(gammaVec[TSVec[0].second]);
    Guess_TS.push_back(TSVec[0].first);

    double t0 = t0Vec[TSVec[0].second];
    double sigma = sigVec[TSVec[0].second];
    double tmin = t0 - nSigmaTrunc_*sigma;
    double tmax = t0 + nSigmaTrunc_*sigma;

    tminVec.push_back(tmin);
    tmaxVec.push_back(tmax);

    //multi-flare overlap removal    
    bool isOverlap=false;
    for(unsigned int i=1; i<t0Vec.size(); i++) {
      t0  = t0Vec[TSVec[i].second];
      sigma = sigVec[TSVec[i].second];
      tmin = t0 - nSigmaTrunc_*sigma;
      tmax = t0 + nSigmaTrunc_*sigma;
      for(unsigned int j=0; j<tminVec.size(); j++) {
        if((tmin >= tminVec[j] && tmax <= tmaxVec[j]) ||                        // 1  2  2  1
        (tmin <= tminVec[j] && tmax >= tmaxVec[j]) ||                           // 2  1  1  2 
        (tmin <= tminVec[j] && tmax <= tmaxVec[j] && tmax >= tminVec[j]) ||     // 2  1  2  1
        (tmin >= tminVec[j] && tmax >= tmaxVec[j] && tmin <= tmaxVec[j])        // 1  2  1  2 
        ) {
          isOverlap=true;
          break;
        }
     } 
     if(!isOverlap) {
       //record parameters of this flare
       Guess_mean.push_back(t0Vec[TSVec[i].second]);
       Guess_rms.push_back(sigVec[TSVec[i].second]);     
       Guess_nsrc.push_back(nsVec[TSVec[i].second]);
       Guess_gamma.push_back(gammaVec[TSVec[i].second]);
       Guess_TS.push_back(TSVec[i].first);                            
       tminVec.push_back(tmin);
       tmaxVec.push_back(tmax);
    }
    else isOverlap=false;
  }
}

    
double MultiGaussAnalysisFn::EvaluateLlh(double *parValueArray) {
    vector<double> parVect;
    for (int i=0; i<nPar; ++i) {
        parVect.push_back(parValueArray[i]);
    }
    double minusLlh = EvalFCN(parVect);
    return -minusLlh;   // that is, max llh = - (minimizer result)
}

double MultiGaussAnalysisFn::EvalFCN(const vector<double>& parVect) const {
  double f;
  const int npar = parVect.size();
  double gin[npar];
  double par[npar];
  int iflag = 0;
  for(int i=0; i<npar; i++) par[i] = parVect[i];
  int npars=npar;
  llhFncMS(npars, gin, f, par, iflag);

  return f;
}

double MultiGaussAnalysisFn::GetProbFromHisto(double teststat) const { 
    int bin = pvalHisto_->FindBin(teststat);
    double ptemp = pvalHisto_->GetBinContent(bin);
    return ptemp;
}

void MultiGaussAnalysisFn::SetNullTestStat(TH1D * inputhisto) {
  
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
