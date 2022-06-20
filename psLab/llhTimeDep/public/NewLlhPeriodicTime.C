#include "llhTimeDep/public/NewLlhPeriodicTime.h"
#include "llhTimeDep/public/TimePdfCollection.h"

#include "TGraph.h"

#include "rootExt/public/FunctionsRoot.h"
#include "rootExt/public/generalfunctions.h"
#include "rootExt/public/log_report.h"
#include "rootExt/public/ModDistanceFn.h"

#include "fluxus/public/FluxFunction.h"

#include "llh/public/BkgSpaceProb.h"
#include "llh/public/EnergyProb.h"
#include "llh/public/I3Event.h"

  // Welcome to LlhPeriodicTime.C!
  
  // This is meant to take a series of (already loaded)
  // I3Events, and a signal location and tests for 
  // compatability with a repeating Gaussian flare (in time).
  
  // If you want to test for a flare from a specific time
  // and duration you can set it by hand in the macro you
  // use to set up the search:
  
  //  llhEnergyFn.SetParDef(2, mean,         0.,  0., 1., true, false);
  //  llhEnergyFn.SetParDef(3, log10(sigma), 0., -8., 4., true, false);
  //
  //         ---> SetParDef(par#, init, step_size, min, max, optFix, optAuto)
  
  // will look for a flare specifically at mean, with width sigma
  // and will not look for the best flare (optFix==true)

void llhPeriodicFnc(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  double srcFrac = par[0]/tdLlhPeriodic->Get_nEvents();   // relative weight of src term
  
  vector<const EnergyProb*> eProbVect = tdLlhPeriodic->GetEProbVect();
  EventPtrList selectedList = tdLlhPeriodic->GetSelectedList();
  vector<double>* eventRatioVect = tdLlhPeriodic->GetEventRatios();
  vector<double> tVect = tdLlhPeriodic->GetTVect();

  if (srcFrac > tdLlhPeriodic->GetOptimizeSrcFracMax()) {
    log_warn("Trying to evaluate srcFrac=%f which is > srcFracMax_=%f",
	     srcFrac,tdLlhPeriodic->GetOptimizeSrcFracMax());
    log_warn("Check your parameter settings and make them consistent.");
    log_warn("(If running optimized, logLambda tolerance could be exceeded.)");
  }

  if (tdLlhPeriodic->GetOptStoreRatios()) { eventRatioVect->clear(); }

  double gamma = par[1];  // Spectral index
 
  if (tdLlhPeriodic->GetUseEnergy() && (gamma >= eProbVect[0]->GetGammaMax() || gamma < eProbVect[0]->GetGammaMin() || TMath::IsNaN(gamma) ) ) {
    f = 1e50;
    return;
  }
  
  //double tmin = tdLlhPeriodic->tmin_;
  //double tmax = tdLlhPeriodic->tmax_;

  double LogLambda = 0.;
  double TimeRatio = 1.;

  double mean = par[2];    // The peak of the Gaussian, in MJD
  if (mean>1.) { mean--; }
  if (mean<0.) { mean++; }
  double sigma = pow(10.,par[3]); // I use the log10 of the sigma, since it can vary over many
                           // orders of magnitude.
    
  if ( TMath::IsNaN(mean) || TMath::IsNaN(sigma) ) { // Sometimes something causes the FCN
                                                     // to test crazy values, here we just
                                                     // say to give a very high answer.
    //cout << "Sigma: " << sigma << " = exp(" << par[3] << ") is causing troubles?" << endl;
    f = 1e50;
    return;
  } 

  TimePdf * timePdf = new PeriodicGaussianTimePdf(0, 1, mean, sigma,1.);
  //timePdf->SetLivetime( 1. );

  int nevs=0; // number of events which contribute S/B > 1

  double lh;
  vector<double> spaceRatioVect   = tdLlhPeriodic->GetProbPairRatios();
  double b=1.;
  
  for (int i=0;i<tdLlhPeriodic->GetNEventsSelected();i++)  { // loop over selected events

    const Event* event = selectedList.GetEvent(i);

    TimeRatio = timePdf->GetPdfValue( tVect[i] ) * tdLlhPeriodic->GetLivetime() / b;
 
    // Perhaps we can speed things up... if the TimeRatio==0, skip the rest of the calculations.
    // include an out for high monitoring levels...
    if (tdLlhPeriodic->GetMonitorLevel() < 4 && TimeRatio < 1e-10) {
      LogLambda += log(1.-srcFrac);
      continue;
    }


    // first comes the event's energy term
    double eRatio=1.;
    if (tdLlhPeriodic->GetUseEnergy()) {
      eRatio = eProbVect[i]->GetEnergyProbGamma(*event, gamma) /
	           eProbVect[i]->GetEnergyProbBkg(*event);
    }
        
    //As implemented, local coord background terms will automagically be included in spaceRatoVect   
    lh = ( srcFrac * ( spaceRatioVect[i] * eRatio * TimeRatio - 1. ) + 1. );
    
    if(lh > 1) { nevs++; }
    LogLambda += log( lh );

    if ( (tdLlhPeriodic->GetMonitorLevel() > 3 && lh > 1.) || tdLlhPeriodic->GetMonitorLevel() > 4 ) { 
                                      // useful monitoring but with massive event-by-event spew
                                      // (~200 lines) going through each FCN call (~100-300)
      cout << lh                      << " " << 
              LogLambda               << " " << 
              spaceRatioVect[i] << " " << 
              eRatio                  << " " << 
              tVect[i]          << " " << 
              fmod(event->GetTime().GetMJD()-tdLlhPeriodic->GetEphemeris(),tdLlhPeriodic->GetPeriod())/tdLlhPeriodic->GetPeriod() << " " << 
              timePdf->GetPdfValue(tVect[i]) << endl;
    }

    if (tdLlhPeriodic->GetOptStoreRatios()) {
       eventRatioVect->push_back(spaceRatioVect[i]*eRatio*TimeRatio);
    }

  } //end loop over events
  
  if (tdLlhPeriodic->GetOptimizeTolerance() > 0.) { // correction for events skipped in optimization
    LogLambda += (tdLlhPeriodic->Get_nEvents()-tdLlhPeriodic->GetNEventsSelected())*log(1.-srcFrac);
  }

  double nn = tdLlhPeriodic->GetWeightPower(); 
  double mn = 1.;
  if (tdLlhPeriodic->JimsTerm_ && sigma>0.1) { mn = TMath::Erf(0.5/sigma); }
    // This (mn) is an attempt at a correction factor as the periodic Gaussian timePDF
    // only uses the PDF value from the most likely phase, so the integration
    // limits for the weighting term will be different.
  
  //The marginalization needs to be < 1. and log < 0. 
  if ( tdLlhPeriodic->JimsTerm_ && log(pow( mn * sigma * 2.0/ tdLlhPeriodic->GetLivetime(), nn )) < 0. ) { 
    
    LogLambda += log( pow( mn * sigma * 2.0/ tdLlhPeriodic->GetLivetime(), nn ) );
  }
 
  tdLlhPeriodic->Set_nevs(nevs);
  tdLlhPeriodic->SetMargValue(log( pow( mn * sigma * 2.0/ tdLlhPeriodic->GetLivetime(), nn ) ));

  // Jim's Spectral index penalty
  // We could expect a gamma between 2.0 and 2.7
  // he applies a Gaussian with sigma 0.2 units in spectrum outside of that range
  if (tdLlhPeriodic->SpectralPenalty_ && tdLlhPeriodic->GetUseEnergy()) { 
    if (gamma < 2.) LogLambda -= (2.-gamma)*(2.-gamma)/0.08;
    else if (gamma > 2.7) LogLambda -= (2.7-gamma)*(2.7-gamma)/0.08;
  }

  f = - LogLambda;    // What Minuit minimizes: -log L
  
  if ( tdLlhPeriodic->GetMonitorLevel()>1 ) { // more useful monitoring, with only
                           // one printout per FCN call. Still lots of spew
                           // but it's pretty managable for one or a few trials.
    printf("LogLambda=%12.6lg : nSrc=%9.4lg : gamma=%5.4lg : mean=%8.8lg : sima=%8.4lg : nEvs=%i : marg_fact=%g : mn=%f\n",
	   LogLambda, srcFrac*tdLlhPeriodic->Get_nEvents(), gamma, mean, sigma, tdLlhPeriodic->Get_nevs(), tdLlhPeriodic->Get_MargWeight(), mn );
  }  
    
  if (timePdf) { delete timePdf; } // otherwise we'd get a memory leak
  
  return;
}



NewLlhPeriodicTime::NewLlhPeriodicTime() :
  i3Set_(NULL),
  useLCBkgProb_(false),
  nullTestStat_(NULL),
  pvalHisto_(NULL),
  //fitfn_(NULL),
  histoForProb_(false),
  optUseEnergy_(true),
  eMaxRatioWarnStatus_(-1),  // default -1 means warn every time
  optStoreRatios_(false),
  icstatWarnLevel_(0),
  nSrcMin_(0.),
  srcFracMax_(0.5), // default is half of total nEvents
  gammaMin_(0.),
  gammaMax_(0.),
  logLambdaBest_(0.),
  nSrcBest_(0.),
  gammaBest_(0.),
  chiSq_(0.),
  chiSqProb_(0.),
  monitorLevel_(0),
  optimizeTolerance_(0.),
  optimizeAngleDeg_(0.),
  nEventsTot_(0),
  margValue_(1.),
  fitfrac_(0.001),
  seedWtMin(1000)
{
  tdLlhPeriodic = this;

  optParAuto_[0] = true;
  optParAuto_[1] = true;
  optParAuto_[2] = true;
  optParAuto_[3] = true;

  // These start when created, so stop immediately
  stopwatch_MaximizeLlh_.Stop();
  stopwatch_optimize_.Stop();
  stopwatch_minuitMigrad_.Stop();
  stopwatch_minuitFCN_.Stop();
}


void NewLlhPeriodicTime::OptimizeEventSelection() {
  const EquatorialDeg *srcEqDeg = 
    ( dynamic_cast<const EquatorialDeg*>(srcCoord_) );

  double decMinDeg = -90.;
  double decMaxDeg =  90.;
  double raRangeDeg = 360.;
  if (optimizeAngleDeg_ > 0.) {
    decMinDeg = srcEqDeg->GetDec() - optimizeAngleDeg_;
    decMaxDeg = srcEqDeg->GetDec() + optimizeAngleDeg_;
    if (decMinDeg < -90.) { decMinDeg = -90.; }
    if (decMaxDeg >  90.) { decMaxDeg =  90.; }
    // scale the ra range according to the dec closest to one of the poles
    double cosFactor = cos(decMinDeg*TMath::DegToRad());
    double cosFactorAlt = cos(decMaxDeg*TMath::DegToRad());
    if (cosFactorAlt < cosFactor) { cosFactor = cosFactorAlt; }
    if (cosFactor>0) {
      raRangeDeg = optimizeAngleDeg_ / cosFactor;
    }
  }

  double threshold = 
    (1./srcFracMax_ - 1.) * ( exp(optimizeTolerance_/nEventsTot_) -1. );

  const EventPtrList* evList = aSet_->GetEventPtrList();

  for (int i=0; i<nEventsTot_; ++i) {
    const I3Event* event = (dynamic_cast<const I3Event*> (evList->GetEvent(i)));
    
    if (!event) { log_fatal("OptimizeEventSelection: failed cast to event.\n");}

    double eventRaDeg = event->GetEquatorialDeg().GetRa();
    // use this function to test whether ra is within range
    if ( ModDistanceFn(eventRaDeg,srcEqDeg->GetRa(),360.) > raRangeDeg ) {
      continue;  // skip this event
    }

    double eventDecDeg = event->GetEquatorialDeg().GetDec();
    if (eventDecDeg < decMinDeg || eventDecDeg > decMaxDeg) {
      continue;  // skip this event
    }

    double sigSpaceProb = event->ProbFrom(*srcCoord_);
    double bkgSpaceProb = event->GetBkgSpaceProbFn()->GetBkgProbDensity(*event);
        
    if ( useLCBkgProb_ ) {
      double lcProb = lcBkgProb_->BackgroundLCProb(*event);
      bkgSpaceProb = bkgSpaceProb*lcProb;
    }
    
    //    if event is in region with b=0 (what to do?)
    if (bkgSpaceProb <=0) { log_fatal("bkgSpaceProb <= 0\n"); }
    double spaceRatio = sigSpaceProb / bkgSpaceProb;

    double eMaxRatio;
    const EnergyProb* eProb(NULL);
    if (optUseEnergy_) { 
      eProb = event->GetEnergyProbFn();
      eMaxRatio = eProb->GetEnergyMaxRatio(*event);
      if (eMaxRatio <= 0.) {
	// Decide what sort of warning to give:
	if (eMaxRatioWarnStatus_ < 0) {          // Always warn
	  log_error("Error: gamma table is zero for this energy.  "
		    "Max Ratio=0.\n");
	  cout << event->GetParams().energyValue << endl;
	} else if (eMaxRatioWarnStatus_ == 0) {  // Warn Once
	  log_error("***\nWARNING!\nAt LEAST one data event has zero for "
		    "its gamma energy pdf.\nAll such events will be "
		    "effectively ELIMINATED when using llh with energy.\n"
		    "Llh is currently set so that you will receive this "
		    "warning ONLY ONE TIME.\n***\n");
	  eMaxRatioWarnStatus_ = 1;
	}
	// else don't warn again
      }
    } else {
      eMaxRatio = 1.;
    }
    
    if (spaceRatio*eMaxRatio > threshold) {  // select this event
      selectedList_.AddEvent(event);
      eVect_.push_back(*event);
      tVect_.push_back( event->GetTime().GetMJD() );
      eProbVect_.push_back(eProb);
      spaceRatioVect_.push_back(spaceRatio);
      if (monitorLevel_ > 2) {
        cout << "(" << eventRaDeg << "," << eventDecDeg << ")  Sr : " << sigSpaceProb << "/" << bkgSpaceProb << " eP : " << eMaxRatio << endl;
      }
    }
  }

  if (monitorLevel_ > 0) {
    printf("Optimizing: %d events selected with maxRatio >%lg"
	   "  out of %d Ntotal.\n",
	   selectedList_.GetSize(), threshold, nEventsTot_);
  }
}



void NewLlhPeriodicTime::PrepareAnalysis() {
  if (!aSet_) { log_fatal("PrepareAnalysis: AnalysisSet was not set.\n"); }
  if (!srcCoord_) { log_fatal("PrepareAnalysis: srcCoord was not set.\n"); }

  const EventPtrList* evList = aSet_->GetEventPtrList();
  if (!evList) { log_fatal("PrepareAnalysis: EventPtrList was not set.\n"); }
  nEventsTot_ = evList->GetSize();
  if (!nEventsTot_) { log_fatal("PrepareAnalysis: EventPtrList was empty.\n"); }
  // or, is there a reason to allow zero events?

  selectedList_.Clear();
  eVect_.clear();
  eProbVect_.clear();
  spaceRatioVect_.clear();
  tVect_.clear();

  eventRatioVect_.clear();

  if (optimizeTolerance_ > 0.) {
    stopwatch_optimize_.Start(false);  //  false = don't reset
    OptimizeEventSelection();
    stopwatch_optimize_.Stop();
  } 
  else 
  { // no optimization
    for (int i=0; i<nEventsTot_; ++i) {
      const I3Event* event = 
	(dynamic_cast<const I3Event*> (evList->GetEvent(i)));
      assert(event);
      double sigSpaceProb = event->ProbFrom(*srcCoord_);
      double bkgSpaceProb = 
	event->GetBkgSpaceProbFn()->GetBkgProbDensity(*event);
      //    if event is in region with b=0 (what to do?)
      if (bkgSpaceProb <=0) { log_fatal("bkgSpaceProb <= 0\n"); }
      double spaceRatio = sigSpaceProb / bkgSpaceProb;
      const EnergyProb* eProb(NULL);
      if (optUseEnergy_) { 
	eProb = event->GetEnergyProbFn();
      }
      eVect_.push_back( *event );
      tVect_.push_back( event->GetTime().GetMJD() );
      selectedList_.AddEvent(event);
      eProbVect_.push_back(eProb);
      spaceRatioVect_.push_back(spaceRatio);
    }
  }
  
  nEventsSelected_ = eVect_.size();
  
  for (int i=0;i<int(tVect_.size());i++) {
    if(tVect_[i]>1.) { // This is because signal events are generated with
                       // a phase [0,1] instead of an MJD. Make sure that
                       // the events ARE generated on the interval [0,1].
      tVect_[i] = fmod(tVect_[i]-t0_,period_)/period_;
    }
  }

  if (optUseEnergy_) {
    gammaMin_ = eProbVect_[0]->GetGammaMin();
    gammaMax_ = eProbVect_[0]->GetGammaMax();
  }
  
}


void NewLlhPeriodicTime::MaximizeLlh()
{
  stopwatch_MaximizeLlh_.Start(false);  //  false = don't reset

  PrepareAnalysis();
  // TO DO: BETTER HANDLING IF NO EVENTS SELECTED:
  assert(eVect_.size());
  //assert(selectedList_.GetSize());

  // sigmamin_ is set up to be the minimimum witdh of the
  // time Gaussian. If you let the width become arbitrarily small
  // the method could find that it can keep shrinking onto one strong event
  // and keep getting a better test statistic, which isn't what we want.
  // Sets a lower bound of the rms between the closest two events in time.

  minuit_ = new TMinuit(4);
  Double_t arglist[1];
  arglist[0] = 0.5; //0.5 for likelihood fit
  Int_t ierflg=0;
  minuit_->SetFCN(llhPeriodicFnc);
  minuit_->mnexcm("SET ERR",arglist,1,ierflg);
  minuit_->SetPrintLevel(-1);

  // DEFINE PARAMETERS

  GetFlareGuess(optUseEnergy_, 1, nsrcGuess_, gammaGuess_, meanGuess_, sigmaGuess_);

  if (optParAuto_[0]) {
    double nSrcMax = srcFracMax_*nEventsTot_;
    nSrcMin_ = 0.; 
    if (!nsrcGuess_) { nsrcGuess_=2.5; }
    minuit_->mnparm(0, "nSrc", nsrcGuess_, 1., nSrcMin_, nSrcMax, ierflg);
  }
 
  if (optParAuto_[1]) {
    if (optUseEnergy_) {
      gammaMin_ = eProbVect_[0]->GetGammaMin();
      gammaMax_ = eProbVect_[0]->GetGammaMax();
      
      double gammaInit = (gammaMin_+gammaMax_)/2.;
      
      gammaInit = gammaGuess_;
      minuit_->mnparm(1, "gamma", gammaInit, 0.2, gammaMin_, gammaMax_, ierflg);
    } else {
      // If not using energy, then FCN will actually replace the
      // energy prob term internally... so these values don't really matter
      // Just fix the gamma so that Minuit doesn't try to minimize w.r.t. it
      minuit_->mnparm(1, "gamma", 0., 0., 0., 0., ierflg);
    }
  }
  
  if (optParAuto_[2]) {
      double lowLimit = -1.;
      double upLimit =  2.;
      double  initStepSize = sigmaGuess_;
      double  initValue = meanGuess_;
      minuit_->mnparm(2, "mean", initValue, initStepSize, lowLimit, upLimit, ierflg);
  }
   
  if (optParAuto_[3]) {
      // Using the log10 of sigma as the fit parameter to make things more smooth.
      double initValueS = log10(sigmaGuess_);
      double lowLimitS = log10(sigmamin_);
      double upLimitS = log10(livetime_);
      double initStepSizeS = 1.;
      minuit_->mnparm(3, "sigma", initValueS, initStepSizeS, lowLimitS, upLimitS, ierflg);
  }

  // MINIMIZE

  stopwatch_minuitMigrad_.Start(false);  //  false = don't reset
  double arglist2[2];
  arglist2[0] = 500;
  arglist2[1] = 0.1;
  minuit_->mnexcm("MIGRAD", arglist2 ,2,minuitOut_);
  stopwatch_minuitMigrad_.Stop();

  // GET RESULTS

  StoreLogLambdaBest();   // Set logLambdaBest_
  nSrcBest_  = GetPar(0);
  gammaBest_ = GetPar(1);
  meanBest_  = GetPar(2);
  
  while(meanBest_<0.) {meanBest_ = meanBest_+1.;}
  while(meanBest_>1.) {meanBest_ = meanBest_-1.;}
  
  sigmaBest_ = pow( 10., GetPar(3) );

  stopwatch_MaximizeLlh_.Stop();
}

///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  END MAXZIMIZELLH()  ///////////////////////
///////////////////////////////////////////////////////////////////////////////

// This is clumsy but seems like the only way to get this info
void NewLlhPeriodicTime::StoreLogLambdaBest() 
{
  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, istat;
  minuit_->mnstat(amin, edm, errdef, nvpar, nparx,istat);
  logLambdaBest_ = -amin;  // that is, max llh = - (minimizer result)

  // The *worst* logLambdaBest should be zero (i.e. null hypothesis).
  // But minimizer will miss exact zero, leading logLambdaBest_ to be slightly
  // negative. Since this will cause probability calculation to choke, we 
  // fix it here
  
  if ( monitorLevel_>1 ) { 
    printf("LogLambdaBest=%12.6lg\n",logLambdaBest_ );
  } 
  
  if (logLambdaBest_ < 0.) { logLambdaBest_ = 0.;}
}

// This function scans over all the events which could have an S/B ratio > 10,
// and tests consecutive doubles, triples... for compatability with and E^-2 flare.
// I tried testing for Gamma = 2 and like 3.5 before, but it doesn't seem to do much.
// This is configurable enough to do whatever you'd like.
// There is a final scan over Gamma at the end to round out the set parameters.

// Oh, I'm using this parameter close_ to set how far into the vector 1,2,3,4,5,10,15,20,25
// to test for flare compatability. You only need to go to 10, I usually go to 15 (close_=6).

void NewLlhPeriodicTime::GetFlareGuess(bool useE, int useSpace, double & Guess_nsrc, double & Guess_gamma, double & Guess_mean, double & Guess_rms) 
{

  if (monitorLevel_ > 1) { cout << "starting initial scan with " << tVect_.size() << " events." << endl; }

  vector<double> tVectortemp;
  vector<double> tVectorclose;

  double llhMax= -100.;
  double llhtemp, rms, avgtime;

  for (int j=0; j<(int)tVect_.size(); j++) {
    if ( (eProbVect_[j]->GetEnergyMaxRatio(eVect_[j])*spaceRatioVect_[j]) > seedWtMin ) {
        if (monitorLevel_ > 2) { cout << tVect_[j] << " " << flush; }
        tVectorclose.push_back( tVect_[j] );
        tVectorclose.push_back( tVect_[j]+1 ); //hopefully this will deal with the period better...
    }
  }
  if (monitorLevel_ > 2) { cout << endl; }
  
  //I didn't sort this vector when changing to the phase, so we should do it here
  sort(tVectorclose.begin(),tVectorclose.end());

  sigmamin_ = 1000.;
  double timed;
  int sp = 1;

  for (unsigned int i=sp; i<tVectorclose.size(); i++) {//This assumes times are phase-ordered.
    timed = tVectorclose[i] - tVectorclose[i-sp];
    if (fabs(timed) < sigmamin_) {sigmamin_ = fabs(timed);}
  }
  

  const int ngamma=2;
  double g[ngamma] = {2.0, 3.5};
//  const int ngamma=1;
//  double g[ngamma] = {2.0};

  int kmax = (int) close_;
  if (10 < kmax) { kmax = 10; }

  int n[] =  {1,2,3,4,5,10,15,20,35,int(tVectorclose.size())/3,int(tVectorclose.size())/2};
  int ns[] = {1,2,2,2,2, 3, 3, 4, 4, 5, 5};

  int kmin=0;
  
    if ( sigmamin_<1e-8 ) { sigmamin_=1e-8; } // this seems to get set to zero sometimes, which is problematic in taking the log later...
  
  if (monitorLevel_ > 0) { cout << tVectorclose.size() << " events with S/B > " << seedWtMin << ". Sigmamin is " << sigmamin_ << endl; }

  for (int mscan=0; mscan<10; mscan++) {
    
    avgtime = mscan/10. + 1./20.;
    rms = 0.2;
  
    llhtemp = EvaluateLlh( 2., 2., avgtime, rms );

    if (llhtemp > llhMax) {
      llhMax = llhtemp;
      Guess_mean = avgtime;
      Guess_rms = rms;
      Guess_nsrc = 2.;
    }
  }

  llhMax = -100.;

  for (int m=1;m<20;m++) {
    
    llhtemp = EvaluateLlh( m/2., 2., Guess_mean, Guess_rms );

    if (llhtemp > llhMax) {
      Guess_nsrc = m/2.;
    }
  } 
  
  while ( tVectorclose.size()*1.0 < n[kmax] ) { kmax--; }

  for (int k=kmin; k<kmax; k++){ //k==1 for pairs, k==2 for triples, etc.
  
    for (unsigned int i=0; i<(tVectorclose.size()-n[k]); i++) {
      for (int j=0;j<=n[k];j++){
        tVectortemp.push_back(tVectorclose[i+j]);
      }
            
   //get mean of the times to seed this gaussian
      double sum=0, sumsq=0;
      llhtemp = 0.;
      for (unsigned int ii=0;ii<tVectortemp.size();ii++) {
        sum += tVectortemp[ii];
      }
      avgtime = sum/tVectortemp.size();

      for (unsigned int ii=0;ii<tVectortemp.size();ii++) {
        sumsq += (tVectortemp[ii] - avgtime)*(tVectortemp[ii] - avgtime);
      }

    //calculate the rms to use as the sigma
      rms = sqrt( sumsq/tVectortemp.size() );

      if (rms < 0.001) {
        tVectortemp.clear();
        continue;
      }

      if (useE) {
        for (int h=0; h<ngamma; h++) {
          llhtemp = EvaluateLlh( ns[k]+1., g[h], avgtime, rms );
          if (llhtemp > llhMax) {
            llhMax = llhtemp;
            Guess_mean = avgtime;
            Guess_rms = rms;
            Guess_nsrc = ns[k]+1.;
          }
        }
        tVectortemp.clear();
      } else { // I have no idea why you'd want to do the analysis without energy, but here you go...
        llhtemp = EvaluateLlh( ns[k]+1., 0., avgtime, rms );
        if (llhtemp > llhMax) {
          llhMax = llhtemp;
          Guess_mean = avgtime;
          Guess_rms = rms;
          Guess_nsrc = ns[k]+1.;
        }
        tVectortemp.clear();
      }
    }
  } //*/
  
  tVectorclose.clear();
  
  double sllhmax = -100;
  double sllhtemp;
  for (double d=1.; d<4.; d+=0.2) { // loop over gamma with 15 steps for best seed
    sllhtemp = EvaluateLlh( Guess_nsrc, d, Guess_mean, Guess_rms );
    if (sllhtemp > sllhmax ) {
      Guess_gamma = d;
      sllhmax = sllhtemp;
    }
  }
  
  if (monitorLevel_ > 0) { printf("Guess Params: %f %f %f %f \n",Guess_nsrc,Guess_gamma,Guess_mean,Guess_rms); }
  
  if (useE || useSpace) { } //nip compiler complaints
  
}


void NewLlhPeriodicTime::SetNullTestStat(TH1D * inputhisto) {
  
  // This reads in a specific TH1D as the null test statistic
  // to use for p-values instead of using a chisquare distribution
  // it also fits an exponential to the upper tail (default set to top 0.01%).
  
  // It takes in the raw test statistic distribution (pdf)
  // and then makes the cdf to work with later.
  
  histoForProb_ = true;
  
  pvalHisto_ = new TH1D();
  nullTestStat_ = new TH1D();
  
  char sc[] = "Scale";
  //double firstfit=0;
  pvalHisto_ = DescendingCumulate(inputhisto,sc);
  nullTestStat_ = (TH1*)inputhisto->Clone("hnew");
  
  //int bins = inputhisto->GetNbinsX();
  //for (int i=1;i<bins;i++) {
  //  if (pvalHisto->GetBinContent(i) < fitfrac_){
  //    firstfit = pvalHisto->GetBinCenter(i);
  //    break;
  //  }
  //}
   
  //double max = 2*inputhisto->GetBinLowEdge(bins) - inputhisto->GetBinLowEdge(bins-1);
  //fitfn_ = new TF1("fitfn_","exp([0]+[1]*x)",firstfit,100);
   
  //pvalHisto->Fit(fitfn_,"QO","",firstfit,max);
  
}

double NewLlhPeriodicTime::GetProbFromHisto(double teststat){ //, bool useFit) {

  // this looks at a null teststat distribution
  // which we loaded, and if it's in the exponential
  // tail, should go with a fit expo->Eval(teststat) (top 1%)
  
  int bin = pvalHisto->FindBin(teststat);
  double ptemp = pvalHisto->GetBinContent(bin);

//  if (ptemp < fitfrac_ && useFit) { //small enough that we should do the histo fit
//    ptemp = fitfn_->Eval(teststat);
//  }

  // in case we don't feel like using a fit and it is higher than the distribution
  //   // just return 1/trials (conservative?)
  if( !pvalHisto->GetBinContent(bin) ) { return 1.0/nullTestStat_->GetEntries(); }

  return ptemp;
  
}


double NewLlhPeriodicTime::EvalFCN(const vector<double>& parVect) const {
  assert(nEventsTot_);
  double f;
  const int npar = parVect.size();
  double gin[npar];
  double par[npar];
  int iflag = 0;
  for(int i=0; i<npar; i++) par[i] = parVect[i];
  int npars=npar;
  llhPeriodicFnc(npars, gin, f, par, iflag);

  return f;
}



double NewLlhPeriodicTime::EvaluateLlh(double nSrc, double gamma, double mean, double sigma) {
  vector<double> parVect;
  parVect.push_back(nSrc);  parVect.push_back(gamma);
  parVect.push_back(mean);  parVect.push_back( log10(sigma) );
  double minusLlh = EvalFCN(parVect);
  return -minusLlh;   // that is, max llh = - (minimizer result)
}

// TGraph* NewLlhPeriodicTime::GetContour(double sigma, int npoints, int pa1, int pa2) {
//   minuit_->SetFCN(MinuitWrapperFCN);
//   SetMinuitWrapperPtr(this);  // point to *our* EvalMinuitFCN

//   minuit_->SetErrorDef(sigma*sigma);
//   TGraph* g = dynamic_cast<TGraph*> (minuit_->Contour(npoints, pa1, pa2));

//   SetMinuitWrapperPtr(NULL);  // reset pointer
//   return g;
// }




// For MultiAnalysis: Translate Global Par values to Individual Par Values:

const vector<double> NewLlhPeriodicTime_ParTranslator::
Translate(int llhIndex, const vector<double>& parGlobal) const {

  //For now assume things are fine with only using the energy weight

  double nSrcGlobal = parGlobal[0];
  double gamma = parGlobal[1];

  vector<double> parLocal(4);

  double weight = 
    LinearInterpolate(gamma, gammaMin_, gammaMax_, srcWeightVect_[llhIndex]);
  // recall: srcWeightVect_[i] is itself a vector<double> 

  // nSrc , must be scaled to relative weight of this source
  parLocal[0] = nSrcGlobal * weight;
  parLocal[1] = gamma;
  parLocal[2] = parGlobal[2];
  parLocal[3] = parGlobal[3];
  return parLocal;
}
  

void NewLlhPeriodicTime_ParTranslator::
SetTranslator(const vector<AnalysisSet*>& aSetVect) {
  int nSets = aSetVect.size();

  srcWeightVect_.clear();
  vector<double> tempVect(nStopsGamma_,0);
  for (int i=0; i<nSets; ++i) {
    srcWeightVect_.push_back( tempVect );
  }

  for (int g = 0; g < nStopsGamma_; ++g) {
    double gamma = gammaMin_ + 
      (gammaMax_-gammaMin_) * (double(g))/(nStopsGamma_-1);
    double sum = 0.;
    for (int i=0; i<nSets; ++i) {
      double nev = 
	aSetVect[i]->GetMeanSrcNevForFluxModel(PowerLawFlux(1,-gamma));
      srcWeightVect_[i][g] = nev;
      sum += nev;
    }
    // for each choice of gamma, normalize weights of different data sets to 1
    cout << "Weights for Gamma=" << gamma << ": ";
    for (int i=0; i<nSets; ++i) {
      srcWeightVect_[i][g] /= sum;
      cout << srcWeightVect_[i][g] << " ";
    }
    cout << endl;
  }    
}

