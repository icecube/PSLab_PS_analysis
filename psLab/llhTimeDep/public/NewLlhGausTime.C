#include "llhTimeDep/public/NewLlhGausTime.h"
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
#include "rootExt/public/randomfunctions.h"

#include <vector>
#include <functional>   // std::greater

  // Welcome to LlhGausTime.C!
  
  // This is meant to take a series of (already loaded)
  // I3Events, and a signal location and tests for 
  // compatability with a single Gaussian flare (in time).
  
  // If you want to test for a flare from a specific time
  // and duration you can set it by hand in the macro you
  // use to set up the search:
  
  //  llhEnergyFn.SetParDef(2, mean, 54561.,54971., 5e6., true, false);
  //  llhEnergyFn.SetParDef(3, log10(sigma), 0., -8., 4., true, false);
  //
  //         ---> SetParDef(par#, init, step_size min, max, optFix, optAuto)
  
  // will look for a flare specifically at mean, with width sigma
  // and will not look for the best flare (optFix==true)


// Driver function to sort the vector elements by  
// first element of pair in descending order

void llhFnc(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  int nFlares = tdLlh->GetNFlareGuess();
  double srcFrac[nFlares]; // relative weight of src term
  double gamma[nFlares];   // Spectral index
  double mean[nFlares];    // The peak of the Gaussian, in MJD
  double sigma[nFlares];   // I use the log10 of the sigma, since it can vary over many
                           // orders of magnitude.

  vector<const EnergyProb*> eProbVect = tdLlh->GetEProbVect();
  for(int iflare=0; iflare<nFlares; iflare++) {
    int i = iflare*4;
    srcFrac[iflare] = par[0+i]/tdLlh->Get_nEvents();
    gamma[iflare]   = par[1+i];
    mean[iflare]    = par[2+i];
    sigma[iflare]   = pow(10.,par[3+i]); 
    if (srcFrac[iflare] > tdLlh->GetOptimizeSrcFracMax()) {
      log_warn("Trying to evaluate srcFrac=%f which is > srcFracMax_=%f",
         srcFrac[iflare],tdLlh->GetOptimizeSrcFracMax());
      log_warn("Check your parameter settings and make them consistent.");
      log_warn("(If running optimized, logLambda tolerance could be exceeded.)");
    }
    if (tdLlh->GetUseEnergy() && (gamma[iflare] >= eProbVect[0]->GetGammaMax() || gamma[iflare] < eProbVect[0]->GetGammaMin() || TMath::IsNaN(gamma[iflare]) ) ) {
      f = 1e50;
      return;
    }
    if ( TMath::IsNaN(mean[iflare]) || TMath::IsNaN(sigma[iflare]) ) { // Sometimes something causes the FCN
                                                                       // to test crazy values, here we just
                                                                       // say to give a very high answer.
      f = 1e50;
      return;
    }
  }

  vector<double>* eventRatioVect;
  vector<double>* spaceWeightVect;
  vector<double>* enWeightVect;
  //vector<double>* enMaxWeightVect;
  vector<double>* timeWeightVect;
  vector<double>* raVect;
  vector<double>* decVect;
  vector<double>* angErrVect;
  vector<double>* eneVect;
  vector<double>* timeVect;
  vector<int>* eventID;
  vector<int>* runID; 
  
  if (tdLlh->GetOptStoreRatios()) {
    eventRatioVect     = tdLlh->GetEventRatios();
    spaceWeightVect    = tdLlh->GetSpatialWeights();
    enWeightVect       = tdLlh->GetEnergyWeights();
    //enMaxWeightVect  = tdLlh->GetEnergyMaxWeights();
    timeWeightVect     = tdLlh->GetTimeWeights();
    raVect             = tdLlh->GetraVect();
    decVect            = tdLlh->GetdecVect(); 
    angErrVect         = tdLlh->GetAngErrVect();
    eneVect            = tdLlh->GeteneVect();
    timeVect           = tdLlh->GettimeVect();
    eventID            = tdLlh->GetEventID();
    runID           = tdLlh->GetRunID();
    eventRatioVect->clear();
    spaceWeightVect->clear();
    enWeightVect->clear();
    //enMaxWeightVect->clear();
    timeWeightVect->clear();
    raVect->clear();
    decVect->clear();
    angErrVect->clear();
    eneVect->clear();
    timeVect->clear();
    eventID->clear();
    runID->clear();
  }
  
  double tmin=0, tmax=0;

  if(tdLlh->GetAnalysisSet()->GetSource()->GetTimePdf()){
    tmin = tdLlh->GetAnalysisSet()->GetSource()->GetTimePdf()->GetTmin();
    tmax = tdLlh->GetAnalysisSet()->GetSource()->GetTimePdf()->GetTmax();
  }
  else if(tdLlh->GetAnalysisSet()->GetSource()->GetTimePdfVect().at(0)){
    tmin = tdLlh->GetAnalysisSet()->GetSource()->GetTimePdfVect().at(0)->GetTmin();
    tmax = tdLlh->GetAnalysisSet()->GetSource()->GetTimePdfVect().at(0)->GetTmax();
  }
  else log_error("llhFnc: no source time pdf.\n");
 
  double LogLambda  = 0.;
  double TimeRatio  = 1.;
  double SpaceRatio = 1.;

  double raSrc  = 0;
  double decSrc = 0.;
  if(tdLlh->IsFitSrc()) {
    int ibin = 4*nFlares;
    raSrc  = par[ibin]; //RA  of the src, to be fitted
    decSrc = par[ibin+1]; //DEC of the src, to be fitted
  }

  if ( tdLlh->GetMonitorLevel() > 4 )
    for(int iflare=0; iflare<nFlares; iflare++) 
      cout << "new GaussianTimePdf( " << tmin << ", " << tmax << ", " << mean[iflare] << ", " << sigma[iflare] << ",1.)" << endl;

  GaussianTimePdf *timePdf[nFlares];
  double sumNorm = 0.;
  double norm[nFlares];
  for(int iflare=0; iflare<nFlares; iflare++) {
    timePdf[iflare] = new GaussianTimePdf(tmin, tmax, mean[iflare], sigma[iflare], 1.);
    if( tdLlh->IsSetNormFromGRL() ){
      timePdf[iflare]->SetNormFromGRL(tdLlh->GetStartMissRunVect(), tdLlh->GetStopMissRunVect());
      timePdf[iflare]->SetLivetime( tdLlh->GetLivetime() );
    }
    else timePdf[iflare]->SetLivetime( tmax - tmin );
    if(tdLlh->GetIsGaussTrunc()) timePdf[iflare]->SetNSigTrunc(tdLlh->GetNSigTrunc());
    norm[iflare] = timePdf[iflare]->GetNorm();
    sumNorm += norm[iflare];
  } 

  int nevs = 0; // number of events which contribute S/B > 1

  double lh;
  vector<double> spaceRatioVect   = tdLlh->GetProbPairRatios();
  vector<double> bkgSpaceProbVect = tdLlh->GetBkgSpaceProbVect();

  for (int i=0;i<tdLlh->GetNEventsSelected() ;i++)  { // loop over selected events
    const Event* event = tdLlh->selectedList_.GetEvent(i);
    lh=0;
    for(int iflare=0; iflare<nFlares; iflare++) {
      if(sumNorm>0) TimeRatio = /*(norm[iflare]/sumNorm) * */timePdf[iflare]->GetPdfValue( event->GetTime().GetMJD() ) * timePdf[iflare]->GetLivetime();
      else TimeRatio=0;
      //Setting space pdf
      if(tdLlh->IsFitSrc()) {
        EquatorialDeg coordSrcED(raSrc,decSrc);
        const Coord *coordSrc = dynamic_cast<const Coord*>(&coordSrcED);
        double sigSpaceProb = event->ProbFrom(*coordSrc);
        double bkgSpaceProb = bkgSpaceProbVect[i];
        if (bkgSpaceProb <=0) { log_fatal("bkgSpaceProb <= 0\n"); }
        SpaceRatio = sigSpaceProb/bkgSpaceProb;
      }

      // Perhaps we can speed things up... if the TimeRatio==0, skip the rest of the calculations.
      // include an out for high monitoring levels...

      if (tdLlh->GetMonitorLevel() < 4 && TimeRatio < 1e-50 && !tdLlh->GetOptStoreRatios()) {
        lh += -srcFrac[iflare];
        continue; 
      }

      // first comes the event's energy term
      double eRatio=1.;
      if (tdLlh->GetUseEnergy()) {
        eRatio = eProbVect[i]->GetEnergyProbGamma(*event, gamma[iflare]) /
                 eProbVect[i]->GetEnergyProbBkg(*event);
      }
      //As implemented, local coord background terms will automagically be included in spaceRatoVect   
      if(!tdLlh->IsFitSrc())  SpaceRatio = spaceRatioVect[i];
      lh += srcFrac[iflare] * ( SpaceRatio * eRatio * TimeRatio - 1. );
      if ( (tdLlh->GetMonitorLevel() > 3 && lh > 1.) || tdLlh->GetMonitorLevel() > 4 ) {
                                      // useful monitoring but with massive event-by-event spew
                                      // (~200 lines) going through each FCN call (~100-300)
        const I3Event* i3ev = (dynamic_cast<const I3Event*>(event));

        cout << "FLARE "<<iflare<<":   " << lh << " " << LogLambda << " " << SpaceRatio << " " <<
        eRatio << " " << event->GetTime().GetMJD() << " " <<
        timePdf[iflare]->GetPdfValue(event->GetTime().GetMJD()) << " " <<
        i3ev->GetParams().runID << " " << i3ev->GetParams().eventID << endl;
       
      }
      if (tdLlh->GetOptStoreRatios()) {
        const I3Event* i3ev = (dynamic_cast<const I3Event*>(event));
        eventRatioVect->push_back(SpaceRatio*eRatio*TimeRatio);
        spaceWeightVect->push_back(SpaceRatio);
        enWeightVect->push_back(eRatio);
        //const EnergyProb* eProb = i3ev->GetEnergyProbFn();
        //double eMaxRatio = eProb->GetEnergyMaxRatio(*i3ev);
        //enMaxWeightVect->push_back(eMaxRatio);
        timeWeightVect->push_back(TimeRatio);
        raVect->push_back(i3ev->GetEquatorialDeg().GetRa());
        decVect->push_back(i3ev->GetEquatorialDeg().GetDec());
        angErrVect->push_back(i3ev->GetSigma());
        eneVect->push_back(i3ev->GetParams().energyValue);
        eventID->push_back(i3ev->GetParams().eventID);
        runID->push_back(i3ev->GetParams().runID);
        timeVect->push_back(i3ev->GetTime().GetMJD());
      }
    }//end loop over flares

    lh=lh+1;
    if(lh > 1) nevs++;
    LogLambda += log( lh );

  }//end loop over events

  if(tdLlh->GetUsePrior()) {
    double testRa  = 0.;
    double testDec = 0.;
    if(tdLlh->IsFitSrc()) {
      int ibin  = 4*nFlares;
      testRa    = par[ibin];
      testDec   = par[ibin+1];
    }
    else {
      const EquatorialDeg *srcEqDeg = ( dynamic_cast<const EquatorialDeg*>(tdLlh->GetSearchCoord()) );
      testRa    = srcEqDeg->GetRa();
      testDec   = srcEqDeg->GetDec();
    }
    LogLambda += log( tdLlh->calcSpatialPrior(testRa,testDec));
  }

  if (tdLlh->GetOptimizeTolerance() > 0.) { // correction for events skipped in optimization
    lh=0;
    for(int iflare=0; iflare<nFlares; iflare++) lh += (-srcFrac[iflare]);
    lh=lh+1;
    LogLambda += (tdLlh->Get_nEvents()-tdLlh->GetNEventsSelected())*log(lh); 
  }

  double margWeight = log(2);
  for(int iflare=0; iflare<nFlares; iflare++) {
    margWeight +=  log(sigma[iflare] * norm[iflare] / ( timePdf[iflare]->GetLivetime() ) );
 
    // Jim's Spectral index penalty
    // We could expect a gamma between 2.0 and 2.7
    // he applies a Gaussian with sigma 0.2 units in spectrum outside of that range
    if (tdLlh->SpectralPenalty_ && tdLlh->GetUseEnergy()) {
      if (gamma[iflare] < 2.) LogLambda -= (2.-gamma[iflare])*(2.-gamma[iflare])/0.08;
      else if (gamma[iflare] > 2.7) LogLambda -= (2.7-gamma[iflare])*(2.7-gamma[iflare])/0.08;
    }
    if ( tdLlh->GetMonitorLevel()>1 ) { // more useful monitoring, with only
                                    // one printout per FCN call. Still lots of spew
                                    // but it's pretty managable for one or a few trials.
      printf("LogLambda=%12.6lg : nSrc=%9.4lg : gamma=%5.4lg : par1=%8.8lg : par2=%8.4lg : nEvs=%i : fact=%g\n",
           LogLambda, srcFrac[iflare]*tdLlh->Get_nEvents(), gamma[iflare], mean[iflare], sigma[iflare], tdLlh->Get_nevs(), tdLlh->Get_margWeight() );
    }
  }

  if ( tdLlh->JimsTerm_ && margWeight < 0) LogLambda += margWeight;

  tdLlh->Set_nevs(nevs);
  tdLlh->SetMargValue(margWeight);
  f = - LogLambda;    // What Minuit minimizes: -log L
  for(int iflare=0; iflare<nFlares; iflare++){
    if (timePdf[iflare]) delete timePdf[iflare];  // otherwise we'd get a memory leak
  }

  return;
}


NewLlhGausTime::NewLlhGausTime() :
  i3Set_(NULL),
  useLCBkgProb_(false),
  fitSrcCoord_(false),
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
  nFlareGuess_(1),
  sigmamin_(10000),
  sigmamax_(-999),
  isSetLowLimSigma_(false),
  isSetUpLimSigma_(false),
  TSthr_(0),
  logLambdaBest_(0.),
  nSigmaTrunc_(4),
  nSrcBest_(0.),
  gammaBest_(0.),
  chiSq_(0.),
  chiSqProb_(0.),
  nEventsTot_(0),
  monitorLevel_(0),
  optimizeTolerance_(0.),
  optimizeAngleDeg_(0.),
  timePdfType_(0),
  doInitGauss_(true),
  optUsePrior_(false),
  setNormFromGRL_(false),
  isGaussTrunc_(false),
  priorCoord_(NULL),
  sigmaX_(1.0),
  sigmaY_(1.0),
  theta_(0.),
  srcCoordUncRA_(0.),
  srcCoordUncDEC_(0.),
  stepSizeSrcCoord_(0.1),
  fitfrac_(0.001),
  seedWtMin(1000)
{ 
  tdLlh = this; //tdLlh defined at the bottom of NewLlhGauss.h

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
  // These start when created, so stop immediately
  stopwatch_MaximizeLlh_.Stop();
  stopwatch_optimize_.Stop();
  stopwatch_minuitMigrad_.Stop();
  stopwatch_minuitFCN_.Stop();
}


void NewLlhGausTime::OptimizeEventSelection() {
  const EquatorialDeg *srcEqDeg = ( dynamic_cast<const EquatorialDeg*>(srcCoord_) );

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

  double threshold = (1./srcFracMax_ - 1.) * ( exp(optimizeTolerance_/nEventsTot_) -1. );
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
        
    double lcProb;
    if ( useLCBkgProb_ ) {
        if (UseFolded_){
            lcProb = lcBkgProb_->BackgroundLCProb_folded(*event);
        }  else {
            lcProb = lcBkgProb_->BackgroundLCProb(*event); 
        }
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
      eProbVect_.push_back(eProb);
      spaceRatioVect_.push_back(spaceRatio);
      bkgSpaceProbVect_.push_back(bkgSpaceProb);
      lcBkgProbVect_.push_back(lcProb);
      sigmaSpaceVect_.push_back(event->GetSigma());
      timeVect_.push_back(event->GetTime().GetMJD());

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



void NewLlhGausTime::PrepareAnalysis() {
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
  bkgSpaceProbVect_.clear();
  lcBkgProbVect_.clear();
  sigmaSpaceVect_.clear();

  eventRatioVect_.clear();
  spaceWeightVect_.clear();
  enWeightVect_.clear();
  //enMaxWeightVect_.clear();
  timeWeightVect_.clear();
  raVect_.clear();
  decVect_.clear();
  angErrVect_.clear();
  eneVect_.clear();
  timeVect_.clear();
  eventID_.clear();

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
      double bkgSpaceProb = event->GetBkgSpaceProbFn()->GetBkgProbDensity(*event);
      //    if event is in region with b=0 (what to do?)
      if (bkgSpaceProb <=0) { log_fatal("bkgSpaceProb <= 0\n"); }
      double spaceRatio = sigSpaceProb / bkgSpaceProb;
      const EnergyProb* eProb(NULL);
      if (optUseEnergy_) { 
	eProb = event->GetEnergyProbFn();
      }
      eVect_.push_back(*event);
      selectedList_.AddEvent(event);
      eProbVect_.push_back(eProb);
      spaceRatioVect_.push_back(spaceRatio);
      bkgSpaceProbVect_.push_back(bkgSpaceProb);
      sigmaSpaceVect_.push_back(event->GetSigma());
    }
  }
  
  nEventsSelected_ = eVect_.size();
}


void NewLlhGausTime::MaximizeLlh()
{
  stopwatch_MaximizeLlh_.Start(false);  //  false = don't reset

  timePdfType_=1; //I should get rid of this, but set it here for now...

  PrepareAnalysis();
  // TO DO: BETTER HANDLING IF NO EVENTS SELECTED:
  assert(eVect_.size());
  //assert(selectedList_.GetSize());

  if(doInitGauss_) GetFlareGuess(optUseEnergy_,1, nsrcGuess_, gammaGuess_, meanGuess_, sigmaGuess_);
  nFlareGuess_ = nsrcGuess_.size();  

  TMinuit *minuit = new TMinuit(4*nFlareGuess_);
  Double_t arglist[1];
  arglist[0] = 0.5; //0.5 for likelihood fit
  Int_t ierflg=0;
  minuit->SetFCN(llhFnc);
  minuit->mnexcm("SET ERR",arglist,1,ierflg);
  minuit->SetPrintLevel(-1);

  for(int i=0; i<nFlareGuess_; i++) {
    if(gammaGuess_[i]<1.) gammaGuess_[i] = 1.1;
    if(gammaGuess_[i]>4.) gammaGuess_[i] = 3.9;
    if (!nsrcGuess_[i]) { nsrcGuess_[i]=0.5; }
    if(!isSetUpLimSigma_) sigmamax_ = (tmax_-tmin_)/2.;
    if(sigmaGuess_[i]>=sigmamax_) sigmaGuess_[i] = sigmamax_-0.1;
    if(sigmaGuess_[i]<=sigmamin_) sigmaGuess_[i] = sigmamin_+0.1;
  }

  // sigmamin_ is set up to be the minimimum witdh of the
  // time Gaussian. If you let the width become arbitrarily small
  // the method could find that it can keep shrinking onto one strong event
  // and keep getting a better test statistic, which isn't what we want.
  // Sets a lower bound of the rms between the closest two events in time.

  if(!isSetLowLimSigma_){
    sigmamin_ = 1000.;
    double timed;
    int sp = 1;

    vector<double> ttt;
    for (int i=0;i<nEventsSelected_;i++)  ttt.push_back( eVect_[i].GetMJD() );
    sort(ttt.begin(),ttt.end());
  
    for (unsigned int i=sp; i<ttt.size(); i++) {//This assumes events are time-ordered.
      timed = ttt.at(i)-ttt.at(i-sp);
      if (monitorLevel_ > 2)  cout << ttt.at(i) << " "; 
      if ( (fabs(timed) < sigmamin_) && (fabs(timed)>0) ) sigmamin_ = fabs(timed);
    }
    sigmamin_ /= sqrt(2.);
    if (monitorLevel_ > 2) { cout << " sigmamin= " << sigmamin_ << endl; }
    ttt.clear();
  }
  
  if(sigmamin_<1E-8) sigmamin_   = 1E-8;

  ///////  ns ///////
  double nSrcMax = srcFracMax_*nEventsTot_;
  nSrcMin_ = 0;
  for(int i=0; i<nFlareGuess_; i++)  {
    minuit->mnparm(0+i*4, Form("nSrc_f%d",i), nsrcGuess_[i], 0.1, nSrcMin_, nSrcMax, ierflg);
    if (!optParAuto_[0])  minuit->FixParameter(0+i*4);
  }
    
  ///////  gamma ///////
  gammaMin_ = eProbVect_[0]->GetGammaMin();
  gammaMax_ = eProbVect_[0]->GetGammaMax();
  gammaMin_=1.;
  gammaMax_=4.;
  for(int i=0; i<nFlareGuess_; i++)  {
    minuit->mnparm(1+i*4, Form("gamma_f%d",i), gammaGuess_[i], 0.1, gammaMin_, gammaMax_, ierflg);
    if (!optParAuto_[1])  minuit->FixParameter(1+i*4);
  }
 
  ///////  mean time ///////
  for(int i=0; i<nFlareGuess_; i++)  {
    double  initStepSize = sigmaGuess_[i];
    double  initValue = meanGuess_[i];
    double lowLimit;
    if(meanGuess_[i]-3*sigmaGuess_[i] > tmin_) lowLimit = meanGuess_[i]-3*sigmaGuess_[i];
    else lowLimit = tmin_;
    double upLimit;
    if(meanGuess_[i]+3*sigmaGuess_[i] < tmax_) upLimit = meanGuess_[i]+3*sigmaGuess_[i];
    else upLimit = tmax_;
    minuit->mnparm(2+i*4, Form("mean_f%d",i), initValue, initStepSize, lowLimit, upLimit, ierflg);
    if (!optParAuto_[2]) minuit->FixParameter(2+i*4);
  }
   
  ///////  sigma time ///////
  // Using the log10 of sigma as the fit parameter to make things more smooth.
  double lowLimitS = log10(sigmamin_);
  double upLimitS = log10(sigmamax_);
  double initStepSizeS = 0.1;
  for(int i=0; i<nFlareGuess_; i++)  {
    double initValueS = log10(sigmaGuess_[i]);
    minuit->mnparm(3+i*4, Form("sigma_f%d",i), initValueS, initStepSizeS, lowLimitS, upLimitS, ierflg);
    if (!optParAuto_[3]) minuit->FixParameter(3+i*4);
  }
 
  ///////  source RA, DEC ///////
  if (optParAuto_[4] && optParAuto_[5]) {
    // here we fit src position coordinates
    double initValueRA  = srcCoordGuess_->GetRa();
    double initValueDEC = srcCoordGuess_->GetDec();
    double lowLimitRA   = initValueRA  - srcCoordUncRA_; 
    double lowLimitDEC  = initValueDEC - srcCoordUncDEC_; 
    double upLimitRA    = initValueRA  + srcCoordUncRA_;
    double upLimitDEC   = initValueDEC + srcCoordUncDEC_;
    minuit->mnparm(4+nFlareGuess_*4, "Ra",  initValueRA,  stepSizeSrcCoord_, lowLimitRA,  upLimitRA, ierflg);
    minuit->mnparm(5+nFlareGuess_*4, "Dec", initValueDEC, stepSizeSrcCoord_, lowLimitDEC, upLimitDEC, ierflg);
  }
  
  // MINIMIZE

  stopwatch_minuitMigrad_.Start(false);  //  false = don't reset
  double arglist2[2];
  arglist2[0] = 500;
  arglist2[1] = 0.1;
  minuit->mnexcm("MIGRAD", arglist2 ,2,minuitOut_);
  
  stopwatch_minuitMigrad_.Stop();

  // GET RESULTS
  if(minuitOut_!=0) {
    int itrial = 0;
    double tmpnsrcGuess  = 0.;
    double tmpgammaGuess = 0.;
    double tmpsigmaGuess = 0.;
    //try to change init ns value
    while(minuitOut_!=0 && itrial<100) {
      for(int i=0; i<nFlareGuess_; i++)  {
        tmpnsrcGuess = random_uniform(0.,2*nsrcGuess_[i]);
        minuit->mnparm(0+i*4, Form("nSrc_f%d",i), tmpnsrcGuess, 0.1, nSrcMin_, nSrcMax, ierflg);
      }
      minuit->mnexcm("MIGRAD", arglist2 ,2,minuitOut_);
      itrial++;
    }
    itrial = 0;
    if(minuitOut_!=0) {
      //reset init ns value to nsrcGuess_
      for(int i=0; i<nFlareGuess_; i++)  minuit->mnparm(0+i*4, Form("nSrc_f%d",i), nsrcGuess_[i], 0.1, nSrcMin_, nSrcMax, ierflg);
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
  if(minuitOut_==0) {
    StoreLogLambdaBest(minuit);   // Set logLambdaBest_
    double par=0., err=0.;
    ClearBestFitParams();
    for(int iflare=0; iflare<nFlareGuess_; iflare++) {
      minuit->GetParameter(0+iflare*4,par,err);
      nSrcBest_.push_back(make_pair(par,err));
      minuit->GetParameter(1+iflare*4,par,err);
      gammaBest_.push_back(make_pair(par,err));
      minuit->GetParameter(2+iflare*4,par,err);
      meanBest_.push_back(make_pair(par,err));
      minuit->GetParameter(3+iflare*4,par,err);
      sigmaBest_.push_back(make_pair(pow(10.,par),pow(10.,err)));
      if(fitSrcCoord_) {
        minuit->GetParameter(4+nFlareGuess_*4,par,err);
        raSrcBest_.push_back(make_pair(par,err));
        minuit->GetParameter(5+nFlareGuess_*4,par,err);
        decSrcBest_.push_back(make_pair(par,err));
      }
    }
    
    chiSq_ = 2.*logLambdaBest_;
    double p_temp;
    chisq_prob(chiSq_, ndof_, &p_temp, &chiSqProb_);
    estProb_ = chiSqProb_ / 2.;  // one-sided chi-sq prob
  }
  stopwatch_MaximizeLlh_.Stop();

  delete minuit;
}

///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  END MAXIMIZELLH()  ///////////////////////
///////////////////////////////////////////////////////////////////////////////

// This is clumsy but seems like the only way to get this info
void NewLlhGausTime::StoreLogLambdaBest(TMinuit *minuit) 
{
  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, istat;
  minuit->mnstat(amin, edm, errdef, nvpar, nparx,istat);  
  logLambdaBest_ = -amin;  // that is, max llh = - (minimizer result)

  // The *worst* logLambdaBest should be zero (i.e. null hypothesis).
  // But minimizer will miss exact zero, leading logLambdaBest_ to be slightly
  // negative. Since this will cause probability calculation to choke, we 
  // fix it here
  if (logLambdaBest_ < 0.) { logLambdaBest_ = 0.;}
}

double NewLlhGausTime::calcSpatialPrior(double testRA, double testDEC) {

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

// This function scans over all the events which could have an S/B ratio > 10,
// and tests consecutive doubles, triples... for compatability with and E^-2 flare.
// I tried testing for Gamma = 2 and like 3.5 before, but it doesn't seem to do much.
// This is configurable enough to do whatever you'd like.
// There is a final scan over Gamma at the end to round out the set parameters.

// Oh, I'm using this parameter close_ to set how far into the vector 1,2,3,4,5,10,15,20,25
// to test for flare compatability. You only need to go to 10, I usually go to 15 (close_=6).

void NewLlhGausTime::GetFlareGuess(bool useE, int useSpace, vector<double> & Guess_nsrc, vector<double> & Guess_gamma, vector<double> & Guess_mean, vector<double> & Guess_sigma) 
{
  nFlareGuess_ = 1;
  Guess_mean.clear();
  Guess_sigma.clear();
  Guess_nsrc.clear();
  Guess_gamma.clear();

  vector<double> tVectortemp;
  vector<double> tVectorclose;

  double llhtemp, rms, avgtime;
  for (int j=0; j<nEventsSelected_; j++) {      
    if ( (eProbVect_[j]->GetEnergyMaxRatio(eVect_[j])*spaceRatioVect_[j]) > seedWtMin ) {
        tVectorclose.push_back(eVect_[j].GetMJD());
    }
  }
  
  sort(tVectorclose.begin(),tVectorclose.end());

  if (monitorLevel_ > 0) { cout << tVectorclose.size() << " events with S/B > " << seedWtMin << ". Sigmamin is " << sigmamin_ << endl; }


//  const int ngamma=2;
//  double g[ngamma] = {2.0, 3.5};
  const int ngamma=1;
  double g[ngamma] = {2.0};

  int kmax = (int) close_;
  if (11 < kmax) { kmax = 11; }

  int n[] =  {1,2,3,4,5,10,15,20,25,30,35,40};
  int ns[] = {1,2,3,4,5, 5, 5, 5, 5, 5, 5, 5};

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
  //for (int k=kmin; k<kmax; k++){ //k==1 for pairs, k==2 for triples, etc.
  for(int k=kmax-1; k>=kmin; k--){
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
      //cout << " GUESS k " << k << " i " << i << " avgtime " << avgtime << " rms " << rms << endl;

      if (useE) {
        for (int h=0; h<ngamma; h++) {
          llhtemp = EvaluateLlh( ns[k]+1., g[h], avgtime, rms );
          if (llhtemp > llhMax) {
            llhMax = llhtemp;
            meanBckup = avgtime;
            rmsBckup = rms;
            nsrcBckup = ns[k]+1.;
            gammaBckup = g[h];
          }
          if(llhtemp>TSthr_) {
            nsVec.push_back(ns[k]+1.);  
            gammaVec.push_back(g[h]);  
            t0Vec.push_back(avgtime);  
            sigVec.push_back(rms);  
            TSVec.push_back(make_pair(llhtemp,ind));
            ind++;
          }
        }
        tVectortemp.clear();
      } else { // I have no idea why you'd want to do the analysis without energy, but here you go...
        llhtemp = EvaluateLlh( ns[k]+1., 0., avgtime, rms );
        if (llhtemp > llhMax) {
          llhMax = llhtemp;
          meanBckup = avgtime;
          rmsBckup = rms;
          nsrcBckup = ns[k]+1.;
          gammaBckup = 0.;
        }
        if(llhtemp>TSthr_) {
          nsVec.push_back(ns[k]+1.);
          gammaVec.push_back(0);
          t0Vec.push_back(avgtime);
          sigVec.push_back(rms);
          TSVec.push_back(make_pair(llhtemp,ind));
          ind++;
        }
        tVectortemp.clear();
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

  tVectorclose.clear();
  
  sort(TSVec.begin(), TSVec.end(), sortinrev);
  tmpTSVec.clear();
  tmpNsVec.clear();
  tmpGammaVec.clear();
  tmpSigVec.clear();
  tmpT0Vec.clear();

  //now select only non-overlapping flares
  Guess_mean.push_back(t0Vec[TSVec[0].second]);                                                                                                                            
  Guess_sigma.push_back(sigVec[TSVec[0].second]);     
  Guess_nsrc.push_back(nsVec[TSVec[0].second]);                                                                                                                                    
  Guess_gamma.push_back(gammaVec[TSVec[0].second]);                                                                                                                                  
  Guess_TS.push_back(TSVec[0].first);                            

  double t0 = t0Vec[TSVec[0].second];
  double sigma = sigVec[TSVec[0].second];
  double tmin = t0 - nSigmaTrunc_*sigma;
  double tmax = t0 + nSigmaTrunc_*sigma;

  tminVec.push_back(tmin);
  tmaxVec.push_back(tmax);

  bool isOverlap=false;
  for(unsigned int i=1; i<t0Vec.size(); i++) {
    t0  = t0Vec[TSVec[i].second];
    sigma = sigVec[TSVec[i].second];
    tmin = t0 - nSigmaTrunc_*sigma;
    tmax = t0 + nSigmaTrunc_*sigma;
    for(unsigned int j=0; j<tminVec.size(); j++) {
      if((tmin >= tminVec[j] && tmax <= tmaxVec[j]) ||                           // 1  2  2  1
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
      Guess_sigma.push_back(sigVec[TSVec[i].second]);     
      Guess_nsrc.push_back(nsVec[TSVec[i].second]);                                                                                                                                    
      Guess_gamma.push_back(gammaVec[TSVec[i].second]);                                                                                                                                  
      Guess_TS.push_back(TSVec[i].first);                            
      tminVec.push_back(tmin);
      tmaxVec.push_back(tmax);
    }
    else isOverlap=false;
  }
 
  ind=0;

  //let's do a finer scan through gamma values and save results in vectors
  for(unsigned int i=0; i<tmpTSVec.size(); i++) {
    llhMax=0.;
    for (double d=1.; d<4.; d+=0.2) { // loop over gamma with 15 steps for best seed                                                                   
      llhtemp = EvaluateLlh( nsVec[i], d, t0Vec[i], sigVec[i] );
      if(llhtemp>llhMax) {
        llhMax = llhtemp;
        Guess_gamma[i] = d;
      }
    }
  }
  
  if (useE || useSpace) { } //nip compiler complaints
 
}


void NewLlhGausTime::SetNullTestStat(TH1D * inputhisto) {
  
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

double NewLlhGausTime::GetProbFromHisto(double teststat){ //, bool useFit) {

  // this looks at a null teststat distribution
  // which we loaded, and if it's in the exponential
  // tail, should go with a fit expo->Eval(teststat) (top 1%)
  
  int bin = pvalHisto->FindBin(teststat);
  double ptemp = pvalHisto->GetBinContent(bin);
  
//  if (ptemp < fitfrac_ && useFit) { //small enough that we should do the fit stuff
//    ptemp = fitfn_->Eval(teststat);
//  }
  return ptemp;
  
  // in case we don't feel like using a fit and it is higher than the distribution
  // just return 1/trials (conservative?)
  if( !pvalHisto->GetBinContent(bin) ) { return 1.0/nullTestStat_->GetEntries(); }
  
  return -1.; //if we get here something is wrong, 
              // and may as well hear about it trying to take the log
}


double NewLlhGausTime::EvalFCN(const vector<double>& parVect) const {

  assert(nEventsTot_);
  double f;
  const int npar = parVect.size();
  double gin[npar];
  double par[npar];
  int iflag = 0;
  for(int i=0; i<npar; i++) par[i] = parVect[i];
  int npars=npar;
  llhFnc(npars, gin, f, par, iflag);

  return f;
}

double NewLlhGausTime::EvaluateLlh(double nSrc, double gamma, double mean, double sigma) {
  vector<double> parVect;
  parVect.push_back(nSrc);
  parVect.push_back(gamma);
  parVect.push_back(mean);
  parVect.push_back( log10(sigma) );
  bool resetUseFitSrc=false;
  if(IsFitSrc()) {
    SetUseFitSrc(false);
    resetUseFitSrc=true;
  }
  double minusLlh = EvalFCN(parVect);
  if(resetUseFitSrc) SetUseFitSrc(true); //this step is necessary when using the llh to fit also the src position.
                                         //In fact, EvaluateLlh is called by GetFlareGauss, which is used to initialized the other llh params;
                                         //The position can be initialized manually by doing a scan of the llh over the interested region;
                                         //hence its initialization does not enter in GetFlareGauss. For this reason, we keep fixed the position when
                                         //initializing the other parameters.
  return -minusLlh;   // that is, max llh = - (minimizer result)
}

// For MultiAnalysis: Translate Global Par values to Individual Par Values:
const vector<double> NewLlhGausTime_ParTranslator::Translate(int llhIndex, const vector<double>& parGlobal) const {
    double nSrcGlobal = parGlobal[0];
    double gamma = parGlobal[1];
    vector<double> parLocal(4);
    double weight =  LinearInterpolate(gamma, gammaMin_, gammaMax_, srcWeightVect_[llhIndex]);
    // recall: srcWeightVect_[i] is itself a vector<double> 
    parLocal[0] = nSrcGlobal * weight;
    parLocal[1] = gamma;
    parLocal[2] = parGlobal[2];
    parLocal[3] = parGlobal[3];
    return parLocal;
}

const vector<double> NewLlhGausTime_ParTranslator::Translate(int llhIndex, const vector<double>& parGlobal, vector<vector<double> > weights) const {
  int nFlares = weights.size();
  vector<double> parLocal;
  for(int iflare=0; iflare<nFlares; iflare++) {
    int ibin = iflare*4;
    double nSrcGlobal = parGlobal[0+ibin];
    double gamma = parGlobal[1+ibin];
    // recall: srcWeightVect_[i] is itself a vector<vector<double> > 
    parLocal.push_back(nSrcGlobal * weights[iflare][llhIndex]);
    parLocal.push_back(gamma);
    parLocal.push_back(parGlobal[2+ibin]);
    parLocal.push_back(parGlobal[3+ibin]);
  }
  return parLocal;
}

vector<vector<double> > NewLlhGausTime_ParTranslator::MakeWeights(vector<NewLlhGausTime*> llhVect, const vector<double>& parGlobal, int nFlares) const {
  vector<vector<double> > weightVect;
  for(int iflare=0; iflare<nFlares; iflare++) {
    int ibin     = iflare*4;
    double gamma = parGlobal[1+ibin];
    double tmean = parGlobal[2+ibin];
    double sigma = pow(10, parGlobal[3+ibin]);

    double zlow, zhigh;
    double weightSum = 0;
    vector<double> tmpWeightVect;
    
    for(int llhIndex = 0; llhIndex<int(llhVect.size()); llhIndex++){
      double weightGamma =  LinearInterpolate(gamma, gammaMin_, gammaMax_, srcWeightVect_[llhIndex]);

      double tmin=0, tmax=0;

      if(llhVect[llhIndex]->GetAnalysisSet()->GetSource()->GetTimePdf()){
        tmin = llhVect[llhIndex]->GetAnalysisSet()->GetSource()->GetTimePdf()->GetTmin();
        tmax = llhVect[llhIndex]->GetAnalysisSet()->GetSource()->GetTimePdf()->GetTmax();
      }
      else if(llhVect[llhIndex]->GetAnalysisSet()->GetSource()->GetTimePdfVect().at(0)){
        tmin = llhVect[llhIndex]->GetAnalysisSet()->GetSource()->GetTimePdfVect().at(0)->GetTmin();
        tmax = llhVect[llhIndex]->GetAnalysisSet()->GetSource()->GetTimePdfVect().at(0)->GetTmax();
      }
      else log_error("ParTranslator: no source time pdf.\n"); 

      zlow  = TMath::Erf( (tmin-tmean)/(sigma*TMath::Sqrt(2.)) );
      zhigh = TMath::Erf( (tmax-tmean)/(sigma*TMath::Sqrt(2.)) );
      double weightTmp = (zhigh-zlow)/2.;
      tmpWeightVect.push_back( weightGamma*weightTmp );
      weightSum += weightGamma*weightTmp;
    }
    for(unsigned int i=0; i<tmpWeightVect.size(); ++i){ tmpWeightVect[i] /= weightSum; }
    weightVect.push_back(tmpWeightVect);
  }
  return weightVect;
}

void NewLlhGausTime_ParTranslator::SetTranslator(const vector<AnalysisSet*>& aSetVect) {
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
      double nev = aSetVect[i]->GetMeanSrcNevForFluxModel(PowerLawFlux(1,-gamma));
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
