#include "llhTimeDep/public/NewLlhBlockTime.h"
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

  // Welcome to NewLlhBlockTime.C!
  
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

void llhBlockFnc(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  double srcFrac = par[0]/tdLlhBlock->Get_nEvents();   // relative weight of src term 

  vector<I3Event> eVect = tdLlhBlock->GetEventVector();
  vector<double> eventRatioVect = tdLlhBlock->GetEventRatios();
  EventPtrList selectedList = tdLlhBlock->GetSelectedList();
  vector<const EnergyProb*> eProbVect = tdLlhBlock->GetEProbVect();
  vector<double> spaceRatioVect = tdLlhBlock->GetSpaceRatioVect();

  if (srcFrac > tdLlhBlock->GetOptimizeSrcFracMax()) {
    log_warn("Trying to evaluate srcFrac=%f which is > srcFracMax_=%f",
             srcFrac, tdLlhBlock->GetOptimizeSrcFracMax());
    log_warn("Check your parameter settings and make them consistent.");
    log_warn("(If running optimized, logLambda tolerance could be exceeded.)");
  }

  if (tdLlhBlock->GetOptStoreRatios()) { eventRatioVect.clear(); }

  double gamma = par[1];  // Spectral index

  if (tdLlhBlock->GetUseEnergy() && (gamma > 4. || gamma < 1. || TMath::IsNaN(gamma) ) ) {
    f = 1e50;
    return;
  }

  double LogLambda = 0.;
  double TimeRatio = 1.;

  //  for blocks

  if ( tdLlhBlock->GetMonitorLevel() > 4 ) {
    cout << "timePdf->SetBlockLevels(" << tdLlhBlock->GetBlocksFile() << ")" << endl;
  }

  double lag = par[2];
  double threshold = par[3];
  //if (lag) { }
  TimePdf * timePdf = new BlockTimePdf1();
  timePdf->SetBlockLevels(tdLlhBlock->GetBlocksFile(), threshold);
                                        //this loads the blocks 
                                        //from a text file    
  timePdf->SetLivetime( tdLlhBlock->GetLivetime() );
  //vector<double> trange = tdLlhBlock->i3Set_->GetEventTimeModulePtr()->GetTimeRange();
  timePdf->CheckTimeBounds(tdLlhBlock->tmin_, tdLlhBlock->tmax_);

  tdLlhBlock->Set_nevs(0); // number of events which contribute S/B > 1

  double lh;

  for (int i=0;i<tdLlhBlock->Get_nEventsSelected();i++)  { // loop over selected events

    const Event* event = selectedList.GetEvent(i);


    TimeRatio = timePdf->GetPdfValue( event->GetTime().GetMJD(), lag ) * tdLlhBlock->GetLivetime();
    //    TimeRatio = timePdf->GetPdfValue( tdLlhBlock->eVect_[i].GetTime().GetMJD(), lag ) * tdLlhBlock->livetime_;

    // Perhaps we can speed things up... if the TimeRatio==0, skip the rest of the calculations.
    // include an out for high monitoring levels...
    if (tdLlhBlock->GetMonitorLevel() < 4 && TimeRatio < 1e-10) {
      LogLambda += log(1.-srcFrac);
      continue;
    }

    // first comes the event's energy term
    double eRatio=1.;
    if (tdLlhBlock->GetUseEnergy()) {
      eRatio = eProbVect[i]->GetEnergyProbGamma(*event, gamma) /
                   eProbVect[i]->GetEnergyProbBkg(*event);
    }

    //As implemented, local coord background terms will automagically be included in spaceRatoVect   
    lh = ( srcFrac * ( spaceRatioVect[i] * eRatio * TimeRatio - 1. ) + 1. );

    if(lh > 1) {tdLlhBlock->Increment_nevs(1);}
    LogLambda += log( lh );

    if ( (tdLlhBlock->GetMonitorLevel() > 3 && lh > 1.) || tdLlhBlock->GetMonitorLevel() > 4 ) {
                                      // useful monitoring but with massive event-by-event spew
                                      // (~200 lines) going through each FCN call (~100-300)
      cout << lh << " " << LogLambda << " " << srcFrac << " " << spaceRatioVect[i] << " " <<
              eRatio << " " << eVect[i].GetTime().GetMJD() << " " <<
              timePdf->GetPdfValue(eVect[i].GetTime().GetMJD(),lag) << " " << tdLlhBlock->GetLivetime() << endl;
    }

    if (tdLlhBlock->GetOptStoreRatios()) {
      eventRatioVect.push_back(spaceRatioVect[i]*eRatio*TimeRatio);
    }

  } //end loop over events

  if (tdLlhBlock->GetOptimizeTolerance() > 0.) { // correction for events skipped in optimization
    LogLambda += (tdLlhBlock->Get_nEvents()-tdLlhBlock->Get_nEventsSelected())*log(1.-srcFrac);
  }

  //if (tdLlhBlock->JimsTerm_) { 
  //  LogLambda += log( BlockTimeAboveThresh(tdLlhBlock->BlocksTimeFile_, threshold)/timePdf->GetLivetime() );
  //}
  //tdLlhBlock->margValue_ = log( BlockTimeAboveThresh(tdLlhBlock->BlocksTimeFile_, threshold)/timePdf->GetLivetime() );
  tdLlhBlock->SetMargValue(1.);

  // Jim's Spectral index penalty
  // We could expect a gamma between 2.0 and 2.7
  // he applies a Gaussian with sigma 0.2 units in spectrum outside of that range
  if (tdLlhBlock->GetSpectralPenalty() && tdLlhBlock->GetUseEnergy()) {
    if (gamma < 2.) LogLambda -= (2.-gamma)*(2.-gamma)/0.08;
    else if (gamma > 2.7) LogLambda -= (2.7-gamma)*(2.7-gamma)/0.08;
  }

  f = - LogLambda;    // What Minuit minimizes: -log L

  if ( tdLlhBlock->GetMonitorLevel()>1 ) { // more useful monitoring, with only
                           // one printout per FCN call. Still lots of spew
                           // but it's pretty managable for one or a few trials.
    printf("LogLambda=%12.6lg : nSrc=%9.4lg : gamma=%5.4lg : par1=%8.3lg : par2=%8.4lg : nEvs=%i : fact=%g\n",
           LogLambda, srcFrac*tdLlhBlock->Get_nEvents(), gamma, lag, threshold, tdLlhBlock->Get_nevs(), tdLlhBlock->Get_margValue() );
  }

  if (timePdf) { delete timePdf; } // otherwise we'd get a memory leak

  return;
}

NewLlhBlockTime::NewLlhBlockTime() :
  i3Set_(NULL),
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
  nEventsTot_(0),
  monitorLevel_(0),
  optimizeTolerance_(0.),
  optimizeAngleDeg_(0.),
  timePdfType_(0),
  fitfrac_(0.001)
{
  blocksMinThreshold_=0.;
  
  tdLlhBlock = this; //tdLlh defined at the bottom of NewLlhBlock.h
 
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


void NewLlhBlockTime::OptimizeEventSelection() {
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
  const EventPtrList* srcList = aSet_->GetEventPtrListSrc();

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
      if (monitorLevel_ > 2) {
        cout << "(" << eventRaDeg << "," << eventDecDeg << ")  Sr : " << sigSpaceProb << "/" << bkgSpaceProb << " eP : " << eMaxRatio << endl;
      }
    }
  }

  int SrcEventsTot_ = srcList->GetSize();


  for (int i=0; i<SrcEventsTot_; ++i) {

    const I3Event* srcEvent = (dynamic_cast<const I3Event*> (srcList->GetEvent(i)));
    if (!srcEvent) { log_fatal("OptimizeEventSelection: failed cast to event.\n");}

    srcVect_.push_back(*srcEvent);

    }


  if (monitorLevel_ > 0) {
    printf("Optimizing: %d events selected with maxRatio >%lg"
	   "  out of %d Ntotal.\n",
	   selectedList_.GetSize(), threshold, nEventsTot_);
  }
}



void NewLlhBlockTime::PrepareAnalysis() {
  if (!aSet_) { log_fatal("PrepareAnalysis: AnalysisSet was not set.\n"); }
  if (!srcCoord_) { log_fatal("PrepareAnalysis: srcCoord was not set.\n"); }

  const EventPtrList* evList = aSet_->GetEventPtrList();

  const EventPtrList* srcList = aSet_->GetEventPtrListSrc();

  if (!evList) { log_fatal("PrepareAnalysis: EventPtrList was not set.\n"); }
  nEventsTot_ = evList->GetSize();
  if (!nEventsTot_) { log_fatal("PrepareAnalysis: EventPtrList was empty.\n"); }
  // or, is there a reason to allow zero events?

  selectedList_.Clear();
  eVect_.clear();
  srcVect_.clear();
  eProbVect_.clear();
  spaceRatioVect_.clear();

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
      eVect_.push_back(*event);
      selectedList_.AddEvent(event);
      eProbVect_.push_back(eProb);
      spaceRatioVect_.push_back(spaceRatio);
    }

    int SrcEventsTot_ = srcList->GetSize();

    cout << "srceventstot no optimization" << SrcEventsTot_ << endl;

    for (int i=0; i<SrcEventsTot_; ++i) {

        const I3Event* srcEvent = (dynamic_cast<const I3Event*> (srcList->GetEvent(i)));

        if (!srcEvent) { log_fatal("OptimizeEventSelection: failed cast to event.\n");}

        srcVect_.push_back(*srcEvent);

    }

  }
  
  nEventsSelected_ = eVect_.size();
  
}


void NewLlhBlockTime::MaximizeLlh()
{
  stopwatch_MaximizeLlh_.Start(false);  //  false = don't reset

  timePdfType_=33; //I should get rid of this, but set it here for now...
  cout << "before prepare analysis " << endl;
  PrepareAnalysis();
  cout << "after prepare analysis " << endl;
  // TO DO: BETTER HANDLING IF NO EVENTS SELECTED:
  assert(eVect_.size());
  assert(selectedList_.GetSize());

  // Wipes out parameters (minuit crashes easily when changing existing pars)
  minuit_ = new TMinuit(4);
  Double_t arglist[1];
  arglist[0] = 0.5; //0.5 for likelihood fit
  Int_t ierflg=0;
  minuit_->SetFCN(llhBlockFnc);
  minuit_->mnexcm("SET ERR",arglist,1,ierflg);
  minuit_->SetPrintLevel(-1); 

  double bestlag=0.;

  double initValue, initStepSize, lowLimit, upLimit;

  // DEFINE PARAMETERS

  if (optParAuto_[0]) {
    double nSrcMax = srcFracMax_*nEventsTot_;
    nSrcMin_ = 0;
    if (!nsrcGuess_) { nsrcGuess_=2.5; }
    minuit_->mnparm(0, "nSrc", nsrcGuess_, 0.1, nSrcMin_, nSrcMax, ierflg);
  }
 
  if (optParAuto_[1]) {
    if (optUseEnergy_) {
      gammaMin_ = eProbVect_[0]->GetGammaMin();
      gammaMax_ = eProbVect_[0]->GetGammaMax();
      
      double gammaInit = (gammaMin_+gammaMax_)/2.;
      
      //gammaInit = gammaGuess_;
      minuit_->mnparm(1, "gamma", gammaInit, 0.1, gammaMin_, gammaMax_, ierflg);
    } else {
      // If not using energy, then FCN will actually replace the
      // energy prob term internally... so these values don't really matter
      // Just fix the gamma so that Minuit doesn't try to minimize w.r.t. it
      minuit_->mnparm(1, "gamma", 0., 0., 0., 0., ierflg);
    }
  }

  if (optParAuto_[2] ) {
      bestlag = SearchForLag();
      initValue = bestlag;
      initStepSize = 2; //Assuming we're using 1-day blocks something  1<x<2 should be good.
      lowLimit = -1.0*laglimit_;
      upLimit = laglimit_;
      minuit_->mnparm(2, "lag", initValue, initStepSize, lowLimit, upLimit, ierflg);
  }
   
  if (optParAuto_[3] ) {
    double initTV, initTU; //initial and upper for the threshold
    SearchBlockSpace(bestlag, initTV, initTU);
    initValue = initTV;
    initStepSize = initTU/5.;
    lowLimit = blocksMinThreshold_;
    if (lowLimit!=0.) log_warn("Setting the lower limit for threshold to %f \n",lowLimit);
    upLimit =  initTU;
    minuit_->mnparm(3, "threshold", initValue, initStepSize, lowLimit, upLimit, ierflg); 
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
  sigmaBest_ = pow( 10., GetPar(3) );

  /* Handled in StoreLogLambdaBest(), I hope...

  // fix, because minimizer may sometimes give slightly negative value instead
  // of exact zero, which will cause probability calcultion to choke...
  if (logLambdaBest_ < 0.) {
    // we can always do better, since LogLambda(ns=0,gamma) = 0.
    logLambdaBest_ = 0.;
    nSrcBest_ = 0.;
    gammaBest_ = 0.;  // i.e., no fit makes sense, if nSrc=0.
  }
  */

  chiSq_ = 2.*logLambdaBest_;

  double p_temp;
  
  chisq_prob(chiSq_, ndof_, &p_temp, &chiSqProb_);
  
  estProb_ = chiSqProb_ / 2.;  // one-sided chi-sq prob

  stopwatch_MaximizeLlh_.Stop();
}

///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  END MAXZIMIZELLH()  ///////////////////////
///////////////////////////////////////////////////////////////////////////////

// This is clumsy but seems like the only way to get this info
void NewLlhBlockTime::StoreLogLambdaBest() 
{
  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, istat;
  minuit_->mnstat(amin, edm, errdef, nvpar, nparx,istat); 
  logLambdaBest_ = -amin;  // that is, max llh = - (minimizer result)

  // The *worst* logLambdaBest should be zero (i.e. null hypothesis).
  // But minimizer will miss exact zero, leading logLambdaBest_ to be slightly
  // negative. Since this will cause probability calculation to choke, we 
  // fix it here
  if (logLambdaBest_ < 0.) { logLambdaBest_ = 0.;}
}


double NewLlhBlockTime::SearchForLag() {

  double llhMax=-100.;
  double lagb=0., llhtemp;

  for (double d=-1.0*laglimit_; d<laglimit_; d=d+laglimit_/10.){
    llhtemp = EvaluateLlh( 2, 2., d, 0);    
    if (llhtemp > llhMax) {
      llhMax = llhtemp;
      lagb = d;
    }
  }  
  
  return lagb;
}

void NewLlhBlockTime::SearchBlockSpace(double lag, double & initV, double & maxT) {

  maxT = (GetHighestBlock(BlocksTimeFile_.c_str()) + GetSecondHighestBlock(BlocksTimeFile_.c_str()))/2.;

  double step = maxT/20.;
  double llhMax=-100.;
  double llhtemp;
  
  for (double d=0;d<maxT;d+=step) {
    llhtemp = EvaluateLlh( 2.5, 2., lag, d);
    if (llhtemp > llhMax) {
      llhMax = llhtemp;
      initV = d;
    }
  }
  
  if (monitorLevel_ > 1) { cout << maxT << " " << initV << endl; }
    
}


void NewLlhBlockTime::SetNullTestStat(TH1D * inputhisto) {
  
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

double NewLlhBlockTime::GetProbFromHisto(double teststat){ //, bool useFit) {

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


double NewLlhBlockTime::EvalFCN(const vector<double>& parVect) const {
  assert(nEventsTot_);
  double f;
  const int npar = parVect.size();
  double gin[npar];
  double par[npar];
  int iflag = 0;
  for(int i=0; i<npar; i++) par[i] = parVect[i];
  int npars=npar;
  llhBlockFnc(npars, gin, f, par, iflag);

  return f; 
}



double NewLlhBlockTime::EvaluateLlh(double nSrc, double gamma, double mean, double sigma) {
  vector<double> parVect;
  parVect.push_back(nSrc);
  parVect.push_back(gamma);
  parVect.push_back(mean);
  parVect.push_back(sigma);
  double minusLlh = EvalFCN(parVect);
  return -minusLlh;   // that is, max llh = - (minimizer result)
}


// TGraph* NewLlhBlockTime::GetContour(double sigma, int npoints, int pa1, int pa2) {
//   minuit_->SetFCN(MinuitWrapperFCN);
//   SetMinuitWrapperPtr(this);  // point to *our* EvalMinuitFCN

//   minuit_->SetErrorDef(sigma*sigma);
//   TGraph* g = dynamic_cast<TGraph*> (minuit_->Contour(npoints, pa1, pa2));

//   SetMinuitWrapperPtr(NULL);  // reset pointer
//   return g;
// }






// For MultiAnalysis: Translate Global Par values to Individual Par Values:

/*

void NewLlhBlockTime_ParTranslator::
SetUpTranslate(const vector<double>& parGlobal) {
  //double nSrcGlobal = parGlobal[0];
  double gamma = parGlobal[1];
  double thresh = parGlobal[3];
  double lag = parGlobal[2];

  double timeweight[10];
  double nsweight[10];
  //double weight[10];
  
  double timesum=0.;
  localweight_.clear();

  TimePdf * timePdf = new BlockTimePdf1();
  timePdf->SetBlockLevels(BlocksTimeFile_, thresh, lag); 
                                        //this loads the blocks 
                                        //from a text file
                                        
  for (unsigned int i=0;i<dataTimeVect_.size();i++){
    nsweight[i] = LinearInterpolate(gamma, gammaMin_,
                                       gammaMax_, srcWeightVect_[i]);

    timePdf->CheckTimeBounds(dataTimeVect_[i][0],dataTimeVect_[i][1]);

    timeweight[i] = timePdf->GetNorm();
    timesum += timeweight[i]*nsweight[i];
  }

  for (unsigned int i=0;i<dataTimeVect_.size();i++){
    localweight_.push_back(timeweight[i]*nsweight[i]/timesum);
  }

  if (timePdf) { delete timePdf; }

} */


const vector<double> NewLlhBlockTime_ParTranslator::
Translate(int llhIndex, const vector<double>& parGlobal) const {

  // Make sure you do "SetUpTranslate()" to reset the time weights
  // and normalization for each new parGlobal!
  
  // this may be a fraught implementation of this:
  //if(llhIndex==0){ SetUpTranslate(parGlobal); }

  double nSrcGlobal = parGlobal[0];
  double gamma = parGlobal[1];
  double threshold = parGlobal[3];
    
  double weight = 
    LinearInterpolate(gamma, gammaMin_, gammaMax_, srcWeightVect_[llhIndex]);
  // recall: srcWeightVect_[i] is itself a vector<double> 
  
  //cout << "Translate! " << threshold << " " << threshMin_ << " " << threshMax_ << " " << blockWeightVect_[llhIndex].size() << " " << srcWeightVect_[llhIndex].size() << endl;

  //cout << "Threshold " << threshold << endl;

  double weight_b = 
    LinearInterpolate(threshold, threshMin_, threshMax_, blockWeightVect_[llhIndex]);
 
  double weight_new = weight*weight_b;
  
  weight_new /= ( (weight*weight_b) + ((1. - weight)*(1. - weight_b)) );
  //cout << nSrcGlobal << " " << weight_new << " " << weight_b << endl;
      // kludge for testing, only works for two datasets.
      //  Still not sure how to do it for many datasets.
  
  vector<double> parLocal(4); 
  
  // nSrc , must be scaled to relative weight of this source
  parLocal[0] = nSrcGlobal * weight_new; //localweight_[llhIndex];
  parLocal[1] = gamma;
  parLocal[2] = parGlobal[2];
  parLocal[3] = parGlobal[3];
  return parLocal;
}
  

void NewLlhBlockTime_ParTranslator::
SetTranslator(const vector<AnalysisSet*>& aSetVect) {
  int nSets = aSetVect.size();

  cout << "in par translator" << endl;

  srcWeightVect_.clear();
  blockWeightVect_.clear();
  dataTMinVect_.clear();
  dataTMaxVect_.clear();
  vector<double> tempVect(nStopsGamma_,0);
  
  cout << "Starting with N sets: " << nSets << endl;
  
  for (int i=0; i<nSets; ++i) {
    srcWeightVect_.push_back( tempVect );
    blockWeightVect_.push_back( tempVect );

    // Is this horrible? I have a feeling this might be horrible.
    // Also, this isn't completely correct, but should be good to 
    // within a second or so if constructed properly.
    dataTMinVect_.push_back( 
        (dynamic_cast<I3Analysis*>(aSetVect[i]))->GetEventTimeModulePtr()->GetTimeMin().GetTime() );


    dataTMaxVect_.push_back( 
        (dynamic_cast<I3Analysis*>(aSetVect[i]))->GetEventTimeModulePtr()->GetTimeMax().GetTime() );

    cout << "dataTMaxVect " << dataTMaxVect_[i] << " dataTMinVect " << dataTMinVect_[i] << endl;
        
  }

  threshMax_ = ( GetHighestBlock(BlocksTimeFile_) + GetSecondHighestBlock(BlocksTimeFile_) ) /2.;
  //threshMax_ = 1.79746e-06;

  cout << "highest second block " << GetSecondHighestBlock(BlocksTimeFile_) << endl;
  cout << "threshMax is: " << threshMax_ << endl;
  
  threshMin_ = 0.;
  double thSum = 0.;
  
  for (int j=0; j < nStopsGamma_; j++) {
    double th = threshMin_ + (threshMax_-threshMin_) * (double(j))/(nStopsGamma_ - 1.);
    thSum = 0.;
    for (int i=0; i<nSets; ++i) {
      double thNorm = BlockNormForThresh(BlocksTimeFile_,th,dataTMinVect_[i],dataTMaxVect_[i]);
      blockWeightVect_[i][j] = thNorm;
      thSum += thNorm;
    } 
    cout << "Weights for Threshold=" << th << ": " << flush;
    for (int i=0; i<nSets; ++i) {
      blockWeightVect_[i][j] /= thSum;
      cout << blockWeightVect_[i][j] << " " << flush;
    }
    cout << endl;
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
    cout << "Weights for Gamma=" << gamma << ": " << flush;
    for (int i=0; i<nSets; ++i) {
      srcWeightVect_[i][g] /= sum;
      cout << srcWeightVect_[i][g] << " " << flush;
    }
    cout << endl;
  }    
}

