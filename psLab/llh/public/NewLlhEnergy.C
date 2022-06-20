#include "llh/public/NewLlhEnergy.h"

#include "TGraph.h"

#include "rootExt/public/generalfunctions.h"
#include "rootExt/public/log_report.h"
#include "rootExt/public/ModDistanceFn.h"

#include "fluxus/public/FluxFunction.h"

#include "llh/public/BkgSpaceProb.h"
#include "llh/public/EnergyProb.h"
#include "llh/public/I3Event.h"

void llhEnergyFnc(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
  double nSrc = par[0];
  double gamma = par[1];

  vector<double>* eventRatioVect = llh->GetEventRatios();
  EventPtrList selectedList = llh->GetSelectedList(); 
  vector<const EnergyProb*> eProbVect = llh->GetEProbVect();
  vector<double> spaceRatioVect = llh->GetSpaceRatioVect();
  if (llh->GetOptStoreRatios()) { eventRatioVect->clear(); }

  double srcFrac = nSrc/llh->Get_nEvents();   // relative weight of src term
  if (srcFrac > llh->GetOptimizeSrcFracMax()) {
    log_warn("Trying to evaluate srcFrac=%f which is > srcFracMax_=%f",
	     srcFrac,llh->GetOptimizeSrcFracMax());
    log_warn("Check your parameter settings and make them consistent.");
    log_warn("(If running optimized, logLambda tolerance could be exceeded.)");
  }

  if (llh->GetUseEnergy() && (gamma >=eProbVect[0]->GetGammaMax() || gamma < eProbVect[0]->GetGammaMin() || TMath::IsNaN(gamma) ) ) {
    f = 1e50;
    return;
  }

  double logLambda=0.;
  for (int i=0; i<selectedList.GetSize(); ++i) {

    double eRatio = 1.;
    if (llh->GetUseEnergy()) {
      const Event* event = selectedList.GetEvent(i);
      eRatio = eProbVect[i]->GetEnergyProbGamma(*event, gamma) /
	eProbVect[i]->GetEnergyProbBkg(*event);
    }
    // CAN USE THIS STOPWATCH TO ISOLATE ONE PIECE OF THE LLH CALC.
    //  stopwatch_minuitFCN_.Start(false);  //  false = don't reset
  
    logLambda += log( srcFrac * ( spaceRatioVect[i] * eRatio - 1) + 1);

    // DON'T FORGET TO STOP AFTERWARDS!
    //    stopwatch_minuitFCN_.Stop();

    if (llh->GetOptStoreRatios()) {
      eventRatioVect->push_back(spaceRatioVect[i]*eRatio);
    }      
  }

  if (llh->GetOptimizeTolerance() > 0.) { 
    // correction for events skipped in optimization
    logLambda += (llh->Get_nEvents() - selectedList.GetSize())*log(1.-srcFrac);
  }

  if (llh->GetMonitorLevel()>1) {
    printf("LogLambda=%12.6lg  :  nSrc=%9.4lg  :  gamma=%5.4lg\n",
	   logLambda, nSrc, gamma);
  }


  f = -logLambda;
  return;   // What Minuit minimizes: -log L
}




NewLlhEnergy::NewLlhEnergy() :
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
  optimizeAngleDeg_(0.)
{
  llh = this; //llh defined at the end of NewLlhEnergy.h

  minuit_ = new TMinuit(2);
  minuit_->SetFCN(llhEnergyFnc);

  optParAuto_[0] = true;
  optParAuto_[1] = true;

  // These start when created, so stop immediately
  stopwatch_MaximizeLlh_.Stop();
  stopwatch_optimize_.Stop();
  stopwatch_minuitMigrad_.Stop();
  stopwatch_minuitFCN_.Stop();
}

NewLlhEnergy::
~NewLlhEnergy()
{
    if (minuit_) {
      delete minuit_;
      minuit_ = NULL;
    }
}


void NewLlhEnergy::OptimizeEventSelection() {
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
    const I3Event* event = 
      (dynamic_cast<const I3Event*> (evList->GetEvent(i)));
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
      eProbVect_.push_back(eProb);
      spaceRatioVect_.push_back(spaceRatio);
    }
  }

  if (monitorLevel_ > 0) {
    printf("Optimizing: %d events selected with maxRatio >%lg"
	   "  out of %d Ntotal.\n",
	   selectedList_.GetSize(), threshold, nEventsTot_);
  }
}



void NewLlhEnergy::PrepareAnalysis() {
  if (!aSet_) { log_fatal("PrepareAnalysis: AnalysisSet was not set.\n"); }
  if (!srcCoord_) { log_fatal("PrepareAnalysis: srcCoord was not set.\n"); }

  const EventPtrList* evList = aSet_->GetEventPtrList();
  if (!evList) { log_fatal("PrepareAnalysis: EventPtrList was not set.\n"); }
  nEventsTot_ = evList->GetSize();
  if (!nEventsTot_) { log_fatal("PrepareAnalysis: EventPtrList was empty.\n"); }
  // or, is there a reason to allow zero events?

  selectedList_.Clear();
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
      selectedList_.AddEvent(event);
      eProbVect_.push_back(eProb);
      spaceRatioVect_.push_back(spaceRatio);
    }
  }
}


void NewLlhEnergy::MaximizeLlh()
{
  stopwatch_MaximizeLlh_.Start(false);  //  false = don't reset

  PrepareAnalysis();
  // TO DO: BETTER HANDLING IF NO EVENTS SELECTED:
  assert(selectedList_.GetSize());

  // Wipes out parameters (minuit crashes easily when changing existing pars)
  Double_t arglist[1];
  arglist[0] = 0.5; //0.5 for likelihood fit
  Int_t ierflg=0;
 
  minuit_->mnexcm("SET ERR",arglist,1,ierflg);
  minuit_->SetPrintLevel(-1);

  // DEFINE PARAMETERS

  if (optParAuto_[0]) {
    double nSrcMax = srcFracMax_*nEventsTot_;
    nSrcMin_ = 0.;
    minuit_->mnparm(0, "nSrc", nSrcMax/2., 0.1, nSrcMin_, nSrcMax, ierflg);
  }

  if (optParAuto_[1] ) {
    if (optUseEnergy_) {
      // Use eProb of *first* event... assume gamma range same for all????
      gammaMin_ = eProbVect_[0]->GetGammaMin();
      gammaMax_ = eProbVect_[0]->GetGammaMax();
      double gammaInit = (gammaMin_+gammaMax_)/2.;

      minuit_->mnparm(1, "gamma", gammaInit, 0.1, gammaMin_, gammaMax_, ierflg);
    } else {
      // If not using energy, then FCN will actually replace the
      // energy prob term internally... so these values don't matter...
      // Here we simply fix it (implied by the range) to zero 
      //so Minuit will not try to minimize
      minuit_->mnparm(1, "gamma", 0., 0., 0., 0., ierflg);
    }
  }

  // MINIMIZE

  stopwatch_minuitMigrad_.Start(false);  //  false = don't reset
  double arglist2[2];
  arglist2[0] = 500;
  arglist2[1] = 0.1;
  minuit_->mnexcm("MIGRAD", arglist2 ,2,minuitOut_);
  stopwatch_minuitMigrad_.Stop();


 //  if (icstat <= icstatWarnLevel_) {
//     log_warn("icstat=%d : ",icstat);
//     if (icstat == 0) { log_warn("covar. matrix not calculated at all.\n"); }
//     if (icstat == 1) { log_warn("covar. matrix is approx. only\n"); }
//     if (icstat == 2) { log_warn("full covar.matrix, but forced pos-def.\n");}
//     if (icstat == 3) { log_warn("full, accurate covariance matrix.\n");}
//   }


  // GET RESULTS


  StoreLogLambdaBest();   // Set logLambdaBest_
  nSrcBest_ = GetPar(0);
  gammaBest_ = GetPar(1);

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
  if (optUseEnergy_) {
    chisq_prob(chiSq_, 2., &p_temp, &chiSqProb_);
  } else {
    chisq_prob(chiSq_, 1., &p_temp, &chiSqProb_);
  }

  estProb_ = chiSqProb_ / 2.;  // one-sided chi-sq prob

  stopwatch_MaximizeLlh_.Stop();
}


// This is clumsy but seems like the only way to get this info
void NewLlhEnergy::StoreLogLambdaBest() 
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


double NewLlhEnergy::EvalFCN(const vector<double>& parVect) const {
  assert(nEventsTot_);
  double f;
  const int npar = parVect.size();
  double gin[npar];
  double par[npar];
  int iflag = 0;
  for(int i=0; i<npar; i++) par[i] = parVect[i];
  int npars=npar;
  llhEnergyFnc(npars, gin, f, par, iflag);

  return f;
}



double NewLlhEnergy::EvaluateLlh(double nSrc, double gamma) {
  vector<double> parVect;
  parVect.push_back(nSrc);
  parVect.push_back(gamma);
  double minusLlh = EvalFCN(parVect);
  return -minusLlh;   // that is, max llh = - (minimizer result)
}




// TGraph* NewLlhEnergy::GetContour(double sigma, int npoints, int pa1, int pa2) {
//   minuit_->SetFCN(MinuitWrapperFCN);
//   SetMinuitWrapperPtr(this);  // point to *our* EvalMinuitFCN

//   minuit_->SetErrorDef(sigma*sigma);
//   TGraph* g = dynamic_cast<TGraph*> (minuit_->Contour(npoints, pa1, pa2));

//   SetMinuitWrapperPtr(NULL);  // reset pointer
//   return g;
// }






// For MultiAnalysis: Translate Global Par values to Individual Par Values:

const vector<double> NewLlhEnergy_ParTranslator::
Translate(int llhIndex, const vector<double>& parGlobal) const {
  double nSrcGlobal = parGlobal[0];
  double gamma = parGlobal[1];

  vector<double> parLocal(2);

  if(TMath::IsNaN(gamma) || TMath::IsNaN(nSrcGlobal))
  {       
      parLocal[0] = 0;
      parLocal[1] = 1;
      return parLocal;
  }



  double weight = 
    LinearInterpolate(gamma, gammaMin_, gammaMax_, srcWeightVect_[llhIndex]);
  // recall: srcWeightVect_[i] is itself a vector<double> 

  // nSrc , must be scaled to relative weight of this source
  parLocal[0] = nSrcGlobal * weight;
  parLocal[1] = gamma;
  return parLocal;
}
  

void NewLlhEnergy_ParTranslator::
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
    cout << "Weights for Gamma=" << gamma << ": " << flush;
    for (int i=0; i<nSets; ++i) {
      srcWeightVect_[i][g] /= sum;
      cout << srcWeightVect_[i][g] << " " << flush;
    }
    cout << " " << endl;
  }    
}

