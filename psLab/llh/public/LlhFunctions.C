#include "llh/public/LlhFunctions.h"

#include "TMath.h"
#include "TMinuit.h"

#include "rootExt/public/generalfunctions.h"

#include "llh/public/BkgSpaceProb.h"


BinnedAnalysis::BinnedAnalysis() : 
  AnalysisLlh(),
  binSize_(0),
  psfFractionInBin_(1.),
  profileLogProb_(NULL),
  counts_(0),
  bkgDensity_(0),
  nBkgMean_(0),
  nSrcMin_(0),
  nSrcBest_(0)
{ }


void BinnedAnalysis::PrepareAnalysis()
{
  if (!aSet_) { log_fatal("PrepareAnalysis: AnalysisSet was not set.\n"); }
  if (!srcCoord_) { log_fatal("PrepareAnalysis: srcCoord was not set.\n"); }

  distanceVector_.clear();
  const EventPtrList* evList = aSet_->GetEventPtrList();
  for (int i=0; i<evList->GetSize(); ++i) {
    const Event* event = evList->GetEvent(i);
    distanceVector_.push_back(event->GetCoord().DistanceTo(*srcCoord_));
  }

  bkgDensity_ = aSet_->BkgNumberDensity(*srcCoord_);
  nBkgMean_ = bkgDensity_ * binSize_ * binSize_ * TMath::Pi();
}


void BinnedAnalysis::DoAnalysis() 
{
  PrepareAnalysis();

  counts_ = 0;

  // FOR OVERALL RESULTS: Only used if user specified profileLogProb_ histogram
  TH1I *hBinCounts = NULL;
  int nBins_hBinCounts = 0;
  if (profileLogProb_) 
  {
    nBins_hBinCounts = profileLogProb_->GetNbinsX();
    hBinCounts = new TH1I("hBinCounts","hBinCounts",
			  nBins_hBinCounts,
			  profileLogProb_->GetXaxis()->GetXmin(),
			  profileLogProb_->GetXaxis()->GetXmax() );
  }


  // MAIN LOOP OVER EVENTS
  for (vector<double>::iterator itDist = distanceVector_.begin(); 
       itDist != distanceVector_.end(); itDist++) {

    if (*itDist <= binSize_) { ++counts_; }  // main result

    if (profileLogProb_) { hBinCounts->Fill(*itDist); }
  }


  // FOR OVERALL RESULTS
  if (profileLogProb_) {  // give user: p-values for each possible bin size
    int sum = 0;    
    // recall, bin zero is underflow bin, so bin one is first bin...
    for (int b=1; b<=nBins_hBinCounts; b++) {
      sum += int(hBinCounts->GetBinContent(b));
      double binCenter = hBinCounts->GetBinCenter(b);
      double binRadius = hBinCounts->GetBinLowEdge(b+1);
      double binArea = binRadius * binRadius * TMath::Pi();
      double nBkg = bkgDensity_ * binArea;
      profileLogProb_->Fill(binCenter, log10(poisson_prob(nBkg, "ge", sum)));
    }

    delete hBinCounts;

    // profileLogProb_ results will accumulate if user analyzes many
    // scrambled sets and keeps using profileLogProb_
  }

  nSrcBest_ = max(counts_-nBkgMean_ , nSrcMin_*psfFractionInBin_ ) 
    / psfFractionInBin_;

  // Best Hyp. / Null Hyp.
  logLambdaBest_ = 
    log(poisson_prob(nBkgMean_+(nSrcBest_*psfFractionInBin_),"eq",counts_)) - 
    log(poisson_prob(nBkgMean_,"eq",counts_));
    ;
}


double BinnedAnalysis::GetPoissonProb() const 
{
  return poisson_prob(nBkgMean_,"ge",counts_); 
}

double BinnedAnalysis::EvaluateLlh(double nSrc) {
  DoAnalysis();
  return 
    log(poisson_prob(nBkgMean_+(nSrc*psfFractionInBin_),"eq", counts_)) -
    log(poisson_prob(nBkgMean_,"eq",counts_));

}



/* 
   Recall likelihood formula:

              Product ( Weight*ProbSource[i] + WeightBkg*ProbBkg[i] ) 
    lambda = -------------------------------------------------------------
              Product ( ProbBkg[i] )

   which can be simplified as seen in the code below
*/



double llhRatio(vector<double>& ProbRatio, double WeightS)
{
  double WeightBkg = 1. - WeightS;
  double LogLambda=0.;

  for (vector<double>::const_iterator itRatio = ProbRatio.begin();
       itRatio != ProbRatio.end(); itRatio++)
  {
    LogLambda += log( WeightS * (*itRatio) + WeightBkg);
  }

  return LogLambda;
}



SimpleLlh::SimpleLlh() : 
  AnalysisLlh(1),
  method_(SIMPLE),
  nSrcInc_(0.1),
  nSrcMin_(0.),
  nSrcMax_(0.),
  opt_Default_SrcRange_(true),
  logLambdaBest_(0.),
  chiSq_(0.),
  chiSqProb_(0.),
  nSrcBest_(0.),
  nEvents_(0)
{
  parDefArray_[0].name = "nSrc";
}

void SimpleLlh::PrepareAnalysis() {
  if (!aSet_) { log_fatal("PrepareAnalysis: AnalysisSet was not set.\n"); }
  if (!srcCoord_) { log_fatal("PrepareAnalysis: srcCoord was not set.\n"); }

  ProbRatio_.clear();
  const EventPtrList* evList = aSet_->GetEventPtrList();
  for (int i=0; i<evList->GetSize(); ++i) {
    const Event* event = evList->GetEvent(i);
    double sigSpaceProb = event->ProbFrom(*srcCoord_);
    double bkgSpaceProb = 
      aSet_->BkgNumberDensity(event->GetCoord())/evList->GetSize();
    assert(bkgSpaceProb); // if event is in region with b=0 (what to do?)
    ProbRatio_.push_back(sigSpaceProb / bkgSpaceProb);
  }
}

void SimpleLlh::MaximizeLlh_Simple()
{
  PrepareAnalysis();

  nEvents_ = ProbRatio_.size();

  if (opt_Default_SrcRange_) {
    nSrcMax_ = nEvents_;
  }


  // Deciding whether to look up or down for best fit:
  //   For WeightS=0, LogLambda always equals 0.
  //   So, look at next highest WeightS, that is, WeightS = WeightSInc;
  //   If this LogLambda > 0, then WeightSBest will positive.

  nSrcBest_ = 0.;
  double nSrc = nSrcInc_;
  double logLambda = llhRatio(ProbRatio_, nSrc/nEvents_);

  if (logLambda>0) { 
    // okay, max LogLambda will be for positive WeightS
    nSrcBest_ = nSrc;
    logLambdaBest_ = logLambda;
    for (nSrc += nSrcInc_; nSrc<=nSrcMax_; nSrc+=nSrcInc_) {
      logLambda = llhRatio(ProbRatio_, nSrc/nEvents_);
      if (logLambda>=logLambdaBest_) {
	logLambdaBest_ = logLambda;
	nSrcBest_ = nSrc;
      } else {
	break;  // we already found the best, LogLambda is going down now
      }
    }

  } else { // max LogLambda must be for negative WeightS
    nSrcBest_ = 0.;
    logLambdaBest_ = 0.;
    for (nSrc = -nSrcInc_; nSrc>=nSrcMin_; nSrc -= nSrcInc_) {
      logLambda = llhRatio(ProbRatio_, nSrc/nEvents_);
      if (logLambda>=logLambdaBest_) {
	logLambdaBest_ = logLambda;
	nSrcBest_ = nSrc;
      } else {
	break;  // we already found the best, LogLambda is going down now
      }
    }
  }

  chiSq_ = 2.*logLambdaBest_;
  double p_temp;
  chisq_prob(chiSq_, 1., &p_temp, &chiSqProb_);

}


double SimpleLlh::EvaluateLlh(double nSrc) {
  PrepareAnalysis();

  nEvents_ = ProbRatio_.size();
  if (nEvents_ == 0) {
    cout << "Error. No events, cannot evaluate Llh.\n";
    return 0.;
  }

  if (method_ == SIMPLE) {
    return llhRatio(ProbRatio_, nSrc/nEvents_);
  } 
  else if (method_ == MINUIT_MIGRAD) {
    int npar = 1;
    double gin[1];
    double f;
    double par[1];
    int iflag = 0;
    par[0] = nSrc;
    EvalMinuitFCN(npar, gin, f, par, iflag);
    return (-f);  // We want Llh, not -Llh
  } 
  else {
    cout << "Error: SimpleLlh Method not recognized.\n";
    return 0.;
  }
}


void SimpleLlh::EvalMinuitFCN (int &npar, double *gin, double &f, 
			       double *par, int iflag)
{  
  double WeightS = par[0] / nEvents_;          // Minuit Input

  double WeightBkg = 1. - WeightS;
  double LogLambda=0.;

  for (vector<double>::const_iterator itRatio = ProbRatio_.begin();
       itRatio != ProbRatio_.end(); itRatio++)
  {
    LogLambda += log( WeightS * (*itRatio) + WeightBkg);
  }

  f = - LogLambda;    // Minuit Output

  if (npar || gin || iflag) { }; // touch these variables just to
  // eliminate warnings during compile
}




void SimpleLlh::MaximizeLlh_MinuitMigrad()
{
  nEvents_ = ProbRatio_.size();

  if (opt_Default_SrcRange_) {
    nSrcMin_ = -1.; // Really, default means min=0, but that is
    // a problem for the minimizer, so we will fix this at the end...
    nSrcMax_ = nEvents_;
  }


  // Define parameters
  double initValue = 0.;
  double initStepSize = .1;
  double lowLimit = nSrcMin_;
  double upLimit =  nSrcMax_;
  minuit_->DefineParameter(0, "Nsrc", initValue, initStepSize, lowLimit, upLimit);
  // minuit_->FixParameter(0);  // if you want to fix something


  // Set the Minuit FCN indirectly, by using a wrapper fcn which 
  // internally will point to 'this' object and its own EvalMinuitFCN

  minuit_->SetFCN(MinuitWrapperFCN);
  SetMinuitWrapperPtr(this);  // point to *our* EvalMinuitFCN
  minuit_->Migrad();          // do the minimization
  SetMinuitWrapperPtr(NULL);  // reset pointer


  double nSrcError;
  minuit_->GetParameter(0, nSrcBest_, nSrcError);

  double fmin, fedm, errdef;
  int nvpar, nparx, icstat;
  minuit_->mnstat(fmin, fedm, errdef, nvpar, nparx, icstat);
  logLambdaBest_ = -fmin;


  // if we really intended for nSrcMin = 0, then we have to fix now,
  // in case best fit was for negative nSrc
  if (opt_Default_SrcRange_) {
    if (nSrcBest_ < 0.) {
      nSrcBest_ = 0.;
      logLambdaBest_ = 0.;
    }
  }

  chiSq_ = 2.*logLambdaBest_;
  double p_temp;
  chisq_prob(chiSq_, 1., &p_temp, &chiSqProb_);
}
