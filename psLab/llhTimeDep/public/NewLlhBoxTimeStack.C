#include "llhTimeDep/public/NewLlhBoxTimeStack.h"
#include "llhTimeDep/public/TimePdfCollection.h"

#include "TGraph.h"

#include "rootExt/public/generalfunctions.h"
#include "rootExt/public/log_report.h"
#include "rootExt/public/ModDistanceFn.h"

#include "fluxus/public/FluxFunction.h"

#include "llh/public/BkgSpaceProb.h"
#include "llh/public/EnergyProb.h"
#include "llh/public/I3Event.h"




void NewLlhBoxStackFnc(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
  double nSrc = par[0];
  double gamma = par[1];
  f = 0;

  if (tdBoxStack->optStoreRatios_) { tdBoxStack->eventRatioVect_.clear(); }

  double srcFrac = nSrc/tdBoxStack->nEventsTot_;   // relative weight of src term
  if (srcFrac > tdBoxStack->srcFracMax_) {
    log_warn("Trying to evaluate srcFrac=%f which is > srcFracMax_=%f",
             srcFrac,tdBoxStack->srcFracMax_);
    log_warn("Check your parameter settings and make them consistent.");
    log_warn("(If running optimized, logLambda tolerance could be exceeded.)");
  }

  ///
  // JD: I use the convention that flux is proportional to E^gamma
  //  SimpleEnergyProb uses E^-gamma... (sorry for confusion)
  ///

  // round a double to nearest 'long'  This avoids round returning a double of
  //  0.99... which (int) might cast down to zero.
  int gammaIndex=(int)lround(10*(-1.*gamma+4)-0.5); // for look-up in Stacking table
  // For now, table size and binning is fixed in macro_StackedWeights.C
  if (gammaIndex<0) gammaIndex=0;
  if (gammaIndex>29) gammaIndex=29;

  double stackGammaMin = -4;
  double stackGammaMax = -1;
  int nGammaBins = 30; // find a global way of setting this...
  // This gamma index corresponds to a bin center of:
  double gammaCenterThisBin = (gammaIndex+0.5)*(stackGammaMax-stackGammaMin)/nGammaBins + stackGammaMin; 

  double logLambda = 0.; 
  double ratioSum=0.;
  double srcWeightSum=0.;
  //int srcIndex=0;

  //BoxTimePdf *timePdf;
  for (int i=0; i<tdBoxStack->selectedList_.GetSize(); ++i) {
    double eRatio = 1.;
    const Event* event = tdBoxStack->selectedList_.GetEvent(i);
    if (tdBoxStack->optUseEnergy_) {
      eRatio = tdBoxStack->eProbVect_[i]->GetEnergyProbGamma(*event, gamma) /
        tdBoxStack->eProbVect_[i]->GetEnergyProbBkg(*event);
    }
    // CAN USE THIS STOPWATCH TO ISOLATE ONE PIECE OF THE LLH CALC.
    //  stopwatch_minuitFCN_.Start(false);  //  false = don't reset

    // Iterate over each source, summing up signal terms and src weights
    ratioSum=0.;
    srcWeightSum=0.;
    //srcIndex=0;
    //double evtTime = event->GetTime().GetMJD();
    for (unsigned srcIndex=0; srcIndex<tdBoxStack->srcCoords_.size(); srcIndex++) {
      double srcWeightInterpValue = 0;
      bool interpDone = 0;      // Nominal value of srcWeightTable before interpolation
      double valueAtThisBin = tdBoxStack->srcWeightsTable_[srcIndex][gammaIndex];
      if ( gammaIndex==0 && -1.*gamma<gammaCenterThisBin ) {
        double lowEdge = -4.;
        srcWeightInterpValue = valueAtThisBin * (-1.*gamma-lowEdge)/(gammaCenterThisBin-lowEdge);
        interpDone = 1;
      }
      if ( gammaIndex == 29 && -1.*gamma>gammaCenterThisBin ){
        double hiEdge = -1;
        srcWeightInterpValue = valueAtThisBin * (-1.*gamma-hiEdge)/(hiEdge-gammaCenterThisBin);
        interpDone = 1;
      }

      if (!interpDone) {
        int dir = 1;
        if (-1*gamma<gammaCenterThisBin) { dir = -1; }
        double valueAtNextBin = tdBoxStack->srcWeightsTable_[srcIndex][gammaIndex+dir];
        double gammaCenterNextBin = (gammaIndex + dir + 0.5)*(stackGammaMax-stackGammaMin)/nGammaBins + stackGammaMin;

        double fractionOffset = fabs((-1*gamma-gammaCenterThisBin)/(gammaCenterNextBin-gammaCenterThisBin));
        srcWeightInterpValue = valueAtThisBin + (valueAtNextBin-valueAtThisBin)*fractionOffset;
      }
      // Interpolation turned out to be important for MINUIT (smooth llh landscape)
//cerr << "i: " << i << " 1: " << srcWeightInterpValue << " 2: " << tdBoxStack->spaceRatioVects_[srcIndex][i] << " 3: " << eRatio << endl;
      assert(tdBoxStack->spaceRatioVects_[i][srcIndex]==tdBoxStack->spaceRatioVects_[i][srcIndex]);
      ratioSum += srcWeightInterpValue *
                tdBoxStack->spaceRatioVects_[i][srcIndex]*eRatio*tdBoxStack->timeRatioVects_[i][srcIndex];
      srcWeightSum += srcWeightInterpValue; // for normalization
      
      //const I3Event* i3ev = (dynamic_cast<const I3Event*>(event));
      //Printf("event = %d, weight = %f, spaceRatio = %e, eRatio = %e, timeRatio = %e, time = %f, energy = %f, eventID = %d, runID = %d", i, srcWeightInterpValue, tdBoxStack->spaceRatioVects_[i][srcIndex], eRatio, tdBoxStack->timeRatioVects_[i][srcIndex], i3ev->GetTime().GetMJD(), i3ev->GetParams().energyValue, i3ev->GetParams().eventID, i3ev->GetParams().runID);

      //delete timePdf;
    }
    
//cerr << "1: " << srcFrac << " 2: " << ratioSum << " 3: " << srcWeightSum << endl;
    //logLambda += log( srcFrac * ( tdBoxStack->spaceRatioVect_[i] * eRatio - 1) + 1);
    logLambda += log( srcFrac * ( ratioSum/srcWeightSum - 1) + 1);
    // DON'T FORGET TO STOP AFTERWARDS!
    //    stopwatch_minuitFCN_.Stop();
    if (tdBoxStack->optStoreRatios_) {
      tdBoxStack->eventRatioVect_.push_back(ratioSum/srcWeightSum);
    }
  }

  if (tdBoxStack->optimizeTolerance_ > 0.) {
    // correction for events skipped in optimization
    log_fatal("This is skipped unless you implemented optimization!\n (in that case, double-check that this code still works)\n");
    logLambda += (tdBoxStack->nEventsTot_ - tdBoxStack->selectedList_.GetSize())*log(1.-srcFrac);
  }
  logLambda += (tdBoxStack->nEventsTot_ - tdBoxStack->selectedList_.GetSize())*log(1.-srcFrac);

  if (tdBoxStack->monitorLevel_>1) {
    printf("LogLambda=%12.6lg  :  nSrc=%9.4lg  :  gamma=%5.4lg\n",
           logLambda, nSrc, gamma);
  }

  f = - logLambda; // What Minuit minimizes: -log L 
  return;
}

NewLlhBoxStack::NewLlhBoxStack()
{
  optUseEnergy_ = true;
  SetEMaxRatioWarnOnlyOnce(-1);  // default -1 means warn every time
  SetOptStoreRatios(false);
  SetStatWarnLevel(0);
  SetnSrcMin(0.);
  srcFracMax_ = 0.5; // default is half of total nEvents
  gammaMin_ = 0.;
  gammaMax_ = 0.;
  logLambdaBest_ = 0.;
  nSrcBest_ = 0.;
  gammaBest_ = 0.;
  chiSq_ = 0.;
  chiSqProb_ = 0.;
  nEventsTot_ = 0;
  monitorLevel_ = 0;
  optimizeTolerance_ = 0.;
  optimizeAngleDeg_ = 0.;

  tdBoxStack = this; 

  optParAuto_[0] = true;
  optParAuto_[1] = true;
  
  // These start when created, so stop immediately
  stopwatch_MaximizeLlh_.Stop();
  stopwatch_optimize_.Stop();
  stopwatch_minuitMigrad_.Stop();
  stopwatch_minuitFCN_.Stop();
}

void NewLlhBoxStack::SetSourceCoords(vector<EquatorialDeg>& sourceCoords) {
  for ( vector<EquatorialDeg>::iterator iter = sourceCoords.begin();
    iter != sourceCoords.end(); iter++ )
  {
    srcCoords_.push_back(*iter);
  }
  return;
}

void NewLlhBoxStack::SetSourceSigmas(vector<double>& sourceSigmas) {
  // In the past, this std::copy has caused seg faults...
  // Trying agin with resize first:
  srcSigmas_.resize(sourceSigmas.size());
  std::copy(sourceSigmas.begin(), sourceSigmas.end(), srcSigmas_.begin());
  return;
}

void NewLlhBoxStack::SetStackedWeightTable(vector<vector<double> > srcWeightsArray) {
  srcWeightsTable_.clear();

 // This compiled but segfaulted:
 //copy(srcWeightsArray.begin(), srcWeightsArray.end(), srcWeightsTable_.begin());

  for ( vector<vector<double> >::iterator iter = srcWeightsArray.begin();
    iter != srcWeightsArray.end(); iter++ )
  {
    srcWeightsTable_.push_back(*iter);
  }

/*
  // Logging, if desired...
  for (vector< vector<double> >::size_type u = 0; u < srcWeightsTable_.size(); u++)
  {
      for (vector<double>::size_type v = 0; v < srcWeightsTable_[u].size(); v++) {
         cout << srcWeightsTable_[u][v] << " ";
      }
      cout << endl;
  }
*/
}

void NewLlhBoxStack::SetSourceTminTmax(vector<double> minTime, vector<double> maxTime) {
  assert( (minTime.size()==maxTime.size()) );
  for(unsigned int i=0; i<minTime.size(); ++i){
    tminSrc_.push_back(minTime[i]);
    tmaxSrc_.push_back(maxTime[i]);
  }
}

void NewLlhBoxStack::SetPeriodTminTmax(double tmin, double tmax) {
    tminPer_ = tmin;
    tmaxPer_ = tmax;
}

void NewLlhBoxStack::OptimizeEventSelection() {
  log_fatal("Optimiziation not implemented for stacking yet!\n");
  return;
}


void NewLlhBoxStack::PrepareAnalysis() {
  if (!aSet_) { log_fatal("PrepareAnalysis: AnalysisSet was not set.\n"); }
  if (srcCoords_.size() == 0) { log_fatal("PrepareAnalysis: srcCoords were not set.\n"); }
  if (srcSigmas_.size() == 0) { log_fatal("PrepareAnalysis: srcSigmas were not set.\n"); }
  if (srcWeightsTable_.size() == 0) { log_fatal("PrepareAnalysis: srcWeightsTable not set.\n"); }

  const EventPtrList* evList = aSet_->GetEventPtrList();
  if (!evList) { log_fatal("PrepareAnalysis: EventPtrList was not set.\n"); }
  nEventsTot_ = evList->GetSize();
  if (!nEventsTot_) { log_fatal("PrepareAnalysis: EventPtrList was empty.\n"); }
  // or, is there a reason to allow zero events?

  selectedList_.Clear();
  spaceRatioVect_.clear();
  eProbVect_.clear();
  spaceRatioVects_.clear();
  timeRatioVects_.clear();
  //spaceRatioVects_.resize( srcCoords_.size() );

  eventRatioVect_.clear();

  if (optimizeTolerance_ > 0.) {
    stopwatch_optimize_.Start(false);  //  false = don't reset
    log_error("Optimization not implemented yet! Setting Tolerance to 0.\n");
    optimizeTolerance_ = 0;
    stopwatch_optimize_.Stop();
  }
  if (optimizeTolerance_ > 0.) {
    OptimizeEventSelection();
  }
  else
  { // no optimization
    vector<double> spaceRatioSrc;
    //spaceRatioSrc.resize(srcCoords_.size());
    vector<double> timeRatioSrc;
    //timeRatioSrc.resize(srcCoords_.size());

    BoxTimePdf *timePdf;   

    double spaceRatioSum;
 
    for (int i=0; i<nEventsTot_; ++i) {
      const I3Event* event = (dynamic_cast<const I3Event*> (evList->GetEvent(i)));
      assert(event);

      double evtTime = event->GetTime().GetMJD();
      bool inSrcTime=false;
      timeRatioSrc.clear();
      timeRatioSrc.resize(srcCoords_.size());
      for(unsigned k=0; k<tminSrc_.size(); ++k){
        if(evtTime >= tdBoxStack->tminSrc_[k] && evtTime <= tdBoxStack->tmaxSrc_[k]){
          timePdf = new BoxTimePdf(tminPer_, tmaxPer_, tminSrc_[k], tmaxSrc_[k], 1., 0.);
          timeRatioSrc[k] = timePdf->GetPdfValue( evtTime  ) * timePdf->GetLivetime();;
          inSrcTime = true;
          delete timePdf;
        }
        else timeRatioSrc[k]=0;
      }
      if(inSrcTime){ //if event is in the time range of at least one source, add it
        int srcIndex=0;
        spaceRatioSrc.clear();
        spaceRatioSrc.resize(srcCoords_.size());
        spaceRatioSum=0;
        for ( vector<EquatorialDeg>::iterator iter = srcCoords_.begin(); iter != srcCoords_.end(); iter++ ) {
          double r = event->GetCoord().DistanceTo(*iter);
          double es = event->GetParams().parafitSigmaDeg; // event sigma
          double ss = srcSigmas_[srcIndex]; // source sigma
          double sigma = sqrt( (es*es + ss*ss) ); // Convolved error regions
          double sigSpaceProb = CircularGaussUnc(r, sigma);
          //double sigSpaceProb = exp(-r*r/(sigma*sigma*2)) / (2.*TMath::Pi()*sigma*sigma );
          double bkgSpaceProb = event->GetBkgSpaceProbFn()->GetBkgProbDensity(*event);
          //    if event is in region with b=0 (what to do?) 
          if (bkgSpaceProb <=0) { log_fatal("bkgSpaceProb <= 0\n"); }
          spaceRatioSrc[srcIndex] = sigSpaceProb / bkgSpaceProb;
          spaceRatioSum += sigSpaceProb / bkgSpaceProb;
          // energy prob
          //const EnergyProb* eProb(NULL);
          //if (optUseEnergy_) { eProb = event->GetEnergyProbFn(); }
          //eProbVect_.push_back(eProb);

          srcIndex++;
        }
        if(spaceRatioSum>1e-5){
            spaceRatioVects_.push_back(spaceRatioSrc);
            timeRatioVects_.push_back(timeRatioSrc);
            // energy prob
            const EnergyProb* eProb(NULL);
            if (optUseEnergy_) { eProb = event->GetEnergyProbFn(); }
            eProbVect_.push_back(eProb);
            selectedList_.AddEvent(event);
        }
      }
    }
  } // End of no optimization option
}

double NewLlhBoxStack::EvalFCN(const vector<double>& parVect) const {
  assert(nEventsTot_);
  double f;
  const int npar = parVect.size();
  double gin[npar];
  double par[npar];
  int iflag = 0;
  for(int i=0; i<npar; i++) par[i] = parVect[i];
  int npars=npar;
  NewLlhBoxStackFnc(npars, gin, f, par, iflag);

  return f;
}

double NewLlhBoxStack::EvaluateLlh(double nSrc, double gamma) {
  vector<double> parVect;
  parVect.push_back(nSrc);
  parVect.push_back(gamma);
  double minusLlh = EvalFCN(parVect);
  return -minusLlh;   // that is, max llh = - (minimizer result)
}
