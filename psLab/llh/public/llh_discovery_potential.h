#ifndef LLH_LLH_DISCOVERY_POTENTIAL_H_
#define LLH_LLH_DISCOVERY_POTENTIAL_H_

#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"


// Forward Declarations (when feasible, more efficient than including headers)
class AnalysisFn;
class AnalysisSet;


class DiscoveryPotential {

 private:
  AnalysisSet *aSet_;
  AnalysisFn *llh_;

  bool optMedianUpperLimit_;
  int nBkgTrials_;

  double detectionTS_;
  bool detectionByTS_;
  
  // rename to detectionProb_ ??
  double detectionSignificance_;
  double detectionPower_;

  bool DoTrial(int ns);

 public:

  int arraySize_;
  int *trials_;
  int *detections_;
  double *fraction_;
  double *fractionError_;
  double *weightedUncertaintySq_;

  int loops_;  // temporary for development

  int nsMaxSampled_;

  TH2D *h_nSrc_logP;
  TH1D *hDiscoveryFraction;
  TGraph *detectionRate_vs_Mean;
  double incMean_;

  double MeanSrcEv_ForDetection_;
  double powerUncertainty_;

  bool monitor_;

  int method_;

  // old method
  int nDetectionsMax_;
  int nTrialsMin_;
  double MinFractionToFinish_;



  DiscoveryPotential () {
    aSet_ = NULL;
    llh_ = NULL;

    // This is set in SetForDiscovery, below
    //    optMedianUpperLimit_ = false;
    
    detectionTS_= 0;
    detectionByTS_ = false;
    
    nBkgTrials_ = 1000;

    trials_ = NULL;
    detections_ = NULL;
    fraction_ = NULL;
    fractionError_ = NULL;
    weightedUncertaintySq_ = NULL;

    arraySize_ = 1000;  // can be changed by user later
    // but I think 1000 should be big enough for most applications
    SetArraySize(arraySize_);

    h_nSrc_logP = NULL;
    hDiscoveryFraction = NULL;
    detectionRate_vs_Mean = NULL;

    nsMaxSampled_ = 0;

    SetForDiscovery();

    monitor_ = false;

    incMean_ = 0.1;  // precision of mean determination

    MeanSrcEv_ForDetection_ = 0.;
    powerUncertainty_ = 1.;

    method_ = 1;
  }


  ~DiscoveryPotential ();


  // size is maximum number of nSrcEvents considered
  // (remember that poisson mean estimation may look beyond nSrcEvents
  // which have actually been simulated)
  void SetArraySize(unsigned int size) {
    if (trials_)        { delete trials_; }
    if (detections_)    { delete detections_; }
    if (fraction_)      { delete fraction_; }
    if (fractionError_) { delete fractionError_; }
    if (weightedUncertaintySq_) { delete weightedUncertaintySq_; }
    trials_        = new int [size];
    detections_    = new int [size];
    fraction_      = new double [size];
    fractionError_ = new double [size];
    weightedUncertaintySq_ = new double [size];
    arraySize_ = size;
  }


  void SetForDiscovery() {
    // Use the one-sided 5-sigma threshold for detection
    SetDetectionSignificance( 2.86651e-07 );
    // this function will also set 'optMedianUpperLimit_' to false

    SetDetectionPower( 0.5 );
    // we're interested in the flux required to
    // trigger a detection (defined above) in this fraction of trials

    // The setting of these variables is more important for upper limits...
    nTrialsMin_ = 20;
    MinFractionToFinish_ = 0.9;
  }

  void SetForUpperLimit(); // deprecated; fatal error will be thrown now.
    
  void SetForMedianUpperLimit(int nBkgTrials) {
    nBkgTrials_ = nBkgTrials;
    detectionSignificance_ = 0.;
    optMedianUpperLimit_ = true;
    // with this option set, the analysis will first begin with bkg-only
    // trials, in order to establish what the median p-value is, which
    // may differ from 50% exactly because the conversion from logLambda
    // to p-value is not always given perfectly by the chi-sq estimate.

    // detectionSignificance_ will be set to the median p-value found...
    // it should hopefully be close to 50% in most cases.

    // As long as optMedianUpperLimit_ is true, this bkg calculation will 
    // always over-ride whatever value was set previously for 
    // detectionSignificance_;

    detectionPower_ = 0.9; // we're interested in the flux required to
    // trigger a detection (i.e. be more significant than our typical
    // background case) in this fraction of trials

    nTrialsMin_ = 1000;  // need good sampling at high nSrcEvent values,
    // which aren't explored as well by first pass search.
    MinFractionToFinish_ = 0.99;  // want this higher than detection fraction,
    // so high nSrcEvents are well-sampled
  }

  bool GetOptMedianUpperLimit() const { return optMedianUpperLimit_; }

  void SetDetectionSignificance(double value) {
    detectionSignificance_ = value;
    optMedianUpperLimit_ = false;
    // explicitly setting detection significance has to override the
    // option of setting it equal to the median p-value of bkg-only trials
  }

  double GetDetectionSignificance() const {
    return detectionSignificance_;
  }

  void SetDetectionByTS (double TS) {
    detectionByTS_ = true;
    detectionTS_ = TS;
    optMedianUpperLimit_ = false;
  }
  
  bool GetDetectionByTS() const {
    return detectionByTS_;
  }
  
  void SetDetectionPower(double value) {
    detectionPower_ = value;
  }

  double GetDetectionPower() const {
    return detectionPower_;
  }


  void SetValues(int nDetectionsMax) {
    nDetectionsMax_ = nDetectionsMax;
  }


  void InitializePlots();

  void EstimatePowerForMean(double mean, double *power, double *uncertainty);

  void FindThresholdMean(double meanInitial, double meanInc,
			 double &meanBest, double &powerBest,
			 double &uncertainty);

  double MedianProbForBkg(int bTrials);

  void AnalyzeDiscoveryPotential();

  void SetAnalysisFn(AnalysisFn* llh) { llh_ = llh;}
  
  void SetAnalysisSet(AnalysisSet* aSet) { aSet_ = aSet;}

  const AnalysisSet* GetAnalysisSet() {return aSet_;}

};




double uncertaintyEstimator(int nSuccess, int nTrials);

#endif // LLH_LLH_DISCOVERY_POTENTIAL_H_
