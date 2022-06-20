#include "llh/public/llh_discovery_potential.h"

#include<vector>

#include "rootExt/public/CountMonitor.h"
#include "rootExt/public/generalfunctions.h"  // for poisson_prob
#include "rootExt/public/log_report.h"
#include "rootExt/public/MakeNewName.h"
#include "rootExt/public/randomfunctions.h"

#include "llh/public/classes.h"


DiscoveryPotential::~DiscoveryPotential () {
  if (trials_)        { delete trials_; }
  if (detections_)    { delete detections_; }
  if (fraction_)      { delete fraction_; }
  if (fractionError_) { delete fractionError_; }
  
  if (h_nSrc_logP)            { delete h_nSrc_logP; }
  if (hDiscoveryFraction)     { delete hDiscoveryFraction; }
  if (detectionRate_vs_Mean)  { delete detectionRate_vs_Mean; }
}


void DiscoveryPotential::InitializePlots() {
  int nsMaxStart = 40;
  
  if (h_nSrc_logP) { delete h_nSrc_logP; }
  h_nSrc_logP = new TH2D(MakeNewName("h_nSrc_logP"),
			 "-Log_{10} Prob  vs. nSrcEvents;-log_{10} Prob",
			 50, 0., 5., nsMaxStart, 0, nsMaxStart);
  h_nSrc_logP->SetCanExtend(TH1::kAllAxes);

  if (hDiscoveryFraction) { delete hDiscoveryFraction; }
  hDiscoveryFraction = new TH1D(MakeNewName("hDiscoveryFraction"),
				"Discovery Fraction",
				nsMaxStart, -0.5, nsMaxStart-0.5);
  hDiscoveryFraction->SetCanExtend(TH1::kAllAxes);
}


void DiscoveryPotential::SetForUpperLimit() { 
  log_fatal("DiscoveryPotential::SetForUpperLimit() has been replaced by\n"
	    "DiscoveryPotential::SetForMedianUpperLimit().\n"
	    "If you really want to fix the significance as in the old\n"
	    "upper limit calculation to 0.5, you can SetDetectionSignificance"
	    "and SetDetectionPower manually.\n");
}



bool DiscoveryPotential::DoTrial(int ns) {

  if (ns>=arraySize_) {
    log_error("Tried to set nSrcEvents above arraySize_=%d",arraySize_);
    log_error("(perhaps detectionSignificance=%lg is set too small?)\n",
	      detectionSignificance_);
    log_fatal("Fatal Error.");
    return false;
  }
   
  aSet_->GenerateDataSet_with_nSrcEvents(ns);
  llh_->MaximizeLlh();
  ++trials_[ns];

  if (ns>nsMaxSampled_) { nsMaxSampled_ = ns;}

  double prob = llh_->GetEstProb();
  h_nSrc_logP->Fill( -log10(prob), ns);

//   bool detected = false;
//   if ( prob < detectionSignificance_) {
//     ++detections_[ns];
//     detected = true;
//   }
  
  
  bool detected = detectionByTS_ ? (llh_-> GetTestStatistic() > detectionTS_) : (prob < detectionSignificance_); //ternary conditional: evaluates to b if value of a is true , otherwise to c
  if (detected){
      ++detections_[ns];
  }

  fraction_[ns] = double(detections_[ns]) / trials_[ns];

  return detected;
}



double DiscoveryPotential::MedianProbForBkg(int bTrials) {

  if ( !llh_ || !aSet_ ) { 
    log_error("MedianProbForBkg: AnalysisFn or AnalysisSet not set!!!\n");
    return -1.;
  }

  vector<double> bkgProbVect;
  CountMonitor cm(10., bTrials);

  for (int i=0; i<bTrials; ++i) {
    if (monitor_) { cm.UpdateCount(); }

    aSet_->GenerateDataSet_with_nSrcEvents(0);
    llh_->MaximizeLlh();
    bkgProbVect.push_back(llh_->GetEstProb());
  }

  sort(bkgProbVect.begin(), bkgProbVect.end());
  int medianElement = int( (bkgProbVect.size()-1) / 2);
  // e.g. if 5 elements, take [2] (third element)
  //      if 6 elements, take [2] (third element)
  double medBkgProb = bkgProbVect[medianElement];

  return medBkgProb;
}




void DiscoveryPotential::AnalyzeDiscoveryPotential() {

  // IF MEDIAN UPPER LIMIT IS TO BE USED:
  if (optMedianUpperLimit_) {
    if (monitor_) {
      cout<<" Calculating Median Bkg. p-value ("<<nBkgTrials_<<" trials):\n ";
    }

    detectionSignificance_ = MedianProbForBkg(nBkgTrials_);

    if (monitor_) {
      cout << " Set Detection Significance to median bkg p-value: ";
      cout << detectionSignificance_ << endl;
      cout << " Detection Power: " << detectionPower_ << endl;
    }
  } 
  else {
    if (monitor_) {
      cout << " Detection Significance: " << GetDetectionSignificance();
      cout << "  Detection Power: " << GetDetectionPower() << endl;
    }
  }


  nsMaxSampled_ = 0;
  MeanSrcEv_ForDetection_ = 0.;
  powerUncertainty_ = 1.;

  if ( !llh_ || !aSet_) { 
    log_error("AnalyzeDiscoveryPotential: AnalysisFn or AnalysisSet not set!!!\n");
    return;
  }

  InitializePlots();

  for (int i=0; i<arraySize_; ++i) {
    trials_[i] = 0;
    detections_[i] = 0;
    fraction_[i] = 0.;
    fractionError_[i] = 1.;  // i.e. completely unknown
  }



  // METHOD 1

  if (method_ == 1) {

    // FIRST PASS: SAMPLE THE DISTRIBUTION, FIND 'ACTIVE' REGION

    int nDetections = 0;
    int loggerLastOut = 0;

    int ns = 0;

    while (nDetections < nDetectionsMax_) {

      bool detected = DoTrial(ns);
      nDetections += detected;

      // Decide whether to try higher or lower ns, depending on this result:
      
      if (detected && ns>0) { --ns; }
      if (!detected)        { ++ns; }


      // Monitoring
      if ( monitor_ && (nDetections*10)/nDetectionsMax_ > loggerLastOut) {
	cout << nDetections << " " << flush;
	++ loggerLastOut;
      }

    }

    // Finish monitoring...
    if (monitor_) {
      cout << endl << "Now Filling in so each nSrcEvents has at least ";
      cout << nTrialsMin_ << " trials...\nnSrcEvents = ";
    }

    // SECOND PASS: FILL IN FOR NSRCEVENTS VALUES WHICH WERE UNDER-SAMPLED

    // note: we will break out of this loop much earlier if statistics are good
    for (ns = 0; ns < arraySize_; ++ns) {

      if (monitor_) { cout << ns << " " << flush;  }

      // increase statistics for this ns bin
      while (trials_[ns] < nTrialsMin_) {
	DoTrial(ns);
      }

      // Decide if we have enough statistics
      if ( fraction_[ns]>MinFractionToFinish_ && ns == nsMaxSampled_ ) {
	break;
      }
    }

    if (monitor_) {  cout << endl; }

    // Determine Poisson Mean nSrcEvents

    {
      bool foundMean = false;
      for (double mean = incMean_; mean <= nsMaxSampled_; mean += incMean_) {

	double power, uncertainty;
	EstimatePowerForMean(mean, &power, &uncertainty);
      
	if (power>=detectionPower_ && !foundMean) {
	  foundMean = true;
	  MeanSrcEv_ForDetection_ = mean;
	  powerUncertainty_ = uncertainty;
	  break;
	}
      }
    
      if (!foundMean) 
	log_warn("Failed to find mean for detectionPower_=%lg (%lg%% detection rate)!!!\n", detectionPower_, detectionPower_*100.);
    }

  } // END METHOD 1



  // METHOD 2


  if (method_ == 2) {

    // FIRST PASS: Do trials until first detection

    int ns = 0;

    // DoTrial returns true or false for whether trial results in detection
    while ( !DoTrial(ns) ) { ++ns; }

//    cout << "after do trial " << endl;
//    cout << " ns " << ns << endl;



    // SECOND PASS: 

    double meanInitial = ns;  // a good starting point
    double meanInc = 1.;
    double meanBest, powerBest, uncertainty;
    FindThresholdMean(meanInitial, meanInc, meanBest, powerBest, uncertainty);

    for (int loop = 0; loop < loops_; ++loop) {

      if (monitor_) {
	cout << " Mean: " << meanBest << "   Power: ";
	cout << powerBest*100 << "% +/-" << uncertainty*100 << "%\n";
      }

      int totalTrials = 100;

      double remainderErrorSq = 1.;

      // Fill in the sampled range with better statistics where 
      // uncertainties are currently contributing the most

      for (ns=0; ns <= nsMaxSampled_; ++ns) {
	double errorContribution = 
	  weightedUncertaintySq_[ns] / (uncertainty*uncertainty);

//    cout << "errorContribution " << errorContribution << endl;

	remainderErrorSq -= errorContribution;

	double approxTrials = totalTrials * errorContribution;
	int nTrials = random_poisson(approxTrials);
	// Reason for poisson sampling here:
	//   If ns is very large, and total trials is small,
	// approxTrials will always be less than one...
	// and taken as an integer would always yield zero new trials.
	// So using random_poisson guarantees that we still generate trials

//    cout << " nTrials " << nTrials << endl;

	for (int i=0; i<nTrials; ++i) {
	  DoTrial(ns);
	}
      }

//    cout << "meanInitial " << meanBest << endl;

      meanInitial = meanBest;
      FindThresholdMean(meanInitial, meanInc, 
			meanBest, powerBest, uncertainty);

      if (fabs(meanBest-meanInitial)<2.*meanInc) { meanInc /= 5.;}
      // start compressing the stepsize as we stabilize and get closer
    }

//    cout << "end " << endl;
//    cout << "meanBest " << meanBest << endl;

    MeanSrcEv_ForDetection_ = meanBest;
    powerUncertainty_ = uncertainty;

  }



  // For nice plots

  for (int ns=0; ns<=nsMaxSampled_; ++ns) {
    fractionError_[ns] = uncertaintyEstimator(detections_[ns], trials_[ns]);

    int iBin = hDiscoveryFraction->FindBin( double(ns) );
    hDiscoveryFraction->SetBinContent( iBin, fraction_[ns] );
    hDiscoveryFraction->SetBinError(   iBin, fractionError_[ns] );
  }

  vector<double> x;
  vector<double> y;
  for (double mean = incMean_; mean <= nsMaxSampled_; mean += incMean_) {
    double power, uncertainty;
    EstimatePowerForMean(mean, &power, &uncertainty);
    x.push_back(mean);
    y.push_back(power);
  }
  if (detectionRate_vs_Mean) { delete detectionRate_vs_Mean; }
  detectionRate_vs_Mean = new TGraph(x.size(), &x[0], &y[0]);

}


// TO DO:: BETTER HANDLING OF REQUESTS WITH MEAN = 0 (POISSON NOT DEFINED)
// is there a solution that returns reasonable values rather than exits?

// Monotonic Assumption:
// We will always sample until last ns has fraction[ns]=1.
// Then we assume that all higher trials have same fraction =1,
// and the total uncertainty for this and all higher trials is
// given by the error for the last ns combined with the weight for the 
// remaining poisson distribution (inclusive of last ns)

void DiscoveryPotential::EstimatePowerForMean(double mean,
				       double *power, double *uncertainty) {
  if (mean<=0.) {
    log_fatal("Can't call for mean=%lf <= 0.\n",mean);
  }

  // First: extend the sampling if last ns does not have fraction = 1;
  // (also, make sure sampling occurs past requested mean)

  int ns = nsMaxSampled_;

  //cout << "nsMaxSampled_ " << ns << endl;

  while (fraction_[ns] < 1. || ns<mean) {
    // output is mainly for testing...
//    cout<<"ns=" << ns << " mean=" << mean << " fraction=" << fraction_[ns];
//    cout << "... raising sampling\n";
    ++ns;
    DoTrial(ns);
  }

  // Now, increment through each value of ns, calculate detection fraction
  // and uncertainty, and add this to the total power, 
  // with poisson weight according to the mean given as input

  *power = 0.;
  *uncertainty = 0.;
  double sumSqError = 0.;

  for (ns = 0; ns<=nsMaxSampled_; ++ns) {

    double weight = poisson_prob(mean, "eq", ns);

    // except for last ns, which carrries weight of all higher terms
    if (ns == nsMaxSampled_) {
      weight = poisson_prob(mean, "ge", ns);
    }

    *power += weight * fraction_[ns];


    double error = uncertaintyEstimator(detections_[ns], trials_[ns]);
    // equals 1 if trials==0 (i.e. skipped, so maximal uncertainty)

    weightedUncertaintySq_[ns] = weight*weight * error*error;
    
    sumSqError += weightedUncertaintySq_[ns];
  }

  *uncertainty = sqrt(sumSqError);

}




void DiscoveryPotential::FindThresholdMean(double meanInitial,
					   double meanInc,
					   double &meanBest,
					   double &powerBest,
					   double &uncertainty)
{
  // starting point is first "best" point, (but check not zero)
  meanBest = meanInitial;
  if (meanBest<=0) { meanBest = meanInc; }

  EstimatePowerForMean(meanBest, &powerBest, &uncertainty);

  // Monotonic Assumption: if we are below threshold, we only
  // search higher means for a higher threshold,
  // and if we are above threshold, we only search lower
  double dir;
  if (powerBest > detectionPower_ && meanBest > meanInc) { dir = -1.;}
  else { dir = 1.;}

  double mean = meanBest + dir*meanInc;
  //cout << "mean in find threshold mean " << mean << endl;
  double power;

  EstimatePowerForMean(mean, &power, &uncertainty);

  while ( fabs(power-detectionPower_) < fabs(powerBest-detectionPower_) ) {
    // we moved closer to threshold, so try yet again
    meanBest = mean;
    powerBest = power;
    mean = meanBest + dir*meanInc;
    if (mean<=0.) 
      {break;}  // we're done, can't go lower
    EstimatePowerForMean(mean, &power, &uncertainty);
  }
  // stop when we stop moving closer

  // do one more time with best settings, so member arrays get filled
  EstimatePowerForMean(meanBest, &powerBest, &uncertainty);  
}
  




// pEst is your estimate for the true binomial p-value
// Note that if pEst=0.5, then typical uncertainty is 1/sqrt(nTrials+1)
// while if pEst=0 or 1 (no or all successes), unc is 1/(nTrials+1)
// (for large nTrials, very low or high success rate is like Poisson stats.)

double uncertaintyEstimator(int nSuccess, int nTrials) {
  if (nTrials<=0) { return 1.; }
  double pEst = double(nSuccess)/nTrials;
  return sqrt(4*nTrials*pEst*(1.-pEst) + 1.)/(nTrials+1.);
}
