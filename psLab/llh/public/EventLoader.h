#ifndef LLH_EVENTLOADER_H_
#define LLH_EVENTLOADER_H_

#include <vector>

#include "TCut.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"

#include "llh/public/CoordEquatorialDeg.h"
#include "llh/public/EventTimeModule.h"
#include "llh/public/I3Event.h"


class EventLoader {
 public:
  vector<I3Event> rawBkgEventVect_;
  vector<double> bkgWeightVect_;

  vector<I3Event> rawSrcEventVect_;
  vector<double> srcMCZenRadVect_;
  vector<double> srcMCAziRadVect_;
  vector<double> srcMCOneWeightVect_;
  
  

  bool monitor_;
  bool applyAngFloor_;
  bool newDataStructure_; //to adapt event loader to both old and new data files
  bool randomizeOnlyTime_;

  // All types
  vector<TString> cutStringVect_;
  TString name_timeMJD_;
  TString name_recoZenith_rad_;
  TString name_recoAzimuth_rad_;
  TString name_sigmaDeg_;
  TString name_energyValue_;
  TString name_runID_;
  TString name_eventID_;

  // Bkg types
  TTree *bkgTree_;
  TString bkgWeight_;
  double minZenithRad_;
  double maxZenithRad_;
  double spreadZenithRad_;

  // See explanation below along with SetLoadMethod functions
  enum LoadMethod {
    EXACT,
    SELECTIVE_SAMPLE,
    POISSON_SAMPLE_PLUS,
    FIXED_UPSAMPLE
  };
  LoadMethod bkgLoadMethod_;

  enum TimeMethod {
    ACTUAL,
    SCRAMBLE
  };

  TimeMethod bkgTimeMethod_;
  EventTimeModule* evTimeModulePtr_;

  // Source types
  TTree *sourceTree_;
  double zenWidthDeg_;
  vector< vector<double> > zenWidthDegRangeVect_;

  TString name_mcRa_rad_;
  TString name_mcDec_rad_;
  TString name_mcZenith_rad_;
  TString name_mcAzimuth_rad_; 
  
  TString name_mcEnergy_GeV_;
  TString name_mcOneWeight_;
  double totalGeneratedEvents_;

  EventLoader() { 
    monitor_ = true;
    applyAngFloor_ = true;
    newDataStructure_= true;
    randomizeOnlyTime_=false;

    name_timeMJD_ = "UNDEFINED";
    name_recoZenith_rad_ = "UNDEFINED";
    name_recoAzimuth_rad_ = "UNDEFINED";
    name_sigmaDeg_ = "UNDEFINED";
    name_energyValue_ = "UNDEFINED";
    name_runID_ = "UNDEFINED";
    name_eventID_ = "UNDEFINED";

    bkgTree_ = NULL;
    bkgWeight_ = "UNDEFINED";
    minZenithRad_ = 0.;
    maxZenithRad_ = 0.;
    spreadZenithRad_ = 0.;
    bkgLoadMethod_ = EXACT;
    bkgTimeMethod_ = SCRAMBLE;
    evTimeModulePtr_ = NULL;

    sourceTree_ = NULL;
    zenWidthDeg_ = 0.;

    name_mcEnergy_GeV_ = "UNDEFINED";
    name_mcOneWeight_ = "UNDEFINED";
    totalGeneratedEvents_ = 0.;
  }

  void SetMonitor(bool monitor) { monitor_ = monitor; }
  void SetApplyAngFloor(bool flag) { applyAngFloor_ = flag; }
  bool GetApplyAngFloor() { return applyAngFloor_; }
  
  void SetNewDataStructure(bool flag) { 
      newDataStructure_ = flag; 
      
      if (newDataStructure_) {
        name_mcRa_rad_ = "UNDEFINED";
        name_mcDec_rad_ = "UNDEFINED";
      }
      else {
        name_mcZenith_rad_ = "UNDEFINED";
        name_mcAzimuth_rad_ = "UNDEFINED";
      }
  }
  
  bool GetDataStructure() { return newDataStructure_; }
  
  void AddCut(const char* cut) { cutStringVect_.push_back(cut); }

  TCut GetCuts() {
    TCut cut;
    for (unsigned int i=0; i<cutStringVect_.size(); i++) {
      cut = cut * cutStringVect_[i];
    }
    return cut;
  }
  
  void ResetCuts() { cutStringVect_.clear(); }

  void SetZenithRangeDeg(double minZenithDeg, double maxZenithDeg,
			  double spreadZenithDeg) {
    minZenithRad_    = minZenithDeg    * TMath::DegToRad();
    maxZenithRad_    = maxZenithDeg    * TMath::DegToRad();
    spreadZenithRad_ = spreadZenithDeg * TMath::DegToRad();
  }

  double GetZenithRangeDegMin() const { 
    return minZenithRad_ * TMath::RadToDeg(); }
  double GetZenithRangeDegMax() const { 
    return maxZenithRad_ * TMath::RadToDeg(); }
  double GetZenithRangeDegSpread() const { 
    return spreadZenithRad_ * TMath::RadToDeg(); }


  void SetName_recoZenith_rad(const char* name) {name_recoZenith_rad_ = name;}
  void SetName_recoAzimuth_rad(const char* name){name_recoAzimuth_rad_ = name;}
  void SetName_sigmaDeg(const char* name) {name_sigmaDeg_ = name;}
  void SetName_energyValue(const char* name) {name_energyValue_ = name;}
  void SetName_runID(const char* name) {name_runID_ = name;}
  void SetName_eventID(const char* name) {name_eventID_ = name;}
  void SetRndOnlyTime(bool flag) {randomizeOnlyTime_ = flag;}

  void SetBkgTree(TTree *tree) {
    bkgTree_ = tree; 
  }

  TTree* GetBkgTree() const { return bkgTree_; }

  // LOAD METHODS:
  // First, always check if event passes cut.  
  // If so, handle as follows:

  // EXACT: No weighting; load each event. (e.g. real data)

  // SELECTIVE_SAMPLE: Load event if a rand.number between (0,1) beats weight
  // (i.e. each event is either taken or not; if w>=1, event is always taken.)

  // POISSON_SAMPLE_PLUS: Generate random number N from Poisson mean = weight
  // if N=0, skip event
  // if N=1, load event
  // if N>1, load event, and make N-1 psuedo-events and load them as well
  // ('PLUS' is to remind you that it's sampling PLUS creating psuedo-events)

  // Comments:
  // IF all weights are much less than 1, then POISSON_SAMPLE_PLUS 
  // will be fairly similar to SELECTIVE_SAMPLE.  
  // (i.e., most events get skipped, a few get loaded once, and virtually none
  // get sampled twice.)
  // This means you have high statistics for the study at hand, and it's
  // probably best to use SELECTIVE_SAMPLE, so that no artifacts from 
  // the psuedo-event generation are introduced.

  // FIXED_UPSAMPLE: Each event is loaded int(weight) times.  After the first 
  // time an event is loaded, the next N-1 times are psuedo-events.
  // The typical use-cases are to upsample real data e.g. weight="20", or
  // to upsample simulation, e.g. weight="20.*BartolFluxWeight".  Once loaded,
  // the generator can make scrambled data sets at different thinning settings
  // (e.g. 1/20., or 1/10.) to simulate 1*livetime, 2*livetime, etc. datasets.

  void SetBkgLoadMethod_Exact() { 
    bkgLoadMethod_ = EXACT; 
    bkgWeight_ = "1.";
  }

  void SetBkgLoadMethod_SelectiveSample(TString weight) {
    bkgLoadMethod_ = SELECTIVE_SAMPLE;
    bkgWeight_ = weight;
  }

  void SetBkgLoadMethod_PoissonSamplePlus(TString weight) {
    bkgLoadMethod_ = POISSON_SAMPLE_PLUS;
    bkgWeight_ = weight;
  }

  void SetBkgLoadMethod_FixedUpsample(TString weight) {
    bkgLoadMethod_ = FIXED_UPSAMPLE;
    bkgWeight_ = weight;
  }

  void SetTimeMethod_Actual(TString name_timeMJD) {
    bkgTimeMethod_ = ACTUAL;
    name_timeMJD_ = name_timeMJD;
  }

  // The Module will keep track of event times; note that it is not
  // part of the I3Analysis object itself, so it should *NOT* be deleted
  void SetTimeMethod_Scramble(EventTimeModule* evTimeModule) {
    bkgTimeMethod_ = SCRAMBLE;
    if (newDataStructure_){
        name_timeMJD_ = "timeMJD";
    }
    else {
        name_timeMJD_ = "0";  // won't actually be used
    }
    evTimeModulePtr_ = evTimeModule;
  }

  // For Now, this function provides some backward compatibility
  void SetTimeMethod_Scramble() {
    bkgTimeMethod_ = SCRAMBLE;
    if (newDataStructure_){
        name_timeMJD_ = "timeMJD";
    }
    else {
        name_timeMJD_ = "0";  // won't actually be used
    }
    evTimeModulePtr_ = NULL;  // will get set when loading starts
  }

  void InstallEvents(vector<I3Event>& eventVect, TTree* tree, bool installSrc);
  void UnInstallBkgEvents() {
      rawBkgEventVect_.clear();
  }

  void LoadBkgEvents(vector<I3Event>& bkgEvents);


  void SetSourceTree(TTree *tree);

  TTree* GetSourceTree() const { return sourceTree_; }

  // select events within +/- zenWidth of source zenith (constant at pole)
  void SetSourceZenWidthDeg(double zenWidthDeg) {zenWidthDeg_ = zenWidthDeg;}

  // new method to add a zenWidth of source in particular zenith range
  void SetSourceZenWidthDegRange(double zenWidthDeg, double zenMin, double zenMax) {
    vector<double> zenWidthRange;

    zenWidthRange.push_back(zenWidthDeg);
    zenWidthRange.push_back(zenMin);
    zenWidthRange.push_back(zenMax);
    
    zenWidthDegRangeVect_.push_back(zenWidthRange);
    
  }

  void LoadSourceEvents(vector<I3Event>& srcCandidateEvents,
			EquatorialDeg srcEqDeg, bool rotate=true);

  double GetZenWidthDeg() { return zenWidthDeg_; }
};


#endif // LLH_EVENTLOADER_H_
