#include "llh/public/I3Analysis.h"

#include "TMath.h"
#include "TFile.h"
#include "rootExt/public/log_report.h"
#include "rootExt/public/randomfunctions.h"

#include "llh/public/CoordEquatorialDeg.h"

#include <fstream>

// PROTECTED FUNCTIONS

void I3Analysis::BindBaseEvents() {
  // wait until all of these ingredients have been set
  if (baseEvents_.size()>0 && bkgSpaceProb_ && eProb_) {
    if(!isRndOnlyTime_){
      for (int i=0; i<int(baseEvents_.size()); ++i) {
        baseEvents_[i].SetAnalysisSet(this);
        baseEvents_[i].SetBkgSpaceProb(bkgSpaceProb_);
        baseEvents_[i].SetEnergyProb(eProb_);
      }
      modifiableBaseEvents_ = baseEvents_;
    }
    else{
      selEvIDVec_.clear();
      selEvIDVec_ = GetSelEvID();      
      vector<I3Event> baseEventsNoSignificant;
      selEvts_.clear();
      for (int i=0; i<int(baseEvents_.size()); ++i) {       
        if(!IsSelEvent(baseEvents_[i])){
          baseEvents_[i].SetAnalysisSet(this);
          baseEvents_[i].SetBkgSpaceProb(bkgSpaceProb_);
          baseEvents_[i].SetEnergyProb(eProb_);
          baseEventsNoSignificant.push_back(baseEvents_[i]);
        }    
        else selEvts_.push_back(baseEvents_[i]);        
      }
      redModifiableBaseEvents_ = baseEventsNoSignificant;
    }
  }
}

// nSrcEvents<0 ==> generate random number based on source mean
// nSrcEvents>=0  ==> specifically generate nSrcEvents
void I3Analysis::GenerateSrcEvents(vector<I3Event>& SrcEventVector, 
				   int nSrcEvents)
{
  if (nSrcEvents != 0 && !srcModule_) {
    log_fatal("Source was not set... cannot generate source events.\n");
  }

  if(injFlare_==0 && isFirstInjection_ ) SrcEventVector.clear();

  int nGeneratedSrcEvents;
  if (nSrcEvents<0) {
    nGeneratedSrcEvents = random_poisson( srcModule_->GetMeanSrcNev() );
  } else {
    nGeneratedSrcEvents = nSrcEvents;
  }

  if (nGeneratedSrcEvents){
    I3SignalGenerator* src = dynamic_cast<I3SignalGenerator*>(srcModule_);
    if(isStacking_) src = dynamic_cast<I3MultiSignalGenerator*>(srcModule_);

    if (!src) {
      log_fatal("ERROR: Expected SourceModule of type I3SignalGenerator.\n");
    }
    if(isStacking_){
      for (int i=0; i<nGeneratedSrcEvents; ++i) {
        SrcEventVector.push_back( src->GenerateEvent() );
        int sz = SrcEventVector.size();
        SrcEventVector[sz-1].SetAnalysisSet(this);
        SrcEventVector[sz-1].SetBkgSpaceProb(bkgSpaceProb_);
        SrcEventVector[sz-1].SetEnergyProb(eProb_);
      }
    }
    else{
      if(src->GetTimePdf()){
        for (int i=0; i<nGeneratedSrcEvents; ++i) {
          SrcEventVector.push_back( src->GenerateEvent() );
          int sz = SrcEventVector.size();
          SrcEventVector[sz-1].SetAnalysisSet(this);
          SrcEventVector[sz-1].SetBkgSpaceProb(bkgSpaceProb_);
          SrcEventVector[sz-1].SetEnergyProb(eProb_);
        }
      }
      else if(src->GetTimePdfVect().at(0)){
        for(int iPdf=0; iPdf<nFlares_; iPdf++){
          injFlare_=iPdf;
          src->SetTimePdf(src->GetTimePdfVect().at(iPdf));
          for (int i=0; i<nGeneratedSrcEvents; ++i) {
            SrcEventVector.push_back( src->GenerateEvent() );
            int sz = SrcEventVector.size();
            SrcEventVector[sz-1].SetAnalysisSet(this);
            SrcEventVector[sz-1].SetBkgSpaceProb(bkgSpaceProb_);
            SrcEventVector[sz-1].SetEnergyProb(eProb_);
          }
          src->ResetTimePdf();
        }
      }
    }
  }
}


// PUBLIC FUNCTIONS


I3Analysis::I3Analysis() :
  randomizeBase_(true), 
  randomizeSrc_(false),
  rndSrcEventTime_(false),
  addSignalToEProb_(true),
  baseThinningProb_(0.),
  isFirstInjection_(true),
  isLastInjection_(true),
  isStacking_(false),
  bkgSpaceProb_(NULL),
  eProb_(NULL),
  evTimeModulePtr_(NULL),
  isRndOnlyTime_(false),
  selEvFileName_((char *)""),
  selEvTreeName_((char *)""),
  selEvNo_(0),
  injFlare_(0),
  nFlares_(1)
{ }


I3Analysis::~I3Analysis()
{ 
  if (bkgSpaceProb_) { delete bkgSpaceProb_; }
  if (srcModule_) { delete srcModule_; }
}


void I3Analysis::SetSource(const SourceModule& src) {
  if (srcModule_) { delete srcModule_; }
  // this 'new' copy must be deleted in ~I3Analysis and whenever reset (here) 
  srcModule_ = src.Clone();
}


// Returns current number density, with fake signal added to data
double I3Analysis::BkgNumberDensity(const Coord& coord) const {
  return bkgSpaceProb_->GetBkgProbDensity(coord) * eList_.GetSize();
}


void I3Analysis::SetBaseEvents(const vector<I3Event> &inputEvents)
{
  baseEvents_ = inputEvents;
  BindBaseEvents();
}


void I3Analysis::SetBkgSpaceProb(const BkgSpaceProb& bsp) {
  if (bkgSpaceProb_) { delete bkgSpaceProb_; }
  // this 'new' copy must be deleted in ~I3Analysis and whenever reset (here) 
  bkgSpaceProb_ = bsp.Clone();
  BindBaseEvents();
}


// Note this will get modified each time the data + injected signal 
// is passed to eProb_ for updating the bkg tables!
// (The problem with making a local copy is derived classes..., we 
// might make a Clone function for EnergyProb's)
void I3Analysis::SetEnergyProb(SimpleEnergyProb &inputEProb) { 
  eProb_ = &inputEProb; 
  BindBaseEvents();
}


// The Module will keep track of event times; note that it is not
// part of the I3Analysis object itself, so it should *NOT* be deleted
void I3Analysis::SetEventTimeModulePtr( EventTimeModule* evTimeModule) {
    evTimeModulePtr_ = evTimeModule;
  }


void I3Analysis::Inject_nSrcEvents(int nSrcEvents) {
  if(isFirstInjection_) sourceEvents_.clear();
  EventPtrList sList;

  if(nSrcEvents == 0) Printf("Warning: trying to inject signla but 0 events defined for injection");
  else
    {
      if(!isRndOnlyTime_)  GenerateSrcEvents(sourceEvents_, nSrcEvents);
      else  GenerateSrcFromTree(sourceEvents_, nSrcEvents);

      if(injFlare_==nFlares_-1){
        for (unsigned int i=(sourceEvents_.size() - nSrcEvents); i < sourceEvents_.size(); ++i) {
                // If desired, scramble RA for this event
                if (randomizeSrc_) {
                  evTimeModulePtr_->RandomizeEvent( sourceEvents_[i] );
                }

                if(rndSrcEventTime_){
                  double mjd = (evTimeModulePtr_->GetRandomTime()).GetMJD();
                  sourceEvents_[i].SetMJD( mjd );
                }

//      double MJD = sourceEvents_[i].GetMJD();
//
//       //opening a file in writing mode
//     ofstream file("/data/user/sbron/psLab/macro_llh/IC86_II-VI_3C_279_ana_MESE/AnaSourcesCheckEventsInjection/3C_279/source_events.txt", ios::out | ios::app);
//
//     if(file)
//     {
//         file << MJD << endl;
//
//         file.close();
//     }
//     else
//             cerr << "Error at opening !" << endl;
// //

        //sList.AddEvent( &(sourceEvents_[i]));
        //sList_.AddEvent( &(sourceEvents_[i]));
        //eList_.AddEvent(&(sourceEvents_[i]));
      }
    }
  }

  if(sourceEvents_.size()!=0 && isLastInjection_) {
    for (unsigned int i=0; i < sourceEvents_.size(); ++i) {
      sList.AddEvent( &(sourceEvents_[i]));
      sList_.AddEvent( &(sourceEvents_[i]));
      eList_.AddEvent(&(sourceEvents_[i]));
    }
  }


  // are there cases where this is really not necessary, and could speed things
  // up if the sorting were skipped?  (and would it matter much?)
  if(isLastInjection_){
    eList_.SortByTime();

    if ( sList.GetSize() > 0 ) {
      bkgSpaceProb_->FixToBasePlusEvents(sList);
    }
    else {
      bkgSpaceProb_->FixToBase();
    }

    if (addSignalToEProb_) {
      eProb_->SetTableBkg(eList_);
    }
    injFlare_=0;
  }
}

void I3Analysis::GenerateDataSet_with_nSrcEvents(int nSrcEvents) {
  eList_.Clear();
  sList_.Clear();
  if(isRndOnlyTime_) modifiableBaseEvents_ = GenerateBkgSample();

  if ( !evTimeModulePtr_ ) {
    //cout << "!evTimeModulePtr" << endl;

    evTimeModulePtr_ = new EventTimeModule();

    // (very minor memory leak here if never deleted, but it only happens
    //  once so not very worrisome)

    // times from generic module with default constructor 
    // are only generated over one sidereal day... this is *supposed* to look
    // wierd if someone is expecting a realistic distribution of times
    // during the year...

    // This maintains backward compatibility for time-independent scripts...
    // eventually we may think to remove this and force _all_ scripts to
    // specify how time is handled.
  }
  // Clear the list (if any) of used times which the module tracks internally
  evTimeModulePtr_->ResetUsedTimes();


  // Get Modifiable version of Base (i.e. background data) events

  for (unsigned int i=0; i < modifiableBaseEvents_.size(); ++i) {

    // Use this to sample only a fraction of base 
    if (baseThinningProb_ > 0. && baseThinningProb_ < 1.) {
      if ( random_uniform(0.,1.) >= baseThinningProb_ ) {
         continue;  // skip if random number is ouside sampling range
      }
    }

    // If desired, randomize time and RA, according to method in evTimeModule
    if (randomizeBase_) {
      evTimeModulePtr_->RandomizeEvent( modifiableBaseEvents_[i], isRndOnlyTime_);
    }

    eList_.AddEvent(&(modifiableBaseEvents_[i]));
  }
  

  // Get Source events
  sourceEvents_.clear();
  EventPtrList sList;
  if(nSrcEvents != 0)
    {
      if(!isRndOnlyTime_) GenerateSrcEvents(sourceEvents_, nSrcEvents);
      else  GenerateSrcFromTree(sourceEvents_, nSrcEvents);
      for (unsigned int i=0; i < sourceEvents_.size(); ++i) {
          // If desired, scramble RA for this event
	    if (randomizeSrc_) {
	      evTimeModulePtr_->RandomizeEvent( sourceEvents_[i] );
	    }

//	double MJD = sourceEvents_[i].GetMJD();
//
//	 //opening a file in writing mode
//     ofstream file("/data/user/sbron/psLab/macro_llh/IC86_II-VI_3C_279_ana_MESE/AnaSourcesCheckEventsInjection/3C_279/source_events.txt", ios::out | ios::app);
//
//     if(file)
//     {
//         file << MJD << endl;
//
//         file.close();
//     }
//     else
//             cerr << "Error at opening !" << endl;
// //

        //sList.AddEvent( &(sourceEvents_[i]));
        //sList_.AddEvent( &(sourceEvents_[i]));
        //eList_.AddEvent(&(sourceEvents_[i])); 
    }
  }
  if(sourceEvents_.size()!=0) {
    for (unsigned int i=0; i < sourceEvents_.size(); ++i) {
      sList.AddEvent( &(sourceEvents_[i]));
      sList_.AddEvent( &(sourceEvents_[i]));
      eList_.AddEvent(&(sourceEvents_[i]));
    }
  }
  // are there cases where this is really not necessary, and could speed things
  // up if the sorting were skipped?  (and would it matter much?)
 
  eList_.SortByTime();
  if ( sList.GetSize() > 0 ) {
    bkgSpaceProb_->FixToBasePlusEvents(sList);
  }
  else {
    bkgSpaceProb_->FixToBase();
  }
  if (addSignalToEProb_) {
    eProb_->SetTableBkg(eList_);
  }  
}

bool I3Analysis::GenerateSrcFromTree(vector<I3Event>& sourceEventsVector, int nSrcEvents) {
  TFile f(selEvFileName_);
  if(f.IsZombie()) {
    cout << "Not existing file or wrong name\n";    
    return false;
  }
 
  TTree* tree = (TTree*)f.Get(selEvTreeName_);

  double time_mjd;
  double ra_deg;
  double dec_deg; 
  double RecoZenith_deg;
  double RecoAzimuth_deg;
  double parafitSigma_deg;
  double energy;
  int run_ID;
  int event_ID;
  double mc_energy;
  double PS_flatSpectrum_rate;
  double src_weight;
  
  tree->SetBranchAddress("time_mjd", &time_mjd);
  tree->SetBranchAddress("ra_deg", &ra_deg);
  tree->SetBranchAddress("dec_deg", &dec_deg);
  tree->SetBranchAddress("RecoZenith_deg", &RecoZenith_deg);
  tree->SetBranchAddress("RecoAzimuth_deg", &RecoAzimuth_deg);
  tree->SetBranchAddress("parafitSigma_deg", &parafitSigma_deg);
  tree->SetBranchAddress("energy", &energy);
  tree->SetBranchAddress("run_ID", &run_ID);
  tree->SetBranchAddress("event_ID", &event_ID);
  tree->SetBranchAddress("mc_energy", &mc_energy);
  tree->SetBranchAddress("PS_flatSpectrum_rate", &PS_flatSpectrum_rate);
  tree->SetBranchAddress("src_weight", &src_weight);

  I3Event evt;
  I3EventParameters evtParams;
  I3MCParameters evtMCParams;
  I3SignalGenerator* src = dynamic_cast<I3SignalGenerator*>(srcModule_);

  int treeIdx;
  for(Int_t i=0; i<nSrcEvents; ++i){
    treeIdx=random_uniform(0.,modifiableSrcIdxVec_.size());
    tree->GetEntry(modifiableSrcIdxVec_[treeIdx]);

    evtParams.recoZenithDeg = RecoZenith_deg;
    evtParams.recoAzimuthDeg = RecoAzimuth_deg;
    evtParams.parafitSigmaDeg = parafitSigma_deg;
    evtParams.energyValue = energy;
    evtParams.runID = run_ID;
    evtParams.eventID = event_ID;

    evtMCParams.mcEnergy = mc_energy;
    evtMCParams.PS_FlatSpectrumRate = PS_flatSpectrum_rate;
    evtMCParams.srcWeight = src_weight;

    evt.SetParams(evtParams);
    evt.SetMCParams(evtMCParams);
    evt.SetRaDeg(ra_deg);
    evt.SetDecDeg(dec_deg);

    if(src->GetTimePdf()){
      do{
        time_mjd = src->GetTimePdf()->GenerateEventTime().GetMJD();
      } while (!IsTimeInGRL(time_mjd));
    }
    else log_fatal("Source was not set... cannot generate source time.\n");

    evt.SetMJD(time_mjd);
    evt.SetAnalysisSet(this);
    evt.SetBkgSpaceProb(bkgSpaceProb_);
    evt.SetEnergyProb(eProb_);
    sourceEventsVector.push_back(evt);
    
    modifiableSrcIdxVec_.erase(modifiableSrcIdxVec_.begin()+treeIdx);
  }
  f.Close();
  return true;
}
 
bool I3Analysis::IsTimeInGRL(double mjd){
  I3SignalGenerator* src = dynamic_cast<I3SignalGenerator*>(srcModule_);
  vector<double> startMissRunVect = src->GetTimePdf()->GetStartMissRunVect();
  vector<double> stopMissRunVect = src->GetTimePdf()->GetStopMissRunVect();

  for(unsigned int i=0; i<startMissRunVect.size(); ++i){
    if(mjd > startMissRunVect[i] && mjd < stopMissRunVect[i]) return false;
  }
  return true;
}

vector<int> I3Analysis::GetSelEvID() {
  vector<int> vector_ID;

  TFile f(selEvFileName_);
  if(f.IsZombie()) log_fatal("File does not exist or has wrong name\n");
  TTree* tree = (TTree*)f.Get(selEvTreeName_);

  selEvNo_ = tree->GetEntries();

  int event_ID;
  tree->SetBranchAddress("event_ID", &event_ID); 
  
  for(int entry=0; entry<selEvNo_; ++entry){
    tree->GetEntry(entry);
    vector_ID.push_back(event_ID);
  }
  f.Close();
  return vector_ID;
}

bool I3Analysis::IsSelEvent(I3Event evt){
  int evtID = evt.GetParams().eventID;
  for(unsigned int i=0; i<selEvIDVec_.size(); ++i){
    if(selEvIDVec_[i] == evtID) return true;
  }
  return false;
}

void I3Analysis::SetSrcID(){
  TFile f(selEvFileName_);
  if(f.IsZombie()) log_fatal("File does not exist or has wrong name\n");
  TTree* tree = (TTree*)f.Get(selEvTreeName_);

  srcIDVec_.clear();
  int event_ID;
  tree->SetBranchAddress("event_ID", &event_ID);

  for(unsigned int i=0; i<srcIdxVec_.size(); ++i){
    tree->GetEntry(srcIdxVec_[i]);
    srcIDVec_.push_back(event_ID);
  }
  f.Close();
}

bool I3Analysis::IsSrcEv(I3Event evt){
  int evtID = evt.GetParams().eventID;
  for(unsigned int i=0; i<srcIDVec_.size(); ++i){
    if(evtID == srcIDVec_[i]) return true;
  }
  return false;
}

vector<I3Event> I3Analysis::GenerateBkgSample(){
  vector<I3Event> bkgEvts = redModifiableBaseEvents_; //contains all events but those
                                                      //corresponding events in source tree
  //Add to background the events in source tree that are not injected as signal
  for(unsigned int j=0; j<selEvts_.size(); ++j){
    if(!IsSrcEv(selEvts_[j])){
      selEvts_[j].SetAnalysisSet(this);
      selEvts_[j].SetBkgSpaceProb(bkgSpaceProb_);
      selEvts_[j].SetEnergyProb(eProb_);
      bkgEvts.push_back(selEvts_[j]);
    }
  }
  return bkgEvts;
}

void I3Analysis::UseRealData() {
  bkgSpaceProb_->FixToBase();
  modifiableBaseEvents_ = baseEvents_;
  
  eList_.Clear();
  for (unsigned int i=0; i<modifiableBaseEvents_.size(); ++i) {
    eList_.AddEvent( &(modifiableBaseEvents_[i]) );
  }

  // are there cases where this is really not necessary, and could speed things
  // up if the sorting were skipped?  (and would it matter much?)
  eList_.SortByTime();

  eProb_->SetTableBkg(eList_);
}
