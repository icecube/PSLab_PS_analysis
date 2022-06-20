#include "llhSN/public/SNAnalysis.h"

#include "TMath.h"

#include "rootExt/public/log_report.h"
#include "rootExt/public/randomfunctions.h"


// PROTECTED FUNCTIONS

void SNAnalysis::BindBaseEvents() {
  tmin_ = 1e6;
  tmax_ = 0.;

  // wait until all of these ingredients have been set
  if (baseEvents_.size()>0 && sProb_) {
    for (int i=0; i<int(baseEvents_.size()); ++i) {
      baseEvents_[i].SetAnalysisSet(this); 
      baseEvents_[i].SetSimpleProb(sProb_);
      if (baseEvents_[i].GetMJD()<tmin_) {tmin_ = baseEvents_[i].GetMJD();}
      if (baseEvents_[i].GetMJD()>tmax_) {tmax_ = baseEvents_[i].GetMJD();}
    }
    modifiableBaseEvents_ = baseEvents_;
  }
}


// PUBLIC FUNCTIONS


SNAnalysis::SNAnalysis() :
  sProb_(NULL),
  evTimeModulePtr_(NULL)
  //srcModule_(NULL)
{ }

void SNAnalysis::SetBaseEvents(vector<SNEvent> &inputEvents)
{
  baseEvents_ = inputEvents;
  BindBaseEvents();
}

void SNAnalysis::SetSimpleProb(SimpleProb &inputProb) { 
  sProb_ = &inputProb; 
  BindBaseEvents();
}


// The Module will keep track of event times; note that it is not
// part of the SNAnalysis object itself, so it should *NOT* be deleted
void SNAnalysis::SetEventTimeModulePtr( EventTimeModule* evTimeModule) {
    evTimeModulePtr_ = evTimeModule;
}


void SNAnalysis::GenerateDataSet() {
  eList_.Clear();
  
  if ( !evTimeModulePtr_ ) {
    evTimeModulePtr_ = new EventTimeModule();

    // (very minor memory leak here if never deleted, but it only happens
    //  once so not very worrisome)

    // This maintains backward compatibility for time-independent scripts...
  }

  // Clear the list (if any) of used times which the module tracks internally
  evTimeModulePtr_->ResetUsedTimes();


  // Get Modifiable version of Base (i.e. background data) events

  for (int i=0; i <int(modifiableBaseEvents_.size()); ++i) {

    // If desired, randomize time, according to method in evTimeModule
    if (randomizeBase_) {
      modifiableBaseEvents_[i].SetTime(evTimeModulePtr_->GetRandomTime());
    }

    eList_.AddEvent(&(modifiableBaseEvents_[i]));
  }
  
  eList_.SortByTime();
}

void SNAnalysis::GenerateDataSet_with_nSrcEvents(int ns) {
  if (ns) { cout << "Sorry, no signal hypothesis" << endl;}
  GenerateDataSet();
}

void SNAnalysis::UseRealData() {

  modifiableBaseEvents_ = baseEvents_;
  
  eList_.Clear();
  for (int i=0; i<int(modifiableBaseEvents_.size()); ++i) {
    eList_.AddEvent( &(modifiableBaseEvents_[i]) );
  }

  eList_.SortByTime();

  sProb_->SetTableBkg(modifiableBaseEvents_);
}
