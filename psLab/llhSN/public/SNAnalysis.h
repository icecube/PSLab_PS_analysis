#ifndef LLH_SNANALYSIS_H_
#define LLH_SNANALYSIS_H_

#include <vector>

#include "llh/public/classes.h"
#include "llh/public/EventTimeModule.h"
#include "llhSN/public/SNEvent.h"
#include "llhSN/public/SimpleProb.h"


class SNAnalysis : public AnalysisSet {
 protected: 

  bool randomizeBase_;
//  bool randomizeSrc_;

  vector<SNEvent> baseEvents_;  // Does Not Change 
  vector<SNEvent> modifiableBaseEvents_; // Scrambles; What eList_ points to
  
  double tmin_;
  double tmax_;

  EventPtrList eList_;
  
  SimpleProb *sProb_; 
  EventTimeModule* evTimeModulePtr_;

  // Protected Functions
  void BindBaseEvents();

 public:

  SNAnalysis();  
  ~SNAnalysis() { };
  
  double BkgNumberDensity(const Coord&) const {
//    if ( Coord ) { }
    return baseEvents_.size()/(tmax_ - tmin_);
  }
  
  const SourceModule* GetSource() const { return srcModule_; } 

  //  *  CONFIGURATION SETTINGS  *  //

  void SetRandomizeBase(bool randomize) { randomizeBase_ = randomize; }
  void SetBaseEvents(vector<SNEvent> &inputEvents);
  
  vector<SNEvent> GetEvents() { return modifiableBaseEvents_; }

  void SetSimpleProb(SimpleProb &inputProb);

  // The Module will keep track of event times; note that it is not
  // part of the SNAnalysis object itself, so it should *NOT* be deleted
  void SetEventTimeModulePtr( EventTimeModule* evTimeModule);

  //  *  OPERATIONS  *  //

  SimpleProb* GetSimpleProbFn() { return sProb_; }
  EventTimeModule* GetEventTimeModulePtr() { return evTimeModulePtr_; }

  void GenerateDataSet_with_nSrcEvents(int ns);
  void GenerateDataSet();

  void UseRealData();

  const EventPtrList* GetEventPtrList() const { return &eList_; }

};


#endif // LLH_SNANALYSIS_H_

