#ifndef LLH_I3ANALYSIS_H_
#define LLH_I3ANALYSIS_H_

#include <vector>

#include "llh/public/BkgSpaceProb.h"
#include "llh/public/classes.h"
#include "llh/public/EventTimeModule.h"
#include "llh/public/I3Event.h"
#include "llh/public/I3SignalGenerator.h"
#include "llh/public/SimpleEnergyProb.h"


class I3Analysis : public AnalysisSet {
 protected: 

  bool randomizeBase_;
  bool randomizeSrc_;
  bool rndSrcEventTime_;
  bool addSignalToEProb_;
  double baseThinningProb_;
  bool isFirstInjection_; //to be used when you want to inject signal in multiple coordinates
  bool isLastInjection_;  //to be used when you want to inject signal in multiple coordinates
  bool isStacking_;

  vector<I3Event> baseEvents_;  // Does Not Change 
  vector<I3Event> modifiableBaseEvents_; // Scrambles; What eList_ points to
  vector<I3Event> redModifiableBaseEvents_; // reduced modifiable base events: all events but those corresponding to entries in selEvTreeName_, if randomizing only times
  vector<I3Event> sourceEvents_; // also what eList_ points to
  EventPtrList eList_;
  //adding list for only source events
  EventPtrList sList_;

  BkgSpaceProb *bkgSpaceProb_;
  SimpleEnergyProb *eProb_; 

  EventTimeModule* evTimeModulePtr_;

  //sel = selected = events in selEvTreeName_
  //src = source = signal events to be injected
  bool isRndOnlyTime_;   //if randomizing only times
  char* selEvFileName_;  //name of file containing tree of events to be injected, if randomizing only times
  char* selEvTreeName_;  //name of tree of events to be injected, if randomizing only times
  int selEvNo_;  //if randomizing only times, entries in selEvTreeName_
  vector<int> selEvIDVec_;  //vectors containg ID of events in source tree, needed is randomizing only times
  vector<I3Event> selEvts_;  //if randomizing only times: vector contains events of
                             //the injection tree that corresponds to evts of this dataset
  vector<int> modifiableSrcIdxVec_;  //if randomizing only times and injecting signal from selEvTreeName_,
                               //this vector stores indeces of events in the tree not yet injected
  vector<int> srcIdxVec_;  //if randomizing only times, keep memory of the indeces
                             //of events that will be injected as signal.
  vector<int> srcIDVec_;  //ID of evts to be injected as signal
  
  // Protected Functions

  void BindBaseEvents();

  // nSrcEvents<0 ==> generate random number based on source mean
  // nSrcEvents>=0  ==> specifically generate nSrcEvents
  void GenerateSrcEvents(vector<I3Event>& SrcEventVector, int nSrcEvents);

  vector<int> GetSelEvID();
  bool IsSelEvent(I3Event evt);
  bool GenerateSrcFromTree(vector<I3Event>& sourceEventsVector, int nSrcEvents);
  bool IsSrcEv(I3Event evt);
  bool IsTimeInGRL(double mjd);
  vector<I3Event> GenerateBkgSample();

  int injFlare_;  //Current flare considered for injection
  int nFlares_; //Total number of flares to be injected

 public:

  I3Analysis(); 
  ~I3Analysis();


  //  *  CONFIGURATION SETTINGS  *  //

  void SetRandomizeBase(bool randomize) { randomizeBase_ = randomize; }
  void SetRandomizeSrc(bool randomize) { randomizeSrc_ = randomize; }
  void SetRndSrcEventTime(bool randomize) { rndSrcEventTime_ = randomize; }
  void SetAddSignalToEProb(bool add) { addSignalToEProb_ = add; }
  void SetBaseThinningProb(double prob) { baseThinningProb_ = prob; }
  void SetIsFirstInjection(bool first) { isFirstInjection_ = first; }
  void SetIsLastInjection(bool last) { isLastInjection_ = last; }
  void SetIsStacking(bool isStack) {isStacking_ = isStack; }

  void SetBaseEvents(const vector<I3Event> &inputEvents);

  // creates a clone copy which must eventually be deleted
  void SetBkgSpaceProb(const BkgSpaceProb& bsp);

  // Note this will get modified each time the data + injected signal 
  // is passed to eProb_ for updating the bkg tables!
  // (The problem with making a local copy is derived classes..., we 
  // might make a Clone function for EnergyProb's)
  void SetEnergyProb(SimpleEnergyProb &inputEProb);


  // creates a clone copy which must eventually be deleted
  void SetSource(const SourceModule& src);

  // The Module will keep track of event times; note that it is not
  // part of the I3Analysis object itself, so it should *NOT* be deleted
  void SetEventTimeModulePtr( EventTimeModule* evTimeModule);

  // -- This is now handled inside of EventTimeModule
  //  void SetTimeRange(double a, double b);
  
  //  *  OPERATIONS  *  //

  double GetBaseThinningProb() const { return baseThinningProb_; }
  const BkgSpaceProb* GetBkgSpaceProbFn() const { return bkgSpaceProb_; }
  const EnergyProb* GetEnergyProbFn() const { return eProb_; }
  const SourceModule* GetSource() const { return srcModule_; }
  EventTimeModule* GetEventTimeModulePtr() { return evTimeModulePtr_; }

  // Returns current number density, with fake signal added to data
  double BkgNumberDensity(const Coord& coord) const;

  void Inject_nSrcEvents(int nSrcEvents);
  // nSrcEvents<0 ==> generate random number based on source mean
  // nSrcEvents>=0  ==> specifically generate nSrcEvents
  void GenerateDataSet_with_nSrcEvents(int nSrcEvents);

  void GenerateDataSet() { GenerateDataSet_with_nSrcEvents(-1); }

  void UseRealData();

  const EventPtrList* GetEventPtrList() const { return &eList_; }
  const EventPtrList* GetEventPtrListSrc() const { return &sList_; }

  bool GetIsRndOnlyTime() { return isRndOnlyTime_; }
  int GetSelEvNo() const { return selEvNo_; }
  int GetNflares() const { return nFlares_; }
  vector<int> GetModifiableSrcIdxVec() const { return modifiableSrcIdxVec_; }
  vector<int> GetSrcIdxVec() const { return srcIdxVec_; }

  void SetRndOnlyTime(bool flag) { isRndOnlyTime_ = flag; }
  void SetSelEvFileName(char* name) { selEvFileName_ = name; }
  void SetSelEvTreeName(char* name) { selEvTreeName_ = name; }
  void SetModifiableSrcIdx(vector<int> vect) {
    modifiableSrcIdxVec_.clear();
    modifiableSrcIdxVec_ = vect;
  }
  void SetSrcID();
  void SetSrcIdx(vector<int> vect){
    srcIdxVec_ = vect;
    SetSrcID();
  }
  void SetNFlares(int nFl){ nFlares_ = nFl; }
  void SetInjFlare(int iFl){ injFlare_ = iFl; }
};


#endif // LLH_I3ANALYSIS_H_

