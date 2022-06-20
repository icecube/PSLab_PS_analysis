#ifndef LLH_MULTIANALYSISSET_H_
#define LLH_MULTIANALYSISSET_H_

#include <vector>

#include "llh/public/classes.h"
#include "llh/public/I3Analysis.h"


class MultiAnalysisSet : public AnalysisSet {
 protected:
  vector<AnalysisSet*> aSetVect_;
  EventPtrList evPtrListMulti_;
  EventPtrList evPtrListSrcMulti_;

  // Execute this after actions which change event lists
  void ConstructEventPtrList();
  void ConstructEventPtrListSrc();
  void ReCollocateEvents();
 public:

  MultiAnalysisSet() { }
  ~MultiAnalysisSet() { }

  //  *  CONFIGURATION SETTINGS  *  //

  void AddAnalysisSet(AnalysisSet* aSet);

  //  *  OPERATIONS  *  //

  vector<AnalysisSet*>& GetAnalysisSetVect() { return aSetVect_; }

  double BkgNumberDensity(const Coord& coord) const;

  //There is not one SourceModule defined for the MultiAnalysisSet
  const SourceModule* GetSource() const { return NULL; }

  // THIS WILL HAVE TO BE IMPLEMENTED FOR TIME-DEP CODE!
  EventTimeModule* GetEventTimeModulePtr() { return NULL; }

//  //test
//  EventTimeModule* GetEventTimeModulePtr(){
//
//    vector<I3Event> evTimeModulePtr_;
//    vector<I3Event> tempVect;
//
//    for (int i=0; i<int(analysisFnVect_.size()); ++i) {
//
//        I3Analysis* i3an = dynamic_cast<I3Analysis*>(analysisFnVect_[i]);
//
//        tempVect.clear();
//	    tempVect = i3an->GetEventVector();
//
//	    for (int j=0; j<int(tempVect.size()); j++) {
//	        evTimeModulePtr_.push_back(tempVect[j]);
//	        }
//
//    }
//
//    return evTimeModulePtr_;
//
//  }

//
//  vector<I3Event> MultiBlockAnalysisFn::GetAllEvents() {
//	    vector<I3Event> allEvents;
//	    vector<I3Event> tempVect;
//
//	    for (int i=0; i<int(analysisFnVect_.size()); ++i) {
//
//	        //cout << " i " << i << endl;
//
//	        NewLlhBlockTime* i3an = dynamic_cast<NewLlhBlockTime*>(analysisFnVect_[i]);
//
//	        if (i3an->Get_nEvents()==0) { i3an->PrepareAnalysis(); }
//	        tempVect.clear();
//	        tempVect = i3an->GetEventVector();
//
//	        for (int j=0; j<int(tempVect.size()); j++) {
//	        allEvents.push_back(tempVect[j]);
//	        }
//	    }
//	    return allEvents;
// }




  double GetMeanSrcNev() const;
  double GetMeanSrcNevTime() const;
  double GetMeanSrcNevStacking() const;

  
  //double GetEnergyQuantile(double prob) const;
  
  //double GetSrcEmin(){ return GetEnergyQuantile(0.05);}
  //double GetSrcEmax(){ return GetEnergyQuantile(0.95);}

  //double GetEnergyMinSet(int iset);
  //double GetEnergyMaxSet(int iset);



  double GetMeanSrcNevForFluxModel(const FluxBase& fluxModel) const;
  double GetFluxTimeIntegralSet(int iset);
  double GetFluxEnergyIntegralSet(int iset);

  double GetFluxScaleForNev(double) const;

  void GenerateDataSet();

  void GenerateDataSet_with_nSrcEvents(int nSrcEvents);
  void GenerateDataSetStacking_with_nSrcEvents(int nSrcEvents);

  void UseRealData();

  const EventPtrList* GetEventPtrList() const { return &evPtrListMulti_; }
  const EventPtrList* GetEventPtrListSrc() const { return &evPtrListSrcMulti_; }

  void SetRndOnlyTime(bool flag){
    for (int i=0; i<int(aSetVect_.size()); ++i) {
      (dynamic_cast<I3Analysis*> (aSetVect_[i]))->SetRndOnlyTime(flag);
    } 
  }

  void SetSelEvFileName(char* filename){
    for (int i=0; i<int(aSetVect_.size()); ++i) {
      (dynamic_cast<I3Analysis*> (aSetVect_[i]))->SetSelEvFileName(filename);
    }
  }

  void SetSelEvTreeName(char* treename){
    for (int i=0; i<int(aSetVect_.size()); ++i) {
      (dynamic_cast<I3Analysis*> (aSetVect_[i]))->SetSelEvTreeName(treename);
    }
  }

  void SetTimePdfVect(vector<TimePdf*> tPdfVect, vector<float> tmin, vector<float> tmax){
    const int nSets = aSetVect_.size();
    for(int set=0; set<nSets; ++set) const_cast<I3PointGenerator*>(dynamic_cast<const I3PointGenerator*>(aSetVect_[set]->GetSource()))->SetTimePdfVect(tPdfVect, tmin[set], tmax[set]);
  } 
  
  void PrepareTimeOnlyScramble(int nSrcEvents, vector<int> nSrcPerSet, int iFlare); 
};


#endif // LLH_MULTIANALYSISSET_H_
