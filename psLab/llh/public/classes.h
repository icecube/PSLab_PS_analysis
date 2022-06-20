#ifndef LLH_CLASSES_H_
#define LLH_CLASSES_H_

#include <iostream>
#include <vector>
#include <algorithm>

#include "rootExt/public/log_report.h"

#include "llh/public/CoordClasses.h"
#include "llh/public/TimePdf.h"
#include "llh/public/Time.h"

// Forward Declarations (when feasible, more efficient than including headers)
class AnalysisSet;
class EventTimeModule;
class FluxBase;



class Event {
 protected:
  Time time_;
  AnalysisSet* aSet_;

 public:
  Event() : time_(0), aSet_(NULL) { }
  virtual ~Event() { }

  // allows for a pointer back to the AnalysisSet which event belongs to
  virtual void SetAnalysisSet(AnalysisSet* aSet) { aSet_ = aSet; }
  virtual const AnalysisSet* GetAnalysisSet() const { return aSet_; }

  virtual void SetTime(const Time& t) { time_ = t; }
  virtual Time GetTime() const { return time_; }

  virtual const Coord& GetCoord() const = 0;
  virtual double ProbFrom(const Coord& coord2) const = 0;
};


bool CompareEventTime(const Event& e1, const Event& e2) {
  return CompareTime(e1.GetTime(), e2.GetTime());
}

bool CompareEventPtrTime(const Event* e1, const Event* e2) {
  return CompareTime(e1->GetTime(), e2->GetTime());
}


class EventPtrList {
 protected:
  vector<const Event*> evPtrVect_;

 public:
  EventPtrList() { }
  virtual ~EventPtrList() { }

  virtual void Clear() { evPtrVect_.clear(); }

  virtual void AddEvent(const Event* evPtr) { evPtrVect_.push_back(evPtr); }

  virtual const Event* GetEvent(int i) const { return evPtrVect_[i]; }

  //  virtual Event* ModifyEvent(int i)
  //  {
  //    return evPtrVect_[i];
  //  }

  virtual int GetSize() const { return evPtrVect_.size(); }

  // What the pointers _point_ to stays constant,
  // but the order of the pointers in the vector will be re-arranged by time
  virtual void SortByTime() {
    std::sort( evPtrVect_.begin(), evPtrVect_.end(), CompareEventPtrTime );
  }
};



class SourceModule {
 public:
  virtual ~SourceModule() { }

  virtual SourceModule* Clone() const = 0;
  // This allows us to copy a derived class object w/o knowing what it is.
  // THIS MUST BE DEFINED SEPARATELY for each derived class.
  // For simple class (not requiring a deep copy) this should suffice:
  // { return new DerivedClass(*this); }

  virtual double GetMeanSrcNev() const = 0;

  virtual double GetMeanSrcNevForFluxModel(const FluxBase& fluxModel) const = 0;

  virtual double GetFluxScaleForNev(double) const = 0;
  
  virtual TimePdf * GetTimePdf() const = 0;
  virtual vector<TimePdf*> GetTimePdfVect() const = 0;
  
  virtual void SetTimeAzBins(int nbins) {
    if (nbins) { }
    cout << "setTimeAzBins not implemented!\n";
  }
  
};



class AnalysisSet {
 protected:
  SourceModule *srcModule_;

 public:
  AnalysisSet() : srcModule_(NULL) { }
  virtual ~AnalysisSet() { }

  virtual double BkgNumberDensity(const Coord&) const = 0;

  // Get Signal Info

  virtual const SourceModule* GetSource() const = 0;

  virtual EventTimeModule* GetEventTimeModulePtr() = 0;

  // This functionality of SourceModule is added to AnalysisSet here, but will
  // be over-ridden by MultiAnalysisSet (which doesn't have just one source)
  virtual double GetMeanSrcNev() const { 
    return srcModule_->GetMeanSrcNev();
  }
  virtual double GetMeanSrcNevForFluxModel(const FluxBase& fluxModel) const {
    return srcModule_->GetMeanSrcNevForFluxModel(fluxModel);
  }
  virtual double GetFluxScaleForNev(double nev) const {
    return srcModule_->GetFluxScaleForNev(nev);
  }
   
  // Operations

  virtual void GenerateDataSet() = 0;

  virtual void GenerateDataSet_with_nSrcEvents(int) = 0;

  virtual void UseRealData() = 0;

  virtual const EventPtrList* GetEventPtrList() const = 0 ;
  virtual const EventPtrList* GetEventPtrListSrc() const = 0 ;
};


// Use this if a given function cannot be implemented in a particular
// derived class, e.g. UseRealData():
// {
//   log_fatal("UseRealData() not implemented\n");
//   return -1.;  // need a return statement, even though log_fatal should exit
// }



class AnalysisFn {
 protected:
  AnalysisSet* aSet_;
  const Coord* srcCoord_;

 public:
  AnalysisFn() : aSet_(NULL), srcCoord_(NULL) { }
  virtual ~AnalysisFn() { }

  virtual void SetAnalysisSet(AnalysisSet* aSet) { aSet_ = aSet; }
  virtual void SetSearchCoord(const Coord& coord) { srcCoord_ = &coord; }
  const Coord* GetSearchCoord() { return srcCoord_; }

  virtual void PrepareAnalysis() = 0;
  virtual void MaximizeLlh() = 0;

  virtual double EvalFCN(const vector<double>&) const {
    log_error("EvalFCN Not implemented in base class\n");
    return 1;
  }
  virtual double EvaluateLlh(double *parValueArray) = 0;

  virtual double GetPar(int i) const = 0;
  virtual double Get_logLambdaBest() const = 0;
  virtual double GetTestStatistic() const = 0;
  virtual double GetEstProb() const = 0;
  virtual double GetSigmaMin() { return 0.; };
  virtual double Get_nEvents() const = 0;
};


#endif // LLH_CLASSES_H_



