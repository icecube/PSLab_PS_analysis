#ifndef LLH_SNEVENT_H_
#define LLH_SNEVENT_H_

#include <iostream>
#include <fstream>
#include "llh/public/classes.h"
#include "llh/public/CoordEquatorialDeg.h"
#include "llh/public/Time.h"

//#include "llh/public/SimpleProb.h"

// Forward Declarations (when feasible, more efficient than including headers)
class SimpleProb;
class SNAnalysis;


class SNEvent : public Event {
 protected:
  Time time_;
//  AnalysisSet* aSet_;

//  SNAnalysis* snSet_;
 
  SimpleProb* sProb_;
  double strength_;
  double nch_; //# of live DOMs for trigger
  double binsize_; //in seconds
  
  EquatorialDeg eq_; //grumble

 public:
//  SNEvent() : time_(0.), snSet_(NULL), sProb_(NULL) { }
  SNEvent() : sProb_(NULL) { }
  
  // FUNCTIONS SPECIFIC FOR SNEVENT

  const Coord& GetCoord() const { return eq_; }
  double ProbFrom(const Coord& coord2) const { return 1.; }
  
  void SetSimpleProb(SimpleProb* ep) { sProb_ = ep; }
  const SimpleProb* GetSimpleProbFn() const { return sProb_; }

  void SetMJD(double mjd) { time_.SetMJD(mjd); }
  double GetMJD() const { return time_.GetMJD(); }
  
  void SetParams(double strength, double nch, double binsize) {
    strength_ = strength;
    nch_ = nch;
    binsize_ = binsize;
  }
  
  double GetStrength() const { return strength_; }
  double GetNch() const { return nch_; }
  double GetBinsize() const { return binsize_; }

//  void SetAnalysisSet(SNAnalysis* snSet) { snSet_ = snSet; }
//  const SNAnalysis* GetAnalysis() const { return snSet_; }

//  void SetTime(const Time& t) { time_ = t; }
//  Time GetTime() const { return time_; }
  
};

class SNEventLoader {
  private:
    vector<SNEvent> baseVector_;
  
  public:
    SNEventLoader() { baseVector_.clear(); }

  void LoadBkgEvents(string filename, vector<SNEvent>& bkgEvents) {     
    bkgEvents.clear();
    cout << "Reading events from file " << filename.c_str() << " ... " << flush;
    ifstream fin;
    fin.open(filename.c_str());
    double mjd, strength, binsize, nch;
    
    SNEvent ev;
    Time t;
    while ( fin >> mjd ) {
      fin >> strength >> nch >> binsize;
      t.SetMJD(mjd);
      ev.SetParams(strength, nch, binsize);
      ev.SetTime(t);
      
      bkgEvents.push_back(ev);
    }
    fin.close();
    cout << "found " << bkgEvents.size() << " events." << endl;
  }
  
};

#endif // LLH_SNEVENT_H_

