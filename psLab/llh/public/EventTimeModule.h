#ifndef LLH_EVENTTIMEMODULE_H_
#define LLH_EVENTTIMEMODULE_H_

#include <string>

#include "llh/public/Time.h"
#include "llh/public/I3Event.h"

class EventTimeModule {
 protected:
  Time timeMin_;
  Time timeMax_;
  
 public:
  EventTimeModule(); 
  // default constructor, will set time range of exactly 1 sidereal day
  // for uniform distribution in right ascension

  virtual ~EventTimeModule() { }

  virtual void SetTimesWithinRange(Time timeMin, Time timeMax);

  Time GetTimeMin() const { return timeMin_; }
  Time GetTimeMax() const { return timeMax_; }

  virtual Time GetRandomTime();
  virtual Time GetGaussianTime(double mean, double sigma);

  /*vector<double> GetTimeRange() { 
    vector<double> a;
    a.push_back(timeMin_.GetMJD());
    a.push_back(timeMax_.GetMJD());
    return a;
  }*/

  // This fn only has to do something special for classes with finite # of times
  virtual void ResetUsedTimes() { }
  virtual void SetTimesFromMJDFile(TString filename) { if (filename) { } }
  //virtual void SetTimesFromMJDFile(string filename, double boundMin, double boundMax) { if (filename.length()) { } }

  // for now, at least, doesn't need to be overloaded
  void RandomizeEvent(I3Event& ev, bool rndOnlyTimes=false);
  void RandomizeEventGaussian(I3Event& ev, double mean, double sigma);
};

#endif // LLH_EVENTTIMEMODULE_H_
