#ifndef LLH_EVENTTIMEMODULEDISCRETE_H_
#define LLH_EVENTTIMEMODULEDISCRETE_H_

//#include <string>
#include <vector>

#include "llh/public/EventTimeModule.h"


// Make a (large) list of possible event times using a text file list

class EventTimeModuleDiscrete : public EventTimeModule {
 protected:
  vector<Time> eventTimeVect_;
  vector<int> usedVect_;
  int nUsed_;

 public:
  EventTimeModuleDiscrete() {
    SetTimesWithinRange(Time(0), Time(-1));
    // over-ride base class settings with nonsense settings, so that
    // function won't work until file is loaded
  }

  /*// over-ride this function... it could be implemented later if need arises
  virtual void SetTimesWithinRange(Time, Time) {
    log_fatal("EventTimesModuleDiscrete::SetTimesWithinRange function\n"
	      "is not defined. (It could be, but one would have to think\n"
	      "about what exactly it would mean and do.\n");
    exit(1);
  } //*/

  void SetTimesFromMJDFile(TString filename);
  //void SetTimesFromMJDFile(string filename, double boundMin, double boundMax);

  virtual Time GetRandomTime();

  virtual void ResetUsedTimes();
};

#endif // LLH_EVENTTIMEMODULEDISCRETE_H_

