#include "llh/public/EventTimeModuleDiscrete.h"

#include <fstream>

#include "rootExt/public/log_report.h"
#include "rootExt/public/randomfunctions.h"

//
//void EventTimeModuleDiscrete::SetTimesFromMJDFile(string filename, double boundMin, double boundMax) {
//  eventTimeVect_.clear();
//  usedVect_.clear();
//  nUsed_ = 0;
//  timeMin_.SetMJD(1e6); // ridiculously in future, so real min will be found
//  timeMax_.SetMJD(0);   // ridiculously in past, so real max will be found
//
//  cout << "Reading times from file '" << filename << "' ... " << flush;
//  ifstream fin;
//  fin.open(filename.c_str());
//  Time t;
//  double mjd;
//  while ( fin >> mjd ) {
//    t.SetMJD(mjd);
//
//    if (mjd >= boundMin && mjd <= boundMax ){
//        eventTimeVect_.push_back(t);
//        }
//    else{
//        continue;
//        }
//    usedVect_.push_back(0);
//    if ( CompareTime(t, timeMin_) ) { timeMin_ = t; } // i.e.  t < timeMin_
//    if ( CompareTime(timeMax_, t) ) { timeMax_ = t; } // i.e.  timeMax_ < t
//  }
//  fin.close();
//
//  cout << eventTimeVect_.size() << " times read.\n";
//  cout << "Earliest MJD: " << timeMin_.GetMJD();
//  cout << "   Latest MJD: " << timeMax_.GetMJD() << endl;
//}

void EventTimeModuleDiscrete::SetTimesFromMJDFile(TString filename) {
  eventTimeVect_.clear();
  usedVect_.clear();
  nUsed_ = 0;
  timeMin_.SetMJD(1e6); // ridiculously in future, so real min will be found
  timeMax_.SetMJD(0);   // ridiculously in past, so real max will be found

  cout << "Reading times from file '" << filename.Data() << "' ... " << flush;
  ifstream fin;
  fin.open(filename.Data());
  Time t;
  double mjd;
  while ( fin >> mjd ) {
    t.SetMJD(mjd);
    eventTimeVect_.push_back(t);
    usedVect_.push_back(0);
    if ( CompareTime(t, timeMin_) ) { timeMin_ = t; } // i.e.  t < timeMin_
    if ( CompareTime(timeMax_, t) ) { timeMax_ = t; } // i.e.  timeMax_ < t
  }
  fin.close();

  cout << eventTimeVect_.size() << " times read.\n";
  cout << "Earliest MJD: " << timeMin_.GetMJD();
  cout << "   Latest MJD: " << timeMax_.GetMJD() << endl;
}



// This function gets inefficient if almost all event times are used,
// and fails completely if they all get used and another is requested

Time EventTimeModuleDiscrete::GetRandomTime() {
  int total = eventTimeVect_.size();

  //cout << "total: " << total << " nUsed_ " << nUsed_ << endl;

  // Do some routine checks
  if (total==0) { 
      log_fatal("Fatal Error: No event times in list.\n");
  }
  if (nUsed_ > 0.9 * total) { 
    log_warn("WARNING: > 90%% of event times from file have been used!\n");
    if (nUsed_ == total) {
      log_fatal("Fatal Error: new event time requested but list used up.\n");
      exit(1);
    }
  }

  while (1) {
    int ran = (int) (random_uniform(0., total));
    if ( !usedVect_[ran] ) {
      usedVect_[ran] = 1;
      ++nUsed_;
      return eventTimeVect_[ran];
    }
  }
}


void EventTimeModuleDiscrete::ResetUsedTimes() {
  vector<int>::iterator it;
  for (it = usedVect_.begin(); it < usedVect_.end(); ++it)  { *it = 0; }
  nUsed_ = 0;
}
