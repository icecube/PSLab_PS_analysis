#include "rootExt/public/CountMonitor.h"

#include <iostream>


CountMonitor::CountMonitor() : 
  monitorPercent_(10.), nTotal_(0), count_(0),
  currentGoalPercent_(0.), currentGoalCount_(0) 	       
{ }

CountMonitor::CountMonitor(double percent, long int ntotal) {
  SetPercentAndTotal(percent, ntotal);
}

void CountMonitor::SetPercentAndTotal(double percent, long int ntotal) {
  monitorPercent_ = percent; 
  nTotal_ = ntotal;
  currentGoalPercent_ = 0.;
  currentGoalCount_ = 0;
  count_ = 0;
  IncrementGoals();
}

// return true if status was printed
bool CountMonitor::UpdateCount() {
  // only give output if monitorPercent_ && nTotal_ are non-zero
  if (!monitorPercent_ || !nTotal_) { return false; }

  ++count_;
  // now check if we've reached the next goal
  if (count_ >= currentGoalCount_) {
    cout << currentGoalPercent_ << "% " << flush;
    if (count_ == nTotal_) {
      cout << endl; // done
    }

    IncrementGoals();
    return true;  // yes, there was output
  }
  
  return false;  // there was no output this time
}


void CountMonitor::IncrementGoals() {
  while (count_ >= currentGoalCount_) {
    currentGoalPercent_ += monitorPercent_;
    currentGoalCount_ = int(nTotal_*currentGoalPercent_/100.);
  }
  // this loop will repeat if there are more percentage intervals
  // than trials, so that the percentages can catch up with the trials

  // now, if goal is too high, set to final goal
  if (currentGoalCount_ > nTotal_) { 
      currentGoalPercent_ = 100.;
      currentGoalCount_ = nTotal_; 
  }
}
