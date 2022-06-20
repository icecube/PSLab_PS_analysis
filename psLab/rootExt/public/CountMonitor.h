#ifndef ROOTEXT_COUNTMONITOR_H_
#define ROOTEXT_COUNTMONITOR_H_

// Simple tool for keeping track of trials and updating user
// example:
//   CountMonitor countMon(10., 700);
//   for (int i=0; i<700; ++i) {
//     countMon.UpdateCount();
//   }

// Specify frequency of monitoring statements, in percentage, 
// i.e. maximum is 100.  (not 1.)

// if percent==0 or ntotal==0, no monitoring status will be given
// when update is called.


class CountMonitor {

public:

  CountMonitor();
  CountMonitor(double percent, long int ntotal);
  virtual ~CountMonitor() { }

  // This also resets the count
  void SetPercentAndTotal(double percent, long int ntotal);

  // return true if status was printed
  bool UpdateCount();

private:
  double monitorPercent_;
  long int nTotal_;

  long int count_;
  double currentGoalPercent_;
  long int currentGoalCount_;

  void IncrementGoals();
};


#endif // ROOTEXT_COUNTMONITOR_H_
