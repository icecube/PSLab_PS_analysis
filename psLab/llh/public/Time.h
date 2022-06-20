#ifndef LLH_TIME_H_
#define LLH_TIME_H_


// A simple class, for now

class Time {
private:
 double mjd_;

public:
 Time() : mjd_(0.) { }
 Time(double mjd) : mjd_(mjd) { }

 void SetMJD(double mjd) { mjd_ = mjd; }
 double GetMJD() const { return mjd_; }

 // These functions are for backward compatibility.
 // But the more explicitly-named functions above should be preferred.
 void SetTime(double mjd) { SetMJD(mjd); }
 double GetTime() const { return GetMJD(); }
}; 


// Comparison Function, for e.g. sorting algorithms

bool CompareTime(const Time& t1, const Time& t2) {
  return t1.GetMJD() < t2.GetMJD();
}



#endif // LLH_TIME_H_
