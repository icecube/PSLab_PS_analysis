#ifndef LLH_TIMEPDF_H_
#define LLH_TIMEPDF_H_

#include <string>

#include "TH1D.h"

#include "llh/public/Time.h"


class TimePdf {
 protected:
   double norm_;
   double tmin_;
   double tmax_;
   double livetimeTotal_;
   TH1D * histo_;
   bool setNormFromGRL_;
   vector<double> startMissRunVect_;
   vector<double> stopMissRunVect_;

 public:
 
  virtual ~TimePdf() {}

  virtual TimePdf* Clone() const = 0;
  // This allows us to copy a derived class object w/o knowing what it is.
  // THIS MUST BE DEFINED SEPARATELY for each derived class.
  // For simple class (not requiring a deep copy) this should suffice:
  // { return new DerivedClass(*this); }
   
  virtual double GetPdfValue(double) = 0;  
  
  double GetPdfValue(Time t)             { return GetPdfValue( t.GetMJD() ); }
  double GetPdfValue(Time t, double lag) { return GetPdfValue( t.GetMJD()-lag ); }

  virtual void fillHisto(int nbins) = 0;
  virtual TH1D * GetHisto() = 0;
    
   virtual Time GenerateEventTime() = 0;
   
   double GetTmin() { return tmin_; }
   double GetTmax() { return tmax_; }
   double GetNorm() {return norm_;}
   double GetLivetime() { return livetimeTotal_; }
   vector<double> GetStartMissRunVect() { return startMissRunVect_; }
   vector<double> GetStopMissRunVect()  { return stopMissRunVect_; }
   
   virtual void SetNorm() = 0; //for periodic tPdf.
   virtual void SetNormFromGRL(vector<double> start, vector<double> stop) {
     if(!start.empty() || !stop.empty()) { }
   };

   virtual void SetBlockLevels(string fname, double cut=0., double offset=0.) {
     //     cout << "No Blocks Used Here!\n" << fname << endl;
     if (fname.length() || cut || offset) { }
   }
   
   void SetLivetime() { //for some ungodly reason this is giving me ~1e-315.
     livetimeTotal_ = tmax_ - tmin_;
   }
   
   void SetLivetime(double x) { //so let's try overloading it...
     livetimeTotal_ = x;
   }
   
   virtual void CheckTimeBounds(double dataTmin, double dataTmax) {
     if (dataTmin || dataTmax) { }
   }
    
};

#endif // LLH_TIMEPDF_H_

