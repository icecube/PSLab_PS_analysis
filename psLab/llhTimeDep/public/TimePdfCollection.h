#ifndef LLHTIMEDEP_TIMEPDFCOLLECTION_H_
#define LLHTIMEDEP_TIMEPDFCOLLECTION_H_

#include <iostream>
#include <string>
#include <vector>

#include "TH1D.h"
#include "TMath.h"

#include "rootExt/public/randomfunctions.h"
#include "llh/public/Time.h"
#include "llh/public/TimePdf.h"


// A simple box-like Pdf.

class BoxTimePdf : public TimePdf {
 private:
	double flarebegin_;
	double flareend_;
	double flareheight_;
	double quietheight_;

 public:

	BoxTimePdf() {
	  tmin_ = 0.;
	  tmax_ = 0.;
	  flarebegin_ = -1.;
	  flareend_ = -1.;
	  flareheight_ = 0.;
	  quietheight_ = 0.;
	  livetimeTotal_ = 0.;
	  norm_ = 0.;
	}
	
    BoxTimePdf(double tmin1, double tmax1, 
                 double flarebegin1, double flareend1, 
                 double flareheight1, double quietheight1) {
                 
	  tmin_ = tmin1;
	  tmax_ = tmax1;  
	  flarebegin_ = flarebegin1;
	  flareend_ = flareend1;
	  flareheight_ = flareheight1;
	  quietheight_ = quietheight1;
	  livetimeTotal_ = tmax1 - tmin1;
	  //norm_ = (flareend1-flarebegin1)*flareheight1 + 
	         //(tmax1-tmin1-(flareend1-flarebegin1))*quietheight1;
          SetNorm();
	}
	
    TimePdf* Clone() const { return new BoxTimePdf(*this); }

	double GetPdfValue(double x) { // don't forget to normalize!
        if (x<=flareend_ && x>=flarebegin_) return flareheight_/norm_;
  	else if(norm_==0) return 0;
        else return quietheight_/norm_;
	}
	
//  double GetPdfValue(Time t)             { return GetPdfValue( t.GetMJD() ); }
//  double GetPdfValue(Time t, double lag) { return GetPdfValue( t.GetMJD()-lag ); }
	
    void fillHisto(int nbins){
      int nbins_ = nbins;
      double bintime = (tmax_ - tmin_) / nbins_;
      double MJD;
      histo_ = new TH1D("histo_","source",nbins_,tmin_,tmax_);
      for (int i=1; i <= nbins_; i++){ // bins go from 1 to n
        MJD = tmin_ + (i-0.5)*bintime;
        histo_->SetBinContent(i,GetPdfValue(MJD));
      }
    }
    
    TH1D * GetHisto() {return histo_;}
	
	Time GenerateEventTime() {
	  double t;
	  if ( quietheight_< 1e-9) {
	    do{
              t = random_uniform(flarebegin_,flareend_);
            } while (t < tmin_ || t > tmax_);
	  } else {
	    t = histo_->GetRandom();
	  }
	    
        return Time(t);
    }
    
  void SetNorm() {
      //Printf("tmin = %.1f, tmax = %.1f, flarebegin = %.1f, flareend = %.1f", tmin_, tmax_, flarebegin_, flareend_);
      if(flarebegin_>tmax_ || flareend_<tmin_) norm_ = (tmax_-tmin_)*quietheight_;
      else{
          double upperTime, lowerTime;
          if(flareend_>tmax_) upperTime = tmax_;
          else upperTime = flareend_;
   
          if(flarebegin_<tmin_) lowerTime = tmin_;
          else lowerTime = flarebegin_;
          norm_ = ( upperTime-lowerTime)*flareheight_ + 
	                         ( tmax_-tmin_-(upperTime-lowerTime))*quietheight_;
      }
  }   
    
  void CheckTimeBounds(double dataTmin, double dataTmax) {
 
    if ( (flareend_ < dataTmin && flarebegin_ < dataTmin) || 
 	     (flareend_ > dataTmax && flarebegin_ > dataTmax)    ) {
	  cout << "Your entire flare seems to be out of the data period!" << endl;
	}
 
    if (tmin_ < dataTmin)        { tmin_       = dataTmin; }
    //if (flarebegin_ < dataTmin)  { flarebegin_ = dataTmin; }
    //if (flareend_ > dataTmax)    { flareend_   = dataTmax; }
    if (tmax_ > dataTmax)        { tmax_       = dataTmax; }
	
    if (flareend_ < flarebegin_) { 
      double a = flarebegin_;
      flarebegin_ = flareend_;
      flareend_ = a;
    }

    livetimeTotal_ = tmax_ - tmin_;
    SetNorm();
  }
};

// A simple 1-D Gaussian time Pdf

class GaussianTimePdf : public TimePdf {
 private:
  double mu_;
  double sigma_;
  double drawnorm_;
  double nSigTrunc_;

 public:

  GaussianTimePdf() {
    tmin_ = 0.;
    tmax_ = 0.;
    mu_ = 0.;
    sigma_ = 1.;
    livetimeTotal_ = 1.;
    norm_ = 1.;
    setNormFromGRL_ = false;
    nSigTrunc_ = 4;
    histo_ = NULL;
  }
    
  GaussianTimePdf(double tmin, double tmax, double mu, double sigma, double drawnorm) {
    tmin_ = tmin;
    tmax_ = tmax;
    mu_ = mu;
    sigma_ = sigma;
    livetimeTotal_ = tmax_ - tmin_;
    drawnorm_ = drawnorm;
    setNormFromGRL_ = false;
    nSigTrunc_ = 4;
    histo_ = NULL;
    SetNorm();
  }

  GaussianTimePdf(double tmin, double tmax, double mu, double sigma, double drawnorm, double nsigtrunc) {
    tmin_ = tmin;
    tmax_ = tmax;
    mu_ = mu;
    sigma_ = sigma;
    livetimeTotal_ = tmax_ - tmin_;
    drawnorm_ = drawnorm;
    setNormFromGRL_ = false;
    nSigTrunc_ = nsigtrunc;
    histo_ = NULL;
    SetNorm();
  }
  
  TimePdf* Clone() const { return new GaussianTimePdf(*this); }
  
  // For the GaussianTimePdf we can just ust the analytical function.
  
  double GetPdfValue(double x) {
    if (x<tmin_ || x>tmax_) return 0;
    if(x<mu_-nSigTrunc_*sigma_ || x>mu_+nSigTrunc_*sigma_) return 0;

    double a=(1./norm_)*(1./(TMath::Sqrt(2.*TMath::Pi())*sigma_))*exp((-(x-mu_)*(x-mu_))/(2*sigma_*sigma_));
     
    return a;
  }


  Time GenerateEventTime(){
    double t = random_gaussian(mu_, sigma_);
    while (t < tmin_ || t > tmax_) {t = random_gaussian(mu_, sigma_);}
    return Time(t);
  }

  void fillHisto(int nbins){ // this is really just for plotting the pdf after the face
     /*int nbins_ = nbins;
     double bintime = (tmax_ - tmin_) / nbins_;
     double MJD;
     if(!histo_) histo_ = new TH1D(Form("histo_%f", mu_),"source",nbins_,tmin_,tmax_);
     for (int i=1; i <= nbins_; i++){ // bins go from 1 to n
       MJD = tmin_ + (1.0*i-0.5)*bintime;
       histo_->SetBinContent(i,GetPdfValue(MJD)*drawnorm_);
     }*/
  }

  TH1D * GetHisto() {return histo_;}
   
  void SetNorm() {
    double zlow  = TMath::Erf( (tmin_-mu_)/(sigma_*TMath::Sqrt(2.)) );
    double zhigh = TMath::Erf( (tmax_-mu_)/(sigma_*TMath::Sqrt(2.)) );
    norm_  = (zhigh-zlow)/2.;
    if(setNormFromGRL_){
      for(unsigned int i=0; i<startMissRunVect_.size(); ++i){
        if(startMissRunVect_[i] > mu_-nSigTrunc_*sigma_ && stopMissRunVect_[i] < mu_+nSigTrunc_*sigma_){
          zlow = TMath::Erf( (startMissRunVect_[i]-mu_)/(sigma_*TMath::Sqrt(2.)) );
          zhigh = TMath::Erf( (stopMissRunVect_[i]-mu_)/(sigma_*TMath::Sqrt(2.)) );
          norm_ -= (zhigh-zlow)/2.;
        }
      }
    }
    if(norm_ < 0) Printf("ERROR: normalization of time PDF is negative (%f). Please check your code", norm_);
  }

  void SetNormFromGRL(vector<double> start, vector<double> stop) {
    setNormFromGRL_ = true;
    
    startMissRunVect_.clear();
    stopMissRunVect_.clear();
    
    for(unsigned int i=0; i<start.size(); ++i){
      startMissRunVect_.push_back(start[i]);
      stopMissRunVect_.push_back(stop[i]);
    }
    SetNorm();
  }

  void CheckTimeBounds(double dataTmin, double dataTmax) { // this is troublesome if applied
    if (tmin_ < dataTmin) tmin_ = dataTmin;                // to the same pdf multpile times
    if (tmax_ > dataTmax) tmax_ = dataTmax;
    livetimeTotal_ = tmax_ - tmin_;

    SetNorm();
  }

  void SetNSigTrunc(double nSig){ nSigTrunc_ = nSig; }
  void SetMean(double mean){ mu_ = mean; }
 
};

// Two Gaussians with some sort of offset (both with same normalization)
// to test what time structure can look like.

class DoubleGaussian : public TimePdf {
  private:   
    double mu_;
    double sigma_;
    double sep_;
    double drawnorm_;

 public:
  DoubleGaussian() {
    tmin_ = 0.;
    tmax_ = 0.;
    mu_ = 0.;
    sigma_ = 1.;
    sep_ = 2.;
    livetimeTotal_ = 1.;
    norm_ = 1.;
  }
    
  DoubleGaussian(double tmin, double tmax, double mu, double sigma, double sep, double dnorm) {
    tmin_ = tmin;
    tmax_ = tmax;
    mu_ = mu;
    sigma_ = sigma;
    sep_ = sep;
    livetimeTotal_ = tmax - tmin;
    drawnorm_ = dnorm;
  }
  
  TimePdf* Clone() const { return new DoubleGaussian(*this); }
  
  double GetPdfValue(double x) {
    if (x<tmin_ ||x>tmax_) return 0;
    double a = max( (0.3989/sigma_)*exp((-(x-mu_)*(x-mu_))/(2*sigma_*sigma_)), 
                    (0.3989/sigma_)*exp((-(x-mu_-sep_)*(x-mu_-sep_))/(2*sigma_*sigma_)) );
    if (a < .0000000001 ) a=0;
    return a/2.;
  }

//  double GetPdfValue(Time t)             { return GetPdfValue( t.GetMJD() ); }
//  double GetPdfValue(Time t, double lag) { return GetPdfValue( t.GetMJD()-lag ); }

  Time GenerateEventTime(){
    double a = random_uniform(0.,2.);
    double t=0.;
    if (a < 1) {
      t = random_gaussian(mu_, sigma_);
    } else {
      t = random_gaussian(mu_+sep_, sigma_);
    }
    return Time(t);
  }

   void fillHisto(int nbins){
    int nbins_=nbins;
    double bintime = (tmax_ - tmin_) / nbins_;
    double MJD;
    histo_ = new TH1D("histo_","source",nbins_,tmin_,tmax_);
    for (int i=1; i <= nbins_; i++){
      MJD = tmin_ + (1.0*i-0.5)*bintime;
      histo_->SetBinContent(i,GetPdfValue(MJD));
    }
   }

   TH1D * GetHisto() {return histo_;}   
   void SetNorm() {norm_ = 1.;}

};

// A special type of Gaussian which wraps around from tmin_ to tmax_.
// It only takes the most probable value of the central Gaussian.

class PeriodicGaussianTimePdf : public TimePdf {
  private:
    double mu_;
    double sigma_;

 public:
  PeriodicGaussianTimePdf() {
    tmin_ = 0.;
    tmax_ = 1.;
    mu_ = 0.;
    sigma_ = 1.;
    livetimeTotal_ = 1.;
    norm_ = 1.;
  }
    
  PeriodicGaussianTimePdf(double tmin, double tmax, double mu, double sigma, double norm) {
    tmin_ = tmin;
    tmax_ = tmax;
    mu_ = mu;
    sigma_ = sigma;
    livetimeTotal_ = tmax_ - tmin_;
    norm_ = norm;
  }
  
  TimePdf* Clone() const { return new PeriodicGaussianTimePdf(*this); }
  
  void SetNorm() { // the normalization is different here than in the 
                   // Gaussian since we only look at 0 to 1 in phase.
    norm_ = 1.;//TMath::Erf((1.-mu_)/sqrt(2.)*sigma_) - TMath::Erf((0.-mu_)/sqrt(2.)*sigma_);
  }

  
  double GetPdfValue(double x) {
    //if (x<tmin ||x>tmax) { return 0; }
    while ( x-mu_ > (tmax_-tmin_)/2. ) { x -= (tmax_-tmin_); }
    while ( x-mu_ < -(tmax_-tmin_)/2. ) { x += (tmax_-tmin_); }
    double a=(0.3989/sigma_)*exp((-(x-mu_)*(x-mu_))/(2.*sigma_*sigma_));
    return a;
  }

//  double GetPdfValue(Time t)             { return GetPdfValue( t.GetMJD() ); }
//  double GetPdfValue(Time t, double lag) { return GetPdfValue( t.GetMJD()-lag ); }

  Time GenerateEventTime(){
    double t = random_gaussian(mu_, sigma_);
    while ( t > tmax_ ) { t -= tmax_-tmin_; }
    while ( t < tmin_ ) { t += tmax_-tmin_; }
    return Time(t);
  }

  void fillHisto(int nbins){
    int nbins_ = nbins;
    double bintime = (tmax_ - tmin_) / nbins_;
    double MJD;
    histo_ = new TH1D("histo_","source",nbins_,tmin_,tmax_);
    for (int i=1; i <= nbins_; i++){ // bins go from 1 to n
      MJD = tmin_ + (1.0*i-0.5)*bintime;
      histo_->SetBinContent(i,GetPdfValue(MJD)*norm_);
    }
  }

  TH1D * GetHisto() {return histo_;}

  void CheckTimeBounds(double dataTmin, double dataTmax) { 
    //cout << "Not used for periodic searches!" << endl;
    if (dataTmin || dataTmax) { }
  }

};


class PeriodicMultiGaussianTimePdf : public TimePdf {
  private:
    double mu_;
    double sigma_;
    double periodobserved_;
    double periodactual_;
    double nofperiods_;
    double perrorshift_;
    

 public:
  PeriodicMultiGaussianTimePdf() {
    tmin_ = 0.;
    tmax_ = 1.;
    mu_ = 0.;
    sigma_ = 1.;
    livetimeTotal_ = 1.;
    norm_ = 1.;
  }
    
  PeriodicMultiGaussianTimePdf(double tmin, double tmax, double mu, double sigma, double norm, double livetimeTotal, double periodactual, double periodobserved) {
    
    
    
    
    tmin_ = tmin;
    tmax_ = tmax;
    mu_ = mu;
    sigma_ = sigma;
    livetimeTotal_ = livetimeTotal;
    norm_ = norm;
    periodobserved_ = periodobserved;
    periodactual_ = periodactual;
    nofperiods_ = livetimeTotal_/periodobserved_;
    perrorshift_ = (periodobserved_-periodactual_)/periodobserved_;
    fillHisto(10000);
    
    
  }
  
  TimePdf* Clone() const { return new PeriodicMultiGaussianTimePdf(*this); }
  
  void SetNorm() { // the normalization is different here than in the 
                   // Gaussian since we only look at 0 to 1 in phase.
    norm_ = 1.;//TMath::Erf((1.-mu_)/sqrt(2.)*sigma_) - TMath::Erf((0.-mu_)/sqrt(2.)*sigma_);
  }

  
  double GetPdfValueinst(double x, double mu) {
    //if (x<tmin ||x>tmax) { return 0; }
    while ( x-mu > (tmax_-tmin_)/2. ) { x -= (tmax_-tmin_); }
    while ( x-mu < -(tmax_-tmin_)/2. ) { x += (tmax_-tmin_); }
    double a=(0.3989/sigma_)*exp((-(x-mu)*(x-mu))/(2.*sigma_*sigma_));
    return a;
  }
  
  double GetPdfValue(double x) {
    //if (x<tmin ||x>tmax) { return 0; }
    double a = 0;
    double muinst_ = mu_ - perrorshift_*nofperiods_/2;
    
    for(int i=0; i<nofperiods_;i++)
    {
      a += GetPdfValueinst(x, muinst_)/nofperiods_;
      muinst_ += perrorshift_;
    }
    return a;
  }  
//  double GetPdfValue(Time t)             { return GetPdfValue( t.GetMJD() ); }
//  double GetPdfValue(Time t, double lag) { return GetPdfValue( t.GetMJD()-lag ); }

  Time GenerateEventTime(){
    
    

    
     double t=histo_->GetRandom();
     
     
     
//     double muinst_ = mu_ - perrorshift_*nofperiods_/2;
//     
//     double r=0;
//     
//     for(int i=0; i<nofperiods_;i++)
//     {
//       r = random_gaussian(muinst_, sigma_)/nofperiods_;
//       while ( r > tmax_ ) { r -= tmax_-tmin_; }
//       while ( r < tmin_ ) { r += tmax_-tmin_; }
//       muinst_ += perrorshift_;
//       t += r;
//     }
     
     while ( t > tmax_ ) { t -= tmax_-tmin_; }
     while ( t < tmin_ ) { t += tmax_-tmin_; }
    return Time(t);
  }

  void fillHisto(int nbins){
    int nbins_=nbins;
    double bintime = (tmax_ - tmin_) / nbins_;
    double MJD;
    histo_ = new TH1D("histo_","source",nbins_,tmin_,tmax_);
    for (int i=1; i <= nbins_; i++){ // bins go from 1 to n
      MJD = tmin_ + (1.0*i-0.5)*bintime;
      histo_->SetBinContent(i,GetPdfValue(MJD)*norm_);
    }
  }

  TH1D * GetHisto() {return histo_;}

  void CheckTimeBounds(double dataTmin, double dataTmax) { 
    //cout << "Not used for periodic searches!" << endl;
    if (dataTmin || dataTmax) { }
  }

};


// This guy is a sine curve with phase to set. I only used it once.

class PeriodicTimePdf : public TimePdf {
  private:
    double period_;
    double phase_;

 public:
  PeriodicTimePdf() {
    norm_   = 1.;
    period_ = 0.;
    phase_  = 0.;
    tmax_   = 0.;
    tmin_   = 0.;
    livetimeTotal_ = 1.;
  }
  
  PeriodicTimePdf(double tmin, double tmax, double period, double phase) {
    period_ = period;
    phase_  = phase;
    tmax_   = tmax;
    tmin_   = tmin;
    livetimeTotal_ = tmax_ - tmin_;
  }
  
  TimePdf* Clone() const { return new PeriodicTimePdf(*this); }
  
  void SetNorm() {
    norm_ = tmax_ - tmin_ + cos(2*acos(-1)*(tmax_)/period_ + phase_) - cos(2*acos(-1)*(tmin_)/period_ + phase_);
  } 

  void fillHisto(int nbins){
    int nbins_=nbins;
    double bintime = (tmax_ - tmin_) / nbins_;
    double MJD;
    histo_ = new TH1D("histo_","source",nbins_,tmin_,tmax_);
    for (int i=1; i <= nbins_; i++){ // bins go from 1 to n
      MJD = tmin_ + (i*1.0-0.5)*bintime;
      histo_->SetBinContent(i,GetPdfValue(MJD));
    }
  }
  
  TH1D * GetHisto() {return histo_;}
  
  Time GenerateEventTime(){
    double t = histo_->GetRandom();
    return Time(t);
  }

  double GetPdfValue(double x) {
    if (x<tmin_ ||x>tmax_) return 0;
    return (cos(2*acos(-1)*x/period_ - 2*acos(-1)*phase_) + 1)/norm_;
  }
  
//  double GetPdfValue(Time t)             { return GetPdfValue( t.GetMJD() ); }
//  double GetPdfValue(Time t, double lag) { return GetPdfValue( t.GetMJD()-lag ); }

  void CheckTimeBounds(double dataTmin, double dataTmax) {
    //cout << "Not used for periodic searches!" << endl;
    if (dataTmin || dataTmax) { }
  }

};



#endif // LLHTIMEDEP_TIMEPDFCOLLECTION_H_
