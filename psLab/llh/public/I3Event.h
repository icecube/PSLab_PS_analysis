#ifndef LLH_I3EVENT_H_
#define LLH_I3EVENT_H_


#include "llh/public/classes.h"
#include "llh/public/CoordEquatorialDeg.h"
#include "llh/public/Time.h"


// Forward Declarations (when feasible, more efficient than including headers)
class EnergyProb;
class BkgSpaceProb;

class I3EventParameters {
 public:
  double recoZenithDeg;
  double recoAzimuthDeg;
  double parafitSigmaDeg;
  //  int parafitNdir;  // haven't used this since IC9
  double energyValue;  // this can be Nchan, or another direct estimator
  int runID;
  int eventID;
};



class I3MCParameters {
 public:
  double mcEnergy;
  double PS_FlatSpectrumRate;
  // The equivalent rate (Hz) for this event, assuming flat E spectrum
  // (i.e.: for  dPhi/dE = 1 / (GeV cm^2 s)   -- constant for all E.
  //  The purpose is to make it easy to convert to the desired flux later.)

  double srcWeight;
  // The srcWeight is determined by scaling the PS_FlatSpectrumRate_:
  //
  // PS_FlatSpectrumRate_ * flux * livetime * pow(mcEnergy_ , spectralIndex)
  //
  // Note only the pow() term changes the relative weighting of events,
  // but including the flux and livetime means that the Sum of all srcWeight 
  // is in fact the actual number of events expected.
};

class I3PointSpreadFunction {
 public:
  virtual double ProbFrom(double r, double sigmaParam) const = 0;
  virtual ~I3PointSpreadFunction() = 0;
};
inline I3PointSpreadFunction::~I3PointSpreadFunction() { }


class PSFCircularGaussUnc : public I3PointSpreadFunction {
 public:
  virtual double ProbFrom(double r, double sigma) const {
    return CircularGaussUnc(r, sigma);
  }
};

PSFCircularGaussUnc psfCGU;  // just need one instance

I3PointSpreadFunction *defaultI3PSFptr = &psfCGU;


class I3Event : public Event {
 protected:
  EquatorialDeg eq_;
  I3EventParameters params_;
  I3MCParameters mcParams_;
  double sigma_;

  BkgSpaceProb* bkgSpaceProb_;
  EnergyProb* eProb_;

  const I3PointSpreadFunction *psf_;

 public:
  I3Event() { psf_ = defaultI3PSFptr; }

  // BASE-CLASS FUNCTIONS

  // for generic access to coord, consistent with base class Event
  const Coord& GetCoord() const  { return eq_; }
  // for direct access to EquatorialDeg, if you know you have an I3Event
  EquatorialDeg GetEquatorialDeg() const { return eq_; }

  double ProbFrom(const Coord &coord2) const {
    double rDeg = eq_.DistanceTo(coord2);
    return psf_->ProbFrom(rDeg, sigma_);
  }

  // FUNCTIONS SPECIFIC FOR I3EVENT

  virtual void SetBkgSpaceProb(BkgSpaceProb* bsp) { bkgSpaceProb_ = bsp; }
  virtual const BkgSpaceProb* GetBkgSpaceProbFn() const { 
    return bkgSpaceProb_; }

  void SetEnergyProb(EnergyProb* ep) { eProb_ = ep; }
  const EnergyProb* GetEnergyProbFn() const { return eProb_; }
  void SetEnergy(double energy) { params_.energyValue = energy; }

  void SetPointSpreadFunction(const I3PointSpreadFunction *psf) { psf_ = psf; }
  const I3PointSpreadFunction* GetPointSpreadFunction() { return psf_; }

  void SetCoord(const EquatorialDeg& eq) {eq_ = eq;}
  void SetRaDeg(double raDeg)   { eq_.SetRaDeg(  raDeg); }
  void SetDecDeg(double decDeg) { eq_.SetDecDeg(decDeg); }

  void SetParams(I3EventParameters p) 
  { params_ = p;
    sigma_ = params_.parafitSigmaDeg;
  }
  void SetAlternateSigma(double sigma)
  { sigma_ =sigma;}
  void SetMCParams(I3MCParameters mp) {mcParams_ = mp;}

  I3EventParameters GetParams() const {return params_;}
  I3MCParameters GetMCParams() const {return mcParams_;}

  void SetMJD(double mjd) { time_.SetMJD(mjd); }
  double GetMJD() const { return time_.GetMJD(); }

  double GetSigma() const { return sigma_; }

};

#endif // LLH_I3EVENT_H_

