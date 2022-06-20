#ifndef __CLASSES_SPHERICAL_HEADER_INCLUDED
#define __CLASSES_SPHERICAL_HEADER_INCLUDED

#include "classes.h"
#include "llh_random.h"
#include "TVector3.h"
#include "TMath.h"


class Spherical : public Coord {
  TVector3 v_;
public:
  Spherical () {};

  void SetCoordsPolarRad (double thetaRad, double phiRad) 
  { v_.SetMagThetaPhi(1.,thetaRad,phiRad); }

  double GetThetaRad() const {return v_.Theta(); }
  double GetPhiRad() const {return v_.Phi(); }

  void SetCoordsPolarDeg (double thetaDeg, double phiDeg) 
  { SetCoordsPolarRad(thetaDeg*TMath::DegToRad(), 
		      phiDeg*TMath::DegToRad() ); }
  double GetThetaDeg() const {return v_.Theta()*TMath::RadToDeg(); }
  double GetPhiDeg() const {return v_.Phi()*TMath::RadToDeg(); }


  double DistanceTo (const Coord& coord2) const
  {
    const Spherical* s2 = 
      dynamic_cast<const Spherical*>(&coord2);
    assert(s2);
    return v_.Angle( s2->v_);
  }


  void RotateX(double angleRad) {
    v_.RotateX(angleRad);
  }
  void RotateY(double angleRad) {
    v_.RotateY(angleRad);
  }
  void RotateZ(double angleRad) {
    v_.RotateZ(angleRad);
  }




};





class SphericalEvent : public Event {
  Spherical s_;
  double sigma_;
public:
  SphericalEvent () {}
  SphericalEvent (Spherical s, double sigma)
    { s_ = s;  sigma_ = sigma;}

  Spherical GetCoord() const {return s_;}
  Spherical* ModifyCoord() {return &s_;}

  void SetSigma(double s) {sigma_ = s;}
  double GetSigma() const {return sigma_;}

  double ProbFrom(const Coord &coord2) const {
    double r = s_.DistanceTo(coord2);
    return CircularGaussUnc(r, sigma_);
  }
};



class SphericalSource : public SourceModule {
 public:
  virtual ~SphericalSource() {}
  virtual SphericalEvent GenerateEvent() const =0;
  virtual double GetMeanSrcNev() const = 0;
};




class SphericalGaussianSource : public SphericalSource {
 private:
  Spherical sourceCoord_;
  double genSigma_;
  double recoSigma_;
  double meanSrcNev_;
  // genSigma determines the Gaussian width according to which the
  //   event coordinates will be randomly assigned.
  // recoSigma corresponds to what the event reconstruction later "thinks"
  //   the event uncertainty is.  So, that's what gets assigned to the
  //   event when it is handed back to the analysis code.
  // Ideally, these sigmas would be the same.  In practice, recoSigma could,
  // for example, underestimate the true uncertainty (spread).

 public:

  void SetParams (Spherical s, double genSigma, 
		  double recoSigma, double meanSrcNev)
  { sourceCoord_ = s; 
    genSigma_ = genSigma; 
    recoSigma_ = recoSigma;
    meanSrcNev_ = meanSrcNev; }

  SphericalEvent GenerateEvent() const
  {
    // for small sigma (<5deg), this method should produce reasonable results
    double x = random_gaussian(0., genSigma_);
    double y = random_gaussian(0., genSigma_);

    double r = sqrt(x*x+y*y);
    double ang = atan2(y,x);

    // treat these as polar coordinates of event near pole
    Spherical s;
    s.SetCoordsPolarRad(r,ang);

    // apply rotation which would take the pole to the source location.
    // this moves the event to the vicinity of the source location.

    s.RotateY(sourceCoord_.GetThetaRad());
    s.RotateZ(sourceCoord_.GetPhiRad());

    return SphericalEvent(s, recoSigma_);
  }

  double GetMeanSrcNev() const {return meanSrcNev_;}

};




class SphericalSet : public AnalysisSet {
protected: 
  int nevTotal_;
  vector<SphericalEvent> eVector_;
  double uniformSigma_;
  SphericalGaussianSource Src_;
  double zenMinRad_;
  double zenMaxRad_;
public:

  void SetNevTotal(int ntot) {nevTotal_ = ntot;}

  void SetUniformSigma(double sigma) {uniformSigma_ = sigma;}
  
  void SetSphericalSource(const SphericalGaussianSource& Src) {Src_ = Src;}

  void SetZenithRangeDeg(double zenMinDeg, double zenMaxDeg) {
    zenMinRad_ = zenMinDeg*TMath::DegToRad();
    zenMaxRad_ = zenMaxDeg*TMath::DegToRad();
  }

  double BkgNumberDensity(const Coord& SourceCoord) {
    if (&SourceCoord) // just so no warnings during compile
      ;
    return nevTotal_ / (2.*TMath::Pi()*(cos(zenMinRad_)-cos(zenMaxRad_)));
  }
    

  void GenerateBkgEvents(int Nevents, 
			 vector<SphericalEvent>& BkgEventVector) {
    BkgEventVector.clear();
    for (int i = 0; i<Nevents; i++) {
      double thetaRad = acos(random_uniform(cos(zenMaxRad_),cos(zenMinRad_)));
      double phiRad = random_uniform(0, 2.*TMath::Pi());
      Spherical s;
      s.SetCoordsPolarRad(thetaRad,phiRad);
      BkgEventVector.push_back(SphericalEvent(s,uniformSigma_));
    }
  }

  void GenerateSrcEvents(vector<SphericalEvent>& SrcEventVector) {
    SrcEventVector.clear();
    int Nev = random_poisson( Src_.GetMeanSrcNev() );

    for (int i=0; i<Nev; ++i) {
      SrcEventVector.push_back( Src_.GenerateEvent() );
    }
  }

  void GenerateDataSet() {
    vector<SphericalEvent> SrcEventVector;
    GenerateSrcEvents(SrcEventVector);

    vector<SphericalEvent> BkgEventVector;
    int NevBkg = nevTotal_ - SrcEventVector.size();
    if (NevBkg < 0) {
      log_warn("More Source events than Total... Throwing away extra source events.");
      NevBkg = 0;
    }
    GenerateBkgEvents(NevBkg, BkgEventVector);

    MergeVectorsRandomly<SphericalEvent>
      (SrcEventVector, BkgEventVector, eVector_);
  }



  void FillProbPairVector(vector<ProbPair>& ProbVector,
			  const Coord& SourceCoord) 
  {
    ProbVector.clear();
    for (vector<SphericalEvent>::iterator ev = eVector_.begin();
       ev != eVector_.end(); ++ev)
    {
      double sProb = ev->ProbFrom(SourceCoord);
      double bProb = BkgNumberDensity(ev->GetCoord()) / nevTotal_;
      ProbVector.push_back(ProbPair(sProb,bProb));
    }
  }

  void FillDistanceVector(vector<double>& DistanceVector,
			  const Coord& SourceCoord)
  {
    DistanceVector.clear();
    for (vector<SphericalEvent>::iterator ev = eVector_.begin();
       ev != eVector_.end(); ++ev)
    {
      DistanceVector.push_back( ev->GetCoord().DistanceTo(SourceCoord) );
    }
  }


};







#endif


/*

// First approach, independent of ROOT, but would require implementing
// vector rotations, etc....

class Spherical : public Coord {
  double thetaRad_;
  double phiRad_;
public:
  Spherical () {};
  void SetCoordsPolarRad (double thetaRad, double phiRad) 
  { thetaRad_ = thetaRad; 
    phiRad_ = phiRad; }
  void SetCoordsPolarDeg (double thetaDeg, double phiDeg) 
  { thetaRad_ = thetaDeg*PI/180.; 
    phiRad_ = phiDeg*PI/180.; }

  double GetThetaRad() const {return thetaRad_; }
  double GetPhiRad() const {return phiRad_; }

  double DistanceTo (const Coord& coord2) const
  {
    const Spherical* s2 = 
      dynamic_cast<const Spherical*>(&coord2);
    assert(s2);
    double thetaRad2 = s2->GetThetaRad();
    double phiRad2 = s2->GetPhiRad();
    double dotprod = sin(thetaRad_) * sin(thetaRad2) *
                     cos(phiRad_-phiRad2) +
                     cos(thetaRad_) * cos(thetaRad2);
    return acos(dotprod);
  }
};

*/
