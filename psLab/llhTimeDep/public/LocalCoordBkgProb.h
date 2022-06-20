#ifndef LLHTIMEDEP_LOCALCOORDBKGPROB_H_
#define LLHTIMEDEP_LOCALCOORDBKGPROB_H_

#include "TH2D.h"
#include "TGraph.h"
#include "llh/public/classes.h"
#include "llh/public/BkgSpaceProb.h"
#include "llh/public/DecBkgProb.h"
#include "llh/public/I3Event.h"

//class LocalCoordBkgProb : public BkgSpaceProb {
class LocalCoordBkgProb {

private:
  //DecBkgProb *decProb_;
  
 public:
  TH2D dataAZ;
  TH2D pdfAZ;
  TGraph *pdfAZ_1D;
  
  int nbAz;
  int nbZn;
  double livetimeTotal_;
  
  double tMin;
  double tMax;
  int nBinsT;
  bool useCosZ;
  
  TH1D dataZT;
  TH1D pdfZT;

  LocalCoordBkgProb() {
    nBinsT = 0;
    nbAz = 1;
    nbZn = 1;
  }
  
//  void Initialize(int nBinsDec, double sigmaSmooth, int nBinsZn=12, int nBinsAz=90, bool uCz=true) {
  void Initialize(int nBinsZn=12, int nBinsAz=90, bool uCz=true) {
    //decProb_->Initialize(nBinsDec, sigmaSmooth);
    nbAz = nBinsAz;
    nbZn = nBinsZn;
    useCosZ = uCz;
  }
  ~LocalCoordBkgProb() { }

  // BASE-CLASS FUNCTIONS

  LocalCoordBkgProb* Clone() const { return new LocalCoordBkgProb(*this); }
  // This allows us to copy a derived class object w/o knowing what it is.
  // THIS MUST BE DEFINED SEPARATELY for each derived class.
  // For simple class (not requiring a deep copy) this should suffice:
  // { return new DerivedClass(*this); }

  /*void FixToBase() { decProb_->FixToBase(); }
  void FixToBasePlusEvents(const EventPtrList& evList) {
    decProb_->FixToBasePlusEvents(evList);
  }

  void SetDecProb(DecBkgProb d) {
    decProb_ = &d;
  } */

  double BackgroundLCProb(const I3Event& ev);
  double BackgroundLCProb_folded(const I3Event& ev);

  double BackgroundLCProbA_mcozZ(double azimuth,double zenith);
  double BackgroundLCProbA_mcozZ_folded(double azimuth,double zenith);

  double BackgroundLCProb(const Event& event) const {
    cout << "is the Event cat breaking?" << endl;
    double blc = BackgroundLCProb( dynamic_cast<const I3Event&>(event) );
    return blc;
  }

  //double GetBkgProbDensity(const Coord& coord) const {
  //  return decProb_->GetBkgProbDensity(coord);
  //}

  // SPECIFIC FUNCTIONS FOR DecBkgProb CLASS

//  void Initialize(int nMaps, double sigmaSmooth);

  void FillLCBkgHisto(const vector<I3Event>& events);
  void FillLCBkgHisto_folded(const vector<I3Event>& events,double minDegZenith=0.,double maxDegZenith=180.);
  void SetTimeParams(double min, double max, int nbins=40);

  /*
  void SetBaseDecMap(const vector<I3Event>& eVect) {
    //decProb_->SetBaseDecMap(eVect);
    FillLCBkgHisto(eVect);
  }

  double GetBkgProbDensity(const Event& event) const {
    const I3Event i3ev = dynamic_cast<const I3Event&>(event);
    double lcprob = BackgroundLCProb( i3ev );
  
    return lcprob;// * decProb_->GetBkgProbDensity(event.GetCoord());
  }

  double GetBkgProbDensity(const I3Event& event) {
    return BackgroundLCProb(event);// * decProb_->GetBkgProbDensity(event.GetCoord());
  } */

  // To Do: make a EventPtrList version of this, if/when ready
};


#endif // LLHTIMEDEP_LOCALCOORDBKGPROB_H_
