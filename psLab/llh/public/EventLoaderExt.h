#ifndef __EVENTLOADEREXT_HEADER_INCLUDED
#define __EVENTLOADEREXT_HEADER_INCLUDED

#include <vector>
#include "TString.h"
#include "TH2D.h"
#include "TStopwatch.h"

#include "llh/public/EventLoader.h"
#include "llh/public/CoordEquatorialDeg.h"
#include "llh/public/I3Event.h"


/* A class for loading source candidate events for spatially extended sources.
 / Implemented source shapes include a Kent distribution (2D Normal on a sphere) 
 / and Beta-profiles for clusters of galaxies.  These "profiles" are de-projected
 / to build a 3D model of the cluster and multiplied (bin-by-bin in 3D) with some
 / CR distribution (model-dependent) to get a neutrino emission profile in 2D.
 / Methods given are only a good approximation for "small-scale" extended sources, 
 / i.e. less than ~10deg in size.
*/

class EventLoaderExt : public EventLoader {

 private:
  EquatorialDeg srcEqDeg_;

  // The spatial PDF
  TH2D *hSourcePdf_;

  // Each event gets rotated to new location and re-sampled for Extended sources
  int nUpSamples_;
  

 public:
  
  EventLoaderExt() { 
    monitor_  = true;

    name_timeMJD_ = "UNDEFINED";
    name_recoZenith_rad_ = "UNDEFINED";
    name_recoAzimuth_rad_ = "UNDEFINED";
    name_sigmaDeg_ = "UNDEFINED";
    name_energyValue_ = "UNDEFINED";
    name_runID_ = "UNDEFINED";
    name_eventID_ = "UNDEFINED";

    bkgTree_ = NULL;
    bkgWeight_ = "UNDEFINED";
    minZenithRad_ = 0.;
    maxZenithRad_ = 0.;
    spreadZenithRad_ = 0.;
    bkgLoadMethod_ = EXACT;
    bkgTimeMethod_ = SCRAMBLE;

    sourceTree_ = NULL;
    zenWidthDeg_ = 0.;
    name_mcEnergy_GeV_ = "UNDEFINED";
    name_mcZenith_rad_ = "UNDEFINED";
    name_mcAzimuth_rad_ = "UNDEFINED";
    name_mcOneWeight_ = "UNDEFINED";
    totalGeneratedEvents_ = 0.;

    nUpSamples_ = 20;
    hSourcePdf_ = NULL;
  }

  // stopwatch
  TStopwatch stopwatch_LoadSourceEvents_;


  // Sets the number of up-samplings for each signal event
  void SetNUpSamples(double nUp) {nUpSamples_ = (int)nUp;}
  // Gets the number of up-samplings for each signal event
  double GetNUpSamples() {return nUpSamples_;}

  // Kent acts like a 2D Gaussian but on a sphere.
  // Requires a central coordinate and a size sigma.
  void LoadSourcePdf_Kent(EquatorialDeg srcEqDeg, double sigmaSrc);

  void SetEqDegRa(double ra){ srcEqDeg_.SetRaDeg(ra); }
  void SetEqDegDec(double dec){ srcEqDeg_.SetDecDeg(dec); }
  void SetEqDeg(EquatorialDeg srcCoord) {
    SetEqDegRa(srcCoord.GetRa());
    SetEqDegDec(srcCoord.GetDec());
  }
  
  // 4 different possibilities to model for Clusters:
  // Central AGN (approx as a point source instead!)
  //void LoadSourcePdf_CG_CentralAGN(EquatorialDeg srcEqDeg, double z, double rCore, bool singleBeta=1, double n1=0, double rc1=0, double beta1=0, double n2=0, double rc2=0, double beta2=0);

  /* Model "A", Model "B", and "I"sobaric Model can be loaded with this function. 
   / These "profiles" are de-projected from x-ray observations to build a 3D model 
   / of the cluster and multiplied (bin-by-bin in 3D) with some CR distribution 
   / (model-dependent) to get a neutrino emission profile in 2D.
   / Models comes from Murase 2008 (arXiv:0805.0104)
  */
  void LoadSourcePdf_CG_ABI(char *model, EquatorialDeg srcEqDeg, double z, double rVirialMpc, bool singleBeta=1, double n1=0, double rc1=0, double beta1=0, double n2=0, double rc2=0, double beta2=0);

  void LoadSourceEvents(vector<I3Event>& srcCandidateEvents) ;
  void LoadSourceEvents(vector<I3Event>& srcCandidateEvents,
                        EquatorialDeg srcEqDeg, bool rotate=true);

  double GetZenWidthDeg() { return zenWidthDeg_; }
};


#endif
