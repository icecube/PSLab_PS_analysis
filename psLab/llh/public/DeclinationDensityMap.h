#ifndef LLH_DECLINATIONDENSITYMAP_H_
#define LLH_DECLINATIONDENSITYMAP_H_

#include <vector>

#include "llh/public/classes.h"
#include "llh/public/CoordClasses.h"
#include "llh/public/I3Event.h"



double ConvertBinToValue(int iBin, 
			 int nBins, double rangeMin, double rangeMax);

int ConvertValueToBin(double x, int nBins, double rangeMin, double rangeMax);


class DeclinationDensityMap {
 private:
  double decSigmaRange_;  
  // this is how far (in sigmas) to keep calculating in declination

  int nBinsDec_;
  vector<double> avDensityAtDec_; 
  double integratedDensitySum_;

  double BinToDecDeg(int iBin) const {
    return ConvertBinToValue(iBin, nBinsDec_, -90., 90.);
  }

  int DecDegToBin(double decDeg) const {
    return ConvertValueToBin(decDeg, nBinsDec_, -90., 90.);
  }

  void Reset() {
    nBinsDec_ = 0;
    avDensityAtDec_.clear();
    integratedDensitySum_ = 0.;
  }

 public:

  DeclinationDensityMap() { 
    nBinsDec_ = 0;
    integratedDensitySum_ = 0.;

    decSigmaRange_ = 5.0001; // default number of sigmas to keep calcluating at
    // ... skip decs which are farther away to save time...
    // 5 sigma range of 2d Gaussian includes 99.9996% of the distribution.
    // the extra 0.0001 is to avoid floating point discrepancies
    // regarding inclusion of the furthest bin when bin size = sigma size
    // (e.g. 1 degree binning and 1 degree sigma)
  }

  virtual ~DeclinationDensityMap() { Reset(); }

  void SetNBinsDec(const int nBinsDec) {
    Reset();
    nBinsDec_ = nBinsDec;
    avDensityAtDec_.assign(nBinsDec_,0.);
  }

  void SetDecSigmaRange(double decSigmaRange) {decSigmaRange_=decSigmaRange; }

  void SetMap(int nBinsDec, const Coord& inputCoord, double sigmaSmooth);


  int GetNBinsDec() const { return nBinsDec_; }

  double GetDecSigmaRange() const {return decSigmaRange_; }

  double GetIntegratedDensitySum() const { return integratedDensitySum_; }

  double GetProbDensityAtDecDeg(double decDeg) const;

  double GetProbDensity(const Coord& coord) const;
    
  DeclinationDensityMap operator+ (const DeclinationDensityMap&);
};





// Make a catalog of dec maps corresponding to event positions
// distributed in evenly spaced increments from -90 to +90 deg declination
//
// This can then be used to quickly generate a summed declination map
// representing the background distribution based on a vector of events


class DecMapCatalog {
 private:
  int nMaps_;
  double sigmaSmooth_;
  vector<DeclinationDensityMap> decMapVect_;

 public:
  // e.g. nMaps=180 means dec. maps for sources at positions 
  // from -89.5 to 89.5 in 1deg intervals.
  // N.B. the DeclinationMaps themselves will have the same bin size as 
  // this interval

  DecMapCatalog() {
    nMaps_ = 0;
    sigmaSmooth_ = 0.;
  }

  ~DecMapCatalog() { }

  void MakeCatalog(int nMaps, double sigmaSmooth);

  const DeclinationDensityMap& GetDecMap(const Coord& inputCoord) const;

  double GetNMaps() const { return nMaps_; }
};




DeclinationDensityMap DecMapFromCatalogAndEventVector
(const DecMapCatalog& catalog, const vector<I3Event>& eVect); 

DeclinationDensityMap DecMapFromCatalogAndEventPtrList
(const DecMapCatalog& catalog, const EventPtrList& evList); 




#endif // LLH_DECLINATIONDENSITYMAP_H_
