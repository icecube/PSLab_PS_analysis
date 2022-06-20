#include "llh/public/DeclinationDensityMap.h"

#include "rootExt/public/log_report.h"

#include "llh/public/CoordEquatorialDeg.h"


// SUPPORT FUNCTIONS


double ConvertBinToValue(int iBin, 
			 int nBins, double rangeMin, double rangeMax) {
  // General error checking
  if (iBin<0 || iBin>=nBins) {
    log_error("Called FindValueFromBin with i=%d outside range (0,%d)\n",
	      iBin,nBins);
  }
  double x = rangeMin + (rangeMax-rangeMin)*(iBin+0.5)/nBins;
  return x;
}


int ConvertValueToBin(double x, int nBins, double rangeMin, double rangeMax) {
  // General error checking
  if (x<rangeMin || x>rangeMax) {
    log_error("Called FindBinFromValue with x=%lg outside range (%lg,%lg)\n",
	      x,rangeMin,rangeMax);
  }

  int iBin = int( nBins * (x-rangeMin)/(rangeMax-rangeMin) );

  // if x==rangeMax exactly, we tolerate this and 
  // treat it as still within range, and set the bin to the last actual bin
  if (iBin == nBins) {
    iBin = nBins-1;
  }

  return iBin;
}


// MAIN CLASS 


void DeclinationDensityMap::SetMap(int nBinsDec, const Coord& inputCoord, 
				   double sigmaSmooth) 
{
  SetNBinsDec(nBinsDec);

  int nBinsRa = 2*nBinsDec_;  // scale ra step size with dec size
  double raStep = 360./nBinsRa;

  const EquatorialDeg* inputEq = 
    dynamic_cast<const EquatorialDeg*>(&inputCoord);
  assert(inputEq);
  double inputDecDeg = inputEq->GetDec();

  for (int iy = 0; iy<nBinsDec_; ++iy) {
    double decDeg = BinToDecDeg(iy);
    double sumAtDec = 0.;

    if (fabs(decDeg - inputDecDeg) < decSigmaRange_*sigmaSmooth) {
      for (double raDeg = raStep/2.; raDeg<360.; raDeg += raStep) {
	EquatorialDeg eqCoord(raDeg,decDeg);
	double rDeg = eqCoord.DistanceTo(inputCoord);
	sumAtDec += CircularGaussUnc(rDeg, sigmaSmooth);
      }
      // note since we're working in deg, the prob is per sq. deg, not sr
    }

    avDensityAtDec_[iy] = sumAtDec / nBinsRa;
    // this is the average "density of events per sq. deg",
    // as a function of declination
    
    integratedDensitySum_ += avDensityAtDec_[iy] *
      (180./nBinsDec_) *  // integrate "dec deg"
      360. * cos( decDeg*TMath::DegToRad() );   // integrate "ra deg"
  }
}


double DeclinationDensityMap::GetProbDensityAtDecDeg(double decDeg) const 
{
  if (nBinsDec_) {
    int iy = DecDegToBin(decDeg);
    return avDensityAtDec_[ iy ] / integratedDensitySum_;
  } else {
    log_error("DeclinationDensityMap not set before getting prob.\n");
    assert(false);
  }
}


double DeclinationDensityMap::GetProbDensity(const Coord& coord) const 
{
  const EquatorialDeg* eq = 
    dynamic_cast<const EquatorialDeg*>(&coord);
  assert(eq);
  return GetProbDensityAtDecDeg(eq->GetDec());
}


DeclinationDensityMap DeclinationDensityMap::operator+ 
(const DeclinationDensityMap &decMap2) 
{
  if (nBinsDec_ != decMap2.nBinsDec_) {
    log_fatal("Cannot add two DeclinationDensityMaps with different array"
	      "sizes.  Stop.\n");
    assert(false);
  }

  DeclinationDensityMap decMapNew;
  decMapNew.SetNBinsDec(nBinsDec_);

  for (int i=0; i<nBinsDec_; ++i) {
    decMapNew.avDensityAtDec_[i] = 
      avDensityAtDec_[i] + decMap2.avDensityAtDec_[i];
  }
  //  bkgNew.decMinDeg_ = decMinDeg_;
  // bkgNew.decMaxDeg_ = decMaxDeg_;
  decMapNew.integratedDensitySum_ = 
    integratedDensitySum_ + decMap2.integratedDensitySum_;

  return decMapNew;
}



// CLASS DECMAPCATALOG //


void DecMapCatalog::MakeCatalog(int nMaps, double sigmaSmooth) 
{
  decMapVect_.clear();
  decMapVect_.assign(nMaps, DeclinationDensityMap() );

  for (int i=0; i<nMaps; ++i) {
    // locate dec. of simulated 'event' at center of ith bin
    double decDeg = ConvertBinToValue(i, nMaps, -90., 90.);

    // make a coordinate
    EquatorialDeg eq(0.,decDeg);

    // make the declination map for an event at this coordinate
    // (maps will have same binning as catalog)
    decMapVect_[i].SetMap(nMaps, eq, sigmaSmooth);
  }

  nMaps_ = nMaps;
  sigmaSmooth_ = sigmaSmooth;
}  



const DeclinationDensityMap& 
DecMapCatalog::GetDecMap(const Coord& inputCoord) const
{
  if (!decMapVect_.size()) {
    log_error("Called DecMapCatlog::GetDecMap before decMapVect was made.\n");
  }

  const EquatorialDeg* eq = 
    dynamic_cast<const EquatorialDeg*>(&inputCoord);
  assert(eq);

  double decDeg = eq->GetDec();

  int i = ConvertValueToBin(decDeg, nMaps_, -90., 90);
        
  return decMapVect_[i];
}



DeclinationDensityMap DecMapFromCatalogAndEventVector
(const DecMapCatalog &catalog, const vector<I3Event>& eVect) 
{
  DeclinationDensityMap outMap;
  outMap.SetNBinsDec(catalog.GetDecMap(EquatorialDeg(0,0)).GetNBinsDec());

  for (unsigned int i=0; i<eVect.size(); ++i)
  {
    outMap = outMap + catalog.GetDecMap( eVect[i].GetCoord() );
  }
  return outMap;
}


DeclinationDensityMap DecMapFromCatalogAndEventPtrList
(const DecMapCatalog &catalog, const EventPtrList& evList) 
{
  DeclinationDensityMap outMap;
  outMap.SetNBinsDec(catalog.GetDecMap(EquatorialDeg(0,0)).GetNBinsDec());

  for (int i=0; i<evList.GetSize(); ++i) {
    outMap = outMap + catalog.GetDecMap( evList.GetEvent(i)->GetCoord() );
  }
  return outMap;
}
