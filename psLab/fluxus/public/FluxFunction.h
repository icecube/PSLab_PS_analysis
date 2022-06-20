#ifndef FLUXUS_FLUXFUNCTION_H_
#define FLUXUS_FLUXFUNCTION_H_


#include "TF1.h"


// For now, units are:

// *  energy will be GeV 
// *  flux will be GeV^-1 cm^-2 s^-1
// *  spectralIndex must include sign, i.e. E^gamma, where gamma = -2



class FluxBase {
 public:
  virtual ~FluxBase() { }
  virtual double GetFlux(double energy) const = 0;
};


// THE WRAPPER FN (makes class accessible with a simple fn call)


/* EXAMPLE:
  PowerLawFlux pf(1e-8, -2);  // flux constant and spectral index
  FluxFunction(0, &pf);
  FluxFunction(1e3);  // get flux calculated at 1 TeV
*/

double FluxFunction(double energy, const FluxBase *myFluxPtr = NULL);




// POWER LAW FLUX


class PowerLawFlux : public FluxBase {
 public:
  PowerLawFlux();
  PowerLawFlux(double fluxConstant, double spectralIndex);
  ~PowerLawFlux() { }

  double GetFluxConstant() const;
  double GetSpectralIndex() const;
  
  double GetFlux(double energy) const;

 private:
  double fluxConstant_;
  double spectralIndex_;
};




//TRUNCATED POWER LAW FLUX (power law only between eMin eMax)

class TruncatedPowerLawFlux : public FluxBase {
 public:
  TruncatedPowerLawFlux();
  TruncatedPowerLawFlux(double fluxConstant, double spectralIndex, double eMin, double eMax);
  ~TruncatedPowerLawFlux() { }
  
  double GetFluxConstant() const;
  double GetSpectralIndex() const;
  double GetEnergyMin() const;
  double GetEnergyMax() const;
  double GetFlux(double energy) const;

 private:
  double fluxConstant_;
  double spectralIndex_;
  double eMin_;
  double eMax_;
};






// FORMULA FLUX


class FormulaFlux : public FluxBase {
 public:
  FormulaFlux();
  FormulaFlux(const char* formula);
  FormulaFlux(const TF1& f);

  // copy constructor (default doesn't work with TF1 data member in class!!!)
  FormulaFlux(const FormulaFlux& rhs);

  // assignment operator
  FormulaFlux& operator=(const FormulaFlux& rhs);

  // destructor
  ~FormulaFlux() { }

  const char* GetFormula() const;
  TF1 GetTF1() const;

  double GetFlux(double energy) const;

 private:
  TF1 f_;
};


#endif // FLUXUS_FLUXFUNCTION_H_
