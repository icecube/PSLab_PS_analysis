#include "fluxus/public/FluxFunction.h"


// THE WRAPPER FN  (makes class accessible with a simple fn call)

// NOTE: when you set the FluxFunction, no energy will be calculated 
// and zero is returned.  This frees the user from having to worry about
// providing a valid energy value which may not be known yet

// default *myFlux = NULL is specified in header
double FluxFunction(double energy, const FluxBase *myFluxPtr) {
  static const FluxBase *myFlux_ = NULL;

  // check if user is specifying a new flux model
  if (myFluxPtr) {
    myFlux_ = myFluxPtr; 
    return 0.;
  }
  // otherwise, must be requesting a flux
  else {
    return myFlux_->GetFlux(energy);
  }
}


// POWER LAW FLUX


PowerLawFlux::PowerLawFlux() : fluxConstant_(0), spectralIndex_(0) { }

PowerLawFlux::PowerLawFlux(double fluxConstant, double spectralIndex) :
  fluxConstant_(fluxConstant),
  spectralIndex_(spectralIndex) { }

double PowerLawFlux::GetFluxConstant() const { return fluxConstant_; }

double PowerLawFlux::GetSpectralIndex() const { return spectralIndex_; }

double PowerLawFlux::GetFlux(double energy) const {
  return fluxConstant_*pow(energy,spectralIndex_);
}

// TRUNCATED POWER LAW FLUX

TruncatedPowerLawFlux::TruncatedPowerLawFlux() : fluxConstant_(0), spectralIndex_(0), eMin_(0), eMax_(0) { }

TruncatedPowerLawFlux::TruncatedPowerLawFlux(double fluxConstant, double spectralIndex, double eMin, double eMax) :
  fluxConstant_(fluxConstant),
  spectralIndex_(spectralIndex),
  eMin_(eMin),
  eMax_(eMax) { }

double TruncatedPowerLawFlux::GetFluxConstant() const { return fluxConstant_; }

double TruncatedPowerLawFlux::GetSpectralIndex() const { return spectralIndex_; }

double TruncatedPowerLawFlux::GetEnergyMin() const {return eMin_;}
double TruncatedPowerLawFlux::GetEnergyMax() const {return eMax_;}

double TruncatedPowerLawFlux::GetFlux(double energy) const {

  if(energy > eMin_ && energy < eMax_)
    return fluxConstant_*pow(energy,spectralIndex_);
  else
    return 0;
}



// FORMULA FLUX


FormulaFlux::FormulaFlux() : f_("","0") { }
FormulaFlux::FormulaFlux(const char* formula) : f_("",formula) { }
FormulaFlux::FormulaFlux(const TF1& f) : f_(f) { }

// copy constructor:
FormulaFlux::FormulaFlux(const FormulaFlux& rhs) : 
  FluxBase(rhs),
  f_(rhs.f_) 
{ }

// assignment operator:
FormulaFlux& FormulaFlux::operator=(const FormulaFlux& rhs) {
  f_ = rhs.f_;
  return *this;
}

const char* FormulaFlux::GetFormula() const { return f_.GetTitle(); }

TF1 FormulaFlux::GetTF1() const { return f_; }

double FormulaFlux::GetFlux(double energy) const { return f_.Eval(energy); }
