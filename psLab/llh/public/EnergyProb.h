#ifndef LLH_ENERGYPROB_H_
#define LLH_ENERGYPROB_H_


// #include "llh/public/classes.h"
// Forward Declarations (when feasible, more efficient than including headers)
class Event;


class EnergyProb {
 public:
  virtual ~EnergyProb() { }

  virtual double GetEnergyProbGamma(const Event& event, double gamma) const = 0;
  virtual double GetEnergyProbBkg(const Event& event) const = 0;
  virtual double GetEnergyMaxRatio(const Event& event) const = 0;

  virtual double GetGammaMin() const = 0;
  virtual double GetGammaMax() const = 0;
};

#endif // LLH_ENERGYPROB_H_
