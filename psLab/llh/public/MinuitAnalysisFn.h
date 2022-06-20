#ifndef LLH_MINUITANALYSISFN_H_
#define LLH_MINUITANALYSISFN_H_

#include <string>

//#include "Minuit2/FCNBase.h"
#include "TMinuit.h"

#include "llh/public/classes.h"


//class MinuitAnalysisFn : public ROOT::Minuit2::FCNBase , public AnalysisFn {
class MinuitAnalysisFn : public AnalysisFn {
 protected:
  TMinuit* minuit_;

 public:
  MinuitAnalysisFn() : minuit_(NULL) { }
  virtual ~MinuitAnalysisFn() { if (minuit_) delete minuit_; }

  // This has to be called by each derived class constructor, pointing to 'this'
  virtual void InitializeMinuitAnalysisFn(void (*fcn)(Int_t &, Double_t *, Double_t &f, Double_t *, Int_t));

  // Here we allow direct access to the minimizer for user configuration
  TMinuit* Minuit() { return minuit_; }

  // ** Pure Virtual Function inherited from FCNBase ** //
  // Up() and operator() 
  // must be defined in derived class to instantiate an object

  // This is related to how "errors" are defined, see Minuit2 for documentation
  virtual double Up() const { return 0.5; }

  // This is what you have to define for your likelihood!
  virtual double operator() (const vector<double>& par) const = 0;


  // Functions inherited from AnalysisFn base class still needing definitions:

  virtual void PrepareAnalysis() = 0;
  virtual void MaximizeLlh() = 0;

  virtual double EvaluateLlh(double *parValueArray) = 0;

  virtual double GetPar(int i) const;
  virtual double Get_logLambdaBest() const;
  virtual double GetTestStatistic() const = 0;
  virtual double GetEstProb() const = 0;
};




// For MultiAnalysis: Translate Global Par values to Individual Par Values:


class ParTranslator {
 public:
  virtual ~ParTranslator() { }

  virtual const vector<double>
    Translate(int llhIndex, const vector<double>& parGlobal) const = 0;
};


#endif // LLH_MINUITANALYSISFN_H_

