#ifndef LLH_LLHFUNCTIONSBASE_H_
#define LLH_LLHFUNCTIONSBASE_H_

#include <string>

#include "TMinuit.h"

#include "llh/public/classes.h"
#include "llh/public/MinuitWrapperClass.h"



class MinuitParDef {
 public:
  // default constructor
  MinuitParDef() :
  name(""), initValue(0), initStepSize(0),
    lowLimit(0), upLimit(0), optFix(false), optAuto(true) { }

  MinuitParDef(string Name, double InitValue, double InitStepSize,
	       double LowLimit, double UpLimit, 
	       bool OptFix = false, bool OptAuto = true) :
  name(Name), initValue(InitValue), initStepSize(InitStepSize),
    lowLimit(LowLimit), upLimit(UpLimit), optFix(OptFix), optAuto(OptAuto) { }

  string name;
  double initValue;
  double initStepSize;
  double lowLimit;
  double upLimit;
  bool optFix;
  bool optAuto;
};



// THIS IS COMPLETELY GENERAL

class LogLikelihood : public MinuitWrapperClass , public AnalysisFn {
// class LogLikelihood : public MinuitWrapperClass {
 protected:

  int nParDefs_;
  MinuitParDef *parDefArray_;
  double *parArray_;
  double *errorArray_;
  TMinuit *minuit_;

  // This passes the parameter definitions into minuit itself...
  // presumably you want to do this just before minimizing.
  void AssignMinuitParDef();

  /*
int i, string name, 
			  double initValue=0., double initStepSize=0., 
			  double lowLimit=0., double upLimit=0.,
			  bool optFix=false, bool optAuto=true) {
    parDefArray_[i].name = name;
    parDefArray_[i].optAuto = optAuto;
    SetMinuitParDef(name, initValue, initStepSize, lowLimit, upLimit, 
		    optFix, optAuto);
  }
  */

 public:

  LogLikelihood();
  LogLikelihood(int npar);
  virtual ~LogLikelihood();

  virtual TMinuit* GetMinuit() { return minuit_; }

  // Defaults: optFix = false, optAuto = false
  // (Typical case for when you, the user, put in your own parameter settings)
  virtual void SetParDef(unsigned int i, double initValue, double initStepSize,
			 double lowLimit, double upLimit, 
			 bool optFix = false, double optAuto = false);


  void SetParAuto(int i, bool optAuto) { parDefArray_[i].optAuto = optAuto; }
  void SetParAuto() { 
    for (int i=0; i<nParDefs_; ++i) { SetParAuto(i,true); }
  }

  virtual void PrepareAnalysis() = 0;
  virtual void MaximizeLlh() = 0;

  virtual double EvaluateLlh(double *parValueArray) = 0;

  virtual double GetPar(int i) const;
  virtual double Get_logLambdaBest() const = 0;
  virtual double GetTestStatistic() const = 0;
  virtual double GetEstProb() const = 0;
};



// THIS IS FOR ANY ANALYSIS SET, STILL PRETTY GENERAL (Cartesian, I3Event...)

class AnalysisLlh : public LogLikelihood {
 protected:

 public:
  AnalysisLlh() { }
  AnalysisLlh(int npar) : LogLikelihood(npar) { }
  virtual ~AnalysisLlh() { }

  virtual void PrepareAnalysis() = 0;
  virtual void MaximizeLlh() = 0;

  virtual double EvaluateLlh(double *parValueArray) = 0;

  virtual double Get_logLambdaBest() const = 0;
  virtual double GetTestStatistic() const = 0;
  virtual double GetEstProb() const = 0;
};


#endif // LLH_LLHFUNCTIONSBASE_H_
