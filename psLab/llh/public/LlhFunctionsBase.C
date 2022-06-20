#include "llh/public/LlhFunctionsBase.h"

#include "rootExt/public/log_report.h"


LogLikelihood::LogLikelihood() :
  nParDefs_(0),
  parDefArray_(NULL),
  parArray_(NULL),
  errorArray_(NULL),
  minuit_(NULL)
{ }


LogLikelihood::LogLikelihood(int npar) :
  nParDefs_(npar),
  parDefArray_(new MinuitParDef[nParDefs_]),
  parArray_(new double[nParDefs_]),
  errorArray_(new double[nParDefs_])
{
  // by default, all parameters are initially set to use automatic definitions
  SetParAuto(); 

  // parDefArray_[].name's  should be set in each derived class constructor

  for (int i=0; i<nParDefs_; ++i) {
    parArray_[i] = 0.;
    errorArray_[i] = 0.;
  }

  // MINUIT

  // evidently this has to be here, rather than in the initialization list above
  minuit_ = new TMinuit(nParDefs_);


  // ** DEFAULT MINUIT SETTINGS.  ** //
  //    Can change by using LogLikelihood::GetMinuit()   //

  // for more info, see http://root.cern.ch/root/html/TMinuit.html

  double arglist[10];  // use for passing arguments to control minuit
  int ierflg = 0;  // use for checking that minuit executed commands

  // Set Print Level (-1 no output; 1 standard output)
  minuit_->SetPrintLevel(-1);

  // Set error Definition (1 for Chi square; 0.5 for negative log likelihood) 
  minuit_->SetErrorDef(0.5); 

  // Set No Warnings 
  minuit_->mnexcm("SET NOW", arglist, 1, ierflg);

  // Minimization strategy (1 standard; 2 try to improve minimum (slower)) 
  arglist[0]=2; 
  minuit_->mnexcm("SET STR", arglist, 1, ierflg); 

  // Call Migrad with 500 iterations maximum 
  minuit_->SetMaxIterations(500); 
  // The MIGRAD algorithm is in general the best minimizer for nearly all 
  // functions. It is a variable-metric method with inexact line search, a 
  // stable metric updating scheme, and checks for positive definiteness. 
  // Its main weakness is that it depends heavily on knowledge of the first 
  // derivatives, and fails miserably if they are very inaccurate. 

}


LogLikelihood::~LogLikelihood() {
  if (minuit_) { delete minuit_; }
  if (parDefArray_) { delete [] parDefArray_; }
  if (parArray_) { delete [] parArray_; }
  if (errorArray_) { delete [] errorArray_; }
}



void LogLikelihood::AssignMinuitParDef() {
  assert(minuit_);
  minuit_->mncler();  // reset all paramters to undefined
  for (int i=0; i<nParDefs_; ++i) {
    // minuit_->DefineParameter(i, "name", initVal, initStepSz, lowLim, upLim);
    minuit_->DefineParameter(i, 
			     parDefArray_[i].name.c_str(),
			     parDefArray_[i].initValue, 
			     parDefArray_[i].initStepSize, 
			     parDefArray_[i].lowLimit, 
			     parDefArray_[i].upLimit);
    if (parDefArray_[i].optFix) {
      minuit_->FixParameter(i);
    }
  }
}


// Defaults in header: optFix = false, optAuto = false
// (Typical case for when you, the user, put in your own parameter settings)
void LogLikelihood::SetParDef(unsigned int i, 
			      double initValue, double initStepSize,
			      double lowLimit, double upLimit, 
			      bool optFix, double optAuto) 
{
  if (int(i) < nParDefs_ ) {
    parDefArray_[i].initValue = initValue;
    parDefArray_[i].initStepSize = initStepSize;
    parDefArray_[i].lowLimit = lowLimit;
    parDefArray_[i].upLimit = upLimit;
    parDefArray_[i].optFix = optFix;
    parDefArray_[i].optAuto = optAuto;
  } else {
    log_error("Minuit Parameter %d not defined.\n",i);
  }
}




double LogLikelihood::GetPar(int i) const {
  if (i>=nParDefs_) {
    log_error("Requested Fit parameter index %d, but only %d defined\n", 
	      int(i), nParDefs_);
    return -1;
  }
  assert(minuit_);
  double currentValue, currentError;
  minuit_->GetParameter(int(i), currentValue, currentError);
  return currentValue;
}
