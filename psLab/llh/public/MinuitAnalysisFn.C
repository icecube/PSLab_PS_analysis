#include "llh/public/MinuitAnalysisFn.h"


void MinuitAnalysisFn::InitializeMinuitAnalysisFn(void (*fcn)(Int_t &, Double_t *, Double_t &f, Double_t *, Int_t))
{ 
  if (minuit_) { return; }  // already initialized.

  minuit_ = new TMinuit();

  minuit_->SetFCN(fcn);

  // default is Migrad
  //  minuit_->CreateMinimizer();

  // Set error Definition (1 for Chi square; 0.5 for negative log likelihood) 
  minuit_->SetErrorDef(0.5);

  // Call Migrad with 500 iterations maximum 
  minuit_->SetMaxIterations(500);
  // The MIGRAD algorithm is in general the best minimizer for nearly all 
  // functions. It is a variable-metric method with inexact line search, a 
  // stable metric updating scheme, and checks for positive definiteness. 
  // Its main weakness is that it depends heavily on knowledge of the first 
  // derivatives, and fails miserably if they are very inaccurate. 

  // Set Print Level (-1 no output; 1 standard output)
  minuit_->SetPrintLevel(-1);

  /* THESE ARE PART OF OLDER SETTINGS... UPDATE AND APPLY FOR MINUIT2 ??

  // Set No Warnings 
  minuit_->mnexcm("SET NOW", arglist, 1, ierflg);

  // Minimization strategy (1 standard; 2 try to improve minimum (slower)) 
  arglist[0]=2; 
  minuit_->mnexcm("SET STR", arglist, 1, ierflg); 
  */ 
} 


double MinuitAnalysisFn::GetPar(int i) const
{
  double par=0., err=0.;
  minuit_->GetParameter(i, par, err); 
  return par;
}

// This is clumsy but seems like the only way to get this info
double MinuitAnalysisFn::Get_logLambdaBest() const 
{
  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, istat;
  minuit_->mnstat(amin, edm, errdef, nvpar, nparx, istat);
  return -amin;  // that is, max llh = - (minimizer result)
}
