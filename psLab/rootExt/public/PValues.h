#ifndef ROOTEXT_PVALUES_HEADER_
#define ROOTEXT_PVALUES_HEADER_

#include "TMath.h"


double TwoSidedSigmaToPValue(double sigma) {
  return 1.-TMath::Erf(sigma/sqrt(2.));
}

double OneSidedSigmaToPValue(double sigma) {
  return (1.-TMath::Erf(sigma/sqrt(2.))) / 2.;
}

double PValueToTwoSidedSigma(double pvalue) {
  return TMath::ErfInverse(1.-pvalue)*sqrt(2.);
}

double PValueToOneSidedSigma(double pvalue) {
  return TMath::ErfInverse(1.-2.*pvalue)*sqrt(2.);
}


#endif // ROOTEXT_PVALUES_HEADER_
