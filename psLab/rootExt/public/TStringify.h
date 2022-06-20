#ifndef ROOTEXT_TSTRINGIFY_H_
#define ROOTEXT_TSTRINGIFY_H_


#include "TString.h"

// Convert a double to TString with high precision if required

TString TStringify(double x) {
  char temp[100];
  // formerly, this was %lg, but that didn't keep enough digits for double
  // precision.  This should work better: up to 16 digits if necessary.
  sprintf(temp,"%.16g",x);
  return TString(temp);
}


#endif // ROOTEXT_TSTRINGIFY_H_INCLUDED
