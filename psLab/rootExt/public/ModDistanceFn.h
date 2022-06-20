#ifndef ROOTEXT_MODDISTANCEFN_H_
#define ROOTEXT_MODDISTANCEFN_H_


double ModDistanceFn(double x1, double x2, double divisor) {
  double a = fmod( x1-x2 , divisor);
  // range is now (-divisor , divisor)
  // change to [-divisor/2. , divisor/2.]
  if (a >  divisor/2.) { a -= divisor; }
  else if (a < -divisor/2.) { a += divisor; }
  return fabs(a);
}

#endif // ROOTEXT_MODDISTANCEFN_H_
