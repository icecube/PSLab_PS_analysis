
#ifndef __BETAPROF_INCLUDED
#define __BETAPROF_INCLUDED

#include "TF1.h"
#include "TH1D.h"

// 
// A macro to return a source pdf based on the Beta Profiles of 
//  Pfrommer and Ensslin (2004).  Ignores weak temperature dependence (Murase 2008).
//
//TH1D* BetaProf(bool singleBeta=1, double n1=0, double rc1=0, double beta1=0, double n2=0, double rc2=0, double beta2=0, double rMax=5000){
TF1* BetaProf(bool singleBeta=1, double n1=0, double rc1=0, double beta1=0, double n2=0, double rc2=0, double beta2=0, double rMaxMpc=5){

  TF1 *eDensity;
  // "x" (or "r") is in units of kpc (for Hubble Constant = 70 km/s/Mpc)
  double xMin = 0;
  double xMax = rMaxMpc*1e3; // R_virial range:[1.8, 3.8] Mpc, 5 Mpc is good default 
  if (singleBeta){ // single-beta profile (non-cool core cluster)
    eDensity = new TF1("eDensity",
      "[0]*(1+x**2/[1]**2)**(-3*[2]/2)+[3]*(1+x**2/[4]**2)**(-3*[5]/2)",
      xMin,xMax);
    eDensity->SetParameters(n1,rc1,beta1,n2,rc2,beta2);
  } else { // Must be double-beta profile (cool core clusters)
    eDensity = new TF1("eDensity",
      "(1*[0]*(1+x**2/[1]**2)**(-3*[2]/2)+[3]*(1+x**2/[4]**2)**(-3*[5]))**0.5",
      xMin,xMax); // This excludes Temp profiles, Lambda-Twiddle = ?
    eDensity->SetParameters(n1,rc1,beta1,n2,rc2,beta2);
  }

  // Default sampling size
  eDensity->SetNpx(100);

  // Must visualize or "GetHistogram()" won't work
  //eDensity->Draw();
  //return (TH1D*)eDensity->GetHistogram();
  return eDensity;

}

#endif
