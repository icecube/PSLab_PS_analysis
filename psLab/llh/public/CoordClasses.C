#include "llh/public/CoordClasses.h"

#include "TMath.h"


double CircularGaussUnc(double r, double sigma) {
  if(r>90) return 0;
  if(sigma<10) return exp(-r*r/(sigma*sigma*2)) / (2.*TMath::Pi()*sigma*sigma);
  else{
    double kDeg = 1/(sigma*sigma);
    double kRad = 1/(sigma*TMath::DegToRad() * sigma*TMath::DegToRad());
    return kDeg/(4*TMath::Pi()*sinh(kRad))*exp(kRad*cos(r*TMath::DegToRad()));
  }
  //return exp(-r*r/(sigma*sigma*2)) / (2.*TMath::Pi()*sigma*sigma);
}
