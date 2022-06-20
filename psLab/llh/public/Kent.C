#ifndef __KENT_INCLUDED
#define __KENT_INCLUDED

#include "TF2.h"
#include "TH2D.h"
#include "TMath.h"

// 
// A macro to return a source pdf based on the Kent (FB5) distribution (Beta=0).
// This pdf represents a gaussian with mean at srcMeanZenDeg, srcMeanAziDeg,
// having a width of srcSigmaDeg (acts just like a flat 2d Gaussian sigma)
//
TH2D* Kent(double srcMeanZenDeg, double srcMeanAziDeg, double srcSigmaDeg){
//TF2* Kent(double srcMeanZenDeg, double srcMeanAziDeg, double srcSigmaDeg){
  
  double degToRad = TMath::DegToRad();

  if (srcSigmaDeg==0) srcSigmaDeg+=0.0001; // Add a little fudge if 0
  if (fabs(srcMeanAziDeg)<0.01 && srcSigmaDeg<0.1) srcMeanAziDeg=0.01; // Add a little fudge if 0, has trouble at azi=0 'boundary' when sigma is smaller than function grid spacing

  int nbins = 100;
  if (fabs(srcMeanAziDeg)<0.01 && srcSigmaDeg<3.0 && 3*srcSigmaDeg>fabs(srcMeanAziDeg))
  {
    cout << "WARNING: Increasing Kent resolution since near azimuth boundary (at 0).\n";
    cout << "       -- May decrease speed of signal generation!\n";
    nbins = 500;
  }
  
  double srcMeanZenRad = srcMeanZenDeg*degToRad; 
  double srcMeanAziRad = srcMeanAziDeg*degToRad;
  double srcSigmaRad   = srcSigmaDeg*degToRad;
  double xmin, ymin, xmax, ymax;

  xmin = srcMeanAziRad - 3*srcSigmaRad/sin(srcMeanZenRad);
  xmax = srcMeanAziRad + 3*srcSigmaRad/sin(srcMeanZenRad);
  if (xmax>2*TMath::Pi() || xmin<=0) 
  {
    xmin = 0; 
    xmax = 2*TMath::Pi();
  }
  ymin = srcMeanZenRad - 3*srcSigmaRad;
  ymax = srcMeanZenRad + 3*srcSigmaRad;
  if (ymax > TMath::Pi()) ymax = TMath::Pi();
  if (ymin < 0) ymin = 0;
  // x = Phi, y = Theta
  TF2 *fSourcePdf;

  // Approx to avoid exceeding float limits!
  if ( srcSigmaDeg < 8 )
  {
    fSourcePdf = new TF2("Kent",
    "[0]*[1]**2*exp( -1/(2*[1]**2)*(acos(sin(y)*sin([2])*cos(x-[3])+cos(y)*cos([2])))**2 )"
    ,xmin,xmax,ymin,ymax); // Kent, sigma ~< 8deg
  }
  // The real deal (Beta=0)
  else
  {
    fSourcePdf = new TF2("Kent",
    "[0]*exp( 1/([1]**2)*(sin(y)*sin([2])*cos(x-[3])+cos(y)*cos([2]) ) )"
    ,xmin,xmax,ymin,ymax); // Kent, sigma >~ 8deg
  }

  fSourcePdf->SetParameters(2.*TMath::Pi(), srcSigmaRad, srcMeanZenRad, srcMeanAziRad);
  fSourcePdf->SetNpx(nbins); // Default binning for histo is 100
  fSourcePdf->SetNpy(nbins);

  // Warning, DO NOT USE IN A LOOP, will probably give seg faults!
  // Usually need to delete anything created with "new"

  //fSourcePdf->Draw("surf2"); // Seems like it must be drawn or GetHistogram fails!!
  TH2D *hSourcePdf = (TH2D*)fSourcePdf->GetHistogram();
  hSourcePdf->Scale(1./hSourcePdf->Integral(0,100,0,100)); // Norm
  //delete fSourcePdf;  // <- causes errors!!
 
 //hSourcePdf->Multiply("1.", 1./hSourcePdf->Integral(0,100,0,100));
  //TH2D *returnPDF = hSourcePdf->;
  //return returnPDF;

  return hSourcePdf; // normalized pdf in zen and azi
  //return fSourcePdf; // causes seg faults: similar to old bug in ROOT?
                       // http://root.cern.ch/root/roottalk/roottalk99/0095.html

}

#endif
