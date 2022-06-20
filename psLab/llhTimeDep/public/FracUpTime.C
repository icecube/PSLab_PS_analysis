#include "llhTimeDep/public/FracUpTime.h"

#include <fstream>

#include "TMath.h"


// This takes in a histogram of the uptimes and limits you're interested in,
// and returns the fraction of uptime for things like fluence scaling.

// First off, uptime between some min to max times

double FracUpTime(TH1D* a, double xmin, double xmax) {
  if (xmax < xmin) {
    double a = xmin; xmin = xmax; xmax=a;
  }
  int start = a->FindBin(xmin);
  int finish = a->FindBin(xmax);
  
  double sc = a->GetBinContent(start);
  double fc = a->GetBinContent(finish);
  
  double sle = a->GetBinLowEdge(start+1);
  double fle = a->GetBinLowEdge(finish);
  
  double result;
  
  if (start == finish) {
    result = sc*(xmax-xmin);
  } else {
    result = sc*(sle-xmin) + fc*(xmax-fle);
    for (int i=start+1;i<finish;i++){
      result += a->GetBinContent(i)*(a->GetBinLowEdge(i+1) - a->GetBinLowEdge(i));
    }
  }
  
  result /= (xmax-xmin);
  return result;
}

// This does the same thing but is specifically for a block function.
// farther down it does the same for a Gaussian.

double BlockFracUpTime(TH1D* a, string fname, double thresh) {
  
  ifstream fin;
  fin.open(fname.c_str());
  double blockbegin, blockdur, blocklev;
  double bigresult=0;
  double totaltime=0;
  while (fin >> blockbegin) {
    fin >> blocklev >> blockdur;
    if (blocklev > thresh) {
      bigresult += FracUpTime(a, blockbegin, blockbegin+blockdur )*blockdur;
      totaltime += blockdur;
    }
  }
  fin.close();
  
  return bigresult/totaltime;
  
}

double FracUpTimeGaus(TH1D* a, double mean, double sigma) {

  // I want to go to +/- 4 sigma, or the edges of the uptime histogram
  double tmin = max( a->GetXaxis()->GetXmin(), mean-4.*sigma );
  double tmax = min( a->GetXaxis()->GetXmax(), mean+4.*sigma );
 
  int start  = a->FindBin(tmin);
  int finish = a->FindBin(tmax);
  
  double sc = a->GetBinContent(start);
  double fc = a->GetBinContent(finish);
  
  double sle = a->GetBinLowEdge(start+1);
  double fle = a->GetBinLowEdge(finish);
  
  double total = TMath::Erf( (tmax-mean)/sigma ) - TMath::Erf( (tmin-mean)/sigma );
  
  double result=0.;
  double h,l;
  
  if (start == finish) {
    result = sc*total;
  } else {
    result = sc*(TMath::Erf( (sle-mean)/sigma ) - TMath::Erf( (tmin-mean)/sigma ) ) + 
               fc*( TMath::Erf( (tmax-mean)/sigma ) - TMath::Erf( (fle-mean)/sigma ) );
    for (int i=start+1;i<finish;i++){
    
      h = a->GetBinLowEdge(i+1);
      l = a->GetBinLowEdge(i);
      result += a->GetBinContent(i)*(TMath::Erf( (h-mean)/sigma ) - TMath::Erf( (l-mean)/sigma ));
      
    }
  }
  
  result /= total;
  return result;
}
