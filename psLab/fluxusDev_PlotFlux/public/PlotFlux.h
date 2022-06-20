#ifndef FLUXUSDEV_PLOTFLUX_PLOTFLUX_H_
#define FLUXUSDEV_PLOTFLUX_PLOTFLUX_H_

#include <vector>

#include "TGraph.h"
#include "TH1D.h"
#include "TString.h"

#include "fluxus/public/FluxFunction.h"


class FluxPlotManager {
 public:
  double yOutputUnits;
  double yPower;
  double yMin;
  double yMax;
  double logEMin;
  double logEMax;
  TH1D hFrame;

  vector<TGraph> graphVect;
  vector<TGraph> plotGraphVect;

  FluxPlotManager() :
    yOutputUnits(1e9), // GeV
    yPower(2),
    yMin(0),
    yMax(0),
    logEMin(2),
    logEMax(9)
  { }

  ~FluxPlotManager() { }

  void Clear() { graphVect.clear(); }

  TGraph* GetGraph(int i) { return &(graphVect[i]); }

  void AddGraph(TGraph g) { graphVect.push_back(g); }

  void Plot();
};


TGraph* GraphFlux(FluxBase &flux, double yPower, 
		  double xmin=1e2, double xmax=1e9, int nPoints = 250);

TGraph* GraphFlux(TString formula, double yPower, 
		  double xmin=1e2, double xmax=1e9, int nPoints = 250) {
  FormulaFlux flux(formula);
  return GraphFlux(flux, yPower, xmin, xmax, nPoints);
}


TGraph* PlotFlux(FluxBase &flux, double yPower, 
		 double xmin=1e2, double xmax=1e9, 
		 double ymin=0, double ymax=0);

TGraph* PlotFlux(TString formula, double yPower, 
		 double xmin=1e2, double xmax=1e9, 
		 double ymin=0, double ymax=0) {
  FormulaFlux flux(formula);
  return PlotFlux(flux, yPower, xmin, xmax, ymin, ymax);
}




#endif // FLUXUSDEV_PLOTFLUX_PLOTFLUX_H_
