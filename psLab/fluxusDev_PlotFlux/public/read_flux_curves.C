#include "fluxusDev_PlotFlux/public/read_flux_curves.h"

#include <vector>

#include "TGraph.h"

#include "rootExt/public/TStringify.h"

#include "fluxus/public/FluxFunction.h"

#include "fluxusDev_PlotFlux/public/PlotFlux.h"


void read_flux_curves(char* filename, FluxPlotManager& fpm,
		      FluxPlotManager& fpmEvents,
		      int color, int width, int style,
		      bool optConnect, int connectPoint) {

  FILE *fp = fopen(filename,"r");
  if (!fp) { cout << "Error: file not found\n"; return; }

  double eMin, eMax, index, flux, nSig;
  char buffer[10000] = {0};

  vector<TGraph> gVect;

  while ( fscanf(fp,"%lg %lg %lg %lg %lg\n",
		 &eMin,&eMax,&index,&flux,&nSig) == 5 ) 
  {
    fscanf(fp,"%[^\n]",buffer);
    cout << buffer << endl;

    TGraph *g;

    if (optConnect) {
      if (connectPoint == -1) {
	g = GraphFlux(buffer, 0, eMin, eMin*1.001, 1);// just one point on left
      }
      else if (connectPoint == 1) {
	g = GraphFlux(buffer, 0, eMax*0.999, eMax, 1);//just one point at right
      } 
      else {
	g = GraphFlux(buffer, 0, eMin, eMax, 1);  // just one point in center
      }
      
      gVect.push_back(*g);
    } else {
      g = GraphFlux(buffer, 0, eMin, eMax);
      g->SetLineColor(color);
      g->SetLineWidth(width);
      g->SetLineStyle(style);
      fpm.AddGraph(*g);
    }

    g = GraphFlux(TStringify(nSig), 0, eMin, eMax);
    g->SetLineColor(color);
    g->SetLineWidth(width);
    g->SetLineStyle(style);
    fpmEvents.AddGraph(*g);

    delete g;
  }

  if (optConnect) {
    double x;
    double y;

    int nStart = 0;
    while (gVect[nStart].GetPoint(0,x,y) , y == 0.) {
      nStart++;
    }

    int nEnd = gVect.size()-1;
    while (gVect[nEnd].GetPoint(0,x,y) , y == 0.) {
      nEnd--;
    }

    cout << nStart << " , " << nEnd << endl;

    TGraph *g = new TGraph(nEnd+1-nStart);
    g->SetLineColor(color);
    g->SetLineWidth(width);
    g->SetLineStyle(style);

    int iOut = 0;
    for (int i=nStart; i<=nEnd; ++i) {
      gVect[i].GetPoint(0, x, y);
      g->SetPoint(iOut,x,y);
      iOut++;
    }

    fpm.AddGraph(*g);
    delete g;
  }
}
