#include <vector>

#include "TSystem.h"
#include "TH2D.h"
#include "TGraph.h"

#include "Projection.h"




TH2D* ProjectHistogram(TH2D* hInput, Projection* pro, 
		       int nBinsX = 1600, int nBinsY = 800)
{
  TH2D *hProjection = new TH2D();
  hProjection->SetBins( nBinsX, pro->GetXMin(), pro->GetXMax(),
			nBinsY, pro->GetYMin(), pro->GetYMax() );

  // iterate through each bin of projection
  for (int iy = 1; iy<=nBinsY; ++iy) {
    double y = hProjection->GetYaxis()->GetBinCenter(iy);

    for (int ix = 1; ix<=nBinsX; ++ix) {
      double x = hProjection->GetXaxis()->GetBinCenter(ix);

      // for x,y in projection, figure out what lon,lat would be:
      pro->ProjectInverse(x,y);

      // only proceed if inverse projection was to a valid lon, lat
      if ( pro->IsValid() ) { 
	int nbin = hInput->FindBin( pro->Lon() , pro->Lat() );
	hProjection->SetBinContent(ix, iy, hInput->GetBinContent(nbin));
      }
    }
  }
  return hProjection;
}

TGraph* ProjectTGraph(TGraph* tgInput, Projection* pro)
{
  TGraph *g = new TGraph(*tgInput);

  // iterate through each point on graph and project it
  for (int i=0; i<tgInput->GetN(); ++i) {
    double lon = tgInput->GetX()[i];
    double lat = tgInput->GetY()[i];
    g->SetPoint(i, pro->X(lon,lat), pro->Y(lon,lat));
  }
  return g;
}


TGraph* ProjectGalacticPlane(Projection* pro,
			     double decMinDeg = -90., double decMaxDeg = 90.)
{
  vector<double> vx, vy;
  const char *filename = 
    gSystem->ExpandPathName("$LAB_MAIN_DIR/skyplots/resources/GalPlane3600.coords");
  
  TGraph *galEquatorial = new TGraph(filename);
  if (!galEquatorial) {
    printf("Galactic coord file %s was not found.\n",filename);
    return NULL;
  }

  for (int i=0; i<galEquatorial->GetN(); ++i) {
    double raDeg, decDeg;
    galEquatorial->GetPoint(i,raDeg,decDeg);

    if (decDeg < decMinDeg) { continue; }
    if (decDeg > decMaxDeg) { continue; }

    vx.push_back(pro->X(raDeg,decDeg));
    vy.push_back(pro->Y(raDeg,decDeg));
  }
  delete galEquatorial;
  TGraph *g = new TGraph(vx.size(), &(vx[0]), &(vy[0]) );
  return g;
}


// You can set the TGraph gModel in advance with color, marker, etc.
// and number of points

TGraph* ProjectGraphLon(double lon, Projection* pro, 
			double latMin, double latMax, 
			const TGraph* gModel = NULL)
{
  const int nDefault = 1600;
  TGraph *g;
  if (gModel) { g = new TGraph(*gModel); }
  else { g = new TGraph(nDefault); }
  for (int i=0; i<g->GetN(); ++i) {
    double lat = latMin + i*(latMax-latMin)/(g->GetN()-1);
    g->SetPoint(i, pro->X(lon,lat), pro->Y(lon,lat));
  }
  return g;
}

TGraph* ProjectGraphLat(double lat, Projection* pro, 
			double lonMin, double lonMax, 
			const TGraph* gModel = NULL)
{
  const int nDefault = 1600;
  TGraph *g;
  if (gModel) { g = new TGraph(*gModel); }
  else { g = new TGraph(nDefault); }
  for (int i=0; i<g->GetN(); ++i) {
    double lon = lonMin + i*(lonMax-lonMin)/(g->GetN()-1);
    g->SetPoint(i, pro->X(lon,lat), pro->Y(lon,lat));
  }
  return g;
}
