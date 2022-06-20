#include "TCanvas.h"
#include "TH1D.h"
#include "TLine.h"


void ColorChooser()
{

  TCanvas *canCol = new TCanvas("canColorChooser","canColorChooser",1000,500);
  canCol->cd();

  int nCols = 129;

  TH1D *hCol = new TH1D("hCol","Color Chooser",nCols,0,nCols+1);
  hCol->SetMinimum(0.);
  hCol->SetMaximum(1.);

  hCol->Draw();

  TLine *line;

  for (int i=1; i<=nCols; ++i) {
    line = new TLine(i,0.05,i,1);
    line->SetLineColor(i);
    line->SetLineWidth(2);
    line->Draw();
  }

}
