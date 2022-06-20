#include "rootExt/public/root_draw_tools.h"

#include "TPad.h"


TLine* vline(Double_t x, Color_t lcolor, Width_t lwidth, Style_t lstyle)
{
  TLine *templine;
  if (gPad->GetLogy()) 
    templine = new TLine(x,pow(10,gPad->GetUymin()),
                         x,pow(10,gPad->GetUymax()) );
  else
    templine = new TLine(x, gPad->GetUymin(), x, gPad->GetUymax() );
  templine->SetLineColor(lcolor);
  templine->SetLineWidth(lwidth);
  templine->SetLineStyle(lstyle);
  return templine;
}


TLine* hline(Double_t y, Color_t lcolor, Width_t lwidth, Style_t lstyle)
{
  TLine *templine;
  if (gPad->GetLogx()) 
    templine = new TLine(pow(10,gPad->GetUxmin()), y,
                         pow(10,gPad->GetUxmax()), y );
  else
    templine = new TLine(gPad->GetUxmin(),y,gPad->GetUxmax(),y );
  templine->SetLineColor(lcolor);
  templine->SetLineWidth(lwidth);
  templine->SetLineStyle(lstyle);
  return templine;
}

