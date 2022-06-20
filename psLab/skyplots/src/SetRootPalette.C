#include <iostream>
#include "TColor.h"
#include "TStyle.h"

void
SetRootPalette(const Int_t pal = 0)
{
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  if (pal == 0) {
    cout << "\n  ROOT Palette Function\n  Usage: SetRootPalette(pal)\n\n"
         << "    o pal =  0: print help\n"
         << "    o pal =  1: rainbow\n"
         << "    o pal =  2: reverse-rainbow\n"
         << "    o pal =  3: amber\n"
         << "    o pal =  4: reverse-amber\n"
         << "    o pal =  5: blue/white\n"
         << "    o pal =  6: white/blue\n"
         << "    o pal =  7: red temperature\n"
         << "    o pal =  8: reverse-red temperature\n"
         << "    o pal =  9: green/white\n"
         << "    o pal = 10: white/green\n"
         << "    o pal = 11: orange/blue\n"
         << "    o pal = 12: blue/orange\n"
         << "    o pal = 13: white/black\n"
         << "    o pal = 14: black/white\n"
         << "    o pal = 15: H.E.S.S.\n"
         << "    o pal = 16: IceCube-22\n"
         << endl;
  }
  else if (pal == 1) {
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
  }
  else if (pal == 2) {
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.51, 1.00, 0.87, 0.00, 0.00 };
    Double_t green[NRGBs] = { 0.00, 0.20, 1.00, 0.81, 0.00 };
    Double_t blue[NRGBs]  = { 0.00, 0.00, 0.12, 1.00, 0.51 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
  }
  else if (pal == 3) {
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.17, 0.39, 0.62, 0.79, 1.00 };
    Double_t green[NRGBs] = { 0.01, 0.02, 0.39, 0.68, 1.00 };
    Double_t blue[NRGBs]  = { 0.00, 0.09, 0.18, 0.09, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
  }
  else if (pal == 4) {
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 1.00, 0.79, 0.62, 0.39, 0.17 };
    Double_t green[NRGBs] = { 1.00, 0.68, 0.39, 0.02, 0.01 };
    Double_t blue[NRGBs]  = { 0.00, 0.09, 0.18, 0.09, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
  }
  else if (pal == 5) {
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.00, 0.38, 1.00 };
    Double_t green[NRGBs] = { 0.00, 0.00, 0.38, 0.76, 1.00 };
    Double_t blue[NRGBs]  = { 0.00, 0.47, 0.83, 1.00, 1.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
  }
  else if (pal == 6) {
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 1.00, 0.38, 0.00, 0.00, 0.00 };
    Double_t green[NRGBs] = { 1.00, 0.76, 0.38, 0.00, 0.00 };
    Double_t blue[NRGBs]  = { 1.00, 1.00, 0.83, 0.47, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
  }
  else if (pal == 7) {
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.50, 0.89, 0.95, 1.00 };
    Double_t green[NRGBs] = { 0.00, 0.00, 0.27, 0.71, 1.00 };
    Double_t blue[NRGBs]  = { 0.00, 0.00, 0.00, 0.40, 1.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
  }
  else if (pal == 8) {
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 1.00, 0.95, 0.89, 0.50, 0.00 };
    Double_t green[NRGBs] = { 1.00, 0.71, 0.27, 0.00, 0.00 };
    Double_t blue[NRGBs]  = { 1.00, 0.40, 0.00, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
  }
  else if (pal == 9) {
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.38, 0.75, 1.00 };
    Double_t green[NRGBs] = { 0.00, 0.35, 0.62, 0.85, 1.00 };
    Double_t blue[NRGBs]  = { 0.00, 0.00, 0.00, 0.47, 1.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
  }
  else if (pal == 10) {
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 1.00, 0.75, 0.38, 0.00, 0.00 };
    Double_t green[NRGBs] = { 1.00, 0.85, 0.62, 0.35, 0.00 };
    Double_t blue[NRGBs]  = { 1.00, 0.47, 0.00, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
  }
  else if (pal == 11) {
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.75, 1.00, 0.24, 0.00, 0.00 };
    Double_t green[NRGBs] = { 0.24, 1.00, 0.75, 0.18, 0.00 };
    Double_t blue[NRGBs]  = { 0.00, 0.62, 1.00, 0.68, 0.12 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
  }
  else if (pal == 12) {
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.24, 1.00, 0.75 };
    Double_t green[NRGBs] = { 0.00, 0.18, 0.75, 1.00, 0.24 };
    Double_t blue[NRGBs]  = { 0.12, 0.68, 1.00, 0.62, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
  }
  else if (pal == 13) {
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t green[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t blue[NRGBs]  = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
  }
  else if (pal == 14) {
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 1.00, 0.84, 0.61, 0.34, 0.00 }; 
    Double_t green[NRGBs] = { 1.00, 0.84, 0.61, 0.34, 0.00 };
    Double_t blue[NRGBs]  = { 1.00, 0.84, 0.61, 0.34, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
  }
  else if (pal == 15) {
    // H.E.S.S. Color Palette
    Double_t stops[6] = { 0.,
			  30./128,
			  60./128, 
			  90./128,
			  120./128,
			  1.00 };
    Double_t red[6]   = { 0.00, 0.00, 0.30, 1.00, 1.00, 1.00 }; 
    Double_t green[6] = { 0.00, 0.00, 0.10, 0.00, 1.00, 1.00 };
    Double_t blue[6]  = { 0.25, 0.55, 0.85, 0.00, 0.00, 1.00 };
    TColor::CreateGradientColorTable(6, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
  }
  else if (pal == 16) {
    // IceCube-22 Color Palette
    const int cgNMax = 100;       // this is to dimension the array
    int cgNum = 4;  // this is actual number of stops we will use
    double cgStops[cgNMax] = { 0., 40./90, 70./90, 1. };
    double cgRed[cgNMax]   = { 1.00, 0.00, 1.00, 1.00 }; 
    double cgGreen[cgNMax] = { 1.00, 0.40, 0.30, 1.00 };
    double cgBlue[cgNMax]  = { 1.00, 1.00, 0.30, 0.30 };
    int nCgContours = 255;
    TColor::CreateGradientColorTable(cgNum, cgStops, 
				     cgRed, cgGreen, cgBlue, nCgContours);
    gStyle->SetNumberContours(nCgContours);
  }
  else if (pal == 17) {
    // IceCube-40 6-months Color Palette (enhanced detail)
    const int cgNMax = 100;       // this is to dimension the array
    int cgNum = 5;  // this is actual number of stops we will use
    double cgStops[cgNMax] = { 0.00, 0.01, 0.45, 0.75, 1.00 };
    double cgRed[cgNMax]   = { 1.00, 0.90, 0.00, 1.00, 1.00 }; 
    double cgGreen[cgNMax] = { 1.00, 0.94, 0.40, 0.30, 1.00 };
    double cgBlue[cgNMax]  = { 1.00, 1.00, 1.00, 0.30, 0.30 };
    int nCgContours = 255;
    TColor::CreateGradientColorTable(cgNum, cgStops, 
				     cgRed, cgGreen, cgBlue, nCgContours);
    gStyle->SetNumberContours(nCgContours);
  }
  else {
    cout << "Palette " << pal << " was not understood.  Try '0' for info\n";
  }
}
