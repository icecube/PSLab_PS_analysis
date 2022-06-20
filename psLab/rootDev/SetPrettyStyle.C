#include "TStyle.h"

void SetPrettyStyle() {

  // Canvas printing details: white bg, no borders.
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(0);

  // Canvas frame printing details: white bg, no borders.
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(0);

  // Plot title details: centered, no bg, no border, nice font.
  gStyle->SetTitleX(0.1);
  gStyle->SetTitleW(0.8);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFillColor(0);

  // Font details for titles and labels.
  gStyle->SetTitleFont(42, "xyz");
  gStyle->SetTitleFont(42, "pad");
  gStyle->SetLabelFont(42, "xyz");
  gStyle->SetLabelFont(42, "pad");

  // Details for stat box.
  gStyle->SetStatColor(0);
  gStyle->SetStatFont(42);
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatX(0.975);
  gStyle->SetStatY(0.9);

  // gStyle->SetOptStat(0);
}
