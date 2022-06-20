
{
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
