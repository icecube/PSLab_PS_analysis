
void PrintHistogramAsArray(TH1 *h) {
  const int nBins = h->GetNbinsX();

  cout << "\n";
  cout << "// TO MAKE A HISTOGRAM:\n";
  cout << "//  NOTE: boundaries = nBins+1\n";
  cout << "//        content = nBins+2 (including under- and overflow bins)\n";

  cout << "int nBins = " << nBins << ";\n";


  cout << "double xBinBoundaries[" << nBins+1 << "]={";
  // recall: bin=0 is underflow, bin=nBins+1 is overflow
  for (int i=1; i<=nBins+1; ++i) {
    cout << h->GetBinLowEdge(i);
    if (i<nBins+1)
      cout << ",";
    else
      cout << "};\n";
  }

  cout << "double yBinContent[" << nBins+2 << "]={";
  // recall: bin=0 is underflow, bin=nBins+1 is overflow
  for (int i=0; i<=nBins+1; ++i) {
    cout << h->GetBinContent(i);
    if (i<nBins+1)
      cout << ",";
    else
      cout << "};\n";
  }

  cout << "TH1D *myHist = new TH1D;\n";
  cout << "myHist->SetBins(nBins, xBinBoundaries);\n";
  cout << "myHist->SetContent(yBinContent);\n\n";


  //
  //
  //


  cout << "// TO MAKE A GRAPH:\n";
  cout << "// NOTE:  x = bin center , y = bin content\n";

  cout << "int nBins = " << nBins << ";\n";

  cout << "double x[" << nBins << "]={";
  // recall: bin=0 is underflow, bin=nBins+1 is overflow
  for (int i=1; i<=nBins; ++i) {
    cout << h->GetBinCenter(i);
    if (i<nBins)
      cout << ",";
    else
      cout << "};\n";
  }

  cout << "double y[" << nBins << "]={";
  // recall: bin=0 is underflow, bin=nBins+1 is overflow
  for (int i=1; i<=nBins; ++i) {
    cout << h->GetBinContent(i);
    if (i<nBins)
      cout << ",";
    else
      cout << "};\n";
  }

  cout << "TGraph *myGraph = new TGraph(nBins,x,y);\n\n";
}
