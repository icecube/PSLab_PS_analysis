
void Hist2DNormalizeEachX(TH2 *h) {
  int nBinsX = h->GetNbinsX();
  int nBinsY = h->GetNbinsY();

  for (int ix=0; ix<=nBinsX+1; ++ix) {
    double sum = 0.;

    // iy=0 bin is underflow
    // iy=nBinsY+1 is overflow
    // Here, we normalize along the y direction including ALL values
    // in underflow and overflow (since 2D hist may have been made too small)

    for (int iy=0; iy<=nBinsY+1; ++iy) {
      sum += h->GetBinContent(ix,iy);
    }

    if (sum != 0.) {
      for (int iy=0; iy<=nBinsY+1; ++iy) {
	double newValue = h->GetBinContent(ix,iy) / sum;
	h->SetBinContent(ix,iy, newValue);
      }
    }

  }
}
