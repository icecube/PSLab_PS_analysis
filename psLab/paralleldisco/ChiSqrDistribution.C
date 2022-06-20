{

TF1 *fchi2 = new TF1("fchi2","[0] * [1] * (1./sqrt(2*TMath::Pi())) * exp(-[1]*x/2)/sqrt([1]*x)",0,50);

//double normFactor = nTrials;
//double normFactor = 1e3;
double normFactor = 3.43e4;
normFactor = normFactor/2.;
// 2 correction because we consider only excess (+nSrc) contributions
//   to chi2 distribution,  not deficits (-nSrc)
double xRange = h0.GetXaxis()->GetXmax() - h0.GetXaxis()->GetXmin();
//double binsPerUnit = h0.GetNbinsX() / xRange;
//normFactor = normFactor / binsPerUnit; // to match histogram plot
double binWidth = xRange / h0.GetNbinsX();
normFactor = normFactor * binWidth; // to match histogram plot

fchi2->SetParameter(0,normFactor);
fchi2->SetParameter(1,2);  // correction for chi = 2 * log lambda

h0.Draw();
fchi2->Draw("same");

}
