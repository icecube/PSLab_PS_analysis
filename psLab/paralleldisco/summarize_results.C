{
        gROOT->SetStyle("Plain");

//TFile *inFile = new TFile("results.root","read");
//TFile *inFile = new TFile("results_GP_4059_Try1.root","read");
//TFile *inFile = new TFile("results_GP_4059.root","read");
TFile *inFile = new TFile("results_fermibubble_fudge1.root","read");
//TFile *inFile = new TFile("results_fermibubble_naoko.root","read");

vector<int> meanNss;

// Worthwhile to automate?
meanNss.push_back(0);
//meanNss.push_back(10);
meanNss.push_back(20);
//meanNss.push_back(30);
meanNss.push_back(40);
//meanNss.push_back(50);
meanNss.push_back(60);
//meanNss.push_back(70);
meanNss.push_back(80);
//meanNss.push_back(90);
meanNss.push_back(100);
meanNss.push_back(120);
meanNss.push_back(140);
meanNss.push_back(160);
meanNss.push_back(180);
meanNss.push_back(200);
/*
meanNss.push_back(300);
meanNss.push_back(350);
meanNss.push_back(400);
meanNss.push_back(450);
meanNss.push_back(500);
meanNss.push_back(550);
*/

const int size = meanNss.size();

TH1F *hTs[size];
char nameBuf[10], cutBuf[100];
int sumOfHistEntries=0;
int nBins=100;
double xMin=0, xMax=30;

for (int i=0; i<size; i++) {
  sprintf(nameBuf,"h%d",meanNss[i]);
  sprintf(cutBuf,"meanNs==%d",meanNss[i]);
  hTs[i] = new TH1F(nameBuf, nameBuf, nBins, xMin, xMax);
  ntuple->Draw("ts>>"+TString(nameBuf),cutBuf,"goff");
  cout << nameBuf << " " << hTs[i]->GetEntries() << endl;
  sumOfHistEntries += hTs[i]->GetEntries();

}


if (ntuple->GetEntries() != sumOfHistEntries)
  { cout << "You left out some entries!!!  Okay, if that's what you want...\n"; }

gStyle->SetOptStat("");

TCanvas *can = new TCanvas("can","can",900,700);
can->Divide(2,2);
can->cd(1)->SetLogy(1);
//can->SetLogy(1);

hTs[0]->SetTitle("Test Statistic Summary");
hTs[0]->SetXTitle("-log_{10} #lambda");
hTs[0]->SetLineWidth(2);
hTs[0]->Draw();
int entries0 = hTs[0]->GetEntries();
int entries = 0;
int color=2;
for (int i=1; i<size; i++) {
  if (color==17) color+=3;
  hTs[i]->SetLineColor(color++);
  entries = hTs[i]->GetEntries();
  hTs[i]->Scale(entries0/entries);
  hTs[i]->Draw("same");
}

TLegend *leg = new TLegend(0.70,0.95-size*0.05,1.0,1.0);
//leg->SetFillColor(kWhite);
char legBuf[100];
//sprintf(legBuf,"ns=%d, mean=%0.2f",meanNss[0],hTs[0]->GetMean());
sprintf(legBuf,"mean ns=%d",meanNss[0]);
leg->AddEntry(hTs[0],legBuf,"l");
for (int i=1; i<size; i++) {
  //sprintf(legBuf,"ns=%d, mean=%0.2f",meanNss[i],hTs[i]->GetMean());
  sprintf(legBuf,"mean ns=%d",meanNss[i]);
  leg->AddEntry(hTs[i],TString(legBuf),"l");
}
leg->Draw();

// 1 d.o.f.
TF1 *fchi2_1 = new TF1("fchi2_1","[0] * [1] * (1./sqrt(2*TMath::Pi())) * exp(-[1]*x/2)/sqrt([1]*x)",0,100);
// 2 d.o.f.
TF1 *fchi2_2 = new TF1("fchi2_2","[0] * [1] * (1./2) * exp(-[1]*x/2)", 0,100);
// 3 d.o.f.
TF1 *fchi2_3 = new TF1("fchi2_3","[0] * [1] * (1./sqrt(2*TMath::Pi())) * exp(-[1]*x/2)*sqrt([1]*x)",0,100);


TF1 *fchi2 = new TF1("fchi2","[0] * [1] * (1./sqrt(2*TMath::Pi())) * exp(-[1]*x/2)/sqrt([1]*x)",0,50); // 2 d.o.f.

double normFactor = h0->GetEntries();
normFactor = normFactor/2.;
// 2 correction because we consider only excess (+nSrc) contributions
//   to chi2 distribution,  not deficits (-nSrc)
double xRange = h0.GetXaxis()->GetXmax() - h0.GetXaxis()->GetXmin();
double binWidth = xRange / h0.GetNbinsX();
normFactor = normFactor * binWidth; // to match histogram plot

fchi2->SetParameter(0,normFactor);
fchi2->SetParameter(1,2);  // correction for chi = 2 * log lambda

//h0.Draw();
fchi2->Draw("same");

// SUMMARIZE FRACTION OF TRIALS DETECTED AT TWO SIGNIFICANCES
double maxMean=0;
for (int i=0; i<size; i++) {
  if (maxMean < meanNss[i]) maxMean=meanNss[i];
}

can->cd(3);
//TCanvas *can3 = new TCanvas("can3","can3",700,500);
//TH1F *hFrac = new TH1F("hFrac","hFrac",100,0,maxMean*1.2);

double frac[size];
double meanNsArray[size];
char bufCut1[100], bufCut2[100];
for (int i=0; i<size; i++) {
  sprintf(bufCut1,"ts>0 && meanNs==%d",meanNss[i]);
  sprintf(bufCut2,"meanNs==%d",meanNss[i]);
  frac[i] = (double)ntuple->GetEntries(bufCut1)/(double)ntuple->GetEntries(bufCut2);
  meanNsArray[i]=meanNss[i];
  //hFrac->SetBinContent(bin, frac); 
}
//hFrac->Draw();
TGraph *gFrac = new TGraph(size, meanNsArray, frac);

gFrac->SetMarkerStyle(20);
gFrac->GetYaxis()->SetRangeUser(0.,1.);
gFrac->SetTitle("Median Sensitivity");
gFrac->GetXaxis()->SetTitle("Mean N Signal Events");
gFrac->GetYaxis()->SetTitle("Frac Above Bkg-only Median Signif");
gFrac->Draw("apl");
double sens=0;
double maxNs=200.;
for (int i=0; i<1000; i++) {
  sens = 0.+(200.-0.)/1000*i;
  if (gFrac->Eval(sens) >= 0.9) break;
}
if (sens >= maxNs-1.) cout << "You reached the maximum range of this sensitivity algorithm!\n";

TLine *line = new TLine(0,0.9,sens,0.9);
TLine *line2 = new TLine(sens,0.9,sens,0.0);
line->SetLineColor(kRed);
line->SetLineWidth(2);
line->Draw();
line2->SetLineColor(kRed);
line2->SetLineWidth(2);
line2->Draw();

can->cd(2)->SetLogy(1);

//gROOT->ProcessLine(".L ../rootExt/public/FunctionsRoot.C");
  TString startDir = gSystem->pwd();

  char *maindirpath = gSystem->ExpandPathName("$LAB_MAIN_DIR");
  gSystem->cd(maindirpath);
  gROOT->ProcessLine(".x llh/loadlibs.C");
  gSystem->cd(startDir);

TH1F *h0Int = (TH1F*)DescendingCumulate(h0);
h0Int->Scale(1./h0Int->GetBinContent(1));
h0Int->SetTitle("Descending Cumulative Prob for Bkg only");
h0Int->SetXTitle("-log_{10} #lambda");
h0Int->GetYaxis()->SetRangeUser(5.e-8,2.);
h0Int->Draw();

TF1 *fexp = new TF1("fexp","[0]*exp([1]*x)",0.4,9.);
h0Int->Fit(fexp,"r","",0.4,5.);

double thresh=0;
for (int i=0; i<1000; i++) {
  thresh = 9.-(9.)/1000*i;
  if (fexp->Eval(thresh) >= 2.87e-7) break;
  //if (fexp->Eval(thresh) >= 1.35e-3) break;
}
cout << "5-sigma (one-sided) significance threshold = " << thresh << endl;
//cout << "3-sigma (one-sided) significance threshold = " << thresh << endl;
TLine *line3 = new TLine(0,2.87e-7,thresh,2.87e-7);
//TLine *line3 = new TLine(0,1.35e-3,thresh,1.35e-3);
TLine *line4 = new TLine(thresh,2.87e-7,thresh,5.e-8);
//TLine *line4 = new TLine(thresh,1.35e-3,thresh,5.e-8);
line3->SetLineColor(kRed);
line3->SetLineWidth(2);
line3->Draw();
line4->SetLineColor(kRed);
line4->SetLineWidth(2);
line4->Draw();
// From fitting 133900 background-only trials, llh threshold for 5-sigma (1-sided)
//  is equal to 11.416

can->cd(4)->SetLogy(0);

char bufCut1[100], bufCut2[100];
for (int i=0; i<size; i++) {
  sprintf(bufCut1,"ts>%f&& meanNs==%d",thresh,meanNss[i]);
  sprintf(bufCut2,"meanNs==%d",meanNss[i]);
  frac[i] = (double)ntuple->GetEntries(bufCut1)/(double)ntuple->GetEntries(bufCut2);
  meanNsArray[i]=meanNss[i];
  //hFrac->SetBinContent(bin, frac); 
}
//hFrac->Draw();
TGraph *gFrac2 = new TGraph(size, meanNsArray, frac);

gFrac2->SetMarkerStyle(20);
gFrac2->GetYaxis()->SetRangeUser(0.,1.);
gFrac2->SetTitle("Discovery Potential");
gFrac2->GetXaxis()->SetTitle("Mean N Signal Events");
gFrac2->GetYaxis()->SetTitle("Frac Above Bkg-only 5#sigma Signif");
//gFrac2->GetYaxis()->SetTitle("Frac Above Bkg-only 3#sigma Signif");
gFrac2->Draw("apl");

double discpot=0;
for (int i=0; i<3000; i++) {
  discpot = 0.+(300.-0.)/3000*i;
  if (gFrac2->Eval(discpot) >= 0.5) break;
}

TLine *line5 = new TLine(0,0.5,discpot,0.5);
TLine *line6 = new TLine(discpot,0.5,discpot,0.0);
line5->SetLineColor(kRed);
line5->SetLineWidth(2);
line5->Draw();
line6->SetLineColor(kRed);
line6->SetLineWidth(2);
line6->Draw();

cout << "Sensitivity = " << sens << " events" << endl;
cout << "Discovery Potential = " << discpot << " events" << endl;

}
