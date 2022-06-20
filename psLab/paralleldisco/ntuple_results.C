#include "Riostream.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TNtuple.h"

void ntuple_results()
{
  // Create an ntuple of Rates.C info, and sum livetime
  gROOT->Reset();

  ifstream in;
  //in.open("results_GP_4059.txt");
  in.open("uhecrResults.txt");

  //TFile *outFile = new TFile("results_GP_4059.root","recreate");
  TFile *outFile = new TFile("uhecr.root","recreate");

  Float_t meanNs=0, testNs=0, ts=0, fitNs=0, fitGamma=0;
  Int_t nlines = 0;
  TNtuple *ntuple = new TNtuple("ntuple","ntuple data","meanNs:testNs:ts:fitNs:fitGamma");

  while (1) {
    in >> meanNs >> testNs >> ts >> fitNs >> fitGamma;
    if (!in.good()) 
      break;
    if (nlines < 5) 
      printf("meanNs=%5f, testNs=%5f, ts=%5f, fitNs=%5f, fitGamma=%5f\n",meanNs, testNs, ts, fitNs, fitGamma);
    ntuple->Fill(meanNs, testNs, ts, fitNs, fitGamma);
    nlines++;
  }
  printf(" found %d points\n",nlines);
  in.close();

  ntuple->Write();
  outFile->Close();

}
