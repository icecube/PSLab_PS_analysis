
{
  //  TFile *f = new TFile("/net/user/cfinley/root/sandbox_scripts/ic40/results/ic40_skymaps_superfine.root");
  //  TH2D *hInput = hAllSkyFine;

  //  TFile *f = new TFile("/home/cfinley/lab/macro_llh/ic40_batch_AllSky_test/hists/trial_00000201_skymaps.root");
  //  TH2D *hInput = hAllSkyFine;

  TFile *f = new TFile("~/results/ic40/results/ic40_allsky.root");
  TH2D *hInput = hAllSkyFine;

  /*
  TFile *f = new TFile("~/results/ic22ps/main/IC22final_results_025_025_hAllSky.root");
  TH2D *hInput = hAllSky;
  */

  //  Projection *pro = new HammerAitoffProjection();

  //  Projection *pro = new FlatProjection();

  Projection *proGal = new EquatorialToGalacticProjection();
  Projection *proHA = new HammerAitoffProjection();
  ChainProjection *chainPro = new ChainProjection();
  chainPro->SetFirstProjection(proGal);
  proHA->SetShiftLonDeg(0.);
  chainPro->SetSecondProjection(proHA);
  Projection *pro = chainPro;


  gROOT->ProcessLine(".x SimpleSkyPlot.C");

  // Cyngus box region
  /*
  ProjectGraphLon(72., proHA, -3, 4.)->Draw(gOpt);
  ProjectGraphLon(83., proHA, -3, 4.)->Draw(gOpt);
  ProjectGraphLat(-3, proHA, 72., 83.)->Draw(gOpt);
  ProjectGraphLat( 4, proHA, 72., 83.)->Draw(gOpt);
  */

  // Another Cyngus box region
  /*
  ProjectGraphLon(73., proHA, -7, 7.)->Draw(gOpt);
  ProjectGraphLon(87., proHA, -7, 7.)->Draw(gOpt);
  ProjectGraphLat(-7, proHA, 73., 87.)->Draw(gOpt);
  ProjectGraphLat( 7, proHA, 73., 87.)->Draw(gOpt);
  */
}
