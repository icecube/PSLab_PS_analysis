{
  /* This version of the plotting script can plots significance map, events, 
   / gal plane (shifting gal center to middle), and lat and long lines.  
   / Simply change the bools at the start to configure the options you want. 
   / Calls helper script: OptSkyPlot.C
  */


  /// CONFIGURE ///
  bool OPT_PROBMAP = false;

  bool OPT_EVENTS = true;
  bool OPT_BIG_EVENTS = true; // OPT_EVENTS must be true for this to matter

  bool OPT_GRID     = true;
  bool OPT_GALPLANE = true;
  bool OPT_COORD    = true;
  bool OPT_GAL_SHIFT_COORD= false;
  ///

  gROOT->ProcessLine(".x $LAB_MAIN_DIR/skyplots/loadlibs.C");

  // IC22
  //TFile *f = new TFile("/net/user/cfinley/root/sandbox_scripts/ic22/results/IC22final_results_025_025_hAllSky.root");
  //TH2D *hInput = hAllSky;

  // IC22 events anywhere?

  //  IC40 6-month results
  //TFile *f = new TFile("/net/user/cfinley/root/sandbox_scripts/ic40/results/ic40_skymaps_superfine.root");
  //TH2D *hInput = hAllSkyFine;

  // IC40 Full sample
  //TFile *f = new TFile("/net/user/jdumm/testlab/macro_llh/ic40_full/AllSkyBasic_FineSkyMap.root");

  // IC40 Fix sample
  TFile *f = new TFile("/net/user/jdumm/testlab/macro_llh/ic40_fix/AllSkyBasic_FineSkyMap.root");
  TH2D *hInput = hAllSkyFine;

  //TFile *f = new TFile("/net/user/jdumm/testlab/macro_llh/ic40_fix/AllSkyFC_Em2results.root");
  //TH2D *hInput = hAllSkyUL;

  TGraph *gEvents = new TGraph("/net/user/jdumm/testlab/macro_llh/ic40_fix/IC40_Fix_final_RaDec_NoPoles_ASCII.txt");

  Projection *pro = new HammerAitoffProjection();

  //Projection *pro = new FlatProjection();

  //Projection *pro = new EquatorialToGalacticProjection();
  //((EquatorialToGalacticProjection*)pro)->SetShiftGalLonDeg(0.);

/*
  Projection *proGal = new EquatorialToGalacticProjection();
  Projection *proHA = new HammerAitoffProjection();
  //Projection *proHA = new FlatProjection();
  ChainProjection *chainPro = new ChainProjection();
  chainPro->SetFirstProjection(proGal);
  proHA->SetShiftLonDeg(0.);
  chainPro->SetSecondProjection(proHA);
  Projection *pro = chainPro;
*/

  gROOT->ProcessLine(".x $LAB_MAIN_DIR/skyplots/OptSkyPlot.C"); // with configurable options
}
