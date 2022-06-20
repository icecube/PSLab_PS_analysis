
// Define globals here if you need them in the macros
// (N.B. anything defined as part of the function below has limited scope)

TString resultsDir;
//double srcDecDeg;
char *srcName;

int nSigMin, nSigMax, nTrials;

long ranSeedGlobal;

bool isUHECR = false;

void rootSubmitMacro_Poisson_UHECR(long ranSeedIn , char *srcNameIn, int nSigMinIn, int nSigMaxIn, int nTrialsIn) {

  nSigMin = nSigMinIn;
  nSigMax = nSigMaxIn;
  nTrials = nTrialsIn;
  ranSeedGlobal = ranSeedIn;

  TString startDir = gSystem->pwd();

  resultsDir = startDir + "/results/";
  srcName = srcNameIn;

  //gROOT->ProcessLine(".x /net/user/jdumm/testlab/llh/loadlibs.C");
  //gROOT->ProcessLine(".x /net/user/jdumm/llh_IC40_v2/llh6_Ext/loadlibs.C");
  char *maindirpath = gSystem->ExpandPathName("$LAB_MAIN_DIR");
  gSystem->cd(maindirpath);
  gROOT->ProcessLine(".x llh/loadlibs.C");

  initialize_ran1(ranSeedGlobal); // Make sure this is not re-done in Track script

  //bool OPT_DISABLE_NUMU = false;
  //bool OPT_DISABLE_NUTAU = true;
  // Can have conditions here as well:
/*
  if (srcName = "GP_FermiGalDiffuse_WarrensSample") { // Warren's BDT sample
    gROOT->ProcessLine(".x /net/user/jdumm/llh_IC40_v2/llh6_Ext/ic40_WarrensSample_llh6_Ext/macro_loadClean_Fix_final_Ext_GP_30Ebins_batch.C");
  } else { // IC40 PS sample
    gROOT->ProcessLine(".x $LAB_MAIN_DIR/macro_llh/ic59/galplane/GP_FermiGalDiffuse/Track_multi_4059.C");
  }
*/
  char *macrodirpath = gSystem->ExpandPathName("$LAB_MAIN_DIR/macro_llh/ic59");
  gSystem->cd(macrodirpath);
  if (srcName == "GP_FermiGalDiffuse") {
    gROOT->ProcessLine(".x $LAB_MAIN_DIR/macro_llh/ic59/galplane/GP_FermiGalDiffuse/Track_multi_4059.C");
  } else if (srcName == "fermibubble") {
    gROOT->ProcessLine(".x $LAB_MAIN_DIR/macro_llh/ic59/fermibubble/Track_multi_4059.C");
  } else if (srcName == "UHECR") {
    gROOT->ProcessLine(".x $LAB_MAIN_DIR/macro_llh/ic59/stacking/UHECR/Track_multi_4059_3dof.C");
    isUHECR = true;
  }

  gSystem->cd(startDir);

  gROOT->ProcessLine(".x macro_Poisson_writer_UHECR.C");

}
