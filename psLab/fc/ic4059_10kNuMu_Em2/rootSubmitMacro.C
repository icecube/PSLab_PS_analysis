
// Define globals here if you need them in the macros
// (N.B. anything defined as part of the function below has limited scope)

TString resultsDir;
double srcDecDeg;


void rootSubmitMacro(long ranSeed , double decDeg) {

  TString startDir = gSystem->pwd();

  resultsDir = startDir + "/results/";
  srcDecDeg = decDeg;

  //gSystem->cd("/net/user/jdumm/testlab");
  gROOT->ProcessLine(".x /net/user/jdumm/testlab/llh/loadlibs.C");

  initialize_ran1(ranSeed);

  //gSystem->cd(startDir);

  //gROOT->ProcessLine(".x macro_load_fromSubmit.C");
  //gROOT->ProcessLine(".x /net/user/jdumm/testlab/macro_llh/ic40_fix/macro_loadClean_Fix_final_batch.C");
  gROOT->ProcessLine(".x /net/user/jdumm/testlab/macro_llh/ic40_fix/macro_loadClean_Fix_WithTaus_final_batch.C");

  gROOT->ProcessLine(".x macro_fc_writer.C");

}
