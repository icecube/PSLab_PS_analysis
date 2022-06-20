
// Define globals here if you need them in the macros
// (N.B. anything defined as part of the function below has limited scope)

TString startDir;
TString resultsDir;
double srcDecDeg;
char *srcSingleName;

int nsrc;
int nSigMin, nSigMax, nTrials;
double sigma;
double decDeg;
long ranSeedI;

double period[7];
double t0[7], p0[7], bfsigma[7];
TString srcName[7];

void rootSubmitMacro_integrated(long ranSeedIn, double decDegIn, int nSigMinIn, int nSigMaxIn, int nTrialsIn) {

  //nsrc    = nsrcIn;
  //sigma   = bfsigma[nsrc]; // UL by default
  //if (sigmaIn) { sigma = sigmaIn; }
  nSigMin = nSigMinIn;
  nSigMax = nSigMaxIn;
  nTrials = nTrialsIn;
  decDeg = decDegIn;

  startDir = gSystem->pwd();

  resultsDir = startDir + "/results/";
  //srcSingleName = srcName[nsrc];

  //gROOT->ProcessLine(".x /net/user/jdumm/testlab/llh/loadlibs.C");
  //gROOT->ProcessLine(".x /net/user/jdumm/llh_IC40_v2/llh6_Ext/loadlibs.C");
  //gSystem->cd("/net/user/mfbaker/lab/IC40");
  //gROOT->ProcessLine(".x loadlibs.C");

  //bool OPT_DISABLE_NUMU =  true;
  //bool OPT_DISABLE_NUTAU = true;
  //gROOT->ProcessLine(".x /net/user/mfbaker/lab/IC40/macro_loadClean_IC40Full_newEbins.C");
  //gROOT->ProcessLine(".x macro_loadClean_Fix_WithTaus_final_batch.C");
  
  ranSeedI = ranSeedIn;
  
  //gSystem->cd(startDir);

  gROOT->ProcessLine(".x macro_fc_writer_integrated.C");

}
