
{
  // note that for macs, you may prefer to run preloadmac.C first,
  // which seems to fix a couple of strange bugs

  #ifndef LABCORE_PRELOAD_C_
  #define LABCORE_PRELOAD_C_

  // This file
  cout << "Preload: " << gSystem->ExpandPathName("$LAB_PRELOAD") << endl;

  #include <vector>

  // For loading/building project code, called by loadlibs.C routines
  gROOT->ProcessLine(".L $LAB_CORE_DIR/BasicBuild.C");

  // For loading whole projects
  gROOT->ProcessLine(".L $LAB_CORE_DIR/BasicLoad.C");

  // Define base, reference directory for includes
  // headers are then defined e.g. rootExt/public/CountMonitor.h
  gSystem->AddIncludePath(" -I$LAB_MAIN_DIR ");

  // Rudimentary style enhancements
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  // Starting value for loadlib scripts, determines if loading will continue
  bool LOADSUCCESS = true; 

  // List of loaded projects
  vector<string> PROJECTS_VECTOR;

  #endif // LABCORE_PRELOAD_C_

}
