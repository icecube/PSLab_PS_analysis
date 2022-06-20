
{

  // THIS PROJECT

  // for now, all code is in the public... we can move compile-only to private
  // later...
  TString projectSrcDir = "fluxusDev_PlotFlux/public";

  // LIST OF LIBRARIES

  vector<string> buildList;
  buildList.push_back("PlotFlux.C");
  buildList.push_back("read_flux_curves.C");


  LOADSUCCESS = BasicBuild(projectSrcDir, buildList);
}
