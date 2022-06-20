
{

  // THIS PROJECT

  TString projectSrcDir = "skyplots/src";

  // LIST OF LIBRARIES

  vector<string> buildList;

  buildList.push_back("BasicPlotFns.h");
  buildList.push_back("SetRootPalette.C");
  buildList.push_back("TransformationFns.h");
  buildList.push_back("Projection.h");
  buildList.push_back("ProjectionCases.C");
  buildList.push_back("SkyPlotFns.C");

  LOADSUCCESS = BasicBuild(projectSrcDir, buildList);
}
