
{

  // THIS PROJECT

  TString projectSrcDir = "plotscripts/public";

  // LIST OF LIBRARIES

  vector<string> buildList;
  buildList.push_back("EffectiveArea.C");

  LOADSUCCESS = BasicBuild(projectSrcDir, buildList);
}
