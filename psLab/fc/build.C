
{

  // THIS PROJECT

  TString projectSrcDir = "fc/public";

  // LIST OF LIBRARIES

  vector<string> buildList;

  buildList.push_back("GeneralFC.C");
  buildList.push_back("GeneralFC_Functions.C");

  LOADSUCCESS = BasicBuild(projectSrcDir, buildList);
}
