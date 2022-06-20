
{

  // For now, only thing this project does is load the astro
  // icetray library and provide public headers for interface

  cout << "Loading libastro\n";
  int fail = gSystem->Load("libastro");
  // 0 means loaded, 1 means already loaded, -1 means failed to load
  if (fail < 0) { 
    LOADSUCCESS = false;
    return;
  }

  TString projectSrcDir = "astro_interface/public";

  vector<string> buildList;

  buildList.push_back("AstroHeader.h");
  LOADSUCCESS = BasicBuild(projectSrcDir, buildList);
}
