
{

  // THIS PROJECT

  // for now, all code is in the public... we can move compile-only to private
  // later...
  projectSrcDir = "fluxus/public";

  // LIST OF LIBRARIES

  buildList.clear();
  buildList.push_back("FluxFunction.C");


  LOADSUCCESS = BasicBuild(projectSrcDir, buildList);
}
