
{
  gROOT->Macro( "$LAB_PRELOAD" );

  string this_project = "llhSN";

  vector<string> project_dependencies;
  project_dependencies.push_back("rootExt");
  project_dependencies.push_back("fluxus");
  project_dependencies.push_back("coord_interface");
  project_dependencies.push_back("llh");
  project_dependencies.push_back("llhTimeDep");

  BasicLoad(this_project, project_dependencies);
}
