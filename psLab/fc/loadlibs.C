
{
  gROOT->Macro( "$LAB_PRELOAD" );

  string this_project = "fc";

  vector<string> project_dependencies;
  project_dependencies.push_back("llh");
  
  BasicLoad(this_project, project_dependencies);
}
