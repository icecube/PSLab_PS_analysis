
{
  gROOT->Macro( "$LAB_PRELOAD" );

  string this_project = "skyplots";

  vector<string> project_dependencies;
  project_dependencies.push_back("coord_interface");

  BasicLoad(this_project, project_dependencies);
}
