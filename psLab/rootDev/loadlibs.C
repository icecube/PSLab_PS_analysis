
{
  gROOT->Macro( "$LAB_PRELOAD" );

  this_project = "rootDev";

  project_dependencies.clear();
  project_dependencies.push_back("rootExt");

  BasicLoad(this_project, project_dependencies);
}
