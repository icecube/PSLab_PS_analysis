
{
  gROOT->Macro( "$LAB_PRELOAD" );

  gSystem->AddIncludePath(" -I$I3_SRC/astro/public ");

  string this_project = "astro_interface";

  vector<string> project_dependencies;
  // none

  // For now, only thing this project does is load the astro icetray library 
  // and provide a common header which refers to the icetray project headers

  BasicLoad(this_project, project_dependencies);
}
