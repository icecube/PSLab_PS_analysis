
{
  gROOT->Macro( "$LAB_PRELOAD" );

  this_project = "llh";

  /*TString LAB_COORDINATE_PROJECT = getenv("LAB_COORDINATE_PROJECT");
  if (LAB_COORDINATE_PROJECT != "astro_interface" &&
      LAB_COORDINATE_PROJECT != "coord_interface") {
    cout << "\nLAB_COORDINATE_PROJECT not specified....\n";
    cout << "Will use coordinate-service as default.\n";
    cout << "(See 'labcore/start_lab.sh' or ";
    cout << "'llh/loadlibs.C' for how to change)\n\n";

    LAB_COORDINATE_PROJECT = "coord_interface";
  }*/

  project_dependencies.clear();
  project_dependencies.push_back("rootExt");
  project_dependencies.push_back("fluxus");
  //project_dependencies.push_back(string(LAB_COORDINATE_PROJECT));
  project_dependencies.push_back("coord_interface");

  BasicLoad(this_project, project_dependencies);
}
