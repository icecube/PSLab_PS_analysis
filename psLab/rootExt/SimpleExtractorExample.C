// SimpleExtractorExample.C
//
// A way to write the variables you want from big TTrees, TNtuples, and TChains
// to a small ntuple in a single root file.
//
// SimpleExtractor uses the same formula parser as TTree->Draw(),
// so the main limitation is every variable has to be (or convert to) 
// a number (float).
// SimpleExtractor cannot extract more complex objects! 


// 
// To run the example:
//
// .x SimpleExtractorLoad.C
// .x SimpleExtractorExample.C
//
//
// (Optional: if you are running a mac, you may get better results if
// you first run:  .x loadmac.C )
//


//
// This example script:
//
//   - Loads and compiles SimpleExtractor
//   - Loads a sample chain of Root Files
//   - Defines alias expressions, to calculate the space angle differences
//   - Writes selected variables (energy, space angle) to a new root file
//
//
// Note: SimpleExtractor works the same for input from 
// TTrees, TChains, and TNtuples (which all inherit from TTree)
//
//
// (Note that SimpleExtractor also has handling for empty trees,
// see the documentation in the code for more info.)
//


{

  // Load input files //

  // there's a * in the one's place, so only the first ten files get loaded
  char *filespec = "/data/sim/IceCube/2007/filtered/level2/neutrino-generator_000651/00000-00999/level2_nugen_numu.000651.00000*_Part00000000_Part00000000.i3.gz.root";

  TChain *chain = new TChain("tree");
  int nLoaded = chain->Add(filespec);
  cout << "\n" << nLoaded << " files loaded.\n\n";



  // We can set aliases for TTree's (and TChains and TNtuples too)
  // to make things clearer

  // MC Truth //
  chain->SetAlias("mcZenRad","I3MCTree_MaxPrimary.I3Particle.dir_.zenith_");
  chain->SetAlias("mcAziRad","I3MCTree_MaxPrimary.I3Particle.dir_.azimuth_");

  // Paraboloid Fit to Reconstructed Track //
  chain->SetAlias("paraZenRad","ParaboloidFit.I3Particle.dir_.zenith_");
  chain->SetAlias("paraAziRad","ParaboloidFit.I3Particle.dir_.azimuth_");

  // Make an alias for the Space Angle Difference //
  chain->SetAlias("paraSpaceAngleDeg",
  "TMath::RadToDeg()*SpaceAngle_rad(mcZenRad,mcAziRad,paraZenRad,paraAziRad)");



  // Make SimpleExtractor object, and specify the things you want to write //


  SimpleExtractor s;

  // We'll extract 3 variables (energy, space angle, zenith)

  // 2 parameter example: (New Name , Old Expression)
  s.AddVar("mcEnergy","I3MCTree_MaxPrimary.I3Particle.energy_");

  // if New Name = Old Expression, you just need 1 parameter to specify:
  s.AddVar("paraSpaceAngleDeg");
  s.AddVar("mcZenRad");
  // These names are just the aliases we defined above



  // Specify a cut on which events get written out, if you like


  s.SetCut("mcZenRad*TMath::RadToDeg() > 90");  // keep only upgoing events


  // Write the new file!


  char *outputFileName  = "MyFile.root";
  char* outputTreeName = "ntuple";
  char *outputTreeTitle = "MyTree";

  s.MakeTree(chain, outputFileName, outputTreeName, outputTreeTitle);

}
