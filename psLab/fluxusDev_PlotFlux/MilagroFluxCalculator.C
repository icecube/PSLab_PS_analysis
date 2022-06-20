

// Photon Flux

// double fluxNorm_g; // TeV^-1 cm^-2 s^-1, at Energy eNorm_g;
// double energyNorm_g;    // TeV
// double alpha_g;    // (positive) spectral index
// double energyCut_g;     // TeV

//    fluxNorm_g = k_g * (energyNorm_g)^{-alpha_g} *
//                 exp( -sqrt( energyNorm_g / energyCut_g ) )






FormulaFlux MilagroFluxCalculator_nu(double fluxNorm_g, 
				     double energyNorm_g, 
				     double alpha_g, 
				     double energyCut_g)
{
  // To calculate k, you might try to use the normalized value in the 
  // full expression:
  //
  //    k_g = fluxNorm_g /
  //    ( pow(energyNorm_g,-alpha_g) * exp(-sqrt(energyNorm_g/energyCut_g)) );
  //
  // But to agree with the flux in Halzen, Kappes, OMurchadha arXiv:0803.0314
  // you have to leave off the exponential term (which would have increased
  // k by about 25% for 300 TeV cutoff;

  double k_g = fluxNorm_g / ( pow(energyNorm_g,-alpha_g) );

  // Now Calculate neutrino parameters from Gamma parameters
  // (see Kelner, Aharonian, and Bugayov, astro-ph/0606058 for exact relation,
  //  and Kappes, Hinton, Stegmann, and Aharonian astro-ph/0607286 for 
  //  the approximation below)

  double k_nu = (0.694 - 0.16*alpha_g) * k_g;
  double alpha_nu = alpha_g;
  double energyCut_nu = 0.59 * energyCut_g;


  // Now translate from TeV to GeV for code

  double k_nu_GeV = k_nu * pow(10, -3 + 3*alpha_nu);
  double energyCut_nu_GeV = energyCut_nu * 1000;


  TString fluxString = "";
  fluxString += k_nu_GeV;
  fluxString += "*pow(x,-";
  fluxString += alpha_nu;
  fluxString += ")*exp(-sqrt(x/";
  fluxString += energyCut_nu_GeV;
  fluxString += "))";

  FormulaFlux flux(fluxString);
  return flux;
}




FormulaFlux MilagroFluxCalculator_gamma(double fluxNorm_g, 
					double energyNorm_g, 
					double alpha_g, 
					double energyCut_g)
{
  // To calculate k, you might try to use the normalized value in the 
  // full expression:
  //
  //    k_g = fluxNorm_g /
  //    ( pow(energyNorm_g,-alpha_g) * exp(-sqrt(energyNorm_g/energyCut_g)) );
  //
  // But to agree with the flux in Halzen, Kappes, OMurchadha arXiv:0803.0314
  // you have to leave off the exponential term (which would have increased
  // k by about 25% for 300 TeV cutoff;

  double k_g = fluxNorm_g / ( pow(energyNorm_g,-alpha_g) );

  // Now translate from TeV to GeV for code

  double k_g_GeV = k_g * pow(10, -3 + 3*alpha_g);
  double energyCut_g_GeV = energyCut_g * 1000;

  TString fluxString = "";
  fluxString += k_g_GeV;
  fluxString += "*pow(x,-";
  fluxString += alpha_g;
  fluxString += ")*exp(-sqrt(x/";
  fluxString += energyCut_g_GeV;
  fluxString += "))";

  FormulaFlux flux(fluxString);
  return flux;
}
