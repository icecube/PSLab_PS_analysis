
{
  TString dir = gSystem->ExpandPathName("$LAB_MAIN_DIR/fluxusDev_PlotFlux/");
  gROOT->ProcessLine(".L "+dir+"MilagroFluxCalculator.C");

  // Photon Flux
  double fluxNorm_g; // TeV^-1 cm^-2 s^-1, at Energy eNorm_g;
  double energyNorm_g;    // TeV
  double alpha_g;    // (positive) spectral index
  double energyCut_g;     // TeV

  //    fluxNorm_g = k_g * (energyNorm_g)^{-alpha_g} *
  //                 exp( -sqrt( energyNorm_g / energyCut_g ) )



  // MGRO J1908+06
  fluxNorm_g = 8.8e-15.;
  energyNorm_g = 20.;
  alpha_g = 2.;
  energyCut_g = 300.;

  FormulaFlux fluxMGRO_J1908 = 
    MilagroFluxCalculator_nu(fluxNorm_g, energyNorm_g, alpha_g, energyCut_g);
  EquatorialDeg srcMGRO_J1908( 287.27 , 6.28 );

  FormulaFlux fluxMGRO_J1908_gamma = 
    MilagroFluxCalculator_gamma(fluxNorm_g,energyNorm_g,alpha_g,energyCut_g);


  // MGRO J1852+01
  fluxNorm_g = 5.7e-14.;
  energyNorm_g = 12.;
  alpha_g = 2.;
  energyCut_g = 300.;

  FormulaFlux fluxMGRO_J1852 = 
    MilagroFluxCalculator_nu(fluxNorm_g, energyNorm_g, alpha_g, energyCut_g);
  EquatorialDeg srcMGRO_J1852( 283.12 , 0.51 );


  // MGRO J2019+37
  fluxNorm_g = 8.7e-15.;
  energyNorm_g = 20.;
  alpha_g = 2.;
  energyCut_g = 300.;

  FormulaFlux fluxMGRO_J2019 = 
    MilagroFluxCalculator_nu(fluxNorm_g, energyNorm_g, alpha_g, energyCut_g);
  EquatorialDeg srcMGRO_J2019( 304.83 , 36.83 );

}
