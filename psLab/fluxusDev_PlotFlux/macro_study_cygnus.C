
{
  gROOT->ProcessLine(".L FluxTools.C");
  gROOT->ProcessLine(".x FluxList.C");

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(0.05);
  
  double sizeFactor = 1.0;
  TCanvas *can = new TCanvas("canFluxCurves","canFluxCurves", 20,20, 
			     600*sizeFactor, 600*sizeFactor);
  can->Divide(1,2,0.001,0.001);

  FluxPlotManager fpm, fpmEvents;

  bool optConnect = true;


  plotLevel = 1;


  // IC40 FULL
  if (plotLevel >= 1) {
    read_flux_curves("results/fluxDecades_ic40_full_Em2_sens_E_dec06.txt", 
    		     fpm, fpmEvents, kRed+1, 3, 1);
    read_flux_curves("results/fluxDecades_ic40_full_Em2_sens_E_dec41.txt", 
    		     fpm, fpmEvents, kBlue+2, 3, 1);
  }


 // VERITAS: M82 - arXiv:0911.0873
  if (0) {
    double phi_int_700GeV = 3.7e-13; // cm-2 s-1
    double gamma = -2.5; // best-fit value
    double factor = IntegralToDifferentialFluxConstant(0.7e12,1e9,gamma);
    double phi_0 = factor * phi_int_700GeV;

    PowerLawFlux pflux_M82_VERITAS(phi_0,gamma);
    g = GraphFlux(pflux_M82_VERITAS, 0, 1e3, 5e3);
    g->SetLineColor(kGray);
    g->SetLineWidth(5);
    fpm.AddGraph(*g);

    // Corresponding IC40 (6 months) upper limit

    double factor = ScaleFlux(1e12,1e9,1e12,1e9,-2.5);
    double phi_0 = factor * 9e-11;
    PowerLawFlux pflux_M82_VERITAS_IC40_half_sen(phi_0,-2.5);
    g = GraphFlux(pflux_M82_VERITAS_IC40_half_sen, 0, 3e2, 1e5);
    g->SetLineColor(kGray);
    g->SetLineWidth(3);
    g->SetLineStyle(9);
    fpm.AddGraph(*g);    
  }

  // CYGNUS X-3

  if (plotLevel >= 2) {
    // Fermi:  Science 326, 1512

    //    double ergToGeV = (1000./1.602.);
    //    double phi_int_100MeV = ergToGeV * 4.0e-10;
    double phi_int_100MeV = 1.19e-6;
    double gamma = -2.7;
    double factor = IntegralToDifferentialFluxConstant(1e8, 1e9, gamma);
    double phi_0 = factor * phi_int_100MeV;
    cout << "Cygnus X-3:  " << phi_0 << endl;
    PowerLawFlux pflux_CygX3_Fermi(phi_0, gamma);
    g = GraphFlux(pflux_CygX3_Fermi, 0, 1e-1, 1e2);
    g->SetLineColor(kBlue+2);
    g->SetLineWidth(3);
    g->SetLineStyle(1);
    fpm.AddGraph(*g);
  }

  if (plotLevel >= 3 ) {
    // extension beyond fermi range    
    g = GraphFlux(pflux_CygX3_Fermi, 0, 1e2, 1e5);
    g->SetLineColor(kBlue+2);
    g->SetLineWidth(3);
    g->SetLineStyle(2);
    fpm.AddGraph(*g); 
   
    // ic40 upper limit

    PowerLawFlux pflux_CygX3_Fermi_IC40_full(6.6e-6,-2.7);
    g = GraphFlux(pflux_CygX3_Fermi_IC40_full, 0, 3e2, 1e5);
    g->SetLineColor(kBlue+2);
    g->SetLineWidth(3);
    g->SetLineStyle(9);
    fpm.AddGraph(*g);    
  }

  if (plotLevel >= 4) {
    gROOT->ProcessLine(".x FluxList_Milagro.C");

    /*
    g = GraphFlux(fluxMGRO_J2019, 0, 1e2, 1e7);
    g->SetLineColor(kOrange-6);
    g->SetLineWidth(3);
    fpm.AddGraph(*g);
    read_flux_curves("results/flux_ic40_MGRO_J2019_sens.txt", 
		     fpm, fpmEvents, kOrange-6,   2, 9);
    read_flux_curves("results/flux_ic40_MGRO_J2019_disc.txt", 
		     fpm, fpmEvents, kOrange-6,   2, 1);
    */

    g = GraphFlux(fluxMGRO_J1908, 0, 3e2, 1e6);
    g->SetLineColor(kRed+1);
    g->SetLineWidth(3);
    g->SetLineStyle(2);
    fpm.AddGraph(*g);
    //    read_flux_curves("results/flux_ic40_MGRO_J1908_sens.txt", 
    //		     fpm, fpmEvents, kRed+1,   3, 9);
    //    read_flux_curves("results/flux_ic40_MGRO_J1908_disc.txt", 
    //		     fpm, fpmEvents, kCyan+2,   2, 1);

    //    g = GraphFlux(fluxMGRO_J1908_gamma, 0, 1e2, 1e7);
    //    g->SetLineColor(kGray);
    //    g->SetLineWidth(3);
    //    fpm.AddGraph(*g);


    g = GraphFlux(fluxMGRO_J1852, 0, 3e2, 1e6);
    g->SetLineColor(kCyan+3);
    g->SetLineWidth(3);
    g->SetLineStyle(2);
    fpm.AddGraph(*g);
    //    read_flux_curves("results/flux_ic40_MGRO_J1852_sens.txt", 
    //		     fpm, fpmEvents, kBlack,   2, 9);
    //    read_flux_curves("results/flux_ic40_MGRO_J1852_disc.txt", 
    //		     fpm, fpmEvents, kBlack,   2, 2);
  }

  // Jan-Patrick's first galactic center sensitivities
  if (0) {
    TGraph *g;
    FormulaFlux f("1000*2.26e-9*pow(x,-2)");
    g = GraphFlux(f,0,7e3,700e3);
    fpm.AddGraph(*g);

    FormulaFlux f("0.001*1e-7*pow(x/1000,-2)*exp(-x/1e4)");
    g = GraphFlux(f,0,1e2,1e6);
    fpm.AddGraph(*g);
  }



  can->cd(1);
  fpmEvents.yPower = 0;
  fpmEvents.Plot();
  gPad->SetLogy(0);
  SetPadTitle(";;Events");
  SetPadYMin(0);



  can->cd(2);
  fpm.logEMin = -1;
  fpm.Plot();
  gPad->SetGrid(1);


  TCanvas *canOut = new TCanvas("canOut","canOut",20,20,900,500);
  canOut->cd();
  gPad->SetRightMargin(0.05);
  fpm.yOutputUnits = 1e12;
  fpm.Plot();
  SetPadYMin(1e-14);
  SetPadYMax(1.5e-8);
  gPad->SetGrid();
  fpm.hFrame.GetYaxis()->SetTitleSize(0.045);
  fpm.hFrame.GetYaxis()->SetLabelSize(0.045);

}
