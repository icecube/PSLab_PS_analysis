
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

  if (1) {
    read_flux_curves("results/fluxDecades2_ic22_Em2_disc_E_dec06.txt", 
		     fpm, fpmEvents, kBlue+1, 1, 1);
    read_flux_curves("results/fluxDecades_ic40_final_Em2_disc_E_dec6.txt", 
		     fpm, fpmEvents, kRed+2, 1, 1);
    read_flux_curves("results/fluxDecades_ic40_02years_Em2_disc_E_dec6.txt", 
    		     fpm, fpmEvents, kRed+1, 1, 1);
    read_flux_curves("results/fluxDecades_ic40_10years_Em2_disc_E_dec6.txt", 
    		     fpm, fpmEvents, kRed, 1, 1);

    read_flux_curves("results/flux_ic40_10years_MGRO_J1852_disc.txt", 
    		     fpm, fpmEvents, kRed, 3, 1);


  }


  if (0) {
    read_flux_curves("../sandbox_scripts/ic22/results/flux_SS433_disc.txt", 
		     fpm, fpmEvents, kRed+1, 1, 1);
  }
  if (0) {
    read_flux_curves("results/fluxCutoff_ic40A5_161days_Em2_disc_E_dec6.txt", 
		     fpm, fpmEvents, kRed+1, 1, 1);
    read_flux_curves("results/fluxCutoff_ic40A5_161days_Em2_disc_E_dec6.txt", 
		     fpm, fpmEvents, kRed+1, 3, 1, optConnect, +1);

    read_flux_curves("results/fluxCutoff_ic40A5_161days_Em2_disc_noE_dec6.txt",
		     fpm, fpmEvents, kBlue+1, 1, 1);
    read_flux_curves("results/fluxCutoff_ic40A5_161days_Em2_disc_noE_dec6.txt",
		     fpm, fpmEvents, kBlue+1, 3, 1, optConnect, +1);
  }

  if (0) {
    //  read_flux_curves("results/flux_ic22_Em2_disc_E_dec06.txt", 
    //		   fpm, fpmEvents, kMagenta+2, 3, 1);
    //  read_flux_curves("results/flux_ic22_Em3_disc_E_dec06.txt", 
    //  		   fpm, fpmEvents, kBlue+1, 3, 1);
    //    read_flux_curves("results/flux_ic22_Em1.5_disc_E_dec06.txt", 
    //  		   fpm, fpmEvents, kRed+1, 3, 1);

  
  //  read_flux_curves("results/fluxDecades1_ic22_Em2_disc_E_dec06.txt",
  //		   fpm, fpmEvents, kCyan+2, 3, 1);

    //  read_flux_curves("results/fluxDecades2_ic22_Em2_disc_E_dec06.txt", 
    //		   fpm, fpmEvents, kCyan+2, 1, 1, optConnect);
  //  read_flux_curves("results/fluxDecades4_ic22_Em2_disc_E_dec06.txt",
  //		   fpm, fpmEvents, kCyan+2, 3, 1);
  }

 if (0) {
    read_flux_curves("results/fluxDecades2_ic40_Em2_disc_E_dec60.txt", 
  	    fpm, fpmEvents, kGreen+3,   3, 1, optConnect);
  read_flux_curves("results/fluxDecades2_ic40_Em2_disc_E_dec30.txt", 
		   fpm, fpmEvents, kGreen-2,   3, 1, optConnect);
  read_flux_curves("results/fluxDecades2_ic40_Em2_disc_E_dec06.txt",
		   fpm, fpmEvents, kCyan+2,  3, 1, optConnect);

  read_flux_curves("results/fluxDecades2_ic40_Em2_disc_E_dec-08.txt", 
		   fpm, fpmEvents, kBlue-4, 3, 1, optConnect);
    read_flux_curves("results/fluxDecades2_ic40_Em2_disc_E_dec-30.txt", 
  		   fpm, fpmEvents, kViolet+2,   3, 1, optConnect);
    read_flux_curves("results/fluxDecades2_ic40_Em2_disc_E_dec-60.txt", 
  		   fpm, fpmEvents, kViolet-2, 3, 1, optConnect);
 }
 if (0) {
  read_flux_curves("results/fluxDecades2_ic40_Em2_disc_E_dec30_WH.txt", 
		   fpm, fpmEvents, kGreen-2,   3, 2, optConnect);
  read_flux_curves("results/fluxDecades2_ic40_Em2_disc_E_dec06_WH.txt",
		   fpm, fpmEvents, kCyan+2,  3, 2, optConnect);
 }
 if (0) {
  read_flux_curves("results/fluxDecades2_ic40_Em2_disc_E_dec30.txt", 
		   fpm, fpmEvents, kGreen-2,   3, 1, optConnect);
  read_flux_curves("results/fluxDecades2_ic40_Em2_disc_E_dec30_thin0.003.txt", 
		   fpm, fpmEvents, kGreen-2,   2, 1, optConnect);
  read_flux_curves("results/fluxDecades2_ic40_Em2_disc_E_dec30_WH.txt", 
		   fpm, fpmEvents, kGreen-2,   3, 2, optConnect);
  read_flux_curves("results/fluxDecades2_ic40_Em2_disc_E_dec30_WH_thin0.003.txt", 
		   fpm, fpmEvents, kGreen-2,   2, 2, optConnect);
 }



  //  read_flux_curves("results/fluxDecades2_ic22_Em2_disc_E_dec06_thin0.1.txt", fpm, fpmEvents, kBlue-4, 3, 1);
  //  read_flux_curves("results/fluxDecades2_ic22_Em2_disc_E_dec06_thin0.01.txt", fpm, fpmEvents, kBlue-7, 3, 1);

  if (0) {
    read_flux_curves("results/fluxDecades2_ic22_Em2_disc_E_dec30.txt", fpm, fpmEvents, kGreen-2, 1, 1, optConnect);
    read_flux_curves("results/fluxDecades2_ic22_Em2_disc_E_dec60.txt", fpm, fpmEvents, kGreen+3, 1, 1, optConnect);
  }

  if (1) {
    // VERITAS: M82 - arXiv:0911.0873

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


  if (1) {
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


    g = GraphFlux(fluxMGRO_J1908, 0, 1e2, 1e7);
    g->SetLineColor(kCyan+2);
    g->SetLineWidth(3);
    fpm.AddGraph(*g);
    read_flux_curves("results/flux_ic40_MGRO_J1908_sens.txt", 
		     fpm, fpmEvents, kCyan+2,   2, 9);
    read_flux_curves("results/flux_ic40_MGRO_J1908_disc.txt", 
		     fpm, fpmEvents, kCyan+2,   2, 1);

    //    g = GraphFlux(fluxMGRO_J1908_gamma, 0, 1e2, 1e7);
    //    g->SetLineColor(kGray);
    //    g->SetLineWidth(3);
    //    fpm.AddGraph(*g);
    */

    g = GraphFlux(fluxMGRO_J1852, 0, 1e2, 1e7);
    g->SetLineColor(kBlack);
    g->SetLineWidth(2);
    fpm.AddGraph(*g);
    read_flux_curves("results/flux_ic40_MGRO_J1852_sens.txt", 
		     fpm, fpmEvents, kBlack,   2, 9);
    //    read_flux_curves("results/flux_ic40_MGRO_J1852_disc.txt", 
    //		     fpm, fpmEvents, kBlack,   2, 2);
  }

  if (1) {
    gROOT->ProcessLine(".x FluxList_RXJ_Blasi.C");
    g = GraphFlux(fluxRXJ_Blasi, 0, 1e2, 1e6);
    g->SetLineColor(kViolet-1);
    g->SetLineWidth(3);
    fpm.AddGraph(*g);
    read_flux_curves("results/flux_ic40_RXJ_Blasi_At_Crab_sens.txt", 
		     fpm, fpmEvents, kViolet-1,   2, 9);
    //    read_flux_curves("results/flux_ic40_RXJ_Blasi_At_Crab_disc.txt", 
    //		     fpm, fpmEvents, kViolet-1,   2, 1);
  }


  if (1) { 
    g = GraphFlux(flux3C279, 0, 1e4, 1e9);
    g->SetLineColor(kBlue-4);
    g->SetLineWidth(3);
    fpm.AddGraph(*g);

    //    read_flux_curves("results/flux_ic22_3C279_sens_E.txt", 
    //	     fpm, fpmEvents, kBlue-4,   1, 1);
      read_flux_curves("results/flux_ic40_3C279_sens_E.txt", 
		     fpm, fpmEvents, kBlue-4,   2, 9);
      //      read_flux_curves("results/flux_ic40_3C279_disc.txt", 
      //		     fpm, fpmEvents, kBlue-4,   2, 1);
  }

  if (1) {
    g = GraphFlux(fluxCrab, 0, 1e2, 6e4);
    g->SetLineColor(kGreen-2);
    g->SetLineWidth(3);
    fpm.AddGraph(*g);

    read_flux_curves("results/flux_ic22_Crab_sens_E.txt", 
		     fpm, fpmEvents, kGreen-2,   1,1);

    read_flux_curves("results/flux_ic40_Crab_sens_E.txt", 
		     fpm, fpmEvents, kGreen-2,   2,9);
    
       //       read_flux_curves("results/flux_ic40_Crab_sens_E_WH.txt", 
       //     fpm, fpmEvents, kGreen-2,  3, 2);
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


  if (0) {
    TGraph *g;
    g = GraphFlux(fluxMGRO_J1908_old, 0, 1e2, 1e7);
    g->SetLineColor(kCyan+2);
    g->SetLineWidth(3);
    fpm.AddGraph(*g);

    read_flux_curves("results/flux_ic22_MGRO_J1908_old_sens_E.txt", 
		     fpm, fpmEvents, kCyan+2,   1, 1);
    read_flux_curves("results/flux_ic40_MGRO_J1908_old_sens_E.txt", 
		     fpm, fpmEvents, kCyan+2,   2, 9);
    read_flux_curves("results/flux_ic40_MGRO_J1908_old_sens_E_WH.txt", 
		     fpm, fpmEvents, kCyan+2,   3, 2);
  }



  can->cd(1);
  fpmEvents.yPower = 0;
  fpmEvents.Plot();
  gPad->SetLogy(0);
  SetPadTitle(";;Events");
  SetPadYMin(0);



  can->cd(2);
  fpm.Plot();
  gPad->SetGrid(1);


  TCanvas *canOut = new TCanvas("canOut","canOut",20,20,900,500);
  canOut->cd();
  gPad->SetRightMargin(0.05);
  fpm.yOutputUnits = 1e12;
  fpm.Plot();
  SetPadYMin(1e-15);
  SetPadYMax(1.5e-7);
  gPad->SetGrid();
  fpm.hFrame.GetYaxis()->SetTitleSize(0.045);
  fpm.hFrame.GetYaxis()->SetLabelSize(0.045);

}
