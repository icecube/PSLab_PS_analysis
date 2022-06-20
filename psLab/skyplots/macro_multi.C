
{
  int startID = 938;
  int stopID = 938;

  for (int id=startID; id<=stopID; ++id) {
    char inputHist[1000];
    char inputResults[1000];
    char outputImage[1000];
    char outputText[1000];
    sprintf(inputHist,
	    "/net/user/cfinley/npx/ic40_v1/"
	    "AllSky_CutA5_6months_final/hists/trial_%08d_skymaps.root", id);
    sprintf(inputResults,
	    "/net/user/cfinley/npx/ic40_v1/"
	    "AllSky_CutA5_6months_final/results/trial_%08d_results.txt", id);

    {
      FILE *fp = fopen(inputResults,"r");
      int ranSeed;
      double llhRatio, log10p, raDeg, decDeg, ns, gamma;
      fscanf(fp,"%d %lf %lf %lf %lf %lf %lf",
	     &ranSeed,
	     &llhRatio,
	     &log10p,
	     &raDeg,
	     &decDeg,
	     &ns,
	     &gamma);
      fclose(fp);

      sprintf(outputText,
	      "scrambled ID: %4d       Hottest Spot:   log #lambda = %5.2f    "
	      "-log_{10}p = %4.2f    r.a. = %6.2f    dec = %6.2f    "
              "nSrc = %4.1f    gamma = %3.1f",
	      -ranSeed, llhRatio, log10p, raDeg, decDeg, ns, gamma);
    }


    sprintf(outputImage,"ic40_6months_scrambled_skymap_%08d.png",id);
    cout << "Generating: " << outputImage << endl;
    cout << outputText << endl;


    TFile *f = new TFile(inputHist);

    TH2D *hInput = hAllSkyFine;

    Projection *pro = new HammerAitoffProjection();
    //  Projection *pro = new FlatProjection();

    gROOT->ProcessLine(".x SimpleSkyPlot.C");

    gPad->SetBottomMargin(0.077);

    TPaveText tp(0.01, 0.01, 0.99, 0.075, "NDC");
    tp.AddText(outputText);
    tp.SetFillColor(0);
    tp.SetBorderSize(0);
    tp.Draw();

    can->Update();
    can->SaveAs(outputImage);

    f->Close();
  }

}
