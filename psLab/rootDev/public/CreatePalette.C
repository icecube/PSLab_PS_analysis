#include "TColor.h"
#include "TROOT.h"
#include "TStyle.h"


void CreatePalette(Int_t pal=0) 
{

  float  r = 0.;
  float  g = 0.;
  float  b = 0.;


  // ################   HESS Color Palette   ###################

  if (pal==0) {

    const int n=128;
    int offset = 1200;
    
    Int_t *palette = new Int_t [n];
    
    double step1 = 30.;
    double step2 = 60.;
    double step3 = 90.;
    double step4 = 120;
    
    for(int i=0; i<n; i++) {
      
      if (i <= step1) {
        r = 0. + 0.*i/step1;
        g = 0.;// + 0.1*i/32.;
        b = 0.25 + 0.3*i/step1;
      } else if (i <= step2) {
        r = 0. + 0.3*(i-step1)/(step2-step1);
        g = 0. + 0.1*(i-step1)/(step2-step1);
        b = 0.55 + 0.3*(i-step1)/(step2-step1);
      } else if (i <= step3) {
        r = 0.3 + 0.7*(i-step2)/(step3-step2);
        g = 0.1 - 0.1*(i-step2)/(step3-step2);
        b = 0.85 - 0.85*(i-step2)/(step3-step2);
      } else if (i <= step4) {
        r = 1.;
        g = (i-step3)/(step4-step3);
        b = 0.;
      } else {
        r = 1.;
        g = 1.;
        b = (i-step4-1)/(n-step4);
      }

      if (! gROOT->GetColor(offset+i)) {
	new TColor(offset+i,r,g,b,"user colour");
      }

      palette[i] = offset+i;
    }

    gStyle->SetPalette(n,palette);
  }    



  if (pal==1) {

    const int n=128;
    int offset = 1600;
    
    Int_t *palette = new Int_t [n];
    
    float  r,g,b;
    double step1 = 30.;
    double step2 = 60.;
    double step3 = 90.;
    double step4 = 120;
    
    for(int i=0; i<n; i++) {

      if (i < 10) {
	r = 1.;
	g = 1.;
	b = 1.;
      } else if (i <= step1) {
        r = 0. + 0.*i/step1;
        g = 0.;// + 0.1*i/32.;
        b = 0.25 + 0.3*i/step1;
      } else if (i <= step2) {
        r = 0. + 0.3*(i-step1)/(step2-step1);
        g = 0. + 0.1*(i-step1)/(step2-step1);
        b = 0.55 + 0.3*(i-step1)/(step2-step1);
      } else if (i <= step3) {
        r = 0.3 + 0.7*(i-step2)/(step3-step2);
        g = 0.1 - 0.1*(i-step2)/(step3-step2);
        b = 0.85 - 0.85*(i-step2)/(step3-step2);
      } else if (i <= step4) {
        r = 1.;
        g = (i-step3)/(step4-step3);
        b = 0.;
      } else {
        r = 1.;
        g = 1.;
        b = (i-step4-1)/(n-step4);
      }

      if (! gROOT->GetColor(offset+i)) {
	new TColor(offset+i,r,g,b,"user colour");
      }

      palette[i] = offset+i;
    }

    gStyle->SetPalette(n,palette);
  }    



  if (pal==3) {
    const int n=90;
    int offset = 2000;
    
    Int_t *palette = new Int_t [n];
    
    double step0 = 0.;
    double step1 = 30.;
    double step2 = 60.;
    double step3 = 90.;
    double frac;

    for(int i=0; i<n; i++) {
      if (i<=4) {
	r = 1.;
	b = 1.;
	g = 1.;
    }else if (i <= step1) {
	frac = (i-step0)/(step1-step0);
        r = 0.0 + 0.0*frac;
        g = 0.0 + 1.0*frac;
        b = 1.0 + 0.0*frac;
      } else if (i <= step2) {
	frac = (i-step1)/(step2-step1);
        r = 0. + 1.*frac;
        g = 1. + 0.*frac;
        b = 1. - 1.*frac;
      } else if (i <= step3) {
	frac = (i-step2)/(step3-step2);
        r = 1. + 0.*frac;
        g = 1. - 1.*frac;
        b = 0. + 0.*frac;
      }
      

      if (! gROOT->GetColor(offset+i)) {
	new TColor(offset+i,r,g,b,"user colour");
      } else {
	TColor *color = gROOT->GetColor(offset+i);
	color->SetRGB(r,g,b);
      }
      
      palette[i] = offset+i;
    }

    gStyle->SetPalette(n,palette);
  }    


  if (pal==4) {
    const int n=90;
    int offset = 2400;
    
    Int_t *palette = new Int_t [n];
    
    double step0 = 0.;
    double step1 = 30.;
    double step2 = 60.;
    double step3 = 90.;
    double frac;

    for(int i=0; i<n; i++) {
      
      if (i <= step1) {
	frac = (i-step0)/(step1-step0);
        r = 0. + 0.*frac;
        g = 0. + 1.*frac;
        b = 1. + 0.*frac;
      } else if (i <= step2) {
	frac = (i-step1)/(step2-step1);
        r = 0. + 1.*frac;
        g = 1. + 0.*frac;
        b = 1. - 1.*frac;
      } else if (i <= step3) {
	frac = (i-step2)/(step3-step2);
        r = 1. + 0.*frac;
        g = 1. - 1.*frac;
        b = 0. + 0.*frac;
      }
      

      if (! gROOT->GetColor(offset+i)) {
	new TColor(offset+i,r,g,b,"user colour");
      }
      
      palette[i] = offset+i;
    }

    gStyle->SetPalette(n,palette);
  }    



  if (pal==5) {
    const int n=90;
    int offset = 2600;
    
    Int_t *palette = new Int_t [n];
    
    double step0 = 0.;
    double step1 = 40.;
    double step2 = 70.;
    double step3 = 90.;
    //    double step4 = 120.;
    double frac;

    for(int i=0; i<n; i++) {
      if (i<=4) {
	r = 1.0;
	b = 1.0;
	g = 1.0;
    }else if (i <= step1) {
	frac = (i-step0)/(step1-step0);
        r = 1.0 - 1.0*frac;
        g = 1.0 - 0.6*frac;
        b = 1.0 + 0.0*frac;
      } else if (i <= step2) {
	frac = (i-step1)/(step2-step1);
        r = 0.0 + 1.0*frac;
        g = 0.4 - 0.1*frac;
        b = 1.0 - 0.7*frac;
      } else if (i <= step3) {
	frac = (i-step2)/(step3-step2);
        r = 1.0 + 0.0*frac;
        g = 0.3 + 0.7*frac;
        b = 0.3 + 0.0*frac;
      } 
      

      if (! gROOT->GetColor(offset+i)) {
	new TColor(offset+i,r,g,b,"user colour");
      } else {
	TColor *color = gROOT->GetColor(offset+i);
	color->SetRGB(r,g,b);
      }
      
      palette[i] = offset+i;
    }

    gStyle->SetPalette(n,palette);
  }    

  // SAME AS ABOVE, BUT NO WHITE AT BOTTOM

  if (pal==6) {
    const int n=90;
    int offset = 2700;
    
    Int_t *palette = new Int_t [n];
    
    double step0 = 0.;
    double step1 = 40.;
    double step2 = 70.;
    double step3 = 90.;
    //    double step4 = 120.;
    double frac;

    for(int i=0; i<n; i++) {
      if (i<0) {
	r = 1.0;
	b = 1.0;
	g = 1.0;
    }else if (i <= step1) {
	frac = (i-step0)/(step1-step0);
        r = 1.0 - 1.0*frac;
        g = 1.0 - 0.6*frac;
        b = 1.0 + 0.0*frac;
      } else if (i <= step2) {
	frac = (i-step1)/(step2-step1);
        r = 0.0 + 1.0*frac;
        g = 0.4 - 0.1*frac;
        b = 1.0 - 0.7*frac;
      } else if (i <= step3) {
	frac = (i-step2)/(step3-step2);
        r = 1.0 + 0.0*frac;
        g = 0.3 + 0.7*frac;
        b = 0.3 + 0.0*frac;
      } 
      

      if (! gROOT->GetColor(offset+i)) {
	new TColor(offset+i,r,g,b,"user colour");
      } else {
	TColor *color = gROOT->GetColor(offset+i);
	color->SetRGB(r,g,b);
      }
      
      palette[i] = offset+i;
    }

    gStyle->SetPalette(n,palette);
  }    


  /* GRAY SCALE */

  if (pal==7) {
    const int n=50;
    int offset = 3000;
    
    Int_t *palette = new Int_t [n];
    
    double step0 = 0.;
    double step1 = 50.;
    double frac;

    for(int i=0; i<n; i++) {
      if (i==0) {
	r = 1.0;
	g = 1.0;
	b = 1.0;
    } else if (i <= step1) {
	frac = (i-step0)/(step1-step0);
        r = 1.0 - 1.0*frac;
        g = 1.0 - 1.0*frac;
        b = 1.0 - 1.0*frac;
    } 
      

      if (! gROOT->GetColor(offset+i)) {
	new TColor(offset+i,r,g,b,"user colour");
      } else {
	TColor *color = gROOT->GetColor(offset+i);
	color->SetRGB(r,g,b);
      }
      
      palette[i] = offset+i;
    }

    gStyle->SetPalette(n,palette);
  }    



}

