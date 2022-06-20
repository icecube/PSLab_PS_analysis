#include "GeneralFC.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// POISSON:  NUMERICAL RECIPES
// #include "llh_general.h"

// POISSON:  ROOT
#include "TMath.h"




void GeneralFC::Load(char *dir, unsigned nSigLoadTotal, double sigMax) {

  // All previously stored data gets cleared

  rawData_.clear();
  matrix_.clear();
  maxes_.clear();

  ill_.clear();
  iul_.clear();


  // Set up 2-d array for raw data

  nSigLoadTotal_ = nSigLoadTotal;
  sigMax_ = sigMax;

  rawData_.resize(nSigLoadTotal_+1);
  for (unsigned i=0; i<rawData_.size(); ++i) { rawData_[i].resize(nBinsX_,0); }

  printf("Reading %d files in %s\n",nSigLoadTotal_+1,dir);
  LoadRawData(dir);


  //Build 2D probability matrix, interpolating between integer nSig

  matrix_.resize(nBinsY_);
  for (unsigned i=0; i<matrix_.size(); ++i) { matrix_[i].resize(nBinsX_,0); }

  maxes_.resize(nBinsX_,0);

  printf("Building probability matrix with interpolation\n");
  for (int i = 0; i < nBinsY_; i++) {
    double meanSig = sigMax_*i/nBinsY_;

    for (int k = 0; k < nSigLoadTotal_+1; k++) {
      for (int j = 0; j < nBinsX_; j++) { 
	matrix_[i][j] += rawData_[k][j] * Poisson(k,meanSig);
      }
    }

    for (int j = 0; j < nBinsX_; j++) {
      if (matrix_[i][j] > maxes_[j]) maxes_[j] = matrix_[i][j];
    }
  }


  // INCORPORATE UNCERTAINTY ON SIGNAL EFFICIENCY, IF APPLICABLE

  if (sigUnc_ > 0.) {

    vector<vector<double> > corr;
    corr.resize(nBinsY_);
    for (unsigned i=0; i<corr.size(); ++i) { corr[i].resize(nBinsY_,0); }

    for (int i = 0; i < nBinsY_; i++) {
      double iiunc2 = i*i*sigUnc_*sigUnc_;
      for (int j = 0; j < nBinsY_; j++) {
	corr[i][j] = exp(-0.5*(i-j)*(i-j)/iiunc2);
      }
    }

    vector<vector<double> > matrix2;
    matrix2.resize(nBinsY_);
    for (unsigned i=0; i<matrix2.size(); ++i) { matrix2[i].resize(nBinsX_,0); }
    
    // Computing new probabilities for uncertainty

    double probsum;
    // i == 0 case
    for (int j = 0; j < nBinsX_; j++) matrix2[0][j] = matrix_[0][j];

    for (int i = 1; i < nBinsY_; i++) {
      for (int j = 0; j < nBinsX_; j++) {
	if (sigUnc_ == 0.) matrix2[i][j] = matrix_[i][j];
	else {
	  probsum = 0.;
	  for (int k = 0; k < nBinsY_; k++) {
	    matrix2[i][j] += matrix_[k][j]*corr[i][k];
	    probsum += corr[i][k];
	  }
        matrix2[i][j] /= probsum;
	}
      }
    }

    // Replace original matrix with new matrix having uncertainties included

    for (int i = 0; i < nBinsY_; i++) {
      for (int j = 0; j < nBinsX_; j++) {
	matrix_[i][j] = matrix2[i][j];
      }
    }
  }

}
    



void GeneralFC::BuildConfidenceBand(double cl) {
  if (matrix_.size() == 0) {
    printf("Error: data has not yet been loaded properly.\n");
    return;
  }

  ill_.clear();
  iul_.clear();

  double probsum, max, lr_upper, lr_lower, lr_avg;
  int ul[nBinsY_], ll[nBinsY_], sel_bin, avg_cnt;
  sel_bin = 0;
    
  //Find highest likelihood ratio
  for (int i = 0; i < nBinsY_; i++) {
    max = probsum = 0.;
    for (int j = 0; j < nBinsX_; j++) {
      lr_avg = 0.;
      avg_cnt = 0;
      if (maxes_[j] != 0. && matrix_[i][j]/maxes_[j] > max) {
	lr_avg += matrix_[i][j]/maxes_[j];
	avg_cnt++;
	for (int k = j-1; k >= 0; k--) {
	  if (maxes_[k] != 0) {
	    lr_avg += matrix_[i][k]/maxes_[k];
	    avg_cnt++;
	    break;
	  }
	}
	for (int k = j+1; k < nBinsX_; k++) {
	  if (maxes_[k] != 0) {
	    lr_avg += matrix_[i][k]/maxes_[k];
	    avg_cnt++;
	    break;
	  }
	}
	
	if (lr_avg/avg_cnt > max) {
	  max = lr_avg/avg_cnt;
	  sel_bin = j;
	}
      }
    }
    

    probsum += matrix_[i][sel_bin];
    ul[i] = sel_bin;
    ll[i] = sel_bin;
    while (probsum < cl) {
      if (maxes_[ll[i]-1] == 0. && ll[i] > 0) {
	ll[i]--;
	continue;
      } else if (maxes_[ul[i]+1] == 0. && ul[i] < nBinsX_-1) {
	ul[i]++;
	continue;
      }
      
      lr_lower = ll[i] == 0 ? 0. : matrix_[i][ll[i]-1]/maxes_[ll[i]-1];
      lr_upper = ul[i] == nBinsX_-1 ? 0. : matrix_[i][ul[i]+1]/maxes_[ul[i]+1];
	
      //Expand side of acceptance region with higher likelihood ratio
      if (lr_lower > lr_upper) {
	ll[i]--;
	probsum += matrix_[i][ll[i]];
      } else if (ul[i] < nBinsX_-1) {
	ul[i]++;
	probsum += matrix_[i][ul[i]];
      } else {
	printf("Error: sig=%f (binY=%d) exceeded X with prob=%f < %f C.L.\n",
	       GetSigFromBinY(i),i,probsum,cl);
	break;
      }
    }

  }

  //Increase coverage to guarantee monotonicity
  for (int i = nBinsY_-1; i > 0; i--) {
    if (ll[i] < ll[i-1]) ll[i-1] = ll[i];
  }
  for (int i = 0; i < nBinsY_-1; i++) {
    if (ul[i+1] < ul[i]) ul[i+1] = ul[i];
  }

  //Invert the band -- add +1 to ul since ul is inclusive
  ill_.resize(nBinsX_);
  iul_.resize(nBinsX_);
  double meanSig = 0.;
  for (int j = 0; j < ll[0]; j++) { iul_[j] = 0; }
  for (int j = 0; j < ul[0]+1; j++) { ill_[j] = 0; }
  for (int i = 1; i < nBinsY_; i++) {
    meanSig = sigMax_*i/nBinsY_;
    for (int j = ll[i-1]; j < ll[i]; j++) { iul_[j] = meanSig; }
    for (int j = ul[i-1]+1; j < ul[i]+1; j++) { ill_[j] = meanSig; }
  }
  for (int j = ll[nBinsY_-1]; j < nBinsX_; j++) { iul_[j] = meanSig; }
  for (int j = ul[nBinsY_-1]+1; j < nBinsX_; j++) { ill_[j] = meanSig; }
    
  if (optPrint_) {
    for (int i = 0; i < nBinsX_; i++) {
      printf("BAND %d %f %f %f\n", i, (upLimX_-lowLimX_)*i/nBinsX_ + lowLimX_, ill_[i], iul_[i]);
    }
  }
}




double GeneralFC::GetUpperLimit(double testStat) {
  return iul_[ GetBinXFromTestStat(testStat) ];
}

double GeneralFC::GetLowerLimit(double testStat) {
  return ill_[ GetBinXFromTestStat(testStat) ];
}


double GeneralFC::GetMedianUpperLimit() {
  return GetUpperLimit( GetMedianOfNullHypothesis() );
}

double GeneralFC::GetMedianLowerLimit() {
  return GetLowerLimit( GetMedianOfNullHypothesis() );
}


double GeneralFC::GetAverageUpperLimit() {
  double avUpLim = 0.;
  for (int i = 0; i < nBinsX_; i++) { avUpLim += iul_[i]*rawData_[0][i]; }
  return avUpLim;
}

double GeneralFC::GetAverageLowerLimit() {
  double avUpLim = 0.;
  for (int i = 0; i < nBinsX_; i++) { avUpLim += ill_[i]*rawData_[0][i]; }
  return avUpLim;
}
    



TH2D* GeneralFC::GetHistRawData() {
  TH2D *hraw = new TH2D("hraw","hraw",
			nBinsX_,lowLimX_,upLimX_,
			nSigLoadTotal_+1,0,nSigLoadTotal_+1);// includes 0 bin
  for (int ix=0; ix<nBinsX_; ++ix) {
    for (int iy=0; iy<nSigLoadTotal_+1; ++iy) {
      double thisVal = rawData_[iy][ix];
      hraw->SetBinContent(ix+1,iy+1,thisVal);
    }
  }
  return hraw;
}

TH2D* GeneralFC::GetHistMatrix() {
  TH2D *hmatrix = new TH2D("hmatrix","hmatrix",
		       nBinsX_,lowLimX_,upLimX_,
		       nBinsY_,0,sigMax_);
  for (int ix=0; ix<nBinsX_; ++ix) {
    for (int iy=0; iy<nBinsY_; ++iy) {
      double thisVal = matrix_[iy][ix];
      hmatrix->SetBinContent(ix+1,iy+1,thisVal); // 0 is underflow bin of TH2D
    }
  }
  return hmatrix;
}

TH1D* GeneralFC::GetHistLowerLimit() {
  TH1D* hll = new TH1D("hll","",nBinsX_,lowLimX_,upLimX_);
  for (int i=0; i< nBinsX_; i++) { hll->SetBinContent(i+1,ill_[i]); }
  return hll;
}

TH1D* GeneralFC::GetHistUpperLimit() {
  TH1D* hul = new TH1D("hul","",nBinsX_,lowLimX_,upLimX_);
  for (int i=0; i< nBinsX_; i++) { hul->SetBinContent(i+1,iul_[i]); }
  return hul;
}








double GeneralFC::GetMedianOfNullHypothesis() {
  int medianBinX = nBinsX_;
  double sum = 0.;
  for (int j = 0; j < nBinsX_; j++) {
    sum += rawData_[0][j];
    if (sum>=0.5) {
      medianBinX = j;
      break;
    }
  }
  if (sum<0.5) { 
    printf("Error: null hypothesis not weighted sufficiently to calculate"
	   " median.\n");
  }
  return GetTestStatFromBinX(medianBinX);
}



void GeneralFC::LoadRawData(char *dir) {

  for (int i = 0; i < nSigLoadTotal_+1; i++) {
    char filename[1000];
    sprintf(filename, "%s/nSig_%03d.txt", dir, i);
    FILE *f = fopen(filename, "r");
    if (!f) {
      fprintf(stderr, "File not found: %s\n", filename);
      exit(1);
    }

    int cnt = 0;
    for (;;) {
      double testStat;
      // fill testStat from file, but break once end of file is reached 
      if ( LoadTestStat(f, testStat) == false) { break; } ;
      rawData_[i][ GetBinXFromTestStat(testStat) ] += 1.;
      cnt++;
    }
    fclose(f);

    for (int j = 0; j < nBinsX_; j++) { rawData_[i][j] /= cnt; }
    if (optPrint_) { printf("Loaded %s: %d entries\n", filename, cnt); }
  }
}



bool GeneralFC::LoadTestStat(FILE *f, double& testStat) {
  double val;
  // Read first double, then skip everything else until end of line
  if (fscanf(f, "%lg%*[^\n]\n", &val) == EOF) { return false; }
  // Perform calculation, if any, here:
  testStat = val;
  return true;
}

bool SignedSqrtFC::LoadTestStat(FILE *f, double& testStat) {
  double val;
  // Read first double, then skip everything else until end of line
  if (fscanf(f, "%lg%*[^\n]\n", &val) == EOF) { return false; }
  // Perform calculation, if any, here:
  testStat = sqrt(2*fabs(val));
  if (testStat < 0.) { testStat *= -1.; }
  return true;
}

 

double GeneralFC::Poisson(int k, double mean) {
  double poisson_p;
  if (mean==0 && k>0) { 
    poisson_p = 0.;
  } else {
    // NUMERICAL RECIPES
    //    poisson_p = poisson_prob(mean,"eq",k);
    // ROOT
    poisson_p = TMath::Poisson(k,mean);
    if (poisson_p != poisson_p) {
      poisson_p = 0.;  // NaN can happen if too small
    }
  }
  return poisson_p;
}


int GeneralFC::GetBinXFromTestStat(double testStat) {
  int binX = int(nBinsX_*(testStat-lowLimX_)/(upLimX_-lowLimX_));
  if (binX<0) {
    printf("Error: testStat=%f < lowerLimit=%f\n",testStat,lowLimX_);
    binX=0;
  }
  if (binX>=nBinsX_) {
    printf("Error: testStat=%f >= upperLimit=%f\n",testStat,upLimX_);
    binX=nBinsX_-1;
  }
  return binX;
}


double GeneralFC::GetTestStatFromBinX(int binX) {
  // add 0.5 to get value of center of bin
  double binValueX = (upLimX_-lowLimX_)*(binX+0.5)/nBinsX_+lowLimX_;
  return binValueX;
}

double GeneralFC::GetSigFromBinY(int binY) {
  // add 0.5 to get value of center of bin
  double binValueY = (sigMax_-0.)*(binY+0.5)/nBinsY_+0.;
  return binValueY;
}
