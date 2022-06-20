#ifndef LLHTIMEDEP_BLOCKLEVEL_H_
#define LLHTIMEDEP_BLOCKLEVEL_H_

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TH1D.h"
#include "TRandom.h"

#include "llh/public/Time.h"
#include "llh/public/TimePdf.h"

#include "llh/public/MultiAnalysisFn.h"


class BlockLevel { // the basic piece for the Maximum Likelihood Blocks time pdf.
 public:
  double BBegin;
  double BLevel;
  double BDuration;
  BlockLevel() {}
  BlockLevel(double BBegin_, double BLevel_, double BDur_) {
    BBegin = BBegin_;
    BLevel = BLevel_;
    BDuration = BDur_;
  }
  void SetBlockLevel(double BBegin_, double BLevel_, double BDur_) {
    BBegin = BBegin_;
    BLevel = BLevel_;
    BDuration = BDur_;
  }
  double GetBlockBegin() { return BBegin; }
  double GetBlockLevel() { return BLevel; }
  double GetBlockDuration() { return BDuration; }
};





// About the Maximum Likelihood Blocks:
// We're using them to denoise lightcurves, it is a
// function which returns a series of times which are
// compatable with a constant flux.

// The output is a text file which has three columns:
// BeginTime(MJD) FluxLevel Duration(day)

// and that just gets read in here and stored.
// The Block pdf can have a threshold on emission and a time lag.

// There are two versions BlockTimePdf which has a simple cut on
// threshold and BlockTimePdf1, which makes the pdf proportional
// to the level above the threshold. I only use BlockTimePdf1 now.

// It still includes a histogram, which by default isn't created
// or set and is only meant for plotting the pdf afterwards.


//class BlockTimePdf : public TimePdf {
//  private:
//    vector<BlockLevel> levels_;
//    string infile;
//
// public:
//
//  BlockTimePdf() {
//    norm_ = 0.;
//    livetimeTotal_ = 0.;
//  }
//
//  BlockTimePdf(BlockLevel a){
//    levels_.push_back(a);
//    norm_ = norm_ + (a.GetBlockLevel() * a.GetBlockDuration());
//    livetimeTotal_ = livetimeTotal_ + a.GetBlockDuration();
//  }
//
//  TimePdf* Clone() const { return new BlockTimePdf(*this); }
//
//  void SetBlockLevels(string fname, double cut=0., double offset=0.) {
//    infile = fname;
//    norm_ = 0;
//    levels_.clear(); // resetting the blocks each time we call this
//    livetimeTotal_ = 0.;
//    ifstream fin;
//    fin.open(fname.c_str());
//    double blockbegin, blockdur, blocklev;
//    while (fin >> blockbegin) {
//      fin >> blocklev >> blockdur;
//      blockbegin += offset;
//      if (blocklev < cut) {blocklev = 0.;}
//      levels_.push_back(BlockLevel(blockbegin, blocklev, blockdur));
//      norm_ = norm_ + (blocklev * blockdur);
//      livetimeTotal_ = livetimeTotal_ + blockdur;
//    }
//    tmin_ = levels_[0].GetBlockBegin();
//    tmax_ = tmin_+livetimeTotal_;
//    fin.close();
//  }
//
//  void SetNorm() { cout << "Norm is set when we load the block file, \n so we're not doing anything here." << endl; }
//
//  double GetPdfValue(double MJD){
//    if (MJD<tmin_ ||MJD>tmax_) return 0;
//    for (unsigned int i=0; i <= levels_.size(); i++) {
//      if (levels_[i].BBegin <= MJD && (levels_[i].BBegin + levels_[i].BDuration) > MJD) {
//        return levels_[i].BLevel/norm_;
//      }
//    }
//    return 0;
//  }
//
////  double GetPdfValue(Time t)             { return GetPdfValue( t.GetMJD() ); }
////  double GetPdfValue(Time t, double lag) { return GetPdfValue( t.GetMJD()-lag ); }
//
//  void fillHisto(int nbins){
//    //this can be different for the block Pdf
//    int nbins_ = levels_.size();
//    double begin[nbins_];
//    for (int i=0; i<nbins_; i++){
//      begin[i] = levels_[i].BBegin;
//    }
//    histo_ = new TH1D("histo_","source",nbins_,begin);
//    double MJD;
//    for (int i=1; i <= nbins_; i++){ // bins go from 1 to n
//      MJD = levels_[1].BBegin;
//      histo_->SetBinContent(i,GetPdfValue(MJD));
//    }
//
// //    int nbins_=nbins;
////    double bintime = (tmax_ - tmin_) / nbins_;
////    double MJD;
////    histo_ = new TH1D("histo_","source",nbins_,tmin_,tmax_);
////    for (int i=1; i <= nbins_; i++){ // bins go from 1 to n
////      MJD = tmin_ + (i-0.5)*bintime;
////      histo_->SetBinContent(i,GetPdfValue(MJD));
////    }
//  }
//
//  TH1D * GetHisto() {return histo_;}
//
//  Time GenerateEventTime(){
//    double t = histo_->GetRandom();
//    return Time(t);
//  }
//
//  void CheckTimeBounds(double dataTmin, double dataTmax) {
//    // a bit of housekeeping above to keep things within [dataTmin, dataTmax]
//    // so the norm is more useful for use with multiple datasets.
//
//    // in this function it won't matter if I made something with a cut
//    // or offset or whatever to begin with, since we're dealing with
//    // the blocks themselves.
//
//    double norm1 = 0.;
//    double lTtemp = 0;
//    vector<BlockLevel> bTemp;
//    double blockbegin, blocklev, blockdur;
//    for (unsigned int i=0;i<levels_.size();i++) {
//      blockbegin = levels_[i].BBegin;
//      blocklev = levels_[i].BDuration;
//      blockdur = levels_[i].BLevel;
//      if ( blockbegin+blockdur <= dataTmin || blockbegin > dataTmax ) { continue; }
//      if ( blockbegin < dataTmin && blockbegin+blockdur > dataTmin ) {
//        blockdur -= dataTmin-blockbegin;
//        blockbegin = dataTmin;
//      }
//      if ( blockbegin + blockdur > dataTmax ) { blockdur = dataTmax - blockbegin; }
//      //a bit of housekeeping above to keep things within [dataTmin, dataTmax] so the norm is more useful
//
//      bTemp.push_back(BlockLevel(blockbegin, blocklev, blockdur));
//      norm1 = norm1 + (blocklev * blockdur);
//      lTtemp = lTtemp + blockdur;
//    }
//
//    levels_ = bTemp;
//    livetimeTotal_ = lTtemp;
//    norm_ = norm1;
//  }
//
//};



class BlockTimePdf1 : public TimePdf {
  private:
    vector<BlockLevel> levels_;
    string infile;
    double threshMax_;

 public:
  BlockTimePdf1() {
    norm_ = 0.;
    livetimeTotal_ = 0.;
  }

  BlockTimePdf1(BlockLevel a){
    levels_.push_back(a);
    norm_ = norm_ + (a.GetBlockLevel() * a.GetBlockDuration());
    livetimeTotal_ = livetimeTotal_ + a.GetBlockDuration();
  }

  TimePdf* Clone() const { return new BlockTimePdf1(*this); }

  void SetBlockLevels(string fname, double cut=0., double offset=0.) {
    infile = fname;
    norm_ = 0;
    levels_.clear(); // resetting the blocks each time we call this
    livetimeTotal_ = 0.;
    ifstream fin;
    fin.open(fname.c_str());
    double blockbegin, blockdur, blocklev;
    double blast=0; double bdlast=0;
    double tMax = -100.;
    while (fin >> blockbegin) {
      fin >> blocklev >> blockdur;
      if (blocklev > tMax) { tMax = blocklev; }
      blockbegin += offset;
      if (blocklev < cut) {
        levels_.push_back(BlockLevel(blockbegin, 0., blockdur));
      } else {
        levels_.push_back(BlockLevel(blockbegin, blocklev-cut, blockdur));
        norm_ = norm_ + ( (blocklev-cut) * blockdur);
      }
      livetimeTotal_ = livetimeTotal_ + blockdur;
      if (blast && blockbegin-(blast+bdlast)>1e-4) {cout << blockbegin << " " << blast << endl;}
      blast = blockbegin;
      bdlast = blockdur;
    }
    tmin_ = levels_[0].GetBlockBegin();
    tmax_ = tmin_+livetimeTotal_;

    threshMax_ = tMax;
    fin.close();
  }

  void SetNorm() { cout << "Norm is set when we load the block file, \n so we're not doing anything here." << endl; }

  double GetPdfValue(double MJD){
    if (MJD<tmin_ ||MJD>tmax_) return 0.;
    if ( norm_<1e-20 ) return 0.; //saying no block is above the threshold
    for (unsigned int i=0; i <= levels_.size(); i++) {
      if (levels_[i].BBegin <= MJD && (levels_[i].BBegin + levels_[i].BDuration) > MJD) {
        return levels_[i].BLevel/norm_;
      }
    }
    return 0.;
  }

//  double GetPdfValue(Time t)             { return GetPdfValue( t.GetMJD() ); }
//  double GetPdfValue(Time t, double lag) { return GetPdfValue( t.GetMJD()-lag ); }

  void fillHisto(int nbins=0){
    // this can be different for the block Pdf using a variable binned histogram
    // which is basically what the blocks are anyway

//     int nbins_ = levels_.size();
//     cout << "nbins " << nbins_ << endl;
//     double begin[nbins_+1];
//     for (int i=0; i<nbins_; i++){
//       begin[i] = levels_[i].BBegin;
//     }
//     begin[nbins_] = tmax_;
//     histo_ = new TH1D("histo_","source",nbins_,begin);
//     double MJD;
//     for (int i=0; i < nbins_; i++){ // bins go from 1 to n
//       MJD = levels_[i].BBegin+levels_[i].BDuration/2.;
//       histo_->SetBinContent(i+1,GetPdfValue(MJD));
//       cout << "MJD: " << MJD << " getpdfvalue " << GetPdfValue(MJD) << endl;
//     }

         if(nbins==0) nbins=10000;
         double tmin = levels_[0].BBegin;
         double tmax = tmin + livetimeTotal_;
         cout << "livetimeTotal " << livetimeTotal_ << endl;

         double bintime = (tmax - tmin) / nbins;
         double MJD;
         histo_ = new TH1D("histo","source",nbins,tmin,tmax);
         for (int i=1; i <= nbins; i++){ // bins go from 1 to n
         MJD = tmin + (i-0.5)*bintime;
         histo_->SetBinContent(i,GetPdfValue(MJD));
       }

  }

  TH1D * GetHisto() {return histo_;}


  Time GenerateEventTime(){
    //gRandom->SetSeed(0);
    //gRandom->SetSeed(get_ran1_seed());
    double t = histo_->GetRandom();
    return Time(t);
  }

  void CheckTimeBounds(double dataTmin, double dataTmax) {
    // a bit of housekeeping above to keep things within [dataTmin, dataTmax]
    // so the norm is more useful for use with multiple datasets.

    // in this function it won't matter if I made something with a cut
    // or offset or whatever to begin with, since we're dealing with
    // the blocks themselves.

    double norm1 = 0.;
    double lTtemp = 0;
    vector<BlockLevel> bTemp;
    double blockbegin = 0.;
    double blocklev   = 0.;
    double blockdur   = 0.;
    for (unsigned int i=0; i<levels_.size(); i++) {
      blockbegin = levels_[i].BBegin;
      blockdur = levels_[i].BDuration;
      blocklev = levels_[i].BLevel;
      if ( blockbegin+blockdur <= dataTmin || blockbegin > dataTmax ) { continue; }
      if ( blockbegin < dataTmin && blockbegin+blockdur > dataTmin ) {
        blockdur -= dataTmin-blockbegin;
        blockbegin = dataTmin;
      }
      if ( blockbegin + blockdur > dataTmax ) {
        blockdur = dataTmax - blockbegin;
      }

      bTemp.push_back(BlockLevel(blockbegin, blocklev, blockdur));
      norm1 = norm1 + (blocklev * blockdur);
      lTtemp = lTtemp + blockdur;
    }

    levels_ = bTemp;
    livetimeTotal_ = lTtemp;
    norm_ = norm1;
  }

};

// I need a few extra functions for working with blocks
// First is a function which tells me the time above threshold
// given a text file with the blocks and the threshold.

double BlockTimeAboveThresh(string fname, double thresh, double dataTmin=0, double dataTmax=1e6) {


  ifstream fin;
  fin.open(fname.c_str());
  double blockbegin, blockdur, blocklev;
  double t=0;

  while (fin >> blockbegin) {
    fin >> blocklev >> blockdur;

    if(blocklev > thresh) {
      if ( blockbegin < dataTmin && blockbegin+blockdur > dataTmin ) {
        blockdur = blockbegin + blockdur - dataTmin;
      }
      if ( blockbegin + blockdur > dataTmax ) {
        blockdur = dataTmax - blockbegin;
      }
      if ( blockbegin+blockdur <= dataTmin || blockbegin > dataTmax ) {
        blockdur=0.;
      }

      t += blockdur;
    }
  }
  fin.close();

  return t;
}


double BlockNormForThresh(string fname, double thresh, double tmin, double tmax) {

  TimePdf * ttt = new BlockTimePdf1();
  ttt->SetBlockLevels(fname,thresh);
  ttt->CheckTimeBounds(tmin,tmax);

  double toReturn = ttt->GetNorm();
  delete ttt;

  return toReturn;

}

// Two quick functions for getting the highest and second highest
// flux level in a MLB series. I noticed that sometimes the analysis
// will try to go to the highest block level and squirm around there
// before returning a bad test statistic, now what I do is to set
// the maximum fit value halfway between the highest and second
// highest value from a block function.

double GetHighestBlock(string fname) {
  ifstream fin;
  fin.open(fname.c_str());
  double blockbegin, blockdur, blocklev;
  double max=-1;
  while (fin >> blockbegin) {
    fin >> blocklev >> blockdur;
    if (blocklev > max) {
      max = blocklev;
    }
  }
  fin.close();
  return max;
}

double GetSecondHighestBlock(string fname) {
  ifstream fin;
  fin.open(fname.c_str());

  double blockbegin, blockdur, blocklev;
  double max=-1;
  double secondmax=-2;

  int NUMlines = 0;

  while (fin >> blockbegin) {
    fin >> blocklev >> blockdur;
    if (blocklev > secondmax) {
      if (blocklev > max) {
        secondmax = max;
        max = blocklev;
      } else {
        secondmax = blocklev;
      }
    }

    //test in order to treat constant lightcurves
    if(secondmax==max){
        cout << "Constant lightcurve: setting max threshold to 0.9*highest block" << endl;
        secondmax = 0.9*max;
    }
   NUMlines++;
  
    }

  
  if(NUMlines==1){
     cout << "File has only one block: setting max threshold to 0.9*highest block." << endl;
     secondmax = 0.9*GetHighestBlock(fname.c_str());

  }

  fin.close();
  return secondmax;
}


double SearchForLag(MultiAnalysisFn maf, double laglimit) {

  //maf.PrepareAnalysis();

  double llhMax=-100.;
  double lagb=0., llhtemp;
  
  double pardef[4];
  
  pardef[0] = 2.;
  pardef[1] = 2.;
  pardef[3] = 0.;

  for (double d=-1.0*laglimit; d<laglimit; d=d+laglimit/10.){
    pardef[2] = d;
    llhtemp = maf.EvaluateLlh( pardef );
    if (llhtemp > llhMax) {
      llhMax = llhtemp;
      lagb = d;
    }
  }  
  
  return lagb;
}


void SearchBlockSpace(MultiAnalysisFn maf, string BlocksTimeFile, double lag, double & initV, double & maxT) {

  maxT = ( GetHighestBlock(BlocksTimeFile) + GetSecondHighestBlock(BlocksTimeFile) )/2.;

  double step = maxT/20.;
  double llhMax=-100.;
  double llhtemp;
  double pardef[4];
  
  pardef[0] = 2.;
  pardef[1] = 2.;
  pardef[2] = lag;
  
  for (double d=maxT/20.; d<maxT; d+=step) {
    pardef[3] = d;    
    llhtemp = maf.EvaluateLlh( pardef );
    if (llhtemp > llhMax) {
      llhMax = llhtemp;
      initV = d;
    }
  }    
}


#endif // LLHTIMEDEP_BLOCKLEVEL_H_
