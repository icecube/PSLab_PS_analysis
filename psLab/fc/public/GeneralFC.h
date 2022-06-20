#include <vector>

#include "TH1D.h"
#include "TH2D.h"



class GeneralFC {

 public:

  GeneralFC() {
    nSigLoadTotal_ = 0;
    sigMax_ = 0.;
    optPrint_ = false;

    // somewhat arbitrary defaults
    SetParams(1100, -10., 100., 510);
    // Originally, these were much larger
    //#define NBINSX 3000
    //#define NBINSY 5000
  }

  // default sigUnc = 0., meaning no signal uncertainty (usual FC)
  void SetParams(int nBinsX, double lowLimX, double upLimX, int nBinsY, 
		 double sigUnc = 0.) {
    nBinsX_ = nBinsX;
    lowLimX_ = lowLimX;
    upLimX_ = upLimX;
    nBinsY_ = nBinsY;
    sigUnc_ = sigUnc;

    // must clear these if parameters get changed:
    rawData_.clear();
    matrix_.clear();
    maxes_.clear();
    ill_.clear();
    iul_.clear();    
  }


  virtual ~GeneralFC() { }

  void Load(char *dir, unsigned nSigLoadTotal, double sigMax);

  void BuildConfidenceBand(double cl);


  double GetUpperLimit(double testStat);
  double GetLowerLimit(double testStat);

  double GetMedianUpperLimit();
  double GetMedianLowerLimit();

  double GetAverageUpperLimit();
  double GetAverageLowerLimit();

  double GetMedianOfNullHypothesis();


  TH2D* GetHistRawData();
  TH2D* GetHistMatrix();
  TH1D* GetHistLowerLimit();
  TH1D* GetHistUpperLimit();


 private:

  vector<vector<double> > rawData_;
  vector<vector<double> > matrix_;
  vector<double> maxes_;
  vector<double> ill_;
  vector<double> iul_;


  int nSigLoadTotal_;
  double sigMax_;
  bool optPrint_;
  double sigUnc_;

  int nBinsX_;
  double lowLimX_;
  double upLimX_;
  int nBinsY_;


  virtual void LoadRawData(char *dir);
  virtual bool LoadTestStat(FILE *f, double& testStat);
  // can replace these load functions to accomodate any kind of data

  double Poisson(int k, double mean);

  int GetBinXFromTestStat(double testStat);
  double GetTestStatFromBinX(int binX);

  double GetSigFromBinY(int binY);
};



// Takes first number on line as val, then testStat = sign(val)*sqrt(2*fabs(val))

class SignedSqrtFC : public GeneralFC {
  virtual bool LoadTestStat(FILE *f, double& testStat);
};
