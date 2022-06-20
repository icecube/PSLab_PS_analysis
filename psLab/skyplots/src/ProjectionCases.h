#ifndef PROJECTIONCASES_H_
#define PROJECTIONCASES_H_

#include "Projection.h"


class ChainProjection : public Projection {
 public:
  ChainProjection() { }
  virtual ~ChainProjection() { }
  void SetFirstProjection(Projection* firstPro) { firstPro_ = firstPro; } 
  void SetSecondProjection(Projection* secondPro) { 
    secondPro_ = secondPro; 
    xMin_ = secondPro->GetXMin();
    xMax_ = secondPro->GetXMax();
    yMin_ = secondPro->GetYMin();
    yMax_ = secondPro->GetYMax();
  }

 private:
  Projection* firstPro_;
  Projection* secondPro_;
  virtual void project_fn() {
    valid_ = firstPro_->Project(lon_ , lat_);
    if (valid_) {
      valid_ = secondPro_->Project( firstPro_->X() , firstPro_->Y() );
      if (valid_) {
	x_ = secondPro_->X();
	y_ = secondPro_->Y();
      }
    }
    return;
  }
  virtual void project_inverse_fn() {
    valid_ = secondPro_->ProjectInverse(x_ , y_ );
    if (valid_) {
      valid_ = firstPro_->ProjectInverse(secondPro_->Lon(),secondPro_->Lat() );
      if (valid_) {
	lon_ = firstPro_->Lon();
	lat_ = firstPro_->Lat();
      }
    }
  }
};
  


class HammerAitoffProjection : public Projection {
 public:
  HammerAitoffProjection();
  virtual ~HammerAitoffProjection() { }
 private:
  virtual void project_fn();
  virtual void project_inverse_fn();
};


class FlatProjection : public Projection {
 public:
  FlatProjection();
  virtual ~FlatProjection() { }
 private:
  virtual void project_fn();
  virtual void project_inverse_fn();

};


class EquatorialToGalacticProjection : public Projection {
 public:
  EquatorialToGalacticProjection();
  virtual ~EquatorialToGalacticProjection() { }
 private:
  virtual void project_fn();
  virtual void project_inverse_fn();
};

class GalacticToEquatorialProjection : public Projection {
 public:
  GalacticToEquatorialProjection();
  virtual ~GalacticToEquatorialProjection() { }
 private:
  virtual void project_fn();
  virtual void project_inverse_fn();
};


#endif // PROJECTIONCASES_H_

