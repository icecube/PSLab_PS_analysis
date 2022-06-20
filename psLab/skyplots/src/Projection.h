#ifndef PROJECTION_H_
#define PROJECTION_H_


class Projection {

 public:
  Projection() :
    lonLatInDeg_(false),
    shiftLonDeg_(0.),
    reverseLon_(false),
    xMin_(0.),
    xMax_(0.),
    yMin_(0.),
    yMax_(0.)
  { 
    ResetLast(); 
  }

  virtual ~Projection() { }

  void SetLonLatInDeg(bool lonLatInDeg) {
    lonLatInDeg_ = lonLatInDeg; 
    ResetLast();
  }

  void SetShiftLonDeg(double shiftDeg) {
    shiftLonDeg_ = shiftDeg;
    ResetLast();
  }

  void SetReverseLon(bool reverseLon) { 
    reverseLon_ = reverseLon;
    ResetLast();
  }


  bool Project(double lon, double lat) {
    if ( lastProject_ && lon == lon_ && lat == lat_) {  // no change
      return valid_; 
    }
    lon_ = lon;
    lat_ = lat;
    project_fn();
    lastProject_ = true;
    lastProjectInverse_ = false;
    return valid_;
  }
 
  bool ProjectInverse(double x, double y) {
    if ( lastProjectInverse_ && x == x_ && y == y_) {   // no change
      return valid_; 
    }
    x_ = x;
    y_ = y;
    project_inverse_fn();
    lastProject_ = false;
    lastProjectInverse_ = true;
    return valid_;
  }

  bool IsValid() const { return valid_; }

  double X() const { return x_; }
  double Y() const { return y_; }
  double Lon() const { return lon_; }
  double Lat() const { return lat_; }

  double X(double lon, double lat) { 
    Project(lon,lat);
    return x_;
  }

  double Y(double lon, double lat) { 
    Project(lon,lat);
    return y_;
  }

  double Lon(double x, double y) {
    ProjectInverse(x,y);
    return lon_;
  }

  double Lat(double x, double y) {
    ProjectInverse(x,y);
    return lat_;
  }


  double GetXMin() const { return xMin_; }
  double GetXMax() const { return xMax_; }
  double GetYMin() const { return yMin_; }
  double GetYMax() const { return yMax_; }


 protected:
  // General projection parameters
  bool lonLatInDeg_;
  double shiftLonDeg_;
  bool reverseLon_;

  // Keep track of whether last operation was Project or ProjectInverse
  bool lastProject_;
  bool lastProjectInverse_;

  // Keep track of last values used:
  double lon_;
  double lat_;
  double x_;
  double y_;

  // Valid_ can provide a check if transformation was in range
  bool valid_;

  // Projection Range
  double xMin_;
  double xMax_;
  double yMin_;
  double yMax_;

  void ResetLast() {
    lastProject_ = false;
    lastProjectInverse_ = false;
    lon_ = 0.; 
    lat_ = 0.; 
    x_ = 0.;
    y_ = 0.;
    valid_ = false;
  }

  virtual void
    project_fn()=0;
  virtual void 
    project_inverse_fn()=0;
};


#endif // PROJECTION_H_
