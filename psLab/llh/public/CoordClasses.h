#ifndef LLH_COORDCLASSES_H_
#define LLH_COORDCLASSES_H_


// Abstract Base Class

class Coord {
public:
  virtual ~Coord() { }
  virtual double DistanceTo (const Coord& coord2) const = 0;
};


double DistanceBetween (const Coord& c1, const Coord& c2) 
{
  return c1.DistanceTo(c2);
}


double CircularGaussUnc(double r, double sigma);


// CARTESIAN COORDINATES //


class Cartesian : public Coord {
  double x_;
  double y_;
public:
  Cartesian () { x_=0.; y_=0.; }
  Cartesian (double x, double y) { x_=x; y_=y; }
  void SetCoords (double x, double y) { x_=x; y_=y; }
  double GetX() const {return x_; }
  double GetY() const {return y_; }

  double DistanceTo (const Coord& coord2) const
  {
    const Cartesian* cart2 = dynamic_cast<const Cartesian*>(&coord2);
    assert(cart2);
    double dx=cart2->GetX() - x_; 
    double dy=cart2->GetY() - y_;
    return sqrt(dx*dx+dy*dy);
  }
};




#endif // LLH_COORDCLASSES_H_
