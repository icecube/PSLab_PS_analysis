#include "ProjectionCases.h"

#include "TMath.h"

#include "TransformationFns.h"

#include "SLAHeader_Galactic.h"


//
// HAMMER-AITOFF
//

HammerAitoffProjection::HammerAitoffProjection() {
  SetLonLatInDeg(true);
  SetShiftLonDeg(-180.);
  SetReverseLon(true);
  xMin_ = - 2 * sqrt(2);
  xMax_ =   2 * sqrt(2);
  yMin_ = - sqrt(2);
  yMax_ =   sqrt(2);
}


void HammerAitoffProjection::project_fn() {
  double lon = lon_;
  double lat = lat_;

  // prepare lon,lat based on ProjectionParameters
  if (lonLatInDeg_) {
    lon *= TMath::DegToRad();
    lat *= TMath::DegToRad();
  }
  lon += shiftLonDeg_ * TMath::DegToRad();
  if (reverseLon_) { lon = -lon; }
  
  lon = Mod_MinusPi_To_Pi(lon);
  if (lat >  TMath::PiOver2()) { lat =  TMath::PiOver2(); }
  if (lat < -TMath::PiOver2()) { lat = -TMath::PiOver2(); }
  
  hammer_aitoff_project(lon, lat, &x_, &y_);
  
  // restrictions applied to lon,lat above mean the result is always valid 
  valid_ = true;
}
    

void HammerAitoffProjection::project_inverse_fn() {
  hammer_aitoff_inverse_project(x_, y_, &lon_, &lat_);
  valid_ = true;
  if (lon_ < -TMath::Pi() || lon_ > TMath::Pi() ) { 
    valid_ = false; 
  }
  // prepare lon,lat based on ProjectionParameters
  if (reverseLon_) { lon_ = -lon_; }
  lon_ -= shiftLonDeg_ * TMath::DegToRad();
  lon_ = Mod_Zero_To_TwoPi(lon_);

  if (lonLatInDeg_) {
    lon_ *= TMath::RadToDeg();
    lat_ *= TMath::RadToDeg();
  }
}


//
// FLAT PROJECTION
//

FlatProjection::FlatProjection() {
  SetLonLatInDeg(true);
  SetShiftLonDeg(-180.);
  SetReverseLon(true);
  // This is somewhat arbitrary
  xMin_ = -1.;
  xMax_ =  1.;
  yMin_ = -1.;
  yMax_ =  1.;
}


void FlatProjection::project_fn() {
  double lon = lon_;
  double lat = lat_;

  // prepare lon,lat based on ProjectionParameters
  if (lonLatInDeg_) {
    lon *= TMath::DegToRad();
    lat *= TMath::DegToRad();
  }
  lon += shiftLonDeg_ * TMath::DegToRad();
  if (reverseLon_) { lon = -lon; }

  lon = Mod_MinusPi_To_Pi(lon);
  if (lat >  TMath::PiOver2()) { lat =  TMath::PiOver2(); }
  if (lat < -TMath::PiOver2()) { lat = -TMath::PiOver2(); }

  x_ = lon / TMath::Pi();
  y_ = lat / TMath::PiOver2();

  // restrictions applied to lon,lat above mean the result is always valid 
  valid_ = true;
}
    

void FlatProjection::project_inverse_fn() {
  lon_ = x_ * TMath::Pi();
  lat_ = y_ * TMath::PiOver2();
  valid_ = true;
  if (lon_ < -TMath::Pi() || lon_ > TMath::Pi() ||
      lat_ < -TMath::PiOver2() || lat_ > TMath::PiOver2() ) { 
    valid_ = false; 
  }

  // prepare lon,lat based on ProjectionParameters
  if (reverseLon_) { lon_ = -lon_; }
  lon_ -= shiftLonDeg_ * TMath::DegToRad();
  lon_ = Mod_Zero_To_TwoPi(lon_);
  if (lonLatInDeg_) {
    lon_ *= TMath::RadToDeg();
    lat_ *= TMath::RadToDeg();
  }
}



//
// EQUATORIAL TO GALACTIC PROJECTION
//

EquatorialToGalacticProjection::EquatorialToGalacticProjection() {
  SetLonLatInDeg(true);
  SetShiftLonDeg(0.);     // NOT VALID FOR COORDINATE TRANSFORMATION
  SetReverseLon(false);   // NOT VALID FOR COORDINATE TRANSFORMATION
  // This is somewhat arbitrary
  xMin_ = 0.;
  xMax_ = 360.;
  yMin_ = -90.;
  yMax_ =  90.;
}


void EquatorialToGalacticProjection::project_fn() {
  double eqLon = lon_;
  double eqLat = lat_;

  // prepare lon,lat based on ProjectionParameters
  if (lonLatInDeg_) {
    eqLon *= TMath::DegToRad();
    eqLat *= TMath::DegToRad();
  }
  // NOT VALID FOR COORDINATE TRANSFORMATION
  //  lon += shiftLonDeg_ * TMath::DegToRad();
  // NOT VALID FOR COORDINATE TRANSFORMATION
  //  if (reverseLon_) { lon = -lon; }

  eqLon = Mod_Zero_To_TwoPi(eqLon);
  if (eqLat >  TMath::PiOver2()) { eqLat =  TMath::PiOver2(); }
  if (eqLat < -TMath::PiOver2()) { eqLat = -TMath::PiOver2(); }

  const double mjd_J2000 = 51544.5;

  double galLon = SLACoordinateTransform::Equatorial2GalacticLongitude(
    eqLon, eqLat, mjd_J2000);
  double galLat = SLACoordinateTransform::Equatorial2GalacticLatitude(
    eqLon, eqLat, mjd_J2000);

  x_ = galLon * TMath::RadToDeg();
  y_ = galLat * TMath::RadToDeg();
  valid_ = true;
}


void EquatorialToGalacticProjection::project_inverse_fn() {
  double galLon = x_ * TMath::DegToRad();
  double galLat = y_ * TMath::DegToRad();
  valid_ = true;

  galLon = Mod_Zero_To_TwoPi(galLon);
  if (galLat >  TMath::PiOver2()) { galLat =  TMath::PiOver2(); }
  if (galLat < -TMath::PiOver2()) { galLat = -TMath::PiOver2(); }

  const double mjd_J2000 = 51544.5;

  lon_ = SLACoordinateTransform::Galactic2RA(
    galLon, galLat, mjd_J2000);
  lat_ = SLACoordinateTransform::Galactic2Dec(
    galLon, galLat, mjd_J2000);

  // prepare lon,lat based on ProjectionParameters
  // NOT VALID FOR COORDINATE TRANSFORMATION
  //  if (reverseLon_) { lon_ = -lon_; }
  // NOT VALID FOR COORDINATE TRANSFORMATION
  //  lon_ -= shiftLonDeg_ * TMath::DegToRad();
  // NOT NEEDED
  //  lon_ = Mod_Zero_To_TwoPi(lon_);
  if (lonLatInDeg_) {
    lon_ *= TMath::RadToDeg();
    lat_ *= TMath::RadToDeg();
  }
}


//
// GALACTIC TO EQUATORIAL PROJECTION
//

GalacticToEquatorialProjection::GalacticToEquatorialProjection() {
  SetLonLatInDeg(true);
//  SetShiftLonDeg(0.);     // NOT VALID FOR COORDINATE TRANSFORMATION
//  SetReverseLon(false);   // NOT VALID FOR COORDINATE TRANSFORMATION
  // This is somewhat arbitrary
  xMin_ = 0.;
  xMax_ = 360.;
  yMin_ = -90.;
  yMax_ =  90.;
}

void GalacticToEquatorialProjection::project_fn() {
  double galLon = lon_ * TMath::DegToRad();
  double galLat = lat_ * TMath::DegToRad();
  valid_ = true;

  galLon = Mod_Zero_To_TwoPi(galLon);
  if (galLat >  TMath::PiOver2()) { galLat =  TMath::PiOver2(); }
  if (galLat < -TMath::PiOver2()) { galLat = -TMath::PiOver2(); }

  const double mjd_J2000 = 51544.5;

  x_ = SLACoordinateTransform::Galactic2RA(
    galLon, galLat, mjd_J2000);
  y_ = SLACoordinateTransform::Galactic2Dec(
    galLon, galLat, mjd_J2000);

  // prepare lon,lat based on ProjectionParameters
  // NOT VALID FOR COORDINATE TRANSFORMATION
  //  if (reverseLon_) { lon_ = -lon_; }
  // NOT VALID FOR COORDINATE TRANSFORMATION
  //  lon_ -= shiftLonDeg_ * TMath::DegToRad();
  // NOT NEEDED
  //  lon_ = Mod_Zero_To_TwoPi(lon_);
  if (lonLatInDeg_) {
    x_ *= TMath::RadToDeg();
    y_ *= TMath::RadToDeg();
  }
}


void GalacticToEquatorialProjection::project_inverse_fn() {
  double eqLon = x_;
  double eqLat = y_;

  // prepare lon,lat based on ProjectionParameters
  if (lonLatInDeg_) {
    eqLon *= TMath::DegToRad();
    eqLat *= TMath::DegToRad();
  }
  // NOT VALID FOR COORDINATE TRANSFORMATION
  //  lon += shiftLonDeg_ * TMath::DegToRad();
  // NOT VALID FOR COORDINATE TRANSFORMATION
  //  if (reverseLon_) { lon = -lon; }

  eqLon = Mod_Zero_To_TwoPi(eqLon);
  if (eqLat >  TMath::PiOver2()) { eqLat =  TMath::PiOver2(); }
  if (eqLat < -TMath::PiOver2()) { eqLat = -TMath::PiOver2(); }

  const double mjd_J2000 = 51544.5;

  double galLon = SLACoordinateTransform::Equatorial2GalacticLongitude(
    eqLon, eqLat, mjd_J2000);
  double galLat = SLACoordinateTransform::Equatorial2GalacticLatitude(
    eqLon, eqLat, mjd_J2000);

  lon_ = galLon * TMath::RadToDeg();
  lat_ = galLat * TMath::RadToDeg();
  valid_ = true;
}

