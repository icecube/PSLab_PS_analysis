#ifndef TRANSFORMATIONFNS_H_
#define TRANSFORMATIONFNS_H_


#include "TMath.h"


double Mod_MinusPi_To_Pi(const double ang) {
  double newAng = fmod(ang, TMath::TwoPi());
  if ( newAng >  TMath::Pi() ) { return newAng - TMath::TwoPi(); }
  if ( newAng < -TMath::Pi() ) { return newAng + TMath::TwoPi(); }
  return newAng;
}

double Mod_Zero_To_TwoPi(const double ang) {
  double newAng = fmod(ang, TMath::TwoPi());
  if ( newAng < 0. ) { return newAng + TMath::TwoPi(); }
  return newAng;
}




//
//  HAMMER-AITOFF
//

// ---------------------------------------------------------------------------
/// \brief Hammer-Aitoff projection: lon,lat --> x,y
// ---------------------------------------------------------------------------
//  lon = [-Pi, +Pi]  ==>  x = [-sqrt(2), +sqrt(2)]
void hammer_aitoff_project(const double lon, const double lat, 
			   double *x, double *y)
{
  double t = sqrt(2./ (1+cos(lat)*cos(lon/2.)));
  *x = 2*t * cos(lat)*sin(lon/2.);
  *y = t * sin(lat);
}

// ---------------------------------------------------------------------------
/// \brief Inverse Hammer-Aitoff projection: x,y --> lon,lat
// ---------------------------------------------------------------------------
// x = [-sqrt(2), +sqrt(2)]  ==>   lon = [-Pi, +Pi]
void hammer_aitoff_inverse_project(const double x, const double y, 
				   double *lon, double *lat)
{
  double z = sqrt(1-(x/4)*(x/4)-(y/2)*(y/2));
  *lon = 2 * atan2(z*x , (2*(2*z*z-1)) );
  *lat = asin(y*z);
}



//
// MOLLWEIDE
//



// Functions associated with map projections...
// ---------------------------------------------------------------------------
/// \brief Transcendental equation for intermediate Mollweide parameter gamma.
// ---------------------------------------------------------------------------
double mollweide_f(double g, double th)
{
   return (2.0 * g + sin(2.0 * g) - TMath::Pi() * sin(th));
}

// ---------------------------------------------------------------------------
/// \brief First derivative of eq. for intermediate Mollweide parameter gamma.
// ---------------------------------------------------------------------------
double mollweide_df(double g, double th)
{
  th = th;  // touch variable, so no warnings during compile
  return (2.0 + 2.0 * cos(2.0 * g));
}



// ---------------------------------------------------------------------------
/// \brief Newton-Raphson root finding algorithm used to solve for gamma.
// ---------------------------------------------------------------------------
double mollweide_root_find(
   double (*f)(double, double),
   double (*df)(double, double),
   double a,
   double b,
   double th
   )
{
   static double precision = 1.0e-4;
   double dx = b - a;
   double x0 = (a + b) / 2.0;
   double x1;

   do {
      x1 = x0 - f(x0, th) / df(x0, th);
      dx = x1 - x0;
      x0 = x1;
   } while (fabs(dx) > precision);

   return x0;
}



// ---------------------------------------------------------------------------
/// \brief Mollweide projection: lon,lat --> x,y
// ---------------------------------------------------------------------------
void mollweide_project(const double lon, const double lat, 
		       double *x, double *y)
{
  double g = mollweide_root_find(mollweide_f, mollweide_df, 
				 -0.5*TMath::Pi(), 0.5*TMath::Pi(), lat);
  *x = 2. * TMath::Sqrt(2.) * lon * TMath::Cos(g) / TMath::Pi();
  *y = TMath::Sqrt(2.) * TMath::Sin(g);
}

// ---------------------------------------------------------------------------
/// \brief Inverse Mollweide projection: x,y --> lon,lat
// ---------------------------------------------------------------------------
void
mollweide_inverse_project(const double x, const double y, 
			  double *lon, double *lat)
{
  double g = TMath::ASin(y / TMath::Sqrt(2.));
  *lon = TMath::Pi() * x / (2.*TMath::Sqrt(2.) * TMath::Cos(g));
  *lat = TMath::ASin((2*g + TMath::Sin(2*g)) / TMath::Pi());
}



#endif //TRANSFORMATIONFNS_H_
