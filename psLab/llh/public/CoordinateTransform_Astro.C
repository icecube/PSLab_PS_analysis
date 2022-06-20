#include "llh/public/CoordinateTransform.h"

#include "astro_interface/public/AstroHeader.h"


void LocalToEq(double zenDeg, double azDeg, double timeMJD,
	       double& raDeg, double& decDeg)
{
  static Astro::IceCubeDetector ice;
  Astro::LocalCoord loc(zenDeg, azDeg, Astro::Deg);
  Astro::Time t(timeMJD);
  // Astro::EqCoord eq = ice.LocalToEquatorial(loc, t);
  Astro::EqCoord eq = ice.LocalToEquatorialJ2000(loc, t);
  raDeg = eq.GetRaDeg();
  decDeg = eq.GetDecDeg();
}

void EqToLocal(double raDeg, double decDeg, double timeMJD,
	       double& zenDeg, double& azDeg)
{
  static Astro::IceCubeDetector ice;
  Astro::EqCoord eq(raDeg, decDeg, Astro::Deg);
  Astro::Time t(timeMJD);
  // Astro::LocalCoord loc = ice.EquatorialToLocal(eq, t);
  Astro::LocalCoord loc = ice.EquatorialJ2000ToLocal(eq, t);
  zenDeg = loc.GetZenithDeg();
  azDeg = loc.GetAzimuthDeg();
}
