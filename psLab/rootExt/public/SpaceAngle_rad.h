#ifndef ROOTEXT_SPACEANGLE_RAD_H_
#define ROOTEXT_SPACEANGLE_RAD_H_


double SpaceAngle_rad(double zenith1_rad, double azimuth1_rad, 
		      double zenith2_rad, double azimuth2_rad)
{
  double dotprod = sin(zenith1_rad) * sin(zenith2_rad) *
                   cos(azimuth1_rad-azimuth2_rad) +
                   cos(zenith1_rad) * cos(zenith2_rad);
  return acos(dotprod);
}

#endif // ROOTEXT_SPACEANGLE_RAD_H_
