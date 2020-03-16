#ifndef _H_INTERPOLATION
#define _H_INTERPOLATION

  double cube_trilinear_interp(double*** values, int ngrid,
    double x, double y, double z);

  double zero_interp(double*** values, int ngrid, 
  	double x, double y, double z);
  
#endif