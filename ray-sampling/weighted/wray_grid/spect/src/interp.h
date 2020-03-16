#ifndef _H_RAY_GRID_INTERPOLATION
#define _H_RAY_GRID_INTERPOLATION

  double cube_trilinear_interp(double*** values, 
			       int ngrid_tfunc, 
			       double x, 
			       double y, 
			       double z);
  
#endif