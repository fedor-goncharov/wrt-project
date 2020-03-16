#ifndef __GRID_INTERPOLATION_H
#define __GRID_INTERPOLATION_H

  double cube_trilinear_interp(double*** values, int ngrid,
    double x, double y, double z);

  double cube_zero_interp(double*** values, int ngrid,
    double x, double y, double z);

#endif
