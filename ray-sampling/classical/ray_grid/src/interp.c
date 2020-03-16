#include <math.h>

double cell_linear_interp(double x, double v1, double v2) {
  return ((1-x)*v1 + x*v2); 
}
  
double cell_bilinear_interp(double x, double y, double v1, double v2, double v3, double v4) {
  double c1 = cell_linear_interp(x, v1, v2),
	 c2 = cell_linear_interp(x, v4, v3);
	 
  return cell_linear_interp(y, c1, c2);
}

double cell_trilinear_interp(double x, double y, double z, 
			     double v1, double v2, double v3, double v4,
			     double v5, double v6, double v7, double v8) {
  double c1 = cell_bilinear_interp(x, y, v1, v2, v3, v4), 
	 c2 = cell_bilinear_interp(x, y, v5, v6, v7, v8);
	 
  return cell_linear_interp(z, c1, c2);
}



double cube_trilinear_interp(double*** values, int ngrid,
    double x, double y, double z) {
  
  if ((fabs(x) > 1) || (fabs(y) > 1) || (fabs(z) > 1))
    return 0;
  
  //get cells numbers
  const double delta = 2.0 / (ngrid - 1);
  int ix = (int)((x + 1.0) / delta), 
      iy = (int)((y + 1.0) / delta), 
      iz = (int)((z + 1.0) / delta);
  
  double cell_x = fmod((x + 1.0), delta) / delta,
	 cell_y = fmod((y + 1.0), delta) / delta,
	 cell_z = fmod((z + 1.0), delta) / delta;
	 
  //control of positive boundary points
  if (ix == (ngrid - 1)) {
     ix -= 1;
     cell_x = 1.0;
  }
  if (iy == (ngrid - 1)) {
     iy -= 1;
     cell_y = 1.0;
  }
  if (iz == (ngrid - 1)) {
     iz -= 1;
     cell_z = 1.0;
  }
 
  //trilinear interpolation at (x,y,z)    
  double interp_value = cell_trilinear_interp(cell_x, cell_y, cell_z, 
			values[ix][iy][iz], values[ix + 1][iy][iz], values[ix + 1][iy + 1][iz], values[ix][iy + 1][iz],
			values[ix][iy][iz + 1], values[ix + 1][iy][iz + 1], values[ix + 1][iy + 1][iz + 1], values[ix][iy + 1][iz + 1]);
  return interp_value;
}

//no interpolation of a test-function
double cube_zero_interp(double*** values, int ngrid, 
  double x, double y, double z) {
  if ((fabs(x) > 1) || (fabs(y) > 1) || (fabs(z) > 1))
    return 0;
  
  //get cells numbers
  const double delta = 2.0 / ngrid - 1;
  int ix = (int)((x + 1.0) / delta), 
      iy = (int)((y + 1.0) / delta), 
      iz = (int)((z + 1.0) / delta);
  
  //control of positive boundary points
  if (ix == ngrid) {
     ix -= 1;
     //cell_x = 1.0;
  }
  if (iy == ngrid) {
     iy -= 1;
     //cell_y = 1.0;
  }
  if (iz == ngrid) {
     iz -= 1;
     //cell_z = 1.0;
  }
  return values[ix][iy][iz];
}