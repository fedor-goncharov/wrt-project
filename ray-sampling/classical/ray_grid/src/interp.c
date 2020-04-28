#include <math.h>

inline double cell_linear_interp(double x, double v1, double v2) __attribute__((always_inline));
double cell_linear_interp(double x, double v1, double v2) {
  return ((1-x)*v1 + x*v2); 
}

inline double cell_bilinear_interp(double x, double y, double v1, double v2, double v3, double v4) __attribute__((always_inline));
double cell_bilinear_interp(double x, double y, double v1, double v2, double v3, double v4) {
  double c1 = cell_linear_interp(x, v1, v2),
	 c2 = cell_linear_interp(x, v4, v3);
	 
  return cell_linear_interp(y, c1, c2);
}


double square_bilinear_interp(double ***values, int npixels, int nslices, double radius, 
    int i_slice, double x, double y) {
  
  if ((fabs(x) > radius) || (fabs(y) > radius))
    return 0;
  
  //get cells numbers
  const double delta = 2.0*radius / (npixels - 1);
  int ix = (int)((x + radius) / delta), 
      iy = (int)((y + radius) / delta);
  
  double cell_x = fmod((x + radius), delta) / delta,
	 cell_y = fmod((y + radius), delta) / delta;
	 
  //control of positive boundary points
  if (ix == (npixels - 1)) {
     ix -= 1;
     cell_x = 1.0;
  }
  if (iy == (npixels - 1)) {
     iy -= 1;
     cell_y = 1.0;
  }
 
  //trilinear interpolation at (x,y,z)    
  double interp_value = cell_bilinear_interp(cell_x, cell_y, 
			values[i_slice][ix][iy], values[i_slice][ix + 1][iy], values[i_slice][ix + 1][iy + 1], values[i_slice][ix][iy + 1]);
  return interp_value;
}

//no interpolation of a test-function
double square_zero_interp(double ***values, int npixels, int nslices, double radius, 
  int i_slice, double x, double y) {
  if ((fabs(x) > radius) || (fabs(y) > radius))
    return 0;
  
  //get cells numbers
  const double delta = 2.0*radius / (npixels - 1);
  int ix = (int)((x + radius) / delta), 
      iy = (int)((y + radius) / delta);
  
  //control of positive boundary points
  if (ix == npixels) {
     ix -= 1;
  }
  if (iy == npixels) {
     iy -= 1;
  }
  return values[i_slice][ix][iy];
}
