#include <math.h>
#include "interp.h"

double ray_transform(double ***function_values, int npixels, int nslices, double radius,  
		     int i_slice, double shift, double phi) {

  double normal_x = cos(phi),
         normal_y = sin(phi);

  double direction_x = cos(phi + M_PI/2), 
         direction_y = sin(phi + M_PI/2);

  //initial point of the ray
  double in_point_x = shift * normal_x, 
         in_point_y = shift * normal_y;

  const double delta = (2.0*radius / (npixels-1))/(2.0);	//half-step of the initial step of the cube grid
  const int nsize = 4 * npixels;	//number of points per slice of the plane

  double integral = 0;
  int i_v;

  double ray_point_x = in_point_x - delta * direction_x * (nsize / 2),
         ray_point_y = in_point_y - delta * direction_y * (nsize / 2);
  for (i_v = 0; i_v < nsize + 1; ++i_v) {
    
       ray_point_x += delta * direction_x;
       ray_point_y += delta * direction_y;

       double value = square_bilinear_interp(function_values, npixels, nslices, radius, i_slice, ray_point_x, ray_point_y);
       integral += value; //1D integration along direction
  }
  return integral * delta;
}
