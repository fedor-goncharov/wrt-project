#include <math.h>

typedef struct {
   double x;
   double y;
} vect2d;

/*
  Code for integration of function of one pixel along lines.
*/

double matrix_element(int ngrid, int npixels, int i_pixel, double delta_grid,
  	                       double phi, double shift, double frequency) {

  const double half_d_grid = 0.5 * delta_grid;

  vect2d normal;
    normal.x = cos(phi);
    normal.y = sin(phi);

  vect2d direction;
    direction.x = cos(phi + M_PI/2);
    direction.y = sin(phi + M_PI/2);

  //initial point of the ray
  vect2d in_point;
    in_point.x = shift * normal.x;
    in_point.y = shift * normal.y;

  int i_pixel_x = i_pixel / (ngrid-1);
  int i_pixel_y = i_pixel % (ngrid-1);

  //coordinates of the pixel's center
  double pixel_x = -1.0 + half_d_grid + delta_grid * i_pixel_x;
  double pixel_y = 1.0  - half_d_grid  - delta_grid * i_pixel_y;


  // pixel is outside of the support domain (unit ball)
  if ((pixel_x * pixel_x + pixel_y * pixel_y) > 1.0) {
     return 0;
  }

  // line does not intersect pixel then integral is zero
  if (fabs(pixel_x*normal.x + pixel_y*normal.y - shift) > (sqrt(2.0)*half_d_grid)) {
    return 0.0;
  }

  double integral = 0.0;
  const double tau_step = half_d_grid; // integration step - pixel half-size
  const int nsteps = 2*ngrid + 1;

  int i_step;
  for (i_step = -nsteps-1; i_step < nsteps + 1; ++i_step) {

     double tau = i_step * tau_step;
     vect2d point;

     point.x = in_point.x + tau * direction.x;
     point.y = in_point.y + tau * direction.y;

     vect2d forward_point;
     forward_point.x = point.x + tau_step * direction.x;
     forward_point.y = point.y + tau_step * direction.y;

     // Simpson's rule for the integration along the line
     double value = 0.0;
     double forward_value = 0.0;

     if ((fabs(point.x) < 1.0) && (fabs(point.y) < 1.0)) {
        if ((fabs(point.x-pixel_x) < half_d_grid) && (fabs(point.x-pixel_x) < half_d_grid)
                     && (fabs(point.y-pixel_y) < half_d_grid) && (fabs(point.y-pixel_y) < half_d_grid)) {
                       value = 1.0;
                     }
     }
     if ((fabs(forward_point.x) < 1.0) && (fabs(forward_point.y) < 1.0)) {
        if ((fabs(forward_point.x-pixel_x) < half_d_grid) && (fabs(forward_point.x-pixel_x) < half_d_grid)
                    && (fabs(forward_point.y-pixel_y) < half_d_grid) && (fabs(forward_point.y-pixel_y) < half_d_grid)) {
                      forward_value = 1.0;
                    }
     }
     double weight = cos(2*M_PI*frequency*tau); // weight = cos(2*pi*frequency * tau)
     integral += weight * 0.5 * (value + forward_value) * tau_step; //1D integration along direction (Simpson's integration rule)
  }
  return integral;
}
