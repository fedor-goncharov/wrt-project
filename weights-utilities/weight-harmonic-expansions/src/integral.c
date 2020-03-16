#include <math.h>
#include <stdbool.h>
#include "interp.h"

#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif 


// exponential ray transform
double exp_ray_transform(double*** att_values,
			 double px,
			 double py,
			 double pz, 
			 double phi,
			 int npixels,
			 double radius, 
			 bool div_form) {
  
  double dir_x = cos((div_form ? phi : phi + M_PI/2.0)),
         dir_y = sin((div_form ? phi : phi + M_PI/2.0)), 
         dir_z = 0.0;
    
  const double half_delta = (2.0*radius / npixels)/2.0;	//half-step of the initial step of the cube grid
  const int nsize = 4 * npixels;	//number of integration steps per ray
  
  double integral = 0;

  int i_v;
  double x = px, y = py, z = pz;
  double step_x = half_delta * dir_x,
         step_y = half_delta * dir_y,
         step_z = half_delta * dir_z;
  
  for (i_v = 0; i_v < (nsize + 1); ++i_v) {
	x += step_x; 
	y += step_y; 
	z += step_z; 	
	double value = cube_trilinear_interp(att_values, npixels, x, y, z);
	integral += value;
  }
  return exp(-integral*half_delta);
}



double wspect_decomposition_real(double* wspect_atpoint, int degree, double* phi, int nphi, bool normalization) {
  double value = 0;
  double norm_factor = normalization ? (1/(2*M_PI)) : 1.0;
  int i_phi;
  for (i_phi = 0; i_phi < nphi; ++i_phi) {
    value += wspect_atpoint[i_phi] * cos(degree * phi[i_phi]);
  }
  return value * (2 * M_PI / nphi) * norm_factor;
}


double wspect_decomposition_imag(double* wspect_atpoint, int degree, double* phi, int nphi, bool normalization) {
  double value = 0;
  double norm_factor = normalization ? (1/(2*M_PI)) : 1.0;
  int i_phi;
  for (i_phi = 0; i_phi < nphi; ++i_phi) {
    value += wspect_atpoint[i_phi] * sin(degree * phi[i_phi]);
  }
  return value * (2 * M_PI / nphi) * norm_factor;
}


