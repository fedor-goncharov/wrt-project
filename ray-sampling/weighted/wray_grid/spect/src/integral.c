#include <math.h>
#include <stdio.h>
#include "interp.h"

typedef struct {
   double x;
   double y; 
   double z; 
} vect3d; 


double exp_beam_transform(double*** att_values,
			 double px,
			 double py,
			 double pz, 
			 double phi,	// phi -- direction of the ray
			 int ngrid) {
  
  vect3d direction;
    direction.x = cos(phi);
    direction.y = sin(phi);
    direction.z = 0;
    
  const double delta = (2.0 / (ngrid - 1))/2.0;	//half-step of the initial step of the cube grid
  const int nstep = 4 * ngrid;	//number of points per slice of the plane
  
  double integral = 0;
  int i_v;
  
  for (i_v = 0; i_v < nstep + 1; ++i_v) {
    double x, y, z;
	  x = px + delta * direction.x * i_v; 
	  y = py + delta * direction.y * i_v; 
	  z = pz + delta * direction.z * i_v; 
	
	  double value = cube_trilinear_interp(att_values, ngrid, x, y, z);
	  integral += value * delta;
  }
  return exp(-integral);
}

double spect_ray_transform(double*** func_values, 
			     double*** att_values,
			     int ngrid, 
			     double zslice,
			     double shift, 
			     double phi) {
  
  vect3d normal;
    normal.x = cos(phi);
    normal.y = sin(phi);
    normal.z = 0;
  
  vect3d direction;
    direction.x = cos(phi + M_PI/2);
    direction.y = sin(phi + M_PI/2);
    direction.z = 0;
    
  //initial point on the ray
  vect3d in_point;
    in_point.x = shift * normal.x;
    in_point.y = shift * normal.y;
    in_point.z = zslice;
  
  const double delta = (2.0 / (ngrid - 1))/2.0;	//half-step of the initial step of the cube grid
  const int nstep = 4 * ngrid;	                //number of points per slice of the plane
  
  double integral = 0;
  int i_v;
  
  for (i_v = -(nstep / 2); i_v < (nstep/2) + 1; ++i_v) {
    vect3d point;
	  point.x = in_point.x + delta * direction.x * i_v; 
	  point.y = in_point.y + delta * direction.y * i_v; 
	  point.z = in_point.z + delta * direction.z * i_v; 
	
	  double value = cube_trilinear_interp(func_values, ngrid, point.x, point.y, point.z);
	  double spect_weight; 

	  if (fabs(value) < 1e-8) { // don't compute weight if value of the function is zero (with double precision) at this point
	    spect_weight = 0;
	  } else {
	    spect_weight = exp_beam_transform(att_values, point.x, point.y, point.z, phi + M_PI/2, ngrid);
	  }
	  integral += spect_weight * value * delta; //1D integration along direction
  }

  return integral;
}

