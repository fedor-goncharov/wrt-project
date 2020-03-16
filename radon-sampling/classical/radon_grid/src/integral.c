#include <math.h>
#include "interp.h"

typedef struct {
   double x;
   double y; 
   double z; 
} vect3d; 

int get_sph_direction(vect3d* vec, double phi, double theta) {
   vec->x = sin(theta) * cos(phi);
   vec->y = sin(theta) * sin(phi);
   vec->z = cos(theta);
   
   return 0;
}

int vec_product(vect3d* vec, vect3d v1, vect3d v2) {
  vec->x = v1.y * v2.z - v1.z * v2.y;
  vec->y = -(v1.x * v2.z - v1.z * v2.x);
  vec->z = v1.x * v2.y - v1.y * v2.x;
  
  return 0;
}

double plane_integral(double*** function_values, int ngrid,
		      double phi, double theta, double shift, double radius) {
  
  vect3d normal;
    get_sph_direction(&normal, phi, theta);
    
  //basis vectors on the plane
  vect3d v1;
    get_sph_direction(&v1, phi, theta + M_PI/2);
  vect3d v2;
    vec_product(&v2, normal, v1);
    
  const double delta = (2.0*radius / ngrid) / 2.0;	//half-step of the initial step of the cube grid
  const int nsize = 4 * ngrid;	        //number of points per slice of the plane

  double integral = 0;
  int i_v1, i_v2;
  double point_x = shift * normal.x,
         point_y = shift * normal.y,
         point_z = shift * normal.z,
         v1_step_x = delta * v1.x,
         v1_step_y = delta * v1.y,
         v1_step_z = delta * v1.z,
         v2_step_x = delta * v2.x,
         v2_step_y = delta * v2.y, 
         v2_step_z = delta * v2.z;
  
  for (i_v1 = -(nsize / 2); i_v1 < (nsize / 2); ++i_v1) {
    point_x += v1_step_x;
    point_y += v1_step_y;
    point_z += v1_step_z;
    
    for (i_v2 = -(nsize / 2); i_v2 < (nsize / 2); ++i_v2) {	
	point_x += v2_step_x;
	point_y += v2_step_y;
	point_z += v2_step_z;
	
        double value = zero_interp(function_values, ngrid,
	                point_x, point_y, point_z);
	
	integral += value; //1D integration along v2
    }
  //1D integration along v1
  }  
  return integral * delta * delta;
}

