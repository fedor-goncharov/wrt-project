#include <math.h>

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

double plane_integral(double (*function_ptr)(double, double, double), int ngrid,
			double phi, double theta, double shift) {
  
  vect3d normal;
    get_sph_direction(&normal, phi, theta);
    
  //basis vectors on the plane
  vect3d v1;
    get_sph_direction(&v1, phi, theta + M_PI/2);
  vect3d v2;
    vec_product(&v2, normal, v1);
    
  const double delta = (2.0 / (ngrid - 1))/2.0;	//half-step of the initial step of the cube grid
  const int nsize = 4 * ngrid;	//number of points per slice of the plane
  double integral = 0;
  int i_v1, i_v2;
  
  for (i_v1 = -(nsize / 2); i_v1 < (nsize/2); ++i_v1) {
    for (i_v2 = -(nsize / 2); i_v2 < (nsize/2); ++i_v2) {
        vect3d point;
	point.x = shift * normal.x + delta * v1.x * i_v1 + delta * v2.x * i_v2;  
	point.y = shift * normal.y + delta * v1.y * i_v1 + delta * v2.y * i_v2;
	point.z = shift * normal.z + delta * v1.z * i_v1 + delta * v2.z * i_v2;
	
	double value = (*function_ptr)(point.x, point.y, point.z);
	
	integral += value * delta * delta; //1D integration along v2
    }
    //1D integration along v1
  }  
  return integral;
}

