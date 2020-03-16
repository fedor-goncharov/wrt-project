#include <math.h>

typedef struct {
   double x;
   double y; 
   double z; 
} vect3d; 


double ray_integral(double (*function_ptr)(double, double, double), 
		        double zslice,
			double shift, double phi, int ngrid) {
  
  vect3d normal;
    normal.x = cos(phi);
    normal.y = sin(phi);
    normal.z = 0;
  
  vect3d direction;
    direction.x = cos(phi + M_PI/2);
    direction.y = sin(phi + M_PI/2);
    direction.z = 0;
    
  //initial point of the ray
  vect3d in_point;
    in_point.x = shift * normal.x;
    in_point.y = shift * normal.y;
    in_point.z = zslice;
    
  const double delta = (2.0 / (ngrid - 1))/2.0;	//half-step of the initial step of the cube grid
  const int nsize = 4 * ngrid;	//number of points per slice of the plane
  
  double integral = 0;
  int i_v;
  
  for (i_v = -(nsize / 2); i_v < (nsize/2) + 1; ++i_v) {
        vect3d point;
	point.x = in_point.x + delta * direction.x * i_v; 
	point.y = in_point.y + delta * direction.y * i_v; 
	point.z = in_point.z + delta * direction.z * i_v; 
	
	double value = (*function_ptr)(point.x, point.y, point.z);
	
	integral += value * delta; //1D integration along direction
  }  
  return integral;
}

