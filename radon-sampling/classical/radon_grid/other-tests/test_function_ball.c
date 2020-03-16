#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// phantom function
// (dcenter_x,y,z) center of the ball
const double dcenter_x = 0.0, dcenter_y = 0.0, dcenter_z = 0.0, 
dradius = 0.45;

double test_function(double x, double y, double z) {
  const double px = x - dcenter_x, py = y - dcenter_y, pz = z - dcenter_z;
  return ((pow(px,2) + pow(py,2) + pow(pz,2)) < pow(dradius,2) ? 1.0 : 0.0);
}

int main(void) {
  
  
  FILE *f;
  f = fopen("centered_ball.dat", "w");
  
  const int ngrid = 129;
  const double delta = 2.0 / (ngrid - 1);
  
  int i, j, k;
  for (i = 0; i < ngrid; ++i) 
    for (j = 0; j < ngrid; ++j)
      for (k = 0; k < ngrid; ++k) {
	  double x = -1.0 + i * delta,
	         y = -1.0 + j * delta,
		 z = -1.0 + k * delta;
		 
	  double value = test_function(x,y,z);
	  fprintf(f, "%lf\n", value);
      }
  fclose(f);
}
