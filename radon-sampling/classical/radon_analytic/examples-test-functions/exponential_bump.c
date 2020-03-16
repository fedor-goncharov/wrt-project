#include <math.h>

#define eps (1e-2)

const double center_x = 0.2, center_y = 0.0, center_z = 0.0;

double test_function(double x, double y, double z) {
  const double xx = x - center_x, yy = y - center_y, zz = z - center_z;
  double value; 
  if ((pow(xx,2) + pow(yy,2) + pow(zz,2)) < (0.5 - eps))
  {
    value = exp((-1)/(0.5 - pow(xx,2) - pow(yy,2) - pow(zz,2)));
  } else {
    value = 0;
  }
  return value;
}
