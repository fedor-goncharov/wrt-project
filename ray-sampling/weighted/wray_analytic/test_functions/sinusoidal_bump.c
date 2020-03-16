#include <math.h>

double sine_bump_template(double x, double y, double z,
			      double center_x, double center_y, double center_z,
			      double amplitude, double rin, double rout) {
  
  double r = sqrt(pow(x-center_x,2) + pow(y-center_y,2) + pow(z-center_z,2));
  if (r > rout)
    return 0;
  
  if ((r < rout) && (r > rin)) {
      const double t = M_PI * (r-rout) / (rout - rin);
      return ( (amplitude/(20*M_PI)) * (-sin(6*t)/3 + 3 * sin(4*t) - 15 * sin(2*t) + 20*t ) );
  }
  if (r < rin)
    return amplitude;
  
  return 0;
}


double test_function(double x, double y, double z) {
  //three balls with different amplitudes and 'transitions'
  const double c1_x = 0.0, c1_y = 0.2, c1_z = 0.0, 
	       a1 = 1.0, rin1 = 0.4, rout1 = 0.6,
	       c2_x = -0.5, c2_y = 0.0, c2_z = 0.0, 
	       a2 = 0.6, rin2 = 0.3, rout2 = 0.45, 
	       c3_x = 0.0, c3_y = 0.2, c3_z = 0.5, 
	       a3 = 0.3, rin3 = 0.4, rout3 = 0.6;
	       
  return sine_bump_template(x,y,z, c1_x, c1_y, c1_z, a1, rin1, rout1) + 
	 sine_bump_template(x,y,z, c2_x, c2_y, c2_z, a2, rin2, rout2) + 
	 sine_bump_template(x,y,z, c3_x, c3_y, c3_z, a3, rin3, rout3);
}