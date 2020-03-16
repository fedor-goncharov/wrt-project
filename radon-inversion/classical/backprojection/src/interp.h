#ifndef __H_INTERPOLATION
#define __H_INTERPOLATION

  //interpolate value of Radon transforms at shift coordinate 's' for fixed direction (phi, theta)
  double radon_interpolate_s(double*** radon_values, double s, int i_phi, int i_theta, 
  	int nphi, int ntheta, int nshift, 
  	double* shift);

#endif