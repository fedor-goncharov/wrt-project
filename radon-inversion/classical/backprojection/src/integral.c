#include <math.h>
#include "interp.h"


// computes adjoint Radon transform at (x,y,z)
double radon_adjoint(const double*** radon_values, double x, double y, double z, 
                     const double* phi, 
                     const double* theta, 
                     const double* theta_weights, 
                     const double* shift, 
                     int nphi, int ntheta, int nshift) {

  int i_phi, i_theta;
  double adjoint_integral = 0;
  const double dphi = 2 * M_PI / nphi;
  const double dshift = 2.0 / (nshift - 1);

  // sphere integral
  for (i_phi = 0; i_phi < nphi; ++i_phi) {
    for (i_theta = 0; i_theta < ntheta; ++i_theta) {
      // adjoint s = <(x,y,z), direciton(phi, theta)>
      double s = x*sin(theta[i_theta])*cos(phi[i_phi]) + y*sin(theta[i_theta])*sin(phi[i_phi]) + z*cos(theta[i_theta]);
      // integral term
      adjoint_integral += dphi * theta_weights[i_theta] * radon_interpolate_s(radon_values, s, i_phi, i_theta, nphi, ntheta, nshift, shift);
    }
  }
  return adjoint_integral;
}
