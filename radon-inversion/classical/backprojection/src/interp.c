#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

double radon_interpolate_s(double*** radon_values, double s, int i_phi, int i_theta, 
    int nphi, int ntheta, int nshift, 
    double* shift) {

  if (fabs(s) > 1) { 
    return 0.0;
  }
  // spline interpolation
  gsl_interp_accel *acc;
  gsl_spline *spline;

  // allocate spline
  acc = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline, nshift);

  gsl_spline_init(spline, shift, radon_values[i_phi][i_theta], nshift);
  double interpolated_value = gsl_spline_eval(spline, s, acc);
  
  // deallocate spline
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);

  return interpolated_value;
}
