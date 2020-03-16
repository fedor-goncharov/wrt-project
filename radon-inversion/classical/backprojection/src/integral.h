#ifndef __H_INTEGRALS
#define __H_INTEGRALS

double radon_adjoint(const double*** radon_values, double x, double y, double z, 
	                 const double* phi, 
	                 const double* theta, 
	                 const double* theta_weights, 
	                 const double* shift, 
	                 int nphi, int ntheta, int nshift);

#endif