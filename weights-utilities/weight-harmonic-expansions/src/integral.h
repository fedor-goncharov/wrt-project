#ifndef __H_RAY_GRID_INTEGRALS
#define __H_RAY_GRID_INTEGRALS

double exp_ray_transform(double*** att_values,
			 double x,
			 double y,
			 double z, 
			 double phi,
			 int npixels,
			 double radius, 
			 bool div_form);


double wspect_decomposition_real(double* wspect_atpoint, int degree, double* phi, int nphi, bool normalization);

double wspect_decomposition_imag(double* wspect_atpoint, int degree, double* phi, int nphi, bool normalization);

#endif
