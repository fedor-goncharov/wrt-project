#ifndef __H_SLICING_REDUCTION
#define __H_SLICING_REDUCTION


double reduction(double*** ray_values, 
		       double* shift, double* phi, double* theta,
		       int nshift, int nphi, int ntheta, int ngrid,
			int i_shift, int i_phi, int i_theta);

#endif
