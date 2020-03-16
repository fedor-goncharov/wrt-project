#ifndef _H_RAY_GRID_INTEGRALS
#define _H_RAY_GRID_INTEGRALS

double spect_ray_transform(double*** func_values, 
		     double*** att_values,
		     int ngrid, 
		     double zslice,
		     double shift, 
		     double phi);

double exp_beam_transform(double*** att_values,
			 double px,
			 double py,
			 double pz, 
			 double phi,
			 int ngrid);


#endif