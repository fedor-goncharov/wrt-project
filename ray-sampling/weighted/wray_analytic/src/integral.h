#ifndef _H_RAY_INTEGRAL
#define _H_RAY_INTEGRAL

double ray_integral3d(double (*function_ptr)(double, double, double), 
			double zslice,
			double shift, double phi, int ngrid);

#endif