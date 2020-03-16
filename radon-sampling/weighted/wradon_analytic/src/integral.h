#ifndef _H_INTEGRALS
#define _H_INTEGRALS

double plane_integral3d(double (*function_ptr)(double, double, double), 
			double (*weight_function_ptr)(double, double, double, double, double),
			int ngrid, 
			double phi, double theta, double shift);

#endif