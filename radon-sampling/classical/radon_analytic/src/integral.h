#ifndef _H_INTEGRALS
#define _H_INTEGRALS

double plane_integral(double (*function_ptr)(double, double, double), int ngrid, 
			double phi, double theta, double shift);

#endif
