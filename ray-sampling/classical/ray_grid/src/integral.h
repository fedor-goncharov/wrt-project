#ifndef __RAY_GRID_INTEGRALS_H
#define __RAY_GRID_INTEGRALS_H

double ray_transform(double ***function_values, int npixels, int nslices, double radius, 
		        int i_slice, 
			double shift, double phi);

#endif
