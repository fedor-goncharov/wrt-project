#ifndef __GRID_INTERPOLATION_H
#define __GRID_INTERPOLATION_H

double square_bilinear_interp(double ***function_values, int npixels, int nslices, double radius, 
    int i_slice,   
    double x, double y);

double square_zero_interp(double ***function_values, int npixels, int nslices, double radius, 
    int i_slice, 
    double x, double y);

#endif
