#ifndef __WRAY_GRID_THREAD_FUNC
#define __WRAY_GRID_THREAD_FUNC


struct spect_threadf_argument {
	int thread_ID;					// id of the thread
	char* output_filename;			// thread's output filename
	int z_slice_min, z_slice_max;	// interval of z-slices [z_slice_min, z_slice_max-1]
	int nphi, nshift, ngrid, nslice; 
	double ***func_values, ***att_values; // values of the test-function on the grid
	double *phi, *shift; 				  // values of angles, shifts
};


void* spect_ray_threadf(void* thread_argument);	// computation of spect ray integrals

#endif