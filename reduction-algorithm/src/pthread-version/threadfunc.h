#ifndef _RAY_RADON_REDUCTION_POSIX_H
#define _RAY_RADON_REDUCTION_POSIX_H


struct reduction_threadf_argument {
	int thread_ID;					// id of the thread
	char* output_filename;			// thread's output filename
	int shift_min, shift_max;	    // interval of shifts [shift_min, sift_max-1]
	int nphi, ntheta, nshift, ngrid;
	double ***ray_transforms_values;// values of the test-function on the grid
	double *phi, *theta, *shift; 	// values of angles, shifts
};

void* reduction_ray_radon_threadf(void* thread_argument);	// computation of spect ray integrals

#endif
