#include <stdio.h>
#include "slicing.h"


struct reduction_threadf_argument {
	int thread_ID;
	char* output_filename;
	int shift_min, shift_max;	// interval of z-slices [z_slice_min, z_slice_max-1]
	int nphi, ntheta, nshift, ngrid;
	double ***ray_transforms_values;
	double *phi, *theta, *shift;
};

void* reduction_ray_radon_threadf(void* thread_argument) {


	struct reduction_threadf_argument* pargument = (struct reduction_threadf_argument*)(thread_argument);

	int thread_ID         = pargument -> thread_ID;
	char *output_filename = pargument -> output_filename;

	int	shift_min   = pargument -> shift_min,
	 	  shift_max   = pargument -> shift_max,
	 	  ngrid       = pargument -> ngrid,
	 	  nphi        = pargument -> nphi,
	 	  ntheta      = pargument -> ntheta,
	 	  nshift      = pargument -> nshift;

	double ***ray_transforms_values  = pargument -> ray_transforms_values,
		  *phi           = pargument -> phi,
		  *theta         = pargument -> theta,
		  *shift         = pargument -> shift;

	printf("  POSIX thread %d. Domain of work (shifts): start %d, end %d, size %d\n", thread_ID, shift_min, shift_max-1, shift_max - shift_min);


  //open ouptut file
  FILE *foutput;
	foutput = fopen(output_filename, "w");

	int i_shift, i_phi, i_theta;
    //iterate plane integrals over shift variable
    for (i_shift = shift_min; i_shift < shift_max; ++i_shift) {
        //iterate plane integrals over phi angle
	    for (i_phi = 0; i_phi < nphi; ++i_phi) {
			//iterate plane integrals over theta angle
			  for (i_theta = 0; i_theta < ntheta; ++i_theta) {

		    	    //reduction of ray integrals to Radon transform over the plane (shift[i_shift], phi[i_phi], theta[i_theta])
		    	 double radon_transform = radon_reduction(ray_transforms_values,
		    		                                        shift, phi, theta,
		    		                                        nshift, nphi, ntheta, ngrid,
		    		                                        i_shift, i_phi, i_theta);
		    	 //write result to file
		    	 fprintf(foutput, "%lf, %lf, %lf, %lf\n", shift[i_shift], phi[i_phi], theta[i_theta], radon_transform);
       }
     }
  }

	fclose(foutput);
	printf("  POSIX thread %d. Job is done.\n", thread_ID);
	return NULL;
}
